"""Experimental mask-only QC. Thresholds are engineering heuristics, not clinical."""
from dataclasses import dataclass, field
import nibabel as nib
import networkx as nx
import numpy as np
from scipy.ndimage import map_coordinates
from .cross_sectional_caliber import estimate_tangents


QC_THRESHOLDS = {
    'minimum_path_points': 5, 'minimum_edge_mm': 1e-8,
    'minimum_extent_coverage': .8, 'maximum_terminal_extent_gap': .15,
    'maximum_mean_turn_deg': 60., 'maximum_tangent_step_deg': 45.,
    'maximum_length_extent_simple': 1.8, 'maximum_length_extent_network': 4.,
    'roundoff_tolerance': 1e-8, 'minimum_occupancy': .5,
    'edge_relative_tolerance': 1e-5, 'edge_absolute_tolerance_mm': 1e-8,
}


@dataclass
class Recovery:
    # Node keys are opaque: refined graph nodes need not be integer voxels.
    graph: nx.Graph | None = None
    path: np.ndarray | None = None
    error: str = ''
    supports_cycles: bool = True
    metadata: dict = field(default_factory=dict)


def smoothness(points):
    p = np.asarray(points)
    edges = np.diff(p, axis=0)
    lengths = np.linalg.norm(edges, axis=1)
    if len(edges) < 2 or np.any(lengths <= 1e-10):
        return {'mean_turn_deg': None, 'max_turn_deg': None, 'cumulative_turn_deg': None, 'mean_curvature_per_mm': None, 'length_chord_ratio': None, 'tangent_step_max_deg': None}
    unit = edges / lengths[:, None]
    turns = np.arccos(np.clip(np.sum(unit[:-1] * unit[1:], axis=1), -1, 1))
    tangent = estimate_tangents(p, 6)
    steps = np.degrees(np.arccos(np.clip(np.sum(tangent[:-1] * tangent[1:], axis=1), -1, 1)))
    chord = np.linalg.norm(p[-1] - p[0])
    return {'mean_turn_deg': float(np.degrees(turns).mean()), 'max_turn_deg': float(np.degrees(turns).max()),
            'cumulative_turn_deg': float(np.degrees(turns).sum()), 'mean_curvature_per_mm': float(np.mean(turns / ((lengths[:-1] + lengths[1:]) / 2))),
            'length_chord_ratio': float(lengths.sum() / chord) if chord > 1e-8 else None,
            'tangent_step_max_deg': float(np.nanmax(steps))}


def occupancy(volume, points):
    vox = nib.affines.apply_affine(np.linalg.inv(volume.affine), np.asarray(points))
    return map_coordinates((np.asarray(volume.data) > .5).astype(float), vox.T, order=1, mode='constant', cval=0, prefilter=False)


def quality_control(volume, recovery, mode='network'):
    if mode not in ('simple', 'network'):
        raise ValueError('unknown_qc_mode')
    checks = []
    t = QC_THRESHOLDS
    eps = t['roundoff_tolerance']
    def add(name, value, threshold, passed, critical=True):
        checks.append({'criterion': name, 'value': value, 'threshold': threshold, 'passed': bool(passed), 'critical': critical})
    g, p = recovery.graph, recovery.path
    exists = g is not None and len(g) > 0 and p is not None and len(p) >= 2
    add('extraction_exists', exists, 'true', exists)
    if exists:
        p = np.asarray(p)
        world = np.array([d['world'] for _, d in g.nodes(data=True)])
        finite = np.isfinite(world).all() and np.isfinite(p).all()
        add('finite_world_coordinates', bool(finite), 'true', finite)
        add('connected_components', nx.number_connected_components(g), '==1', nx.is_connected(g))
        add('minimum_path_points', len(p), '>=5', len(p) >= t['minimum_path_points'])
        repeated = len(p) - len(np.unique(p, axis=0))
        add('repeated_path_points', repeated, '==0', repeated == 0)
        lengths = np.array([d.get('length_mm', np.nan) for _, _, d in g.edges(data=True)])
        positive = bool(len(lengths) and np.isfinite(lengths).all() and np.all(lengths > t['minimum_edge_mm']))
        add('positive_finite_edges', positive, 'true; >1e-8 mm', positive)
        cycles = g.number_of_edges() - len(g) + nx.number_connected_components(g)
        add('simple_mode_cycles', cycles, '==0 in simple mode', mode != 'simple' or cycles == 0)
        if finite and positive and repeated == 0:
            lookup = {tuple(d['world']): n for n, d in g.nodes(data=True)}
            nodes = [lookup.get(tuple(point)) for point in p]
            ordered = all(n is not None for n in nodes) and all(g.has_edge(a, b) for a, b in zip(nodes[:-1], nodes[1:]))
            consistent = all(np.isclose(data['length_mm'], np.linalg.norm(np.asarray(g.nodes[a]['world']) - g.nodes[b]['world']), rtol=t['edge_relative_tolerance'], atol=t['edge_absolute_tolerance_mm']) for a, b, data in g.edges(data=True))
            add('ordered_graph_path', ordered, 'consecutive graph edges', ordered)
            add('physical_edge_consistency', consistent, 'allclose Euclidean world distance', consistent)
            foreground = nib.affines.apply_affine(volume.affine, np.argwhere(volume.data > .5))
            center = foreground.mean(axis=0)
            _, _, axes = np.linalg.svd(foreground - center, full_matrices=False)
            axis = axes[0]
            mask_u, graph_u = (foreground - center) @ axis, (world - center) @ axis
            extent = float(np.ptp(mask_u))
            coverage = float(np.ptp(graph_u) / extent) if extent else 0.
            gap = float(max(graph_u.min() - mask_u.min(), mask_u.max() - graph_u.max(), 0) / extent) if extent else 1.
            add('mask_extent_coverage', coverage, '>=0.80', coverage >= t['minimum_extent_coverage'] - eps)
            add('terminal_extent_gap', gap, '<=0.15 of mask PCA extent', gap <= t['maximum_terminal_extent_gap'] + eps)
            diagnostic = smoothness(p)
            add('mean_turn_deg', diagnostic['mean_turn_deg'], '<=60 degrees', diagnostic['mean_turn_deg'] is not None and diagnostic['mean_turn_deg'] <= t['maximum_mean_turn_deg'] + eps)
            add('tangent_step_max_deg', diagnostic['tangent_step_max_deg'], '<=45 degrees', diagnostic['tangent_step_max_deg'] is not None and diagnostic['tangent_step_max_deg'] <= t['maximum_tangent_step_deg'] + eps)
            ratio = float(np.linalg.norm(np.diff(p, axis=0), axis=1).sum() / extent) if extent else float('inf')
            add('length_extent_ratio', ratio, '<=1.8 simple; <=4 network', ratio <= t['maximum_length_extent_simple' if mode == 'simple' else 'maximum_length_extent_network'] + eps)
            # Diagnostic for raw voxel edges as well as refined paths; a sampled
            # interpolated-mask test is not proof of continuous lumen containment.
            outside = int(np.sum(occupancy(volume, p) < t['minimum_occupancy'] - eps))
            add('outside_path_samples', outside, '==0', outside == 0)
        add('cycle_capability', recovery.supports_cycles, 'reported; not a selection criterion', recovery.supports_cycles, critical=False)
    failed = [c['criterion'] for c in checks if c['critical'] and not c['passed']]
    return {'qc_pass': not failed, 'qc_fail_count': len(failed), 'qc_fail_reasons': ';'.join(failed), 'checks': checks}
