"""Experimental branch-wise physical refinement; no analytic truth dependencies."""
from dataclasses import dataclass, asdict
import networkx as nx
import nibabel as nib
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree
from .centerline_qc import occupancy, Recovery
from .cross_sectional_caliber import estimate_tangents


@dataclass(frozen=True)
class RefinementParameters:
    recenter: bool = True
    spacing_mm: float = .5
    smoothing_scale_mm: float = 1.0
    refine_endpoints: bool = True
    slab_sigma_mm: float = .5
    junction_guard_mm: float = 3.0
    max_displacement_voxels: float = 2.0


def distance(points):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]


def resample(points, spacing):
    if spacing <= 0 or not np.isfinite(spacing):
        raise ValueError('invalid_resampling_spacing')
    points = np.asarray(points, dtype=float)
    d = distance(points)
    if len(points) < 2 or not np.isfinite(points).all() or np.any(np.diff(d) <= 1e-10):
        raise ValueError('invalid_path_for_resampling')
    # Equal intervals no greater than requested spacing; preserve both ends.
    positions = np.linspace(0, d[-1], int(np.ceil(d[-1] / spacing)) + 1)
    return np.column_stack([np.interp(positions, d, points[:, i]) for i in range(3)])


def chains(graph):
    """Split graph at degree != 2, retaining junction identity. Cycles explicit."""
    if not nx.is_tree(graph):
        raise ValueError('refinement_requires_tree; supports_cycles=false')
    visited = set()
    for node in graph:
        if graph.degree[node] == 2:
            continue
        for neighbor in graph.neighbors(node):
            edge = frozenset((node, neighbor))
            if edge in visited:
                continue
            visited.add(edge)
            chain = [node, neighbor]
            previous, current = node, neighbor
            while graph.degree[current] == 2:
                nxt = next(n for n in graph.neighbors(current) if n != previous)
                visited.add(frozenset((current, nxt)))
                chain.append(nxt)
                previous, current = current, nxt
            yield chain


def principal_points(graph):
    endpoints = [n for n in graph if graph.degree[n] == 1]
    best, chosen = -1, None
    for i, a in enumerate(endpoints):
        lengths, paths = nx.single_source_dijkstra(graph, a, weight='length_mm')
        for b in endpoints[i + 1:]:
            if lengths[b] > best:
                best, chosen = lengths[b], paths[b]
    if chosen is None:
        raise ValueError('no_principal_tree_path')
    return np.array([graph.nodes[n]['world'] for n in chosen])


class Refiner:
    def __init__(self, volume):
        self.volume = volume
        voxels = np.argwhere(volume.data > .5)
        self.foreground = nib.affines.apply_affine(volume.affine, voxels)
        self.tree = cKDTree(self.foreground)
        edt = distance_transform_edt(np.pad(volume.data > .5, 1), sampling=volume.spacing)[1:-1, 1:-1, 1:-1]
        self.radii = edt[tuple(voxels.T)]

    def constrain(self, old, candidate, limit):
        delta = candidate - old
        norm = np.linalg.norm(delta)
        if norm > limit:
            candidate = old + delta * limit / norm
        # Backtrack toward the already valid point if a proposal leaves the mask.
        for factor in (1., .5, .25, .125, 0.):
            point = old + factor * (candidate - old)
            if occupancy(self.volume, point[None])[0] >= .5 - 1e-8:
                return point
        raise ValueError('raw_point_outside_interpolated_lumen')

    def refine_chain(self, points, terminals, parameters):
        raw = resample(points, parameters.spacing_mm)
        p = raw.copy()
        d = distance(p)
        fixed = np.zeros(len(p), dtype=bool)
        if not terminals[0]:
            fixed |= d <= parameters.junction_guard_mm
        if not terminals[1]:
            fixed |= d[-1] - d <= parameters.junction_guard_mm
        limit = parameters.max_displacement_voxels * max(self.volume.spacing)
        if parameters.recenter:
            tangents = estimate_tangents(p, 6.)
            for i, (point, tangent) in enumerate(zip(p, tangents)):
                if fixed[i]:
                    continue
                _, nearest = self.tree.query(point)
                radius = 2 * self.radii[nearest] + max(self.volume.spacing)
                ids = self.tree.query_ball_point(point, radius)
                offsets = self.foreground[ids] - point
                axial = offsets @ tangent
                weights = self.radii[ids] * np.exp(-.5 * (axial / parameters.slab_sigma_mm) ** 2)
                centroid = np.average(offsets, axis=0, weights=weights)
                transverse = centroid - (centroid @ tangent) * tangent
                p[i] = self.constrain(point, point + transverse, limit)
        if parameters.smoothing_scale_mm > 0:
            s = distance(p)
            smoothed = p.copy()
            for i in range(len(p)):
                if fixed[i]:
                    continue
                offset = s - s[i]
                chosen = np.abs(offset) <= 2 * parameters.smoothing_scale_mm + 1e-10
                if chosen.sum() < 3:
                    continue
                x = offset[chosen]
                w = np.exp(-.5 * (x / parameters.smoothing_scale_mm) ** 2)
                design = np.column_stack((np.ones(len(x)), x, x * x))
                coefficient = np.linalg.lstsq(design * np.sqrt(w[:, None]), p[chosen] * np.sqrt(w[:, None]), rcond=None)[0]
                smoothed[i] = self.constrain(raw[i], coefficient[0], limit)
            p = smoothed
        if parameters.refine_endpoints:
            tangents = estimate_tangents(p, 6.)
            for i, sign, terminal in ((0, -1, terminals[0]), (-1, 1, terminals[1])):
                if not terminal:
                    continue
                origin = p[i].copy()
                # Stop at the first mask boundary, not at a later foreground island.
                for step in np.arange(.1, 2 * max(self.volume.spacing) + .001, .1):
                    candidate = origin + sign * step * tangents[i]
                    if occupancy(self.volume, candidate[None])[0] < .5:
                        break
                    p[i] = candidate
        # Reparameterize after geometric changes; junction positions stay fixed.
        p = resample(p, parameters.spacing_mm)
        if np.any(occupancy(self.volume, p) < .5 - 1e-8):
            raise ValueError('resampled_point_outside_lumen')
        # Test segments densely, as endpoints alone can miss excursions.
        dense = resample(p, min(self.volume.spacing) / 4)
        if np.any(occupancy(self.volume, dense) < .5 - 1e-8):
            raise ValueError('refined_segment_outside_lumen')
        return p, raw

    def refine(self, recovery, parameters=RefinementParameters()):
        if not all(np.isfinite(x) for x in (parameters.spacing_mm, parameters.smoothing_scale_mm, parameters.slab_sigma_mm, parameters.junction_guard_mm, parameters.max_displacement_voxels)):
            raise ValueError('nonfinite_refinement_parameters')
        if parameters.smoothing_scale_mm < 0 or parameters.slab_sigma_mm <= 0 or parameters.junction_guard_mm < 0 or parameters.max_displacement_voxels <= 0:
            raise ValueError('invalid_refinement_parameters')
        graph = nx.Graph()
        displacements = []
        raw_graph = recovery.graph
        for chain_id, chain in enumerate(chains(raw_graph)):
            points = np.array([raw_graph.nodes[n]['world'] for n in chain])
            terminals = (raw_graph.degree[chain[0]] == 1, raw_graph.degree[chain[-1]] == 1)
            refined, raw = self.refine_chain(points, terminals, parameters)
            # Measure correspondence by normalized raw arc fraction; endpoint
            # extension changes length, so this is not nearest-neighbor motion.
            u = distance(refined) / distance(refined)[-1]
            raw_u = distance(raw) / distance(raw)[-1]
            counterpart = np.column_stack([np.interp(u, raw_u, raw[:, j]) for j in range(3)])
            displacements.extend(np.linalg.norm(refined - counterpart, axis=1).tolist())
            keys = [('anchor', chain[0])] + [('sample', chain_id, i) for i in range(1, len(refined) - 1)] + [('anchor', chain[-1])]
            for key, point in zip(keys, refined):
                graph.add_node(key, world=tuple(point))
            for a, b in zip(keys[:-1], keys[1:]):
                length = np.linalg.norm(np.asarray(graph.nodes[a]['world']) - graph.nodes[b]['world'])
                if length <= 1e-8:
                    raise ValueError('duplicate_refined_points')
                graph.add_edge(a, b, length_mm=float(length))
        if not nx.is_tree(graph):
            raise ValueError('refinement_topology_changed')
        metadata = {'parameters': asdict(parameters), 'mean_displacement_mm': float(np.mean(displacements)),
                    'max_displacement_mm': float(np.max(displacements)), 'points_outside_lumen_count': 0,
                    'displacements_mm': displacements}
        return Recovery(graph, principal_points(graph), supports_cycles=False, metadata=metadata)
