"""Unchanged baseline adapter and coordinate-based, coverage-aware scoring."""
import networkx as nx
import numpy as np
from scipy.optimize import linear_sum_assignment

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.cross_sectional_caliber import estimate_tangents
from neurovasc.geometry.geodesic_centerline import extract_geodesic_tree
from neurovasc.graph.vascular_graph import build_vascular_graph
from .truth import project


def lee(volume):
    return build_vascular_graph(extract_centerline(volume))


ESTIMATORS = {'lee': lee, 'alternative': extract_geodesic_tree}


def topology(graph):
    junctions = [n for n, d in graph.degree if d >= 3]
    groups = list(nx.connected_components(graph.subgraph(junctions)))
    components = nx.number_connected_components(graph)
    return {'connected_components': components, 'endpoint_count': sum(d == 1 for _, d in graph.degree),
            'junction_count': len(groups), 'junction_voxel_count': len(junctions),
            'cycle_count': graph.number_of_edges() - graph.number_of_nodes() + components}, groups


def validate_path(graph, path):
    points = np.asarray(path.world_points)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 2 or not np.isfinite(points).all():
        raise ValueError('malformed_path_coordinates')
    for _, data in graph.nodes(data=True):
        if np.shape(data.get('world')) != (3,) or not np.isfinite(data['world']).all():
            raise ValueError('malformed_graph_coordinates')
    for a, b, data in graph.edges(data=True):
        actual = np.linalg.norm(np.asarray(graph.nodes[a]['world']) - graph.nodes[b]['world'])
        if not actual > 0 or not np.isclose(data.get('length_mm', np.nan), actual):
            raise ValueError('malformed_graph_edge_length')
    nodes = [tuple(p) for p in path.voxel_points]
    if len(nodes) != len(points) or len(set(nodes)) != len(nodes):
        raise ValueError('malformed_path_nodes')
    if any(n not in graph for n in nodes) or any(not graph.has_edge(a, b) for a, b in zip(nodes[:-1], nodes[1:])):
        raise ValueError('discontinuous_path')
    if not np.allclose(points, [graph.nodes[n]['world'] for n in nodes]):
        raise ValueError('path_graph_coordinate_mismatch')
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    if np.any(steps <= 0) or np.shape(path.distance_mm) != (len(points),) or not np.allclose(path.distance_mm, np.r_[0, np.cumsum(steps)]):
        raise ValueError('malformed_path_distance')
    return True


def angular_errors(estimated, truth):
    estimated, truth = np.asarray(estimated), np.asarray(truth)
    norm = np.linalg.norm(estimated, axis=1) * np.linalg.norm(truth, axis=1)
    cosine = np.divide(np.abs(np.sum(estimated * truth, axis=1)), norm, out=np.full(len(norm), np.nan), where=norm > 0)
    return np.degrees(np.arccos(np.clip(cosine, 0, 1)))


def score_points(vessel, points, branch_name=None):
    mapped, tangents, position, labels = project(vessel, points, branch_name)
    route_names = None
    if vessel.geometry == 'bifurcation' and not branch_name:
        # Match the two endpoint branches first. Near the junction, assigning
        # samples to the unused third branch would artificially inflate coverage.
        route_names = [name for name in vessel.branches if name in (labels[0], labels[-1])]
        if len(route_names) == 2:
            mapped, tangents, position, labels = project(vessel, points, route_names)
    errors = np.linalg.norm(points - mapped, axis=1)
    angles = angular_errors(estimate_tangents(points, window_mm=6.0), tangents)
    if vessel.geometry == 'bifurcation':
        # The analytic tangent is undefined at the junction. Exclude the full
        # half-window around it so branch-local fits don't straddle the split.
        angles[np.linalg.norm(mapped, axis=1) <= 3.0] = np.nan
    finite = angles[np.isfinite(angles)]
    length = float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())
    true_length = float(np.linalg.norm(np.diff(vessel.branches[branch_name], axis=0))) if branch_name else vessel.length_mm
    coverage = float(np.ptp(position))
    monotonic = bool(np.all(np.diff(position) >= -1e-8) or np.all(np.diff(position) <= 1e-8))
    start, end = float(position[0]), float(position[-1])
    if vessel.geometry == 'bifurcation' and not branch_name:
        if len(route_names) == 2:
            first, second = route_names
            first_length = np.linalg.norm(vessel.branches[first][1] - vessel.branches[first][0])
            second_length = np.linalg.norm(vessel.branches[second][1] - vessel.branches[second][0])
            route = np.where(np.array(labels) == first, (1 - position) * first_length, first_length + position * second_length) / true_length
            coverage = float(np.ptp(route))
            start, end = float(route[0]), float(route[-1])
            monotonic = bool(np.all(np.diff(route) >= -1e-8) or np.all(np.diff(route) <= 1e-8))
        else:
            coverage, start, end, monotonic = 0.0, None, None, False
    metrics = {'true_length_mm': true_length, 'recovered_length_mm': length,
               'length_error_mm': abs(length - true_length), 'signed_length_error_mm': length - true_length,
               'length_error_percent': 100 * abs(length - true_length) / true_length,
               'coverage_fraction': coverage, 'normalized_start': start, 'normalized_end': end,
               'monotonic_order': monotonic, 'mean_position_error_mm': float(errors.mean()),
               'median_position_error_mm': float(np.median(errors)), 'max_position_error_mm': float(errors.max()),
               'position_rmse_mm': float(np.sqrt(np.mean(errors ** 2))), 'tangent_valid_count': len(finite),
               **{f'{stat}_tangent_error_deg': float(fn(finite)) if len(finite) else None for stat, fn in [('mean', np.mean), ('median', np.median), ('max', np.max)]}}
    samples = [{'sample_index': i, 'x_mm': float(p[0]), 'y_mm': float(p[1]), 'z_mm': float(p[2]),
                'truth_x_mm': float(q[0]), 'truth_y_mm': float(q[1]), 'truth_z_mm': float(q[2]),
                'analytic_position': float(u), 'truth_branch': label, 'position_error_mm': float(e),
                'tangent_error_deg': float(a) if np.isfinite(a) else None}
               for i, (p, q, u, label, e, a) in enumerate(zip(points, mapped, position, labels, errors, angles))]
    return metrics, samples


def branch_scores(vessel, graph, groups):
    names = list(vessel.branches)
    endpoints = [n for n, d in graph.degree if d == 1]
    matches = {}
    if endpoints:
        costs = np.array([[np.linalg.norm(np.asarray(graph.nodes[n]['world']) - vessel.branches[name][1]) for name in names] for n in endpoints])
        a, b = linear_sum_assignment(costs)
        matches = {names[j]: (endpoints[i], float(costs[i, j])) for i, j in zip(a, b)}
    rows, samples = [], []
    for name in names:
        row = {'branch': name, 'branch_path_found': False, 'branch_recovered': False, 'failure_reason': 'missing_endpoint_or_unique_junction'}
        if name in matches and len(groups) == 1:
            endpoint, endpoint_error = matches[name]
            _, paths = nx.multi_source_dijkstra(graph, groups[0], weight='length_mm')
            points = np.array([graph.nodes[n]['world'] for n in paths[endpoint]])
            if len(points) >= 2:
                metrics, profile = score_points(vessel, points, name)
                quality = bool(metrics['coverage_fraction'] >= 0.8 and endpoint_error <= 2 * max(vessel.volume.spacing))
                row.update(metrics, branch_path_found=True, branch_recovered=quality, endpoint_error_mm=endpoint_error,
                           failure_reason='' if quality else 'branch_coverage_or_endpoint_tolerance_not_met')
                samples.extend({**p, 'path_kind': name} for p in profile)
        rows.append(row)
    return rows, samples


def evaluate(vessel, estimator_name):
    row = {'extraction_success': False, 'status': 'failed', 'failure_reason': '', 'path_found': False,
           'recovered_sample_count': 0, 'path_continuity': False, 'correct_principal_path': False,
           'true_length_mm': vessel.length_mm}
    graph = path = None
    try:
        vascular = ESTIMATORS[estimator_name](vessel.volume)
        graph = vascular.graph
        stats, groups = topology(graph)
        row.update(stats)
        path = extract_main_path(vascular)
        row['path_found'] = True
        validate_path(graph, path)
        if stats['connected_components'] != 1:
            raise ValueError('disconnected_recovered_graph')
    except (ValueError, nx.NetworkXException) as error:
        row['failure_reason'] = str(error)
        branches = [{'branch': name, 'branch_recovered': False, 'failure_reason': str(error)} for name in vessel.branches] if vessel.geometry == 'bifurcation' else []
        return row, [], branches, graph, path
    metrics, samples = score_points(vessel, path.world_points)
    row.update(metrics, extraction_success=True, status='valid', recovered_sample_count=len(path.world_points), graph_sample_count=len(graph), path_continuity=True)
    expected_endpoints = 3 if vessel.geometry == 'bifurcation' else 2
    topology_correct = stats['endpoint_count'] == expected_endpoints and stats['junction_count'] == (expected_endpoints - 2) and stats['cycle_count'] == 0
    row['topology_correct'] = topology_correct
    branches = []
    if vessel.geometry == 'bifurcation':
        branches, branch_profiles = branch_scores(vessel, graph, groups)
        labels = [samples[0]['truth_branch'], samples[-1]['truth_branch']]
        correct = topology_correct and 'parent' in labels and len(set(labels)) == 2 and metrics['monotonic_order'] and metrics['coverage_fraction'] >= .8
        row['junction_error_mm'] = float(np.linalg.norm(np.mean([graph.nodes[n]['world'] for n in groups[0]], axis=0))) if len(groups) == 1 else None
        row['network_true_length_mm'] = sum(np.linalg.norm(segment[1] - segment[0]) for segment in vessel.branches.values())
        row['network_recovered_length_mm'] = sum(d['length_mm'] for _, _, d in graph.edges(data=True))
        row['all_branches_recovered'] = all(b['branch_recovered'] for b in branches)
    else:
        branch_profiles = []
        correct = topology_correct and metrics['monotonic_order'] and min(metrics['normalized_start'], metrics['normalized_end']) <= .1 and max(metrics['normalized_start'], metrics['normalized_end']) >= .9
    row['correct_principal_path'] = bool(correct)
    return row, [{**p, 'path_kind': 'main'} for p in samples] + branch_profiles, branches, graph, path
