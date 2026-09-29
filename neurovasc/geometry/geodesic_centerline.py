"""Experimental clearance-weighted geodesic tree; no skeletonization or truth.

Full 26-neighbor foreground graph with physical edge lengths. Disconnected
masks are rejected rather than silently selecting a favorable component.
"""
from dataclasses import dataclass
from itertools import product

import networkx as nx
import nibabel as nib
import numpy as np
from scipy.ndimage import distance_transform_edt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra

from neurovasc.graph.vascular_graph import VascularGraph
from neurovasc.io.volume import Volume


@dataclass(frozen=True)
class GeodesicParameters:
    radius_power: float = 2.0
    terminal_band_voxels: float = 1.0
    minimum_branch_length_mm: float = 3.0
    branch_radius_factor: float = 2.5
    maximum_branches: int = 8
    maximum_foreground_voxels: int = 200_000


def _foreground_graph(volume):
    mask = np.asarray(volume.data) > 0.5
    if mask.ndim != 3 or not mask.any():
        raise ValueError('empty_or_non3d_mask')
    affine = np.asarray(volume.affine, dtype=float)
    if affine.shape != (4, 4) or not np.isfinite(affine).all():
        raise ValueError('invalid_affine')
    axes = affine[:3, :3]
    lengths = np.linalg.norm(axes, axis=0)
    if np.any(lengths <= 0) or not np.allclose(lengths, volume.spacing):
        raise ValueError('affine_spacing_mismatch')
    if not np.allclose(axes.T @ axes, np.diag(lengths ** 2), atol=1e-8):
        raise ValueError('sheared_affine_not_supported')
    voxels = np.argwhere(mask)
    ids = np.full(mask.shape, -1, dtype=np.int32)
    ids[tuple(voxels.T)] = np.arange(len(voxels))
    sources, destinations, weights = [], [], []
    shape = np.asarray(mask.shape)
    for offset in product((-1, 0, 1), repeat=3):
        if offset <= (0, 0, 0):
            continue
        neighbors = voxels + offset
        inside = np.all((neighbors >= 0) & (neighbors < shape), axis=1)
        a = np.flatnonzero(inside)
        b = ids[tuple(neighbors[inside].T)]
        keep = b >= 0
        a, b = a[keep], b[keep]
        length = np.linalg.norm(np.asarray(volume.spacing) * offset)
        sources.extend((a, b))
        destinations.extend((b, a))
        weights.extend((np.full(len(a), length), np.full(len(a), length)))
    graph = csr_matrix((np.concatenate(weights), (np.concatenate(sources), np.concatenate(destinations))), shape=(len(voxels), len(voxels)))
    graph.sort_indices()
    if connected_components(graph, directed=False, return_labels=False) != 1:
        raise ValueError('disconnected_mask')
    return voxels, nib.affines.apply_affine(affine, voxels), graph


def _terminal(distances, physical_graph, points, radii, band_mm):
    """Recenter the farthest geodesic shell's connected terminal patch.

    Ties pick the first voxel in lexicographic index order. A patch centroid is
    weighted by physical EDT radius, then snapped to the nearest patch voxel.
    """
    farthest = int(np.argmax(distances))
    candidates = np.flatnonzero(distances >= distances[farthest] - band_mm)
    _, labels = connected_components(physical_graph[candidates][:, candidates], directed=False)
    component = labels[np.flatnonzero(candidates == farthest)[0]]
    candidates = candidates[labels == component]
    center = np.average(points[candidates], axis=0, weights=radii[candidates])
    squared = np.sum((points[candidates] - center) ** 2, axis=1)
    # Resolve numerical ties in voxel order, including after rigid affine rotation.
    chosen = np.flatnonzero(np.isclose(squared, squared.min(), rtol=1e-12, atol=1e-12))[0]
    return int(candidates[chosen])


def _trace(predecessors, start, targets):
    chain, visited = [], set()
    node = int(start)
    while node not in targets:
        if node < 0 or node in visited:
            raise ValueError('invalid_geodesic_predecessor_chain')
        visited.add(node)
        chain.append(node)
        node = int(predecessors[node])
    return chain + [node]


def extract_geodesic_tree(volume: Volume, parameters=GeodesicParameters()) -> VascularGraph:
    if not all(np.isfinite(x) and x > 0 for x in (parameters.radius_power, parameters.terminal_band_voxels, parameters.minimum_branch_length_mm, parameters.branch_radius_factor)):
        raise ValueError('invalid_geodesic_parameters')
    if any(not isinstance(x, int) or isinstance(x, bool) or x < 1 for x in (parameters.maximum_branches, parameters.maximum_foreground_voxels)):
        raise ValueError('invalid_geodesic_limits')
    count = int(np.count_nonzero(np.asarray(volume.data) > 0.5))
    if count > parameters.maximum_foreground_voxels:
        raise ValueError('foreground_graph_size_limit')
    voxels, points, physical = _foreground_graph(volume)
    if len(voxels) < 2:
        raise ValueError('insufficient_foreground_voxels')
    # Padding ensures the EDT always sees external background, even at array faces.
    radius_map = distance_transform_edt(np.pad(np.asarray(volume.data) > 0.5, 1), sampling=volume.spacing)[1:-1, 1:-1, 1:-1]
    radii = radius_map[tuple(voxels.T)]
    coo = physical.tocoo()
    cost = coo.data * (radii.max() / np.sqrt(radii[coo.row] * radii[coo.col])) ** parameters.radius_power
    weighted = csr_matrix((cost, (coo.row, coo.col)), shape=physical.shape)
    seed = int(np.argmax(radii))
    band = max(volume.spacing) * parameters.terminal_band_voxels
    distance = dijkstra(physical, directed=False, indices=seed)
    a = _terminal(distance, physical, points, radii, band)
    distance = dijkstra(physical, directed=False, indices=a)
    b = _terminal(distance, physical, points, radii, band)
    if a == b:
        raise ValueError('coincident_endpoints')
    _, predecessors = dijkstra(weighted, directed=False, indices=a, return_predecessors=True)
    chain = _trace(predecessors, b, {a})
    tree = nx.Graph()

    def add_chain(chain):
        for node in chain:
            key = tuple(map(int, voxels[node]))
            tree.add_node(key, voxel=key, world=tuple(map(float, points[node])), foreground_id=node)
        for x, y in zip(chain[:-1], chain[1:]):
            tree.add_edge(tuple(voxels[x]), tuple(voxels[y]), length_mm=float(np.linalg.norm(points[x] - points[y])))

    add_chain(chain)
    selected = set(chain)
    threshold = max(parameters.minimum_branch_length_mm, parameters.branch_radius_factor * radii.max())
    for branch_index in range(parameters.maximum_branches + 1):
        residual = dijkstra(physical, directed=False, indices=sorted(selected), min_only=True)
        if residual.max() <= threshold:
            break
        if branch_index == parameters.maximum_branches:
            raise ValueError('branch_limit_exceeded')
        endpoint = _terminal(residual, physical, points, radii, band)
        _, predecessors, _ = dijkstra(weighted, directed=False, indices=sorted(selected), min_only=True, return_predecessors=True)
        branch = _trace(predecessors, endpoint, selected)
        if len(branch) < 2:
            raise ValueError('branch_growth_stalled')
        add_chain(branch)
        selected.update(branch)
    if not nx.is_tree(tree):
        raise ValueError('invalid_recovered_tree')
    return VascularGraph(tree)
