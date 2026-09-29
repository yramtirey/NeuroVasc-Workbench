from dataclasses import dataclass
from collections import deque

import networkx as nx
import numpy as np

from neurovasc.graph.vascular_graph import VascularGraph


@dataclass
class Branch:
    branch_id: str
    start_node: str
    end_node: str
    voxel_path: list[tuple[int, int, int]]
    world_path: np.ndarray
    length_mm: float
    chord_length_mm: float
    tortuosity: float


@dataclass
class BranchGraph:
    logical_graph: nx.Graph
    branches: list[Branch]

    @property
    def n_logical_nodes(self) -> int:
        return self.logical_graph.number_of_nodes()

    @property
    def n_logical_edges(self) -> int:
        return self.logical_graph.number_of_edges()

    @property
    def n_branches(self) -> int:
        return len(self.branches)

    @property
    def junction_nodes(self) -> list[str]:
        return [
            node
            for node, data in self.logical_graph.nodes(data=True)
            if data["kind"] == "junction"
        ]

    @property
    def endpoint_nodes(self) -> list[str]:
        return [
            node
            for node, data in self.logical_graph.nodes(data=True)
            if data["kind"] == "endpoint"
        ]


def _connected_components_from_subset(
    graph: nx.Graph,
    subset_nodes: set,
) -> list[set]:
    """
    Connected components induced by a subset of nodes.
    """
    subgraph = graph.subgraph(subset_nodes)
    return [set(component) for component in nx.connected_components(subgraph)]


def _cluster_centroid_world(
    graph: nx.Graph,
    cluster: set,
) -> tuple[float, float, float]:
    points = np.array(
        [graph.nodes[node]["world"] for node in cluster],
        dtype=float,
    )
    centroid = points.mean(axis=0)
    return tuple(float(x) for x in centroid)


def build_branch_graph(
    vascular_graph: VascularGraph,
) -> BranchGraph:
    """
    Convert the raw voxel-level vascular graph into a branch-level graph.

    Logical nodes are:
      - endpoint voxels
      - clusters of touching junction voxels

    Logical edges are centerline paths between logical nodes.
    """

    raw_graph = vascular_graph.graph

    endpoint_set = set(vascular_graph.endpoints)
    junction_set = set(vascular_graph.branchpoints)

    # Cluster adjacent junction voxels into logical junctions
    junction_clusters = _connected_components_from_subset(
        raw_graph,
        junction_set,
    )

    logical_graph = nx.Graph()

    voxel_to_logical: dict[tuple[int, int, int], str] = {}

    # Add logical junction nodes
    for idx, cluster in enumerate(junction_clusters, start=1):
        node_id = f"J{idx}"

        logical_graph.add_node(
            node_id,
            kind="junction",
            voxels=sorted(cluster),
            world=_cluster_centroid_world(raw_graph, cluster),
        )

        for voxel in cluster:
            voxel_to_logical[voxel] = node_id

    # Add logical endpoint nodes
    for idx, voxel in enumerate(sorted(endpoint_set), start=1):
        node_id = f"E{idx}"

        logical_graph.add_node(
            node_id,
            kind="endpoint",
            voxels=[voxel],
            world=raw_graph.nodes[voxel]["world"],
        )

        voxel_to_logical[voxel] = node_id

    logical_voxels = set(voxel_to_logical.keys())

    branches: list[Branch] = []
    visited_edges = set()

    branch_counter = 1

    for start_voxel in logical_voxels:
        for neighbor in raw_graph.neighbors(start_voxel):
            edge_key = frozenset((start_voxel, neighbor))

            if edge_key in visited_edges:
                continue

            path = [start_voxel]
            prev = start_voxel
            curr = neighbor

            visited_edges.add(edge_key)

            while True:
                path.append(curr)

                if curr in logical_voxels:
                    break

                next_candidates = [
                    n for n in raw_graph.neighbors(curr)
                    if n != prev
                ]

                if len(next_candidates) == 0:
                    break

                if len(next_candidates) > 1:
                    # This should be rare after junction clustering,
                    # but if it happens, stop the branch here.
                    break

                nxt = next_candidates[0]
                visited_edges.add(frozenset((curr, nxt)))
                prev, curr = curr, nxt

            end_voxel = path[-1]

            if end_voxel == start_voxel:
                continue

            if end_voxel not in logical_voxels:
                continue

            start_node = voxel_to_logical[start_voxel]
            end_node = voxel_to_logical[end_voxel]

            if start_node == end_node:
                continue

            world_path = np.array(
                [raw_graph.nodes[v]["world"] for v in path],
                dtype=float,
            )

            segment_lengths = np.linalg.norm(
                np.diff(world_path, axis=0),
                axis=1,
            )

            length_mm = float(segment_lengths.sum())

            chord_length_mm = float(
                np.linalg.norm(world_path[-1] - world_path[0])
            )

            tortuosity = (
                length_mm / chord_length_mm
                if chord_length_mm > 0
                else 1.0
            )

            branch_id = f"B{branch_counter:03d}"
            branch_counter += 1

            branch = Branch(
                branch_id=branch_id,
                start_node=start_node,
                end_node=end_node,
                voxel_path=path,
                world_path=world_path,
                length_mm=length_mm,
                chord_length_mm=chord_length_mm,
                tortuosity=tortuosity,
            )

            branches.append(branch)

            logical_graph.add_edge(
                start_node,
                end_node,
                branch_id=branch.branch_id,
                length_mm=branch.length_mm,
                chord_length_mm=branch.chord_length_mm,
                tortuosity=branch.tortuosity,
            )

    return BranchGraph(
        logical_graph=logical_graph,
        branches=branches,
    )