from dataclasses import dataclass

import networkx as nx
import numpy as np

from neurovasc.graph.vascular_graph import VascularGraph


@dataclass
class MainCenterlinePath:
    voxel_points: np.ndarray
    world_points: np.ndarray
    distance_mm: np.ndarray

    @property
    def length_mm(self) -> float:
        if len(self.distance_mm) == 0:
            return 0.0

        return float(self.distance_mm[-1])


def _component_length(
    graph: nx.Graph,
) -> float:
    return float(
        sum(
            data["length_mm"]
            for _, _, data in graph.edges(data=True)
        )
    )


def extract_main_path(
    vascular_graph: VascularGraph,
) -> MainCenterlinePath:
    """
    Extract the primary ordered path through a vascular graph.

    The largest connected component is selected, then the
    longest weighted path between terminal endpoints is used
    as the vessel's main centerline.
    """

    graph = vascular_graph.graph

    if graph.number_of_nodes() == 0:
        raise ValueError(
            "Cannot extract main path from an empty graph."
        )

    #
    # Find connected components.
    #
    components = [
        graph.subgraph(component).copy()
        for component in nx.connected_components(graph)
    ]

    #
    # Prefer the component with the greatest physical length.
    #
    main_component = max(
        components,
        key=_component_length,
    )

    endpoints = [
        node
        for node, degree
        in main_component.degree()
        if degree == 1
    ]

    #
    # Normal vessel case:
    # find longest shortest path between endpoints.
    #
    if len(endpoints) >= 2:

        best_start = None
        best_end = None
        best_length = -1.0

        for i, start in enumerate(endpoints):

            lengths = nx.single_source_dijkstra_path_length(
                main_component,
                start,
                weight="length_mm",
            )

            for end in endpoints[i + 1:]:

                if end not in lengths:
                    continue

                length = lengths[end]

                if length > best_length:
                    best_length = length
                    best_start = start
                    best_end = end

        path_nodes = nx.shortest_path(
            main_component,
            best_start,
            best_end,
            weight="length_mm",
        )

    else:
        #
        # Fallback for loop-like structures with no clear
        # terminal endpoints.
        #
        start = next(
            iter(main_component.nodes)
        )

        distances = nx.single_source_dijkstra_path_length(
            main_component,
            start,
            weight="length_mm",
        )

        farthest_a = max(
            distances,
            key=distances.get,
        )

        distances = nx.single_source_dijkstra_path_length(
            main_component,
            farthest_a,
            weight="length_mm",
        )

        farthest_b = max(
            distances,
            key=distances.get,
        )

        path_nodes = nx.shortest_path(
            main_component,
            farthest_a,
            farthest_b,
            weight="length_mm",
        )

    voxel_points = np.asarray(
        path_nodes,
        dtype=int,
    )

    world_points = np.asarray(
        [
            main_component.nodes[node]["world"]
            for node in path_nodes
        ],
        dtype=float,
    )

    segment_lengths = np.linalg.norm(
        np.diff(
            world_points,
            axis=0,
        ),
        axis=1,
    )

    distance_mm = np.concatenate(
        (
            [0.0],
            np.cumsum(segment_lengths),
        )
    )

    return MainCenterlinePath(
        voxel_points=voxel_points,
        world_points=world_points,
        distance_mm=distance_mm,
    )