from dataclasses import dataclass
from itertools import product

import networkx as nx
import numpy as np

from neurovasc.geometry.centerline import Centerline


NEIGHBOR_OFFSETS = [
    offset
    for offset in product(
        (-1, 0, 1),
        repeat=3,
    )
    if offset != (0, 0, 0)
]


@dataclass
class VascularGraph:
    """Graph representation of a vascular centerline."""

    graph: nx.Graph

    @property
    def n_nodes(self) -> int:
        return self.graph.number_of_nodes()

    @property
    def n_edges(self) -> int:
        return self.graph.number_of_edges()

    @property
    def endpoints(self) -> list:
        """Centerline voxels with only one neighbor."""
        return [
            node
            for node, degree in self.graph.degree()
            if degree == 1
        ]

    @property
    def branchpoints(self) -> list:
        """Centerline voxels connected to three or more neighbors."""
        return [
            node
            for node, degree in self.graph.degree()
            if degree >= 3
        ]

    @property
    def n_components(self) -> int:
        return nx.number_connected_components(
            self.graph
        )

    @property
    def total_length_mm(self) -> float:
        return float(
            sum(
                data["length_mm"]
                for _, _, data
                in self.graph.edges(data=True)
            )
        )


def build_vascular_graph(
    centerline: Centerline,
) -> VascularGraph:
    """
    Convert a 3D centerline into a NetworkX graph.

    Each skeleton voxel becomes a graph node.
    Neighboring skeleton voxels become graph edges.
    """

    graph = nx.Graph()

    voxel_coordinates = [
        tuple(int(value) for value in point)
        for point in centerline.voxel_points
    ]

    coordinate_set = set(
        voxel_coordinates
    )

    world_lookup = {
        coordinate: np.asarray(
            world,
            dtype=float,
        )
        for coordinate, world in zip(
            voxel_coordinates,
            centerline.world_points,
        )
    }

    #
    # Add nodes
    #
    for coordinate in voxel_coordinates:

        world = world_lookup[coordinate]

        graph.add_node(
            coordinate,
            voxel=coordinate,
            world=tuple(
                float(value)
                for value in world
            ),
        )

    #
    # Add 26-connected edges
    #
    for coordinate in voxel_coordinates:

        world_a = world_lookup[
            coordinate
        ]

        for offset in NEIGHBOR_OFFSETS:

            neighbor = tuple(
                coordinate[i] + offset[i]
                for i in range(3)
            )

            if neighbor not in coordinate_set:
                continue

            # Avoid adding every edge twice.
            if coordinate >= neighbor:
                continue

            world_b = world_lookup[
                neighbor
            ]

            distance_mm = float(
                np.linalg.norm(
                    world_b - world_a
                )
            )

            graph.add_edge(
                coordinate,
                neighbor,
                length_mm=distance_mm,
            )

    return VascularGraph(
        graph=graph
    )