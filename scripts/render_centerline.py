import argparse

import numpy as np
import pyvista as pv

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import load_nifti


def graph_to_polydata(vascular_graph) -> pv.PolyData:
    """Convert the NetworkX centerline graph into PyVista line geometry."""

    graph = vascular_graph.graph

    nodes = list(graph.nodes())

    node_to_index = {
        node: index
        for index, node in enumerate(nodes)
    }

    points = np.array(
        [
            graph.nodes[node]["world"]
            for node in nodes
        ],
        dtype=float,
    )

    lines = []

    for node_a, node_b in graph.edges():

        lines.extend(
            [
                2,
                node_to_index[node_a],
                node_to_index[node_b],
            ]
        )

    polydata = pv.PolyData(points)

    polydata.lines = np.asarray(
        lines,
        dtype=np.int64,
    )

    return polydata


def node_cloud(
    vascular_graph,
    nodes,
) -> pv.PolyData:
    """Create a PyVista point cloud from selected graph nodes."""

    graph = vascular_graph.graph

    if len(nodes) == 0:
        return pv.PolyData()

    points = np.array(
        [
            graph.nodes[node]["world"]
            for node in nodes
        ],
        dtype=float,
    )

    return pv.PolyData(points)


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Render vessel surface and centerline "
            "with NeuroVasc Workbench."
        )
    )

    parser.add_argument(
        "path",
        help="Path to vascular NIfTI mask.",
    )

    args = parser.parse_args()

    print()
    print("NeuroVasc Workbench")
    print("=" * 55)

    volume = load_nifti(args.path)

    print("Generating vessel mesh...")

    vessel = mask_to_mesh(volume)

    print("Extracting centerline...")

    centerline = extract_centerline(volume)

    print("Building vascular graph...")

    vascular_graph = build_vascular_graph(
        centerline
    )

    print()
    print("Centerline Geometry")
    print("-" * 55)

    print(
        f"Centerline points: "
        f"{centerline.n_points:,}"
    )

    print(
        f"Graph nodes:       "
        f"{vascular_graph.n_nodes:,}"
    )

    print(
        f"Graph edges:       "
        f"{vascular_graph.n_edges:,}"
    )

    print(
        f"Endpoints:         "
        f"{len(vascular_graph.endpoints):,}"
    )

    print(
        f"Junction voxels:   "
        f"{len(vascular_graph.branchpoints):,}"
    )

    print(
        f"Total length:      "
        f"{vascular_graph.total_length_mm:.2f} mm"
    )

    print("-" * 55)

    centerline_mesh = graph_to_polydata(
        vascular_graph
    )

    endpoint_points = node_cloud(
        vascular_graph,
        vascular_graph.endpoints,
    )

    junction_points = node_cloud(
        vascular_graph,
        vascular_graph.branchpoints,
    )

    plotter = pv.Plotter()

    # Vessel surface
    plotter.add_mesh(
        vessel.mesh,
        opacity=0.25,
        smooth_shading=True,
    )

    # Centerline
    plotter.add_mesh(
        centerline_mesh,
        line_width=5,
    )

    # Endpoints
    if endpoint_points.n_points > 0:
        plotter.add_mesh(
            endpoint_points,
            point_size=18,
            render_points_as_spheres=True,
        )

    # Junctions
    if junction_points.n_points > 0:
        plotter.add_mesh(
            junction_points,
            point_size=14,
            render_points_as_spheres=True,
        )

    plotter.add_axes()
    plotter.show_grid()

    plotter.camera_position = "iso"

    plotter.show()


if __name__ == "__main__":
    main()