import argparse

import numpy as np
import pyvista as pv

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.graph.branch_graph import build_branch_graph
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import load_nifti


def make_branch_polyline(
    world_path: np.ndarray,
) -> pv.PolyData:
    """
    Convert an ordered branch centerline path into a PyVista polyline.
    """

    if len(world_path) < 2:
        return pv.PolyData()

    return pv.lines_from_points(
        world_path,
        close=False,
    )


def midpoint(
    points: np.ndarray,
) -> np.ndarray:
    """
    Return an approximate center point along a branch path.
    """

    return points[len(points) // 2]


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Render branch-level vascular topology "
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
    print("=" * 60)

    #
    # Load volume
    #
    volume = load_nifti(
        args.path
    )

    print("Generating vessel surface...")

    vessel_mesh = mask_to_mesh(
        volume
    )

    print("Extracting centerline...")

    centerline = extract_centerline(
        volume
    )

    print("Building voxel graph...")

    vascular_graph = build_vascular_graph(
        centerline
    )

    print("Extracting logical branches...")

    branch_graph = build_branch_graph(
        vascular_graph
    )

    print()
    print("Branch Visualization")
    print("-" * 60)

    print(
        f"Branches:          "
        f"{branch_graph.n_branches}"
    )

    print(
        f"Logical junctions: "
        f"{len(branch_graph.junction_nodes)}"
    )

    print(
        f"Logical endpoints: "
        f"{len(branch_graph.endpoint_nodes)}"
    )

    print("-" * 60)

    plotter = pv.Plotter()

    #
    # Transparent vessel shell
    #
    plotter.add_mesh(
        vessel_mesh.mesh,
        opacity=0.15,
        smooth_shading=True,
    )

    #
    # Give subsequently added branches different colors.
    #
    plotter.set_color_cycler(
        "default"
    )

    label_points = []
    label_names = []

    #
    # Draw each logical branch individually
    #
    for branch in branch_graph.branches:

        polyline = make_branch_polyline(
            branch.world_path
        )

        if polyline.n_points == 0:
            continue

        #
        # Tube makes centerlines easier to see in 3D.
        #
        tube = polyline.tube(
            radius=0.35
        )

        plotter.add_mesh(
            tube,
            smooth_shading=True,
        )

        label_points.append(
            midpoint(
                branch.world_path
            )
        )

        label_names.append(
            branch.branch_id
        )

    #
    # Add branch IDs
    #
    if label_points:

        plotter.add_point_labels(
            np.asarray(
                label_points
            ),
            label_names,
            font_size=12,
            point_size=5,
            always_visible=True,
        )

    #
    # Logical endpoints
    #
    endpoint_positions = []

    for node_id in branch_graph.endpoint_nodes:

        node_data = (
            branch_graph
            .logical_graph
            .nodes[node_id]
        )

        endpoint_positions.append(
            node_data["world"]
        )

    if endpoint_positions:

        endpoints = pv.PolyData(
            np.asarray(
                endpoint_positions,
                dtype=float,
            )
        )

        plotter.add_mesh(
            endpoints,
            point_size=18,
            render_points_as_spheres=True,
        )

    #
    # Logical junctions
    #
    junction_positions = []

    for node_id in branch_graph.junction_nodes:

        node_data = (
            branch_graph
            .logical_graph
            .nodes[node_id]
        )

        junction_positions.append(
            node_data["world"]
        )

    if junction_positions:

        junctions = pv.PolyData(
            np.asarray(
                junction_positions,
                dtype=float,
            )
        )

        plotter.add_mesh(
            junctions,
            point_size=24,
            render_points_as_spheres=True,
        )

    plotter.add_axes()
    plotter.show_grid()

    plotter.camera_position = "iso"

    print(
        "Opening interactive branch viewer..."
    )

    plotter.show()


if __name__ == "__main__":
    main()