import argparse

import pyvista as pv

from neurovasc.geometry.mesh import mask_to_mesh
from neurovasc.io.volume import load_nifti


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Render a vascular segmentation "
            "with NeuroVasc Workbench."
        )
    )

    parser.add_argument(
        "path",
        help="Path to a vascular .nii or .nii.gz mask",
    )

    args = parser.parse_args()

    print()
    print("NeuroVasc Workbench")
    print("=" * 55)

    volume = load_nifti(args.path)

    print(f"Loading: {volume.path.name}")
    print(f"Shape: {volume.shape}")
    print(f"Spacing: {volume.spacing} mm")

    print()
    print("Extracting vascular surface...")

    vessel = mask_to_mesh(volume)

    print()
    print("Vascular Geometry")
    print("-" * 55)

    print(f"Vertices:      {vessel.n_vertices:,}")
    print(f"Faces:         {vessel.n_faces:,}")

    print(
        f"Surface area:  "
        f"{vessel.surface_area_mm2:,.2f} mm²"
    )

    print(
        f"Mesh volume:   "
        f"{vessel.volume_mm3:,.2f} mm³"
    )

    print("-" * 55)

    plotter = pv.Plotter()

    plotter.add_mesh(
        vessel.mesh,
        smooth_shading=True,
        show_edges=False,
    )

    plotter.add_axes()

    plotter.show_grid()

    plotter.camera_position = "iso"

    plotter.show()


if __name__ == "__main__":
    main()