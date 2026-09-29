import argparse

import numpy as np

from neurovasc.io.volume import load_nifti


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect a 3D imaging volume "
            "with NeuroVasc Workbench."
        )
    )

    parser.add_argument(
        "path",
        help="Path to a .nii or .nii.gz file",
    )

    args = parser.parse_args()

    volume = load_nifti(args.path)

    data = volume.data

    print()
    print("NeuroVasc Workbench")
    print("=" * 55)

    print(f"File:            {volume.path.name}")
    print(f"Shape:           {volume.shape}")
    print(f"Voxel spacing:   {volume.spacing} mm")

    size = volume.physical_size_mm

    print(
        "Physical size:   "
        f"{size[0]:.2f} × "
        f"{size[1]:.2f} × "
        f"{size[2]:.2f} mm"
    )

    print(
        f"Voxel volume:    "
        f"{volume.voxel_volume_mm3:.6f} mm³"
    )

    print("-" * 55)

    print(f"Data type:       {data.dtype}")
    print(f"Minimum:         {np.min(data):.3f}")
    print(f"Maximum:         {np.max(data):.3f}")
    print(f"Mean:            {np.mean(data):.3f}")
    print(f"Median:          {np.median(data):.3f}")
    print(f"Nonzero voxels:  {np.count_nonzero(data):,}")

    print("-" * 55)
    print("Affine transformation:")
    print(volume.affine)

    print("=" * 55)


if __name__ == "__main__":
    main()