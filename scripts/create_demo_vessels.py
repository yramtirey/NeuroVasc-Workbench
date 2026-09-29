from pathlib import Path

import nibabel as nib
import numpy as np


def add_vessel(
    mask: np.ndarray,
    start: tuple[float, float, float],
    end: tuple[float, float, float],
    radius: float,
) -> None:
    """
    Draw a cylindrical vessel between two 3D points.
    """

    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)

    x, y, z = np.indices(mask.shape)

    points = np.stack(
        (x, y, z),
        axis=-1,
    )

    direction = end - start

    length_squared = np.dot(
        direction,
        direction,
    )

    if length_squared == 0:
        raise ValueError(
            "Vessel start and end points cannot be identical."
        )

    relative = points - start

    t = np.sum(
        relative * direction,
        axis=-1,
    ) / length_squared

    t = np.clip(
        t,
        0.0,
        1.0,
    )

    projection = (
        start
        + t[..., None] * direction
    )

    distance = np.linalg.norm(
        points - projection,
        axis=-1,
    )

    mask[distance <= radius] = 1


def main() -> None:
    shape = (128, 128, 128)

    mask = np.zeros(
        shape,
        dtype=np.uint8,
    )

    # Main vessel trunk
    add_vessel(
        mask,
        start=(64, 64, 12),
        end=(64, 64, 70),
        radius=5,
    )

    # Left branch
    add_vessel(
        mask,
        start=(64, 64, 64),
        end=(35, 45, 100),
        radius=4,
    )

    # Right branch
    add_vessel(
        mask,
        start=(64, 64, 64),
        end=(93, 45, 100),
        radius=4,
    )

    # Secondary left branch
    add_vessel(
        mask,
        start=(44, 51, 88),
        end=(25, 70, 116),
        radius=3,
    )

    # Secondary right branch
    add_vessel(
        mask,
        start=(84, 51, 88),
        end=(103, 70, 116),
        radius=3,
    )

    # Small anterior-style branch
    add_vessel(
        mask,
        start=(64, 64, 67),
        end=(64, 92, 96),
        radius=3,
    )

    # 0.5 mm isotropic voxels
    affine = np.diag(
        [0.5, 0.5, 0.5, 1.0]
    )

    image = nib.Nifti1Image(
        mask,
        affine,
    )

    output = Path(
        "data/samples/demo_vessels.nii.gz"
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    nib.save(
        image,
        output,
    )

    print()
    print("NeuroVasc vascular demo created ✓")
    print(f"Saved: {output}")
    print(f"Shape: {mask.shape}")
    print(f"Vessel voxels: {mask.sum():,}")


if __name__ == "__main__":
    main()