from pathlib import Path

import nibabel as nib
import numpy as np


def main() -> None:
    shape = (128, 128, 96)

    data = np.zeros(
        shape,
        dtype=np.uint8,
    )

    # Create a simple spherical structure.
    x, y, z = np.indices(shape)

    center = np.array(
        [
            shape[0] / 2,
            shape[1] / 2,
            shape[2] / 2,
        ]
    )

    radius = 20

    distance = np.sqrt(
        (x - center[0]) ** 2
        + (y - center[1]) ** 2
        + (z - center[2]) ** 2
    )

    data[distance <= radius] = 1

    # 0.5 x 0.5 x 0.8 mm voxel spacing
    affine = np.diag(
        [0.5, 0.5, 0.8, 1.0]
    )

    image = nib.Nifti1Image(
        data,
        affine,
    )

    output = Path(
        "data/samples/demo_volume.nii.gz"
    )

    output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    nib.save(
        image,
        output,
    )

    print(
        f"Created demo volume: {output}"
    )


if __name__ == "__main__":
    main()