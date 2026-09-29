import argparse
import json
from pathlib import Path

import numpy as np

from neurovasc.io.volume import load_nifti
from neurovasc.report.case_report import (
    analyze_case,
)
from neurovasc.segmentation.labels import (
    TOPCOW_LABELS,
)


def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "Generate a NeuroVasc case-level "
            "vascular analysis."
        )
    )

    parser.add_argument(
        "path",
    )

    parser.add_argument(
        "--case-id",
        default="case_001",
    )

    parser.add_argument(
        "--trim-mm",
        type=float,
        default=2.0,
    )

    args = parser.parse_args()

    volume = load_nifti(
        args.path
    )

    print()
    print("NeuroVasc Workbench")
    print("=" * 72)

    print(
        f"Case:      {args.case_id}"
    )

    print(
        f"Volume:    {volume.path.name}"
    )

    print(
        f"Shape:     {volume.shape}"
    )

    print(
        f"Spacing:   {volume.spacing} mm"
    )

    labels = sorted(
        int(value)
        for value in np.unique(volume.data)
        if int(value) != 0
    )

    print()
    print(
        f"Anatomical labels present: "
        f"{len(labels)}"
    )

    print(
        ", ".join(
            TOPCOW_LABELS.get(
                label,
                str(label),
            )
            for label in labels
        )
    )

    print()
    print(
        "Running vessel-level analysis..."
    )

    table, warnings = analyze_case(
        volume,
        trim_mm=args.trim_mm,
    )

    output_dir = (
        Path("outputs")
        / "cases"
        / args.case_id
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    #
    # CSV
    #
    csv_path = (
        output_dir
        / "vessel_metrics.csv"
    )

    table.to_csv(
        csv_path,
        index=False,
    )

    #
    # JSON
    #
    json_path = (
        output_dir
        / "case_report.json"
    )

    report = {
        "case_id":
            args.case_id,

        "source_file":
            volume.path.name,

        "shape":
            list(volume.shape),

        "spacing_mm":
            list(volume.spacing),

        "labels_present":
            labels,

        "vessel_count":
            len(table),

        "vessels":
            table.to_dict(
                orient="records"
            ),

        "warnings":
            warnings,

        "measurement_note": (
            "Measurements are derived from "
            "segmentation geometry and are not "
            "clinical diagnoses."
        ),
    }

    with open(
        json_path,
        "w",
    ) as file:

        json.dump(
            report,
            file,
            indent=2,
        )

    #
    # Terminal report
    #
    display_columns = [
        "abbreviation",
        "main_path_length_mm",
        "mean_diameter_mm",
        "minimum_diameter_mm",
        "candidate_reduction_percent",
        "centerline_components",
    ]

    print()
    print(
        table[
            display_columns
        ].to_string(
            index=False,
            float_format=lambda value:
                f"{value:.2f}",
        )
    )

    if warnings:

        print()
        print("QC warnings")
        print("-" * 72)

        for warning in warnings:
            print(
                f"• {warning}"
            )

    print()
    print(
        f"CSV  → {csv_path}"
    )

    print(
        f"JSON → {json_path}"
    )

    print("=" * 72)


if __name__ == "__main__":
    main()