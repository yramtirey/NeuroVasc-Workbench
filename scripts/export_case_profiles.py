"""Export UI profiles using the same analysis functions as the case report."""
import argparse
import json
import math
from pathlib import Path

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.diameter_profile import build_diameter_profile
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.profile_analysis import analyze_profile
from neurovasc.graph.vascular_graph import build_vascular_graph
from neurovasc.io.volume import load_nifti
from neurovasc.report.case_report import make_label_volume


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", help="Original case segmentation NIfTI")
    parser.add_argument("--case-id", default="case_001")
    parser.add_argument("--trim-mm", type=float, default=2.0)
    args = parser.parse_args()
    case_dir = Path(__file__).resolve().parents[1] / "outputs" / "cases" / args.case_id
    report = json.loads((case_dir / "case_report.json").read_text())
    volume = load_nifti(args.path)
    if volume.path.name != report["source_file"]:
        raise ValueError("Segmentation must match the case report source file")
    profiles_dir = case_dir / "profiles"
    profiles_dir.mkdir(exist_ok=True)

    for vessel in report["vessels"]:
        label = vessel["label"]
        output = profiles_dir / f"label_{label}.json"
        if output.exists():
            print(f"Preserving existing {output.name}", flush=True)
            continue
        label_volume = make_label_volume(volume, label)
        graph = build_vascular_graph(extract_centerline(label_volume))
        profile = build_diameter_profile(label_volume, extract_main_path(graph), trim_mm=args.trim_mm)
        analysis = analyze_profile(profile)
        metrics = {
            "reference_diameter_mm": analysis.reference_diameter_mm,
            "minimum_diameter_mm": analysis.minimum_diameter_mm,
            "minimum_location_mm": analysis.minimum_distance_mm,
            "candidate_reduction_percent": analysis.maximum_reduction_percent,
        }
        for key, value in metrics.items():
            if not math.isclose(value, vessel[key], rel_tol=1e-6, abs_tol=1e-6):
                raise ValueError(f"Label {label}: {key} differs from the case report; check source and trim settings")
        samples = [
            dict(zip(
                ("distance_mm", "raw_diameter_mm", "smooth_diameter_mm", "x_mm", "y_mm", "z_mm"),
                map(float, (distance, raw, smooth, *point)),
            ))
            for distance, raw, smooth, point in zip(
                analysis.distance_mm, analysis.raw_diameter_mm,
                analysis.smooth_diameter_mm, profile.world_points, strict=True,
            )
        ]
        payload = {
            "case_id": report["case_id"], "label": label,
            "abbreviation": vessel["abbreviation"], "name": vessel["name"],
            **metrics, "samples": samples, "note": report["measurement_note"],
        }
        # Exclusive creation protects existing analysis outputs.
        with output.open("x") as file:
            file.write(json.dumps(payload, indent=2, allow_nan=False) + "\n")
        print(f"Exported {output.name}: {len(samples)} samples", flush=True)


if __name__ == "__main__":
    main()
