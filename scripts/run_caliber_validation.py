"""Run deterministic synthetic validation; never read or write real case data."""
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))  # Support python3 scripts/run_caliber_validation.py.
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "validation" / "outputs" / ".matplotlib"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np

from validation.synthetic_caliber.evaluate import evaluate
from validation.synthetic_caliber.geometries import (
    straight_cylinder, oblique_cylinder, curved_vessel, tapered_vessel,
)

SPACINGS = ((0.5, 0.5, 0.5), (0.46875, 0.46875, 0.8))


def experiment_matrix():
    for spacing in SPACINGS:
        for factory, diameters in (
            (straight_cylinder, (1, 2, 3, 4, 5)),
            (oblique_cylinder, (2, 3, 4)),
            (curved_vessel, (2, 3, 4)),
        ):
            for diameter in diameters:
                yield factory(diameter_mm=diameter, spacing=spacing)
        for proximal, distal in ((4, 2), (3, 1.5)):
            yield tapered_vessel(proximal, distal, spacing=spacing)


def write_csv(path, rows):
    columns = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def aggregate(rows):
    if not rows:
        return {"runs": 0, "mae_mm": None, "rmse_mm": None, "signed_bias_mm": None}
    return {
        "runs": len(rows), "mae_mm": float(np.mean([r["mae_mm"] for r in rows])),
        "rmse_mm": float(np.sqrt(np.mean([r["rmse_mm"] ** 2 for r in rows]))),
        "signed_bias_mm": float(np.mean([r["signed_bias_mm"] for r in rows])),
    }


def figures(output, rows, samples):
    directory = output / "figures"
    directory.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "figure.figsize": (6.4, 4.4), "savefig.dpi": 180})

    def save(name, xlabel, ylabel, title):
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(alpha=0.2)
        plt.legend(fontsize=8)
        plt.tight_layout()
        plt.savefig(directory / name)
        plt.close()

    for metric, name, ylabel in (
        ("estimated_mean_mm", "estimated_vs_true.png", "Estimated mean diameter (mm)"),
        ("mae_mm", "absolute_error_vs_true.png", "Mean absolute sample error (mm)"),
        ("percent_error", "percent_error_vs_true.png", "Mean absolute percentage error (%)"),
    ):
        plt.figure()
        for geometry in ("straight", "oblique", "curved"):
            for sz, marker, style, label in ((0.5, "o", "-", "isotropic"), (0.8, "s", "--", "anisotropic")):
                selected = sorted([r for r in rows if r["geometry"] == geometry and r["spacing_z_mm"] == sz and r["trim_mm"] == 2], key=lambda r: r["true_diameter_mm"])
                plt.plot([r["true_diameter_mm"] for r in selected], [r[metric] for r in selected], marker=marker, linestyle=style, label=f"{geometry}, {label}")
        if metric == "estimated_mean_mm":
            plt.plot([1, 5], [1, 5], ":", label="Identity y=x")
        save(name, "True diameter (mm)", ylabel, "Constant-caliber vessels · 2 mm endpoint trimming")

    for row in rows:
        if row["trim_mm"] != 2 or not (row["geometry"] == "tapered" or (row["geometry"] == "curved" and row["true_diameter_mm"] == 3)):
            continue
        selected = [s for s in samples if s["run_id"] == row["run_id"]]
        axis = "normalized_position" if row["geometry"] == "tapered" else "distance_mm"
        selected.sort(key=lambda s: s[axis])
        plt.figure()
        for key, label, style in (("true_diameter_mm", "Analytic truth", "--"), ("raw_diameter_mm", "Production raw EDT", "-"), ("smooth_diameter_mm", "Production smoothed", ":")):
            plt.plot([s[axis] for s in selected], [s[key] for s in selected], style, label=label)
        save(f'{row["run_id"]}.png', "Analytic normalized position" if axis == "normalized_position" else "Distance along recovered path (mm)", "Diameter (mm)", f'{row["geometry"].title()} · {row["proximal_diameter_mm"]} → {row["distal_diameter_mm"]} mm · {row["spacing_z_mm"]} mm z-spacing')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "outputs" / "caliber_v1")
    args = parser.parse_args()
    output = args.output.resolve()
    if not output.is_relative_to(ROOT / "validation" / "outputs"):
        parser.error("Output must stay under validation/outputs to protect case outputs")
    output.mkdir(parents=True, exist_ok=True)
    (output / "masks").mkdir(exist_ok=True)
    rows, samples, geometries = [], [], []
    for index, vessel in enumerate(experiment_matrix(), 1):
        geometry_id = f'{index:02d}_{vessel.geometry}_{vessel.proximal_diameter_mm:g}_{vessel.distal_diameter_mm:g}_z{vessel.volume.spacing[2]:g}'
        geometries.append({"geometry_id": geometry_id, **vessel.metadata()})
        nib.save(nib.Nifti1Image(vessel.volume.data, vessel.volume.affine), output / "masks" / f"{geometry_id}.nii.gz")
        pair = []
        for trim in (0.0, 2.0):
            run_id = f"{geometry_id}_trim{trim:g}"
            try:
                row, profile = evaluate(vessel, trim)
                samples.extend({"run_id": run_id, **sample} for sample in profile)
            except (ValueError, RuntimeError) as error:
                row = {
                    "geometry": vessel.geometry, "trim_mm": trim,
                    "spacing_x_mm": vessel.volume.spacing[0], "spacing_y_mm": vessel.volume.spacing[1], "spacing_z_mm": vessel.volume.spacing[2],
                    "proximal_diameter_mm": vessel.proximal_diameter_mm, "distal_diameter_mm": vessel.distal_diameter_mm,
                    "status": "failed", "error": str(error), "sample_count": 0,
                }
            row.update(run_id=run_id, geometry_id=geometry_id)
            pair.append(row)
            rows.append(row)
            print(f'{run_id}: {row["status"]}, MAE={row.get("mae_mm")}', flush=True)
        for row in pair:
            row["trim_applied"] = row["trim_mm"] > 0 and row["sample_count"] < pair[0]["sample_count"] if all(r["status"] == "ok" for r in pair) else None
            row["samples_removed"] = pair[0]["sample_count"] - row["sample_count"] if row["status"] == "ok" else None

    successful = [r for r in rows if r["status"] == "ok"]
    for row in successful:
        baseline = next((r for r in successful if r["geometry"] == "straight" and r["true_diameter_mm"] == row["true_diameter_mm"] and r["spacing_z_mm"] == row["spacing_z_mm"] and r["trim_mm"] == row["trim_mm"]), None)
        row["mae_delta_vs_straight_mm"] = row["mae_mm"] - baseline["mae_mm"] if baseline else None
    all_stats = aggregate(successful)
    worst = max(successful, key=lambda r: r["percent_error"]) if successful else None
    summary = {
        "version": "caliber_v1", "number_of_geometries": len(geometries), "number_of_experiments": len(rows),
        "successful_runs": len(successful), "failed_runs": len(rows) - len(successful),
        "overall_mae_mm": all_stats["mae_mm"], "overall_rmse_mm": all_stats["rmse_mm"],
        "overall_signed_bias_mm": all_stats["signed_bias_mm"],
        "by_trim": {str(trim): aggregate([r for r in successful if r["trim_mm"] == trim]) for trim in (0.0, 2.0)},
        "isotropic_mae_mm": aggregate([r for r in successful if r["spacing_z_mm"] == 0.5])["mae_mm"],
        "anisotropic_mae_mm": aggregate([r for r in successful if r["spacing_z_mm"] == 0.8])["mae_mm"],
        "worst_case_geometry": worst["geometry"] if worst else None,
        "worst_case_percent_error": worst["percent_error"] if worst else None,
        "worst_case_run": worst,
        "by_geometry_spacing_trim": [
            {"geometry": g, "spacing_z_mm": sz, "trim_mm": trim, **aggregate([r for r in successful if r["geometry"] == g and r["spacing_z_mm"] == sz and r["trim_mm"] == trim])}
            for g in ("straight", "oblique", "curved", "tapered") for sz in (0.5, 0.8) for trim in (0.0, 2.0)
        ],
        "runtime": {"python": platform.python_version(), "packages": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "scikit-image", "networkx", "nibabel", "matplotlib")}},
        "source_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for pattern in ("neurovasc/geometry/*.py", "neurovasc/graph/*.py", "neurovasc/io/*.py", "validation/synthetic_caliber/*.py", "scripts/run_caliber_validation.py")
            for path in sorted(ROOT.glob(pattern))
        },
        "notes": [
            "Geometric/synthetic validation of a segmentation-derived estimator; not clinical validation.",
            "Overall MAE is the equal-run mean of sample MAEs; overall RMSE is sqrt(mean of run MSEs). Both trim settings count equally; runs are paired, not independent replicates.",
            "percent_error is mean absolute sample-relative error times 100; absolute_error_mm is absolute mean signed error (not MAE).",
            "Raw production EDT is primary; smoothing is retained for secondary comparison and figures.",
            "Failures are recorded and excluded from aggregates, never treated as zero error.",
            "One grid phase and orientation per geometry; no noise, partial volume, or segmentation uncertainty.",
            "Primary comparison figures use trim=2 mm; all untrimmed and trimmed profiles are saved in profiles.csv.",
        ],
    }
    write_csv(output / "results.csv", rows)
    write_csv(output / "profiles.csv", samples)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    (output / "geometries.json").write_text(json.dumps(geometries, indent=2, allow_nan=False) + "\n")
    figures(output, successful, samples)
    (output / "manifest.txt").write_text("\n".join(sorted(
        {str(p.relative_to(output)) for p in output.rglob("*") if p.is_file()} | {"manifest.txt"}
    )) + "\n")
    print(json.dumps({k: summary[k] for k in ("number_of_experiments", "failed_runs", "overall_mae_mm", "overall_rmse_mm", "isotropic_mae_mm", "anisotropic_mae_mm", "worst_case_percent_error")}, indent=2))
    if summary["failed_runs"]:
        raise SystemExit("Some experiments failed; inspect results.csv")


if __name__ == "__main__":
    main()
