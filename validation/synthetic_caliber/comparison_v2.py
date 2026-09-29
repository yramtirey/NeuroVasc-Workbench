"""Matched adapters and explicit missingness accounting for caliber v2."""
from dataclasses import dataclass
from typing import Callable

import numpy as np

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.main_path import extract_main_path, MainCenterlinePath
from neurovasc.geometry.diameter_profile import build_diameter_profile, DiameterProfile
from neurovasc.geometry.cross_sectional_caliber import measure_cross_sectional_caliber
from neurovasc.graph.vascular_graph import build_vascular_graph
from .geometries import SyntheticVessel
from .metrics import comparison_metrics


@dataclass
class Context:
    vessel: SyntheticVessel
    path: MainCenterlinePath
    edt_profile: DiameterProfile


@dataclass
class Estimates:
    diameter_mm: np.ndarray
    area_mm2: np.ndarray
    tangents: np.ndarray
    reasons: list[str]


def prepare(vessel):
    path = extract_main_path(build_vascular_graph(extract_centerline(vessel.volume)))
    return Context(vessel, path, build_diameter_profile(vessel.volume, path, trim_mm=0))


def edt_adapter(context):
    d = context.edt_profile.diameter_mm
    return Estimates(d, np.full(len(d), np.nan), np.full((len(d), 3), np.nan), ["ok"] * len(d))


def cross_sectional_adapter(context):
    result = measure_cross_sectional_caliber(context.vessel.volume, context.path, tangent_window_mm=6.0)
    return Estimates(result.diameter_mm, result.area_mm2, result.tangent_vectors, result.failure_reasons)


ESTIMATORS: dict[str, Callable[[Context], Estimates]] = {
    "edt": edt_adapter, "cross_sectional": cross_sectional_adapter,
}


def trimmed_indices(context, trim_mm):
    # Ask production code which samples survive; do not copy its trimming rules.
    profile = build_diameter_profile(context.vessel.volume, context.path, trim_mm=trim_mm)
    lookup = {tuple(point): i for i, point in enumerate(context.path.world_points)}
    return np.array([lookup[tuple(point)] for point in profile.world_points], dtype=int), profile.distance_mm


def score(estimated, truth):
    estimated, truth = np.asarray(estimated, dtype=float), np.asarray(truth, dtype=float)
    if estimated.ndim != 1 or estimated.shape != truth.shape:
        raise ValueError("Expected matching one-dimensional profiles")
    valid = np.isfinite(estimated) & (estimated > 0)
    metrics = comparison_metrics(estimated[valid], truth[valid]) if valid.any() else {
        key: None for key in comparison_metrics(np.array([1.0]), np.array([1.0]))
    }
    metrics.update(
        sample_count=len(estimated), valid_sample_count=int(valid.sum()),
        invalid_sample_count=int((~valid).sum()),
        valid_fraction=float(valid.mean()) if len(valid) else None,
        status="ok" if valid.all() and len(valid) else "partial" if valid.any() else "invalid",
    )
    return metrics


def matched_metrics(edt, cross, truth, tie_tolerance_mm=0.01):
    """Compare on identical valid support; missing samples remain in coverage."""
    edt, cross, truth = map(np.asarray, (edt, cross, truth))
    if not np.isfinite(tie_tolerance_mm) or tie_tolerance_mm < 0:
        raise ValueError("Tie tolerance must be finite and nonnegative")
    if edt.shape != cross.shape or edt.shape != truth.shape:
        raise ValueError("Matched profiles must have identical shapes")
    valid = np.isfinite(edt) & np.isfinite(cross) & (edt > 0) & (cross > 0)
    result = {
        "common_valid_samples": int(valid.sum()), "total_samples": len(truth),
        "common_valid_fraction": float(valid.mean()) if len(truth) else None,
        "tie_tolerance_mm": tie_tolerance_mm,
        "comparison_status": "complete" if valid.all() and len(valid) else "partial" if valid.any() else "unscorable",
        "winner": "unscorable",
    }
    if valid.any():
        a, b = comparison_metrics(edt[valid], truth[valid]), comparison_metrics(cross[valid], truth[valid])
        for short, key in (("mae", "mae_mm"), ("rmse", "rmse_mm"), ("bias", "signed_bias_mm")):
            result[f"edt_{short}"] = a[key]
            result[f"cross_sectional_{short}"] = b[key]
            result[f"cross_sectional_{short}_minus_edt_{short}"] = b[key] - a[key]
        delta = b["mae_mm"] - a["mae_mm"]
        result["winner"] = "tie" if abs(delta) <= tie_tolerance_mm else "cross_sectional" if delta < 0 else "edt"
    return result
