"""Adapters call the production estimator unchanged; alternatives can be injected."""
from collections.abc import Callable

from neurovasc.geometry.centerline import extract_centerline
from neurovasc.geometry.diameter_profile import DiameterProfile, build_diameter_profile
from neurovasc.geometry.main_path import extract_main_path
from neurovasc.geometry.profile_analysis import analyze_profile
from neurovasc.graph.vascular_graph import build_vascular_graph
from .geometries import SyntheticVessel
from .metrics import comparison_metrics

Estimator = Callable[[SyntheticVessel, float], DiameterProfile]


def production_estimator(vessel: SyntheticVessel, trim_mm: float) -> DiameterProfile:
    centerline = extract_centerline(vessel.volume)
    path = extract_main_path(build_vascular_graph(centerline))
    return build_diameter_profile(vessel.volume, path, trim_mm=trim_mm)


def evaluate(vessel: SyntheticVessel, trim_mm: float = 2.0,
             estimator: Estimator = production_estimator, estimator_name="production_edt"):
    profile = estimator(vessel, trim_mm)
    analysis = analyze_profile(profile)
    truth = vessel.true_diameter(profile.world_points)
    u = vessel.normalized_position(profile.world_points)
    row = {
        "geometry": vessel.geometry, "estimator": estimator_name,
        "spacing_x_mm": vessel.volume.spacing[0],
        "spacing_y_mm": vessel.volume.spacing[1],
        "spacing_z_mm": vessel.volume.spacing[2],
        "true_diameter_mm": vessel.proximal_diameter_mm if vessel.geometry != "tapered" else None,
        "proximal_diameter_mm": vessel.proximal_diameter_mm,
        "distal_diameter_mm": vessel.distal_diameter_mm,
        "trim_mm": trim_mm, "status": "ok", "error": None,
        **comparison_metrics(profile.diameter_mm, truth),
        "smooth_mae_mm": comparison_metrics(analysis.smooth_diameter_mm, truth)["mae_mm"],
        "normalized_position_min": float(u.min()), "normalized_position_max": float(u.max()),
        "profile_length_mm": profile.length_mm,
    }
    samples = [
        {
            "distance_mm": float(distance), "normalized_position": float(position),
            "true_diameter_mm": float(true), "raw_diameter_mm": float(raw),
            "smooth_diameter_mm": float(smooth), "error_mm": float(raw - true),
            "x_mm": float(point[0]), "y_mm": float(point[1]), "z_mm": float(point[2]),
        }
        for distance, position, true, raw, smooth, point in zip(
            profile.distance_mm, u, truth, profile.diameter_mm,
            analysis.smooth_diameter_mm, profile.world_points, strict=True,
        )
    ]
    return row, samples
