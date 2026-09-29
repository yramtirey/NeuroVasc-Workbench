from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d

from neurovasc.geometry.diameter_profile import DiameterProfile


@dataclass
class ProfileAnalysis:
    distance_mm: np.ndarray
    raw_diameter_mm: np.ndarray
    smooth_diameter_mm: np.ndarray
    reference_diameter_mm: float
    relative_diameter: np.ndarray
    reduction_percent: np.ndarray
    minimum_index: int

    @property
    def minimum_distance_mm(self) -> float:
        return float(
            self.distance_mm[self.minimum_index]
        )

    @property
    def minimum_diameter_mm(self) -> float:
        return float(
            self.smooth_diameter_mm[self.minimum_index]
        )

    @property
    def maximum_reduction_percent(self) -> float:
        return float(
            self.reduction_percent[self.minimum_index]
        )


def analyze_profile(
    profile: DiameterProfile,
    sigma: float = 1.25,
) -> ProfileAnalysis:
    """
    Smooth a vessel caliber profile and quantify reductions
    relative to a robust reference caliber.

    This is a geometric QC metric, not a diagnostic stenosis score.
    """

    raw = np.asarray(
        profile.diameter_mm,
        dtype=float,
    )

    distance = np.asarray(
        profile.distance_mm,
        dtype=float,
    )

    if len(raw) < 3:
        raise ValueError(
            "Diameter profile contains less than 3 samples."
        )

    smooth = gaussian_filter1d(
        raw,
        sigma=sigma,
        mode="nearest",
    )

    # Robust reference that is less influenced by isolated
    # unusually large or small samples.
    reference = float(
        np.percentile(
            smooth,
            75,
        )
    )

    relative = (
        smooth / reference
    )

    reduction = (
        1.0 - relative
    ) * 100.0

    minimum_index = int(
        np.argmin(smooth)
    )

    return ProfileAnalysis(
        distance_mm=distance,
        raw_diameter_mm=raw,
        smooth_diameter_mm=smooth,
        reference_diameter_mm=reference,
        relative_diameter=relative,
        reduction_percent=reduction,
        minimum_index=minimum_index,
    )