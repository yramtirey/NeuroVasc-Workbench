"""Sample-wise errors; no normalization by recovered path length."""
import numpy as np


def comparison_metrics(estimated, truth):
    estimated, truth = np.asarray(estimated, dtype=float), np.asarray(truth, dtype=float)
    if estimated.ndim != 1 or estimated.shape != truth.shape or not estimated.size:
        raise ValueError("Expected matching nonempty 1D profiles")
    if not np.all(np.isfinite(estimated)) or not np.all(np.isfinite(truth)) or np.any(truth <= 0):
        raise ValueError("Profiles must be finite with positive ground truth")
    errors = estimated - truth
    bias = float(errors.mean())
    return {
        "estimated_mean_mm": float(estimated.mean()),
        "estimated_median_mm": float(np.median(estimated)),
        "estimated_min_mm": float(estimated.min()),
        "estimated_max_mm": float(estimated.max()),
        "mae_mm": float(np.abs(errors).mean()),
        "rmse_mm": float(np.sqrt(np.mean(errors ** 2))),
        "max_absolute_error_mm": float(np.abs(errors).max()),
        "absolute_error_mm": abs(bias),
        "signed_bias_mm": bias,
        "percent_error": float(np.mean(np.abs(errors) / truth) * 100),
        "signed_percent_error": float(np.mean(errors / truth) * 100),
        "correlation": float(np.corrcoef(estimated, truth)[0, 1])
            if np.ptp(truth) > 1e-10 and np.ptp(estimated) > 1e-10 else None,
        "sample_count": int(estimated.size),
    }
