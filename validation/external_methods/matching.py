"""Physical-coordinate support, one-to-one endpoints and explicit missingness."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import cKDTree


def errors(estimate, truth):
    estimate, truth = np.asarray(estimate, float), np.asarray(truth, float)
    if estimate.shape != truth.shape:
        raise ValueError('estimate_truth_shape_mismatch')
    good = np.isfinite(estimate) & np.isfinite(truth) & (truth > 0)
    if not good.any():
        return dict(sample_count=0, estimated_mean_mm=None, mae_mm=None, rmse_mm=None, bias_mm=None, percent_error=None)
    delta = estimate[good] - truth[good]
    return dict(sample_count=int(good.sum()), estimated_mean_mm=float(estimate[good].mean()),
                mae_mm=float(np.abs(delta).mean()), rmse_mm=float(np.sqrt(np.mean(delta**2))),
                bias_mm=float(delta.mean()), percent_error=float(100*np.mean(np.abs(delta)/truth[good])))


def endpoints(found, truth):
    found, truth = np.asarray(found, float), np.asarray(truth, float)
    if not len(found) or not len(truth):
        return dict(endpoint_error_mm=None, matched_endpoint_count=0,
                    missing_endpoint_count=len(truth), extra_endpoint_count=len(found))
    costs = np.linalg.norm(found[:,None]-truth[None], axis=2)
    i,j = linear_sum_assignment(costs)
    return dict(endpoint_error_mm=float(costs[i,j].mean()), matched_endpoint_count=len(i),
                missing_endpoint_count=len(truth)-len(i), extra_endpoint_count=len(found)-len(i))


def sample_on_truth(vessel, points, values, stations, tolerance_mm):
    """Interpolate projected physical stations, never extrapolate or bridge NaNs.

    Only single tubes: ordered projected positions must be monotone (up to
    duplicate projections). A proximity gate and 2-voxel maximum bracket gap
    prevent a short/absent section from inheriting distant values.
    """
    from validation.synthetic_centerline.truth import project
    points, values, stations = np.asarray(points,float), np.asarray(values,float), np.asarray(stations,float)
    if len(points) != len(values) or len(points) < 2:
        return np.full(len(stations), np.nan)
    _,_,u,_ = project(vessel, points)
    _,_,target,_ = project(vessel, stations)
    if np.all(np.diff(u) <= 1e-7):
        u,values,points = u[::-1],values[::-1],points[::-1]
    if np.any(np.diff(u) < -1e-7):
        return np.full(len(stations), np.nan)
    u,idx = np.unique(u, return_index=True); values=values[idx]; points=points[idx]
    if len(u)<2: return np.full(len(stations),np.nan)
    result = np.interp(target,u,values,left=np.nan,right=np.nan)
    high=np.clip(np.searchsorted(u,target),1,len(u)-1); low=high-1
    gap=(u[high]-u[low])*vessel.length_mm
    spatial=cKDTree(points).query(stations)[0]
    result[(gap>2*tolerance_mm+1e-7)|(spatial>tolerance_mm)] = np.nan
    return result


def matched(a, b, truth):
    a,b,truth = map(lambda x:np.asarray(x,float),(a,b,truth))
    if a.shape != b.shape or a.shape != truth.shape:
        raise ValueError('common_station_shape_mismatch')
    support=np.isfinite(a)&np.isfinite(b)&np.isfinite(truth)
    left,right=errors(a[support],truth[support]),errors(b[support],truth[support])
    return dict(common_sample_count=int(support.sum()), requested_sample_count=len(truth),
                common_support_fraction=float(support.mean()) if len(truth) else None,
                left_valid_count=int(np.isfinite(a).sum()),right_valid_count=int(np.isfinite(b).sum()),
                left_mae_mm=left['mae_mm'],right_mae_mm=right['mae_mm'],
                right_minus_left_mae_mm=right['mae_mm']-left['mae_mm'] if support.any() else None,
                left_rmse_mm=left['rmse_mm'],right_rmse_mm=right['rmse_mm'],
                status='matched' if support.any() else 'unscorable')
