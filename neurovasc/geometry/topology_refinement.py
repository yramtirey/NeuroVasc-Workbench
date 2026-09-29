"""Experimental geometry refinement with fixed logical topology and explicit rollback."""
from copy import deepcopy
from dataclasses import dataclass, asdict
import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
from .centerline_refinement import Refiner, RefinementParameters, resample, distance
from .centerline_qc import occupancy
from .cross_sectional_caliber import estimate_tangents
from neurovasc.graph.junction_refinement import classify_junctions, junction_center, classify_spurs
from neurovasc.graph.topology_invariants import materialize, check_invariants, RefinementResult


@dataclass(frozen=True)
class NetworkParameters:
    junction_cluster_radius_mm: float = 1.0
    resample_spacing_mm: float = .5
    smoothing_scale_mm: float = 1.0
    spur_handling: str = 'conservative'
    endpoint_mode: str = 'raw'
    recenter: bool = True


CANDIDATE = NetworkParameters()


def clean_path(points):
    points = np.asarray(points, float)
    keep = np.r_[True, np.linalg.norm(np.diff(points, axis=0), axis=1)>1e-10]
    points = points[keep]
    if len(points)<2:
        raise ValueError('collapsed_branch_geometry')
    return points


def endpoint_position(network, node, refiner, mode):
    original = np.array(network.logical.nodes[node]['world'], float)
    if mode=='raw': return original
    branch = next(b for b in network.branches.values() if node in (b['start'],b['end']))
    points = branch['points'] if branch['start']==node else branch['points'][::-1]
    sampled = resample(clean_path(points), .5)
    inward = estimate_tangents(sampled, 4.)[0]
    _, nearest = refiner.tree.query(original)
    ids = refiner.tree.query_ball_point(original, 2*refiner.radii[nearest]+max(refiner.volume.spacing))
    offsets = refiner.foreground[ids]-original
    axial = offsets@inward
    weights = refiner.radii[ids]*np.exp(-.5*(axial/.5)**2)
    shift = np.average(offsets, axis=0, weights=weights)
    shift -= (shift@inward)*inward
    center = refiner.constrain(original, original+shift, max(refiner.volume.spacing))
    if mode=='recentered': return center
    limit = 2*max(refiner.volume.spacing)
    if mode=='tangent':
        return refiner.constrain(center, center-limit*inward, limit)
    # Boundary mode stops at the first failed occupancy sample, never jumps gaps.
    chosen = center
    for step in np.arange(.05, limit+.001, .05):
        point = center-step*inward
        if occupancy(refiner.volume, point[None])[0]<.5: break
        chosen=point
    return chosen


def geometric_checks(volume, network, before):
    step = min(volume.spacing)/8
    components = list(nx.connected_components(network.logical))
    samples = []
    for component in components:
        paths = [resample(clean_path(b['points']), step) for b in network.branches.values() if b['start'] in component]
        points = np.concatenate(paths)
        if np.any(occupancy(volume, points)<.5-1e-8):
            raise ValueError('dense_segment_outside_lumen')
        samples.append(points)
    # Distinct components never gain graph edges. Also reject geometric contact
    # at the sampling resolution, even though coordinate proximity adds no edge.
    for i in range(len(samples)):
        for j in range(i):
            if cKDTree(samples[j]).query(samples[i])[0].min() <= step:
                raise ValueError('cross_component_geometric_contact')
    return True


def _refine_local(volume, original, parameters=CANDIDATE, stage='full', refiner=None):
    values = (parameters.junction_cluster_radius_mm, parameters.resample_spacing_mm, parameters.smoothing_scale_mm)
    if not np.isfinite(values).all() or values[0]<=0 or values[1]<=0 or values[2]<0:
        raise ValueError('invalid_network_refinement_parameters')
    if parameters.endpoint_mode not in ('raw','recentered','tangent','boundary') or parameters.spur_handling not in ('off','conservative') or stage not in ('junction','full'):
        raise ValueError('invalid_network_refinement_mode')
    refiner = refiner or Refiner(volume)
    diagnostics = classify_junctions(original, parameters.junction_cluster_radius_mm)
    spurs = classify_spurs(original, refiner, parameters.junction_cluster_radius_mm, parameters.spur_handling=='conservative')
    proposal = deepcopy(original)
    lengths = []
    try:
        for node, data in proposal.logical.nodes(data=True):
            if data['node_kind']=='junction':
                data['world']=junction_center(refiner, data['world'], parameters.junction_cluster_radius_mm).tolist()
            elif data['node_kind']=='endpoint' and stage=='full':
                data['world']=endpoint_position(original, node, refiner, parameters.endpoint_mode).tolist()
        for key, branch in proposal.branches.items():
            source = clean_path(original.branches[key]['points'])
            sampled = resample(source, parameters.resample_spacing_mm)
            points = source.copy()
            points[0]=proposal.logical.nodes[branch['start']]['world']
            points[-1]=proposal.logical.nodes[branch['end']]['world']
            points = clean_path(points)
            if stage=='full' and len(points)>=3:
                config=RefinementParameters(recenter=parameters.recenter, spacing_mm=parameters.resample_spacing_mm,
                    smoothing_scale_mm=parameters.smoothing_scale_mm, refine_endpoints=False, junction_guard_mm=0.)
                points,_=refiner.refine_chain(points, (False,False), config)
            # Both ends attach exactly to the same immutable logical identities.
            points[0]=proposal.logical.nodes[branch['start']]['world']
            points[-1]=proposal.logical.nodes[branch['end']]['world']
            branch['points']=clean_path(points)
            branch['length_mm']=float(distance(branch['points'])[-1])
            branch['tangents']=estimate_tangents(branch['points'], 4.)
            proposal.logical[branch['start']][branch['end']][key]['length_mm']=branch['length_mm']
            lengths.append(dict(branch_id=key, raw_branch_length_mm=float(distance(source)[-1]),
                                resampled_length_mm=float(distance(sampled)[-1]), refined_length_mm=branch['length_mm']))
        proposal.raw=materialize(proposal)
        invariants=check_invariants(original, proposal)
        if not invariants['topology_invariant_pass']:
            raise ValueError('topology_invariant_violation')
        geometric_checks(volume, proposal, original)
        proposal.metadata.update(refinement_parameters=asdict(parameters), production_eligible=False)
        return RefinementResult(original,proposal,True,'',invariants,diagnostics,spurs,lengths,asdict(parameters))
    except ValueError as error:
        # Explicitly expose rejection; metrics for this retained result must be
        # labelled, never presented as a successful refined network.
        invariants=check_invariants(original,proposal)
        invariants['proposal_accepted']=False
        return RefinementResult(original,original,False,str(error),invariants,diagnostics,spurs,lengths,asdict(parameters))


def refine_network(volume, original, parameters=CANDIDATE, stage='full', refiner=None):
    """Operate in orthonormal voxel-axis mm coordinates, then restore world frame.

    This removes floating boundary/tie differences under rigid affine relabeling.
    It does not rotate/resample the mask or change physical distances. Rounding
    at 1e-12 mm stabilizes transformations far below the geometric tolerances.
    """
    spacing=np.asarray(volume.spacing)
    axes=np.asarray(volume.affine)[:3,:3]/spacing[None,:]
    if not np.allclose(axes.T@axes,np.eye(3),atol=1e-8):
        raise ValueError('orthogonal_affine_required')
    origin=np.asarray(volume.affine)[:3,3]
    def convert(network, forward):
        result=deepcopy(network)
        def points(p):
            p=np.asarray(p,float)
            return np.round((p-origin)@axes,12) if forward else p@axes.T+origin
        for _,data in result.logical.nodes(data=True):data['world']=points(data['world']).tolist()
        for _,data in result.raw.nodes(data=True):data['world']=tuple(points(data['world']))
        for branch in result.branches.values():
            branch['points']=points(branch['points'])
            if 'tangents' in branch:branch['tangents']=branch['tangents']@axes if forward else branch['tangents']@axes.T
        return result
    if refiner is not None and hasattr(refiner,'_topology_local_frame'):
        local_volume,local_refiner=refiner._topology_local_frame
    else:
        from copy import copy
        local_volume=copy(volume)
        local_volume.affine=np.diag([*spacing,1.])
        local_refiner=Refiner(local_volume)
        if refiner is not None:refiner._topology_local_frame=(local_volume,local_refiner)
    local=convert(original,True)
    result=_refine_local(local_volume,local,parameters,stage,local_refiner)
    result.original=original
    if result.accepted:
        result.network=convert(result.network,False)
        result.invariants=check_invariants(original,result.network)
        if not result.invariants['topology_invariant_pass']:
            result.accepted=False;result.failure_reason='world_frame_invariant_violation';result.network=original
    else:result.network=original
    return result
