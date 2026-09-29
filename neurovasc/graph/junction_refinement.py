"""Conservative, truth-independent junction classification and shared coordinates."""
import networkx as nx
import numpy as np


def classify_junctions(network, radius_mm):
    """Describe local incidence; do not contract an ambiguous real connector."""
    records = []
    simple = nx.Graph(network.logical)
    bridges = {frozenset(edge) for edge in nx.bridges(simple)}
    for node, data in network.logical.nodes(data=True):
        if data['node_kind'] != 'junction':
            continue
        incident = [b for b in network.branches.values() if node in (b['start'], b['end'])]
        cycle_ids = sorted({c for b in incident for c in b['cycle_membership']})
        for branch in incident:
            other = branch['end'] if branch['start'] == node else branch['start']
            if other == node or network.logical.nodes[other]['node_kind'] != 'junction':
                continue
            if node > other or branch['length_mm'] > 2*radius_mm:
                continue
            connector = not branch['cycle_membership'] and frozenset((node, other)) in bridges
            records.append(dict(node_id=node, other_node=other, branch_id=branch['branch_id'],
                                junction_degree=network.logical.degree[node], cycle_ids=cycle_ids,
                                classification='ambiguous_split_or_true_connector' if connector else 'cycle_passage',
                                decision='retain_distinct_junctions', merge_performed=False,
                                reason='Incidence cannot distinguish a split from a real short connector; no truth-based merge'))
        records.append(dict(node_id=node, other_node=None, branch_id=None,
                            junction_degree=network.logical.degree[node], cycle_ids=cycle_ids,
                            classification='high_degree_junction' if network.logical.degree[node] >= 4 else 'branch_junction',
                            decision='shared_geometry_center_only', merge_performed=False,
                            reason='Preserve identities and all incident chains'))
    return records


def junction_center(refiner, coordinate, radius_mm):
    point = np.asarray(coordinate, float)
    ids = refiner.tree.query_ball_point(point, radius_mm)
    if not ids:
        return point
    offsets = refiner.foreground[ids]-point
    weights = refiner.radii[ids]**2 * np.exp(-np.sum(offsets**2, axis=1)/(2*radius_mm**2))
    candidate = np.average(refiner.foreground[ids], axis=0, weights=weights)
    return refiner.constrain(point, candidate, min(radius_mm, max(refiner.volume.spacing)))


def classify_spurs(network, refiner, radius_mm, enabled=True):
    records = []
    for key, branch in network.branches.items():
        terminals = [n for n in (branch['start'], branch['end']) if network.logical.nodes[n]['node_kind']=='endpoint']
        if not terminals:
            continue
        endpoint = np.asarray(network.logical.nodes[terminals[0]]['world'])
        other = branch['end'] if branch['start']==terminals[0] else branch['start']
        junction = network.logical.nodes[other]
        _, i = refiner.tree.query(endpoint)
        clearance = float(refiner.radii[i])
        _, j = refiner.tree.query(junction['world'])
        reasons = []
        if branch['length_mm'] < 2*max(refiner.volume.spacing): reasons.append('sub_two_voxel_length')
        if not branch['cycle_membership']: reasons.append('no_cycle_membership')
        if junction['node_kind']=='junction' and np.linalg.norm(endpoint-np.asarray(junction['world'])) < refiner.radii[j]: reasons.append('inside_junction_clearance_region')
        if clearance < .75*refiner.radii[j]: reasons.append('low_relative_clearance')
        candidate = enabled and len(reasons)>=3
        records.append(dict(branch_id=key, spur_candidate=candidate, spur_score=len(reasons)/4,
                            spur_reasons=';'.join(reasons), spur_removed=False,
                            topology_effect_if_removed='branch and endpoint identities lost; rank/components may remain unchanged',
                            decision='mark_only_identity_contract' if candidate else 'retain',
                            clearance_mm=clearance, length_mm=branch['length_mm']))
    return records
