"""Strict identity-aware invariants for experimental geometry-only refinement."""
from dataclasses import dataclass
from copy import deepcopy
import networkx as nx
import numpy as np
from .topology_network import TopologyNetwork, cycle_rank


def adjacency(network):
    return {key: tuple(sorted((b['start'], b['end']))) for key, b in network.branches.items()}


def materialize(network):
    """Embed branch paths without inferring connectivity from coordinate proximity."""
    graph = nx.Graph()
    for node, data in network.logical.nodes(data=True):
        graph.add_node(('node', node), world=tuple(data['world']))
    for key, branch in network.branches.items():
        points = np.asarray(branch['points'])
        ids = [('node', branch['start'])] + [('sample', key, i) for i in range(1, len(points)-1)] + [('node', branch['end'])]
        if len(points) < 3:
            # Every parallel branch needs a distinct interior identity.
            points = np.vstack((points[0], (points[0]+points[-1])/2, points[-1]))
            ids = [ids[0], ('sample', key, 1), ids[-1]]
        for node, point in zip(ids, points):
            graph.add_node(node, world=tuple(point))
        for a, b in zip(ids[:-1], ids[1:]):
            length = float(np.linalg.norm(np.asarray(graph.nodes[a]['world'])-graph.nodes[b]['world']))
            if a == b or length <= 1e-10:
                raise ValueError('collapsed_branch_geometry')
            graph.add_edge(a, b, length_mm=length, branch_id=key)
    return graph


def check_invariants(before, after):
    result = {}
    for name, network in (('before', before), ('after', after)):
        for metric, value in network.counts().items():
            if metric in ('components', 'cycle_rank', 'branches', 'endpoints', 'junctions'):
                result[f'{metric}_{name}'] = value
        result[f'logical_components_{name}'] = nx.number_connected_components(network.logical)
    result['branch_identity_preserved'] = set(before.branches) == set(after.branches)
    result['node_identity_preserved'] = set(before.logical) == set(after.logical)
    result['adjacency_preserved'] = adjacency(before) == adjacency(after)
    result['logical_edges_preserved'] = sorted(before.logical.edges(keys=True)) == sorted(after.logical.edges(keys=True))
    result['cycle_membership_preserved'] = before.cycles == after.cycles and all(
        key in after.branches and b['cycle_membership'] == after.branches[key]['cycle_membership']
        for key, b in before.branches.items())
    result['node_kinds_preserved'] = all(n in after.logical and d['node_kind'] == after.logical.nodes[n]['node_kind'] for n, d in before.logical.nodes(data=True))
    result['branch_attachments_exact'] = all(
        np.allclose(b['points'][0], after.logical.nodes[b['start']]['world'], atol=1e-10, rtol=0) and
        np.allclose(b['points'][-1], after.logical.nodes[b['end']]['world'], atol=1e-10, rtol=0)
        for b in after.branches.values())
    result['finite_geometry'] = all(np.isfinite(b['points']).all() for b in after.branches.values())
    result['embedded_rank_consistent'] = cycle_rank(after.raw) == cycle_rank(after.logical)
    try:
        expected = materialize(after)
        result['embedded_edges_preserved'] = (set(expected) == set(after.raw) and
            {frozenset(e) for e in expected.edges} == {frozenset(e) for e in after.raw.edges} and
            all(np.allclose(d['world'], after.raw.nodes[n]['world'], atol=1e-10, rtol=0) for n,d in expected.nodes(data=True)))
    except (ValueError, KeyError):
        result['embedded_edges_preserved'] = False
    conserved = all(result[f'{k}_before'] == result[f'{k}_after'] for k in ('components', 'logical_components', 'cycle_rank', 'branches', 'endpoints', 'junctions'))
    result['topology_invariant_pass'] = conserved and all(v for k,v in result.items() if isinstance(v, (bool, np.bool_)))
    return result


@dataclass
class RefinementResult:
    original: TopologyNetwork
    network: TopologyNetwork
    accepted: bool
    failure_reason: str
    invariants: dict
    junction_diagnostics: list
    spur_diagnostics: list
    lengths: list
    parameters: dict
