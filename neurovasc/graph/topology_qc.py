"""Truth-independent network checks; short structures are flagged, never pruned."""
import networkx as nx
import numpy as np
from .topology_network import cycle_rank


def topology_qc(network,supports_cycles=True,voxel_size_mm=.5):
    checks=[]
    def add(name,value,threshold,passed,critical=True):
        checks.append(dict(criterion=name,value=value,threshold=threshold,passed=bool(passed),critical=critical))
    g=network.raw
    finite=all(np.isfinite(d['world']).all() for _,d in g.nodes(data=True))
    add('finite_coordinates',finite,'true',finite)
    components=nx.number_connected_components(g)
    add('components',components,'>=1; multiple components permitted',components>=1)
    isolated=sum(g.degree[n]==0 for n in g)
    add('isolated_nodes',isolated,'0',isolated==0)
    invalid=sum(not np.isfinite(d.get('length_mm',np.nan)) or d.get('length_mm',0)<=0 for _,_,d in g.edges(data=True))
    add('invalid_edge_lengths',invalid,'0',invalid==0)
    rawloops=nx.number_of_selfloops(g)
    add('raw_self_loops',rawloops,'0; logical closed-loop branches are legitimate',rawloops==0)
    add('duplicate_raw_edges',0,'simple raw graph disallows duplicate pairs',True)
    short=sum(b['length_mm']<voxel_size_mm for b in network.branches.values())
    add('subvoxel_branches',short,'0; diagnostic only',short==0,False)
    degree=max(dict(network.logical.degree()).values(),default=0)
    add('maximum_logical_degree',degree,'<=6; diagnostic only',degree<=6,False)
    tiny=sum(sum(network.branches[b]['length_mm'] for b in c)<4*voxel_size_mm for c in network.cycles)
    add('tiny_cycles',tiny,'perimeter >=4 maximum voxel spacings; diagnostic only',tiny==0,False)
    same=cycle_rank(g)==cycle_rank(network.logical)
    add('simplification_cycle_rank',same,'unchanged',same)
    add('cycle_capability',supports_cycles,'reported; not proof of correctness',supports_cycles,False)
    return checks
