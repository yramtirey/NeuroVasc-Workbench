"""Lossless cycle-rank branch reduction using a MultiGraph, including closed loops."""
from dataclasses import dataclass
import networkx as nx
import numpy as np


def cycle_rank(g):
    return g.number_of_edges()-g.number_of_nodes()+nx.number_connected_components(g)


@dataclass
class TopologyNetwork:
    raw: nx.Graph
    logical: nx.MultiGraph
    branches: dict
    cycles: list
    metadata: dict

    def counts(self):
        return {'components':nx.number_connected_components(self.raw),
                'endpoints':sum(d['node_kind']=='endpoint' for _,d in self.logical.nodes(data=True)),
                'junctions':sum(d['node_kind']=='junction' for _,d in self.logical.nodes(data=True)),
                'branches':len(self.branches),'cycle_rank':cycle_rank(self.logical),
                'nodes':len(self.raw),'edges':self.raw.number_of_edges()}


def fundamental_cycles(logical):
    forest=nx.Graph();forest.add_nodes_from(logical)
    cycles=[]
    for a,b,key in sorted(logical.edges(keys=True),key=lambda edge:edge[2]):
        if a==b:cycles.append([key])
        elif nx.has_path(forest,a,b):
            route=nx.shortest_path(forest,a,b)
            cycles.append(sorted([key]+[forest[x][y]['branch_id'] for x,y in zip(route[:-1],route[1:])]))
        else:forest.add_edge(a,b,branch_id=key)
    return cycles


def reduce_network(raw):
    if not len(raw):raise ValueError('empty_recovered_network')
    logical=nx.MultiGraph();mapping={};branches={};visited=set()
    junctions=[n for n in raw if raw.degree[n]>=3]
    clusters=list(nx.connected_components(raw.subgraph(junctions)))
    for i,cluster in enumerate(clusters):
        key=f'J{i}'
        logical.add_node(key,node_kind='junction',world=np.mean([raw.nodes[n]['world'] for n in cluster],axis=0).tolist(),raw_nodes=list(cluster))
        for n in cluster:mapping[n]=key
    for n in raw:
        if raw.degree[n]<=1:
            key=f'E{len(logical)}'
            logical.add_node(key,node_kind='endpoint' if raw.degree[n]==1 else 'regular',world=list(raw.nodes[n]['world']),raw_nodes=[n]);mapping[n]=key
    # A degree-two closed component still needs an anchor and a self-loop branch.
    for component in nx.connected_components(raw):
        if not any(n in mapping for n in component):
            n=next(n for n in raw if n in component);key=f'R{len(logical)}'
            logical.add_node(key,node_kind='regular',world=list(raw.nodes[n]['world']),raw_nodes=[n]);mapping[n]=key
    def add_branch(a,b,path,kind='chain'):
        key=f'B{len(branches):04d}'
        points=np.array([raw.nodes[n]['world'] for n in path])
        length=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())
        branches[key]={'branch_id':key,'start':a,'end':b,'points':points,'length_mm':length,'kind':kind}
        logical.add_edge(a,b,key=key,branch_id=key,length_mm=length)
    # Contract only a spanning tree of each logical junction cluster. Every
    # remaining internal edge is retained as a diagnostic self-loop cycle.
    contracted=0
    for cluster in clusters:
        sub=raw.subgraph(cluster)
        tree=nx.minimum_spanning_tree(sub,weight='length_mm')
        tree_edges={frozenset(e) for e in tree.edges}
        for a,b in sub.edges:
            visited.add(frozenset((a,b)))
            if frozenset((a,b)) in tree_edges:contracted+=1
            else:
                route=nx.shortest_path(tree,b,a)
                add_branch(mapping[a],mapping[a],[a]+route,'junction_internal_cycle')
    for start in mapping:
        for neighbor in raw.neighbors(start):
            edge=frozenset((start,neighbor))
            if edge in visited:continue
            path=[start,neighbor];visited.add(edge);previous,current=start,neighbor
            while current not in mapping:
                candidates=[n for n in raw.neighbors(current) if n!=previous]
                if len(candidates)!=1:raise RuntimeError('unclassified_branch_node')
                nxt=candidates[0];visited.add(frozenset((current,nxt)));path.append(nxt)
                previous,current=current,nxt
            add_branch(mapping[start],mapping[current],path)
    if len(visited)!=raw.number_of_edges():raise RuntimeError('unvisited_network_edge')
    if cycle_rank(raw)!=cycle_rank(logical):raise RuntimeError('branch_reduction_changed_cycle_rank')
    cycles=fundamental_cycles(logical)
    for key,branch in branches.items():
        branch['cycle_membership']=[i for i,c in enumerate(cycles) if key in c]
        branch['adjacent_branches']=sorted({k for n in (branch['start'],branch['end']) for _,_,k in logical.edges(n,keys=True) if k!=key})
    return TopologyNetwork(raw,logical,branches,cycles,{'contracted_junction_tree_edges':contracted,'cycles_pruned':[],
                                                    'raw_cycle_rank':cycle_rank(raw),'logical_cycle_rank':cycle_rank(logical)})
