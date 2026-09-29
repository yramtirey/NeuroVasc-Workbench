"""Experimental distance-ordered (6 foreground, 26 background) thinning.

No Lee, geodesic tree or analytic-truth calls. Retains the full digital graph.
"""
from functools import lru_cache
import heapq
from itertools import product
import networkx as nx
import nibabel as nib
import numpy as np
from scipy.ndimage import distance_transform_edt, generate_binary_structure, label

FACES=np.array([(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)])
OFFSETS=np.array([p for p in product((-1,0,1),repeat=3) if p!=(0,0,0)])
SIX=generate_binary_structure(3,1)
FULL=np.ones((3,3,3),bool)
N18=generate_binary_structure(3,2);N18[1,1,1]=False


@lru_cache(maxsize=100_000)
def simple_point(bits):
    block=np.unpackbits(np.frombuffer(bits,dtype=np.uint8))[:27].reshape(3,3,3).astype(bool)
    # Preserve terminal voxels using the entire immediate neighborhood.
    if block.sum()<=2:
        return False
    foreground=block & N18
    components,_=label(foreground,SIX)
    touching={int(components[tuple(1+offset)]) for offset in FACES if foreground[tuple(1+offset)]}
    if len(touching)!=1:
        return False
    background=~block
    background[1,1,1]=False
    return label(background,FULL)[1]==1


def extract_topology_network(volume, maximum_voxels=150_000):
    original=np.asarray(volume.data)>.5
    if original.ndim!=3 or not original.any():raise ValueError('empty_or_non3d_mask')
    if original.sum()>maximum_voxels:raise ValueError('topology_size_limit')
    axes=np.asarray(volume.affine)[:3,:3]
    if not np.isfinite(volume.affine).all() or not np.allclose(axes.T@axes,np.diag(np.asarray(volume.spacing)**2),atol=1e-8):
        raise ValueError('orthogonal_affine_with_matching_spacing_required')
    mask=np.pad(original,1)
    radius=distance_transform_edt(mask,sampling=volume.spacing)
    heap=[];queued=set()
    def enqueue(p):
        key=tuple(int(v) for v in p)
        if key not in queued and mask[key] and any(not mask[tuple(p+d)] for d in FACES):
            heapq.heappush(heap,(float(radius[key]),key));queued.add(key)
    for p in np.argwhere(mask):enqueue(p)
    deleted=0
    while heap:
        _,key=heapq.heappop(heap);queued.remove(key)
        if not mask[key]:continue
        x,y,z=key
        block=mask[x-1:x+2,y-1:y+2,z-1:z+2]
        if not simple_point(np.packbits(block).tobytes()):continue
        mask[key]=False;deleted+=1
        for d in OFFSETS:enqueue(np.array(key)+d)
    graph=nx.Graph()
    for p in np.argwhere(mask):
        voxel=tuple(int(v) for v in p-1)
        world=nib.affines.apply_affine(volume.affine,voxel)
        graph.add_node(voxel,world=tuple(world),voxel=voxel,clearance_mm=float(radius[tuple(p)]))
    for a in graph:
        for d in ((1,0,0),(0,1,0),(0,0,1)):
            b=tuple(np.array(a)+d)
            if b in graph:
                graph.add_edge(a,b,length_mm=float(np.linalg.norm(np.asarray(graph.nodes[a]['world'])-graph.nodes[b]['world'])))
    before=label(original,SIX)[1];after=nx.number_connected_components(graph)
    if before!=after:raise RuntimeError('thinning_component_invariant_failed')
    for n in graph:graph.nodes[n]['node_kind']='endpoint' if graph.degree[n]==1 else 'junction' if graph.degree[n]>=3 else 'regular'
    graph.graph.update(supports_cycles=True,production_eligible=False,foreground_connectivity=6,background_connectivity=26,
                       deleted_voxels=deleted,input_components_6=before,pruned_cycles=[],method='distance_ordered_simple_point_thinning')
    return graph
