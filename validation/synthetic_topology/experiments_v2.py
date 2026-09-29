"""Exact archived cohort plus independent harder analytic junction networks."""
from pathlib import Path
import networkx as nx
import nibabel as nib
import numpy as np
from neurovasc.io.volume import Volume
from validation.synthetic_topology.geometries import Phantom, matrix

SPACINGS=((.5,.5,.5),(.46875,.46875,.8))


def junction_phantom(kind, spacing=SPACINGS[0], phase=0.):
    graph=nx.MultiGraph();branches={}
    def node(key,point): graph.add_node(key,world=list(map(float,point)))
    def branch(a,b,vertices=None,diameter=2.,bridge=False):
        vertices=np.array(vertices if vertices is not None else [graph.nodes[a]['world'],graph.nodes[b]['world']],float)
        lengths=np.linalg.norm(np.diff(vertices,axis=0),axis=1)
        points=np.concatenate([p+np.linspace(0,1,int(np.ceil(length/.08)),endpoint=False)[:,None]*(q-p) for p,q,length in zip(vertices[:-1],vertices[1:],lengths)]+[vertices[-1:]])
        key=f'T{len(branches):03d}'
        branches[key]=dict(branch_id=key,start=a,end=b,points=points,length_mm=float(lengths.sum()),diameter_mm=diameter,vertices=vertices.tolist(),shape='polyline',communicating_bridge=bridge)
        graph.add_edge(a,b,key=key,length_mm=float(lengths.sum()))
    if kind=='four_way':
        node('J',[0,0,0])
        for i,p in enumerate(([-8,0,0],[8,0,0],[0,0,-8],[0,0,8])):
            node(f'E{i}',p);branch('J',f'E{i}')
    elif kind=='double_junction':
        node('L',[-1.5,0,0]);node('R',[1.5,0,0]);branch('L','R',bridge=True)
        for i,(a,p) in enumerate((('L',[-8,0,-6]),('L',[-8,0,6]),('R',[8,0,-6]),('R',[8,0,6]))):
            node(f'E{i}',p);branch(a,f'E{i}')
    elif kind=='asymmetric_figure_eight':
        node('J',[0,0,0])
        for scale,side in ((1.,-1),(.65,1)):
            v=np.array([[0,0,0],[side*6,0,6],[side*12,0,0],[side*6,0,-6],[0,0,0]])*scale
            branch('J','J',v,diameter=2.)
    elif kind=='short_communicator':
        node('L',[-.75,0,0]);node('R',[.75,0,0])
        for side in (-1,1):
            branch('L','R',[[-.75,0,0],[-5,0,side*4],[-5,0,side*8],[5,0,side*8],[5,0,side*4],[.75,0,0]],diameter=1.5)
        branch('L','R',diameter=1.,bridge=True)
    elif kind=='near_junction_branch':
        node('J',[0,0,0]);node('K',[3,0,0]);branch('J','K')
        for i,(a,p,d) in enumerate((('J',[-8,0,0],2.),('J',[0,0,8],2.),('K',[10,0,0],2.),('K',[3,0,-2],1.))):
            node(f'E{i}',p);branch(a,f'E{i}',diameter=d)
    else: raise ValueError('unknown_junction_phantom')
    for n in graph: graph.nodes[n]['node_kind']='endpoint' if graph.degree[n]==1 else 'junction' if graph.degree[n]>=3 else 'regular'
    spacing=np.array(spacing);points=np.concatenate([b['points'] for b in branches.values()])
    half=np.ceil((np.max(np.abs(points),axis=0)+4)/spacing).astype(int)
    affine=np.diag([*spacing,1.]);affine[:3,3]=-half*spacing+phase*spacing
    shape=tuple(2*half+1);world=nib.affines.apply_affine(affine,np.moveaxis(np.indices(shape),0,-1))
    mask=np.zeros(shape,bool)
    for b in branches.values():
        vertices=np.array(b['vertices'])
        for p,q in zip(vertices[:-1],vertices[1:]):
            delta=q-p;t=np.clip((world-p)@delta/(delta@delta),0,1)
            mask|=np.linalg.norm(world-p-t[...,None]*delta,axis=-1)<=b['diameter_mm']/2+1e-10
    volume=Volume(mask.astype(np.uint8),affine,tuple(spacing),Path('synthetic_junction.nii.gz'))
    return Phantom(volume,graph,branches,kind,phase,'canonical')


def matrix_v2(archive=None):
    for eid,phantom in matrix():
        if archive is not None:
            image=nib.load(Path(archive)/'masks'/f'{eid}.nii.gz')
            if not np.array_equal(image.dataobj,phantom.volume.data) or not np.allclose(image.affine,phantom.volume.affine,atol=1e-6):
                raise ValueError('v1_mask_or_affine_mismatch:'+eid)
        yield eid,phantom,'core_v1'
    for spacing in SPACINGS:
        for phase in (0.,.25,.5):
            for kind in ('four_way','double_junction','asymmetric_figure_eight','short_communicator','near_junction_branch'):
                yield f'J_{kind}_z{spacing[2]}_phase{phase}',junction_phantom(kind,spacing,phase),'junction_added'
