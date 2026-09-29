"""Analytic multigraph tubes, voxelized in physical mm; truth is evaluation-only."""
from dataclasses import dataclass
from pathlib import Path
import networkx as nx
import nibabel as nib
import numpy as np
from scipy.spatial import cKDTree
from neurovasc.io.volume import Volume
from neurovasc.graph.topology_network import fundamental_cycles
from validation.synthetic_centerline.geometries import y_bifurcation


@dataclass
class Phantom:
    volume: Volume
    graph: nx.MultiGraph
    branches: dict
    geometry: str
    phase: float
    orientation: str
    bridge_diameter_mm: float | None = None
    minimum_separation_mm: float | None = None

    def metadata(self):
        return {'geometry':self.geometry,'phase_voxels':[self.phase]*3,'orientation':self.orientation,
                'spacing_mm':list(self.volume.spacing),'affine':self.volume.affine.tolist(),'shape':list(self.volume.shape),
                'bridge_diameter_mm':self.bridge_diameter_mm,'minimum_separation_mm':self.minimum_separation_mm,
                'nodes':{n:d for n,d in self.graph.nodes(data=True)},
                'branches':{k:{kk:vv.tolist() if isinstance(vv,np.ndarray) else vv for kk,vv in b.items()} for k,b in self.branches.items()},
                'true_components':nx.number_connected_components(self.graph),'true_cycle_rank':len(self.graph.edges)-len(self.graph)+nx.number_connected_components(self.graph),
                'true_endpoints':sum(d['node_kind']=='endpoint' for _,d in self.graph.nodes(data=True)),
                'true_junctions':sum(d['node_kind']=='junction' for _,d in self.graph.nodes(data=True)),
                'true_branch_count':len(self.branches),'true_length_mm':sum(b['length_mm'] for b in self.branches.values()),
                'cycle_basis':fundamental_cycles(self.graph),'voxelization':'binary tube inclusion at voxel centers; circles use exact finite-arc distance, lines finite-segment distance; round caps except reused Y'}


def phantom(kind,spacing=(.5,.5,.5),phase=0.,bridge_diameter=1.5,separation=.4,rotation=False):
    g=nx.MultiGraph();branches={}
    def node(name,p):g.add_node(name,world=list(map(float,p)))
    def line(a,b,diameter=3.,bridge=False):
        p,q=np.array(g.nodes[a]['world']),np.array(g.nodes[b]['world'])
        length=np.linalg.norm(q-p);t=np.linspace(0,1,int(np.ceil(length/.08))+1)
        add(a,b,p+t[:,None]*(q-p),diameter,length,'line',bridge)
    def add(a,b,points,diameter,length,shape,bridge=False,**extra):
        key=f'T{len(branches):03d}'
        branches[key]=dict(branch_id=key,start=a,end=b,points=points,diameter_mm=diameter,length_mm=float(length),shape=shape,communicating_bridge=bridge,**extra)
        g.add_edge(a,b,key=key,length_mm=float(length))
    def arc(a,b,center,r,start,end,diameter=3.):
        theta=np.linspace(start,end,int(np.ceil(r*(end-start)/.08))+1)
        points=np.column_stack((center[0]+r*np.cos(theta),np.zeros(len(theta)),center[2]+r*np.sin(theta)))
        add(a,b,points,diameter,r*(end-start),'arc',center=list(center),radius_mm=r,start_angle=start,end_angle=end)
    if kind in ('simple_ring','ring_branches','ring_multiple'):
        count={'simple_ring':1,'ring_branches':2,'ring_multiple':4}[kind]
        for i in range(count):
            theta=2*np.pi*i/count;node(f'J{i}',[8*np.cos(theta),0,8*np.sin(theta)])
        for i in range(count):arc(f'J{i}',f'J{(i+1)%count}',(0,0,0),8,2*np.pi*i/count,2*np.pi*(i+1)/count)
        if count>1:
            for i in range(count):node(f'E{i}',np.array(g.nodes[f'J{i}']['world'])*1.75);line(f'J{i}',f'E{i}')
    elif kind=='figure_eight':
        # Two diamond loops sharing one four-way junction; straight tubes avoid
        # tangentially touching circle boundaries merging into a broad lens.
        node('J',[0,0,0])
        for side in (-1,1):
            vertices=np.array([[0,0,0],[side*6,0,6],[side*12,0,0],[side*6,0,-6],[0,0,0]],float)
            length=np.linalg.norm(np.diff(vertices,axis=0),axis=1).sum()
            points=np.concatenate([a+np.linspace(0,1,110,endpoint=False)[:,None]*(b-a) for a,b in zip(vertices[:-1],vertices[1:])]+[vertices[-1:]])
            add('J','J',points,3.,length,'polyline',vertices=vertices.tolist())
    elif kind in ('communicating_bridge','small_bridge'):
        # Two parallel rails joined at their ends, with a short middle bridge.
        node('L',[-3,0,0]);node('R',[3,0,0])
        for side in (-1,1):
            vertices=np.array([[-3,0,0],[-3,0,side*8],[3,0,side*8],[3,0,0]],float)
            points=np.concatenate([a+np.linspace(0,1,int(np.ceil(np.linalg.norm(z-a)/.08)),endpoint=False)[:,None]*(z-a) for a,z in zip(vertices[:-1],vertices[1:])]+[vertices[-1:]])
            add('L','R',points,3.,22.,'polyline',vertices=vertices.tolist())
        line('L','R',bridge_diameter if kind=='small_bridge' else 3.,bridge=True)
    elif kind=='incomplete_loop':
        angle=.35
        node('A',[8*np.cos(angle),0,8*np.sin(angle)]);node('B',[8*np.cos(2*np.pi-angle),0,8*np.sin(2*np.pi-angle)])
        arc('A','B',(0,0,0),8,angle,2*np.pi-angle)
    elif kind=='y':
        y=y_bifurcation(spacing)
        node('J',[0,0,0])
        for name,segment in y.branches.items():node(name,segment[1]);line('J',name,y.diameters[name])
    elif kind=='near_touch':
        a=(2.+separation)/(2*np.sqrt(2))
        for i,sign in enumerate((-1,1)):
            node(f'A{i}',[sign*a,sign*a,-12]);node(f'B{i}',[sign*a,sign*a,12]);line(f'A{i}',f'B{i}',2.)
    else:raise ValueError('unknown_phantom')
    for n in g:g.nodes[n]['node_kind']='endpoint' if g.degree[n]==1 else 'junction' if g.degree[n]>=3 else 'regular'
    transform=np.eye(3)
    if rotation:
        # Physical orientation changes relative to the voxel lattice, not merely
        # a relabeling affine. All analytic branch points rotate consistently.
        ay,az=.43,.31
        ry=np.array([[np.cos(ay),0,np.sin(ay)],[0,1,0],[-np.sin(ay),0,np.cos(ay)]])
        rz=np.array([[np.cos(az),-np.sin(az),0],[np.sin(az),np.cos(az),0],[0,0,1]])
        transform=rz@ry
        for n in g:g.nodes[n]['world']=(transform@np.array(g.nodes[n]['world'])).tolist()
    all_points=np.concatenate([b['points']@transform.T for b in branches.values()])
    maxradius=max(b['diameter_mm'] for b in branches.values())/2
    spacing=np.array(spacing)
    half=np.ceil((np.max(np.abs(all_points),axis=0)+maxradius+3*spacing.max())/spacing).astype(int)
    affine=np.diag([*spacing,1.]);affine[:3,3]=-half*spacing+phase*spacing
    shape=tuple(2*half+1)
    world=nib.affines.apply_affine(affine,np.moveaxis(np.indices(shape),0,-1))
    canonical=world@transform
    mask=np.zeros(shape,bool)
    def segment_distance(points,a,b):
        d=b-a;t=np.clip((points-a)@d/(d@d),0,1)
        return np.linalg.norm(points-a-t[...,None]*d,axis=-1)
    for b in branches.values():
        if b['shape']=='arc':
            center=np.array(b['center']);rel=canonical-center
            theta=np.arctan2(rel[...,2],rel[...,0])
            theta=(theta-b['start_angle'])%(2*np.pi)+b['start_angle']
            full=b['end_angle']-b['start_angle']>=2*np.pi-1e-9
            radial=np.hypot(np.hypot(rel[...,0],rel[...,2])-b['radius_mm'],rel[...,1])
            ends=np.minimum(np.linalg.norm(canonical-b['points'][0],axis=-1),np.linalg.norm(canonical-b['points'][-1],axis=-1))
            dist=np.where(full|(theta<=b['end_angle']),radial,ends)
        else:
            vertices=np.array(b.get('vertices',[b['points'][0],b['points'][-1]]))
            dist=np.minimum.reduce([segment_distance(canonical,a,z) for a,z in zip(vertices[:-1],vertices[1:])])
        mask|=dist<=b['diameter_mm']/2+1e-10
        b['points']=b['points']@transform.T
    if kind=='y' and not rotation:
        # Reuse the original flat-capped Y inclusion for every phase.
        mask[:]=False
        for b in branches.values():
            a,z=b['points'][[0,-1]];d=(z-a)/b['length_mm'];along=(world-a)@d
            radial=np.linalg.norm(world-a-along[...,None]*d,axis=-1)
            mask|=(along>=-1e-10)&(along<=b['length_mm']+1e-10)&(radial<=b['diameter_mm']/2+1e-10)
    volume=Volume(mask.astype(np.uint8),affine,tuple(spacing),Path('synthetic_topology.nii.gz'))
    return Phantom(volume,g,branches,kind,phase,'rotated' if rotation else 'canonical',bridge_diameter if kind=='small_bridge' else None,separation if kind=='near_touch' else None)


def matrix():
    specs=[(k,{}) for k in ('simple_ring','ring_branches','figure_eight','communicating_bridge','incomplete_loop','y','ring_multiple')]
    specs += [('near_touch',{'separation':s}) for s in (.4,1.2)]
    specs += [('small_bridge',{'bridge_diameter':d}) for d in (1.,1.5,2.)]
    index=0
    for spacing in ((.5,.5,.5),(.46875,.46875,.8)):
        for phase in (0.,.25,.5):
            for kind,kwargs in specs:
                index+=1;v=phantom(kind,spacing,phase,**kwargs)
                yield f'{index:03d}_{kind}_z{spacing[2]}_phase{phase}'+(''.join(f'_{k}{val}' for k,val in kwargs.items())),v
        for kind in ('simple_ring','communicating_bridge'):
            index+=1;yield f'{index:03d}_{kind}_z{spacing[2]}_rotated',phantom(kind,spacing,rotation=True)
