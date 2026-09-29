"""Frozen v1 matrix plus analytic ring networks; truth never enters selection."""
from dataclasses import dataclass
from pathlib import Path
import nibabel as nib
import numpy as np
from neurovasc.io.volume import Volume
from validation.synthetic_caliber.experiments_v2 import Experiment
from .geometries import matrix as matrix_v1


@dataclass
class LoopVessel:
    volume: Volume
    geometry: str
    branches: dict
    radius_mm: float = 8.
    proximal_diameter_mm: float = 3.
    distal_diameter_mm: float = 3.

    @property
    def length_mm(self):
        return 2 * np.pi * self.radius_mm + sum(np.linalg.norm(s[1] - s[0]) for s in self.branches.values())

    def metadata(self):
        return {'geometry': self.geometry, 'ring_radius_mm': self.radius_mm, 'diameter_mm': self.proximal_diameter_mm,
                'branches': {k:v.tolist() for k,v in self.branches.items()}, 'true_network_length_mm': self.length_mm,
                'true_cycles': 1, 'true_endpoints': len(self.branches), 'true_junctions': len(self.branches)}


def loop_vessel(spacing=(.5,.5,.5), attached=False):
    spacing=np.asarray(spacing)
    radius, tube_radius=8.,1.5
    branches={'left':np.array([[-8.,0,0],[-16.,0,0]]),'right':np.array([[8.,0,0],[16.,0,0]])} if attached else {}
    half=np.ceil((np.array([18.,4.,10.])+3*spacing.max())/spacing).astype(int)
    affine=np.diag([*spacing,1.])
    affine[:3,3]=-half*spacing
    shape=tuple(2*half+1)
    points=nib.affines.apply_affine(affine,np.moveaxis(np.indices(shape),0,-1))
    mask=(np.hypot(np.hypot(points[...,0],points[...,2])-radius,points[...,1])<=tube_radius+1e-10)
    for segment in branches.values():
        direction=(segment[1]-segment[0])/8
        along=(points-segment[0])@direction
        radial=np.linalg.norm(points-segment[0]-along[...,None]*direction,axis=-1)
        mask|=(along>=0)&(along<=8)&(radial<=tube_radius+1e-10)
    return LoopVessel(Volume(mask.astype(np.uint8),affine,tuple(spacing),Path('loop.nii.gz')),'loop_branches' if attached else 'loop',branches)


def matrix_v2():
    yield from matrix_v1()
    index=61
    for spacing in ((.5,.5,.5),(.46875,.46875,.8)):
        for attached in (False,True):
            vessel=loop_vessel(spacing,attached)
            yield f'{index}_{vessel.geometry}_z{spacing[2]}',Experiment(vessel,cohort='loops',orientation_id='ring_xz')
            index+=1


def loop_truth(vessel, points):
    points=np.asarray(points)
    theta=np.arctan2(points[:,2],points[:,0])
    ring=np.column_stack((vessel.radius_mm*np.cos(theta),np.zeros(len(points)),vessel.radius_mm*np.sin(theta)))
    tangent=np.column_stack((-np.sin(theta),np.zeros(len(points)),np.cos(theta)))
    distances=[np.linalg.norm(points-ring,axis=1)]
    candidates=[ring]
    tangents=[tangent]
    from .truth import segment_projection
    for segment in vessel.branches.values():
        mapped,t,unused=segment_projection(points,segment)
        candidates.append(mapped); tangents.append(t); distances.append(np.linalg.norm(points-mapped,axis=1))
    selected=np.argmin(distances,axis=0)
    return np.array([candidates[k][i] for i,k in enumerate(selected)]),np.array([tangents[k][i] for i,k in enumerate(selected)])


def loop_samples(vessel, step=.1):
    theta=np.linspace(0,2*np.pi,int(np.ceil(2*np.pi*vessel.radius_mm/step)),endpoint=False)
    ring=np.column_stack((vessel.radius_mm*np.cos(theta),np.zeros(len(theta)),vessel.radius_mm*np.sin(theta)))
    result=[ring]
    result.extend(segment[0]+np.linspace(0,1,int(np.ceil(8/step))+1)[:,None]*(segment[1]-segment[0]) for segment in vessel.branches.values())
    return result
