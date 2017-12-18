#!/usr/bin/env python

import numpy as np
from ase.io import read, write
from ase import Atoms

def posFromXdatcar(inf='XDATCAR', direct=True, return_vel=False, dt=1.0):
    '''
    extract coordinates from XDATCAR and create POSCAR files.
    
    Input arguments:
        inf:        location of XDATCAR
        direct:     coordinates in fractional or cartesian 
        return_vel: whether to calculate velocity or not
        dt:         time step of MD
    '''

    inp = open(inf).readlines()

    ElementNames = inp[5].split()
    ElementNumbers = np.array([int(x) for x in inp[6].split()], dtype=int)
    cell = np.array([line.split() for line in inp[2:5]], dtype=float)
    Natoms = ElementNumbers.sum()

    # construct ASE atoms object without positions.
    ChemFormula = ''.join(['%s%d' % (xx, nn) for xx, nn in zip(ElementNames, ElementNumbers)])
    geo = Atoms(ChemFormula, positions=np.zeros((Natoms, 3)), cell=cell, pbc=True)

    # the coordinates of atoms at each time step
    positions = [line.split() for line in inp[8:] if line.strip() and 'config'
                 not in line]
    positions = np.array(positions, dtype=float).reshape((-1, Natoms, 3))
    if not direct:
        pos = np.dot(positions, cell)
    else:
        pos = positions

    # the velocity of atoms at each time step
    if return_vel:
        # a simple way to estimate the velocities from positons.
        vel = np.diff(positions, axis=0) / dt
        # periodic boundary condition.
        vel[vel > 0.5] -= 1.0
        vel[vel <-0.5] += 1.0
        vel = np.dot(vel, cell)
    else:
        vel = None

    return geo, pos, vel


############################################################
geo, positions, velocities = posFromXdatcar('XDATCAR', direct=False)
niter = positions.shape[0]
natom = positions.shape[1]

# extract the last NSW configurations from the trajectory
NSW=2000
start = niter - NSW
end = niter
for ii in range(start, end):
    outF = 'POSCAR_%04d' % (ii + 1 - start)

    geo.set_scaled_positions(positions[ii,...])
    # geo = Atoms(ChemFormula, scaled_positions=positions[ii,...], cell=cell,
    #         pbc=True)

    write(outF, geo, format='vasp', vasp5=True, direct=True)
