#!/usr/bin/env python

from ase.io import read

xx = read('POSCAR', format='vasp')

zmax = xx.positions[:,-1].max()
zmin = xx.positions[:,-1].min()

print (zmax - zmin) / (2. * xx.cell[-1,-1])
