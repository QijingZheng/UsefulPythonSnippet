#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import Voronoi

from ase.io import read

icell = np.linalg.inv(
    read('POSCAR', format='vasp').cell
).T
# icell = np.linalg.inv([[5, 5, 0], [0, 5, 5], [5, 0, 5]]).T
Grid = np.tensordot(
    icell,
    np.mgrid[-1:2, -1:2, -1:2],
    axes=[0, 0]
).reshape((3, -1)).T

vor = Voronoi(Grid)
bz_ridges = []
for r in vor.ridge_dict:
    if r[0] == 13 or r[1] == 13:
        bz_ridges.append([vor.vertices[i] for i in vor.ridge_dict[r]])

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for ridge in bz_ridges:
    x, y, z = np.array(ridge).T
    ax.plot(x, y, z, color='k')

plt.show()
