#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

band = np.loadtxt('kband.dat')
kptbd = np.loadtxt('kpt_bound.dat')

nbands = band.shape[1] - 4
xmin, xmax = band[:,0].min(), band[:,0].max()
x = band[:,0]

for bd in kptbd:
    ax.axvline(x=bd, color='r', lw=1.0, ls=':')

for i in range(nbands):
    plt.plot(x, band[:,i+4], '-', lw=1.0, color='k', ms=2,
            alpha=0.4)
    plt.plot(x, band[:,i+4], 'o', color='r', ms=2,
            alpha=0.4)


ax.axhline(y=0, color='b', ls='--', lw=1.0)
# ax.get_xaxis().set_visible(False)
ax.set_xlim((xmin, xmax))
ax.set_ylim((-2, 3))

pos = [0,] + list(kptbd) + [1,]
ax.set_xticks(pos)
print pos
ax.set_xticklabels(['M',r'$\Gamma$','X' ,'M' ], pos, 
                   fontsize=15)                   
ax.set_ylabel('Energy [eV]', fontsize=14)

for sp in ax.spines:
    ax.spines[sp].set_linewidth(1.5)

# ax.tick_params(axis='both', which='both',
#       top='off', bottom='off',
#       left='off', right='off',
#       labelleft='off', labelbottom='off')
# ax.minorticks_off()
# plt.show()
plt.savefig('kaka.png', dpi=360)
