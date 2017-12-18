#!/usr/bin/env python

import re
import numpy as np
import matplotlib as mpl
mpl.use('agg')

import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftfreq, fftshift
from os.path import isfile

def extract_from_vasp_outcar (outFile = 'OUTCAR', whichK=1):
    # outFile = 'OUTCAR'

    OUTCAR = [line for line in open(outFile, 'r') if line.strip()]
    for line in OUTCAR:
        if 'NBANDS' in line and 'NKPTS' in line:
            NBANDS = int(line.split()[-1])
            NKPTS  = int(line.split()[ 3])
            break
    print 'Number of Bands: %d' % NBANDS
    print 'Number of Kpoints: %d' % NKPTS

    # where_scf_ends = [ii for ii, line in enumerate(OUTCAR)
    #                 if 'aborting loop because EDIFF' in line]
    # where_Efermi_starts = []

    where_scf_ends = []
    where_Efermi_starts = [ii for ii, line in enumerate(OUTCAR)
                           if 'E-fermi :' in line]
    # for ii in where_Efermi_starts:
    #     start = ii
    #     while 'E-fermi' not in OUTCAR[start]:
    #         start -= 1
    #     where_Efermi_starts += [start]

    Niters = len(where_Efermi_starts)
    # TDKSEN = np.zeres((Niters, NBANDS))
    TDKSEN = []
    for it, ii in enumerate(where_Efermi_starts):
        start = ii + 1
        end   = start + (NBANDS + 2) * NKPTS
        tmp = [ line.split()[1] for line in OUTCAR[start:end]
                if not re.search('[a-zA-Z]', line) ]
        TDKSEN += [tmp]
    TDKSEN = np.asarray(TDKSEN, dtype=float).reshape((-1, NKPTS, NBANDS))

    assert TDKSEN.shape[0] == Niters
    # np.savetxt('tden.dat', TDKSEN, fmt='%10.4f')

    return TDKSEN[:,whichK-1,:]

################################################################################
if isfile('tden.npy'):
    TDKS = np.load('tden.npy')
else:
    TDKS = extract_from_vasp_outcar()
    np.save('tden.npy', TDKS)

NSW   = TDKS.shape[0]
NBAND = TDKS.shape[1]
dt = 1.0
omega = 33.35640951981521 * 1E3 * fftfreq(NSW, dt)

vbmIndex = 324

TIME  = np.arange(NSW) * dt
TDKS -= np.average(TDKS[:, vbmIndex-1]) 
################################################################################

fig = plt.figure()
fig.set_size_inches((6.0, 4.0))

ax = plt.subplot(111)

for ii in range(NBAND):
    ax.plot(TIME, TDKS[:,ii], ls='-', lw=1.0, color='b', alpha=0.7)

# ax.set_xlim(TIME.max() - 2000, TIME.max())
ax.set_ylim(-1.0, 2.0)

plt.tight_layout()
plt.savefig('tdks.png', dpi=360)

from subprocess import call
call(['feh', '-xdF', 'tdks.png'])

# fig, axes = plt.subplots(nrows=2, ncols=1,
#                          sharex=False,
#                          sharey=False)
# plt.subplots_adjust(left=0.15, right=0.95,
#                     bottom=0.15, top=0.95,
#                     wspace=0.10, hspace=0.30)
# fig.set_size_inches((6, 4))
#
# Et = TDKS[:,324 - 1]
# power = np.abs(fft(Et - np.average(Et)))
#
# axes[0].plot(np.arange(NSW), TDKS[:,323], ls='-', lw=1.5, color='r')
# axes[0].plot(np.arange(NSW), TDKS[:,322], ls='-', lw=1.5, color='r')
# axes[0].plot(np.arange(NSW), TDKS[:,321], ls='-', lw=1.5, color='r')
# axes[1].plot(omega, power, ls='-', lw=1.5, color='k')
#
# axes[1].set_xlim(0, 1000)
# plt.show()
