#!/usr/bin/env python

import numpy as np

def load_vibmodes_from_outcar(inf='OUTCAR', include_imag=False):
    '''
    Read vibration eigenvectors and eigenvalues from OUTCAR.
    '''

    out = [line for line in open(inf) if line.strip()]
    ln  = len(out)
    for line in out:
        if "NIONS =" in line:
            nions = int(line.split()[-1])
            break

    THz_index = []
    for ii in range(ln-1,0,-1):
        if '2PiTHz' in out[ii]:
            THz_index.append(ii)
        if 'Eigenvectors and eigenvalues' in out[ii]:
            i_index = ii + 2
            break
    j_index = THz_index[0] + nions + 2

    real_freq = [False if 'f/i' in line else True
                 for line in out[i_index:j_index] 
                 if '2PiTHz' in line]

    omegas = [line.split()[ -4] for line in out[i_index:j_index]
              if '2PiTHz' in line]
    modes  = [line.split()[3:6] for line in out[i_index:j_index]
              if ('dx' not in line) and ('2PiTHz' not in line)]

    omegas = np.array(omegas)
    modes  = np.array(modes).reshape((-1, nions, 3))

    if include_imag:
        omegas = omegas[real_freq]
        modes  = modes[real_freq]

    return omegas, modes

