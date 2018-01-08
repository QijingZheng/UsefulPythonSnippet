#!/usr/bin/env python

import re, os, sys
import numpy as np

################################################################################    
def WeightFromPro(infile='PROCAR', lsorbit=False, use_last_column=True):
    """
    Contribution of selected atoms to the each KS orbital
    """

    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub('[^0-9]', ' ', FileContents[1]).split()]

    # Weights = np.asarray([line.split()[-1] for line in FileContents
    #                       if not re.search('[a-zA-Z]', line)], dtype=float)
    Weights = np.asarray([line.split()[1:] for line in FileContents
                          if not re.search('[a-zA-Z]', line)], dtype=float)

    kpt_weight = np.asarray([line.split()[-1] for line in FileContents if 'weight' in line], dtype=float)
    
    energies = np.asarray([line.split()[-4] for line in FileContents
                            if 'occ.' in line], dtype=float)

    nlmax = Weights.shape[-1]
    nspin = Weights.shape[0] / (nkpts * nbands * nions)
    nspin /= 4 if lsorbit else 1

    if lsorbit:
        Weights.resize(nspin, nkpts, nbands, 4, nions, nlmax)
        Weights = Weights[:,:,:,0,:,:]
    else:
        Weights.resize(nspin, nkpts, nbands, nions, nlmax)

    if use_last_column:
        Weights = Weights[..., -1]
    else:
        Weights = Weights[..., :-1]

    kpt_weight.resize(nspin, nkpts)
    energies.resize(nspin, nkpts, nbands)
    
    return energies, kpt_weight, Weights

############################################################    
# select the spin index, starting from 0
whichSpin   = 0
# select the k-point index, starting from 0
whichKpoint = 0
# select the atom index, starting from 0
m = np.arange(54)
w = np.arange(54) + 54

if __name__ == '__main__':
    rundir = sys.argv[1]
    proF   = rundir + '/PROCAR'
    outF   = rundir + '/enw.dat'

    en, kptw, wht = WeightFromPro(proF)
    # select the spin and the k-point
    wht = wht[whichSpin, whichKpoint]

    # sum over atoms
    wht_m = np.sum(wht[:, m], axis=-1)
    wht_w = np.sum(wht[:, w], axis=-1)

    np.savetxt(outF, np.vstack([en[whichSpin, whichKpoint, :], wht_m, wht_w]), fmt='%8.4f')
