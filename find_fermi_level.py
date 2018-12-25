#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import time


def find_fermi_level(band_energies, kpt_weight,
                     nelect, occ=2.0, sigma=0.01, nedos=100,
                     soc_band=False,
                     nmax=1000):
    '''
    '''

    if band_energies.ndim == 2:
        band_energies = band_energies[None, :]

    nspin, nkpts, nbnds = band_energies.shape

    if nspin == 1 and (not soc_band):
        occ = 2.0
    else:
        occ = 1.0

    if nbnds > nedos:
        nedos = nbnds * 5

    kpt_weight = np.asarray(kpt_weight, dtype=float)
    assert kpt_weight.shape == (nkpts,)
    kpt_weight /= np.sum(kpt_weight)

    emin = band_energies.min()
    emax = band_energies.max()
    e0 = np.linspace(emin, emax, nedos)
    de = e0[1] - e0[0]

    # find the approximated Fermi level
    nelect_lt_en = np.array([
        np.sum(occ * (band_energies <= en) * kpt_weight[None, :, None])
        for en in e0
    ])
    ne_tmp = nelect_lt_en[nedos/2]
    if (np.abs(ne_tmp - nelect) < 0.05):
        i_fermi = nedos / 2
        i_lower = i_fermi - 1
        i_upper = i_fermi + 1
    elif (ne_tmp > nelect):
        for ii in range(nedos/2-1, -1, -1):
            ne_tmp = nelect_lt_en[ii]
            if ne_tmp < nelect:
                i_fermi = ii
                i_lower = i_fermi
                i_upper = i_fermi + 1
                break
    else:
        for ii in range(nedos/2+1, nedos):
            ne_tmp = nelect_lt_en[ii]
            if ne_tmp > nelect:
                i_fermi = ii
                i_lower = i_fermi - 1
                i_upper = i_fermi
                break

    ############################################################
    # Below is the algorithm used by VASP, much slower
    ############################################################
    # find the approximated Fermi level
    # x = (e0[None, None, None, :] - band_energies[:, :, :, None]) / sigma
    # x = x.clip(-100, 100)
    # dos = 1./sigma * np.exp(x) / (np.exp(x) + 1)**2 * \
    #       kpt_weight[None, :, None, None] * de
    # ddos = np.sum(dos, axis=(0,1,2))
    #
    # nelect_from_dos_int = np.sum(ddos[:nedos/2])
    # if (np.abs(nelect_from_dos_int - nelect) < 0.05):
    #     i_fermi = nedos / 2 - 1
    #     i_lower = i_fermi - 1
    #     i_upper = i_fermi + 1
    # elif (nelect_from_dos_int > nelect):
    #     for ii in range(nedos/2, -1, -1):
    #         nelect_from_dos_int = np.sum(ddos[:ii])
    #         if nelect_from_dos_int < nelect:
    #             i_fermi = ii
    #             i_lower = i_fermi
    #             i_upper = i_fermi + 1
    #             break
    # else:
    #     for ii in range(nedos/2, nedos):
    #         nelect_from_dos_int = np.sum(ddos[:ii])
    #         if nelect_from_dos_int > nelect:
    #             i_fermi = ii
    #             i_lower = i_fermi - 1
    #             i_upper = i_fermi
    #             break

    # Locate the exact Fermi level using bisectioning
    e_lower = e0[i_lower]
    e_upper = e0[i_upper]
    lower_B = False
    upper_B = False
    for ii in range(nmax):
        e_fermi = (e_lower + e_upper) / 2.

        F_nk = occ / (np.exp(band_energies - e_fermi) + 1)
        N = np.sum(F_nk * kpt_weight[None, :, None])
        # print ii, e_lower, e_upper, N

        if (np.abs(N - nelect) < 1E-10):
            break
        if (np.abs(e_upper - e_lower / (np.abs(e_fermi) + 1E-10)) < 1E-14):
            raise ValueError("Cannot reach the specified precision!")

        if (N > nelect):
            if not lower_B: e_lower -= de
            upper_B = True
            e_upper = e_fermi
        else:
            if not upper_B: e_upper += de
            lower_B = True
            e_lower = e_fermi

    if (ii == nmax - 1):
        raise ValueError("Cannot reach the specified precision!")

    return e_fermi


if __name__ == "__main__":
    kw = np.array([
        0.01234568, 0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136,
        0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136,
        0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136,
        0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136,
        0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136,
        0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136,
        0.02469136, 0.02469136, 0.02469136, 0.02469136, 0.02469136
    ])

    e_nk = np.load('e_nk.npy')

    print find_fermi_level(e_nk, kw, nelect=54, sigma=0.01, soc_band=True)
