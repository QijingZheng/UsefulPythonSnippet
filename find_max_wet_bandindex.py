#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import os

import matplotlib as mpl
mpl.use('agg')
mpl.rcParams['axes.unicode_minus'] = False

import matplotlib.pyplot as plt
############################################################
# time step in fs
dt   = 1.0
bmin = 109
bmax = 155

# randomly selected from the initial 1000 configs
rand_sele = np.arange(1, 1000 + 1, dtype=int)
np.random.shuffle(rand_sele)
rand_sele = rand_sele[:100]
inicon    = open('INICON', 'w')
############################################################
if os.path.isfile('all_wht.npy'):
    Wht = np.load('all_wht.npy')
    Enr = np.load('all_en.npy')

    nsw    = Enr.shape[0]
    nband  = Enr.shape[1]
    bndInd = np.arange(1, nband + 1, dtype=int)

    selected_states = []
    for tt in range(nsw):
        W = Wht[tt,bmin-1:bmax]
        E = Enr[tt,bmin-1:bmax]
        I = bndInd[bmin-1:bmax]

        sorted_index = np.argsort(W)
        selected_states.append([
            tt + 1,
            I[sorted_index[-1]],
            W[sorted_index[-1]]
        ])
        if (tt + 1) in rand_sele:
            inicon.write('{:5d} {:5d}\n'.format(tt+1, I[sorted_index[-1]]))

    with open('time_index_wht.dat', 'w') as out:
        for item in selected_states:
            out.write('{:5d} {:5d} {:8.4f}\n'.format(item[0], item[1], item[2]))

    ############################################################
    fig = plt.figure()
    fig.set_size_inches(4.8, 3.0)

    ########################################
    ax      = plt.subplot()
    nband   = Enr.shape[1]
    T, dump = np.mgrid[0:nsw:dt, 0:nband]
    sFac    = 8

    ax.scatter( T, Enr, s=Wht / Wht.max() * sFac, color='red', lw=0.0, zorder=1)
    for ib in range(nband):
        ax.plot(T[:,ib], Enr[:,ib], lw=0.5, color='k', alpha=0.5)

    newT = np.arange(nsw, dtype=int)
    markS = [xx[1] - 1 for xx in selected_states]
    ax.scatter(newT, Enr[newT, markS],
               color='b', lw=0.0,
               marker='^', s=1,
               zorder=2)

    ax.axhline(y=np.average(Enr[:,bmin-1]), lw=0.5, color='b', ls='--')
    ax.axhline(y=np.average(Enr[:,bmax-1]), lw=0.5, color='b', ls='--')

    ax.set_xlim(0, nsw)
    ax.set_ylim(0.0, 8.0)

    ax.set_xlabel('Time [fs]',   fontsize='small', labelpad=5)
    ax.set_ylabel('Energy [eV]', fontsize='small', labelpad=8)
    ax.tick_params(which='both', labelsize='x-small')

    ########################################
    plt.tight_layout(pad=0.2)
    plt.savefig('selected_state.png', dpi=360)
