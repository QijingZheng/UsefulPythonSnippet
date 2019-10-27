#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, argparse
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms, FixScaled, FixedPlane, FixedLine

def parse_cml_args(cml):
    '''
    CML parser.
    '''
    arg = argparse.ArgumentParser(add_help=True)

    arg.add_argument('-i', dest='poscar', action='store', type=str,
                     default='POSCAR',
                     help='The POSCAR based on which to make supercell.')
    arg.add_argument('-o', dest='out', action='store', type=str,
                     default='out.vasp',
                     help='Default supercell filename.')
    arg.add_argument('-n', dest='new_sym_order', action='store',
                     type=str, default=None, nargs='*',
                     help='New order of the chemical symbols.')
    arg.add_argument('-x', dest='xyz', action='store', type=str,
                     default='z', choices=['x', 'y', 'z'],
                     help='Sort with increasing x/y/z-coordinates.')
    arg.add_argument('-s', '--size', dest='size', action='store', type=int,
                     default=[1, 1, 1], nargs=3,
                     help='The supercell size.')

    return arg.parse_args(cml)

def mk_supercell(cml):
    arg = parse_cml_args(cml)

    pc = read(arg.poscar)
    sc = pc * arg.size
    org_atom_index = np.arange(len(sc), dtype=int)

    # First sort accoring to the x/y/z-coordinates
    if arg.xyz:
        new_atom_index = np.argsort(
                np.round(sc.positions, 4)[:, 'xyz'.index(arg.xyz)]
                )
        sc = sc[new_atom_index]

    # New order of chemical symbols
    if arg.new_sym_order:
        assert set(arg.new_sym_order) == set(org_chem_symbols)
        chem_sym_order = arg.new_sym_order
    else:
        # Just stick to the alphabetical order
        chem_sym_order = sorted(set(org_chem_symbols))

    # Re-arrange the atoms according to the new order of chemical symbols
    new_atom_index = [ii for ss in chem_sym_order
                         for ii in org_atom_index[sc.symbols == ss]]
    sc = sc[new_atom_index]

    write(arg.out, sc, vasp5=True, direct=True)

if __name__ == '__main__':
    mk_supercell(sys.argv[1:])
