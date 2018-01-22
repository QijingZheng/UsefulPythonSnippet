#!/usr/bin/env python
import os

################################################################################
# VDW-D2 method of Grimme, parameters extracted from vdwforcefield.F in VASP
# 5.4.1.
vdwd2_elements = [
    'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', 'Na', 'Mg', 'Al',
    'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe',
    'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ',
    'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
    'I ', 'Xe', 'X ', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
    'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir',
    'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn'
]

vdw_c6 = [
    0.14, 0.08, 1.61, 1.61, 3.13, 1.75, 1.23, 0.70, 0.75, 0.63, 5.71, 5.71, 10.79,
    9.23, 7.84, 5.57, 5.07, 4.61, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8, 10.8,
    10.8, 10.8, 10.8, 10.8, 16.99, 17.10, 16.37, 12.64, 12.47, 12.01, 24.67, 24.67,
    24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 24.67, 37.32,
    38.71, 38.44, 31.74, 31.50, 29.99, 29.99, 315.275, 226.994, 176.252, 140.68,
    140.68, 140.68, 140.68, 140.68, 140.68, 140.68, 140.68, 140.68, 140.68, 140.68,
    140.68, 140.68, 140.68, 105.112, 81.24, 81.24, 81.24, 81.24, 81.24, 81.24,
    81.24, 57.364, 57.254, 63.162, 63.540, 55.283, 57.171, 56.64
]

vdw_r0 = [
    1.001, 1.012, 0.825, 1.408, 1.485, 1.452, 1.397, 1.342, 1.287, 1.243, 1.144,
    1.364, 1.639, 1.716, 1.705, 1.683, 1.639, 1.595, 1.485, 1.474, 1.562, 1.562,
    1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.650, 1.727, 1.760,
    1.771, 1.749, 1.727, 1.628, 1.606, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639,
    1.639, 1.639, 1.639, 1.639, 1.672, 1.804, 1.881, 1.892, 1.892, 1.881, 1.881,
    1.802, 1.762, 1.720, 1.753, 1.753, 1.753, 1.753, 1.753, 1.753, 1.753, 1.753,
    1.753, 1.753, 1.753, 1.753, 1.753, 1.753, 1.788, 1.772, 1.772, 1.772, 1.772,
    1.772, 1.772, 1.772, 1.758, 1.989, 1.944, 1.898, 2.005, 1.991, 1.924
]
################################################################################

def get_elements_from_poscar(inp='POSCAR'):
    if os.path.isfile(inp):
        pos = open(inp).readlines()
        if len(pos) < 6:
            return None
        else:
            ele_list       = pos[5].split()
            alpha_or_digit = [xx.isalpha() for xx in ele_list]

            if all(alpha_or_digit):
                return ele_list
            else:
                return None

def get_elements_from_potcar(inp='POTCAR'):
    if os.path.isfile(inp):
        ele_list = [line.split('=')[1].split(':')[0].strip() for line in open(inp) 
                    if 'VRHFIN' in line]

        if ele_list:
            return ele_list
        else:
            return None

if __name__ == '__main__':
    posL = get_elements_from_poscar('POSCAR')
    potL = get_elements_from_potcar('POTCAR')

    if posL and potL:
        assert len(posL) <= len(potL) and all([n1.upper() == n2.upper() for n1, n2 in zip(posL, potL)]), \
               "\nElements in POSCAR and POTCAR does NOT match!\n  POSCAR:{}\n  POTCAR:{}".format(
                       ' '.join(['%3s' % xx for xx in posL]), ' '.join(['%3s' % xx for xx in potL])
                       )

    ############################################################
    L = posL if posL else potL
    if L:
        E      = [xx.strip().upper()  for xx in vdwd2_elements]
        EINDEX = [E.index(xx.upper()) for xx in L]
        c6     = [vdw_c6[ii] for ii in EINDEX]
        r0     = [vdw_r0[ii] for ii in EINDEX]

        ############################################################
        print "LVDW = .TRUE."
        print "# %s"        % ' '.join(["%2s"   % xx for xx in L ])
        print "VDW_C6 = %s" % ' '.join(["%7.3f" % xx for xx in c6])
        print "VDW_R0 = %s" % ' '.join(["%7.3f" % xx for xx in r0])
        ############################################################
