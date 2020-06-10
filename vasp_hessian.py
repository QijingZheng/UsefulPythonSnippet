#!/usr/bin/env python

import numpy as np
import xml.etree.ElementTree as ET


def read_hessian_vasprun(inf='vasprun.xml', dynm=True):
    '''
    Read the dynamical matrix from vasprun.xml, multiply the masses if
    necessary.
    '''

    atom_types = []
    masses = []
    num_atoms = []
    sel = None

    for ev, elem in ET.iterparse(inf):
        if elem.tag == 'array':
            if 'name' in elem.attrib:
                if elem.attrib['name'] == 'atomtypes':
                    for rc in elem.findall('./set/rc'):
                        atom_info = [x.text for x in rc.findall('./c')]
                        num_atoms.append(int(atom_info[0]))
                        atom_types.append(atom_info[1].strip())
                        masses.append(float(atom_info[2]))

        if elem.tag == 'varray':
            if elem.attrib['name'] == 'hessian':
                # the dynamical matrix, divided by SQRT(m_i * m_j)
                fc = np.array(
                    [
                        np.fromstring(v.text, sep=' ')
                        for v in elem.findall('./v')
                    ]
                )
            if elem.attrib['name'] == 'selective':
                sel = np.array(
                    [
                        True if x == "T" else False
                        for v in elem.findall('./v')
                        for x in v.text.strip().split()
                    ]
                )

    if dynm:
        return num_atoms, masses, sel, fc
    else:
        # the force matrix, multiply by SQRT(m_i * m_j)
        M = np.concatenate([3 * n * [m] for n, m in zip(num_atoms, masses)])
        if sel is not None:
            M = M[sel]
        Mi, Mj = np.meshgrid(M, M)
        SQRT_M = np.sqrt(Mi * Mj)
        return num_atoms, masses, sel, fc * SQRT_M


if __name__ == '__main__':
    mass_silver    = 107.868
    mass_deuterium = 2.014
    mass_hydrogen  = 1.000
    mass_carbon    = 12.011
    new_masses     = [mass_silver, mass_carbon, mass_deuterium]
    # new_masses     = [mass_silver, mass_carbon, mass_hydrogen]

    EVTOJ          = 1.60217733E-19
    AMTOKG         = 1.6605402E-27
    PLANK          = 6.626075E-34
    factor         = -np.sqrt(EVTOJ * 1E20 / AMTOKG)

    omegas_nelects = []
    modes_nelects  = []
    nelects  = range(5)
    for ne in nelects:
        vasprun = 'vib_c{:d}_2/vasprun.xml'.format(ne)
        print(vasprun)
        num_atoms, matoms, sel, dynm = read_hessian_vasprun(vasprun, dynm=False)
        natoms = np.sum(num_atoms)
        # ndof = natoms * 3 if sel is None else sel.sum()

        M = np.concatenate([3 * n * [m] for n, m in zip(num_atoms, new_masses)])
        if sel is not None:
            M = M[sel]

        Mi, Mj = np.meshgrid(M, M)
        SQRT_M = np.sqrt(Mi * Mj)
        dynm /= SQRT_M

        e, v = np.linalg.eigh(dynm)
        isort = np.argsort(e)
        e = e[isort]
        if sel is None:
            v_ext = v[:,isort].T.reshape((-1, natoms, 3))
        else:
            neig = len(e)
            v_ext = np.zeros((neig, natoms*3))
            v_ext[:,sel] = v[:,isort].T
            v_ext = v_ext.reshape((neig, natoms, 3))

        omegas_nelects.append(np.sign(e) * np.sqrt(np.abs(e)) * factor * 1000 *
                              PLANK / EVTOJ/ 2 / np.pi)
        modes_nelects.append(v_ext)

    np.save('all_omegas_deu.npy', omegas_nelects)
    np.save('all_modes_deu.npy',  modes_nelects)

    # np.save('all_omegas.npy', omegas_nelects)
    # np.save('all_modes.npy',  modes_nelects)
