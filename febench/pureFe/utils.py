from ase.build import bcc100, bcc110, bcc111, bulk, make_supercell
from ase.lattice.cubic import BodyCenteredCubic
from ase.constraints import FixAtoms
from ase.io import write
import numpy as np
import sys

from febench.pureFe.matscipy_metrics import *

import pandas as pd

from ase.units import GPa, J, m
EvAToJm = (m ** 2)/J

def write_fe_base(config, a):
    fe_unit = BodyCenteredCubic(directions=np.diag([1,1,1]), size=(1,1,1), symbol='Fe', pbc=True, latticeconstant=a)
    atoms = make_supercell(fe_unit, np.diag(config["carbon"]["supercell"]))
    write(f'{config["cwd"]}/POSCAR_base', atoms, format='vasp')

def constrain_slab(slab):
    z = slab.positions[:,2].copy()
    index = [atom.index for atom in slab if atom.position[2] < np.percentile(z,62) and atom.position[2]>np.percentile(z,38)]
    fix_atom = FixAtoms(indices=index)
    slab.set_constraint(fix_atom)
    return slab


def get_slab(hkl, a0, vacuum=10, size=(4,4,20)):
    if hkl=='100':
        slab=bcc100('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    elif hkl=='110':
        slab=bcc110('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    elif hkl=='111':
        slab=bcc111('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    else:
        print('unknown surface miler index; check config.yaml')
        sys.exit(1)

    slab = constrain_slab(slab)
    return slab


def get_Cij(C_least_squares):
    df = pd.DataFrame(columns=['j1','j2','j3','j4','j5','j6'], index=['i1', 'i2','i3','i4','i5','i6'])
    for row in range(6):
        col_list = []
        for col in range(6):
            col_list.append(C_least_squares[row][col]/GPa)
        df.loc[f'i{int(row+1)}'] = col_list
    return df

def get_elastic_constants(C_least_squares):
    df = pd.DataFrame()
    df['bulk_voigt'] = get_voigt_bulk_modulus(C_least_squares)/GPa
    try:
        df['bulk_reuss'] = get_reuss_bulk_modulus(np.linalg.inv(C_least_squares/GPa))
        df['bulk_vrh'] = get_vrh_bulk_modulus(C_least_squares)/GPa
    except:
        df['bulk_reuss'] = np.nan
        df['bulk_vrh'] = np.nan
    df['shear_voigt'] = get_voigt_shear_modulus(C_least_squares)/GPa
    try:
        df['shear_reuss'] = get_reuss_shear_modulus(np.linalg.inv(C_least_squares/GPa))
        df['shear_vrh'] = get_vrh_shear_modulus(C_least_squares)/GPa
    except:
        df['shear_reuss']=np.nan
        df['shear_vrh']=np.nan
    return df 

def write_csv(file, atoms, idx='pre', delimiter=','):
    try:
        conv = atoms.info['conv']
    except:
        conv = '-'
    file.write(f"{idx}{delimiter}{atoms.info['e_fr_energy']}{delimiter}{atoms.info['volume']}{delimiter}{atoms.info['surface_area']}{delimiter}{len(atoms)}{delimiter}{atoms.info['a']}{delimiter}{atoms.info['b']}{delimiter}{atoms.info['c']}{delimiter}{atoms.info['alpha']}{delimiter}{atoms.info['beta']}{delimiter}{atoms.info['gamma']}{delimiter}{conv}\n")



