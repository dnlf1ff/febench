from ase.build import bcc100, bcc110, bcc111, bulk, make_supercell
from ase.lattice.cubic import BodyCenteredCubic
from ase.constraints import FixAtoms
from ase.io import write
import numpy as np
import sys

from ase.units import J, m
EvAToJm = (m ** 2)/J

def write_fe_base(config, a):
    fe_unit = BodyCenteredCubic(directions=np.diag([1,1,1]), size=(1,1,1), symbol='Fe', pbc=True, latticeconstant=a)
    atoms = make_supercell(fe_unit, np.diag(config["carbon"]["supercell"]))
    write(f'{config["cwd"]}/POSCAR_base', atoms, format='vasp')

def get_slab(hkl, a0, vacuum=20, size=(4,4,20)):
    if hkl=='100':
        slab=bcc100('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    elif hkl=='110':
        slab=bcc110('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    elif hkl=='111':
        slab=bcc111('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    else:
        print('unknown surface miler index; check config.yaml')
        sys.exit(1)

    return slab


def write_csv(file, atoms, idx='pre', delimiter=','):
    try:
        conv = atoms.info['conv']
    except:
        conv = '-'
    file.write(f"{idx}{delimiter}{atoms.info['e_fr_energy']}{delimiter}{atoms.info['surface_area']}{delimiter}{len(atoms)}{delimiter}{atoms.info['a']}{delimiter}{atoms.info['b']}{delimiter}{atoms.info['c']}{delimiter}{atoms.info['alpha']}{delimiter}{atoms.info['beta']}{delimiter}{atoms.info['gamma']}{delimiter}{conv}\n")



