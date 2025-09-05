from ase.build import make_supercell
from ase.lattice.cubic import BodyCenteredCubic
from ase.io import write
import numpy as np

from ase.units import J, m
EvAToJm = (m ** 2)/J

def write_fe_base(config, a):
    fe_unit = BodyCenteredCubic(directions=np.diag([1,1,1]), size=(1,1,1), symbol='Fe', pbc=True, latticeconstant=a)
    atoms = make_supercell(fe_unit, np.diag(config["carbon"]["supercell"]))
    write(f'{config["cwd"]}/POSCAR_base', atoms, format='vasp')


def write_csv(file, atoms, idx='pre', delimiter=',',conv=None):
    if conv is not None:
        try:
            conv=f"{atoms.info['opt_fa']},{atoms.info['opt_step']},{atoms.info['force_cnv']}"
        except:
            conv = '-,-,-'
    file.write(f"{idx}{delimiter}{atoms.info['e_fr_energy']}{delimiter}{delimiter}{len(atoms)}{delimiter}{atoms.info['a']}{delimiter}{atoms.info['b']}{delimiter}{atoms.info['c']}{delimiter}{atoms.info['alpha']}{delimiter}{atoms.info['beta']}{delimiter}{atoms.info['gamma']}{delimiter}{conv}\n")



