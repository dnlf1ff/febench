from ase import Atoms
from ase.io import write, read
from ase.build import make_supercell 
from ase.lattice.cubic import BodyCenteredCubic
import yaml
import numpy as np
from febench.util.relax import aar_from_config
import pandas as pd
import sys

def find_vac_idx(atoms, a, vac_pos):
    x = vac_pos[0] * a
    y = vac_pos[1] * a
    z = vac_pos[2] * a

    pos = atoms.positions.copy()
    x_pos = pos[:,0]
    y_pos = pos[:,1]
    z_pos = pos[:,2]

    x_indices = np.where(x_pos==x)[0]
    y_indices = np.where(y_pos==y)[0]
    z_indices = np.where(z_pos==z)[0]

    index = set(x_indices) & set(y_indices) & set(z_indices)
    return int(list(index)[0])


def write_poscar_from_config(config, a, label, n_carbon, n_vac, carbon_pos, vac_pos):
    # Fe(n)C, Fe(n-q)Vac(q), Fe(n-q)C(p)Vac(q)


    struct_dir = f'{config["carbon"]["save"]}/structure'
    # Fe(n)
    base = read(f'{struct_dir}/POSCAR_base', format='vasp')

    # for Fe(n)C
    base_carbon = base.copy()
    print(n_carbon, n_vac)

    if n_carbon == 1 and n_vac == 1:
        print(f'1 carbon, 1 vacancy')
        vac_idx = find_vac_idx(base, a, vac_pos)
        # Fe(n-q)Vac(q)
        del base[vac_idx]
        write(f'{struct_dir}/POSCAR_{label}_vac', base)

        # Fe(n-q)C(p)Vac(q)
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos) 

        # Fe(n)C
        base_carbon.append('C')
        base_carbon.positions[-1] = a * np.array(carbon_pos) 
        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        write(f'{struct_dir}/POSCAR_{label}_C_1', base_carbon, format='vasp')
        return

    if n_carbon == 1 and n_vac == 2:
        print(f'1 carbon, 2 vacanceis')
        # Fe(n-q)Vac(q)
        vac_idx_1 = find_vac_idx(base, a, vac_pos[0])
        del base[vac_idx_1]
        vac_idx_2 = find_vac_idx(base, a, vac_pos[1])
        del base[vac_idx_2]
        write(f'{struct_dir}/POSCAR_{label}_vac', base)

        # Fe(n-q)C(p)Vac(q)
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos) 

        # Fe(n)C
        base_carbon.append('C')
        base_carbon.positions[-1] = a * np.array(carbon_pos) 
        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        write(f'{struct_dir}/POSCAR_{label}_C_1', base_carbon, format='vasp')
        return


    if n_carbon == 2 and n_vac == 0:
        print(f'2 carbons, 0 vacanceis')
        # two types of Fe(n)C ...
        base_carbon_2 = base_carbon.copy()

        carbon_pos_1 = carbon_pos[0]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_1) 

        # Fe(n)C type 1...
        base_carbon.append('C')
        base_carbon.positions[-1] = a * np.array(carbon_pos_1) 

        # Fe(n-q)C(p)Vac(q) q=0, p=2
        carbon_pos_2 = carbon_pos[1]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_2) 

        # Fe(n)C type 2...
        base_carbon_2.append('C')
        base_carbon_2.positions[-1] = a * np.array(carbon_pos_2) 
        write(f'{struct_dir}/POSCAR_{label}_C_2', base_carbon_2, format='vasp')
        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        write(f'{struct_dir}/POSCAR_{label}_C_1', base_carbon, format='vasp')
        return


    if n_carbon == 2 and n_vac == 1:
        print(f'2 carbons, 1 vacancy')
        # two types of Fe(n)C ...
        base_carbon_2 = base_carbon.copy()

        # Fe(n-q)Vac(q)
        vac_idx = find_vac_idx(base, a, vac_pos)
        del base[vac_idx]
        write(f'{struct_dir}/POSCAR_{label}_vac', base)

        carbon_pos_1 = carbon_pos[0]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_1) 

        # Fe(n)C type 1...
        base_carbon.append('C')
        base_carbon.positions[-1] = a * np.array(carbon_pos_1) 

        # Fe(n-q)C(p)Vac(q) q=1, p=2
        carbon_pos_2 = carbon_pos[1]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_2) 

        # Fe(n)C type 2...
        base_carbon_2.append('C')
        base_carbon_2.positions[-1] = a * np.array(carbon_pos_2) 
        write(f'{struct_dir}/POSCAR_{label}_C_2', base_carbon_2, format='vasp')
        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        write(f'{struct_dir}/POSCAR_{label}_C_1', base_carbon, format='vasp')
        return


    else:
        print(f'unknown carbon configuration')
        return

