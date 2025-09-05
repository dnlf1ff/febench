from ase.io import write, read
import numpy as np
import sys

def find_vac_idx(first=True):
    if first:
        return 1
    else:
        return 8

def write_FeC_poscar(config, a):
    carbon_pos = [0.5, 0.5, 0]
    struct_dir = f'{config["carbon"]["save"]}/structure'
    # Fe(n)
    base = read(f'{config["cwd"]}/POSCAR_base', format='vasp')
    # Fe(n)C
    base.append('C')
    base.positions[-1] = a * np.array(carbon_pos) 
    write(f'{struct_dir}/POSCAR_C', base, format='vasp')
    return

def write_poscar_from_config(config, a, label, n_carbon, n_vac, carbon_pos, vac_pos):
    # Fe(n)C, Fe(n-q)Vac(q), Fe(n-q)C(p)Vac(q)

    struct_dir = f'{config["carbon"]["save"]}/structure'
    base_atoms = read(f'{config["cwd"]}/POSCAR_base', format='vasp')
    base = base_atoms.copy()

    if n_carbon == 1 and n_vac == 1:
        vac_idx = find_vac_idx()
        # Fe(n-q)Vac(q)
        print(f'removing {vac_idx}th Fe atom in {base.positions[vac_idx]}')
        del base[vac_idx]

        # Fe(n-q)C(p)Vac(q)
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos) 

        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        return

    if n_carbon == 1 and n_vac == 2:
        vac_idx_1 = find_vac_idx()
        print(f'removing {vac_idx_1}th Fe atom in {base.positions[vac_idx_1]}')
        del base[vac_idx_1]
        vac_idx_2 = find_vac_idx(first=False)
        print(f'removing {vac_idx_2}th Fe atom in {base.positions[vac_idx_2]}')
        del base[vac_idx_2]

        # Fe(n-q)C(p)Vac(q)
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos) 

        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        return


    if n_carbon == 2 and n_vac == 0:
        carbon_pos_1 = carbon_pos[0]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_1) 

        # Fe(n-q)C(p)Vac(q) q=0, p=2
        carbon_pos_2 = carbon_pos[1]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_2) 

        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        return


    if n_carbon == 2 and n_vac == 1:
        # Fe(n-q)Vac(q)
        vac_idx = find_vac_idx()
        del base[vac_idx]

        carbon_pos_1 = carbon_pos[0]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_1) 

        # Fe(n-q)C(p)Vac(q) q=1, p=2
        carbon_pos_2 = carbon_pos[1]
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos_2) 

        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        return

    else:
        raise NotImplementedError

