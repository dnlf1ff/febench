from ase.io import write, read
import pandas as pd
import numpy as np
import os.path as osp

# TODO:  csv file config
def find_idx(atoms, idx_pos, a):
    x = idx_pos[0] * a
    y = idx_pos[1] * a
    z = idx_pos[2] * a

    pos = atoms.positions.copy()
    x_pos = pos[:,0]
    y_pos = pos[:,1]
    z_pos = pos[:,2]

    x_indices = np.where(x_pos==x)[0]
    y_indices = np.where(y_pos==y)[0]
    z_indices = np.where(z_pos==z)[0]

    index = set(x_indices) & set(y_indices) & set(z_indices)
    return int(list(index)[0])

def write_atoms(triplet, a, base, idx, pos, sol_type, sol, struct_dir):
    filename=f'{struct_dir}/POSCAR_{triplet}_{sol_type}_{sol}'
    if osp.isfile(filename):
        return
    base_atoms = base.copy()
    del base_atoms[idx]
    if sol == 'Vac':
        write(filename, base_atoms)
        return
    else:
        base_atoms.append(sol)
        base_atoms.positions[-1] = np.array(pos) * a
        write(filename, base_atoms)
        return

def write_poscar_from_config(config, a, triplet):
    df_ref = pd.read_csv(f'triplet_{triplet}.csv')
    struct_dir = f'{config["triplet"]["save"]}/structure'
    base = read(f'{config["cwd"]}/POSCAR_base', format='vasp')

    pos_dict = config['triplet']['position'][triplet]
    A_pos = pos_dict[0]
    A_idx = find_idx(base, A_pos, a)
    B_pos = pos_dict[1]
    B_idx = find_idx(base, B_pos, a)
    C_pos = pos_dict[2]
    C_idx = find_idx(base, C_pos, a)

    A_sols = df_ref['A']
    B_sols = df_ref['B']
    C_sols = df_ref['C']

    for A, B, C in zip(A_sols, B_sols, C_sols):
        A_base = base.copy()
        B_base = base.copy()
        C_base = base.copy()

        write_atoms(triplet, a, A_base, A_idx, A_pos, sol_type='A', sol=A, struct_dir=struct_dir)
        write_atoms(triplet, a, B_base, B_idx, B_pos, sol_type='B', sol=B, struct_dir=struct_dir)
        write_atoms(triplet, a, C_base, C_idx, C_pos, sol_type='C', sol=C, struct_dir=struct_dir)

        tot_base = base.copy()

        A_idx_tot = find_idx(tot_base, A_pos, a)
        del tot_base[A_idx_tot]
        if A != 'Vac':
            tot_base.append(A)
            tot_base.positions[-1] = np.array(A_pos) * a

        B_idx_tot = find_idx(tot_base, B_pos, a)
        del tot_base[B_idx_tot]
        if B != 'Vac':
            tot_base.append(B)
            tot_base.positions[-1] = np.array(B_pos) * a

        C_idx_tot = find_idx(tot_base, C_pos, a)
        del tot_base[C_idx_tot]
        if C != 'Vac':
            tot_base.append(C)
            tot_base.positions[-1] = np.array(C_pos) * a

        write(f'{struct_dir}/POSCAR_{triplet}_ABC_{A}{B}{C}', tot_base)
