from ase.io import write, read
import numpy as np
import sys
from ase.build import make_supercell

def find_nn_idx(atoms, nn_pos, a, config):
    x = nn_pos[0] * a
    y = nn_pos[1] * a
    z = nn_pos[2] * a

    pos = atoms.positions.copy()
    x_pos = pos[:,0]
    y_pos = pos[:,1]
    z_pos = pos[:,2]

    x_indices = np.where(x_pos==x)[0]
    y_indices = np.where(y_pos==y)[0]
    z_indices = np.where(z_pos==z)[0]

    index = set(x_indices) & set(y_indices) & set(z_indices)
    atoms_index = int(list(index)[0])
    return atoms_index

def write_poscar_from_config(config, solute, a):
    struct_dir = f'{config["tm"]["save"]}/structure'
    pos_dict = config['tm']['position']
    base_pos = pos_dict['base']

    # Fe(n-1)M
    base = read(f'{config["cwd"]}/POSCAR_base', format='vasp')
    base_idx = find_nn_idx(base, base_pos, a, config)
    del base[base_idx]
    base.append(solute)
    base.positions[-1] = np.array(base_pos) * a
    write(f'{struct_dir}/POSCAR_{solute}', base, format='vasp')

    nn1_pos = pos_dict['1nn']

    base_copy = base.copy()
    nn_idx = find_nn_idx(base_copy, nn1_pos, a, config)
    del base_copy[nn_idx]

    # Fe(n-2)M(2)
    base_copy.append(solute)
    base_copy.positions[-1] = np.array(nn1_pos) * a
    write(f'{struct_dir}/POSCAR_{solute}_{solute}_1nn', base_copy, format='vasp')

