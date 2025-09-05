from ase.io import write, read
import numpy as np
import sys
from ase.build import make_supercell

def write_poscar_from_config(config, solute, a):
    struct_dir = f'{config["tm"]["save"]}/structure'

    # Fe(n-1)M
    base = read(f'{config["cwd"]}/POSCAR_base', format='vasp')
    base_idx = 32

    print(f'removing {base_idx}th Fe atom in {base.positions[base_idx]}')
    del base[base_idx]

    base.append(solute)
    base.positions[-1] = np.array(base_pos) * a
    write(f'{struct_dir}/POSCAR_{solute}', base, format='vasp')

    nn_idx = 1
    nn_pos = [0.5, 0.5, 0.5]
    base_copy = base.copy()
    print(f'removing {nn_idx}th Fe atom in {base.positions[nn_idx]}')
    del base_copy[nn_idx]

    # Fe(n-2)M(2)
    base_copy.append(solute)
    base_copy.positions[-1] = np.array(nn1_pos) * a
    write(f'{struct_dir}/POSCAR_{solute}_{solute}_1nn', base_copy, format='vasp')

