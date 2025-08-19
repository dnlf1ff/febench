from ase.io import write, read
import numpy as np

def find_nn_idx(atoms, nn_pos, a):
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
    return int(list(index)[0])

def write_poscar_base(config, a):
    struct_dir = f'{config["pair"]["save"]}/structure'
    pos_dict = config['pair']['position']
    base_pos = pos_dict['base']

    for sol in config['pair']['solute'] + ['Vac']:
        base = read(f'{config["cwd"]}/POSCAR_base', format='vasp')
        base_idx = find_nn_idx(base, base_pos, a)
        del base[base_idx]
        if sol != 'Vac':
            base.append(sol)
            base.positions[-1] = np.array(base_pos) * a
        write(f'{struct_dir}/POSCAR_{sol}', base, format='vasp')


def write_poscar_from_config(config, a):
    struct_dir = f'{config["pair"]["save"]}/structure'
    pos_dict = config['pair']['position']
    nn1_pos = pos_dict['1nn']
    nn2_pos = pos_dict['2nn']
    nn3_pos = pos_dict['3nn']
    nn4_pos = pos_dict['4nn']
    nn5_pos = pos_dict['5nn']

    for sol1 in config['pair']['solute'] + ["Vac"]:
        base = read(f'{struct_dir}/POSCAR_{sol1}', format='vasp')

        for i, nn_pos in enumerate([nn1_pos, nn2_pos, nn3_pos, nn4_pos, nn5_pos]):
            for sol2 in config['pair']['solute'] + ["Vac"]:
                base_copy = base.copy()
                nn_idx = find_nn_idx(base_copy, nn_pos, a)
                del base_copy[nn_idx]

                if sol2 != 'Vac':
                    base_copy.append(sol2)
                    base_copy.positions[-1] = np.array(nn_pos) * a
                write(f'{struct_dir}/POSCAR_{sol1}_{sol2}_{int(i+1)}nn', base_copy, format='vasp')

