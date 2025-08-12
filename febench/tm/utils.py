from ase.io import write, read
import numpy as np

def find_nn_idx(atoms, nn_pos, a, config, cont=False):
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
    try:
        atoms_index = int(list(index)[0])

    except:
        if cont:
            sys.exit()
        else:
            atoms = BodyCenteredCubic(directions=np.diag([1,1,1]), size=(1,1,1),
                    symbol='Fe', pbc=True, latticeconstant=a)
            atoms = make_supercell(atoms,np.diag(config['carbon']['supercell']))
            write(f'{config["carbon"]["save"]}/POSCAR_base',atoms, format='vasp')
        
            atoms_index = find_vac_idx(atoms, a, nn_pos, config, cont=True)
    return atoms_index

def write_poscar_from_config(config, solute, a):
    struct_dir = f'{config["tm"]["save"]}/structure'
    pos_dict = config['tm']['position']
    base_pos = pos_dict['base']
    nn1_pos = pos_dict['1nn']
    nn2_pos = pos_dict['2nn']
    nn3_pos = pos_dict['3nn']
    nn4_pos = pos_dict['4nn']
    nn5_pos = pos_dict['5nn']

    # Fe(n-1)M
    base = read(f'{config["cwd"]}/POSCAR_base', format='vasp')
    base_idx = find_nn_idx(base, base_pos, a, config)
    del base[base_idx]
    base.append(solute)
    base.positions[-1] = np.array(base_pos) * a
    write(f'{struct_dir}/POSCAR_{solute}', base, format='vasp')

    for i, nn_pos in enumerate([nn1_pos, nn2_pos, nn3_pos, nn4_pos, nn5_pos]):
        base_copy = base.copy()
        nn_idx = find_nn_idx(base_copy, nn_pos, a, config)
        del base_copy[nn_idx]

        # Fe(n-2)MVac
        write(f'{struct_dir}/POSCAR_{solute}_Vac_{int(i+1)}nn', base_copy, format='vasp')

        # Fe(n-2)M(2)
        base_copy.append(solute)
        base_copy.positions[-1] = np.array(nn_pos) * a
        write(f'{struct_dir}/POSCAR_{solute}_{solute}_{int(i+1)}nn', base_copy, format='vasp')


