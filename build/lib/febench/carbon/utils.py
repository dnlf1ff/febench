from ase.io import write, read
import numpy as np
import sys

from ase.lattice.cubic import BodyCenteredCubic
from ase.build import make_supercell

def find_vac_idx(atoms, a, vac_pos, config, cont=False):
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
        
            atoms_index = find_vac_idx(atoms, a, vac_pos, config, cont=True)
    return atoms_index

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
        vac_idx = find_vac_idx(base, a, vac_pos, config)
        # Fe(n-q)Vac(q)
        del base[vac_idx]

        # Fe(n-q)C(p)Vac(q)
        base.append('C')
        base.positions[-1] = a * np.array(carbon_pos) 

        write(f'{struct_dir}/POSCAR_{label}', base, format='vasp')
        return

    if n_carbon == 1 and n_vac == 2:
        vac_idx_1 = find_vac_idx(base, a, vac_pos[0], config)
        del base[vac_idx_1]
        vac_idx_2 = find_vac_idx(base, a, vac_pos[1], config)
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
        vac_idx = find_vac_idx(base, a, vac_pos, config)
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

