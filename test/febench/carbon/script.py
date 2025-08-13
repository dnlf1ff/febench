from ase import Atoms
from ase.io import write, read
from ase.build import make_supercell 
from ase.lattice.cubic import BodyCenteredCubic
import yaml
import gc
import numpy as np
from febench.util.parse_args import parse_base_args
from febench.util.utils import dumpYAML 
from febench.util.parse_config import parse_config_yaml 
from febench.util.relax import aar_from_config
from febench.carbon.utils import write_poscar_from_config

import pandas as pd

def fe_base(config):
    bulk_df = pd.read_csv(f'{config["pureFe"]["save"]}/bulk.csv')
    a = bulk_df['a'][1]/config["pureFe"]["bulk"]["supercell"][0]
    struct_dir = f'{config["carbon"]["save"]}/structure'

    fe_unit = BodyCenteredCubic(directions=np.diag([1,1,1]), size=(1,1,1), symbol='Fe', pbc=True, latticeconstant=a)
    atoms = make_supercell(fe_unit, np.diag(config["carbon"]["supercell"]))
    write(f'{struct_dir}/POSCAR_base', atoms, format='vasp')

    return a

def process_carbon(config, calc):
    a = fe_base(config)

    save_dir = config["carbon"]["save"]
    struct_dir = f'{config["carbon"]["save"]}/structure'
    log_dir = f'{config["carbon"]["save"]}/log'

    labels = config["carbon"]["label"]
    
    with open(config["carbon"]["config"], 'r') as f:
        carbon_config = yaml.load(f, Loader=yaml.FullLoader)
    df = pd.DataFrame(columns = ['E_FeCVac', 'n_atom', 'n_carbon', 'n_vacancy', 'E_FeC_1', 'E_FeC_2', 'E_FeVac'], index=labels)

    for label in labels:
        print(label)
        carbon_args = {
                'a': a,
                'label': label,
                'n_carbon': carbon_config[label]['n_carbon'],
                'n_vac': carbon_config[label]['n_vacancy'],
                'carbon_pos': carbon_config[label]['carbon'],
                'vac_pos': carbon_config[label]['vacancy'],
                }

        print(carbon_args)


        write_poscar_from_config(config, **carbon_args)

        # calc #Fe(n-q)C(p)Vac(q) 
        atoms = read(f'{struct_dir}/POSCAR_{label}', format='vasp')

        ase_atom_relaxer = aar_from_config(config, calc,opt=config["carbon"]["opt"], logfile = f'{log_dir}/{label}_FeCVac.log')
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.info['conv'] = conv
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{label}', atoms, format='vasp')

        del  ase_atom_relaxer
        gc.collect()

        # calc Fe(n-q)Vac(q)
        try:
            atoms_vac = read(f'{struct_dir}/POSCAR_{label}_vac', format='vasp')

            ase_atom_relaxer = aar_from_config(config, calc,opt=config["carbon"]["opt"], logfile = f'{log_dir}/{label}_FeVac.log')
            atoms_vac, conv = ase_atom_relaxer.relax_atoms(atoms_vac)
            atoms_vac = ase_atom_relaxer.update_atoms(atoms_vac)
            atoms_vac.info['conv'] = conv
            atoms_vac.calc = None
            write(f'{struct_dir}/CONTCAR_{label}_vac', atoms_vac, format='vasp')

            del  ase_atom_relaxer
            gc.collect()

        except:
            # config w/o vacancy
            atoms_vac = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')

 
        # calc Fe(n)C(p)
        atoms_C_1 = read(f'{struct_dir}/POSCAR_{label}_C_1', format='vasp')

        ase_atom_relaxer = aar_from_config(config, calc,opt=config["carbon"]["opt"], logfile = f'{log_dir}/{label}_FeC_1.log')
        atoms_C_1, conv = ase_atom_relaxer.relax_atoms(atoms_C_1)
        atoms_C_1 = ase_atom_relaxer.update_atoms(atoms_C_1)
        atoms_C_1.info['conv'] = conv
        atoms_C_1.calc = None
        write(f'{struct_dir}/CONTCAR_{label}_C_1', atoms_C_1, format='vasp')
        del  ase_atom_relaxer
        gc.collect()

        try:
            # there are two types of Fe(n)C for configurations with two Carbon interstitals
            atoms_C_2 = read(f'{struct_dir}/POSCAR_{label}_C_2', format='vasp')

            ase_atom_relaxer = aar_from_config(config, calc,opt=config["carbon"]["opt"], logfile = f'{log_dir}/{label}_FeC_2.log')
            atoms_C_2, conv = ase_atom_relaxer.relax_atoms(atoms_C_2)
            atoms_C_2 = ase_atom_relaxer.update_atoms(atoms_C_2)
            atoms_C_2.info['conv'] = conv
            atoms_C_2.calc = None
            write(f'{struct_dir}/CONTCAR_{label}_C_2', atoms_C_2, format='vasp')

        except Exception as e:
            # in the case p=1 at Fe(n-q)C(p)Vac(q)
            print(f'Exception {e} occured while attempting to calculate Fe(n)C config 2 for label {label}')
            atoms_C_2 = atoms_C_1

        df.loc[label] = [atoms.info['e_fr_energy'], len(atoms), carbon_args['n_carbon'],carbon_args['n_vac'], atoms_C_1.info['e_fr_energy'], atoms_C_2.info['e_fr_energy'], atoms_vac.info['e_fr_energy']]
        df.to_csv(f'{save_dir}/carbon.csv', index_label = 'config')

        del atoms, atoms_vac, atoms_C_1, atoms_C_2


def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    args = parse_base_args(argv)
    config_dir = args.config
    calc = args.calc
    modal = args.modal

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config['calculator']['prefix'] = calc
    config['calculator']['modal'] = modal
    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/config_carbon.yaml')
    calc = calc_from_config(config)

    process_carbon(config, calc)



if __name__ == '__main__':
    main()
