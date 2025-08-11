from ase import Atoms
from ase.io import write, read
import yaml
import gc
import numpy as np
from febench.util.parse_args import parse_base_args
from febench.util.utils import dumpYAML 
from febench.util.parse_config import parse_config_yaml 
from febench.util.relax import aar_from_config
from febench.carbon.utils import write_poscar_from_config, write_FeC_poscar

import torch
from tqdm import tqdm
import warnings

def process_carbon(config, calc):
    save_dir = config["carbon"]["save"]
    struct_dir = f'{config["carbon"]["save"]}/structure'
    log_dir = f'{config["carbon"]["save"]}/log'
    
    atoms_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    a = atoms_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]
    E_Fe = atoms_bulk.info['e_fr_energy']
    n_Fe = len(atoms_bulk)

    atoms_Vac = read(f'{config["pureFe"]["save"]}/structure/Vac_opt.extxyz')
    E_FeVac = atoms_Vac.info['e_fr_energy']
    n_FeVac = len(atoms_Vac)

    write_FeC_poscar(config, a)
    
    with open(config["carbon"]["config"], 'r') as f:
        carbon_config = yaml.load(f, Loader=yaml.FullLoader)

    csv_file = open(f'{save_dir}/carbon.csv', 'w', buffering = 1)
    csv_file.write('config,E_bind,E_FeVac,n_FeVac,E_FeC,n_FeC,E_Fe,n_Fe,E_FeCVac,n_FeCVac,n_carbon,n_vacancy\n')

    atoms_C = read(f'{struct_dir}/POSCAR_C', format='vasp')
    ase_atom_relaxer = aar_from_config(config, calc,opt=config["carbon"]["opt"], logfile = f'{log_dir}/FeC.log')
    atoms_C, conv = ase_atom_relaxer.relax_atoms(atoms_C)
    atoms_C = ase_atom_relaxer.update_atoms(atoms_C)
    atoms_C.info['conv'] = conv
    atoms_C.calc = None
    write(f'{struct_dir}/CONTCAR_C', atoms_C, format='vasp')
    write(f'{struct_dir}/FeC.extxyz', atoms_C, format='extxyz')

    E_FeC = atoms_C.info['e_fr_energy']
    n_FeC = len(atoms_C)

    del atoms_bulk, atoms_Vac, atoms_C, ase_atom_relaxer
    gc.collect()

    labels = config["carbon"]["label"]
    for idx, label in enumerate(tqdm(labels, desc='processing carbon in iron ...')):
        carbon_args = {
                'a': a,
                'label': label,
                'n_carbon': carbon_config[label]['n_carbon'],
                'n_vac': carbon_config[label]['n_vacancy'],
                'carbon_pos': carbon_config[label]['carbon'],
                'vac_pos': carbon_config[label]['vacancy'],
                }

        write_poscar_from_config(config, **carbon_args)

        # calc #Fe(n-q)C(p)Vac(q) 
        atoms = read(f'{struct_dir}/POSCAR_{label}', format='vasp')

        ase_atom_relaxer = aar_from_config(config, calc,opt=config["carbon"]["opt"], logfile = f'{log_dir}/{label}_FeCVac.log')
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.info['conv'] = conv
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{label}', atoms, format='vasp')

        if not conv:
            warnings.warn(f'{idx+1}-th structure, label({label}) did not converge in {config["opt"]["ortho"]["steps"]}steps\n')  

        E_FeCVac = atoms.info['e_fr_energy']
        n_FeCVac = len(atoms)

        del  ase_atom_relaxer, atoms
        gc.collect()

        # equation 3
        n_Vac, n_C = carbon_args["n_vac"], carbon_args["n_carbon"]
        E_bind = n_Vac * E_FeVac + n_C * E_FeC - (n_C + n_Vac - 1) * E_Fe - E_FeCVac

        csv_file.write(f'{label},{E_bind},{E_FeVac},{n_FeVac},{E_FeC},{n_FeC},{E_Fe},{n_Fe},{E_FeCVac},{n_FeCVac},{n_C},{n_Vac}\n')

    torch.cuda.empty_cache()
    csv_file.close()
    write(f'{struct_dir}/FeCVac_opt.extxyz',[read(f'{struct_dir}/CONTCAR_{label}') for label in labels])
    write(f'{struct_dir}/FeCVac.extxyz',[read(f'{struct_dir}/POSCAR_{label}') for label in labels])

def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    args = parse_base_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/febench_carbon_config.yaml')
    calc = calc_from_config(config)

    process_carbon(config, calc)


if __name__ == '__main__':
    main()
