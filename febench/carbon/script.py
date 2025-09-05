from ase.io import write, read
import yaml
import gc
import numpy as np

from febench.util.utils import dumpYAML 

from febench.util.parser import parse_args
from febench.util.parser import parse_config 
from febench.util.relax import aar_from_config
from febench.carbon.utils import write_poscar_from_config, write_FeC_poscar

import torch
from tqdm import tqdm
import warnings
import os

def process_carbon(config, calc):
    save_dir = config["carbon"]["save"]
    struct_dir = f'{config["carbon"]["save"]}/structure'
    log_dir = f'{config["carbon"]["save"]}/log'

    Fe_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    a = Fe_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]
    E_Fe = Fe_bulk.info['e_fr_energy']
    n_Fe = len(Fe_bulk)

    Fe_Vac = read(f'{config["pureFe"]["save"]}/structure/Vac_opt.extxyz')
    E_FeVac = Fe_Vac.info['e_fr_energy']
    n_FeVac = len(Fe_Vac)

    carbon_config = config['carbon_config']

    if config['carbon']['cont']:
        csv_file = open(f'{save_dir}/carbon.csv', 'a', buffering = 1)
    else:
        csv_file = open(f'{save_dir}/carbon.csv', 'w', buffering = 1)
        csv_file.write('config,E_bind,E_FeVac,n_FeVac,E_FeC,n_FeC,E_Fe,n_Fe,E_FeCVac,n_FeCVac,n_carbon,n_vacancy,opt_fa,opt_step,force_cnv\n')

    if config['carbon']['cont']:
        FeC = read(f'{struct_dir}/FeC.extxyz', format='extxyz')
    else:
        write_FeC_poscar(config, a)
        FeC = read(f'{struct_dir}/POSCAR_C', format='vasp')
        ase_relaxer = aar_from_config(config, calc, logfile = f'{log_dir}/FeC.log', opt_type='carbon')
        FeC = ase_relaxer.relax_atoms(FeC)
        FeC = ase_relaxer.update_atoms(FeC)
        FeC.calc = None
        write(f'{struct_dir}/CONTCAR_C', FeC, format='vasp')
        write(f'{struct_dir}/FeC.extxyz', FeC, format='extxyz')
        del ase_relaxer

    E_FeC = FeC.info['e_fr_energy']
    n_FeC = len(FeC)
    csv_file.write(f'# FeC: energy:{E_FeC}, opt_fa: {FeC.info["opt_fa"]}, opt_steps: {FeC.info["opt_step"]}, force_cnv: {FeC.info["force_cnv"]}\n')

    del Fe_bulk, Fe_Vac, FeC
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

        atoms = read(f'{struct_dir}/POSCAR_{label}', format='vasp')
        ase_relaxer = aar_from_config(config, calc, logfile = f'{log_dir}/{label}_relax.log',)
        atoms = ase_relaxer.relax_atoms(atoms)
        atoms = ase_relaxer.update_atoms(atoms)
        conv=f"{atoms.info['opt_fa']},{atoms.info['opt_step']},{atoms.info['force_cnv']}"
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{label}', atoms, format='vasp')

        E_FeCVac = atoms.info['e_fr_energy']
        n_FeCVac = len(atoms)

        # equation 3
        n_Vac, n_C = carbon_args["n_vac"], carbon_args["n_carbon"]
        E_bind = n_Vac * E_FeVac + n_C * E_FeC - (n_C + n_Vac - 1) * E_Fe - E_FeCVac

        del  ase_relaxer, atoms, carbon_args
        gc.collect()

 
        csv_file.write(f'{label},{E_bind},{E_FeVac},{n_FeVac},{E_FeC},{n_FeC},{E_Fe},{n_Fe},{E_FeCVac},{n_FeCVac},{n_C},{n_Vac},{conv}\n')

    torch.cuda.empty_cache()
    csv_file.close()
    write(f'{struct_dir}/FeCVac_opt.extxyz',[read(f'{struct_dir}/CONTCAR_{label}') for label in labels])
    write(f'{struct_dir}/FeCVac.extxyz',[read(f'{struct_dir}/POSCAR_{label}') for label in labels])

def main(argv: list[str] | None=None) -> None:
    from febench.calculator.loader import load_calc
    args = parse_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = parse_config(config)
    dumpYAML(config, f'{config["cwd"]}/febench_carbon_config.yaml')
    calc = load_calc(config)

    process_carbon(config, calc)


if __name__ == '__main__':
    main()
