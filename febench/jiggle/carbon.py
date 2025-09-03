from ase import Atoms
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

from febench.calculator.loader import load_calc

def main(argv: list[str] | None=None) -> None:
    args = parse_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = parse_config(config)
    dumpYAML(config, f'{config["cwd"]}/febench_carbon_jiggle.yaml')
    calc = load_calc(config)

    conf = config["carbon"]["jiggle"]
    
    save_dir = config["carbon"]["save"]
    struct_dir = f'{config["carbon"]["save"]}/structure'

    ptb_dir = conf['save']
    ptb_struct = f"{conf['save']}/structure"

    Fe_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    a = Fe_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]
    E_Fe = Fe_bulk.info['e_fr_energy']
    n_Fe = len(Fe_bulk)

    Fe_Vac = read(f'{config["pureFe"]["save"]}/structure/Vac_opt.extxyz')
    E_FeVac = Fe_Vac.info['e_fr_energy']
    n_FeVac = len(Fe_Vac)

    Fe_C = read(f'{struct_dir}/FeC.extxyz', format='extxyz')
    E_FeC = Fe_C.info['e_fr_energy']
    n_FeC = len(Fe_C)

    del Fe_bulk, Fe_Vac, Fe_C
    gc.collect()


    carbon_config = config['carbon_config']

    if conf['cont']:
        csv_file = open(f'{ptb_dir}/carbon_ptb.csv', 'a', buffering = 1)
    else:
        csv_file = open(f'{ptb_dir}/carbon.csv', 'w', buffering = 1)
        csv_file.write('config,E_bind,E_FeVac,n_FeVac,E_FeC,n_FeC,E_Fe,n_Fe,E_FeCVac,n_FeCVac,n_carbon,n_vacancy,conv\n')


    label, nums = conf["label"], conf["num"]
    carbon_args = config["carbon_config"][label]
    for num in tqdm(range(nums), desc='jiggle jiggle ...'):
        logfile=f'{ptb_dir}/{num+1}_relax.log'
        trajfile=f'{ptb_dir}/{num+1}_traj.traj'
        atoms = read(f'{ptb_struct}/POSCAR_{label}_{num+1}', format='vasp')
        ase_relaxer = aar_from_config(config, calc, logfile = logfile, trajfile=trajfile, opt_type='tm')
        atoms, conv = ase_relaxer.relax_atoms(atoms)
        atoms = ase_relaxer.update_atoms(atoms)
        atoms.info['conv'] = conv
        atoms.calc = None
        write(f'{ptb_struct}/CONTCAR_{label}_{num+1}', atoms, format='vasp')

        E_FeCVac = atoms.info['e_fr_energy']
        n_FeCVac = len(atoms)

        # equation 3
        n_Vac, n_C = carbon_args["n_vac"], carbon_args["n_carbon"]
        E_bind = n_Vac * E_FeVac + n_C * E_FeC - (n_C + n_Vac - 1) * E_Fe - E_FeCVac

        del  ase_relaxer, atoms, carbon_args
        gc.collect()
 
        csv_file.write(f'{label}_{num}+1,{E_bind},{E_FeVac},{n_FeVac},{E_FeC},{n_FeC},{E_Fe},{n_Fe},{E_FeCVac},{n_FeCVac},{n_C},{n_Vac},{conv}\n')

    torch.cuda.empty_cache()
    csv_file.close()

if __name__ == '__main__':
    main()
