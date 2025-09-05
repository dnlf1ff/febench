from ase.io import read, write
from ase.build import make_supercell
import os, sys, gc
import numpy as np
from contextlib import redirect_stdout, redirect_stderr
from tqdm import tqdm

from febench.util.relax import aar_from_config
from febench.util.utils import dumpYAML
from febench.util.parser import parse_args
from febench.util.parser import parse_config

from febench.pureFe.utils import write_fe_base, write_csv, EvAToJm


def process_bulk(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'
    csv_file = None

    input_atoms = read(config["data"]["input"], **config["data"]["load_args"])
    atoms = make_supercell(input_atoms, np.diag(config["pureFe"]["bulk"]["supercell"]))
    ase_relaxer = aar_from_config(config, calc, logfile=f'{log_dir}/bulk_relax.log', opt_type='bulk')
   
    if config['pureFe']['cont']:
        csv_file = open(f'{save_dir}/bulk.csv', 'a', buffering=1)
    else:
        csv_file = open(f'{save_dir}/bulk.csv', 'w', buffering=1)
        csv_file.write('idx,energy,surface_area,natom,a,b,c,alpha,beta,gamma,opt_fa,opt_step,force_cnv\n')

    atoms = ase_relaxer.update_atoms(atoms)
    atoms.calc = None

    write_csv(csv_file, atoms, idx='bulk-pre')
    input_list = [atoms]
    for idx, atoms in enumerate(tqdm(input_list, desc = 'relaxing bulk structure ...')):
        atoms.calc = calc
        atoms = ase_relaxer.relax_atoms(atoms)
        atoms = ase_relaxer.update_atoms(atoms)

        conv=f"{atoms.info['opt_fa']},{atoms.info['opt_step']},{atoms.info['force_cnv']}"
        atoms.calc = None
    
        write(f"{struct_dir}/bulk_opt.extxyz", atoms, format='extxyz')
        write(f"{struct_dir}/CONTCAR_bulk", atoms, format='vasp')
        write_csv(csv_file, atoms, idx='bulk-post', conv=conv)

        a0 = atoms.info['a']/config['pureFe']['bulk']['supercell'][0]
        write_fe_base(config, a0)

        csv_file.close()

        del atoms, csv_file
        gc.collect()

def process_vacancy(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    csv_file = open(f'{save_dir}/bulk.csv', 'a', buffering=1)

    atoms=read(f'{struct_dir}/CONTCAR_bulk', **config["data"]["load_args"])
    ase_relaxer = aar_from_config(config, calc,  logfile=f'{log_dir}/Vac_relax.log', opt_type='vacancy')

    del atoms[0]
    atoms = ase_relaxer.update_atoms(atoms)
    atoms.calc = None
    write(f"{struct_dir}/POSCAR_Vac", atoms, format='vasp')
    write_csv(csv_file, atoms, idx='vac-pre')

    input_list = [atoms]
    for idx, atoms in enumerate(tqdm(input_list, desc = 'relaxing sturcture with one vacancy ...')):
        atoms.calc = calc
        atoms = ase_relaxer.relax_atoms(atoms)
        atoms = ase_relaxer.update_atoms(atoms)
        atoms.calc = None
        write(f"{struct_dir}/Vac_opt.extxyz", atoms, format='extxyz')
        write(f"{struct_dir}/CONTCAR_Vac", atoms, format='vasp')
        write_csv(csv_file, atoms, idx='vac-post')
        csv_file.close()

        del atoms, csv_file
        gc.collect()


def main(argv: list[str] | None=None) -> None:
    from febench.calculator.loader import load_calc
    import yaml

    args = parse_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    
    config = parse_config(config)
    dumpYAML(config, f'{config["cwd"]}/config_pure.yaml')
    calc = load_calc(config)
    
    if config['pureFe']['bulk']['run']:
        process_bulk(config, calc)

    if config['pureFe']['vacancy']['run']:
        process_vacancy(config, calc)

if __name__ == '__main__':
    main()
