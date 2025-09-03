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
    ase_atom_relaxer = aar_from_config(config, calc, opt=config["pureFe"]["bulk"]["opt"], logfile=f'{log_dir}/bulk_relax.log')
   
    if config['pureFe']['cont']:
        csv_file = open(f'{save_dir}/bulk.csv', 'a', buffering=1)
    else:
        csv_file = open(f'{save_dir}/bulk.csv', 'w', buffering=1)
        csv_file.write('idx,energy,surface_area,natom,a,b,c,alpha,beta,gamma,conv\n')

    atoms = ase_atom_relaxer.update_atoms(atoms)
    atoms.calc = None

    write_csv(csv_file, atoms, idx='bulk-pre')
    input_list = [atoms]
    for idx, atoms in enumerate(tqdm(input_list, desc = 'relaxing bulk structure ...')):
        atoms.calc = calc
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.info['conv'] = conv
        atoms.calc = None
    
        write(f"{struct_dir}/bulk_opt.extxyz", atoms, format='extxyz')
        write(f"{struct_dir}/CONTCAR_bulk", atoms, format='vasp')
        write_csv(csv_file, atoms, idx='bulk-post')

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
    ase_atom_relaxer = aar_from_config(config, calc, opt=config["pureFe"]["vacancy"]["opt"], logfile=f'{log_dir}/Vac_relax.log')

    del atoms[0]
    atoms = ase_atom_relaxer.update_atoms(atoms)
    atoms.calc = None
    write(f"{struct_dir}/POSCAR_Vac", atoms, format='vasp')
    write_csv(csv_file, atoms, idx='vac-pre')

    input_list = [atoms]
    for idx, atoms in enumerate(tqdm(input_list, desc = 'relaxing sturcture with one vacancy ...')):
        atoms.calc = calc
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.calc = None
        atoms.info['conv'] = conv
        write(f"{struct_dir}/Vac_opt.extxyz", atoms, format='extxyz')
        write(f"{struct_dir}/CONTCAR_Vac", atoms, format='vasp')
        write_csv(csv_file, atoms, idx='vac-post')
        csv_file.close()

        del atoms, csv_file
        gc.collect()


def main(argv: list[str] | None=None) -> None:
    from febench.calculator import load_calc
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
