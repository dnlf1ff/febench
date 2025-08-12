from ase.io import read, write
from ase.build import make_supercell
from ase.optimize import FIRE
import os, sys, gc
import numpy as np
import pandas as pd
from contextlib import redirect_stdout, redirect_stderr
from tqdm import tqdm

from febench.util.relax import aar_from_config
from febench.util.utils import dumpYAML
from febench.util.parse_args import parse_base_args
from febench.util.parse_config import parse_config_yaml

from febench.pureFe.utils import write_fe_base, get_slab, write_csv, EvAToJm


def process_bulk(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

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

        del input_atoms, atoms, csv_file
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

def process_surfaces(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    csv_file = open(f"{save_dir}/bulk.csv", "a", buffering = 1)

    atoms = read(f'{struct_dir}/bulk_opt.extxyz', format='extxyz')
    a0 = atoms.info['a']/config['pureFe']['bulk']['supercell'][0]
 
    for idx, hkl in enumerate(tqdm(config["pureFe"]["surface"]["hkl"], desc ='relaxing surfaces ...')):
        ase_atom_relaxer = aar_from_config(config, calc, opt=config["pureFe"]["surface"]["opt"], logfile=f'{log_dir}/surface_{hkl}.log')
        slab=get_slab(hkl, a0=a0)
        slab = ase_atom_relaxer.update_atoms(slab)
        slab.calc = None
        write(f"{struct_dir}/surface_{hkl}.extxyz", slab, format='extxyz')
        write(f"{struct_dir}/POSCAR_surface_{hkl}", slab, format='vasp')
        write_csv(csv_file, slab, idx=f'{hkl}-pre')

        slab.calc = calc
        slab, conv = ase_atom_relaxer.relax_atoms(slab)
        slab = ase_atom_relaxer.update_atoms(slab)
        slab.calc = None
        slab.info['conv'] = conv
        write(f"{struct_dir}/surface_{hkl}_opt.extxyz", slab, format='extxyz')
        write(f"{struct_dir}/CONTCAR_surface_{hkl}", slab, format='vasp')
        write_csv(csv_file, slab, idx=f'{hkl}-post')

        del slab, ase_atom_relaxer
        gc.collect()
    csv_file.close()
    del csv_file
    gc.collect()

def post_process(config):
    # ugly and explicit .. only for pureFe
    save_dir = config["pureFe"]["save"]
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    atoms_bulk = read(f'{struct_dir}/bulk_opt.extxyz', format='extxyz')
    a = atoms_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]
    E_bulk = atoms_bulk.info['e_fr_energy']
    n_bulk = len(atoms_bulk)

    atoms_Vac = read(f'{struct_dir}/Vac_opt.extxyz', format='extxyz') 
    E_Vac = atoms_Vac.info['e_fr_energy']
    E_Vac_f = E_Vac - (n_bulk - 1) * E_bulk / n_bulk

    property_vals = [a, E_Vac_f]

    for hkl in config['pureFe']['surface']['hkl']:
        atoms_hkl = read(f'{struct_dir}/surface_{hkl}_opt.extxyz', format='extxyz')
        E_hkl = atoms_hkl.info['e_fr_energy']
        A_hkl = atoms_hkl.info['surface_area']
        E_ref = E_bulk * 320/128
        E_hkl_f = (E_hkl - E_ref)/ (2*A_hkl) * EvAToJm
        property_vals.append(E_hkl_f)

    property_keys = ['a0','E_vac_f', 'E_100_f', 'E_110_f', 'E_111_f']
    df = pd.DataFrame(index=property_keys)

    df[config['prefix']] = property_vals

    df.index.name = 'property'
    df.to_csv(f'{save_dir}/pureFe_properties.csv', index_label='property')

def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    import yaml

    args = parse_base_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    
    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/config_pure.yaml')
    calc = calc_from_config(config)
    
    if config['pureFe']['bulk']['run']:
        process_bulk(config, calc)

    if config['pureFe']['vacancy']['run']:
        process_vacancy(config, calc)

    if config['pureFe']['surface']['run']:
        process_surfaces(config, calc)

    if config['pureFe']['post']['run']:
        post_process(config)

if __name__ == '__main__':
    main()
