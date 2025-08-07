from ase.io import read, write
from ase.build import make_supercell
from ase.optimize import FIRE
import os, sys, gc
import numpy as np
import pandas as pd
from contextlib import redirect_stdout, redirect_stderr

from febench.util.relax import aar_from_config
from febench.util.utils import dumpYAML
from febench.util.parse_args import parse_base_args
from febench.util.parse_config import parse_config_yaml

from febench.pureFe.utils import get_slab, get_surface_ref_bulk, get_Cij, get_elastic_constants, write_csv

from matscipy.elasticity import fit_elastic_constants

def process_bulk(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    input_atoms = read(config["data"]["input"], **config["data"]["load_args"])
    atoms = make_supercell(input_atoms, np.diag(config["pureFe"]["bulk"]["supercell"]))
    ase_atom_relaxer = aar_from_config(config, calc, opt=config["pureFe"]["bulk"]["opt"], logfile=f'{log_dir}/bulk_relax.log')
    
    csv_file = open(f'{save_dir}/bulk.csv', 'w', buffering=1)
    csv_file.write('idx,energy,volume,surface_area,natom,a,b,c,alpha,beta,gamma,conv\n')

    atoms = ase_atom_relaxer.update_atoms(atoms)
    atoms.calc = None

    write_csv(csv_file, atoms, idx='bulk-pre')

    print('running bulk relaxation')

    atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
    atoms = ase_atom_relaxer.update_atoms(atoms)
    atoms.info['conv'] = conv
    atoms.calc = None
    
    write(f"{struct_dir}/bulk_opt.extxyz", atoms, format='extxyz')
    write(f"{struct_dir}/CONTCAR_bulk", atoms, format='vasp')
    write_csv(csv_file, atoms, idx='bulk-post')

    print('relaxation for 4 by 4 by 4 cell -- done')
    csv_file.close()

    del input_atoms, atoms, csv_file
    gc.collect()

def process_vacancy(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    print('running pure Iron sturcture with a vacancy')
    csv_file = open(f'{save_dir}/bulk.csv', 'a', buffering=1)

    atoms=read(f'{struct_dir}/CONTCAR_bulk', **config["data"]["load_args"])
    ase_atom_relaxer = aar_from_config(config, calc, opt=config["pureFe"]["vacancy"]["opt"], logfile=f'{log_dir}/vacancy_relax.log')

    del atoms[0]
    atoms = ase_atom_relaxer.update_atoms(atoms)
    atoms.calc = None
    write(f"{struct_dir}/POSCAR_vac", atoms, format='vasp')
    write_csv(csv_file, atoms, idx='vac-pre')

    atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
    atoms = ase_atom_relaxer.update_atoms(atoms)
    atoms.calc = None
    atoms.info['conv'] = conv
    write(f"{struct_dir}/vac_opt.extxyz", atoms, format='extxyz')
    write(f"{struct_dir}/CONTCAR_vac", atoms, format='vasp')
    write_csv(csv_file, atoms, idx='vac-post')
    csv_file.close()

    del atoms, csv_file
    gc.collect()

def process_surfaces(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    csv_file = open(f"{save_dir}/bulk.csv", "a", buffering = 1)

    df = pd.read_csv(f"{save_dir}/bulk.csv")
    a0 = df['a'][1]/config['pureFe']['bulk']['supercell'][0]

    for hkl in config["pureFe"]["surface"]["hkl"]:
        ase_atom_relaxer = aar_from_config(config, calc, opt=config["pureFe"]["surface"]["opt"], logfile=f'{log_dir}/surface_{hkl}.log')
        slab=get_slab(hkl, a0=a0)

        slab = ase_atom_relaxer.update_atoms(slab)
        slab.calc = None
        write(f"{struct_dir}/POSCAR_surface_{hkl}", slab, format='vasp')
        write_csv(csv_file, slab, idx=f'{hkl}-pre')

        slab, conv = ase_atom_relaxer.relax_atoms(slab)
        slab = ase_atom_relaxer.update_atoms(slab)
        slab.calc = None
        slab.info['conv'] = conv
        write(f"{struct_dir}/surface_{hkl}_opt.extxyz", slab, format='extxyz')
        write(f"{struct_dir}/CONTCAR_surface_{hkl}", slab, format='vasp')
        write_csv(csv_file, slab, idx=f'{hkl}-post')

        del slab, ase_atom_relaxer
        gc.collect()
        print(f'surface energy for {hkl} -- done')
    csv_file.close()
    del csv_file
    gc.collect()

def process_stiffness(config, calc):
    save_dir = config["pureFe"]["save"]
    log_dir = f'{config["pureFe"]["save"]}/log'
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    atoms = read(f'{struct_dir}/CONTCAR_bulk', format='vasp')
    atoms.calc = calc

    with open(f'{log_dir}/elastic_stdout.x', 'w') as f, redirect_stdout(f), redirect_stderr(f):
        C_least_squares, _ = fit_elastic_constants(atoms,optimizer=FIRE, **config['opt']['stiffness'], logfile=f'{log_dir}/elastic.log')

    os.chdir(config['root'])

    tensor_df = get_Cij(C_least_squares)
    tensor_df.to_csv(f'{save_dir}/Cij.csv', index_label='row-column')

def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    import yaml
    args = parse_base_args(argv)
    config_dir = args.config
    calc_type = args.calc_type
    calc = args.calc
    modal = args.modal
    potential_path = args.potential_path
    potential_ext = args.potential_ext

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config['calculator']['calc_type'] = calc_type
    config['calculator']['prefix'] = calc
    config['calculator']['path'] = potential_path
    config['calculator']['extension'] = potential_ext

    if modal.lower() != 'null':
        config['calculator']['modal'] = modal
    
    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/config_pure.yaml')
    calc = calc_from_config(config)
    
    if config['pureFe']['bulk']['run']:
        process_bulk(config, calc)

    if config['pureFe']['vacancy']['run']:
        process_vacancy(config, calc)

    if config['pureFe']['surface']['run']:
        process_surfaces(config, calc)

    if config['pureFe']['stiffness']['run']:
        process_stiffness(config, calc)


if __name__ == '__main__':
    main()
