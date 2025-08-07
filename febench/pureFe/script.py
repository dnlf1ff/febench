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

from febench.pureFe.utils import write_fe_base, get_slab, get_Cij, get_elastic_constants, write_csv,EvAToJm

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

    a0 = atoms.info['a']/config['pureFe']['bulk']['supercell'][0]
    write_fe_base(config, a0)

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

def post_process(config):
    # ugly and explicit .. only for pureFe
    save_dir = config["pureFe"]["save"]

    df = pd.read_csv(f'{save_dir}/bulk.csv', index_col='idx')
    a = df.loc['bulk-post']['a']/config['pureFe']['bulk']['supercell'][0]
    E_bulk = df.loc['bulk-post']['energy']
    n_bulk = df.loc['bulk-post']['natom']
    E_vac = df.loc['vac-post']['energy']
    E_vac_f = E_vac - (n_bulk - 1) * E_bulk / n_bulk

    E_100 = df.loc['100-post']['energy']
    A_100 = df.loc['100-post']['surface_area']
    E_110 = df.loc['110-post']['energy']
    A_110 = df.loc['110-post']['surface_area']
    E_111 = df.loc['111-post']['energy']
    A_111 = df.loc['111-post']['surface_area']
    E_ref = E_bulk * 320/128

    E_100_f = (E_100 - E_ref)/ (2*A_100) * EvAToJm
    E_110_f = (E_110 - E_ref)/ (2*A_110) * EvAToJm
    E_111_f = (E_111 - E_ref)/ (2*A_111) * EvAToJm

    df_Cij= pd.read_csv(f'{save_dir}/Cij.csv', index_col='row-column')

    # B = (C11 + 2*C12)/3
    # C_prime = (C11-C12)/2
    B = (df_Cij.loc['i1']['j1'] + 2 * df_Cij.loc['i1']['j2'])/3
    C_prime = (df_Cij.loc['i1']['j1'] - df_Cij.loc['i1']['j2'])/2
    C44 = df_Cij.loc['i4']['j4']

    property_keys = ['a0','E_vac_f', 'E_100_f', 'E_110_f', 'E_111_f', 'B', 'C_prime', 'C44']
    df_pureFe = pd.DataFrame(index=property_keys)

    pureFe_props = [a, E_vac_f, E_100_f, E_110_f, E_111_f, B, C_prime, C44]
    df_pureFe[config['calculator']['prefix']] = pureFe_props

    df_pureFe.index.name = 'property'
    df_pureFe.to_csv(f'{save_dir}/pureFe_properties.csv', index_label='property')

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

    if config['pureFe']['post']['run']:
        post_process(config)

if __name__ == '__main__':
    main()
