from ase.io import write, read
from ase.build import make_supercell 
import yaml

from febench.util.parse_args import parse_base_args
from febench.util.utils import dumpYAML 
from febench.util.parse_config import parse_config_yaml 

from febench.util.relax import aar_from_config
from febench.triplet.utils import write_poscar_from_config
import pandas as pd
import numpy as np
import gc
import torch

def process_triplet(config, calc):
    save_dir = config["triplet"]["save"]
    struct_dir = f'{config["triplet"]["save"]}/structure'
    log_dir = f'{config["triplet"]["save"]}/log'

    atoms_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    E_Fe = atoms_bulk.info['e_fr_energy']
    a = atoms_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]

    del atoms_bulk
    gc.collect()

    for triplet in config['triplet']['label']:
        ref_df = pd.read_csv(f'triplet_{triplet}.csv')
        csv_file = open(f'{save_dir}/triplet_{triplet}_E_bind.csv', 'w', buffering = 1)
        csv_file.write('A,B,C,E_bind,E_Fe,E_A,E_B,E_C,E_ABC\n')

        write_poscar_from_config(config, a, triplet)

        for A, B, C in zip(ref_df['A'], ref_df['B'], ref_df['C']):
            A_atoms = read(f'{struct_dir}/POSCAR_{triplet}_A_{A}', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["triplet"]["opt"], logfile = f'{log_dir}/{triplet}_{A}{B}{C}_{A}.log')

            A_atoms, conv = ase_atom_relaxer.relax_atoms(A_atoms)
            A_atoms = ase_atom_relaxer.update_atoms(A_atoms)
            A_atoms.info['conv'] = conv
            A_atoms.calc = None
            write(f'{struct_dir}/CONTCAR_{triplet}_{A}{B}{C}_A_{A}', A_atoms, format='vasp')
            E_A = A_atoms.info['e_fr_energy']

            del A_atoms, ase_atom_relaxer
            gc.collect()

            B_atoms = read(f'{struct_dir}/POSCAR_{triplet}_B_{B}', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["triplet"]["opt"], logfile = f'{log_dir}/{triplet}_{A}{B}{C}_{B}.log')

            B_atoms, conv = ase_atom_relaxer.relax_atoms(B_atoms)
            B_atoms = ase_atom_relaxer.update_atoms(B_atoms)
            B_atoms.info['conv'] = conv
            B_atoms.calc = None
            write(f'{struct_dir}/CONTCAR_{triplet}_{A}{B}{C}_B_{B}', B_atoms, format='vasp')
            E_B = B_atoms.info['e_fr_energy']

            del B_atoms, ase_atom_relaxer
            gc.collect()

            C_atoms = read(f'{struct_dir}/POSCAR_{triplet}_C_{C}', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["triplet"]["opt"], logfile = f'{log_dir}/{triplet}_{A}{B}{C}_{C}.log')

            C_atoms, conv = ase_atom_relaxer.relax_atoms(C_atoms)
            C_atoms = ase_atom_relaxer.update_atoms(C_atoms)
            C_atoms.info['conv'] = conv
            C_atoms.calc = None
            write(f'{struct_dir}/CONTCAR_{triplet}_{A}{B}{C}_C_{C}', C_atoms, format='vasp')
            E_C = C_atoms.info['e_fr_energy']

            del C_atoms, ase_atom_relaxer
            gc.collect()

            ABC_atoms = read(f'{struct_dir}/POSCAR_{triplet}_ABC_{A}{B}{C}', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc, opt=config["triplet"]["opt"], logfile = f'{log_dir}/{triplet}_{A}{B}{C}.log')

            ABC_atoms, conv = ase_atom_relaxer.relax_atoms(ABC_atoms)
            ABC_atoms = ase_atom_relaxer.update_atoms(ABC_atoms)
            ABC_atoms.info['conv'] = conv
            ABC_atoms.calc = None
            write(f'{struct_dir}/CONTCAR_{triplet}_{A}{B}{C}_ABC_{A}{B}{C}', ABC_atoms, format='vasp')
            E_ABC = ABC_atoms.info['e_fr_energy']

            del ABC_atoms, ase_atom_relaxer
            gc.collect()
            
            E_bind = E_ABC - E_A - E_B - E_C + 2 * E_Fe
            csv_file.write(f'{A},{B},{C},{E_bind},{E_Fe},{E_A},{E_B},{E_C},{E_ABC}\n')

        torch.cuda.empty_cache()
        csv_file.close()

def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    args = parse_base_args(argv)
    config_dir = args.config
    calc_type = args.calc_type
    calc = args.calc
    modal = args.modal
    potential_path = args.potential_path
    potential_ext = args.potential_ext

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    if modal.lower() != 'null':
        config['calculator']['modal'] = modal
    
    config['calculator']['calc_type'] = calc_type
    config['calculator']['prefix'] = calc
    config['calculator']['path'] = potential_path
    config['calculator']['extension'] = potential_ext
    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/febench_triplet_config.yaml')
    calc = calc_from_config(config)

    process_triplet(config, calc)



if __name__ == '__main__':
    main()
