from ase.io import write, read
from ase.build import make_supercell 
import yaml

from febench.util.parse_args import parse_base_args
from febench.util.utils import dumpYAML 
from febench.util.parse_config import parse_config_yaml 

from febench.util.relax import aar_from_config
from febench.pair.utils import write_poscar_from_config, write_poscar_base
import pandas as pd
import numpy as np
import gc
import torch

def process_pair(config, calc):
    save_dir = config["pair"]["save"]
    struct_dir = f'{config["pair"]["save"]}/structure'
    log_dir = f'{config["pair"]["save"]}/log'

    atoms_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    E_Fe = atoms_bulk.info['e_fr_energy']
    a = atoms_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]

    atoms_vac = read(f'{config["pureFe"]["save"]}/structure/vac_opt.extxyz')
    E_FeVac = atoms_vac.info['e_fr_energy']

    del atoms_bulk, atoms_vac
    gc.collect()

    write_poscar_base(config, a)
    write_poscar_from_config(config, a)

    csv_file = open(f'{save_dir}/pair.csv', 'w', buffering = 1)
    csv_file.write('sol1,sol2,nn,E_bind,E_Fe,E_FeM_1,E_FeM_2,E_tot\n')

    sols = config["pair"]["solute"] + ['Vac']
    for sol in sols:
        if sol != 'Vac':
            atoms = read(f'{struct_dir}/POSCAR_{sol}', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["pair"]["opt"], logfile = f'{log_dir}/A_{sol}.log')
            atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
            atoms = ase_atom_relaxer.update_atoms(atoms)
            atoms.info['conv'] = conv
            atoms.calc = None
            write(f'{struct_dir}/CONTCAR_A_{sol}', atoms, format='vasp')

            E_FeM_1 = atoms.info['e_fr_energy']

            del  atoms, ase_atom_relaxer
            gc.collect()

        else:
            E_FeM_1 = E_FeVac

        for sol2 in sols:
            if sol2 != 'Vac':
                atoms = read(f'{struct_dir}/POSCAR_{sol2}', format='vasp')
                ase_atom_relaxer = aar_from_config(config, calc,opt=config["pair"]["opt"], logfile = f'{log_dir}/B_{sol2}.log')
                atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
                atoms = ase_atom_relaxer.update_atoms(atoms)
                atoms.info['conv'] = conv
                atoms.calc = None
                write(f'{struct_dir}/CONTCAR_B_{sol2}', atoms, format='vasp')

                E_FeM_2 = atoms.info['e_fr_energy']
            else:
                E_FeM_2 = E_FeVac

            for i in range(5): 
                atoms = read(f'{struct_dir}/POSCAR_{sol}_{sol2}_{i+1}nn', format='vasp')
                ase_atom_relaxer = aar_from_config(config, calc,opt=config["pair"]["opt"], logfile = f'{log_dir}/{sol}_{sol2}_{i+1}nn.log')
                atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
                atoms = ase_atom_relaxer.update_atoms(atoms)
                atoms.info['conv'] = conv
                atoms.calc = None
                E_tot = atoms.info['e_fr_energy']
                write(f'{struct_dir}/CONTCAR_{sol}_{sol2}_{i+1}nn', atoms, format='vasp')

                E_b = E_FeM_1 + E_FeM_2 - E_Fe - E_tot

                csv_file.write(f'{sol},{sol2},{i+1},{E_b},{E_Fe},{E_FeM_1},{E_FeM_2},{E_tot}\n')
            del  atoms, ase_atom_relaxer
            gc.collect()
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
    dumpYAML(config, f'{config["cwd"]}/febench_pair_config.yaml')
    calc = calc_from_config(config)

    process_pair(config, calc)



if __name__ == '__main__':
    main()
