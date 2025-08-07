from ase import Atoms
from ase.io import write, read
from ase.build import make_supercell 
from ase.lattice.cubic import BodyCenteredCubic
import yaml

from febench.util.parse_args import parse_base_args
from febench.util.utils import dumpYAML 
from febench.util.parse_config import parse_config_yaml 

from febench.util.relax import aar_from_config
from febench.tm.utils import write_poscar_from_config
import pandas as pd
import numpy as np
import gc
import torch

def process_tm(config, calc):
    save_dir = config["tm"]["save"]
    struct_dir = f'{config["tm"]["save"]}/structure'
    log_dir = f'{config["tm"]["save"]}/log'

    atoms_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    E_Fe = atoms_bulk.info['e_fr_energy']
    a = atoms_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]

    atoms_vac = read(f'{config["pureFe"]["save"]}/structure/vac_opt.extxyz')
    E_FeVac = atoms_vac.info['e_fr_energy']

    del atoms_bulk, atoms_vac
    gc.collect()

    csv_file = open(f'{save_dir}/tm.csv', 'w', buffering = 1)
    csv_file.write('solute,E_Fe,E_FeVac,E_FeM,' + ','.join(f'E_FeMM_{i+1}nn' for i in range(5)) + ','.join(f'E_FeMVac_{i+1}nn' for i in range(5))+'\n')

    ss_file = open(f'{save_dir}/ss.csv', 'w', buffering = 1)
    ss_file.write('solute,' + ','.join(f'E_ss_{i+1}nn' for i in range(5)) + '\n')

    sv_file = open(f'{save_dir}/sv.csv', 'w', buffering = 1)
    sv_file.write('solute,' + ','.join(f'E_sv_{i+1}nn' for i in range(5)) + '\n')

    sols = config["tm"]["solute"]
    for sol in sols:
        write_poscar_from_config(config, sol, a)

        # calc Fe(n-1)M
        atoms = read(f'{struct_dir}/POSCAR_{sol}', format='vasp')

        ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_relax.log')
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.info['conv'] = conv
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{sol}', atoms, format='vasp')
        write(f'{struct_dir}/{sol}_opt.extxyz', atoms, format='extxyz')

        E_FeM = atoms.info['e_fr_energy']

        del  atoms, ase_atom_relaxer
        gc.collect()
        
        E_vac = f'' # E_FeMVac_1nn ... E_FeMVac_5nn
        E_M = f'' # E_FeMM_1nn ... E_FeMM_5nn

        E_ss= f'' # equation 4
        E_sv = f'' # equation 5

        for i in range(5): 
            # calc Fe(n-2)MVac
            atoms_vac = read(f'{struct_dir}/POSCAR_{sol}_vac_{i+1}nn', format='vasp')

            ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_vac_{i+1}nn.log')
            atoms_vac, conv = ase_atom_relaxer.relax_atoms(atoms_vac)
            atoms_vac = ase_atom_relaxer.update_atoms(atoms_vac)
            atoms_vac.info['conv'] = conv
            atoms_vac.calc = None
            write(f'{struct_dir}/CONTCAR_{sol}_vac_{i+1}nn', atoms_vac, format='vasp')

            E_vac += f',{atoms_vac.info["e_fr_energy"]}'
            E_sv += f',{E_FeVac + E_FeM - E_Fe - atoms_vac.info["e_fr_energy"]}'
            del  atoms_vac, ase_atom_relaxer
            gc.collect()

            # calc Fe(n-2)M(2)
            atoms_sol = read(f'{struct_dir}/POSCAR_{sol}_{sol}_{i+1}nn', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_{sol}_{i+1}nn.log')
            atoms_sol, conv = ase_atom_relaxer.relax_atoms(atoms_sol)
            atoms_sol = ase_atom_relaxer.update_atoms(atoms_sol)
            atoms_sol.info['conv'] = conv
            atoms_sol.calc = None
            write(f'{struct_dir}/CONTCAR_{sol}_{sol}_{i+1}nn', atoms_sol, format='vasp')

            E_M += f',{atoms_sol.info["e_fr_energy"]}'
            E_ss += f',{2*E_FeM - E_Fe - atoms_sol.info["e_fr_energy"]}'
            del  atoms_sol, ase_atom_relaxer
            gc.collect()

        csv_file.write(f'{sol},{E_Fe},{E_FeVac},{E_FeM}{E_M}{E_vac}\n')
        sv_file.write(f'{sol}{E_sv}\n')
        ss_file.write(f'{sol}{E_ss}\n')
        torch.cuda.empty_cache()

        write(f'{struct_dir}/{sol}_{sol}.extxyz', [read(f'{struct_dir}/CONTCAR_{sol}_{sol}_{i+1}nn') for i in range(5)]) 
        write(f'{struct_dir}/{sol}_vac.extxyz', [read(f'{struct_dir}/CONTCAR_{sol}_vac_{i+1}nn') for i in range(5)]) 
    csv_file.close()
    sv_file.close()
    ss_file.close()

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
    dumpYAML(config, f'{config["cwd"]}/febench_tm_config.yaml')
    calc = calc_from_config(config)

    process_tm(config, calc)



if __name__ == '__main__':
    main()
