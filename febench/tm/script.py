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
import numpy as np
import gc
import torch

from tqdm import tqdm
import warnings


def process_tm(config, calc):
    save_dir = config["tm"]["save"]
    struct_dir = f'{config["tm"]["save"]}/structure'
    log_dir = f'{config["tm"]["save"]}/log'

    atoms_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    E_Fe = atoms_bulk.info['e_fr_energy']
    a = atoms_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]

    atoms_Vac = read(f'{config["pureFe"]["save"]}/structure/Vac_opt.extxyz')
    E_FeVac = atoms_Vac.info['e_fr_energy']

    del atoms_bulk, atoms_Vac
    gc.collect()

    if config['tm']['cont']:
        csv_file = open(f'{save_dir}/tm.csv', 'a', buffering = 1)
        tm_file = open(f'{save_dir}/tm_E_bind.csv', 'a', buffering = 1)
    else:
        csv_file = open(f'{save_dir}/tm.csv', 'w', buffering = 1)
        csv_file.write('solute,E_Fe,E_FeVac,E_FeM,FeM_conv,' + ','.join(f'E_FeMM_{i+1}nn' for i in range(5)) + ','.join(f'E_FeMVac_{i+1}nn' for i in range(5))+'\n')

        tm_file = open(f'{save_dir}/tm_E_bind.csv', 'w', buffering = 1)
        tm_file.write('sol_1,sol_2,nn,E_bind,FeMM(Vac)_conv\n')

    sols = config["tm"]["solute"]
    for idx, sol in enumerate(tqdm(sols, desc='processing transition metals ...')):
        write_poscar_from_config(config, sol, a)

        # calc Fe(n-1)M
        atoms = read(f'{struct_dir}/POSCAR_{sol}', format='vasp')

        ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_relax.log')
        atoms, FeM_conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.info['conv'] = FeM_conv
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{sol}', atoms, format='vasp')
        write(f'{struct_dir}/{sol}_opt.extxyz', atoms, format='extxyz')

        E_FeM = atoms.info['e_fr_energy']

        del  atoms, ase_atom_relaxer
        gc.collect()

        E_FeMM = ''

        for i in range(5): 
            # calc Fe(n-2)MVac
            # calc Fe(n-2)M(2)
            atoms = read(f'{struct_dir}/POSCAR_{sol}_{sol}_{i+1}nn', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_{sol}_{i+1}nn_relax.log')
            atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
            atoms = ase_atom_relaxer.update_atoms(atoms)
            if not conv:
                warnings.warn(f'{idx+1}th structure, i.e. {i+1}th nn of {sol}-{sol}, did not converge in {config["opt"]["ortho"]["steps"]}steps\n')  


            atoms.info['conv'] = conv
            atoms.calc = None
            write(f'{struct_dir}/CONTCAR_{sol}_{sol}_{i+1}nn', atoms, format='vasp')

            E_FeMM += f',{atoms.info["e_fr_energy"]}'
            E_bind = 2 * E_FeM - E_Fe - atoms.info['e_fr_energy']
            tm_file.write(f'{sol},{sol},{i+1},{E_bind},{conv}\n')
            del  atoms, ase_atom_relaxer
            gc.collect()

   
        E_FeMVac =''
        for i in range(5): 
            # calc Fe(n-2)MVac
            atoms = read(f'{struct_dir}/POSCAR_{sol}_Vac_{i+1}nn', format='vasp')

            ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_Vac_{i+1}nn_relax.log')
            atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
            atoms = ase_atom_relaxer.update_atoms(atoms)
            atoms.info['conv'] = conv

            if not conv:
                warnings.warn(f'{idx+1}th structure, i.e. {i+1}th nn of {sol}-Vac, did not converge in {config["opt"]["ortho"]["steps"]}steps\n')  
            atoms.calc = None
            write(f'{struct_dir}/CONTCAR_{sol}_Vac_{i+1}nn', atoms, format='vasp')

            E_FeMVac += f',{atoms.info["e_fr_energy"]}'
            E_bind = E_FeVac + E_FeM - E_Fe - atoms.info['e_fr_energy']
            tm_file.write(f'{sol},Vac,{i+1},{E_bind},{conv}\n')
            del  atoms, ase_atom_relaxer
            gc.collect()

        csv_file.write(f'{sol},{E_Fe},{E_FeVac},{E_FeM},{conv}{E_FeMM}{E_FeMVac}\n')
        torch.cuda.empty_cache()

        write(f'{struct_dir}/{sol}_{sol}.extxyz', [read(f'{struct_dir}/CONTCAR_{sol}_{sol}_{i+1}nn') for i in range(5)]) 
        write(f'{struct_dir}/{sol}_Vac.extxyz', [read(f'{struct_dir}/CONTCAR_{sol}_Vac_{i+1}nn') for i in range(5)]) 
    csv_file.close()

def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    args = parse_base_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/febench_tm_config.yaml')
    calc = calc_from_config(config)

    process_tm(config, calc)



if __name__ == '__main__':
    main()
