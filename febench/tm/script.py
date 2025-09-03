from ase import Atoms
from ase.io import write, read
from ase.build import make_supercell 
from ase.lattice.cubic import BodyCenteredCubic
import yaml

from febench.util.parser import parse_args
from febench.util.utils import dumpYAML 
from febench.util.parser import parse_config 

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

    dir_args = {'save_dir': save_dir, 'struct_dir': struct_dir,
                'log_dir': log_dir}
    Fe_bulk = read(f'{config["pureFe"]["save"]}/structure/bulk_opt.extxyz')
    E_Fe = Fe_bulk.info['e_fr_energy']
    a = Fe_bulk.info['a']/config['pureFe']['bulk']['supercell'][0]

    Fe_Vac = read(f'{config["pureFe"]["save"]}/structure/Vac_opt.extxyz')
    E_FeVac = Fe_Vac.info['e_fr_energy']

    del Fe_bulk, Fe_Vac
    gc.collect()

    if config['tm']['cont']:
        tm_file = open(f'{save_dir}/tm_E_bind.csv', 'a', buffering = 1)
    else:
        tm_file = open(f'{save_dir}/tm_E_bind.csv', 'w', buffering = 1)
        tm_file.write('sol,E_bind,conv\n')

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

        atoms = read(f'{struct_dir}/POSCAR_{sol}_{sol}_1nn', format='vasp')
        ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_{sol}_1nn_relax.log')
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)

        if not conv:
            warnings.warn(f'1nn of {sol}-{sol}, did not converge in {config["opt"]["ortho"]["steps"]}steps\n')  

        atoms.info['conv'] = conv
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{sol}_{sol}_1nn', atoms, format='vasp')

        E_FeMM = atoms.info["e_fr_energy"]
        E_bind = 2 * E_FeM - E_Fe - atoms.info['e_fr_energy']
        tm_file.write(f'{sol},1,{E_bind},{conv}\n')
        del  atoms, ase_atom_relaxer
        gc.collect()

        torch.cuda.empty_cache()


def main(argv: list[str] | None=None) -> None:
    from febench.calculator.loader import load_calc
    args = parse_args(argv)
    config_dir = args.config

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config = parse_config(config)
    dumpYAML(config, f'{config["cwd"]}/febench_tm_config.yaml')
    calc = load_calc(config)

    process_tm(config, calc)



if __name__ == '__main__':
    main()
