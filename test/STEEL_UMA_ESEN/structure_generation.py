from ase.io import read, write
import os, sys, gc
import numpy as np
from ase.lattice.cubic import BodyCenteredCubic
from ase.build import bcc100, bcc110, bcc111, bulk, make_supercell
from febench.util.parse_args import parse_base_args
from febench.util.parse_config import parse_config_yaml
import yaml

def get_slab(hkl, a0, vacuum=20, size=(4,4,20)):
    if hkl=='100':
        slab=bcc100('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    elif hkl=='110':
        slab=bcc110('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    elif hkl=='111':
        slab=bcc111('Fe', size=size, a=a0, vacuum=vacuum, orthogonal=True)
    else:
        raise NotImplementedError
    return slab
 
def generate_vacancy(config):
    save_dir = config["pureFe"]["save"]
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    atoms=read(f'{struct_dir}/CONTCAR_bulk', **config["data"]["load_args"])
    del atoms[0]
    write(f"{struct_dir}/POSCAR_Vac", atoms, format='vasp')


def generate_surfaces(config):
    save_dir = config["pureFe"]["save"]
    struct_dir = f'{config["pureFe"]["save"]}/structure'

    atoms = read(f'{struct_dir}/bulk_opt.extxyz', format='extxyz')
    a0 = atoms.info['a']/config['pureFe']['bulk']['supercell'][0]
    
    for hkl in ['100', '110', '111']:
        slab=get_slab(hkl, a0=a0)
        write(f"{struct_dir}/POSCAR_surface_{hkl}", slab, format='vasp')


def generate_carbon_configs(config):
    struct_dir = f'{config["pureFe"]["save"]}/structure'
    atoms = read(f'{struct_dir}/bulk_opt.extxyz', format='extxyz')
    a0 = atoms.info['a']/config['pureFe']['bulk']['supercell'][0]
    from febench.carbon.utils import write_FeC_poscar
    from febench.carbon.utils import write_poscar_from_config as write_carbon_poscars

    write_FeC_poscar(config, a0)
    carbon_config = config['carbon_config']
    for label in config['carbon']['label']:
        carbon_args = {
            'a': a0,
            'label': label,
            'n_carbon': carbon_config[label]['n_carbon'],
            'n_vac': carbon_config[label]['n_vacancy'],
            'carbon_pos': carbon_config[label]['carbon'],
            'vac_pos': carbon_config[label]['vacancy'],
            }

        write_carbon_poscars(config, **carbon_args)

def generate_tm_configs(config): 
    struct_dir = f'{config["pureFe"]["save"]}/structure'
    atoms = read(f'{struct_dir}/bulk_opt.extxyz', format='extxyz')
    a0 = atoms.info['a']/config['pureFe']['bulk']['supercell'][0]
    from febench.tm.utils import write_poscar_from_config as write_tm_poscars
    for solute in config['tm']['solute']:
        write_tm_poscars(config, solute, a0)

def main(argv: list[str] | None=None) -> None:
    args = parse_base_args(argv)
    config_dir = args.config 

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
   
    config = parse_config_yaml(config)
    generate_vacancy(config)
    generate_surfaces(config)
    generate_carbon_configs(config)
    generate_tm_configs(config)

if __name__ == '__main__':
    main()

