import sys, yaml, gc
from ase.io import read, write

from febench.util.calc import calc_from_config
from febench.util.utils import dumpYAML
from febench.util.parse_config import parse_config_yaml
from febench.util.parse_args import parse_base_args
from febench.pureFe.script import *
from febench.carbon.script import process_carbon
from febench.tm.script import process_tm


import pandas as pd
import torch
import warnings

def main(argv: list[str] | None=None) -> None:
    args = parse_base_args(argv)

    # config.yaml file to read
    config_dir = args.config 


    # see util/calc.py for details
    # mlip type; default is sevennet
    calc_type = args.calc_type 
   
    # will read f'{potential_path}/{calc}.{potential_ext}' as a ASE calc. binary
    # dir named calc will be automatically generated, being the working directory

    # omni, omat, chgTot, o50~  etc
    # customizable if you save potential files somewhere else with a coded name)
    calc = args.calc 

    # dirname of potential file aka what you get when you os.path.dirname($YOUR_POT_FILE)
    # default value is . 
    potential_path = args.potential_path

    # extention format of the potential file default value is pth. some potentials use pt
    potential_ext = args.potential_ext

    # moodal of the potential this will be passed to the ASE.calculator object
    # default value is null which will be ignored
    # for ompa and other models that require a modal key value, pass it in the bash script
    modal = args.modal

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
   
    config['calculator']['calc_type'] = calc_type
    config['calculator']['prefix'] = calc
    config['calculator']['modal'] = modal
    config['calculator']['path'] = potential_path
    config['calculator']['extension'] = potential_ext

    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/config.yaml')

    calc = calc_from_config(config)

    print('processing calculations for pure Iron ...')
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

    if config['carbon']['run']:
        process_carbon(config, calc)

    if config['tm']['run']:
        process_tm(config, calc)
    

if __name__ == '__main__':
    main()
