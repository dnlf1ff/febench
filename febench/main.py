import sys, yaml, gc
from ase.io import read, write

from febench.util.calc import calc_from_config
from febench.util.utils import dumpYAML
from febench.util.parse_config import parse_config_yaml
from febench.util.parse_args import parse_base_args
from febench.pureFe.script import *
from febench.carbon.script import process_carbon
from febench.tm.script import process_tm
from febench.pair.script import process_pair
from febench.triplet.script import process_triplet


import pandas as pd
import torch
import warnings

def main(argv: list[str] | None=None) -> None:
    args = parse_base_args(argv)

    # config.yaml file to read
    config_dir = args.config 

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
   
    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/config.yaml')

    calc = calc_from_config(config)

    print('processing calculations for pure Iron ...')
    if config['pureFe']['run']:
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

    if config['pair']['run']:
        process_pair(config, calc)

    if config['triplet']['run']:
        process_triplet(config, calc)
    

if __name__ == '__main__':
    main()
