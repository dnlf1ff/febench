import sys, yaml, gc
from ase.io import read, write

from febench.util.utils import dumpYAML
from febench.util.parser import parse_config
from febench.util.parser import parse_args

from febench.calculator.loader import load_calc 

from febench.pureFe.script import *
from febench.carbon.script import process_carbon
from febench.tm.script import process_tm

import pandas as pd
import torch
import warnings

def main(argv: list[str] | None=None) -> None:
    args = parse_args(argv)

    # config.yaml file to read
    config_dir = args.config 

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)
   
    config = parse_config(config)
    dumpYAML(config, f'{config["cwd"]}/config_parsed.yaml')

    calc = load_calc(config)

    print('processing calculations for pure Iron ...')
    if config['pureFe']['run']:
        if config['pureFe']['bulk']['run']:
            process_bulk(config, calc)

        if config['pureFe']['vacancy']['run']:
            process_vacancy(config, calc)

    if config['carbon']['run']:
        process_carbon(config, calc)

    if config['tm']['run']:
        process_tm(config, calc)

if __name__ == '__main__':
    main()
