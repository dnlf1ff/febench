import os
import warnings
from febench.util.parse_calc_config import check_calc_config
from febench.util.parse_args import parse_base_args


def overwrite_default(config, argv: list[str] | None=None):
    args = parse_base_args(argv)
    config['prefix'] =args.prefix
    config['calculator']['calc_type'] =args.calc_type
    config['calculator']['modal'] =args.modal
    config['calculator']['model'] =args.model
    config['calculator']['dirname'] =args.dirname
    config['calculator']['functional'] =args.functional
    config['calculator']['dispersion'] =args.dispersion
    return config


def check_data_config(config):
    assert os.path.exists(config['data']['input']), 'no input data found'

def check_pure_config(config):
    config_pure = config['pureFe']
    os.makedirs(config_pure['save'], exist_ok = True)
    os.makedirs(f"{config_pure['save']}/structure", exist_ok = True)
    os.makedirs(f"{config_pure['save']}/log", exist_ok = True)

def check_carbon_config(config):
    config_carbon = config['carbon']
    os.makedirs(config_carbon['save'], exist_ok = True)
    os.makedirs(f"{config_carbon['save']}/structure", exist_ok = True)
    os.makedirs(f"{config_carbon['save']}/log", exist_ok = True)

def check_tm_config(config):
    config_tm = config['tm']
    os.makedirs(config_tm['save'], exist_ok = True)
    os.makedirs(f"{config_tm['save']}/structure", exist_ok = True)
    os.makedirs(f"{config_tm['save']}/log", exist_ok = True)

def update_config_dirs(config):
    prefix = config['prefix']
    config['cwd'] = cwd = f"./mlip/{prefix}"

    os.makedirs(cwd, exist_ok=True)
    os.makedirs('output', exist_ok=True)
    tasks = ['pureFe', 'carbon', 'tm']

    for opt_type in config['opt'].keys():
        config['opt'][opt_type]['logfile'] = cwd + f"/{config['opt'][opt_type]['logfile']}"

    for task in tasks:
        if (save_path := config[task].get('save', None)) is not None:
            config[task]['save'] = f"{cwd}/{save_path}"
    return config

def parse_config_yaml(config, argv: list[str] | None=None):
    config = overwrite_default(config, argv)
    config = check_calc_config(config)
    config = update_config_dirs(config)

    check_data_config(config)
    config = check_calc_config(config)

    if config['pureFe']['run']:
        check_pure_config(config)
    if config['carbon']['run']:
        check_carbon_config(config)
    if config['tm']['run']:
        check_tm_config(config)

    config['root'] = os.path.abspath(os.getcwd()) # short stopper
    config['cwd'] = os.path.join(os.path.abspath(os.getcwd()), config['cwd'])
    config['output'] = os.path.join(config['root'], 'output')

    return config
