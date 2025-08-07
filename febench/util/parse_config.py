import os

def check_calc_config(config):
    config_calc = config['calculator']
    assert config_calc['calc_type'].lower() in ['sevennet', 'sevennet-batch', 'custom']
    assert isinstance(config_calc['path'], str)
    assert os.path.isfile(potential := f'{config_calc["path"]}/{config_calc["prefix"]}.{config_calc["extension"]}'), f'MLIP potential not found at {potential}'


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
    config_SnV = config['tm']
    os.makedirs(config_SnV['save'], exist_ok = True)
    os.makedirs(f"{config_SnV['save']}/structure", exist_ok = True)
    os.makedirs(f"{config_SnV['save']}/log", exist_ok = True)

def update_config_dirs(config):
    prefix = config['calculator']['prefix']

    if config['calculator'].get('modal', None) in ['mpa', 'omat24']: # for now
        config['calculator']['calc_args']['modal'] = config['calculator']['modal']
        prefix = f"{config['calculator']['prefix']}/{config['calculator']['modal']}"

    config['data']['cwd'] = (cwd := f"./mlip/{prefix}")

    os.makedirs(cwd, exist_ok=True)
    os.makedirs('output', exist_ok=True)
    tasks = ['pureFe', 'carbon', 'tm']

    for opt_type in config['opt'].keys():
        if opt_type == 'stiffness':
            continue
        config['opt'][opt_type]['logfile'] = cwd + f"/{config['opt'][opt_type]['logfile']}"

    for task in tasks:
        if (save_path := config[task].get('save', None)) is not None:
            config[task]['save'] = f"{cwd}/{save_path}"
    return config

def parse_config_yaml(config):
    config = update_config_dirs(config)

    check_data_config(config)
    check_calc_config(config)
    check_pure_config(config)
    check_carbon_config(config)
    check_tm_config(config)

    config['root'] = os.path.abspath(os.getcwd()) # short stopper
    config['cwd'] = os.path.join(os.path.abspath(os.getcwd()), config['data']['cwd'])
    config['output'] = os.path.join(config['root'], 'output')

    return config
