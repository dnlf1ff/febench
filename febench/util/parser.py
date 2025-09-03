import os
import warnings
import argparse

def parse_args(argv: list[str]| None=None):
    parser = argparse.ArgumentParser(description= "cli tool")

    parser.add_argument('--config', type=str, default='./config.yaml', 
                        help='config yaml file directory')

    parser.add_argument('--calc', type=str, default='omni',
                        help='ompa, mace, orb, esen, dpa, uma, omni')

    parser.add_argument('--modal', type=str, default='omat',
                        help='mpa, omat24, matpes_pbe')

    return parser.parse_args(argv)

def overwrite_default(config, argv: list[str] | None=None):
    args = parse_args(argv)
    config['calculator']['calc'] =args.calc.lower()
    config['calculator']['model'] =args.calc.lower()

    if args.calc == 'test':
        config['prefix'] = args.calc.lower()
        config['calculator']['modal'] = 'na'
    else:
        config['prefix'] = f'{args.calc.lower()}/{args.modal.lower()}'
        config['calculator']['modal'] = args.modal.lower()

    return config

def check_data_config(config):
    assert os.path.exists(config['data']['input']), 'no input data found'

def check_calc_config(config):
    conf = config['calculator']
    calc, modal = conf['calc'], conf['modal']
    assert calc in ['ompa', 'mace', 'orb', 'esen', 'dpa', 'uma', 'omni', 'test']
    if calc in ['mace', 'orb']:
        assert modal in ['mpa', 'omat']
    elif calc in ['esen']:
        assert modal in ['oam', 'omat']
    elif calc in ['dpa']:
        assert modal in ['omat', 'mp']
    elif calc in ['uma']:
        assert modal in ['omat', 'omc']
    elif calc in ['ompa']:
        assert modal in ['mpa', 'omat24']
    elif calc in ['omni']:
        assert modal in ['mpa', 'omat24', 'matpes_pbe']
    elif calc == 'test':
        print("TEST MODE: 7net-0 will be automatically loaded")

    if modal == 'omc':
        config['calculator']['dispersion']: 'true'
        config['calculator']['functional']: 'PBE'
    else:
        config['calculator']['dispersion']: 'false'
        config['calculator']['functional']: 'na'

    return config

def check_pure_config(config):
    conf = config['pureFe']
    os.makedirs(conf['save'], exist_ok = True)
    os.makedirs(f"{conf['save']}/structure", exist_ok = True)
    os.makedirs(f"{conf['save']}/log", exist_ok = True)

def check_carbon_config(config):
    conf = config['carbon']
    os.makedirs(conf['save'], exist_ok = True)
    os.makedirs(f"{conf['save']}/structure", exist_ok = True)
    os.makedirs(f"{conf['save']}/log", exist_ok = True)

def check_tm_config(config):
    conf = config['tm']
    os.makedirs(conf['save'], exist_ok = True)
    os.makedirs(f"{conf['save']}/structure", exist_ok = True)
    os.makedirs(f"{conf['save']}/log", exist_ok = True)

def update_config_dirs(config):
    prefix = config['prefix']
    config['cwd'] = cwd = f"./{prefix}"

    os.makedirs(cwd, exist_ok=True)
    tasks = ['pureFe', 'carbon', 'tm']

    for task in tasks:
        if (save_path := config[task].get('save', None)) is not None:
            config[task]['save'] = f"{cwd}/{save_path}"
    return config

def parse_config(config, argv: list[str] | None=None):
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

    config['cwd'] = os.path.join(os.path.abspath(os.getcwd()), config['cwd'])
    return config
