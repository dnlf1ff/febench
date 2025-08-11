import os
import warnings

def check_calc_config(config):
    config_calc = config['calculator']
    calc_types = ['sevennet', '7net', '7net-mf', 'sevennet-mf', 'dpa', 'esen', 'uma', 'orb']
    model = config_calc['model'].lower()
    modal = config_calc['modal'].lower()
    dirname = config_calc['dirname']

    assert config_calc['calc_type'].lower() in calc_types
    assert os.path.exists(dirname), f'directory for potential file {dirname} does not exists'

    if config_calc['calc_type'].lower() in ['sevennet-mf', '7net-mf']:
        assert os.path.isfile(potential := f'{dirname}/{model}'), f'MLIP potential not found at {potential}'

        config['calculator']['calc_args']['model'] = potential
        assert modal in ['mpa', 'omat24'], 'SevenNet-MF potentials require "modal" as a keyword'
        config['calculator']['calc_args']['modal'] = modal

    elif config_calc['calc_type'].lower() in ['sevennet', '7net']:
        assert os.path.isfile(potential := f'{dirname}/{model}'), f'MLIP potential not found at {potential}'
        config['calculator']['calc_args']['model'] = potential
        if modal in ['mpa', 'omat24']:
            warnings.warn('SevenNet-0 and SevenNet-omat does not take "modal" as a keyword ...')
            config['calculator']['calc_args']['modal'] = None

    elif config_calc['calc_type'].lower() == 'dpa':
        from febench.util.parse_calc import DPA_MODELS, DPA_MODALS 
        assert model in DPA_MODELS.keys(),f'unknown DPA model {model}'
        model = DPA_MODELS[model]
        assert os.path.isfile(potential := f'{dirname}/{model}'), f'MLIP potential not found at {potential}'
        config['calculator']['calc_args']['model'] = potential
        assert modal in DPA_MODALS.keys(), f'unknown DPA modal {modal}'
        config['calculator']['calc_args']['modal'] = DPA_MODALS[modal]

    elif config_calc['calc_type'].lower() == 'esen':
        from febench.util.parse_calc import ESEN_MODELS
        assert model in ESEN_MODELS.keys(), f'unknown eSEN model {model}'
        model = ESEN_MODELS[model]
        assert os.path.isfile(potential := f'{dirname}/{model}'), f'MLIP potential not found at {potential}'
        config['calculator']['calc_args']['model'] = potential
 
    elif config_calc['calc_type'].lower() == 'orb':
        from febench.util.parse_calc import ORB_MODELS
        assert model in ORB_MODELS.keys(),f'unknown ORB model {model}'
        model = ORB_MODELS[model]
        assert os.path.isfile(potential := f'{dirname}/{model}'), f'MLIP potential not found at {potential}'
        config['calculator']['calc_args']['model'] = potential

    elif config_calc['calc_type'].lower() == 'uma':
        from febench.util.parse_calc import UMA_MODELS, UMA_MODALS
        if potential.endswith('pt'):
            assert model.lower() in UMA_MODELS.values(), f'unknown UMA model {model}'
            potential = f"{dirname}/../{model}"
        else:
            assert model.lower() in UMA_MODELS.keys(), f'unknown UMA model {model}'
            potential = f"{dirname}/{model}"

        assert modal in UMA_MODALS, f'unknown UMA modal {modal}'
        assert os.path.isfile(potential), f'MLIP potential not found at {potential}'
 
        config['calculator']['calc_args']['model'] = potential
        config['calculator']['calc_args']['modal'] = config_calc['modal']

    else:
        print('smting is wrong ...')
    

    return config

