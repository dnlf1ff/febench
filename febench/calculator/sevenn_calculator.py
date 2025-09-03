"""
Modified based on Jinmu Yu's code
"""

from types import NotImplementedType
import warnings
from sevenn.calculator import SevenNetCalculator

CALC_DCT = {
    'ompa': '/data2/shared_data/cps/o50mp7a_ft/checkpoint_2.pth',
    'omni': '/data2/shared_data/cps/omni/oob_v6/checkpoint_1.pth',
    }

FUNC_DCT = {
    'mpa': 'PBE',
    'omat24': 'PBE',
    'matpes_pbe': 'PBE',
    'spice': 'wB97M',
    'qcml': 'PBE0',
    'oc20': 'RPBE',
    'oc22': 'PBE',
    'mp_r2scan': 'r2SCAN',
    'matpes_r2scan': 'r2SCAN',
}

def return_calc(config, dispersion=None):
    conf = config['calculator']
    model, modal = conf['model'], conf['modal']
    model_path = CALC_DCT[model]
    calc_kwargs = {
        'model': model_path,
        'modal': modal,
    }

    functional = FUNC_DCT.get(modal, None)
    print(f"[SevenNet] model={model}, modal={modal}")
    print(f"[SevenNet] potential path: {model_path}")

    calc = SevenNetCalculator(**calc_kwargs)
    functional_name = functional_name if dispersion else None

    return calc
