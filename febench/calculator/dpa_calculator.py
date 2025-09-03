"""
Modified based on Jinmu Yu's code
"""

from types import NotImplementedType
import warnings

from ase.calculators.mixing import MixedCalculator
from deepmd.calculator import DP

DPA_MODELS ={
    # Matbench
    # https://www.aissquare.com/models/detail?pageType=models&name=DPA-3.1-3M&id=343
    'dpa31-openlam': f'./DPA3/DPA-3.1-3M.pt',
    }

DPA_MODALS ={
    # is read as **calc_kwargs only for dpa31-openlam
    'mpa': 'MP_traj_v024_alldata_mixu',
    'mp': 'MP_traj_v024_alldata_mixu',
    'alex2d': 'Alex2D',
    'omat24': 'Omat24',
    'omat': 'Omat24',
    'oc22': 'OC22',
    }

# https://www.aissquare.com/models/detail?pageType=models&name=DPA-2.3.1-v3.0.0rc0&id=287#data-used-for-pretraining
def return_calc(config):
    conf = config['calculator']

    calc_kwargs = {
            'model': conf['model'],
            'device': conf['device'],
            'head': conf['modal']
            }

    if conf['model'].endswith('pt'):
        return DP(**calc_kwargs)

    elif conf['model'].endswith('pth'):
        calc_kwargs.pop('head', None)
        return DP(**calc_kwargs)
