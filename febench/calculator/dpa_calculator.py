"""
Modified based on Jinmu Yu's code
"""

from types import NotImplementedType
import warnings

from deepmd.calculator import DP

# https://www.aissquare.com/models/detail?pageType=models&name=DPA-2.3.1-v3.0.0rc0&id=287#data-used-for-pretraining
HEAD = f'/data2/shared_data/cps'
# https://www.aissquare.com/models/detail?pageType=models&name=DPA-3.1-3M&id=343
# 'dpa31-openlam': f'./DPA3/DPA-3.1-3M.pt',

model_name = 'dpa31-openlam'

model_path = f'{HEAD}/DPA3/DPA-3.1-3M.pt'

MODAL_DCT ={
    'mp': 'MP_traj_v024_alldata_mixu',
    'omat': 'Omat24',
    }

def return_calc(config):
    conf = config['calculator']
    model, modal = conf['model'], conf['modal']

    calc_kwargs = {
            'model': model_path,
            'device': 'cuda',
            'head': MODAL_DCT[modal], 
            }

    print(f"[DPA] model={model_name}, modal={MODAL_DCT[modal]}")
    print(f"[DPA] potential directory = {model_path}")

    calc = DP(**calc_kwargs)
    return calc
