"""
Modified based on Jinmu Yu's code
"""
import os

from fairchem.core import OCPCalculator
root_dir='/data2/shared_data/cps'

ESEN_MODELS = {
    'esen'      : "./eSEN/esen_30m_oam.pt",
    'esen_mp'   : "./eSEN/esen_30m_mptrj.pt",
    'esen_omat' : "./eSEN/esen_30m_omat.pt",
    'eqV2_31M_oms'  : "./eqV2/eqV2_31M_omat_mp_salex.pt",
    'eqV2_153M_oms' : "./eqV2/eqV2_153M_omat_mp_salex.pt",
    'eqV2_dens_31M_mp'  : "./eqV2/eqV2_dens_31M_mp.pt",
    'eqV2_dens_153M_mp' : "./eqV2/eqV2_dens_153M_mp.pt",
    }
 

def load_esen_calc(model):
    model = f'{root_dir}/{ESEN_MODELS[model]}'
    assert os.path.isfile(model)
    calc_kwargs = {
            'checkpoint_path': model,
            'cpu': False
            }

    return OCPCalculator(**calc_kwargs)

