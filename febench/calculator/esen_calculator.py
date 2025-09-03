"""
Modified based on Jinmu Yu's code
"""

from fairchem.core import OCPCalculator

HEAD = '/data2/shared_data/cps'

CALC_DCT = {
    'oam': f"{HEAD}/eSEN/esen_30m_oam.pt",
    'omat': f"{HEAD}/eSEN/esen_30m_omat.pt",
    }
 

def return_calc(config):
    conf = config['calculator']
    model, modal = conf['model'], conf['modal']
    model_path = CALC_DCT[modal]

    calc_kwargs = {
            'checkpoint_path': model_path,
            'cpu': False
            }

    print(f"[eSEN] model={model}, modal={modal}")
    print(f"[eSEN] potential path: {model_path}")

    calc = OCPCalculator(**calc_kwargs)
    return calc

