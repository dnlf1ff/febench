"""
Modified based on Jinmu Yu's code
"""

from types import NotImplementedType
import warnings

from ase.calculators.mixing import MixedCalculator
from orb_models.forcefield import pretrained
from orb_models.forcefield.calculator import ORBCalculator
 

CALC_DCT = {
    "omat": "orb_v3_conservative_inf_omat",
    "mpa": "orb_v3_conservative_inf_mpa",
    }


def return_calc(config, dispersion=None):
    conf = config['calculator']
    model, modal = conf['model'], conf['modal']
    model_name = CALC_DCT[modal]
    device = 'cuda'

    print(f"[ORB] model={model_name}")
 
    orbff = getattr(pretrained, model_name)(
            precision="float32-highest",
            device=device,
        )

    calc = ORBCalculator(orbff, device=device)
    return calc
 
