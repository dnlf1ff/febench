"""
Modified based on Jinmu Yu's code
"""

from types import NotImplementedType
import warnings

from fairchem.core.units.mlip_unit import MLIPPredictUnit
from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
from fairchem.core.calculate.pretrained_mlip import get_isolated_atomic_energies
from fairchem.core import FAIRChemCalculator

from ase.calculators.mixing import MixedCalculator

HEAD = '/data2/shared_data/cps'

CALC_DCT = {
    'uma-m-1p1': f'{HEAD}/UMA/uma-m-1p1.pt',
    }

MODAL_DCT = {
    'omat': 'omat', # 'PBE',
    'omc':  'omc' # None,  # PBE-D3
    }

FUNC_DCT = {
    'omat': 'PBE',
    'omc':  'pbe', # 'PBE-D3' # None
    }
  
# https://www.aissquare.com/models/detail?pageType=models&name=DPA-2.3.1-v3.0.0rc0&id=287#data-used-for-pretraining
 
def return_calc(config):
    conf = config['calculator']
    model = conf['model']
    modal = task = conf['modal']

    model_name = list(CALC_DCT.keys())[0]
    model_path =  CALC_DCT[model_name]

    print(f"[UMA] model={model_name}, modal(task_name)={modal}")
    print(f"[UMA] potential path: {model_path}")

    mlip_unit_kwargs = {
        'inference_model_path': model_path,
        'device': 'cuda',
        'inference_settings': InferenceSettings(
            external_graph_gen=False,
            ),
        'atom_refs': get_isolated_atomic_energies(model_name=model_name),
        }

    mlip_predict_unit = MLIPPredictUnit(**mlip_unit_kwargs)
    calc_uma = FAIRChemCalculator(
            predict_unit=mlip_predict_unit,
            task_name=MODAL_DCT[modal]
            )

    if modal == 'omc':
        from sevenn.calculator import D3Calculator
        functional = FUNC_DCT[modal]
        calc_d3 = D3Calculator(functional_name = functional)
        calc_mix = MixedCalculator(calc_uma, calc_d3, +1, -1)
        print('WARNING: EXCLUDING D3 CONTRIBUTION')
        print(f"[UMA] functional: {functional}")
        print('WARNING: CALCULATING WITH UMA-omc WITH D3 CONTRIBUTION')
        return calc_mix

    else:
        return calc_uma

if __name__  == '__main__':
    import sys
    model, modal = sys.argv[1], sys.argv[2]
    dispersion = False
    functional = 'PBE'

    if modal == 'omc':
        dispersion = True

    config = {'calculator': {'model': model, 'modal': modal, 'dispersion': dispersion, 'functional': functional}}

    return_calc(config)
