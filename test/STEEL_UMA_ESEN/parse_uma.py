"""
Modified based on Jinmu Yu's code
"""
import os

from fairchem.core.units.mlip_unit import MLIPPredictUnit
from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
from fairchem.core import FAIRChemCalculator
import warnings

root_dir = '/data2/shared_data/cps'

UMA_MODELS = {
    'uma_sm': './UMA/uma_sm.pt',
    'uma-s-1p1': './UMA/uma-s-1p1.pt',
    'uma-m-1p1': './UMA/uma-m-1p1.pt',
    }

UMA_MODALS = {
    'omat': 'omat', # 'PBE',
    'oc20': 'oc20', # 'RPBE',
    'omol': 'omol', # None,  # wB97M-V
    'odac': 'odac', # None,  # PBE-D3
    'omc':  'omc' # None,  # PBE-D3
    }

UMA_FUNCTIONALS = {
    'omat': 'PBE',
    'oc20': 'RPBE',
    'omol': 'wB97M-V', # None
    'odac': 'PBE', # 'PBE-D3' # None
    'omc':  'PBE', # 'PBE-D3' # None
    }
  
def load_uma_calc(modal, dispersion=False, functional=None, model='uma-m-1p1', device='cuda'):

    model = f'{root_dir}/{UMA_MODELS[model.lower()]}'
    assert modal in UMA_MODALS.keys()
    assert os.path.isfile(model)
    modal = UMA_MODALS[modal.lower()]

    mlip_unit_kwargs = {
        'inference_model_path': model,
        'device': device,
        'inference_settings': InferenceSettings(
            external_graph_gen=False,
            ),
        }

    mlip_predict_unit = MLIPPredictUnit(**mlip_unit_kwargs)
    calc_uma = FAIRChemCalculator(predict_unit=mlip_predict_unit, task_name=modal)

    if modal in ['omc']:
        assert dispersion is True

    if dispersion:
        functional = UMA_FUNCTIONALS[modal]
        from ase.calculators.mixing import MixedCalculator
        from sevenn.calculator import D3Calculator
        calc_d3 = D3Calculator(functional_name = functional)
        calc = MixedCalculator(calc_uma, calc_d3, +1, -1)
        warnings.warn('Excluding D3 contribution')
        return calc
    else:
        return calc_uma


