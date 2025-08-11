"""
Modified based on Jinmu Yu's code
"""

from types import NotImplementedType

DPA_MODELS ={
    # Matbench
    'dpa31-mp': f'./DPA3/dpa-3.1-mptrj.pth',
    'dpa-mp': f'./DPA3/dpa-3.1-mptrj.pth',
    'dpa31_mp': f'./DPA3/dpa-3.1-mptrj.pth',
    'dpa_mp': f'./DPA3/dpa-3.1-mptrj.pth',
    'mp': f'./DPA3/dpa-3.1-mptrj.pth',
    'dpa31-ft': f'./DPA3/dpa-3.1-3m-ft.pth', 
    'dpa31_ft': f'./DPA3/dpa-3.1-3m-ft.pth', 
    'ft': f'./DPA3/dpa-3.1-3m-ft.pth', 
    # https://www.aissquare.com/models/detail?pageType=models&name=DPA-3.1-3M&id=343
    'dpa31-openlam': f'./DPA3/DPA-3.1-3M.pt',
    'dpa-openlam': f'./DPA3/DPA-3.1-3M.pt',
    'dpa31_openlam': f'./DPA3/DPA-3.1-3M.pt',
    'dpa_openlam': f'./DPA3/DPA-3.1-3M.pt',
    'openlam': f'./DPA3/DPA-3.1-3M.pt',
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

ESEN_MODELS = {
    'esen'      : "./eSEN/esen_30m_oam.pt",
    'esen_oam'      : "./eSEN/esen_30m_oam.pt",
    'esen-oam'      : "./eSEN/esen_30m_oam.pt",
    'oam'      : "./eSEN/esen_30m_oam.pt",
    'esen_mp'   : "./eSEN/esen_30m_mptrj.pt",
    'esen-mp'   : "./eSEN/esen_30m_mptrj.pt",
    'esen_mpa'   : "./eSEN/esen_30m_mptrj.pt",
    'esen-mpa'   : "./eSEN/esen_30m_mptrj.pt",
    'mpa'   : "./eSEN/esen_30m_mptrj.pt",
    'mp'   : "./eSEN/esen_30m_mptrj.pt",
    'esen_omat' : "./eSEN/esen_30m_omat.pt",
    'esen_omat24' : "./eSEN/esen_30m_omat.pt",
    'esen-omat' : "./eSEN/esen_30m_omat.pt",
    'esen-omat24' : "./eSEN/esen_30m_omat.pt",
    'omat' : "./eSEN/esen_30m_omat.pt",
    'omat24' : "./eSEN/esen_30m_omat.pt",
    'eqV2_31M_oms'  : "./eqV2/eqV2_31M_omat_mp_salex.pt",
    'eqV2_153M_oms' : "./eqV2/eqV2_153M_omat_mp_salex.pt",
    'eqV2_dens_31M_mp'  : "./eqV2/eqV2_dens_31M_mp.pt",
    'eqV2_dens_153M_mp' : "./eqV2/eqV2_dens_153M_mp.pt",
    }
 

ORB_MODELS = {
    "orb_omat": "orb_v3_conservative_inf_omat",
    "orb_omat24": "orb_v3_conservative_inf_omat",
    "orb-omat": "orb_v3_conservative_inf_omat",
    "orb-omat24": "orb_v3_conservative_inf_omat",
    "omat": "orb_v3_conservative_inf_omat",
    "omat24": "orb_v3_conservative_inf_omat",
    "mpa": "orb_v3_conservative_inf_omat",
    "mp": "orb_v3_conservative_inf_omat",
    "orb-mpa": "orb_v3_conservative_inf_mpa",
    "orb-mpa": "orb_v3_conservative_inf_mpa",
    }


UMA_MODELS = {
    'uma_sm': './UMA/uma_sm.pt',
    'uma-s-1p1': './UMA/uma-s-1p1.pt',
    'uma-m-1p1': './UMA/uma-m-1p1.pt',
    }

UMA_MODALS = {
    'omat': 'omat', # 'PBE',
    'omat24': 'omat', # 'PBE',
    'oc20': 'oc20', # 'RPBE',
    'omol': 'omol', # None,  # wB97M-V
    'odac': 'odac', # None,  # PBE-D3
    'omc':  'omc' # None,  # PBE-D3
    }
  
# https://www.aissquare.com/models/detail?pageType=models&name=DPA-2.3.1-v3.0.0rc0&id=287#data-used-for-pretraining
def load_dpa_calc(config):
    from deepmd.calculator import DP
    calc_args = config['calculator']['calc_args']

    calc_kwargs = {
            'model': calc_args['model'],
            'device': calc_args['device'],
            'head': calc_args['modal']
            }

    if calc_args['model'].endswith('pt'):
        return DP(**calc_kwargs)

    elif calc_args['model'].endswith('pth'):
        calc_kwargs.pop('head', None)
        return DP(**calc_kwargs)

def load_esen_calc(config):
    from fairchem.core import OCPCalculator
    calc_args = config['calculator']['calc_args']
    calc_kwargs = {
            'checkpoint_path': calc_args['model'],
            'cpu': False
            }

    return OCPCalculator(**calc_kwargs)

def load_orb_calc(config):
    from orb_models.forcefield import pretrained
    from orb_models.forcefield.calculator import ORBCalculator
    calc_args = config['calculator']['calc_args']

    orbff = getattr(pretrained, ORB_MODELS[calc_args['modal'].lower()])(
            precision="float32-highest",
            device=calc_args['device'],
        )

    return ORBCalculator(orbff, device=calc_args['device'])
 
def load_uma_calc(config):
    from fairchem.core.units.mlip_unit import MLIPPredictUnit
    from fairchem.core.units.mlip_unit.api.inference import InferenceSettings
    from fairchem.core import FAIRChemCalculator

    calc_args = config['calculator']['calc_args']
    mlip_unit_kwargs = {
        'inference_model_path': calc_args["model"],
        'device': calc_args['device'],
        'inference_settings': InferenceSettings(
            external_graph_gen=False,
            ),
        }

    mlip_predict_unit = MLIPPredictUnit(**mlip_unit_kwargs)
    return FAIRChemCalculator(predict_unit=mlip_predict_unit, 
                              task_name=UMA_MODALS[calc_args['modal'].lower()])

def load_calc(config):
    calc_type = config['calculator']['calc_type']

    if calc_type.lower() == 'dpa':
        calc = load_dpa_calc(config)

    elif calc_type.lower() == 'esen':
        calc = load_esen_calc(config)

    elif calc_type.lower() == 'orb':
        calc = load_orb_calc(config)

    elif calc_type.lower() == 'uma':
        calc = load_uma_calc(config)

    else:
        calc = None
        raise NotImplementedError
    return calc
