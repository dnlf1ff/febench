from mace.calculators import MACECalculator

HEAD = '/data2/shared_data/cps/MACE'

MAP_DCT = {
        'mpa': 'mace-mpa-0-medium',
        'omat': 'mace-omat-0-medium'
        }

CALC_DCT = {
    'mace-mpa-0-medium': f'{HEAD}/mace-mpa-0-medium/mace-mpa-0-medium.model',
    'mace-omat-0-medium': f'{HEAD}/mace-omat-0-medium/mace-omat-0-medium.model',
}

FUNC_DCT = {
    'mace-mpa-0-medium': 'PBE',
    'mace-omat-0-medium': 'PBE',
}


def return_calc(config, dispersion=None):
    conf = config['calculator']
    model= modal = conf['modal']
    model_name = MAP_DCT[model]
    model_path = CALC_DCT[model_name]

    calc_kwargs = {
            'model_paths': model_path,
            'device': 'cuda',
        }
    functional_name = FUNC_DCT[model_name]

    print(f"[MACE] model={model_name}")
    print(f"[MACE] potential directory = {model_path}")
    functional_name = functional_name if dispersion else None
    calc_mlp = MACECalculator(**calc_kwargs)
    return calc_mlp
