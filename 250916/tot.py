import pandas as pd
import sys
head = '/data2_1/jinvk/25_Fe'

model_dct = {
    'ompa': {'modal': ['mpa', 'omat24'], 'label': 'ompa'},
    'mace': {'modal': ['mpa', 'omat'], 'label': 'MACE'},
    'orb': {'modal': ['mpa', 'omat'], 'label': 'ORB'},
    'esen': {'modal': ['oam', 'omat'], 'label': 'eSEN'},
    'dpa': {'modal': ['mp', 'omat'], 'label': 'DPA'},
    'uma': {'modal': ['omat'], 'label': 'UMA'},
    'omni': {'modal': ['mpa', 'omat24', 'matpes_pbe'], 'label': 'omni'},
    'grace': {'modal': ['oam', 'omat'], 'label': 'GRACE'},
    'nequip': {'modal': ['oam'], 'label': 'NequIP'},
        }

run_params = {'date': '250916', 'fmax': ['1e-2', '0.5e-3'], 'cell': ['ucf', 'frech']}
carbon_dft = [0.47, -0.01, 0.68, 0.47, 0.34, 0.58, 0.12, 0.16, 0.13, -0.09, -0.09, -0.65, -1.67, 1.07, 1.5] 
tm_dft = [-0.043, -0.239, 0.246, 0.079, -0.277, -0.381, 0.02, -0.246, -0.233]
carbon_config = ['a1', 'a2', 'b1', 'b2', 'b3', 'c1', 'c2', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']
tm_config = ['Co', 'Cr', 'Cu', 'Mn', 'Mo', 'Nb', 'Ni', 'Ti', 'V']
import sys

def get_columns():
    columns = ['dft']
    for model, dct in model_dct.items():
        for modal in dct['modal']:
            label = f'{model}-{modal}'
            columns.append(label)
    return columns
 
def get_vals(fmax, cell):
    prefix = f'{run_params["date"]}_fmax_{fmax}_{cell}'
    df = pd.DataFrame(index=carbon_config+tm_config)
    df['dft'] = carbon_dft + tm_dft
    for model, dct in model_dct.items():
        for modal in dct['modal']:
            try:
                df_carbon = pd.read_csv(f'{head}/{prefix}/{model}/{modal}/carbon/carbon.csv', comment='#')
                df_tm = pd.read_csv(f'{head}/{prefix}/{model}/{modal}/tm/tm_E_bind.csv')
                df[f'{model}-{modal}']  = df_carbon['config'].tolist() + df_tm['E_Fe'].tolist() 
                # TODO check if this is necessary .. can't we just do reverse-calc?
            except:
                print(f'ERROR:{model}-{modal} {fmax} {cell}')
    df.to_csv(f'tot_aug_fmax_{fmax}_{cell}.csv')


def get_errors(fmax, cell):
    suffix = f'aug_fmax_{fmax}_{cell}'
    _df = pd.read_csv(f'tot_{suffix}.csv', index_col=0)
    _true = _df['dft']
    _error = pd.DataFrame(index=carbon_config+tm_config)
    for col in _df.columns:
        _error[col] = _df[col] - _true

    _error.to_csv(f'tot_error_fmax_{fmax}_{cell}.csv')
    
def get_abs(fmax, cell):
    _error = pd.read_csv(f'tot_error_fmax_{fmax}_{cell}.csv', index_col=0)
    for col in _error.columns:
        _error[col] = _error[col].apply(lambda x: abs(x))
    _error.to_csv(f'tot_abs_fmax_{fmax}_{cell}.csv')

def get_mae(fmax, cell):
    suffix = f'abs_fmax_{fmax}_{cell}'
    _df = pd.read_csv(f'tot_{suffix}.csv', index_col=0)
    _mae = []
    for _col in _df.columns:
        _mae.append(_df[_col].mean())
    return _mae

if __name__ == '__main__':
    _df = pd.DataFrame()
    mlips = get_columns()
    _df['mlip'] = mlips
    for fmax in run_params['fmax']:
        for cell in run_params['cell']:
            get_vals(fmax, cell)
            get_errors(fmax, cell)
            get_abs(fmax, cell)
            label = f'fmax_{fmax}_{cell}'
            _mae_list = get_mae(fmax, cell)
            _df[label] = _mae_list

    _df.to_csv('tot_mae.csv', index=False)
