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
 
def save_carbon(fmax, cell):
    prefix = f'{run_params["date"]}_fmax_{fmax}_{cell}'
    df = pd.DataFrame(index=carbon_config)
    df['dft'] = carbon_dft
    for model, dct in model_dct.items():
        for modal in dct['modal']:
            try:
                df_carbon = pd.read_csv(f'{head}/{prefix}/{model}/{modal}/carbon/carbon.csv', comment='#')
                df[f'{model}-{modal}']  = df_carbon['config'].tolist() #TODO
            except:
                print(f'ERROR:{model}-{modal} {fmax} {cell}')
    df.to_csv(f'carbon_aug_fmax_{fmax}_{cell}.csv')

def save_tm(fmax, cell):
    prefix = f'{run_params["date"]}_fmax_{fmax}_{cell}'
    df = pd.DataFrame(index=tm_config)
    df['dft'] = tm_dft
    for model, dct in model_dct.items():
        for modal in dct['modal']:
            try:
                df_tm = pd.read_csv(f'{head}/{prefix}/{model}/{modal}/tm/tm_E_bind.csv', ) # comment=#
                df[f'{model}-{modal}'] = df_tm['E_Fe'].tolist() # need to fix this
            except Exception as exc:
                print(f'ERROR:{model}-{modal} {fmax} {cell}')
    df.to_csv(f'tm_aug_fmax_{fmax}_{cell}.csv')


def get_carbon_errors(fmax, cell):
    suffix = f'aug_fmax_{fmax}_{cell}'
    carbon_df = pd.read_csv(f'carbon_{suffix}.csv', index_col=0)
    carbon_true = carbon_df['dft']
    print(carbon_true)
    carbon_error = pd.DataFrame(index=carbon_config)
    for col in carbon_df.columns:
        print(carbon_df[col])
        carbon_error[col] = carbon_df[col] - carbon_true

    carbon_error.to_csv(f'carbon_error_fmax_{fmax}_{cell}.csv')
    
def get_tm_errors(fmax, cell):
    suffix = f'aug_fmax_{fmax}_{cell}'
    tm_df = pd.read_csv(f'tm_{suffix}.csv', index_col=0)
    tm_true = tm_df['dft']
    print(tm_true)
    tm_error = pd.DataFrame(index=tm_config, index_col=0)
    for col in tm_df.columns:
        print(tm_df[col])
        tm_error[col] = tm_df[col] - tm_true

    tm_error.to_csv(f'tm_error_fmax_{fmax}_{cell}.csv')
 

def get_abs(fmax, cell):
    prefix = f'250916_fmax_{fmax}_{cell}'
    carbon_error = pd.read_csv(f'{prefix}/carbon_error_fmax_{fmax}_{cell}.csv', index_col=0)
    for col in carbon_error.columns:
        carbon_error[col] = carbon_error[col].apply(lambda x: abs(x))
    carbon_error.to_csv(f'{prefix}/carbon_abs_fmax_{fmax}_{cell}.csv')

    tm_error = pd.read_csv(f'{prefix}/tm_error_fmax_{fmax}_{cell}.csv', index_col=0)
    for col in tm_error.columns:
        tm_error[col] = tm_error[col].apply(lambda x: abs(x))
    tm_error.to_csv(f'{prefix}/tm_abs_fmax_{fmax}_{cell}.csv')

def get_mae(fmax, cell):
    suffix = f'abs_fmax_{fmax}_{cell}'
    prefix = f'250916_fmax_{fmax}_{cell}'
    carbon_df = pd.read_csv(f'{prefix}/carbon_{suffix}.csv', index_col=0)
    tm_df = pd.read_csv(f'{prefix}/tm_{suffix}.csv', index_col=0)
    carbon_mae = []
    tm_mae = []
    for carbon_col in carbon_df.columns:
        carbon_mae.append(carbon_df[carbon_col].mean())
    for tm_col in tm_df.columns:
        tm_mae.append(tm_df[tm_col].mean())

    return carbon_mae, tm_mae


if __name__ == '__main__':
    carbon_df = pd.DataFrame()
    tm_df = pd.DataFrame()
    mlips = get_columns()
    carbon_df['mlip'] = mlips
    tm_df['mlip'] = mlips
    for fmax in run_params['fmax']:
        for cell in run_params['cell']:
            get_abs(fmax, cell)
            label = f'fmax_{fmax}_{cell}'
            carbon_mae_list, tm_mae_list = get_mae(fmax, cell)
            carbon_df[label] = carbon_mae_list
            tm_df[label] = tm_mae_list

    carbon_df.to_csv('carbon_mae.csv', index=False)
    tm_df.to_csv('tm_mae.csv', index=False)
