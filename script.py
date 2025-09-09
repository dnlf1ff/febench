import pandas as pd
import numpy as np
HEAD='/Users/jinvk/250904_clsd'
TAIL='/Users/jinvk/febench'
MLIP=['ompa', 'mace', 'orb', 'esen', 'dpa', 'uma', 'omni']
MLIP_CONFIG={
        'ompa': {'label': 'ompa', 'modal': ['mpa', 'omat24'],},
        'mace': {'label': 'MACE', 'modal': ['mpa', 'omat'],},
        'orb': {'label': 'ORB', 'modal': ['mpa', 'omat'],},
        'esen': {'label': 'eSEN', 'modal': ['omat', 'oam'],},
        'dpa': {'label': 'DPA', 'modal': ['mp', 'omat'],},
        'uma': {'label': 'UMA', 'modal': ['omc', 'omat'],},
        'omni': {'label': 'omni', 'modal': ['mpa', 'omat24', 'matpes_pbe'],},
        }
MLIP_LIST=[]
for mlip in MLIP:
    for modal in MLIP_CONFIG[mlip]['modal']:
        MLIP_LIST.append(f"{mlip}-{modal}")

CARBON=['a1', 'a2', 'b1', 'b2', 'b3', 'c1', 'c2', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']

TM=['Co', 'Cr', 'Cu', 'Mn', 'Mo', 'Nb', 'Ni', 'Ti', 'V']

CONFIG = ['Fe128', 'Fe127Vac', 'Fe128C', 'a1', 'a2', 'b1', 'b2', 'b3', 'c1', 'c2', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'Fe127Co', 'Fe126Co2', 'Fe127Cr', 'Fe126Cr2', 'Fe127Cu', 'Fe126Cu2', 'Fe127Mn', 'Fe126Mn2', 'Fe127Mo', 'Fe126Mo2', 'Fe127Nb', 'Fe126Nb2', 'Fe127Ni', 'Fe126Ni2','Fe127Ti', 'Fe126Ti2', 'Fe127V', 'Fe126V2']

df = pd.DataFrame(index=CONFIG)
from ase.io import read
def get_e_fr_energy():
    for mlip in MLIP:
        for modal in MLIP_CONFIG[mlip]['modal']:
            label = MLIP_CONFIG[mlip]['label']
            mlip_path = f"{HEAD}/{mlip}/{modal}"
            energy_list = []
            energy_list.append(read(f"{mlip_path}/pureFe/structure/bulk_opt.extxyz").info['e_fr_energy'])
            energy_list.append(read(f"{mlip_path}/pureFe/structure/Vac_opt.extxyz").info['e_fr_energy'])
            energy_list.append(read(f"{mlip_path}/carbon/structure/FeC.extxyz").info['e_fr_energy'])
            carbon_df = pd.read_csv(f"{mlip_path}/carbon/carbon.csv", comment='#')
            energy_list.extend(carbon_df['E_FeCVac'])
            for tm in TM:
                energy_list.append(read(f"{mlip_path}/tm/structure/{tm}_opt.extxyz").info['e_fr_energy'])
                logfile = open(f"{mlip_path}/tm/log/{tm}_{tm}_1nn_relax.log", 'r')
                line = logfile.readlines()[-1]
                energy_list.append(float(line.split()[-2]))
                logfile.close()
            df[f'{mlip}-{modal}'] = energy_list

    df.to_csv(f"{HEAD}/e_fr_energy.csv")

def parse_logfile(logfile):
    logfile = open(logfile, 'r')
    line = logfile.readlines() [-1]
    step = int(line.split()[1])
    energy = float(line.split()[-2])
    fmax = float(line.split()[-1])
    logfile.close()
    return None

def get_carbon():
    df = pd.DataFrame(index=CARBON)
    for mlip in MLIP:
        for modal in MLIP_CONFIG[mlip]['modal']:
            label = MLIP_CONFIG[mlip]['label']
            mlip_path = f"{HEAD}/{mlip}/{modal}"
            carbon_df = pd.read_csv(f"{mlip_path}/carbon/carbon.csv", comment='#')
            df[f'{mlip}-{modal}'] = carbon_df['E_bind'].values
    df.to_csv(f"{TAIL}/carbon.csv")


def get_opt_info():
    for mlip in MLIP:
        for modal in MLIP_CONFIG[mlip]['modal']:
            label = MLIP_CONFIG[mlip]['label']
            mlip_path = f"{HEAD}/{mlip}/{modal}"
    pass
 
if __name__ == "__main__":
#    get_e_fr_energy()
    get_carbon()
