
from types import NotImplementedType
import warnings
from sevenn.calculator import SevenNetCalculator

HEAD = '/data2/shared_data'

model_path = f'{HEAD}/pretrained/7net_chgTot/checkpoint_best.pth'

def return_calc(config, dispersion=None):
    calc_kwargs = {
        'model': model_path,
        'device': 'cuda',
    }

    print(f"[TEST] model=7net-0")
    print(f"[SevenNet] potential path: {model_path}")

    calc = SevenNetCalculator(**calc_kwargs)

    return calc
