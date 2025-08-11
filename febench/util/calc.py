import numpy as np
from tqdm import tqdm
from sevenn.calculator import SevenNetCalculator
from ase.calculators.singlepoint import SinglePointCalculator

"""
modified based on Jaesun Kim's code
"""

def calc_from_py(config):
    import importlib.util
    from pathlib import Path

    file_path = Path(__file__).resolve().parent
    spec = importlib.util.spec_from_file_location(f'load_calc', f'{file_path}/parse_calc.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    calc = module.load_calc(config)
    # print(calc)
    return calc


def calc_from_config(config):
    calc_config = config['calculator']
    calc_type = calc_config['calc_type'].lower()

    if calc_type in ['sevennet', 'sevennet-mf', '7net', '7net-mf']:
        calc_kwargs = {'model': calc_config['calc_args']['model'],
                   'modal': calc_config['calc_args']['modal'],
                   'device': calc_config['calc_args']['device']
                   }
        if calc_type in ['sevennet-mf', '7net-mf']:
            return SevenNetCalculator(**calc_kwargs)

        else:
            calc_kwargs.pop('modal', None)
            return SevenNetCalculator(**calc_kwargs)


    else:
        return calc_from_py(config)


def single_point_calculate(atoms, calc):
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    stress = atoms.get_stress()

    calc_results = {"energy": energy, "forces": forces, "stress": stress}
    calculator = SinglePointCalculator(atoms, **calc_results)
    new_atoms = calculator.get_atoms()

    return new_atoms


def single_point_calculate_list(atoms_list, calc, desc=None):
    calculated = []
    for atoms in tqdm(atoms_list, desc=desc, leave=False):
        calculated.append(single_point_calculate(atoms, calc))

    return calculated
