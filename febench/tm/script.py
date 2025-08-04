from ase import Atoms
from ase.io import write, read
from ase.build import make_supercell 
from ase.lattice.cubic import BodyCenteredCubic
import yaml

from febench.util.parse_args import parse_base_args
from febench.util.utils import dumpYAML 
from febench.util.parse_config import parse_config_yaml 

from febench.util.relax import aar_from_config
from febench.tm.utils import write_poscar_from_config
import pandas as pd
import numpy as np
import gc

def fe_base(config):
    bulk_df = pd.read_csv(f'{config["pureFe"]["save"]}/bulk.csv')
    a = bulk_df['a'][1]/config["pureFe"]["bulk"]["supercell"][1]
    struct_dir = f'{config["tm"]["save"]}/structure'

    fe_unit = BodyCenteredCubic(directions=np.diag([1,1,1]), size=(1,1,1), symbol='Fe', pbc=True, latticeconstant=a)
    atoms = make_supercell(fe_unit, np.diag(config["tm"]["supercell"]))
    write(f'{struct_dir}/POSCAR_base', atoms, format='vasp')
    return a

def process_tm(config, calc):
    save_dir = config["tm"]["save"]
    struct_dir = f'{config["tm"]["save"]}/structure'
    log_dir = f'{config["tm"]["save"]}/log'

    a = fe_base(config)

    sols = config["tm"]["solute"]
    df = pd.DataFrame(columns = [['E_FeM']]+[f'E_FeMM_{int(i+1)}nn' for i in range(5)]+[f'E_FeMVac_{int(i+1)}nn' for i in range(5)], index = sols)

    for sol in sols:
        sol_list = []
        vac_list = []
        write_poscar_from_config(config, sol, a)

        # calc Fe(n-1)M
        atoms = read(f'{struct_dir}/POSCAR_{sol}', format='vasp')

        ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_relax.log')
        atoms, conv = ase_atom_relaxer.relax_atoms(atoms)
        atoms = ase_atom_relaxer.update_atoms(atoms)
        atoms.info['conv'] = conv
        atoms.calc = None
        write(f'{struct_dir}/CONTCAR_{sol}', atoms, format='vasp')
        write(f'{struct_dir}/{sol}_opt.extxyz', atoms, format='extxyz')

        del  ase_atom_relaxer
        gc.collect()

        for i in range(5): 
            # calc Fe(n-2)MVac
            atoms_vac = read(f'{struct_dir}/POSCAR_{sol}_vac_{int(i+1)}nn', format='vasp')

            ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_vac_{int(i+1)}nn.log')
            atoms_vac, conv = ase_atom_relaxer.relax_atoms(atoms_vac)
            atoms_vac = ase_atom_relaxer.update_atoms(atoms_vac)
            atoms_vac.info['conv'] = conv
            atoms_vac.calc = None
            write(f'{struct_dir}/CONTCAR_{sol}_vac_{int(i+1)}nn', atoms_vac, format='vasp')
            vac_list.append(atoms_vac)

            del  ase_atom_relaxer
            gc.collect()

            # calc Fe(n-2)M(2)
            atoms_sol = read(f'{struct_dir}/POSCAR_{sol}_{sol}_{int(i+1)}nn', format='vasp')
            ase_atom_relaxer = aar_from_config(config, calc,opt=config["tm"]["opt"], logfile = f'{log_dir}/{sol}_{sol}_{int(i+1)}nn.log')
            atoms_sol, conv = ase_atom_relaxer.relax_atoms(atoms_sol)
            atoms_sol = ase_atom_relaxer.update_atoms(atoms_sol)
            atoms_sol.info['conv'] = conv
            atoms_sol.calc = None
            write(f'{struct_dir}/CONTCAR_{sol}_{sol}_{int(i+1)}nn', atoms_sol, format='vasp')
            sol_list.append(atoms_sol)

            del  ase_atom_relaxer
            gc.collect()

        write(f'{struct_dir}/{sol}_{sol}.extxyz', sol_list)
        write(f'{struct_dir}/{sol}_vac.extxyz', vac_list)
 
        df.loc[sol] = [atoms.info['e_fr_energy']] +  [atom.info['e_fr_energy'] for atom in sol_list] + [atom.info['e_fr_energy'] for atom in vac_list]
        df.to_csv(f'{save_dir}/tm.csv', index_label = 'solute')
        del atoms, sol_list, vac_list
        gc.collect()


def main(argv: list[str] | None=None) -> None:
    from febench.util.calc import calc_from_config
    args = parse_base_args(argv)
    config_dir = args.config
    calc = args.calc
    modal = args.modal

    with open(config_dir, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    config['calculator']['prefix'] = calc
    config['calculator']['modal'] = modal
    config = parse_config_yaml(config)
    dumpYAML(config, f'{config["cwd"]}/config_tm.yaml')
    calc = calc_from_config(config)

    process_tm(config, calc)



if __name__ == '__main__':
    main()
