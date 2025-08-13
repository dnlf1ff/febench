from ase.io import write, read
from ase.filters import UnitCellFilter
from ase.optimize import FIRE
from ase.constraints import FixAtoms
import gc, torch, os
import numpy as np

def relax_pureFe(calc, debug=False):
    save_dir = './pureFe'
    log_dir = './pureFe/log'
    struct_dir = './pureFe/structure'
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(struct_dir, exist_ok=True)

    if debug:
        print(os.path.abspath(save_dir))
        print(os.path.abspath(log_dir))
        print(os.path.abspath(struct_dir))

    print('Iniciating relaxation of BCC Fe with one vacancy')
    atoms=read(f'{struct_dir}/POSCAR_Vac', format='vasp')
    if not debug:
        cell_filter = UnitCellFilter(atoms, mask=[1,1,1,0,0,0])
        fire = FIRE(cell_filter, logfile=f'{log_dir}/Vac_relax.log')
        fire.run(fmax=0.0001, steps=100000)
        write(f'{struct_dir}/CONTCAR_Vac', atoms, format='vasp')
        del atoms, cell_filter, fire
        gc.collect()
    print('relaxation of BCC Fe with one vacancy - completed\n')

    for hkl in ['100', '110', '111']:
        print(f'Iniciating relaxation of {hkl} surface structure of BCC Fe\n')
        atoms = read(f'{struct_dir}/POSCAR_surface_{hkl}', format='vasp')
        z =atoms.positions[:,2].copy()
        indices = [atom.index for atom in atoms if atom.position[2] < np.percentile(z,62) and atom.position[2]>np.percentile(z,38)]
        atoms.set_constraint(FixAtoms(indices=indices))
        atoms.set_pbc([True, True, True])

        if not debug:
            fire = FIRE(atoms, logfile=f'{log_dir}/surface_{hkl}_relax.log')
            fire.run(fmax=0.0001, steps=100000)
            write(f'{struct_dir}/CONTCAR_surface_{hkl}', atoms, format='vasp')
            del atoms, fire
            gc.collect()
        print(f'relaxation of {hkl} surface structure of BCC Fe - completed\n')


def relax_carbons(calc, debug=False):
    print('Iniciating relaxation of BCC Fe with Carbons')
    save_dir = './carbon'
    struct_dir = './carbon/structure'
    log_dir = './carbon/log'
    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(struct_dir, exist_ok=True)

    if debug:
        print(os.path.abspath(save_dir))
        print(os.path.abspath(log_dir))
        print(os.path.abspath(struct_dir))


    atoms = read(f'{struct_dir}/POSCAR_C', format='vasp')
    if not debug:
        cell_filter = UnitCellFilter(atoms, mask=[1,1,1,0,0,0])
        fire = FIRE(cell_filter, logfile=f'{log_dir}/FeC_relax.log')
        fire.run(fmax=0.0001, steps=100000)
        write(f'{struct_dir}/CONTCAR_FeC', atoms, format='vasp')
        del atoms, cell_filter, fire
        gc.collect()

    labels = ['a1', 'a2', 'b1', 'b2', 'b3', 'c1', 'c2', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']
    for idx, label in enumerate(labels):
        atoms = read(f'{struct_dir}/POSCAR_{label}', format='vasp')
        if not debug:
            cell_filter = UnitCellFilter(atoms, mask=[1,1,1,0,0,0])
            fire = FIRE(cell_filter, logfile=f'{log_dir}/{label}_relax.log')
            fire.run(fmax=0.0001, steps=100000)
            write(f'{struct_dir}/CONTCAR_{label}', atoms, format='vasp')
            del atoms, fire, cell_filter
            gc.collect()
        print(f'relaxation of BCC Fe with Carbon config {label} - completed\n')
    print('Relaxation of BCC Fe with Carbons - completed')

def relax_tms(calc, debug=False):
    print('Iniciating relaxation of BCC Fe with transition metals\n')
    save_dir = './tm'
    struct_dir = './tm/structure'
    log_dir = './tm/log'

    os.makedirs(save_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    os.makedirs(struct_dir, exist_ok=True)

    if debug:
        print(os.path.abspath(save_dir))
        print(os.path.abspath(log_dir))
        print(os.path.abspath(struct_dir))


    solutes = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au']
    n_solute = len(solutes)

    for j, sol in enumerate(solutes):
        print(f'Iniciating relaxation of BCC Fe with solute {sol} {j+1}/{n_solute}')
        atoms = read(f'{struct_dir}/POSCAR_{sol}', format='vasp')
        if not debug:
            cell_filter = UnitCellFilter(atoms, mask=[1,1,1,0,0,0])
            fire = FIRE(cell_filter, logfile=f'{log_dir}/{sol}_relax.log')
            fire.run(fmax=0.001, steps=100000)
            write(f'{struct_dir}/CONTCAR_{sol}', atoms, format='vasp')
            del atoms, cell_filter, fire
            gc.collect()

        for i in range(5): 
            print(f'\nIniciating relaxation of BCC Fe with {sol}-{sol} {i+1}nn/5nn')
            atoms = read(f'{struct_dir}/POSCAR_{sol}_{sol}_{i+1}nn', format='vasp')
            if not debug:
                cell_filter = UnitCellFilter(atoms, mask=[1,1,1,0,0,0])
                fire = FIRE(cell_filter, logfile=f'{log_dir}/{sol}_{sol}_{i+1}nn_relax.log')
                fire.run(fmax=0.001, steps=100000)
                write(f'{struct_dir}/CONTCAR_{sol}_{sol}_{i+1}nn', atoms, format='vasp')
                del atoms, cell_filter, fire
                gc.collect()
 
            print(f'Iniciating relaxation of BCC Fe with {sol}-Vac {i+1}nn/5nn')
            atoms = read(f'{struct_dir}/POSCAR_{sol}_Vac_{i+1}nn', format='vasp')
            if not debug:
                cell_filter = UnitCellFilter(atoms, mask=[1,1,1,0,0,0])
                fire = FIRE(cell_filter, logfile=f'{log_dir}/{sol}_Vac_{i+1}nn_relax.log')
                fire.run(fmax=0.001, steps=100000)
                write(f'{struct_dir}/CONTCAR_{sol}_Vac_{i+1}nn', atoms, format='vasp')
                del atoms, cell_filter, fire
                gc.collect()
 
        print(f'Relaxation of BCC Fe with solute {sol} {j+1}/{n_solute} -- completed\n')
    print('Relaxation of BCC Fe with transition metals - completed \n')
