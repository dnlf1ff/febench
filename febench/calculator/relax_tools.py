import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from ase.atoms import Atoms
from ase.optimize import FIRE
from ase.io import read, write
from ase.constraints import FixAtoms, FixInternals
from ase.calculators.singlepoint import SinglePointCalculator
from sklearn.metrics import mean_absolute_error


def read_extxyz_with_selective(atoms_path, **kwargs):
    atoms = read(atoms_path, **kwargs)
    if 'fix' in atoms.arrays:
        fix_mask = atoms.arrays['fix']
    else:
        fix_mask = [False] * len(atoms)
    fix_indices = [i for i, x in enumerate(fix_mask) if x]
    atoms.set_constraint(FixAtoms(indices=fix_indices))
    return atoms


def write_extxyz_with_selective(filename, atoms):
    fix_mask = [False] * len(atoms)
    if atoms.constraints is not None:
        constraints = atoms.constraints
        for c in constraints:
            if isinstance(c, FixAtoms):
                for index in c.todict()['kwargs']['indices']:
                    fix_mask[index] = True
            else:
                print("WARNING: not FixAtoms", flush=True)
        atoms.arrays['fix'] = np.array(fix_mask)
    write(filename, atoms)


def load_atoms(atoms_path, infer_path):
    atomslist = []
    extxyzs = os.listdir(atoms_path)
    extxyzs = [f for f in extxyzs if f.endswith('.extxyz')]
    for extxyz in extxyzs:
        label = extxyz.replace('.extxyz', '')
        if os.path.exists(f'{infer_path}/{label}.extxyz'):
            print(f"Skipping {label}, already exists in {infer_path}", flush=True)
            continue
        atoms = read_extxyz_with_selective(os.path.join(atoms_path, extxyz))
        atoms.info['label'] = label
        atomslist.append(atoms)
    print(f"Read {len(atomslist)} atoms from {atoms_path}")
    return atomslist


def attach_calc_with_d3(atoms, calc_mlp, functional_name, subtract_dispersion=False):
    from sevenn.calculator import D3Calculator

    if subtract_dispersion:
        from ase.calculators.mixing import MixedCalculator
        calc_d3 = D3Calculator(functional_name=functional_name)
        calc = MixedCalculator(calc_mlp, calc_d3, +1, -1)
        atoms.calc = calc
        print("WARNING: YOU ARE SUBTRACTING D3 CONTRIBUTION.")

    else:
        from ase.calculators.mixing import SumCalculator
        calc_d3 = D3Calculator(functional_name=functional_name)
        calc = SumCalculator([calc_mlp, calc_d3])
        atoms.calc = calc
    

def relax_atoms(
        atomslist,
        infer_path,
        calc_mlp,
        dispersion: bool,
        subtract_dispersion=False,
        functional_name=None,
        fmax=0.02,
        steps=1000,
        cell_relax=False,
        fix_symmetry=False,
        symprec=1e-5,
        ):
    
    if dispersion:
        assert functional_name is not None, \
            "If dispersion is True, functional_name must be provided."
        functional_name = functional_name.lower()
        print(f"Using {functional_name} for D3", flush=True)
        from sevenn.calculator import D3Calculator
        from ase.calculators.mixing import SumCalculator
    else:
        print("No dispersion correction applied.", flush=True)
        functional_name = None
        
    for atoms in atomslist:
        # remove 'charge' or 'spin' info
        if 'charge' in atoms.info:
            del atoms.info['charge']
        if 'spin' in atoms.info:
            del atoms.info['spin']
        
        # get initial positions for rmsd calculation
        try:
            movable_indices = [i for i in range(len(atoms)) if atoms.arrays['fix'][i] == False]
        except KeyError:
            movable_indices = [i for i in range(len(atoms))]
        movable_pos_ini = atoms.get_positions()[movable_indices]

        # set calculator
        if dispersion:
            attach_calc_with_d3(atoms, calc_mlp, functional_name, subtract_dispersion)
        else:
            atoms.calc = calc_mlp
        
        # set symmetry constraints
        if fix_symmetry:
            from ase.constraints import FixSymmetry
            from ase.spacegroup.symmetrize import check_symmetry
            initial_spg_data = check_symmetry(atoms, symprec, verbose=True)
            atoms.set_constraint(FixSymmetry(atoms, symprec=symprec, verbose=True))
            print(f"Fixing symmetry with symprec={symprec}", flush=True)
            # print all constraints
            if atoms.constraints is not None:
                for c in atoms.constraints:
                    print(f"Constraint: {c}", flush=True)
        
        # When having single atoms, just update potential energy and forces and escape the loop
        if len(atoms) == 1:
            print(f"Single atom {atoms.info['label']}, bypassing optimizer...", flush=True)
            atoms.set_pbc(False) # for uma
            atoms.calc = SinglePointCalculator(
                atoms,
                energy=atoms.get_potential_energy(),
                forces=np.array([[0, 0, 0]]),
                stress=np.array([0, 0, 0, 0, 0, 0]),
            )
            # Dummy info
            converged = True
            atoms.info['converged'] = converged
            atoms.info['rmsd'] = 0.0
            atoms.info['max_displacement'] = 0.0
            label = atoms.info['label']
            write_extxyz_with_selective(f'{infer_path}/{label}.extxyz', atoms)
            print(f"Saved {label}.extxyz")
            continue

        # optimize w/ or w/o cell relaxation
        if cell_relax:
            from ase.filters import FrechetCellFilter
            print("Relaxing cell shape using FrechetCellFilter...", flush=True)
            fcf = FrechetCellFilter(atoms)
            opt = FIRE(fcf, logfile=f'{infer_path}/log')
        else:
            print("Relaxing atomic positions only...", flush=True)
            opt = FIRE(atoms, logfile=f'{infer_path}/log')
        converged = opt.run(fmax=fmax, steps=steps)

        if steps > 0:
            print(f"{atoms.info['label']} converged: {converged}", flush=True)
            atoms.info['converged'] = converged
            if fix_symmetry:
                final_spg_data = check_symmetry(atoms, symprec, verbose=True)
                atoms.info['initial_spg'] = initial_spg_data['international']
                atoms.info['final_spg'] = final_spg_data['international']
                if initial_spg_data['international'] != final_spg_data['international']:
                    print(f"Warning: Initial and final space groups differ: "
                          f"{initial_spg_data['international']} vs {final_spg_data['international']}", flush=True)

        # get rmsd and max displacement
        movable_pos_fin = atoms.get_positions()[movable_indices]
        displacements = np.linalg.norm(movable_pos_fin - movable_pos_ini, axis=1)
        rmsd = np.sqrt(np.average(np.square(displacements)))
        max_displacement = np.max(np.abs(displacements))
        atoms.info['rmsd'] = rmsd
        atoms.info['max_displacement'] = max_displacement
    
        label = atoms.info['label']
        write_extxyz_with_selective(f'{infer_path}/{label}.extxyz', atoms)
        print(f"Saved {label}.extxyz")


def get_normal_vector(atoms, indices):
    assert len(indices) == 3, "Indices must be a list of 3 integers"
    v1 = atoms[indices[1]].position - atoms[indices[0]].position
    v2 = atoms[indices[2]].position - atoms[indices[0]].position
    normal = np.cross(v1, v2)
    return normal


def _dihedral_init(atoms: Atoms):

    dihedral_indices = atoms.info["dihedrals"]
    dihedral_indices = [int(i) for i in dihedral_indices.split("-")]

    dihedral_angle = atoms.get_dihedral(*dihedral_indices)
    dihedral_constraint = FixInternals(
        dihedrals_deg=[[dihedral_angle, dihedral_indices]], mic=True
    )
    atoms.set_constraint(dihedral_constraint)
    print(f"Initialized dihedral constraints for indices {dihedral_indices}", flush=True)


def relax_atoms_torsion(
        atomslist,
        infer_path,
        calc_mlp,
        dispersion: bool,
        functional_name=None,
        fmax=0.01,
        steps=1000,
        ):

    if dispersion:
        assert functional_name is not None, \
            "If dispersion is True, functional_name must be provided."
        functional_name = functional_name.lower()
        print(f"Using {functional_name} for D3", flush=True)
        from sevenn.calculator import D3Calculator
        from ase.calculators.mixing import SumCalculator
    else:
        print("No dispersion correction applied.", flush=True)
        functional_name = None

    for atoms in atomslist:
        # append 10 Ang of vacuum if cell is not present
        print("Add vacuum", flush=True)
        atoms.center(vacuum=30.0)

        # set PBC as False
        if atoms.pbc is None:
            print("Setting PBC to False", flush=True)
            atoms.set_pbc(False)

        # remove 'charge' or 'spin' info
        if 'charge' in atoms.info:
            del atoms.info['charge']
        if 'spin' in atoms.info:
            del atoms.info['spin']

        # Initialize dihedral constraints if specified
        assert "dihedrals" in atoms.info, \
            "Atoms info must contain 'dihedrals' key with indices."
        _dihedral_init(atoms)

        # set calculator
        if dispersion:
            calc_d3 = D3Calculator(functional_name=functional_name)
            calc = SumCalculator([calc_mlp, calc_d3])
            atoms.calc = calc
        else:
            atoms.calc = calc_mlp

        pos_ini = atoms.get_positions()
        print("Relaxing atomic positions with dihedral constraint...", flush=True)
        opt = FIRE(atoms, logfile=f'{infer_path}/log')
        converged = opt.run(fmax=fmax, steps=steps)
        if steps > 0:
            print(f"{atoms.info['label']} converged: {converged}", flush=True)
            atoms.info['converged'] = converged

        # get rmsd and max displacement
        pos_fin = atoms.get_positions()
        displacements = np.linalg.norm(pos_fin - pos_ini, axis=1)
        rmsd = np.sqrt(np.average(np.square(displacements)))
        max_displacement = np.max(np.abs(displacements))
        atoms.info['rmsd'] = rmsd
        atoms.info['max_displacement'] = max_displacement
        
        label = atoms.info['label']
        write_extxyz_with_selective(f'{infer_path}/{label}.extxyz', atoms)
        print(f"Saved {label}.extxyz")


def infer(
    atoms_path: str,
    infer_path: str,
    calc_mlp,
    **kwargs
):
    atomslist = load_atoms(atoms_path, infer_path)
    relax_atoms(
        atomslist=atomslist,
        infer_path=infer_path,
        calc_mlp=calc_mlp,
        **kwargs
    )
    
def infer_torsion(
    atoms_path: str,
    infer_path: str,
    calc_mlp,
    **kwargs
):
    atomslist = load_atoms(atoms_path, infer_path)
    relax_atoms_torsion(
        atomslist=atomslist,
        infer_path=infer_path,
        calc_mlp=calc_mlp,
        **kwargs
    )

def oneshot_atomslist(
        atoms_path,
        infer_path,
        calc_mlp,
        dispersion: bool,
        functional_name=None,
        ):
    
    from tqdm import tqdm
    
    atomslist = read(atoms_path, ":")
    if os.path.exists(f'{infer_path}/total.extxyz'):
        print(f"Skipping {infer_path}/total.extxyz, already exists.", flush=True)
        return

    if dispersion:
        assert functional_name is not None, \
            "If dispersion is True, functional_name must be provided."
        functional_name = functional_name.lower()
        print(f"Using {functional_name} for D3", flush=True)
        from sevenn.calculator import D3Calculator
        from ase.calculators.mixing import SumCalculator
    else:
        print("No dispersion correction applied.", flush=True)
        functional_name = None

    atomslist_calc = []   
    for atoms in tqdm(atomslist):
        # remove 'charge' or 'spin' info
        if 'charge' in atoms.info:
            del atoms.info['charge']
        if 'spin' in atoms.info:
            del atoms.info['spin']
        
        # set calculator        
        if dispersion:
            calc_d3 = D3Calculator(functional_name=functional_name)
            calc = SumCalculator([calc_mlp, calc_d3])
            atoms.calc = calc
        else:
            atoms.calc = calc_mlp

        atoms.calc = SinglePointCalculator(
            atoms,
            energy=atoms.get_potential_energy(),
            forces=atoms.get_forces(),
            stress=atoms.get_stress()
        )   
        atomslist_calc.append(atoms)

    write(f'{infer_path}/total.extxyz',
          atomslist_calc, format='extxyz')


def collect_results(infer_path, dft=False, results_filename='indiv_results.csv'):
    data = []
    extxyz_files = [f for f in os.listdir(infer_path) if f.endswith('.extxyz')]

    for filename in extxyz_files:
        name = os.path.splitext(filename)[0]
        file_path = os.path.join(infer_path, filename)
        try:
            atoms = read(file_path)
            predicted_energy = atoms.get_potential_energy()

        except Exception as e:
            print(f"Error reading {file_path}: {e}", flush=True)
            predicted_energy = None
        
        if dft:
            data.append({
                'name': name,
                'dft_energy': predicted_energy,
                })
        else:
            data.append({
                'name': name,
                'pred_energy': predicted_energy,
                'converged': atoms.info.get('converged', None),
                'rmsd': atoms.info.get('rmsd', None),
                'max_displacement': atoms.info.get('max_displacement', None),
                })

    df = pd.DataFrame(data)
    output_csv_path = os.path.join(infer_path, results_filename)
    df.to_csv(output_csv_path, index=False)
    print(f"Predicted energies saved to {output_csv_path}", flush=True)


def collect_oneshot_results(infer_path, filename='total.extxyz'):
    import numpy as np
    import pandas as pd
    from ase.io import read
    from ase.units import bar

    atomslist = read(f'{infer_path}/{filename}', ":")
    data = []

    # to per_graph.csv
    num_atoms = []
    energies_pred = []
    stresses_pred = []

    # to per_atom.csv
    stct_id = []
    atom_id = []
    atomic_numbers = []
    forces_pred = []

    for stct, atoms in enumerate(atomslist):
        # to per_graph.csv
        num_atoms.append(len(atoms))
        energies_pred.append(atoms.get_potential_energy())
        stresses_pred.append(atoms.get_stress(voigt=True) / (-1000 * bar))
        
        # to per_atom.csv
        stct_id.extend([stct] * len(atoms))
        atom_id.extend(range(len(atoms)))

        numbers = atoms.get_atomic_numbers()
        atomic_numbers.extend(numbers)

        forces = atoms.get_forces()
        forces_pred.append(forces)

    num_atoms = np.array(num_atoms)
    energies_pred = np.array(energies_pred)
    stresses_pred = np.array(stresses_pred)

    stct_id = np.array(stct_id)
    atom_id = np.array(atom_id)
    forces_pred = np.vstack(forces_pred)

    # save labels, energies_ref, energies_pred
    df = pd.DataFrame({
        'inferred_total_energy': energies_pred,
        'inferred_stress_xx': stresses_pred[:,0],
        'inferred_stress_yy': stresses_pred[:,1],
        'inferred_stress_zz': stresses_pred[:,2],
        'inferred_stress_xy': stresses_pred[:,5],
        'inferred_stress_yz': stresses_pred[:,3],
        'inferred_stress_zx': stresses_pred[:,4],
    })
    df.to_csv(f'{infer_path}/per_graph.csv', index=False)

    # save forces_ref, forces_pred
    df = pd.DataFrame({
        'stct_id': stct_id,
        'atom_id': atom_id,
        'atomic_numbers': atomic_numbers,
        'inferred_force_x': forces_pred[:, 0],
        'inferred_force_y': forces_pred[:, 1],
        'inferred_force_z': forces_pred[:, 2],
    })
    df.to_csv(f'{infer_path}/per_atom.csv', index=False)

    
def calculate_reaction_energies(
        rxn_sheet_csv,
        infer_path,
        dft=False,
        results_filename='indiv_results.csv',
        rxn_filename='rxn_energies.csv',
        ):
    indiv_results_csv = f'{infer_path}/{results_filename}'
    output_csv = f'{infer_path}/{rxn_filename}'
    
    assert os.path.exists(indiv_results_csv), \
        f"File not found: {indiv_results_csv}"

    energy_df = pd.read_csv(indiv_results_csv)
    energy_df = energy_df.set_index('name')
    rxn_df = pd.read_csv(rxn_sheet_csv, header=None)

    rxn_names = []
    pred_rxn_energies = []
    for i, row in rxn_df.iterrows():
        # e.g. [(1, 'IF5_hf_1'), (-1, 'IF5_1'), (-1, 'hf')]
        # if it has nan values after specific index, terminate there
        rxn_name = row[0]
        terms = list(zip(row[1::2], row[2::2]))
        for j, (coeff, name) in enumerate(terms):
            if pd.isna(coeff) or pd.isna(name):
                terms = terms[:j]
                break
        
        coeff_energy_pred = 0
        has_none = False

        for coeff, name in terms:
            try:
                if dft:
                    e_pred = energy_df.at[name, 'dft_energy']
                else:
                    e_pred = energy_df.at[name, 'pred_energy']
                print(f"Processing {name}: {e_pred}", flush=True)
                if pd.isna(e_pred):
                    print(f"Warning: {name} has NaN value in energy_df.")
                    has_none = True
                    break
                coeff_energy_pred += coeff * e_pred
            except (KeyError, ValueError):
                print(f"Warning: {name} not found in energy_df or has NaN value.", flush=True)
                has_none = True
                break

        rxn_names.append(rxn_name)
        pred_rxn_energies.append(None if has_none else coeff_energy_pred)

    if dft:
        output_df = pd.DataFrame({
            'rxn_name': rxn_names,
            'dft_rxn_energy': pred_rxn_energies
        })

    else:
        output_df = pd.DataFrame({
            'rxn_name': rxn_names,
            'pred_rxn_energy': pred_rxn_energies
        })
    output_df.to_csv(output_csv, index=False)
