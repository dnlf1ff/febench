from ase.filters import UnitCellFilter, FrechetCellFilter
from ase.optimize import FIRE
import numpy as np
from ase import Atoms
from ase.io import Trajectory

OPT_DICT = {'fire': FIRE}
FILTER_DICT = {'frechet': FrechetCellFilter, 'unitcell': UnitCellFilter}


"""
modified based on Jaesun Kim's code
"""

class AseAtomRelax:
    def __init__(
        self,
        calc,
        optimizer,
        cell_filter,
        mask,
        fmax=0.0001,
        steps=10000,
        logfile='ase_relaxer.log',
        trajfile='ase_traj.traj',

    ):
        self.calc = calc
        self.optimizer = optimizer
        self.cell_filter = cell_filter
        self.mask = mask
        self.fmax = fmax
        self.steps = steps
        self.logfile = logfile
        self.trajfile = trajfile

    def update_atoms(self, atoms):
        atoms = atoms.copy()
        atoms.calc = self.calc

        try:
            atoms.info['e_fr_energy'] = atoms.get_potential_energy(force_consistent=True)
        except:
            atoms.info['e_fr_energy'] = atoms.get_potential_energy()
        atoms.info['e_0_energy'] = atoms.get_potential_energy()
        atoms.info['force'] = atoms.get_forces()
        atoms.info['volume'] = atoms.get_volume()
        atoms.info['a'] = atoms.cell.lengths()[0]
        atoms.info['b'] = atoms.cell.lengths()[1]
        atoms.info['c'] = atoms.cell.lengths()[2]
        atoms.info['alpha'] = atoms.cell.angles()[0]
        atoms.info['beta'] = atoms.cell.angles()[1]
        atoms.info['gamma'] = atoms.cell.angles()[2]
        return atoms

    def relax_atoms(self, atoms):
        atoms = atoms.copy()
        atoms.calc = self.calc

        cell_filter = self.cell_filter(atoms, mask=self.mask)
        opt = self.optimizer(cell_filter, logfile=self.logfile)
        traj = Trajectory(filename=self.trajfile, mode='w', atoms=atoms)
        opt.attach(traj.write, interval=20)
        opt.run(fmax=self.fmax, steps=self.steps)
        conv = check_atoms_conv(atoms.get_forces())
        traj.close()
        return atoms, conv

def aar_from_config(config, calc, logfile, trajfile, opt_type='carbon'):
    arr_args = config['opt'][opt_type].copy()
    opt = OPT_DICT['fire']
    cell_filter = FILTER_DICT['unitcell']

    arr_args['calc'] = calc
    arr_args['optimizer'] = opt
    arr_args['cell_filter'] = cell_filter
    arr_args['logfile'] = logfile 
    arr_args['trajfile'] = trajfile 

    return AseAtomRelax(**arr_args)

def check_atoms_conv(forces: np.ndarray) -> bool:
    conv = True
    for i in range(forces.shape[-1]):
        if np.any(forces[:,i]) < 0:
            conv = False
    return conv


