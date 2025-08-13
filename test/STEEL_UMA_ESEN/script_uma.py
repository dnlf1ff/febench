from parse_uma import load_uma_calc
from relax_metrics import relax_pureFe, relax_carbons, relax_tms

import os

basedir = os.path.abspath(os.getcwd())

print('calculating with uma_omat')
os.chdir(f'{basedir}/./uma_omat')
calc = load_uma_calc(modal='omat')
relax_pureFe(calc)
relax_carbons(calc)
relax_tms(calc)

os.chdir(basedir)

print('calculating with uma_omc')
os.chdir(f'{basedir}/./uma_omc')
calc = load_uma_calc('omc', dispersion=True)
relax_pureFe(calc)
relax_carbons(calc)
relax_tms(calc)


