from parse_esen import load_esen_calc
from relax_metrics import relax_pureFe, relax_carbons, relax_tms

import os

basedir = os.path.abspath(os.getcwd())
print('calculating with eSEN-omat')
os.chdir(f'{basedir}/./esen_omat')
calc = load_esen_calc('esen_omat')

relax_pureFe(calc)
relax_carbons(calc)
relax_tms(calc)

os.chdir(basedir)

print('calculating with eSEN')
os.chdir(f'{basedir}/./esen')
calc = load_esen_calc('esen')

relax_pureFe(calc)
relax_carbons(calc)
relax_tms(calc)


