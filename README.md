### Description
python package for reproducing the following paper: 

S. E. Restrepo, N. K. Mohandas, M. Sluiter, and A. T. Paxton, Applicability of universal machine learning interatomic potentials to the simulation of steels, Model. Simul. Mater. Sci. Eng. (2025).

### Installation
pip install -e .

### Usage
for command line interface .. use 
febench --config ./config.yaml --calc omni --modal omat24 --carbon_config ./carbon_config.yaml
indicated arguments above are the default argument of the parser

to run pure iron, carbon in iron, and transition metal solute/vacancy interactions
you can also tweak the config.yaml's 'run' key's value to false to skip a certain type of calculation

other command lines available
febench (same with febench-run)

febench-pure to only run pure iron calc.

febench-carbon to only run carbon in iron calc.

febench-tm to only run transition metal solute/vacancy interactions

