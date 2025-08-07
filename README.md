### Description
python package for reproducing the following paper: 

S. E. Restrepo, N. K. Mohandas, M. Sluiter, and A. T. Paxton, Applicability of universal machine learning interatomic potentials to the simulation of steels, Model. Simul. Mater. Sci. Eng. (2025).

### Installation
pip install -e .

### Usage
run the following  command to run the module as a whole
febench --config ./config.yaml --modal omat24 --calc ompa --potential_path . --potential_ext pth

this reads the binary file {args.potential_path}/{args.calc}.{args.potential_ext}
as an ASE calculator object

go to febench/util/parse_args.py for default arguments pass
or febench/main.py for usage on scripts

### config.yaml
holds calc. parameters concerning each task
primary forcus is to reimplement the aforementioned paper above

### other command lines available
febench (same with febench-run)

febench-pure to only run pure iron calc.

febench-carbon to only run carbon in iron calc.

febench-tm to only run transition metal solute/vacancy interactions

or tweak config.yaml's 'run' key-value

