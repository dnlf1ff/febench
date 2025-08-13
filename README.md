#### Branch including code for calculating the components ($$C_{ij}$$) of the stiffness tensor($$C$$) of pure iron[1]
<br>
as matscipy requires numpy > 3.0.0 || numpy < 2.0.0 which conflicts with u-MLIPs requiring numpy >= 2.3.0, this metric is singled out
<br>
[1] M. de Jong, W. Chen, T. Angsten, A. Jain, R. Notestine, A. Gamst, M. Sluiter, C.K. Ande, S. van der Zwaag, J.J. Plata, C. Toher, S. Curtarolo, G. Ceder, K.A. Persson, M. Asta, Charting the complete elastic properties of inorganic crystalline compounds, Sci. Data 2 (2015) 150009. https://doi.org/10.1038/sdata.2015.9.


### Description
python package for reproducing the following paper: 

S. E. Restrepo, N. K. Mohandas, M. Sluiter, and A. T. Paxton, Applicability of universal machine learning interatomic potentials to the simulation of steels, Model. Simul. Mater. Sci. Eng. (2025).

### Usage
run the following  command to run the module as a whole <br>
febench --config ./config.yaml --cwd ompa_omat --modal omat24 --model checkpoint.pth --calc_type 7net-mf --potential_path $PATH_TO_POTENTIAL --dispersion false --functional PBE <br><br>


for general usage of argument passing, go to febench/util/parse_args.py & febench/util/parse_config.py<br>
for detailed calculator configurations see febench/util/parse_calc.py for details<br>

<br>
#### config.yaml
holds calc. parameters concerning each task
primary forcus is to reimplement the aforementioned paper above

#### other command lines
febench (same with febench-run)

febench-pure to only run pure iron calc.

febench-carbon to only run carbon in iron calc.

febench-tm to only run transition metal solute/vacancy interactions

or tweak config.yaml's 'run' key-value

### Installation
1. create virtual environment (venv works as well) <br>
micromamba create -n febench python=3.11 <br>
micromamba activate febench <br>

2. install required libraries
``` 
pip install matscipy ase torch tqdm pandas numpy
# pip install sevenn
```
3. install package  <br>
    git clone --branch matscipy git@github.com:dnlf1ff/febench.git <br>
    cd ./febench <br>
    pip install .  <br>
