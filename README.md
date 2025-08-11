### Description
python package for reproducing the following paper: 

S. E. Restrepo, N. K. Mohandas, M. Sluiter, and A. T. Paxton, Applicability of universal machine learning interatomic potentials to the simulation of steels, Model. Simul. Mater. Sci. Eng. (2025).

### Installation
pip install .

### Usage
run the following  command to run the module as a whole <br>
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

### environment
1. create virtual environment (venv works as well)
micromamba create -n febench python=3.11
micromamba activate febench

2. install required libraries for mlips

    <details><<summary style="background-color:white;color:green;font-weight:normal;width:220px;">SevenNet</summary>
    > pip install torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121<br>
> pip install torch==2.5.1 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121<br>
> pip install torch-scatter -f https://data.pyg.org/whl/torch-2.5.1+cu121.html <br>
> pip install sevenn ase matscipy <br>
    </details>

3. install package 
> git clone git@github.com:dnlf1ff/febench.git
> cd ./febench
> pip install . 
