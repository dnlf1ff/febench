### Description
python package for reproducing the following paper: 

S. E. Restrepo, N. K. Mohandas, M. Sluiter, and A. T. Paxton, Applicability of universal machine learning interatomic potentials to the simulation of steels, Model. Simul. Mater. Sci. Eng. (2025).

### Usage
run the following  command to run the module as a whole <br>
febench --config ./config.yaml --cwd ompa_omat --modal omat24 --model checkpoint.pth --calc_type 7net-mf --potential_path $PATH_TO_POTENTIAL --dispersion false --functional PBE <br><br>


this reads the binary file {args.potential_path}/{args.model} as an ASE calculator object <br>
set --dispersion to true in order to exclude D3 calc. <br>

for general usage of argument passing, go to febench/util/parse_args.py & febench/util/parse_config.py<br>
for detailed calculator configurations see febench/util/parse_calc.py for details<br>
<details><summary style="background-color:white;color:green;font-weight:normal;width:220px;">notes</summary>
     DPA, ORB, UMA requires numpy >= 2.3.0 which conflicts with febench/pureFe/script.process_stiffness <br>
     .. which is a rather trivial prob. as our target property is the binding energy b/w solutes and vacancies <br>
     Skip the stiffness tensor part .. what matters in febench/pureFe is the relaxation of bulk + vacancy <br>
</details>
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

2. install required libraries for mlips

   <details><summary style="background-color:white;color:green;font-weight:normal;width:220px;">SevenNet</summary>
    pip install torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121'<br>
    pip install torch==2.5.1 torchvision torchaudio --index-url https://download.pytorch.org/whl/cu121<br>
    pip install torch-scatter -f https://data.pyg.org/whl/torch-2.5.1+cu121.html <br>
    pip install sevenn ase matscipy <br>
    </details>

   <details><summary style="background-color:white;color:green;font-weight:normal;width:220px;">DPA31</summary>
    pip install deepmd-kit[torch]<br>
    pip install sevenn <br>
    pip install numpy==2.3.1 <br>
    </details>

   <details><summary style="background-color:white;color:green;font-weight:normal;width:220px;">eSEN</summary>
    pip install fairchem-core==1.10.0 <br>
    pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.4.1+cu121.html <br>
    pip install sevenn <br>
    </details>

   <details><summary style="background-color:white;color:green;font-weight:normal;width:220px;">ORB</summary>
    pip install orb-models <br>
    pip install sevenn <br>
    pip install numpy==2.3.1<br>
    </details>

   <details><summary style="background-color:white;color:green;font-weight:normal;width:220px;">UMA</summary>
    pip install fairchem-core <br>
    pip install sevenn <br>
    pip install numpy==2.2.6 <br>
    </details>


3. install package  <br>
    git clone git@github.com:dnlf1ff/febench.git <br>
    cd ./febench <br>
    pip install .  <br>
