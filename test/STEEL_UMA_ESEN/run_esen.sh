#!/bin/bash
#SBATCH --job-name=esen
#SBATCH --output=esen.x
#SBATCH --error=esen.x
#SBATCH --nodes=1   
#SBATCH --ntasks=1
#SBATCH --time=3:00:00 
#SBATCH --partition=gpu2
#SBATCH --gres=gpu:1

export FLASH=1

echo "SLURM_NTASKS: $SLURM_NTASKS"

# source ~/.bash_profile

if [ -z "$SLURM_NTASKS" ] || [ "$SLURM_NTASKS" -le 0 ]; then
	echo "Error: SLURM_NTASKS is not set or is less than or equal to 0"
	exit 1
fi

python script_esen.py
