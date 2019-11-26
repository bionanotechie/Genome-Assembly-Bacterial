#!/bin/bash
#SBATCH --job-name=falc
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=80G
#SBATCH --mail-user=sunny.sarker@uconn.edu
#SBATCH -o falcon_%j.out
#SBATCH -e falcon_%j.err

#module load Miniconda/Miniconda3

. ~/miniconda3/etc/profile.d/conda.sh
conda activate denovo_py3



fc_run fc_run.cfg

#fc_unzip.py fc_unzip.cfg
