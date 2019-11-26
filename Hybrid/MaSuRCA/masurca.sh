#!/bin/bash
#SBATCH --job-name=masurca
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mem=450G
#SBATCH --mail-user=sunny.sarker@uconn.edu
#SBATCH -o masurca_%j.out
#SBATCH -e masurca_%j.err

module load singularity/3.1.1
module load MaSuRCA/3.3.3


#masurca config.txt

./assemble.sh
