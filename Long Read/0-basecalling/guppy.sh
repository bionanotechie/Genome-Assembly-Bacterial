#!/bin/bash
#SBATCH --job-name=guppy_test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user= username@uconn.edu
#SBATCH -o guppy_%j.out
#SBATCH -e guppy_%j.err

module load guppy/2.3.1

guppy basecaller -i /path_to_reads --cpu_threads_per_caller 36 --flowcell FLO-MIN06 --kit SQK-LSK109 --qscore_filtering 
