#!/bin/bash
#SBATCH --job-name=purge_haplotigs_step2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=jacqueline.guillemin@uconn.edu
#SBATCH -o purge_halpotigs_step2_%j.out
#SBATCH -e purge_haplotigs_step2_%j.err

module load purge_haplotigs/1.0

purge_haplotigs contigcov -i physcomitrellopsis_africana_genome.reads.sorted.bam.gencov -l 3 -m 57 -h 195

