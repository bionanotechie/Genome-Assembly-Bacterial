#!/bin/bash
#SBATCH --job-name=quast_long
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sunny.sarker@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
module load quast/5.0.2

# flye statistics
quast.py /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/flye_assembly/flye_assembly_initial/assembly.fasta -o Flye

# shasta statistics
quast.py /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/test_shasta_assembly/ShastaRun_pafricana_rmv_contam_minreadlen_500/Assembly.fasta -o Shasta

#falcon statistics

#masurca stats
