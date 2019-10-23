#!/bin/bash
#SBATCH --job-name=purge_haplotigs_step3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=jacqueline.guillemin@uconn.edu
#SBATCH -o purge_haplotigs_step3_%j.out
#SBATCH -e purge_haplotigs_step3_%j.err

module load purge_haplotigs/1.0

purge_haplotigs purge -g /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/flye_assembly/assembly.fasta -c coverage_stats.csv -a 60 -d -b /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/physcomitrellopsis_africana_genome.reads.sorted.bam
