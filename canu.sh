#!/bin/bash

#SBATCH --job-name=canu_limulus_assembly

#SBATCH -N 1

#SBATCH -n 1

#SBATCH -c 24

#SBATCH --partition=general

#SBATCH --qos=general

#SBATCH --mail-type=END

#SBATCH --mem=200G

#SBATCH --mail-user=sunny.sarker@uconn.edu

#SBATCH -o %x_%j.out

#SBATCH -e %x_%j.err


#module load gnuplot/5.2.2

module load canu/1.8

canu -p limulus_polyphemus -d limulus_assembly_canu genomeSize=1g -nanopore-raw /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/5074_test_LSK109_30JAN19-reads/5074_test_LSK109_30JAN19-reads-pass.fastq gridOptions="--partition=general --qos=general --mem-per-cpu=8032m --cpus-per-task=24" canuIteration=1

