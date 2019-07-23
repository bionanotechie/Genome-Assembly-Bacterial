#!/bin/bash
#SBATCH --job-name=flye_assembly_tutorial
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 36
#SBATCH --partition=himem
#SBATCH --qos=himem
#SBATCH --mail-type=END
#SBATCH --mem=500G
#SBATCH --mail-user=jacqueline.guillemin@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
module load flye/2.4.2
flye --nano-raw /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/5074_test_LSK109_30JAN19-reads/5074_test_LSK109_30JAN19-reads-pass.fastq
 --genome-size 1g --out-dir /labs/Wegrzyn/Long_read_assembly_tutorial --threads 36
