#!/bin/bash
#SBATCH --job-name=busco_polished_tutorial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=jacqueline.guillemin@uconn.edu
#SBATCH -o busco_%j.out
#SBATCH -e busco_%j.err

module load busco/3.0.2b
module unload augustus
export PATH=/home/CAM/jguillemin/augustus/bin:/home/CAM/jguillemin/augustus/scripts:$PATH
export AUGUSTUS_CONFIG_PATH=$HOME/augustus/config
###if your run crashes uncomment the following:
##module unload blast/2.7.1
##module load blast/2.2.31

run_BUSCO.py -i /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/polished_assembly/physcomitrellopsis_africana_polished_3kb_assembly.fasta -l /labs/Wegrzyn/Moss/Physcomitrium/viridi/ -o physcomitrellopsis_africana_polished_busco_tutorial -m geno -c 1
