#!/bin/bash
#SBATCH --job-name=nanopolish_tutorial
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=60G
#SBATCH --mail-user=jacqueline.guillemin@uconn.edu
#SBATCH -o nanopolish_%j.out
#SBATCH -e nanopolish_%j.err

module load nanopolish/0.11.1
module load python/3.6.3
module load biopython/1.70
module load parallel/20180122

python /labs/Wegrzyn/genome_assembly_tutorial/nanopolish/nanopolish_makerange.py /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/flye_assembly/assembly.fasta | parallel --results /labs/Wegrzyn/genome_assembly_tutorial/nanopolish/results -p 9 \
nanopolish variants --consensus -o /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/nanopolish_contigs/physcomitrellopsis_africana_assembly_polished. {1}.vcf -w {1} -r /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/5074_test_LSK109_30JAN19-read/5074_test_LSK109_30JAN19-reads-pass.fastq.index -b /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/nanopolish_scripts/physcomitrellopsis_africana_genome.reads.sorted.bam -g /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/Flye_assembly/assembly.fasta -t --min-canidate-frequency 0.1
