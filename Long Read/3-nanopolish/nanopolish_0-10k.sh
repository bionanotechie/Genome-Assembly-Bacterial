#!/bin/bash
#SBATCH --job-name=nanopolish_tutorial
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=250G
#SBATCH --mail-user=sunny.sarker@uconn.edu
#SBATCH -o nanopolish_%j.out
#SBATCH -e nanopolish_%j.err

module load nanopolish/0.11.1
module load python/3.6.3
module load biopython/1.70
module load parallel/20180122

python nanopolish_makerange.py flye_moss_rm3kb_0-10k.fasta | parallel --results nanopolish.results -P 8 --tmpdir nanopolish_contigs_1_tmp/ \ nanopolish variants --consensus -o vcf_out/nanopolish_flye_1.{1}.vcf -w {1} -r /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/5074_test_LSK109_30JAN19-reads/5074_test_LSK109_30JAN19-reads-pass.fastq  -b /labs/Wegrzyn/Moss/Physcomitrellopsis_africana/Physcomitrellopsis_africana_Genome/RawData_Nanopore_5074/5074_test_LSK109_30JAN19/nanopolish_scripts/physcomitrellopsis_africana_genome.reads.sorted.bam -g /labs/Wegrzyn/genome_assembly_tutorial/nanopolish_flye/flye_moss_rm3kb_0-10k.fasta -t 8  --min-candidate-frequency 0.1
