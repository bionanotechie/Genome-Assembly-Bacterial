
#!/bin/bash
#SBATCH --job-name=Sample_QC_job
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err


#gunzip *.gz

module load sickle/1.33

sickle pe -f /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Sample_R1.fastq -r /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Sample_R2.fastq -t sanger -o Sample_1.fastq -p Sample_2.fastq -s Sample_s.fastq -q 30 -l 45
