#!/bin/bash
#SBATCH --job-name=Sample_quast
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

hostname
module load quast/5.0.2

# SOAPdenovo statistics
quast.py /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SOAP/graph_Sample_*.scafSeq -o SOAP

# SPAdes statistics
quast.py /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SPAdes/scaffolds.fasta -o SPAdes

#MaSuRCA statistics
quast.py /Path_to_tutorial/Assembly_Tutorial/samples/CA/final.genome.scf.fasta -o MaSuRCA
