#!/bin/bash
#SBATCH --job-name=Sample_assembly
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=30G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=firs.name@uconn.edu
#SBATCH -o %x_%A.out
#SBATCH -e %x_%A.err

# SOAPdenovo assembly


module load SOAP-denovo/2.04

cd /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SOAP


SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 31 -R -o graph_Sample_31 1>ass31.log 2>ass31.err
SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 35 -R -o graph_Sample_35 1>ass35.log 2>ass35.err
SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 41 -R -o graph_Sample_41 1>ass41.log 2>ass41.err
module unload SOAP-denovo/2.04



# SPAdes assembly


module load SPAdes/3.13.0


spades.py --careful -o SPAdes -1 /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_1.fastq -2 /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_2.fastq -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_s.fastq
module unload SPAdes/3.13.0



# MaSuRCA assembly (without sickle)


cd /Path_to_tutorial/Assembly_Tutorial/Assembly/MaSuRCA



#run MaSuRCA


module load MaSuRCA/3.2.4


masurca config.txt


bash assemble.sh
