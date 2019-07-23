#!/bin/bash
#SBATCH --job-name=bowtie2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=END
#SBATCH --mem=50G
#SBATCH --mail-user=your.email@uconn.edu
#SBATCH -o bowtie2_%j.out
#SBATCH -e bowtie2_%j.err
cat Sample_1.fastq Sample_2.fastq > genome.fastq
mkdir index
cd index
module load bowtie2/2.3.3.1

#Masurca
bowtie2-build /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/samples/CA/final.genome.scf.fasta masurca.index

bowtie2 -x masurca.index -U /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/genome.fastq -S masurca.bowtie2.sam

#SPAdes
bowtie2-build /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/Assembly/SPAdes_out/scaffolds.fasta SPAdes.index

bowtie2 -x SPAdes.index -U /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/genome.fastq -S SPAdes.bowtie2.sam

#SOAP31
bowtie2-build /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/Assembly/SOAP/graph_Sample_31.scafSeq SOAP31.index

bowtie2 -x SOAP31.index -U /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/genome.fastq -S SOAP31.bowtie2.sam

#SOAP35
bowtie2-build /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/Assembly/SOAP/graph_Sample_35.scafSeq SOAP35.index

bowtie2 -x SOAP35.index -U /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/genome.fastq -S SOAP35.bowtie2.sam

#SOAP41
bowtie2-build /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/Assembly/SOAP/graph_Sample_41.scafSeq SOAP41.index

bowtie2 -x SOAP41.index -U /home/CAM/jguillemin/Assembly_tutorial/Assembly_tutorial/genome.fastq -S SOAP41.bowtie2.sam
