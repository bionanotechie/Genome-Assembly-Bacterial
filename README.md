# Genome Assembly  
  
This repository is a usable, publicly available tutorial. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   
  
## Table of Contents  
1. [Short Read Genome Assembly](#short)  
   - [Overview](#over)
   - [Copy the Assembly Directory to your account node](#copy)
   - [Quality Control with Sickle](#sickle)
   - [Assembly with SOAPdenovo, SPAdes, and MaSuRCA](#short-assemble)
     - [Assembly with SOAPdenovo](#soap)
     - [Assembly with SPAdes](#spades)
     - [Assembly with MaSuRCA](#ma)
   - [Assembly Statistics with QUAST](#quast)
2. [Long Read Genome Assembly](#long)
   - [Base Calling with Guppy](#gup)
   - [Assembly with Flye and Falcon](#ff)
   - [Checking completeness with BUSCOMP](#bus)
   - [Polishing with Nanopolish](#nano)
   - [Organizing with Purge Haplotigs](#ph)
3. [Hybrid Assembly](#ha)


<a name="short"></a>
# Short Read Genome Assembly

<a name="over"></a>
## Overview

This tutorial will teach you how to use open source quality control, genome assembly, and assembly assessment tools to complete a high quality de novo assembly, this means that you do not have any data to compare back to. Moving through the tutorial you will take pair end short read data from a bacterial species and perform assemblies via MaSuRCA, SOAPDenovo, and SPAdes. With these assemblies completed we will use a program called QUAST to assess the quality of the data produced. Once finished with the short read data we will move on to performing a long read assembly.
In this tutorial, you will work with common genome assemblers, such as 
<a href="https://github.com/aquaskyline/SOAPdenovo2">SOAPdenovo</a>, 
<a href="https://github.com/ablab/spades">SPAdes</a>, 
<a href="https://github.com/alekseyzimin/masurca">MaSuRCA</a>, 
and quality assessment tool <a href="https://github.com/ablab/quast">QUAST</a>

In order to do this tutorial you will need to need to run submission scripts, the structure of this is explained 
<a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/">here</a>.

If you would like to include your email to be notified when the jobs are done enter the script doc and add your email to the user mail line.

<a name="copy"></a>
## Step 1: Copy the Assembly Directory to your account node 

Run the following command to download the directory:

```
cp -avr /UCHC/PublicShare/Tutorials/Assembly_Tutorial your_directory
```
Wait for download to finish (3.6 gb file will take a bit).

Make sure to not be in the head node in order for the download to be quick and secure.

Enter the directory you created.


<a name="sickle"></a>
## Step 2: Quality Control with Sickle
Sickle takes raw reads and outputs data with the 3’ and 5’ ends trimmed to assure that the quality of the read is high enough for assembly, it will also trim low quality reads. 

The flags meanings are as follows:

- The pe flag stands for pair-end, sickle also has the ability to perform quality control on single end reads. 
- the -f flag designates the input file containing the forward reads.
- -r is the input file containing the reverse reads.
- -o is the output file containing the trimmed forward reads.
- -p is the output file containing the trimmed reverse reads.
- -s is the output file containing trimmed singles. 
- -q flag represents the minimum quality of the output.
- -l is the minimum read length.
- -t is the type of read you are inputting.

### Running Sickle
```
Module load sickle/1.33
sickle pe -f /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Sample_R1.fastq -r /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Sample_R2.fastq -t sanger -o Sample_1.fastq -p Sample_2.fastq -s Sample_s.fastq -q 30 -l 45
```

The commands are located in Sample_QC.sh in Quality Control.

Run the shell script file with sbatch.

The files created appear as:
```
File_with_tutorial/
|-- Sample_1.fastq
|-- Sample_2.fastq
|-- Sample_s.fastq
```

<a name="short-assemble"></a>
## Step 3: Assembly with SOAPdenovo, SPAdes, and MaSuRCA

Run Sample_assembly.sh in the Assembly folder to perform SOAPdenovo, SPAdes, and MaSuRCA at once.
```
sbatch Sample_assembly.sh
```
You should expect an output of an .out and .err file. Check these to assure your assembly ran properly.

Here is an explanation of each step within Sample_assembly.sh:
<a name="soap"></a>
### **Assembly with SOAPdenovo:**

This is a de novo assembler, this assembler, like MaSuRCA which we will be encountering later, requires a config file to run through the data. The configuration file can be found here https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/SOAPdenovo%20config%20file

**Running SOAPdenovo:**

Run SOAPdenovo with the following commands:
```
module load SOAP-denovo/2.04

cd /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SOAP


SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 31 -R -o graph_Sample_31 1>ass31.log 2>ass31.err
SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 35 -R -o graph_Sample_35 1>ass35.log 2>ass35.err
SOAPdenovo-63mer all -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/Sample.config -K 41 -R -o graph_Sample_41 1>ass41.log 2>ass41.err

module unload SOAP-denovo/2.04
```

The meaning of the flags within the command have these meanings: 
- -s is the path to the config file
- -K is the size of the k-mer, a k-mer is a set of nucleotides, k is the number of nucleotides in that set. (biostar)
- -o is the output file
- 1 is for the assembly log and 2 is for the assembly errors. 

For this assembly we use the reads that have been run through Sickle for quality control.
<a name="spades"></a>
### **Assembly with SPAdes:**

Instead of manually selecting k-mers, SPAdes automatically selects k-mers based off the maximum read length data of your input. This is a called a de Bruijn graph based assembler, meaning that it assigns (k-1)-mers to nodes and every possible matching prefix and suffix of these nodes are connected with a line(Compeau). 

**Running SPAdes:**

Use the spades.py with the following parameters:
```
module load SPAdes/3.13.0

spades.py --careful -o SPAdes -1 /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_1.fastq -2 /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_2.fastq -s /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_s.fastq

module unload SPAdes/3.13.0
```
The meanings of the flags are:

- --careful to reduce mismatches in the contigs 
- -o for the output folder
- -1 for the location of forward reads file
- -2 for the location of the reverse reads file
- -s for the path to the singles reads 

***If desired, a list of kmers can be specified with the -k flag which will override automatic kmer selection.
<a name="ma"></a>
### **Assembly with MaSuRCA:**

This assembler is a combination of a De Bruijn graph and an Overlap-Layout-Consensus model. The Overlap-Layout-Consensus model consists of three steps, Overlap, which is the process of overlapping matching sequences in the data, this forms a long branched line. Layout, which is the process of picking the least branched line in from the overlap sequence created earlier, the final product here is called a contig. Consensus is the process of lining up all the contigs and picking out the most similar nucleotide line up in this set of sequences (OIRC). This assembly DOES NOT require a preprocessing step, such as Sickle, you will only input the raw data. For this assembly you will have another configuration file which can be found here https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/MaSuRCA%20Config%20file

(**Note**: Replace path to tutorial with the folder location)

**Running MaSuRCA:**

```
cd /Path_to_Tutorial/Assembly_Tutorial/Assembly/MaSuRCA

#run MaSuRCA

module load MaSuRCA/3.2.4

masurca config.txt

bash assemble.sh
```
<a name="quast"></a>
## Step 4: Assembly Statistics with QUAST

The final step for short read data is to analyze the quality of the assemblies. We will be using the program QUAST which will give us the number of contigs, total length and N50 value; the data we are most interested in. A good assembly would have small number of contigs, a total length that makes sense for the specific species, and a large N50 value. 

**Running QUAST:**

We run QUAST with the following commands:

```
module load quast/5.0.2

#SOAPdenovo statistics
quast.py /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SOAP/graph_Sample_*.scafSeq -o SOAP

#SPAdes statistics
quast.py /UCHC/PublicShare/Tutorials/Assembly_Tutorial/Assembly/SPAdes/scaffolds.fasta -o SPAdes

#MaSuRCA statistics
quast.py /UCHC/PublicShare/Tutorials/Assembly_Tutorial/samples/CA -o MaSuRCA
```

The following commands are located in the Sample_quast.sh file in the QUAST folder.

After running QUAST you will be able to access output files through two different processes. 
The first process would be to use an application like Cyberduck and pull the files from the transfer server to you home for viewing. 
The second would be to use pscp through the windows command prompt.

The statistics that are outputted via QUAST should follow this pattern.

|Info                    | MaSuRCA   | SOAPdenovo  |           |           | SPAdes    |
| -------------          | --------- | ----------  | --------- | --------- | --------- |
|                        |           |K-mer 31     |K-mer 35   |K-mer 41   |           |
| -------------          | --------- | ----------  | --------- | --------- | --------- |
|# contigs (>=0bp)       |116        |1507         |1905       |1486       |78         |
|# contigs (>= 1000bp)   |113        |249          |220        |198        |52         |
|# contigs (>= 5000bp)   |85         |             |           |           |           |
|# contigs (>=10000bp)   |73         |             |           |           |           |
|# contigs(>=25000bp)    |45         |             |           |           |           |
|# contigs (>=50000bp)   |18         |             |           |           |           |
|Total length (>=0bp)    |2871471    |3743924      |3764281    |3630629    |2885291    |
|Total length (>=1000bp) |2869528    |3554783      |3525490    |3426820    |2875160    |           
|Total length (>=5000bp) |2782331    |             |           |           |           |
|Total length (>=10000bp)|2696889    |             |           |           |           |
|Total length (>=25000bp)|2263199    |             |           |           |           |
|Total length (>=50000bp)|1389271    |             |           |           |           |
|# contigs               |115        |276          |246        |214        |59         |      
|Largest contig          |162425     |103125       |86844      |99593      |255551     |
|Total length            |2871084    |3574101      |3543834    |3438095    |2880184    |
|GC (%)                  |32.63      |32.44        |32.46      |32.46      |32.65      |
|N50                     |40374      |26176        |27766      |36169      |147660     |
|N75                     |27444      |14642        |16356      |16752      |54782      |
|L50                     |20         |44           |42         |33         |8          |
|L75                     |41         |91           |84         |69         |16         |
|# N's per 100 kbp       |0          |26547.43     |25459.35   |22602.08   |20.48      |

According to our requirements regarding n50 and contigs it would appear that the best assembly perfromed was via SPAdes.

## Step 5: Genomescope for Assessing Genome size

## Step 6: Read Alignment with Bowtie2
Bowtie2 takes read sequences and aligns them with long reference sequences. Since this is de novo assembly you will take the data from the assemblies you have and align them back to the raw read data. You want to use unpaired data. You will find the outputted data in the .err file, it should look like this:

|MaSuRCA                                   |
|------------------------------------------|
|894184 reads; of these:                   | 
|894184 (100.00%) were unpaired; of these: |
|125031 (13.98%) aligned 0 times           |
|748206 (83.67%) aligned exactly 1 time    |
|20947 (2.34%) aligned >1 times            |
|86.02% overall alignment rate             |

|SPAdes                                    |
|------------------------------------------|
|894184 reads; of these:                   | 
|894184 (100.00%) were unpaired; of these: |
|93720 (10.48%) aligned 0 times            |
|784462 (87.73%) aligned exactly 1 time    |
|16002 (1.79%) aligned >1 times            |
|89.52% overall alignment rate             |

|SOAP31                                    |
|------------------------------------------|
|894184 reads; of these:                   | 
|894184 (100.00%) were unpaired; of these: |
|456710 (51.08%) aligned 0 times           |
|437262 (48.90%) aligned exactly 1 time    |
|212 (0.02%) aligned >1 times              |
|48.92% overall alignment rate             |

|SOAP35                                    |
|------------------------------------------|
|894184 reads; of these:                   |
|894184 (100.00%) were unpaired; of these: |
|440065 (49.21%) aligned 0 times           |
|453682 (50.74%) aligned exactly 1 time    |
|437 (0.05%) aligned >1 times              |
|50.79% overall alignment rate             |

|SOAP41                                    |
|------------------------------------------|
|894184 reads; of these:                   |
|894184 (100.00%) were unpaired; of these: |
|406420 (45.45%) aligned 0 times           |
|486826 (54.44%) aligned exactly 1 time    |
|938 (0.10%) aligned >1 times              |
|54.55% overall alignment rate             |

## Step 7: Busco Evaluation
<a name="long"></a>
# Long Read Genome Assembly

<a name="gup"></a>
## Step 1: Base Calling with Guppy

<a name="ff"></a>
## Step 2: Assembly with Flye and Falcon

<a name="bus"></a>
## Step 3: Checking completeness with BUSCOMP

<a name="nano"></a>
## Step 4: Polishing with Nanopolish 

<a name="ph"></a>
## Step 5: Organizing with Purge Haplotigs

<a name="ha"></a>
# Hybrid Assembly 
