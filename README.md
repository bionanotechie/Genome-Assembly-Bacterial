# Genome Assembly Tutorial

## Table of Contents  
1. [Overview](#over)
2. [Short Read Genome Assembly](#short)  
   - [Copy the Assembly Directory to your account node](#copy)
   - [Quality Control with Sickle](#sickle)
   - [Assembly with SOAPdenovo, SPAdes, and MaSuRCA](#short-assemble)
     - [Assembly with SOAPdenovo](#soap)
     - [Assembly with SPAdes](#spades)
     - [Assembly with MaSuRCA](#ma)
   - [Assembly Statistics with QUAST](#quast)
   - [Assessing Genome size](#genome)
   - [Read Alignment with Bowtie2](#bow)
   - [Busco Evaluation](#bus)
3. [Long Read Genome Assembly](#long)
   - [Base Calling with Guppy](#gup)
   - [Assembly with Flye, Canu and Falcon](#ff)
     - [Assembly with Flye](#flye)
     - [Assembly with Canu](#canu)
     - [Assembly with Falcon](#falcon)
   - [Checking completeness with BUSCO](#bus)
   - [Polishing with Nanopolish](#nano)
   - [Organizing with Purge Haplotigs](#ph)
4. [Hybrid Assembly](#ha)

<a name="over"></a>
## Overview
  
This repository is a usable, publicly available tutorial. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.   

This tutorial will teach you how to use open source quality control, genome assembly, and assembly assessment tools to complete a high quality de novo assembly, this means that you do not have any data to compare back to. Moving through the tutorial you will take pair end short read data from a bacterial species and perform assemblies via various commonly used genome assmeblers. With these assemblies completed we will then need to assess the quality of the data produced. Once finished with the short read data we will move on to performing a long read assembly.

<a name="short"></a>
# Short Read Genome Assembly

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

This is a de novo assembler, this assembler, like MaSuRCA which we will be encountering later, requires a config file to run through the data. The configuration file can be found [here](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Short%20Read/sample_soap_config.cfg):
```
#maximal read length
max_rd_len=250
[LIB]
#average insert size
avg_ins=550
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#use only first 250 bps of each read
rd_len_cutoff=250
#in which order the reads are used while scaffolding
rank=1
# cutoff of pair number for a reliable connection (at least 3 for short insert size)
pair_num_cutoff=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
# path to genes
q1=/UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_1.fastq
q2=/UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_2.fastq
q=/UCHC/PublicShare/Tutorials/Assembly_Tutorial/Quality_Control/Sample_s.fastq
```
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

This assembler is a combination of a De Bruijn graph and an Overlap-Layout-Consensus model. The Overlap-Layout-Consensus model consists of three steps, Overlap, which is the process of overlapping matching sequences in the data, this forms a long branched line. Layout, which is the process of picking the least branched line in from the overlap sequence created earlier, the final product here is called a contig. Consensus is the process of lining up all the contigs and picking out the most similar nucleotide line up in this set of sequences (OIRC). This assembly DOES NOT require a preprocessing step, such as Sickle, you will only input the raw data. For this assembly you will have another configuration file which can be found [here](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Short%20Read/sample_masurca_config.cfg): 
```
# example configuration file 

# DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
# innies, i.e. --->.<---, and JUMP are assumed to be outties
# <---.--->. If there are any jump libraries that are innies, such as
# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
# are optional for PE libraries and mandatory for JUMP libraries. Any
# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
# converted into Celera Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
#Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
#if single-end, do not specify <reverse_reads>
#MUST HAVE Illumina paired end reads to use MaSuRCA
PE= pe 500 50  /FULL_PATH/frag_1.fastq  /FULL_PATH/frag_2.fastq
#Illumina mate pair reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#pacbio OR nanopore reads must be in a single fasta or fastq file with absolute path, can be gzipped
#if you have both types of reads supply them both as NANOPORE type
#PACBIO=/FULL_PATH/pacbio.fa
#NANOPORE=/FULL_PATH/nanopore.fa
#Other reads (Sanger, 454, etc) one frg file, concatenate your frg files into one if you have many
#OTHER=/FULL_PATH/file.frg
#synteny-assisted assembly, concatenate all reference genomes into one reference.fa; works for Illumina-only data
#REFERENCE=/FULL_PATH/nanopore.fa
END

PARAMETERS
#PLEASE READ all comments to essential parameters below, and set the parameters according to your project
#set this to 1 if your Illumina jumping library reads are shorter than 100bp
EXTEND_JUMP_READS=0
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)
USE_LINKING_MATES = 0
#specifies whether to run the assembly on the grid
USE_GRID=0
#specifies grid engine to use SGE or SLURM
GRID_ENGINE=SGE
#specifies queue (for SGE) or partition (for SLURM) to use when running on the grid MANDATORY
GRID_QUEUE=all.q
#batch size in the amount of long read sequence for each batch on the grid
GRID_BATCH_SIZE=500000000
#use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads
#can increase this to 30 or 35 if your reads are short (N50<7000bp)
LHE_COVERAGE=25
#set to 0 (default) to do two passes of mega-reads for slower, but higher quality assembly, otherwise set to 1
MEGA_READS_ONE_PASS=0
#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#CABOG ASSEMBLY ONLY: set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#CABOG ASSEMBLY ONLY: whether to attempt to close gaps in scaffolds with Illumina  or long read data
CLOSE_GAPS=1
#auto-detected number of cpus to use, set this to the number of CPUs/threads per node you will be using
NUM_THREADS = 16
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*20
JF_SIZE = 200000000
#ILLUMINA ONLY. Set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>=8Gbp) genomes from Illumina-only data
SOAP_ASSEMBLY=0
#Hybrid Illumina paired end + Nanopore/PacBio assembly ONLY.  Set this to 1 to use Flye assembler for final assembly of corrected mega-reads.  A lot faster than CABOG, at the expense of some contiguity. Works well even when MEGA_READS_ONE_PASS is set to 1.  DO NOT use if you have less than 15x coverage by long reads.
FLYE_ASSEMBLY=0
END 
```

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

<a name="genome"></a>
## Step 5: Assessing Genome size

You can learn how to asses the genome size by refering to this [tutorial](https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/)

Using SPAdes data you can expect the output to appear as such:


|Info                                 |             |
|-------------------------------------|-------------|
|Total k-mers                         |132587548    |
|Genome size estimation               |3233843      |
|Single copy region                   |2870240      |
|Single copy region/Genome estimation |0.8875633    |

![Image of Jelly](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Screen%20Shot%202019-07-30%20at%2011.01.39%20AM.png)

![Image of Jelly](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Screen%20Shot%202019-07-30%20at%2011.01.51%20AM.png)

![Image of Jelly](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Screen%20Shot%202019-07-30%20at%2011.02.16%20AM.png)

<a name="bow"></a>
## Step 6: Read Alignment with Bowtie2
Bowtie2 takes read sequences and aligns them with long reference sequences. Since this is de novo assembly you will take the data from the assemblies you have and align them back to the raw read data. You want to use unpaired data. 

The steps to perform a read alignment will be within the bowtie2.sh command. These steps include at first creating a combined fastq file  from your raw reads, to do this you use the command
``` 
Cat Sample_1.fastq Sample_2.fastq > genome.fastq
```
After this you will want to create directory where the all the indexes from your assemblies can be saved. After moving to that directory load bowtie2
```
Module load bowtie2/ 2.3.3.1
```
To create indexes from the assemblies you created you will need to use the command 
```
bowtie2-build PATH_TO_SCAFFOLD NAME_OF_SCAFFOLD.index
```
After this you will align the reads back to the combined raw read fastq file
```
bowtie2 -x NAME_OF_SCAFFOLD.index -U PATH_TO_COMBINED_RAWREAD.fastq -S NAME_OF_SCAFFOLD.bowtie2.sam
```


You will find the outputted data in the .err file, see the outputted results below. 

**Running Bowtie2**

You can run Bowtie2 by running [short_read_bowtie2.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Short%20Read/short_read_bowtie2.sh).

### Bowtie2 Results:
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

<a name="bus"></a>
## Step 7: BUSCO Evaluation
BUSCO stands for Benchmarking Universal Single-Copy Orthologs. This program assists with checking assemblies, annotations, and transcriptomes to see if they appear complete. It does this by taking an orthologous gene set of your species of interest and comparing it back to the genome of interest, taking into consideration possibly evolutionary changes.

BUSCO uses four flags, which are as follows for genome evaluation:
-i which will be your sequence file in fasta format
-l which is the lineage of your species of interest, in the tutorials case the bacterial lineage, these lineages can be obtained from the BUSCO website.
-o which will be the name of you would like for your output file
-m which will be the mode you plan to use BUSCO in, for this tutorial it is geno, which stands for genome

When preparing to run BUSCO you first need to have Augustus in your home directory, to do this use the command:
```
-cp avr /isg/shared/apps/augustus/3.2.3 your_directory
```

To run BUSCO use the command [short_read_busco.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Short%20Read/short_read_busco.sh)


Using the SPAdes data the BUSCO results look like:
```
C:98.6%[S:98.6%,D:0.0%],F:0.0%,M:1.4%,n:148

	146	Complete BUSCOs (C)
	146	Complete and single-copy BUSCOs (S)
	0	Complete and duplicated BUSCOs (D)
	0	Fragmented BUSCOs (F)
	2	Missing BUSCOs (M)
	148	Total BUSCO groups searched
  ```
  
<a name="long"></a>
# Long Read Genome Assembly
For long read assembly there is an additional step that is not used for short read data called **base calling**. This is performed first before any long read assembly. The process involves taking the data that the sequencer outputs, which appears as a squiggle line, and applying a base to the hills and valleys of the squiggle. For this step we will use the basecaller [Guppy](https://github.com/rrwick/Basecalling-comparison/blob/master/basecalling_scripts/guppy_basecalling.sh). To complete the assemblies, we will use the assemblers [Flye](https://github.com/fenderglass/Flye), [Canu](https://canu.readthedocs.io/en/latest/), and [Falcon](https://pb-falcon.readthedocs.io/en/latest/). [Purge Haplotags](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) will be used to assure that the contigs that are assembled are not being combined with the Haplotig of that sequence. After this, the assembly will be polished via [Nanopolish](https://github.com/nanoporetech/nanopolish). Assessment of the quality of the genome assembled will be completed through the same program as with the short read data, QUAST.


<a name="gup"></a>
## Step 1: Base Calling with Guppy
As mentioned earlier this is required as the data that is outputted via the sequencer is different from those outputted via a short read sequencer. Guppy runs via a neural network and can work off GPUs and CPUs.

**Running Guppy:**

To run guppy, run guppy.sh located in _______.

The output files of guppy will be:
- guppy_basecaller_log-2019-06-26_14-07-05.log 
- sequencing_telemetry.js 
- sequencing_summary.txt

<a name="ff"></a>
## Step 2: Assembly with Flye, Canu and Falcon
In this step we will run all of the Long read genome assemblers on the basecalled data.

<a name="flye"></a>
### Assembly with Flye
This assembler takes data from Pacbio or Oxford Nanopore technologies sequencers and outputs polished contigs. It will repeat graph, that is similar in appearance to the De Bruijn graph. The manner in which this graph is assembled reveals the repeats in the genome allowing for the most accurate assembly. 

**Running Flye**

To run Flye run [flye.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/flye.sh) located in _______.

The flags in flye are:
- --pacbio-raw - fasta file
- --out-dir - location of result
- --genome-size
- --threads

<a name="canu"></a>
### Assembly with Canu
Canu assembles the with the same data as Flye, in a three step process, these steps are: Correction, trimming, and unitig construction.


**Running Canu**

To run Flye run canu.sh located in _______.

The flags in Canu are:
-p which applies to the prefix of the intermediate and output file names. 
-d specifies the directory
-s imports parameters from the specification file. 


<a name="falcon"></a>
### Assembly with Falcon

Falcon is another de novo assembler which is used for PacBio Long Read data. 

Falcon has 3 inputs:
- your PacBio data in fasta format (can be one or many files), 
- a text file telling FALCON where to find your fasta files,
- and a configuration file

**Running Falcon**

To run Flye run falcon.sh located in _______.


<a name="bus"></a>

## Step___: Checking for errors with Centrifuge

## Step 3: Checking completeness with BUSCO
BUSCO was discussed earlier during the short read tutorial, here we will use it to assess the genome before and after polishing. Which was described earlier during the short read assembly.
You can run BUSCO with [long_read_BUSCO.sh] (https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Long_read_busco.sh)

BUSCO with flye before polishing:
```
C:56.1%[S:48.4%,D:7.7%],F:14.4%,M:29.5%,n:430

	241	Complete BUSCOs (C)
	208	Complete and single-copy BUSCOs (S)
	33	Complete and duplicated BUSCOs (D)
	62	Fragmented BUSCOs (F)
	127	Missing BUSCOs (M)
	430	Total BUSCO groups searched
```

BUSCO with flye after polishing:
```
C:82.8%[S:71.4%,D:11.4%],F:4.7%,M:12.5%,n:430

	356	Complete BUSCOs (C)
	307	Complete and single-copy BUSCOs (S)
	49	Complete and duplicated BUSCOs (D)
	20	Fragmented BUSCOs (F)
	54	Missing BUSCOs (M)
	430	Total BUSCO groups searched
```
<a name="nano"></a>
## Step 4: Polishing with Nanopolish 
Nanopolish is used to strengthen consensus data from yoru assembly.It will take the assembly you have created and align it, break it into segments, and then a consensus algorithm can run through the segments to polish them.

<a name="ph"></a>
## Step 5: Organizing with Purge Haplotigs
Purge Haplotigs assures that there is not a combination of sequences between contigs and haplotigs. It uses a system that uses the mapped reads that you assembled and Minimap2 to assess which contigs should be kept in the assembly.

## 

<a name="ha"></a>
# Hybrid Assembly 
To perform a hybrid assembly it requires you have both short and long read data to complete the genome.

## Step 1:

## Step 2: Assembly with Kraken and MaSuRCA
