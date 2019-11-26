# Genome Assembly Tutorial

This repository is a usable, publicly available tutorial. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**. 

## Table of Contents  
1. [Overview](#over)
2. [Short Read Genome Assembly](#short)  
   - [Copy the Assembly Directory to your account node](#copy)
   - [Quality Control with Sickle](#sickle)
   - [Assembly with SOAPdenovo, SPAdes, and MaSuRCA](#short-assemble)
     - [Assembly with SOAPdenovo](#soap)
     - [Assembly with SPAdes](#spades)
     - [Assembly with MaSuRCA](#ma)
   - [Assessing Genome size](#genome)
   - [Assembly Statistics with QUAST](#quast)
   - [Read Alignment with Bowtie2](#bow)
   - [Busco Evaluation](#bus)
3. [Long Read Genome Assembly](#long)
   - [Base Calling with Guppy](#gup)
   - [Assembly with Flye Shasta, and Falcon](#ff)
     - [Assembly with Flye](#flye)
     - [Assembly with Shasta](#shas)
     - [Assembly with Falcon](#falcon)
   - [Checking completeness with BUSCO](#bus2)
   - [Polishing with Nanopolish](#nano)
   - [Organizing with Purge Haplotigs](#ph)
4. [Hybrid Assembly](#ha)
   - [Preprocessing with CCS](#ccs)
   - [Assembly with MaSuRCA](#mas)

<a name="over"></a>
## Overview  

This tutorial will teach you how to use open source quality control, genome assembly, and assembly assessment tools to complete a high quality de novo assembly which is commonly utilized when you dont have a reference genome. Moving through the tutorial, you will take pair end short read data from a bacterial species and perform assemblies via various commonly used genome assmeblers. With these assemblies completed we will then need to assess the quality of the data produced. Once finished with the short read data we will move on to performing a long read assembly with long read nanopore data using basecalling and commonly used long read assemblers. Finally, we will then move on to Hybrid PacBio data. 

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
- -K is the size of the k-mer, a k-mer is a set of nucleotides, k is the number of nucleotides in that set. It is a crucial parameter in most de Brujin Graph assemblers and assemblers work with the highest accuracy if the k-mer size estimation is accurate. 
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

<a name="genome"></a>
## Step 4: Assessing Genome size

You can learn how to asses the genome size by refering to this [tutorial](https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/)

Using the raw data you can expect the output to appear as such:


|Info                                 |             |
|-------------------------------------|-------------|
|Total k-mers                         |132587548    |
|Genome size estimation               |3233843      |
|Single copy region                   |2870240      |
|Single copy region/Genome estimation |0.8875633    |

![Image of Jelly](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Screen%20Shot%202019-07-30%20at%2011.01.39%20AM.png)

![Image of Jelly](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Screen%20Shot%202019-07-30%20at%2011.01.51%20AM.png)

![Image of Jelly](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Screen%20Shot%202019-07-30%20at%2011.02.16%20AM.png)

<a name="quast"></a>
## Step 5: Assembly Statistics with QUAST

The final step for short read data is to analyze the quality of the assemblies. We will be using the program QUAST which will give us the number of contigs, total length and N50 value; the data we are most interested in. A good assembly would have small number of contigs, a total length that makes sense for the specific species, and a large N50 value. N50 is a measure to describe the quality of assembled genomes fragmented in contigs of different length. The N50 is the minimum contig length needed to cover 50% of the genome.

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


<a name="bow"></a>
## Step 6: Read Alignment with Bowtie2
Bowtie2 is a tool you would use for comparitive genomics via alignment. Alignment is the process where we discover how and where the read sequences are similar to a reference sequence. An alignment is a way of lining up the characters in the read with some characters from the reference to reveal how they are similar. 

Bowtie2 in our case takes read sequences and aligns them with long reference sequences. Since this is de novo assembly you will take the data from the assemblies you have and align them back to the raw read data. You want to use unpaired data. 

You will find the outputted data in the .err file, see the outputted results below. 

**Running Bowtie2**

You can run Bowtie2 by running [short_read_bowtie2.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Short%20Read/short_read_bowtie2.sh).

### Bowtie2 Results:
|MaSuRCA                                   |
|------------------------------------------|
|447092 reads; of these:                   | 
|447092 (100.00%) were paired; of these:   |
|189344 (42.35%) aligned 0 times           |
|255081 (57.05%) aligned exactly 1 time    |
|2667 (0.60%) aligned  >1 times            |
|88.51% overall alignment rate             |

|SPAdes                                    |
|------------------------------------------|
|447092 reads; of these:                   | 
|447092 (100.00%) were paired; of these    |
|93720 (10.48%) aligned 0 times            |
|267852 (59.91%) aligned exactly 1 time    |
|3681 (0.82%) aligned  >1 times            |
|89.52% overall alignment rate             |

|SOAP31                                    |
|------------------------------------------|
|447092 reads; of these:                   | 
|447092 (100.00%) were paired; of these:   |
|338885 (75.80%) aligned 0 times           |
|108205 (24.20%) aligned exactly 1 time    |
|2 (0.00%) aligned >1 times                |
|48.96% overall alignment rate             |

|SOAP35                                    |
|------------------------------------------|
|447092 reads; of these:                   |
|447092 (100.00%) were paired; of these:   |
|332703 (74.41%) aligned 0 times           |
|114376 (25.58%) aligned exactly 1 time    |
|13 (0.00%) aligned  >1 times              |
|50.82% overall alignment rate             

|SOAP41                                    |
|------------------------------------------|
|447092 reads; of these:                   |
|447092 (100.00%) were paired; of these:   |
|318539 (71.25%) aligned 0 times           |
|128499 (28.74%) aligned exactly 1 time    |
|54 (0.01%) aligned  >1 times              |
|54.57% overall alignment rate             |

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
  Using the MaSuRCA data the BUSCO results look like:
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
For long read assembly there is an additional step that is not used for short read data called **base calling**. This is performed first before any long read assembly. The process involves taking the data that the sequencer outputs, which appears as a squiggle line, and applying a base to the hills and valleys of the squiggle. For this step we will use the basecaller [Guppy](https://github.com/rrwick/Basecalling-comparison/blob/master/basecalling_scripts/guppy_basecalling.sh). To complete the assemblies, we will use the assemblers [Flye](https://github.com/fenderglass/Flye), [Shasta](https://chanzuckerberg.github.io/shasta/), and [Falcon](https://pb-falcon.readthedocs.io/en/latest/). [Purge Haplotags](https://bitbucket.org/mroachawri/purge_haplotigs/src/master/) will be used to assure that the contigs that are assembled are not being combined with the Haplotig of that sequence. After this, the assembly will be polished via [Nanopolish](https://github.com/nanoporetech/nanopolish). Assessment of the quality of the genome assembled will be completed through the same program as with the short read data, QUAST.


<a name="gup"></a>
## Step 1: Base Calling with Guppy
Oxford Nanopore long-read sequencing works by doing the following: pass a single strand of DNA through a membrane with a nanopore and apply a voltage difference across the membrane. By doing so the nucleotides present in the pore will affect the pore’s electrical resistance so current measurements over time indicate the sequence of DNA bases passing through the pore. This  current signal is the raw data gathered by an ONT sequencer. 

Basecalling for ONT devices is the process of translating this raw signal into a DNA sequence. Guppy runs via a neural network and can work off GPUs and CPUs.

**Running Guppy:**

To run guppy, run guppy.sh located in the long read directory.

The output files of guppy will be:
- guppy_basecaller_log-2019-06-26_14-07-05.log 
- sequencing_telemetry.js 
- sequencing_summary.txt

<a name="ff"></a>
## Step 2: Assembly with Flye, Shasta and Falcon 
In this step we will run all of the Long read genome assemblers on the basecalled data.

<a name="flye"></a>
### Assembly with Flye
This assembler takes data from Pacbio or Oxford Nanopore technologies sequencers and outputs polished contigs. It will repeat graph, that is similar in appearance to the De Bruijn graph. The manner in which this graph is assembled reveals the repeats in the genome allowing for the most accurate assembly. 

**Running Flye**

To run Flye run [flye.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/flye.sh) located in the long read assembly folder.

The flags in flye are:
- --pacbio-raw - fasta file
- --out-dir - location of result
- --genome-size
- --threads

<a name="canu"></a>
### Assembly with Shasta
Similar to Flye, the Shasta long read assemblers purpose is to rapidly produce an accurate assembled sequence using Oxford Nanopoore sequencing data.

The computational methods used in the Shasta assembler are:
- A run-length representation of the read sequence which makes the assembly process more resilient to errors in homopolymer repeat counts (very commmon in nanopore data).
- Using a fixed subset of kmers in some phases of the computation a representation of the read sequence based on markers.


**Running Shasta**

To run Shasta run [shasta.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/shasta.sh) located in the long read assembly folder.

The flags in Shasta are:
- -p which applies to the prefix of the intermediate and output file names. 
- -d specifies the directory
- -s imports parameters from the specification file. 

<a name="falcon"></a>
### Assembly with Falcon

Falcon is another de novo assembler which is used for PacBio Long Read data. 

Falcon has 3 inputs:
- your PacBio data in fasta format (can be one or many files), 
- a text file telling FALCON where to find your fasta files,
- and a configuration file

**Running Falcon**

To run Flye run falcon.sh located in the hybrid assembly directory.

<a name="bus2"></a>
## Step 3: Checking completeness with BUSCO
BUSCO was discussed earlier during the short read tutorial, here we will use it to assess the genome before and after polishing. Which was described earlier during the short read assembly.
You can run BUSCO with [long_read_BUSCO.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Long_read_busco.sh)

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

BUSCO with flye after Purge Haplotigs:
``` 
C:81.1%[S:72.3%,D:8.8%],F:5.1%,M:13.8%,n:430

	349    Complete BUSCOs (C)
	311    Complete and single-copy BUSCOs (S)
	38     Complete and duplicated BUSCOs (D)
	22     Fragmented BUSCOs (F)
	59     Missing BUSCOs (M)
	430    Total BUSCO groups searched
```
<a name="nano"></a>
## Step 4: Polishing with Nanopolish 
Nanopolish is used to strengthen consensus data from your assembly.It will take the assembly you have created and align it, break it into segments, and then a consensus algorithm can run through the segments to polish them.

The original purpose of nanopolish was to improve the consensus accuracy of an assembly of Oxford Nanopore Technology sequencing reads. 

**Running Nanopolish**
To run nanopolish run the [nanopolish0-10kb.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/Long%20Read/3-nanopolish/nanopolish.sh) file located in the long read folder inside folder 3.

In our script we first run nanopolish_makerange.py in order to split the draft of larger genomes so that the algorithm can run in parallel on each part. In our case we first run the divide_genome.py script on our full genome assembly (whichever you want to access, we used the flye output). then we run nanopolish_0-10kb.sh and just run nanopolish on that portion.

Here are the following meanings of the parameters:
- -r is the input of the raw reads
- -b is the bam file of the raw reads
- -g is the partial genome you want to examine
- -o is the output (a vcf file)

<a name="ph"></a>
## Step 5: Organizing with Purge Haplotigs
Purge Haplotigs is a pipeline to help with curating genome assemblies. It assures that there is not a combination of sequences between contigs and haplotigs. It uses a system that uses the mapped reads that you assembled and Minimap2 to assess which contigs should be kept in the assembly.

We use this because some parts of a genome may have a very high degree of heterozygosity which causes contigs for both haplotypes of that part of the genome to be assembled as separate primary contigs, rather than as a contig with a haplotig which may cause an error in analysis.

What the purge halpotigs algortihm does is identify pairs of contigs that are syntenic and group them as haplotig. The pipeline uses mapped read coverage and blast/lastz alignments to determine which contigs to keep for the assembly. Dotplots are produced for all flagged contig matches to help the user analyze any remaining ambiguous contigs. 

This part must be done in seperate steps as the parameters in each part depend on the results of the previous steps.

**Running Purge Haplotigs**

To run purge haplotigs, you mut run each script seperately. you must run [purge_haplotigs_step1.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/purge_haplotigs_step1.sh), [purge_haplotigs_step2.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/purge_haplotigs_step2.sh), and [purge_haplotigs_step3.sh](https://github.com/CBC-UCONN/Genome-Assembly-Bacterial/blob/master/purge_haplotigs_step3.sh) in the long read folder.


<a name="ha"></a>
# Hybrid Assembly 

For the long read assembly section, we have been working with long read Nanopore data. In this section we will be working with hybrid assemblers which will be compatible with long read PacBio Data and Short read Illumina data. Nanopore and PacBio are currently both the main long read sequencing technologies but the major differences in them are that PacBio reads a molecule multiple times to generate high-quality consensus data while Nanopore can only sequence a molecule twice. As a result, PacBio generates data with lower error rates compared to Oxford Nanopore. 

To perform a hybrid assembly it requires you have both short and long read data to complete the genome. Hybrid assembly uses short read data to resolve ambiguity in the long read data as it is assembled. For this tutorial we are using data from a boxelder genome. We will begin with preprocessing with PacBio Circular Consensus Sequence analysis application (CCS), run the Falcon, Kraken, and MaSuRCA assemblers,

<a name="ccs"></a>
## Step 1:Preprocessing with PacBio CCS
CCS takes multiple subreads of the same molecule and combines them using a statistical model to produce one accurate consensus sequence (HiFi read), with base quality values. For more information, refer to the [PacBio Website](https://www.pacb.com/smrt-science/smrt-sequencing/smrt-sequencing-modes/)

For our data, CCS was already run on the boxelder data. For more information on CCS please refer to [this PacBio CCS Tutorial](https://www.pacb.com/videos/tutorial-circular-consensus-sequence-analysis-application/) And if you have a UConn PacBio account please refer to [this tutorial](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/pacbio-v7/).






<a name="fkm"></a>
## Step 2: Assembly with MaSuRCA


