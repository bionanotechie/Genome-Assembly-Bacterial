#! /usr/bin/env python

import sys

from Bio import SeqIO

fasta_file = "assembly.fasta"

fasta1 = "flye_moss_rm3kb_0-10k.fasta"

fasta2 = "flye_moss_rm3kb_10k-20k.fasta"

fasta3 = "flye_moss_rm3kb_20k-30k.fasta"

fasta4 = "flye_moss_rm3kb_30k-40k.fasta"

fasta5 = "flye_moss_rm3kb_40k-47k.fasta"

with open(fasta1, "w") as a, open(fasta2, "w") as b, open(fasta3, "w") as c, open(fasta4, "w") as d, open(fasta5, "w") as e:

   for seq in SeqIO.parse(fasta_file, 'fasta'):

       if seq.id[0] is 'c':

           if int(seq.id[7:]) > 0 and int(seq.id[7:]) <= 10000:

               SeqIO.write(seq, a, "fasta")

           elif int(seq.id[7:]) > 10000 and int(seq.id[7:]) <= 20000:

               SeqIO.write(seq, b, "fasta")

           elif int(seq.id[7:]) > 20000 and int(seq.id[7:]) <= 30000:

               SeqIO.write(seq, c, "fasta")

           elif int(seq.id[7:]) > 30000 and int(seq.id[7:]) <= 40000:

               SeqIO.write(seq, d, "fasta")

           else:

               SeqIO.write(seq, e, "fasta")

       else:

           if int(seq.id[9:]) > 0 and int(seq.id[9:]) <= 10000:

               SeqIO.write(seq, a, "fasta")

           elif int(seq.id[9:]) > 10000 and int(seq.id[9:]) <= 20000:

               SeqIO.write(seq, b, "fasta")

           elif int(seq.id[9:]) > 20000 and int(seq.id[9:]) <= 30000:

               SeqIO.write(seq, c, "fasta")

           elif int(seq.id[9:]) > 30000 and int(seq.id[9:]) <= 40000:

               SeqIO.write(seq, d, "fasta")

           else:

               SeqIO.write(seq, e, "fasta")
