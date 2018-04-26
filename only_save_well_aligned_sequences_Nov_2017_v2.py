#!/usr/bin/env python

#Remove sequences of fly lines that are mal-aligned (Ns at all positions) V2 (with more strict filter)
#Nov 14th 2017
#Joo Hyun Im (ji72)

import sys
import csv

from collections import Counter
from Bio import SeqIO

#file_name = "/workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/Dmel_FBgn0000146.fa_149_selected.fa_aligned.best.fas"
file_name = sys.argv[1]
gene_id = file_name.split('/')[6]
gene_id = gene_id.split('_')[1]

#Ideally, every sequence should have less than 20% of their sequence covered with '-' 
#So the count should be 150 for Dmel and 15 for Dsim

count = 1
with open(file_name,'r') as data:
    for record in SeqIO.parse(data, "fasta"):
        sequence = str(record.seq) #individual sequence
        if sequence.count('-') < len(sequence)*0.2: #good alignment only
            count = count + 1

#Change the number of cutoff depending on the species
#Dmel_cutoff = 149, Dsim_cutoff = 14

cutoff = 149

if count > cutoff: #if the number of lines without alignment problem is more than half
    with open(file_name,'r') as data:
        for record in SeqIO.parse(data, "fasta"):
            sequence = str(record.seq) #individual sequence
    
            #good alignment only
            if sequence.count('-') < len(sequence)*0.2:
                print str('>') + str(record.id)
                print sequence

else:
    sys.exit("This gene, " + gene_id +  ", has a potential alignment problem.")
