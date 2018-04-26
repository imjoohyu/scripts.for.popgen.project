#!/usr/bin/env python

#Randomly select 14 Dsim sequences V2 (python)
#Nov 11th, 2017
#Joo Hyun Im (ji72)

#read-in command from UNIX:
#python randomly_select_Dsim_14.py [file]

import sys
from Bio import SeqIO
import glob
import csv
import random

#1. Read in a sequence file (will be aligned aftewards)
#name_original = "/workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_FBgn0261514.fa"
name_original = sys.argv[1] #feed it in by unix loop
gene_id = name_original.split('/')[5]
gene_id = gene_id.split('_')[1]

print "gene ID: "
print gene_id

#2. Randomly pick 14 Dsim sequences and the top Dsim sequence (6/20/2017, 10:55am)
#Get the position/line number for each line
seq_num = []; i=0
with open(name_original,'r+') as data_original:
	for row in data_original:
		row_str = str(row)
		line_number=i
		if '>mel' in str(row): #Exclude the outgroup
			print "Dmel is excluded"
		elif '>' in str(row):
			seq_num.append(line_number)
		i+=1

#Randomly choose 14 Dsim positions and add the Dmel outgroup position to the list
seq_num_rand = random.sample(seq_num, 14)
seq_num_rand = [0] + seq_num_rand


#Add the position for the sequences
seq_num_rand_double = [j+1 for j in seq_num_rand] #get the sequence
seq_num_rand = seq_num_rand + seq_num_rand_double
seq_num_rand.sort()


#Get the sequences
selected_sequences = []; i=0
with open(name_original,'r+') as data_original:
	for row in data_original:
		if i in seq_num_rand:
			selected_sequences.append(row) 
		i+=1



#3. Save the data
output_address = '/workdir/ji72/DH/Dsim_prepped_for_DH_Nov_2017/Dsim_' + gene_id + '_14_selected.fa'
f = open(output_address, 'wb')
writer = csv.writer(f, delimiter='\n', quoting = csv.QUOTE_NONE, escapechar='\n')
writer.writerow(selected_sequences)
f.close()
