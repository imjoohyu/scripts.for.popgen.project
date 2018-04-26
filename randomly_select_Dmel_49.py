#!/usr/bin/env python

#Randomly select 49 Dmel sequences (python)
#June 19th, 2017
#Joo Hyun Im (ji72)

#read-in command from UNIX:
#python randomly_select_Dmel_49.py [file]

import sys
from Bio import SeqIO
import glob
import csv
import random

#1. Read in an aligned sequence file
#name_original = "Dsim_FBgn0000146.fa_aligned.best.fas"
#name_original = "Dsim_FBgn0261514.fa_aligned.best.fas_aligned_confirmed.fasta"
#name_original = "/workdir/ji72/DH/Dmel_prepped_for_DH_done/Dmel_FBgn0261514.fa_aligned_confirmed_ordered.fasta"
name_original = sys.argv[1] #feed it in by unix loop
gene_id = name_original.split('/')[5]
gene_id = gene_id.split('_')[1]

#2. Randomly pick 49 sequences -- (6/19/2017, 6:33pm)
#Pick the first 50 sequences (Dmel sequences were already "randomized" after alignment)
# selected_sequences = []; i = 0; count_list = []; count = None
# with open(name_original,'r+') as data_original:
# 	for row in data_original:
# 		if i<445: #first 50
# 			row_str = str(row)
# 			line_number = i
# 			if '>' in str(row): #only choose the header
# 				selected_sequences.append(row_str)
# 				count =i+3 #header, two newlines, and then sequence
# 				count_list.append(count)
# 			elif count in count_list and line_number == count:
# 				selected_sequences.append(row_str)
# 			i+=1

#Truly randomly pick 49 Dmel sequences and the top Dsim sequence (6/20/2017, 10:55am)
#Get the position/line number for each line
seq_num = []; i=0
with open(name_original,'r+') as data_original:
	for row in data_original:
		row_str = str(row)
		line_number=i
		if '>sim' in str(row): #Exclude the outgroup
			print "Dsim is excluded"
		elif '>' in str(row):
			seq_num.append(line_number)
		i+=1

#Randomly choose 49 Dmel positions and add the Dsim outgroup position to the list
seq_num_rand = random.sample(seq_num, 49)
seq_num_rand = [0] + seq_num_rand


#Add the position for the sequences
seq_num_rand_double = [j+3 for j in seq_num_rand] #get the sequence
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
output_address = '/workdir/ji72/DH/Dmel_prepped_for_DH_done/Dmel_' + gene_id + '_aligned_confirmed_ordered_49_selected.fasta'
f = open(output_address, 'wb')
writer = csv.writer(f, delimiter='\n', quoting = csv.QUOTE_NONE, escapechar='\n')
writer.writerow(selected_sequences)
f.close()
