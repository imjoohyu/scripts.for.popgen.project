#!/usr/bin/env python

#Prepare the aligned sequence to have the outgroup sequence on top (python)
#June 6th, 2017
#Joo Hyun Im (ji72)

#read-in command from UNIX:
#python put_outgroup_on_top_DH.py [file]

import sys
from Bio import SeqIO
import glob
import csv


#0. Tell me which species is the main species
#indicator = 'Dmel'
indicator = 'Dsim'


#1. Read in an aligned sequence file
#name_original = "Dsim_FBgn0261514.fa_aligned.best.fas_aligned_confirmed.fasta"
name_original = sys.argv[1] #feed it in by unix loop
gene_id = name_original.split('/')[6]
gene_id = gene_id.split('_')[1]

print gene_id

#2. Remove the outgroup sequence sequence

#Get the number for the outgroup

outgroup = 'mel'
if indicator == 'Dmel':
	outgroup = 'sim'
if indicator == 'Dsim':
	outgroup = 'mel'

with open(name_original, 'r+') as data_original:
	for num, line in enumerate(data_original, 1):
		has_match = line.find(outgroup)
		#print has_match
		#print num
		if has_match !=-1:
			outgroup_num = num


lines = [outgroup_num-1,outgroup_num]
i=0; outgroups_removed = [] 
with open(name_original,'r+') as data_original:
	for row in data_original:
		if not i in lines:
			outgroups_removed.append(row) 
		i+=1


#Save the outgroup separately
lines = [outgroup_num-1,outgroup_num]; outgroup_only = []; i=0 
with open(name_original,'r+') as data_original:
	for row in data_original:
		if i in lines:
			outgroup_only.append(row) 
		i+=1


#3. Add them together
merged_list = outgroup_only + outgroups_removed


#4. Save the sequence
if indicator == 'Dmel':
	output_address = '/workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/Dmel_' + gene_id + '_149_selected_aligned_confirmed_ordered.fasta' #for matching a Dsim to Dmels
if indicator == 'Dsim':
	output_address = '/workdir/ji72/DH/Dsim_prepped_for_DH_Nov_2017/Dsim_14_selected/Dmel_' + gene_id + '_14_selected_aligned_confirmed_ordered.fasta' #for matching a Dmel to Dsims

f = open(output_address, 'wb')
writer = csv.writer(f, delimiter='\n', quoting = csv.QUOTE_NONE, escapechar='\n')
writer.writerow(merged_list)
f.close()


