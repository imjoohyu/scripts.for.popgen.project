#!/usr/bin/env python

#Prepare the aligned sequence to have the outgroup sequence on top (python)
#Nov 21st, 2017
#Joo Hyun Im (ji72)

#read-in command from UNIX:
#python put_outgroup_on_top_MK_Nov_2017.py [file]

import sys
from Bio import SeqIO
import glob
import csv
from collections import Counter

#0. Tell me which species is the main species
indicator = 'Dmel'
#indicator = 'Dsim'


#1. Read in an aligned sequence file
#name_original = "Dmel_FBgn0261514_149_selected.fa_aligned.best.fas_confirmed.fa"
name_original = sys.argv[1] #feed it in by unix loop
gene_id = name_original.split('/')[5]
gene_id = gene_id.split('_')[1]

#print gene_id

#2. Remove the outgroup sequence sequence

#Get the number for the outgroup

outgroup = 'mel'
if indicator == 'Dmel':
	outgroup = 'sim'
if indicator == 'Dsim':
	outgroup = 'mel'
outgroup_2 = 'yak'


with open(name_original, 'r+') as data_original:
	for num, line in enumerate(data_original, 1):
		has_match = line.find(outgroup)
		has_match_yak = line.find(outgroup_2)
		#print has_match
		#print num
		if has_match !=-1:
			outgroup_num = num
		if has_match_yak !=-1:
			outgroup_yak_num = num


#Save the non-outgroup sequences
lines = [outgroup_num-1,outgroup_num, outgroup_yak_num-1, outgroup_yak_num]
i=0; outgroups_removed = [] 
with open(name_original,'r+') as data_original:
	for row in data_original:
		if not i in lines:
			outgroups_removed.append(row) 
		i+=1


#Save the outgroup separately
lines = [outgroup_num-1, outgroup_num, outgroup_yak_num-1, outgroup_yak_num]; outgroup_only = []; i=0 
with open(name_original,'r+') as data_original:
	for row in data_original:
		if i in lines:
			outgroup_only.append(row) 
		i+=1

#remove multiple ';'s from >yak and >sim
outgroup_only_name_fixed = []
for row in outgroup_only:
	num_of_semicolon = row.count(';')
	#print num_of_semicolon
	if num_of_semicolon > 1:
		row = row.split(';')[0]
		row = row + ';'
		outgroup_only_name_fixed.append(row)
	if not num_of_semicolon > 1:
		outgroup_only_name_fixed.append(row)


#3. Add them together
merged_list = outgroup_only + outgroups_removed

#4. Remove the space ('\n') from each line and save the sequence
for row in merged_list:
	row_stripped = row.split('\n')[0]
	print row_stripped

