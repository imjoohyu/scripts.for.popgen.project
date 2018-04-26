#!/usr/bin/env python

#Match the Dmel and Dsim sequences for DH run(python)
#June 6th, 2017
#Joo Hyun Im (ji72)

#read-in command from UNIX:
#python matching.py FBgn0261514_nimA_Chr2L_13957898_13962981.fa.out_pasted_excluding_unusable_lines.fa

import sys
from Bio import SeqIO
import glob
import csv


#0. Tell me which species is the first species
indicator = 'Dmel'
#indicator = 'Dsim'

#1. Pick a file in the Dmel folder and read the gene ID off of it (split the file name by '_' and take the first part [0]

##for instance: name_in_the_first_folder = 'FBgn0261514_nimA_Chr2L_13957898_13962981.fa.out_pasted_excluding_unusable_lines.fa'
#or '/workdir/ji72/Nexus_CDS_pasted_Ns_removed/FBgn0261514_nimA_Chr2L_13957898_13962981.fa.out_pasted_excluding_unusable_lines.fa' in workdir (CBSU)
name_in_the_first_folder = sys.argv[1] #feed it in by unix loop
gene_id_in_the_first_folder = name_in_the_first_folder.split('/')[4]
gene_id_in_the_first_folder = gene_id_in_the_first_folder.split('_')[0]
data_in_the_first_folder = open(name_in_the_first_folder,'r')

#data_in_the_first_folder.read() #checking that the data gets read in


#2. Look for a file that contains that gene ID in the second folder.
#3. Open the file from the second folder and copy the second sequence (the first sequence the sim.line00_reference_haploid.fa sometimes has gaps, so it's safer to use the second sequence sim.line01_MD221_haploid.fa).

if indicator == 'Dmel':
	second_folder_address = '/workdir/ji72/Rogers_CDS_pasted_Ns_removed/' + gene_id_in_the_first_folder + '*' #for matching Dsim to Dmels
if indicator == 'Dsim':
	second_folder_address = '/workdir/ji72/Nexus_CDS_pasted_Ns_removed/' + gene_id_in_the_first_folder + '*' #for matching Dmel to Dsims
#second_folder_address = '/home/ji72/popgen/Dsim_Rogers/Rogers_CDS_pasted_Ns_removed/' + gene_id_in_the_first_folder + '*'


name_in_the_second_folder = None
for name_in_the_second_folder in glob.glob(second_folder_address):
    print name_in_the_second_folder

#If the matching case is not available in the second folder:
if name_in_the_second_folder == None:
	sys.exit("This gene, " + gene_id_in_the_first_folder +  ", does not have the equivalent gene in the second folder.")

#Pick an outgroup 
lines = [2,3]; i=0; chosen_data_in_the_second_folder = [] #Pick line 3 and 4
with open(name_in_the_second_folder,'r') as data_in_the_second_folder:
	reader = csv.reader(data_in_the_second_folder, delimiter='\n')
	for row in reader:
		if i in lines:
			#print row
			chosen_data_in_the_second_folder.append(row) 
		i+=1	

#change 
chosen_header = ''.join(map(str, chosen_data_in_the_second_folder[0]))
chosen_seq = ''.join(map(str, chosen_data_in_the_second_folder[1]))



#4. Concatenate that sequence to the first file.
put_together_both_sequences = [chosen_data_in_the_second_folder[0],chosen_data_in_the_second_folder[1]]
with open(name_in_the_first_folder,'r') as data_in_the_first_folder:
	reader = csv.reader(data_in_the_first_folder, delimiter='\n')
	for row in reader:
		put_together_both_sequences.append(row)


#5. Change the name of each sequence to be '>mel;' and then change the name of the first sequence to be '>sim;'. Make the output as a new file.
i=0; first_i = 0
for line in put_together_both_sequences:
	if '>' in str(line): #if the line has '>', replace it with the new name
		line_name = str(line).split('>')[1]
		line_name = line_name.split("']")[0]
		if i == first_i:
			if indicator == 'Dmel':
				line_name_changed = '>sim; ' + line_name #for matching a Dsim to Dmels
			if indicator == 'Dsim':
				line_name_changed = '>mel; ' + line_name #for matching a Dmel to Dsims
		else:
			if indicator == 'Dmel':
				line_name_changed = '>mel; ' + line_name #for matching a Dsim to Dmels
			if indicator == 'Dsim':
				line_name_changed = '>sim; ' + line_name #for matching a Dmel to Dsims
			
		put_together_both_sequences[i] = [line_name_changed]
	i+=1



#6. Save the sequence
if indicator == 'Dmel':
	output_address = '/workdir/ji72/Dmel_prepped_for_DH/Dmel_' + gene_id_in_the_first_folder + '.fa' #for matching a Dsim to Dmels
if indicator == 'Dsim':
	output_address = '/workdir/ji72/Dsim_prepped_for_DH/Dsim_' + gene_id_in_the_first_folder + '.fa' #for matching a Dmel to Dsims
print output_address

f = open(output_address, 'wb')
writer = csv.writer(f, delimiter='\n', quoting = csv.QUOTE_NONE, escapechar='\n')
writer.writerows(put_together_both_sequences)
f.close()


