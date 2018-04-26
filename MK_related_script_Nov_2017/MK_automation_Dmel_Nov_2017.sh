#! /bin/bash

#MK_automation_Dmel_Nov_2017.sh
#Prepare data to run MK test
#Nov 19th, 2017

#Starting materials:
#FBgn00XXXXX_Dsim_Dyak_matched_to_Dmels.fa which contains one Dsim sequence, one Dyak sequence and 20 Dmel sequences

#Make sure the following are present in the /workdir folder before starting this run:
#remove_lines_with_missing_data_Nov_2017_MK.py
#randomly_select_Dmel_149_Nov_2017_MK.py
#Move prank to the /workdir: cp -r /home/ji72/prank-msa/ /workdir/ji72/ | cd /workdir/ji72/ | export PATH=$PATH:/workdir/ji72/prank-msa/src
#only_save_well_aligned_sequences_Nov_2017_MK.py -- change the setting to the right species


#1. Remove lines from each gene sequence that have missing data

#Input: FBgn0261514_Dsim_Dyak_matched_to_Dmels.fa

for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/*.fa; do
     python remove_lines_with_missing_data_Nov_2017.py $f >`echo $f`_filtered.fasta
done

#Output: FBgn0261514_Dsim_Dyak_matched_to_Dmels.fa_filtered.fasta

#Check if your outgroup has been lost:
for f in ./*_filtered.fasta; do
	grep -o '>sim' $f | wc -l
done > mel_lost.txt

for f in ./*_filtered.fasta; do
	grep -o '>yak' $f | wc -l
done > yak_lost.txt

ls > Dmel_matched_to_Dsims_Nov_2017_gene_list.txt

#or simply
for f in ./*_filtered.fasta; do
	echo $f
	grep -o '>sim' $f | wc -l
done

#Check how many lines are left for each sample and find the lower bound
for f in ./*_filtered.fasta; do
	echo $f
	grep -o '>' $f | wc -l
done

#149 (actually 149 + 1 outgroup + 1 outgroup)



#2. Randomly select 149 sequences.

for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/*_filtered.fasta; do
	python randomly_select_Dmel_149_Nov_2017_MK.py $f
done

mkdir /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK
mv /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/*_selected.fa /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK

#Output: Dsim_FBgn0261241_14_selected.fa

#Check sim and yak were still included after selection
for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_selected.fa; do
	echo $f
	grep -o '>sim' $f | wc -l
done

for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_selected.fa; do
	echo $f
	grep -o '>yak' $f | wc -l
done

#Check if 151 have been selected
for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_selected.fa; do
	echo $f
	grep -o '>' $f | wc -l
done


#Make sure that screen is on


#3. Align the sequences.

export PATH=$PATH:/workdir/ji72/prank-msa/src
for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_149_selected.fa; do
	prank -d=$f `echo -o=$f`_aligned -F -iterate=4
echo "Done aligning $f"
done

#Output: Dsim_FBgnXXXXXXX.fa_aligned.best.fas

#When done, send me an email indicating that the job is done:
echo "Your Dsim sequences are aligned" | mailx -s "Your Dsim sequences are aligned" joohyun.im@gmail.com


#4. Remove sequences of fly lines that are complelety unusable (Ns at all positions)

for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*.fas; do
	python only_save_well_aligned_sequences_Nov_2017_MK.py $f >`echo $f`_confirmed.fa
done

#Output: 
#In theory, only the ones that pass the filter will be made into new documents.

#Remove anything without contents
ls -lSr

#Check the number of samples (should be 151)
for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_confirmed.fa; do
	echo $f
	grep -o '>' $f | wc -l
done


#5. Put the outgroups on top
echo "Let's re-order the sequences in order to run the DH test"

for f in /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_confirmed.fa; do
	python put_outgroup_on_top_MK_Nov_2017.py $f > `echo $f`_ordered.fa
done


mkdir /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK_run
cp /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK/*_ordered.fa /workdir/ji72/Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK_run


#6. Run MK test
perl MK_Yasir_Jeff_Vedanayagam_David_Begun.pl -outfile Dmels_Nov_2017_MK_results.txt -pol pol -outgroup sim -outgroup yak -ingroup mel -dir ./Dsim_matched_to_Dmels_Nov_2017/Dmel_149_selected_MK_run/
#perl MK_Yasir_Jeff_Vedanayagam_David_Begun.pl -outfile Dsims_Nov_2017_MK_results.txt -pol pol -outgroup mel -outgroup yak -ingroup sim -dir ./Dsim_14_selected_MK_prepped/


#When done, send me an email indicating that the job is done:
echo "Your MK test is done for Dsim" | mailx -s "Your MK run" joohyun.im@gmail.com