#! /bin/bash

#MK_automation_Dsim_Nov_2017.sh
#Prepare data to run MK test
#Nov 11th, 2017

#Starting materials:
#FBgn0261514_Dsim_Dyak_matched_to_Dmels.fa which contains one Dmel sequence, one Dyak sequence and 20 Dsim sequences

#Make sure the following are present in the /workdir folder before starting this run:
#remove_lines_with_missing_data_Nov_2017.py
#randomly_select_Dmel_14_Nov_2017.py
#Move prank to the /workdir: cp -r /home/ji72/prank-msa/ /workdir/ji72/ | cd /workdir/ji72/ | export PATH=$PATH:/workdir/ji72/prank-msa/src
#only_save_well_aligned_sequences_Nov_2017.py -- change the setting to the right species


#1. Remove lines from each gene sequence that have missing data

#Input: FBgn0033452_Dmel_Dyak_matched_to_Dsims.fa

for f in /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/*.fa; do
     python remove_lines_with_missing_data_Nov_2017.py $f >`echo $f`_filtered.fasta
done

#Output: FBgn0033452_Dmel_Dyak_matched_to_Dsims.fa_filtered.fasta

#Check if your outgroup has been lost:
for f in ./*_filtered.fasta; do
	grep -o '>mel' $f | wc -l
done > mel_lost.txt

for f in ./*_filtered.fasta; do
	grep -o '>yak' $f | wc -l
done > yak_lost.txt

ls > Dmel_matched_to_Dsims_Nov_2017_gene_list.txt

#or simply
for f in ./*_filtered.fasta; do
	echo $f
	grep -o '>mel' $f | wc -l
done

#Check how many lines are left for each sample and find the lower bound
for f in ./*_filtered.fasta; do
	echo $f
	grep -o '>' $f | wc -l
done

#14 looks like a good number (actually 14 + 1 outgroup + 1 outgroup)

#2. Randomly select 14 sequences.

for f in /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/*_filtered.fasta; do
	python randomly_select_Dsim_14_Nov_2017_MK.py $f
done

mkdir /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK
mv /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/*_selected.fa /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK

#Output: Dsim_FBgn0261241_14_selected.fa

#Check mel and yak were still included after selection
for f in /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK/*_selected.fa; do
	echo $f
	grep -o '>mel' $f | wc -l
done


#3. Align the sequences.

export PATH=$PATH:/workdir/ji72/prank-msa/src
for f in /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK/*_14_selected.fa; do
	prank -d=$f `echo -o=$f`_aligned -F -iterate=4
echo "Done aligning $f"
done

#Output: Dsim_FBgnXXXXXXX.fa_aligned.best.fas

#When done, send me an email indicating that the job is done:
echo "Your Dsim sequences are aligned" | mailx -s "Your Dsim sequences are aligned" joohyun.im@gmail.com


#4. Remove sequences of fly lines that are complelety unusable (Ns at all positions)

for f in /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK/*.fas; do
	python only_save_well_aligned_sequences_Nov_2017_MK.py $f >`echo $f`_confirmed.fa
done

#Output: 
#In theory, only the ones that pass the filter will be made into new documents.

for f in ./*_confirmed.fa; do
	echo $f
	grep -o '>' $f | wc -l
done


#5. Put the outgroups on top
echo "Let's re-order the sequences in order to run the DH test"

for f in /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK/*_confirmed.fa; do
	python put_outgroup_on_top_MK_Nov_2017.py $f > `echo $f`_ordered.fa
done

mkdir /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK_run
cp /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK/*_ordered.fa /workdir/ji72/Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK_run


#6. Run MK test
perl MK_Yasir_Jeff_Vedanayagam_David_Begun.pl -outfile Dsims_Nov_2017_MK_results.txt -pol pol -outgroup mel -outgroup yak -ingroup sim -dir ./Dmel_matched_to_Dsims_Nov_2017/Dsim_14_selected_MK_run/
#perl MK_Yasir_Jeff_Vedanayagam_David_Begun.pl -outfile Dsims_Nov_2017_MK_results.txt -pol pol -outgroup mel -outgroup yak -ingroup sim -dir ./Dsim_14_selected_MK_prepped/


#When done, send me an email indicating that the job is done:
echo "Your MK test is done for Dsim" | mailx -s "Your MK run" joohyun.im@gmail.com