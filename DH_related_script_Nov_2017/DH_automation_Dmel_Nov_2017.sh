#DH_automation_Dmel_Nov_2017.sh
#Run the DH test
#Nov 11th, 2017

#Starting materials:
#Dmel_FBgnXXXXX.fa which contains one Dsim sequence and 180+ Dmel sequences
#Dsim_FBgnXXXXX.fa which contains one Dmel sequence and 20 Dsim sequences

#Make sure the following are present in the /workdir folder before starting this run:
#randomly_select_Dmel_149_Nov_2017.py
#Move prank to the /workdir: cp -r /home/ji72/prank-msa/ /workdir/ji72/ | cd /workdir/ji72/ | export PATH=$PATH:/workdir/ji72/prank-msa/src
#only_save_well_aligned_sequences_Nov_2017.py -- change the setting to the right species
#put_outgroup_on_top_DH_Nov_2017.py -- change the setting to the right species
#Move DH to the /workdir


#1. Remove lines from each gene sequence that have missing data

for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/*.fa; do
     python remove_lines_with_missing_data_Nov_2017.py $f >`echo $f`_filtered.fasta
done

#Output: Dmel_FBgn0003654.fa_filtered.fasta
#Remove the following: FBgn0003654, FBgn0030956, FBgn0031939, FBgn0243514

#2. Randomly select 149 sequences.

for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/*.fasta; do
	python randomly_select_Dmel_149_Nov_2017.py $f
done

mkdir /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected
mv /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/*_149_selected.fa /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected

#Output: Dmel_FBgn0261514.fa_149_selected.fa

#2. Align the sequences.

for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/*_149_selected.fa; do
	prank -d=$f `echo -o=$f`_aligned -F -iterate=4
echo "Done aligning $f"
done

#Output: Dmel_FBgnXXXXXXX.fa_aligned.best.fas

#When done, send me an email indicating that the job is done:
echo "Your Dmel sequences are aligned" | mailx -s "Your Dmel sequences are aligned" joohyun.im@gmail.com


#3. Remove sequences of fly lines that are complelety unusable (Ns at all positions)

for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/*.fas; do
	python only_save_well_aligned_sequences_Nov_2017_2.py $f >`echo $f`_confirmed.fa
done #> only_save_well_aligned_sequences_Nov_2017_output.txt

#Output: Dmel_FBgn0000146.fa_149_selected.fa_aligned.best.fas_confirmed.fa

#4. Order the sequences correctly in order to run the DH test (outgroup first)
echo "Let's re-order the sequences in order to run the DH test"


for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/*_confirmed.fa; do
	python put_outgroup_on_top_DH_Nov_2017.py $f
done

#Output: Dmel_FBgn0000146.fa_149_selected_aligned_confirmed_ordered.fasta

#When done, send me an email indicating that the job is done:
echo "Your Dmel sequences are ready for DH run" | mailx -s "Your DH run is about to start for Dmel" joohyun.im@gmail.com


#Check if outgroups are present
for f in ./*.fasta; do
	echo $f
	grep -o '>sim' $f | wc -l
done

#Check if all 150 lines are present
for f in ./*.fasta; do
	echo $f
	grep -o '>' $f | wc -l
done

#5. Run DH test
echo "Let's run the DH test"

#Take the aligned, confirmed and ordered data as input

#Method 1: using a for loop
#for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/*.fasta; do
#	seq_num=$(grep -o '>' $f | wc -l)
#	java -cp dh.jar dh.ReadFastaDHEW $f $seq_num 1000 > `echo $f`_Dmel_DHEW_results.txt
#done
#echo "Your DHEW is done. Now extracting information..." | mailx -s "Your DHEW run - Part I" joohyun.im@gmail.com

#Method 2: parallelize the commands
perl /programs/bin/perlscripts/perl_fork_univ.pl DH_Dmel_parallel.txt 510


#6. Extract the relevant information
echo "Let's extract the summary statistics"

for f in /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/*.txt; do
	File_name=$(echo $f)
	S_seg=$(grep "Number of segregating sites: *" $f)
	Watt_theta=$(grep "Watterson's theta_W = *" $f)
	TajD=$(grep "Tajima's D = *" $f)
	TajD_pval=$(grep "p{D_neutral <= D_obs} = *" $f)
	nFWH=$(grep "Fay and Wu's H = *" $f)
	nFWH_pval=$(grep "p{H_neutral <= H_obs} = *" $f)
	EW=$(grep "The value of the Ewens-Watterson (EW) test statistic = *" $f)
	EW_pval=$(grep "p{EW_neutral >= EW_obs} = *" $f)
	DH_pval=$(grep "p-value of the DH test = *" $f)
	HEW_pval=$(grep "p-value of the HEW test = *" $f)
	DHEW_pval=$(grep "p-value of the DHEW test = *" $f)

	destdir=`echo $f`_Dmel_DHEW_extracted.txt

		#Save the variables into a file
	echo -e "${File_name}\t${S_seg}\t${Watt_theta}\t${TajD}\t${TajD_pval}\t${nFWH}\t${nFWH_pval}\t${EW}\t${EW_pval}\t${
DH_pval}\t${HEW_p
val}\t${DHEW_pval}" > "$destdir"
	echo "We've extracted the summary statistics for this gene"
done

cat /workdir/ji72/DH/Dmel_prepped_for_DH_Nov_2017/Dmel_149_selected/*_Dmel_DHEW_extracted.txt > /workdir/ji72/DH/Dmel_DHEW_extracted_compiled_Nov_2017.txt


#When done, send me an email indicating that the job is done:
echo "Your DHEW is done" | mailx -s "Your DHEW run" joohyun.im@gmail.com