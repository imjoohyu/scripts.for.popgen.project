#!/bin/bash

#DH_automation_Jul_2017.sh
#Run the DH test
#July 16th, 2017

#Starting materials:
#Dmel_FBgnXXXXX.fa which contains one Dsim sequence and 180+ Dmel sequences
#Dsim_FBgnXXXXX.fa which contains one Dmel sequence and 20 Dsim sequences

#Make sure the following are present in the /workdir folder before starting this run:
#randomly_select_Dmel_49_Jul_2017.py
#Move prank to the /workdir: cp -r /home/ji72/prank-msa/ /workdir/ji72/ | cd /workdir/ji72/ | export PATH=$PATH:/workdir/ji72/prank-msa/src
#only_save_well_aligned_sequences_Jul_2017.py -- change the setting to the right species
#put_outgroup_on_top_DH_Jul_2017.py -- change the setting to the right species
#Move DH to the /workdir



#1. If working with Dmel sequences, randomly select 49 sequences.

if [ $1 == 'Dmel']
then
	for f in /workdir/ji72/Dmel_prepped_for_DH/*.fa; do
		python randomly_select_Dmel_49_Jul_2017.py $f
	done

mkdir /workdir/ji72/Dmel_prepped_for_DH/Dmel_49_selected
mv /workdir/ji72/Dmel_prepped_for_DH/*_49_selected.fa /workdir/ji72/DH/Dmel_prepped_for_DH/Dmel_49_selected

fi


#2. Align the sequences.

if [ $1 == 'Dmel']
then
	for f in /workdir/ji72/Dmel_prepped_for_DH/Dmel_49_selected/*_49_selected.fa; do
		prank -d=$f `echo -o=$f`_aligned -F -iterate=4
	echo "Done aligning $f"
	done
	mv /workdir/ji72/Dmel_prepped_for_DH/Dmel_49_selected/*.fas /workdir/ji72/Dmel_prepped_for_DH/
fi

if [ $1 == 'Dsim']
then
	for f in /workdir/ji72/Dsim_prepped_for_DH/*.fa; do
		prank -d=$f `echo -o=$f`_aligned -F -iterate=4
	echo "Done aligning $f"
	done
fi

#Output: Dmel_FBgnXXXXXXX.fa_aligned.best.fas


#3. Remove sequences of fly lines that are complelety unusable (Ns at all positions)

if [ $1 == 'Dmel']
then
	for f in /workdir/ji72/Dmel_prepped_for_DH/*.fas; do
		python only_save_well_aligned_sequences_Jul_2017.py $f >`echo $f`_confirmed.fasta
	done
fi

if [ $1 == 'Dsim']
then
	for f in /workdir/ji72/Dsim_prepped_for_DH/*.fas; do
		python only_save_well_aligned_sequences_Jul_2017.py $f >`echo $f`_confirmed.fasta
	done
fi


#4. Order the sequences correctly in order to run the DH test (outgroup first)
echo "Let's re-order the sequences in order to run the DH test"

if [ $1 == 'Dmel']
then
	for f in /workdir/ji72/Dmel_prepped_for_DH/*.fasta; do
		python put_outgroup_on_top_DH_Jul_2017.py $f 

	done
fi

if [ $1 == 'Dsim']
then
	for f in /workdir/ji72/Dsim_prepped_for_DH/*.fasta; do
		python put_outgroup_on_top_DH_Jul_2017.py $f 
	done
fi

#When done, send me an email indicating that the job is done:
echo "Your sequences are ready for DH run" | mailx -s "Your DH run is about to start" joohyun.im@gmail.com


#5. Move the files to run DH test
if [ $1 == 'Dmel']
then
	mkdir /workdir/ji72/Dmel_prepped_for_DH_done
	mv /workdir/ji72/Dmel_prepped_for_DH/Dmel_*_aligned_confirmed_ordered.fasta /workdir/ji72/Dmel_prepped_for_DH_done
fi

if [ $1 == 'Dsim']
then
	mkdir /workdir/ji72/Dsim_prepped_for_DH_done
	mv /workdir/ji72/Dsim_prepped_for_DH/Dsim_*_aligned_confirmed_ordered.fasta /workdir/ji72/Dsim_prepped_for_DH_done
fi


#6. Run DH test
echo "Let's run the DH test"

#Take the aligned, confirmed and ordered data

if [ $1 == 'Dmel']
then
	for f in /workdir/ji72/Dmel_prepped_for_DH_done/*_ordered.fasta; do
		seq_num=$(grep -o '>' $f | wc -l)
		java -cp dh.jar dh.ReadFastaDHEW $f $seq_num > `echo $f`_Dmel_DHEW_results.txt
	done
fi

if [ $1 == 'Dsim']
then
	for f in /workdir/ji72/Dsim_prepped_for_DH_done/*_ordered.fasta; do
		seq_num=$(grep -o '>' $f | wc -l)
		java -cp dh.jar dh.ReadFastaDHEW $f $seq_num > `echo $f`_Dsim_DHEW_results.txt
	done
fi


#7. Extract the relevant information
echo "Let's extract the summary statistics"

if [ $1 == 'Dmel']
then
	for f in /workdir/ji72/Dmel_prepped_for_DH_done/*.txt; do
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
		echo -e "${File_name}\t${S_seg}\t${Watt_theta}\t${TajD}\t${TajD_pval}\t${nFWH}\t${nFWH_pval}\t${EW}\t${EW_pval}\t${DH_pval}\t${HEW_p
val}\t${DHEW_pval}" > "$destdir"
		echo "We've extracted the summary statistics for this gene"
	done

	cat /workdir/ji72/Dmel_prepped_for_DH_done/*_Dmel_DHEW_extracted.txt > /workdir/ji72/Dmel_DHEW_extracted_compiled.txt

fi

if [ $1 == 'Dsim']
then
	for f in /workdir/ji72/Dsim_prepped_for_DH_done/*.txt; do
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

		destdir=`echo $f`_Dsim_DHEW_extracted.txt

		#Save the variables into a file
		echo -e "${File_name}\t${S_seg}\t${Watt_theta}\t${TajD}\t${TajD_pval}\t${nFWH}\t${nFWH_pval}\t${EW}\t${EW_pval}\t${DH_pval}\t${HEW_p
val}\t${DHEW_pval}" > "$destdir"
		echo "We've extracted the summary statistics for this gene"
	done

	cat /workdir/ji72/Dsim_prepped_for_DH_done/*_Dsim_DHEW_extracted.txt > /workdir/ji72/Dsim_DHEW_extracted_compiled.txt

fi


#When done, send me an email indicating that the job is done:
echo "Your DHEW is done" | mailx -s "Your DHEW run" joohyun.im@gmail.com