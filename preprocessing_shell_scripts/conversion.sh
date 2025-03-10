#!/usr/bin/env bash

# Original script written by Fatemeh, adjusted by Patrick for Philipp, then adjusted by Philipp to allow for complex data
#This script converts DICOM data to Nifti :
#Before running this script you need to check these parameter:

    is_create_separate_files_for_real_and_imag=true

    maindir="SPECIFY"
	indir="dicoms"
	#outdirname="derivatives"
	#outdirname="niftis"
	outdir=$maindir #Could also choose something else
	
	#file_postfix=""
	# ... This script also generates a text file of $indir/$patient/DicVsNifti.txt  
	
	##########################################################################################

#mkdir -p $maindir/$outdirname
	

for dcm in $maindir/$indir/*/ ;do
	dcm2niix -b y -z y -x n -m n -f output -s n -v n ${dcm}
	cd ${dcm}
	for jsonname in *.json ; do
		name=${jsonname%.json}
		#echo "The name is:"
		#echo ${name}
		protocol=`jq ".SeriesDescription" ${name}.json`
        	protocol=${protocol:1:size-1}
        	sernum=`jq ".SeriesNumber" ${name}.json`
        	
			if [ "$is_create_separate_files_for_real_and_imag" = true ]; then
				if [[ "$name" == *"real"* ]];then
					file_postfix="_real" # If name contains "real", add that to the output file ending
					to_delete_string="real" # In the later step in the loop, onlz files containing this string will be deleted
				elif [[ "$name" == *"imaginary"* ]];then
					file_postfix="_imag" # If name contains "imaginary", add that to the output file ending
					to_delete_string="imaginary"
				else
					file_postfix=""
					to_delete_string=""
  				fi
			else
				file_postfix=""
				to_delete_string=""
			fi

        	echo -e "$(basename ${dcm}), file ${name} ----> ${sernum}-${protocol}${file_postfix}" >> $maindir/$outdirname/DicVsNifti.txt

		cp ${dcm}/${name}.json "$outdir/${sernum}-${protocol}${file_postfix}.json"
		cp ${dcm}/${name}.nii.gz "$outdir/${sernum}-${protocol}${file_postfix}.nii.gz"
        	
			#echo "Deleting the following (new version)"
			#find . -type f \( -name "*.json" -a -name "*${to_delete_string}*" \)
			find . -type f \( -name "*.json" -a -name "*${to_delete_string}*" \) -delete
			find . -type f \( -name "*.nii.gz" -a -name "*${to_delete_string}*" \) -delete
			#echo "Deleting the following (old version)"
        	#find . -name "*.json" -type f
        	#find . -name "*.json" -type f -delete
        	#find . -name "*.nii.gz" -type f -delete
	    	
			index=$((index+1))
        done
	
done
	
    	
