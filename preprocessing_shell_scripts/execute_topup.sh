
#This script performs topup distortion correction on FASE data and brings it into a form usable by quantiphyse:

    is_data_complex=true

    codedir="SPECIFY"
	indir="SPECIFY"
	outdir=$indir/derivatives
	#mask_registeration='yesn'
	#FASE_analysis='yesn'
	#FASE_IS_analysis='yes'
	FASE_image='FASE'

    select_only_second_volume_patient_list=()
    select_only_second_volume_tau_list=()
    select_only_third_volume_patient_list=()
    select_only_third_volume_tau_list=()

	tau=("00" "03" "06" "09" "12" "15" "18" "21" "24")
	shifts=("0.0" "0.006" "0.012" "0.018" "0.024" "0.030" "0.036" "0.042" "0.048")

    principal_ped="_LR" # "_AP", "_PA"
    secondary_ped="_RL" # "_PA", _AP
    postfix="_P2" # "_P2", "_P4"

    patient_list=('volunteer02') #('volunteer01') #s('volunteer01_session2')


if [ "$is_data_complex" = true ]; then
    complex_postfix_list=("_real" "_imag")
else
    complex_postfix_list=("")
fi

for patient in "${patient_list[@]}" ;do

	echo "... Dealing with patient  : ${patient}"

	FASEdir=$outdir/${patient}/FASE
	mkdir -p ${FASEdir}
	cd ${FASEdir}

	for img in "${tau[@]}" ;do

		echo "... The tau value is : ${img}"

        if [ "$img" = "00" ]; then 
            # Extract acquisition parameters
            echo "... Writing acquisition parameters into a text file"
            echo ${indir}/${patient}/*-${FASE_image}_${img}${principal_ped}${postfix}${complex_postfix_list[0]}.json
            echo acqparams${principal_ped}${postfix}.txt
            python3 $codedir/python_command_line_scripts/x.create_acqparam_file.py --paramsin1 ${indir}/${patient}/*-${FASE_image}_${img}${principal_ped}${postfix}${complex_postfix_list[0]}.json --paramsout acqparams${principal_ped}${postfix}.txt

            # For fase 00 image, process reverse phase encoding direction as well (needed for topup)
            peds_to_be_processed=("${principal_ped}" "${secondary_ped}")
        else
            # Otherwise, process only principal phase encoding direction
            peds_to_be_processed=("${principal_ped}")
        fi

        for ped in "${peds_to_be_processed[@]}" ;do
            
            echo "... Dealing with phase encoding direction: ${ped}"

            # If complex values are provided, do calculation of mean etc for both real and imaginary. Otherwise just for the given image.
            for complex_postfix in "${complex_postfix_list[@]}" ;do

                echo "... Dealing with complex postfix: ${complex_postfix}"



                cp ${indir}/${patient}/*-${FASE_image}_${img}${ped}${postfix}${complex_postfix}.json ./FASE.json
                cp $indir/${patient}/*-${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz ./${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz

                echo "... remove dummy images:"

                # Determine which volumes should be used based on user input
                # I apologize for the ugly code, I am not good at shell scripts
                is_use_only_second_volume=false
                is_use_only_third_volume=false
                for i in "${!select_only_second_volume_patient_list[@]}"; do
                    if [ "$patient" = "${select_only_second_volume_patient_list[i]}" ] && [ "$img" = "${select_only_second_volume_tau_list[i]}" ]; then
                        is_use_only_second_volume=true
                    fi
                done
                for i in "${!select_only_third_volume_patient_list[@]}"; do
                    if [ "$patient" = "${select_only_third_volume_patient_list[i]}" ] && [ "$img" = "${select_only_third_volume_tau_list[i]}" ]; then
                        is_use_only_third_volume=true
                    fi
                done

                if [ "$is_use_only_second_volume" = true ]; then
                    echo "... only second image volume used for this tau"
                    fslroi ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz 1 1
                elif [ "$is_use_only_third_volume" = true ]; then
                    echo "... only third image volume used for this tau"
                    fslroi ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz 2 1
                else
                    #fslroi ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz 1 2
                    fslroi ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz 1 3 # 3 for volunteer
                fi

                echo "... Tmean of volumes"
                fslmaths ${FASE_image}_${img}${ped}${postfix}${complex_postfix}.nii.gz -Tmean ${FASE_image}_${img}${ped}${postfix}${complex_postfix}_mean.nii.gz
            done

            # If the images are real and imaginary, conver them to magnitude and phase here. 
            # The name will then be such that the rest of the script works.
            if [ "$is_data_complex" = true ]; then
                # Calculate magnitude images
                fslmaths ${FASE_image}_${img}${ped}${postfix}_real_mean.nii.gz -sqr ${FASE_image}_${img}${ped}${postfix}_real_mean_squared.nii.gz
                fslmaths ${FASE_image}_${img}${ped}${postfix}_imag_mean.nii.gz -sqr ${FASE_image}_${img}${ped}${postfix}_imag_mean_squared.nii.gz
                fslmaths ${FASE_image}_${img}${ped}${postfix}_real_mean_squared.nii.gz -add ${FASE_image}_${img}${ped}${postfix}_imag_mean_squared.nii.gz ${FASE_image}_${img}${ped}${postfix}_sum_mean_squared.nii.gz
                fslmaths ${FASE_image}_${img}${ped}${postfix}_sum_mean_squared.nii.gz -sqrt ${FASE_image}_${img}${ped}${postfix}_mean.nii.gz
                
                # Calculate phase images
                fslmaths ${FASE_image}_${img}${ped}${postfix}_imag_mean.nii.gz -div ${FASE_image}_${img}${ped}${postfix}_real_mean.nii.gz ${FASE_image}_${img}${ped}${postfix}_imag_real_ratio_mean.nii.gz
                fslmaths ${FASE_image}_${img}${ped}${postfix}_imag_real_ratio_mean.nii.gz -atan ${FASE_image}_${img}${ped}${postfix}_phase_mean.nii.gz

                # Delete helper images
                rm ${FASE_image}_${img}${ped}${postfix}_real_mean_squared.nii.gz
                rm ${FASE_image}_${img}${ped}${postfix}_imag_mean_squared.nii.gz
                rm ${FASE_image}_${img}${ped}${postfix}_sum_mean_squared.nii.gz
                rm ${FASE_image}_${img}${ped}${postfix}_imag_real_ratio_mean.nii.gz
            fi

            # Create brain mask, is used for registration of images.
            # The mask is eroded in z-direction in order to exclude the lower parts of the brain that show high off-resonances and thus tau-dependent signal
            echo "... Extract brain mask from image"
            bet FASE_${img}${ped}${postfix}_mean.nii.gz FASE_${img}${ped}${postfix}_brain.nii.gz -f 0.6 -m
            rm FASE_${img}${ped}${postfix}_brain.nii.gz
            fslmaths FASE_${img}${ped}${postfix}_brain_mask.nii.gz -kernel boxv3 1 1 5 -ero FASE_${img}${ped}${postfix}_brain_mask_ero.nii.gz

            #echo "... Intensity correction with biasfield"
            #python3 $indir/code/oef/x.biasfield.py --img FASE_${img}_mean.nii.gz --imgout FASE_${img}_mean_bf.nii.gz
            # in case you do not want to use bias field correction, run this instead and comment the line above: 
            cp FASE_${img}${ped}${postfix}_mean.nii.gz FASE_${img}${ped}${postfix}_mean_bf.nii.gz

            if [ "$ped" = "${secondary_ped}" ]; then 
                echo "... For reverse phase encoding direction image, do not do registration:"
                cp FASE_${img}${ped}${postfix}_mean_bf.nii.gz FASE_${img}${ped}${postfix}_mean_bf_reg.nii.gz

                cp FASE_${img}${ped}${postfix}_real_mean.nii.gz FASE_${img}${ped}${postfix}_real_mean_reg.nii.gz
                cp FASE_${img}${ped}${postfix}_imag_mean.nii.gz FASE_${img}${ped}${postfix}_imag_mean_reg.nii.gz
            elif [ "$img" = "00" ]; then 
                echo "... For tau=00 image, do not do registration:"
                cp FASE_${img}${ped}${postfix}_mean_bf.nii.gz FASE_${img}${ped}${postfix}_mean_bf_reg.nii.gz

                cp FASE_${img}${ped}${postfix}_real_mean.nii.gz FASE_${img}${ped}${postfix}_real_mean_reg.nii.gz
                cp FASE_${img}${ped}${postfix}_imag_mean.nii.gz FASE_${img}${ped}${postfix}_imag_mean_reg.nii.gz
            else
                echo "... For other tau images with normal phase encoding direction, do registration:"
                flirt -dof 6 -in FASE_${img}${ped}${postfix}_mean_bf.nii.gz -ref FASE_00${ped}${postfix}_mean_bf.nii.gz -out FASE_${img}${ped}${postfix}_mean_bf_reg.nii.gz -refweight FASE_00${ped}${postfix}_brain_mask_ero.nii.gz -inweight FASE_${img}${ped}${postfix}_brain_mask_ero.nii.gz -omat reg_mat_FASE_${img}${ped}${postfix}.mat -cost corratio -searchrx -10 10 -searchry -10 10 -searchrz -10 10 -interp trilinear -noresample

                # Apply transformation to the complex images
                flirt -in FASE_${img}${ped}${postfix}_real_mean.nii.gz -ref FASE_00${ped}${postfix}_mean_bf.nii.gz -out FASE_${img}${ped}${postfix}_real_mean_reg.nii.gz -init reg_mat_FASE_${img}${ped}${postfix}.mat -applyxfm -noresample
                flirt -in FASE_${img}${ped}${postfix}_imag_mean.nii.gz -ref FASE_00${ped}${postfix}_mean_bf.nii.gz -out FASE_${img}${ped}${postfix}_imag_mean_reg.nii.gz -init reg_mat_FASE_${img}${ped}${postfix}.mat -applyxfm -noresample

                # If no registration is wanted, use this":
                #echo "... For other tau images with normal phase encoding direction, chose to omit registration:"
                #cp FASE_${img}_mean_bf.nii.gz FASE_${img}_mean_bf_reg.nii.gz
            fi
            
            # The following lines were used by Fatmeh, I decided to go a slightly different route with the registration, as seen above.
            #echo "... Registering each volume to the first volume with elastix"
            #elastix -f FASE_00_mean_bf.nii.gz -m FASE_${img}_mean_bf.nii.gz -p $indir/code/oef/x.elastix_parameters_rigid.txt -out .
            #mv result.0.nii.gz FASE_${img}_mean_bf_reg.nii.gz 
            #rm IterationInfo.0*
            #rm TransformParameters.0.txt
            #rm elastix.log

            #In case elastix does not work you can use flirt instead: 
            #flirt -dof 6 -in FASE_${img}_mean_bf.nii.gz -ref FASE_00_mean_bf.nii.gz -out FASE_${img}_mean_bf_reg.nii.gz 
            
            #In case do not want to use elastix nor flirt, then just use the following line:  
            #cp FASE_${img}_mean_bf.nii.gz FASE_${img}_mean_bf_reg.nii.gz
        done
	done

    echo "... Preparing for topup correction"

    echo "... rename two FASE volumes to standardPED and revPED"

    cp FASE_00${principal_ped}${postfix}_mean_bf.nii.gz FASE_standardPED.nii.gz
    cp FASE_00${secondary_ped}${postfix}_mean_bf.nii.gz FASE_reversePED.nii.gz

    # Remove secondary phase encoding direction mean files so they do not get included in the merging on accident
    rm FASE_00${secondary_ped}${postfix}_mean_bf*
    rm FASE_00${secondary_ped}${postfix}_real_mean*
    rm FASE_00${secondary_ped}${postfix}_imag_mean*

    echo "... merging all tau values"
    fslmerge -t FASE${principal_ped}${postfix}_raw.nii.gz `ls FASE_*${principal_ped}${postfix}_mean_bf.nii.gz`
    fslmerge -t FASE${principal_ped}${postfix}_raw_real.nii.gz `ls FASE_*${principal_ped}${postfix}_real_mean_reg.nii.gz`
    fslmerge -t FASE${principal_ped}${postfix}_raw_imag.nii.gz `ls FASE_*${principal_ped}${postfix}_imag_mean_reg.nii.gz`
    #fslmerge -t FASE${principal_ped}${postfix}_raw_real.nii.gz `ls FASE_*${principal_ped}${postfix}_real_mean.nii.gz`
    #fslmerge -t FASE${principal_ped}${postfix}_raw_imag.nii.gz `ls FASE_*${principal_ped}${postfix}_imag_mean.nii.gz`

    # this is needed for topup correction
    fslmerge -t standardPED_reversePED.nii.gz FASE_standardPED.nii.gz FASE_reversePED.nii.gz

    #rm FASE_*_mean*

    echo "... Looking if the number of slices are even for topup correction"

    info=(`fslinfo FASE_standardPED.nii.gz`)
    if [ $((info[7]%2)) -eq 0 ];then 
            echo "perfect"
    else
        fslroi standardPED_reversePED.nii.gz standardPED_reversePED.nii.gz 0 ${info[3]} 0 ${info[5]} 0 $(( info[7] - 1 ))
        fslroi FASE${principal_ped}${postfix}_raw.nii.gz FASE${principal_ped}${postfix}_raw.nii.gz 0 ${info[3]} 0 ${info[5]} 0 $(( info[7] - 1 ))
        fi

    echo "... Applying topup"
    topup --imain=standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --config=b02b0.cnf --out=topup_standardPED_reversePED --fout=topup_B0field${principal_ped}${postfix}
    # Apply topup with jacobian method to all images of the default phase encoding direction
    applytopup --imain=FASE${principal_ped}${postfix}_raw.nii.gz --topup=topup_standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --inindex=1 --method=jac --interp=spline --out=FASE${principal_ped}${postfix}
    # Apply topup with jacobian method to the real and imaginary images
    applytopup --imain=FASE${principal_ped}${postfix}_raw_real.nii.gz --topup=topup_standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --inindex=1 --method=jac --interp=spline --out=FASE${principal_ped}${postfix}_real
    applytopup --imain=FASE${principal_ped}${postfix}_raw_imag.nii.gz --topup=topup_standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --inindex=1 --method=jac --interp=spline --out=FASE${principal_ped}${postfix}_imag
    # Apply topup with least squares method to the tau=0 images.
    #applytopup --imain=FASE_AP,FASE_PA --topup=topup_AP_PA --datain=acqparams.txt --inindex=1,2 --method=lsr --interp=spline --out=topup_lsq_AP_PA
    # Apply topup with least squares method to all images (will only work if files exist).
    applytopup --imain=FASE${principal_ped}${postfix}_raw,FASE${secondary_ped}${postfix}_raw --topup=topup_standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --inindex=1,2 --method=lsr --interp=spline --out=FASE${principal_ped}${postfix}_topuplsqr
    # Do also for the real and imaginary images (no idea if this is a good plan or not)
    applytopup --imain=FASE${principal_ped}${postfix}_raw_real,FASE${secondary_ped}${postfix}_raw_real --topup=topup_standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --inindex=1,2 --method=lsr --interp=spline --out=FASE${principal_ped}${postfix}_topuplsqr_real
    applytopup --imain=FASE${principal_ped}${postfix}_raw_imag,FASE${secondary_ped}${postfix}_raw_imag --topup=topup_standardPED_reversePED --datain=acqparams${principal_ped}${postfix}.txt --inindex=1,2 --method=lsr --interp=spline --out=FASE${principal_ped}${postfix}_topuplsqr_imag

    # Final steps to get normalized data which can be used in quantiphyse
	fslroi FASE${principal_ped}${postfix}.nii.gz FASE${principal_ped}${postfix}_0_aftertopup.nii.gz 0 1
	bet FASE${principal_ped}${postfix}_0_aftertopup.nii.gz FASE${principal_ped}${postfix}_0_aftertopup_brain.nii.gz -f 0.2 -m
	fslmaths FASE${principal_ped}${postfix}_0_aftertopup_brain_mask.nii.gz -ero FASE${principal_ped}${postfix}_0_aftertopup_brain_mask_ero.nii.gz
	fslmaths FASE${principal_ped}${postfix}.nii.gz -div FASE${principal_ped}${postfix}_0_aftertopup.nii.gz FASE${principal_ped}${postfix}_norm.nii.gz
	fslmaths FASE${principal_ped}${postfix}.nii.gz -div FASE${principal_ped}${postfix}_0_aftertopup.nii.gz -mul FASE${principal_ped}${postfix}_0_aftertopup_brain_mask_ero.nii.gz FASE${principal_ped}${postfix}_norm.nii.gz

    # Also get the normalized image for the leastsquares topup image
	fslroi FASE${principal_ped}${postfix}_topuplsqr.nii.gz FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup.nii.gz 0 1
	bet FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup.nii.gz FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup_brain.nii.gz -f 0.2 -m
	fslmaths FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup_brain_mask.nii.gz -ero FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup_brain_mask_ero.nii.gz
	fslmaths FASE${principal_ped}${postfix}_topuplsqr.nii.gz -div FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup.nii.gz FASE${principal_ped}${postfix}_topuplsqr_norm.nii.gz
	fslmaths FASE${principal_ped}${postfix}_topuplsqr.nii.gz -div FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup.nii.gz -mul FASE${principal_ped}${postfix}_topuplsqr_0_aftertopup_brain_mask_ero.nii.gz FASE${principal_ped}${postfix}_topuplsqr_norm.nii.gz
done