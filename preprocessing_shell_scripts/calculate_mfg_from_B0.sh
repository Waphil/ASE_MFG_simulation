#This script performs a correction of 

  codedir="SPECIFY"
	indir="SPECIFY"
	outdir=$indir/derivatives

    output_folder='mfg_from_b0'
    segmentations_folder='segmentations'
    fase_folder='FASE'

    is_do_b0_unwrap=true

    b0_name_input="B0 Map whole brain for P2"
    b0_name="b0map"
    
    t1w_name="T1w_image"

    phase_encoding_direction="_RL" #"_PA" # Refers only to the FASE image that is the target for registration
    postfix="_P4" # "_P2"# Refers only to the FASE image that is the target for registration

    fase_image_name="FASE${phase_encoding_direction}${postfix}_0_aftertopup"
    fase_brain_mask_name="FASE${phase_encoding_direction}${postfix}_0_aftertopup_brain_mask_ero"


    patient_list=('volunteer02') #('volunteer01') #('volunteer01_session2')



for patient in "${patient_list[@]}" ;do

	echo "... Dealing with patient  : ${patient}"

    fasedir=$outdir/${patient}/${fase_folder}
    segdir=$outdir/${patient}/${segmentations_folder}
	mfgdir=$outdir/${patient}/${output_folder}
	mkdir -p ${mfgdir}
	cd ${mfgdir}

    # Copy B0 map
    cp ${indir}/${patient}/*-"${b0_name_input}.json" ./${b0_name}.json
    cp $indir/${patient}/*-"${b0_name_input}.nii.gz" ./${b0_name}.nii.gz
    # Copy t1w images and mask
    cp $segdir/"${t1w_name}.json" ./${t1w_name}.json
    cp $segdir/"${t1w_name}.nii.gz" ./${t1w_name}.nii.gz
    cp $segdir/${t1w_name}_bet_mask.nii.gz ./${t1w_name}_bet_mask.nii.gz
    # Copy the FASE image and mask
    cp $fasedir/"${fase_image_name}.nii.gz" ./${fase_image_name}.nii.gz
    cp $fasedir/"${fase_brain_mask_name}.nii.gz" ./${fase_brain_mask_name}.nii.gz

    # Copy the flirt transformation matrix
    #cp $segdir/reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat ./reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat
    # Define affine matrix corresponding to identity transform
	echo -e "1. 0. 0. 0. \n0. 1. 0. 0.\n0. 0. 1. 0.\n0. 0. 0. 1." > identity.mat

    # Resample T1-w image into the b0 map so that we have a "magnitude image" for PRELUDE
    # IMPORTANT: Somehow this code does not work as intended, but it is easy to resample the T1-w image to the B0 map in FSLEYES
    #flirt -in ${t1w_name}.nii.gz -ref ${b0_name}.nii.gz -out ${t1w_name}_resampledtob0.nii.gz -init identity.mat -applyxfm -interp trilinear
    #flirt -in ${t1w_name}.nii.gz -ref ${b0_name}.nii.gz -out ${t1w_name}_resampledtob0.nii.gz -interp trilinear

    if [ "$is_do_b0_unwrap" = true ]; then
        # Rescale b0 map from hz to radians
        fslmaths ${b0_name}.nii.gz -div 198 -mul 6.283 ${b0_name}_rescaled.nii.gz

        # Apply PRELUDE to unwrap b0 image
        prelude -a ${t1w_name}_resampledtob0.nii.gz -p ${b0_name}_rescaled.nii.gz -o ${b0_name}_unwrapped_rescaled.nii.gz
        #prelude -p ${b0_name}.nii.gz -o ${b0_name}_unwrapped.nii.gz

        # Rescale unwrapped B0 back to Hz
        fslmaths ${b0_name}_unwrapped_rescaled.nii.gz -div 6.283 -mul 198 ${b0_name}_unwrapped.nii.gz
    else
        cp ${b0_name}.nii.gz ${b0_name}_unwrapped.nii.gz
    fi

    # Transform unwrapped B0 into FASE coordinates
    flirt -dof 6 -in ${t1w_name}_resampledtob0.nii.gz -ref ${fase_image_name}.nii.gz -out ${t1w_name}_resampledtob0_registered${phase_encoding_direction}${postfix}.nii.gz -refweight ${fase_brain_mask_name}.nii.gz -omat reg_mat_t1w_resampled_to${phase_encoding_direction}${postfix}.mat -cost corratio -searchrx -10 10 -searchry -10 10 -searchrz -10 10 -interp trilinear -noresample

    flirt -in ${b0_name}_unwrapped.nii.gz -ref ${fase_image_name}.nii.gz -out ${b0_name}_unwrapped_reg${phase_encoding_direction}${postfix}.nii.gz -init reg_mat_t1w_resampled_to${phase_encoding_direction}${postfix}.mat -applyxfm -interp trilinear

    # Use python script to calculate magnetic field gradients (in Hz/mm)
    python3 $codedir/python_command_line_scripts/x.calculate_mfg.py --inpath ${b0_name}_unwrapped_reg${phase_encoding_direction}${postfix}.nii.gz    
done