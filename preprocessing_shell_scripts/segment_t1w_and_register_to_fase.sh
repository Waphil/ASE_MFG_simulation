#This script performs brain tissue segmentation on the t1-w images and registers the masks to the distortion corrected fase images.


    codedir="SPECIFY" # Where python code is stored
	indir="SPECIFY"
	outdir=$indir/derivatives
	FASE_image='FASE'
    output_folder='segmentations'
    
    t1w_name_input="Ax FSPGR  BRAVO"
    t1w_name="T1w_image"

    phase_encoding_direction="_PA" #"_PA" # Refers only to the FASE image that is the target for registration
    postfix="_P2" # Refers only to the FASE image that is the target for registration

    fase_image_name="FASE${phase_encoding_direction}${postfix}_0_aftertopup"
    fase_brain_mask_name="FASE${phase_encoding_direction}${postfix}_0_aftertopup_brain_mask_ero"

    is_skip_bet=true
    is_skip_fast=true

    patient_list=('volunteer02') #('volunteer01') ('volunteer01_session2')



for patient in "${patient_list[@]}" ;do

	echo "... Dealing with patient  : ${patient}"

    segdir=$outdir/${patient}/${output_folder}
	FASEdir=$outdir/${patient}/FASE
	mkdir -p ${segdir}
	cd ${segdir}

    # Copy t1w images
    cp ${indir}/${patient}/*-"${t1w_name_input}.json" ./${t1w_name}.json
    cp $indir/${patient}/*-"${t1w_name_input}.nii.gz" ./${t1w_name}.nii.gz
    # Copy the FASE image and mask
    cp $FASEdir/"${fase_image_name}.nii.gz" ./${fase_image_name}.nii.gz
    cp $FASEdir/"${fase_brain_mask_name}.nii.gz" ./${fase_brain_mask_name}.nii.gz

    if [ "$is_skip_bet" = false ]; then
	    echo "... Extracting the brain mask"
        # Get brain mask
        #bet ${t1w_name}.nii.gz ${t1w_name}_brain.nii.gz -f 0.7 -m  # Old version with FSL bet, not so good
        hd-bet -i ${t1w_name}.nii.gz -device cpu -mode fast -tta 0 # New version with HD-BET, good results
    fi

    if [ "$is_skip_bet" = false ]; then
	    echo "... Segmenting the brain tissues"
        # Get brain tissue segmentations using FSL FAST
        #fast ${t1w_name}_brain.nii.gz # This name is correct for FSL bet
        fast ${t1w_name}_bet.nii.gz # This name is correct for HD-BET
    fi

	echo "... Registering the brain tissues to the FASE image"
    # Do linear registration of the brain tissue masks to the FASE image
    flirt -dof 6 -in ${t1w_name}.nii.gz -ref ${fase_image_name}.nii.gz -out ${t1w_name}_registered${phase_encoding_direction}${postfix}.nii.gz -refweight ${fase_brain_mask_name}.nii.gz -inweight ${t1w_name}_bet_mask.nii.gz -omat reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat -cost corratio -searchrx -10 10 -searchry -10 10 -searchrz -10 10 -interp trilinear -noresample

    # Apply transformation to the segmented partial volume maps
    flirt -in ${t1w_name}_bet_pve_0.nii.gz -ref ${fase_image_name}.nii.gz -out ${t1w_name}_bet_pve_0_reg${phase_encoding_direction}${postfix}.nii.gz -init reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat -applyxfm -interp trilinear
    flirt -in ${t1w_name}_bet_pve_1.nii.gz -ref ${fase_image_name}.nii.gz -out ${t1w_name}_bet_pve_1_reg${phase_encoding_direction}${postfix}.nii.gz -init reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat -applyxfm -interp trilinear
    flirt -in ${t1w_name}_bet_pve_2.nii.gz -ref ${fase_image_name}.nii.gz -out ${t1w_name}_bet_pve_2_reg${phase_encoding_direction}${postfix}.nii.gz -init reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat -applyxfm -interp trilinear
    flirt -in ${t1w_name}_bet_mask.nii.gz -ref ${fase_image_name}.nii.gz -out ${t1w_name}_bet_mask_reg${phase_encoding_direction}${postfix}.nii.gz -init reg_mat_t1w_to${phase_encoding_direction}${postfix}.mat -applyxfm -interp trilinear

done