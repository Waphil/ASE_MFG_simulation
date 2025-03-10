# This script plots the slice visualizations and scatter plots of MFGs, logarithmic signal ratios and echo shifts
# for acquired data in nifti format. Used to create figures 7 & 8 in paper.


import os
import numpy as np
import scipy as sp
import nibabel as nib
import sys
sys.path.append('./')
from utilities.constants_and_helpers import convert_mfg_from_Hz_per_mm_to_microtesla_per_m
from visualization.scatter_plot import scatter_plot_2d_with_histograms
from visualization.image_show_plot import plot_slice_from_3d_image



def process_patient(input_settings_dict, plot_settings_dict=None):
    #####################################
    # Read settings from provided dicts #
    #####################################

    input_folder = input_settings_dict["input_folder"]
    patient_name = input_settings_dict["patient_name"]
    ped_1_file_label = input_settings_dict["ped_1_file_label"]
    ped_2_file_label = input_settings_dict["ped_2_file_label"]
    parallel_imaging_file_label = input_settings_dict["parallel_imaging_file_label"]
    ped_1_true = input_settings_dict["ped_1_true"]
    ped_2_true = input_settings_dict["ped_2_true"]
    mfg_input_folder = input_settings_dict["mfg_input_folder"]
    mfg_multiplication_factor_list = input_settings_dict["mfg_multiplication_factor_list"]
    display_slice_index = input_settings_dict["display_slice_index"]
    tau_list = input_settings_dict["tau_list"]
    evaluate_tau_index_list = input_settings_dict["evaluate_tau_index_list"]
    is_use_topup_image = input_settings_dict["is_use_topup_image"]

    if evaluate_tau_index_list is None:
        evaluate_tau_index_list = np.arange(len(tau_list))
    if mfg_multiplication_factor_list is None:
        mfg_multiplication_factor_list = [1., 1., 1.]

    if plot_settings_dict is None:
        is_show_only_brain_on_plotted_slices = False
        is_skip_plots = True
        is_save_plots = False
    else:
        is_skip_plots = False
        log_mag_ratio_range = plot_settings_dict["log_mag_ratio_range"]
        tau_eff_difference_range = plot_settings_dict["tau_eff_difference_range"]
        gy_range = plot_settings_dict["gy_range"]
        b0_range = plot_settings_dict["b0_range"]
        is_show_only_brain_on_plotted_slices = plot_settings_dict["is_show_only_brain_on_plotted_slices"]
        scatter_plot_grid_kwargs = plot_settings_dict["scatter_plot_grid_kwargs"]
        stairs_kwargs = plot_settings_dict["stairs_kwargs"]
        figsize_slice_plot = plot_settings_dict["figsize_slice_plot"]
        figsize_scatter_plot = plot_settings_dict["figsize_scatter_plot"]
        fontsize_scatter_plot_labels = plot_settings_dict["fontsize_scatter_plot_labels"]
        tau_sub_label = plot_settings_dict["tau_sub_label"]
        tau_sub_unit = plot_settings_dict["tau_sub_unit"]
        log_mag_ratio_label = plot_settings_dict["log_mag_ratio_label"]
        log_mag_ratio_unit = plot_settings_dict["log_mag_ratio_unit"]
        ped_mfg_label = plot_settings_dict["ped_mfg_label"]
        ped_mfg_unit = plot_settings_dict["ped_mfg_unit"]

        plot_folder_path = plot_settings_dict.get("plot_folder_path", None)
        if plot_folder_path is not None:
            # If data for a plot location is given, determine the appropriate folders and names here.
            is_save_plots = True
            plot_name_component_log_sig_rat = plot_settings_dict["plot_name_component_log_sig_rat"]
            plot_name_component_echo_shift = plot_settings_dict["plot_name_component_echo_shift"]
            plot_name_component_gy = plot_settings_dict["plot_name_component_gy"]
            plot_name_component_read_mfg = plot_settings_dict["plot_name_component_read_mfg"]
            plot_name_component_phase_mfg = plot_settings_dict["plot_name_component_phase_mfg"]
            plot_name_component_slice_mfg = plot_settings_dict["plot_name_component_slice_mfg"]
            plot_name_component_b0 = plot_settings_dict["plot_name_component_b0"]

            plot_patient_folder_path = os.path.join(plot_folder_path, patient_name)

            plot_name_component_true_ped = ped_1_true.lower()
            plot_name_component_parallel_imaging = "" if parallel_imaging_file_label == "_P2" else parallel_imaging_file_label.lower()
            plot_name_component_slice = f"slice{display_slice_index}"

            if not os.path.exists(plot_patient_folder_path):
                os.mkdir(plot_patient_folder_path)
        else:
            is_save_plots = False
    
    print("-------------------------------------------------------------------")
    print(f"Evaluating patient from folder: {input_folder}")

    #############
    # Load data #
    #############

    if not is_use_topup_image:
        magnitude_image_file_path_1 = os.path.join(input_folder, f"FASE{ped_1_file_label}{parallel_imaging_file_label}_raw.nii.gz")
        magnitude_image_file_path_2 = os.path.join(input_folder, f"FASE{ped_2_file_label}{parallel_imaging_file_label}_raw.nii.gz")
        echo_shift_time_file_path_1 = os.path.join(input_folder, f"FASE{ped_1_file_label}{parallel_imaging_file_label}_raw_echoshift.nii.gz")
        echo_shift_time_file_path_2 = os.path.join(input_folder, f"FASE{ped_2_file_label}{parallel_imaging_file_label}_raw_echoshift.nii.gz")
        
        log_magnitude_ratio_output_file_path = os.path.join(input_folder, f"FASE_{ped_1_file_label}{parallel_imaging_file_label}_raw_phaseencode_log_magnitude_ratio.nii.gz")
    else:
        magnitude_image_file_path_1 = os.path.join(input_folder, f"FASE{ped_1_file_label}{parallel_imaging_file_label}.nii.gz")
        magnitude_image_file_path_2 = os.path.join(input_folder, f"FASE{ped_2_file_label}{parallel_imaging_file_label}.nii.gz")
        echo_shift_time_file_path_1 = os.path.join(input_folder, f"FASE{ped_1_file_label}{parallel_imaging_file_label}_echoshift.nii.gz")
        echo_shift_time_file_path_2 = os.path.join(input_folder, f"FASE{ped_2_file_label}{parallel_imaging_file_label}_echoshift.nii.gz")

        log_magnitude_ratio_output_file_path = os.path.join(input_folder, f"FASE{ped_1_file_label}{parallel_imaging_file_label}_phaseencode_log_magnitude_ratio.nii.gz")

    brain_mask_path = os.path.join(input_folder, f"FASE{ped_1_file_label}{parallel_imaging_file_label}_0_aftertopup_brain_mask_ero.nii.gz")

    b0_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}.nii.gz")
    
    # Load in MFGs from B0 map
    if ped_1_true in ["AP", "PA"]:
        mfg_read_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_0.nii.gz")
        mfg_phase_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_1.nii.gz")
        mfg_slice_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_2.nii.gz")
    else:
        mfg_read_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_1.nii.gz")
        mfg_phase_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_0.nii.gz")
        mfg_slice_path = os.path.join(mfg_input_folder, f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_2.nii.gz")

    magnitude_image_1 = nib.load(magnitude_image_file_path_1)
    magnitude_image_2 = nib.load(magnitude_image_file_path_2)
    echo_shift_image_1 = nib.load(echo_shift_time_file_path_1)
    echo_shift_image_2 = nib.load(echo_shift_time_file_path_2)
    brain_mask_image = nib.load(brain_mask_path)
    b0_image = nib.load(b0_path)
    mfg_read_image = nib.load(mfg_read_path)
    mfg_phase_image = nib.load(mfg_phase_path)
    mfg_slice_image = nib.load(mfg_slice_path)

    ############################
    # Convert & Calculate data #
    ############################

    magnitude_arr_1 = magnitude_image_1.get_fdata()
    magnitude_arr_2 = magnitude_image_2.get_fdata()
    echo_shift_arr_1 = echo_shift_image_1.get_fdata()
    echo_shift_arr_2 = echo_shift_image_2.get_fdata()
    brain_mask_arr = brain_mask_image.get_fdata() > 0
    b0_arr = b0_image.get_fdata()
    mfg_read_arr = convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg_read_image.get_fdata())*mfg_multiplication_factor_list[0]
    mfg_phase_arr = convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg_phase_image.get_fdata())*mfg_multiplication_factor_list[1]
    mfg_slice_arr = convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg_slice_image.get_fdata())*mfg_multiplication_factor_list[2]

    log_magnitude_ratio = np.log(magnitude_arr_1/magnitude_arr_2)
    echo_shift_difference = echo_shift_arr_1 - echo_shift_arr_2

    masked_mfg_read_arr = mfg_read_arr[brain_mask_arr]
    masked_mfg_phase_arr = mfg_phase_arr[brain_mask_arr]
    masked_mfg_slice_arr = mfg_slice_arr[brain_mask_arr]


    plot_brain_mask_arr = brain_mask_arr if is_show_only_brain_on_plotted_slices else None

    # Print median t_eff across brain, because it is interesting.
    median_echo_shift_1 = np.nanmedian(echo_shift_arr_1[brain_mask_arr], axis=0)
    median_echo_shift_2 = np.nanmedian(echo_shift_arr_2[brain_mask_arr], axis=0)
    print(f"For tau values: {tau_list}")
    print(f"Median echo shift values in {ped_1_true} direction [ms]: {list(median_echo_shift_1)}")
    print(f"Median echo shift values in {ped_2_true} direction [ms]: {list(median_echo_shift_2)}")

    ###########################################################
    # Plot single slice to ensure data sensibility, if chosen #
    ###########################################################

    if not is_skip_plots:
        for evaluate_tau_index in evaluate_tau_index_list:
            if is_save_plots:
                plot_echo_shift_diff_slice_path = os.path.join(
                    plot_patient_folder_path, 
                    f"{plot_name_component_echo_shift}_tau{tau_list[evaluate_tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}_{plot_name_component_slice}.pdf"
                    )
                plot_logsigrat_slice_path = os.path.join(
                    plot_patient_folder_path, 
                    f"{plot_name_component_log_sig_rat}_tau{tau_list[evaluate_tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}_{plot_name_component_slice}.pdf"
                    )
            else:
                plot_echo_shift_diff_slice_path = None
                plot_logsigrat_slice_path = None

            plot_slice_from_3d_image(echo_shift_difference[:,:,:,evaluate_tau_index], display_slice_index, plot_brain_mask_arr,
                                    title=f"{tau_sub_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[evaluate_tau_index]} \\;ms$",
                                    val_range=tau_eff_difference_range, cmap="RdBu",
                                    is_show_colorbar=True, colorbar_label=tau_sub_unit, figsize=figsize_slice_plot, is_show=False,
                                    save_path=plot_echo_shift_diff_slice_path)
            plot_slice_from_3d_image(log_magnitude_ratio[:,:,:,evaluate_tau_index], display_slice_index, plot_brain_mask_arr,
                                    title=f"{log_mag_ratio_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[evaluate_tau_index]} \\;ms$",
                                    val_range=log_mag_ratio_range, cmap="PuOr",
                                    is_show_colorbar=True, colorbar_label=log_mag_ratio_unit, figsize=figsize_slice_plot, is_show=False,
                                    save_path=plot_logsigrat_slice_path)
        #show_slice(t0_to_tau_slope[:,:,display_slice_index], (-0.3, 0.3), "PRGn", f"slope of delta t0 to tau")

        if is_save_plots:
            plot_b0_slice_path = os.path.join(
                plot_patient_folder_path, 
                f"{plot_name_component_b0}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}_{plot_name_component_slice}.pdf"
                )
            plot_read_mfg_slice_path = os.path.join(
                plot_patient_folder_path, 
                f"{plot_name_component_read_mfg}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}_{plot_name_component_slice}.pdf"
                )
            plot_phase_mfg_slice_path = os.path.join(
                plot_patient_folder_path, 
                f"{plot_name_component_phase_mfg}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}_{plot_name_component_slice}.pdf"
                )
            plot_slice_mfg_slice_path = os.path.join(
                plot_patient_folder_path, 
                f"{plot_name_component_slice_mfg}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}_{plot_name_component_slice}.pdf"
                )
        else:
            plot_b0_slice_path = None
            plot_read_mfg_slice_path = None
            plot_phase_mfg_slice_path = None
            plot_slice_mfg_slice_path = None

        plot_slice_from_3d_image(b0_arr, display_slice_index, plot_brain_mask_arr,
                                title=f"$B_0$ field map (off-resonance)",
                                val_range=b0_range, cmap="PRGn",
                                is_show_colorbar=True, colorbar_label="$Hz$", figsize=figsize_slice_plot,
                                is_show=False, save_path=plot_b0_slice_path)
        plot_slice_from_3d_image(mfg_read_arr, display_slice_index, plot_brain_mask_arr,
                                title=f"MFG in read direction\nPrimary PED = {ped_1_true}",
                                val_range=gy_range, cmap="PiYG",
                                is_show_colorbar=True, colorbar_label=ped_mfg_unit, figsize=figsize_slice_plot,
                                is_show=False, save_path=plot_read_mfg_slice_path)
        plot_slice_from_3d_image(mfg_phase_arr, display_slice_index, plot_brain_mask_arr,
                                title=f"MFG in phase direction\nPrimary PED = {ped_1_true}",
                                val_range=gy_range, cmap="PiYG",
                                is_show_colorbar=True, colorbar_label=ped_mfg_unit, figsize=figsize_slice_plot,
                                is_show=False, save_path=plot_phase_mfg_slice_path)
        plot_slice_from_3d_image(mfg_slice_arr, display_slice_index, plot_brain_mask_arr,
                                title=f"MFG in slice direction\nPrimary PED = {ped_1_true}",
                                val_range=gy_range, cmap="PiYG",
                                is_show_colorbar=True, colorbar_label=ped_mfg_unit, figsize=figsize_slice_plot,
                                is_show=(not is_save_plots), save_path=plot_slice_mfg_slice_path)

    #########################################################################################
    # For each tau, analyze voxelwise relationship between variables and plot it, if chosen #
    #########################################################################################

    for tau_index in evaluate_tau_index_list:
        print(f"Evauating tau {tau_list[tau_index]}")

        # Select the data for one tau
        tau_log_magnitude_ratio = log_magnitude_ratio[..., tau_index]
        tau_echo_shift_difference = echo_shift_difference[..., tau_index]

        masked_log_magnitude_ratio = tau_log_magnitude_ratio[brain_mask_arr]
        masked_echo_shift_difference = tau_echo_shift_difference[brain_mask_arr]

        # Remove nan values from analysis
        masked_echo_shift_difference = masked_echo_shift_difference[~np.isnan(masked_log_magnitude_ratio)]
        masked_mfg_phase_nonnan = masked_mfg_phase_arr[~np.isnan(masked_log_magnitude_ratio)]
        masked_log_magnitude_ratio = masked_log_magnitude_ratio[~np.isnan(masked_log_magnitude_ratio)]


        pearson_result_echo_shift_diff_to_log_magnitude = sp.stats.pearsonr(
            masked_echo_shift_difference, masked_log_magnitude_ratio
            )
        pearson_result_echo_shift_diff_to_mfg = sp.stats.pearsonr(
            masked_echo_shift_difference, masked_mfg_phase_nonnan
            )
        pearson_result_mfg_to_log_magnitude = sp.stats.pearsonr(
            masked_mfg_phase_nonnan, masked_log_magnitude_ratio
            )
        
        pearson_corr_echo_shift_diff_to_log_magnitude = pearson_result_echo_shift_diff_to_log_magnitude.statistic
        pearson_corr_ci_echo_shift_diff_to_log_magnitude = pearson_result_echo_shift_diff_to_log_magnitude.confidence_interval(confidence_level=.95)
        pearson_corr_cilow_echo_shift_diff_to_log_magnitude = pearson_corr_ci_echo_shift_diff_to_log_magnitude.low
        pearson_corr_cihigh_echo_shift_diff_to_log_magnitude = pearson_corr_ci_echo_shift_diff_to_log_magnitude.high
        pearson_corr_echo_shift_diff_to_mfg = pearson_result_echo_shift_diff_to_mfg.statistic
        pearson_corr_ci_echo_shift_diff_to_mfg = pearson_result_echo_shift_diff_to_mfg.confidence_interval(confidence_level=.95)
        pearson_corr_cilow_echo_shift_diff_to_mfg = pearson_corr_ci_echo_shift_diff_to_mfg.low
        pearson_corr_cihigh_echo_shift_diff_to_mfg = pearson_corr_ci_echo_shift_diff_to_mfg.high
        pearson_corr_mfg_to_log_magnitude = pearson_result_mfg_to_log_magnitude.statistic
        pearson_corr_ci_mfg_to_log_magnitude = pearson_result_mfg_to_log_magnitude.confidence_interval(confidence_level=.95)
        pearson_corr_cilow_mfg_to_log_magnitude = pearson_corr_ci_mfg_to_log_magnitude.low
        pearson_corr_cihigh_mfg_to_log_magnitude = pearson_corr_ci_mfg_to_log_magnitude.high

        print(f"Pearson correlation for log_mag_rat to echo_shift_diff [95% CI]: {pearson_corr_echo_shift_diff_to_log_magnitude:.4f} "
              f"[{pearson_corr_cilow_echo_shift_diff_to_log_magnitude:.4f}, {pearson_corr_cihigh_echo_shift_diff_to_log_magnitude:.4f}]")
        print(f"Pearson correlation for mfg to echo_shift_diff [95% CI]: {pearson_corr_echo_shift_diff_to_mfg:.4f} "
              f"[{pearson_corr_cilow_echo_shift_diff_to_mfg:.4f}, {pearson_corr_cihigh_echo_shift_diff_to_mfg:.4f}]")
        print(f"Pearson correlation for log_mag_rat to mfg [95% CI]: {pearson_corr_mfg_to_log_magnitude:.4f} "
              f"[{pearson_corr_cilow_mfg_to_log_magnitude:.4f}, {pearson_corr_cihigh_mfg_to_log_magnitude:.4f}]")

        if not is_skip_plots:
            if is_save_plots:
                plot_corr_echo_shift_diff_to_log_magnitude_path = os.path.join(
                    plot_patient_folder_path, 
                    f"correlation_{plot_name_component_log_sig_rat}_to_{plot_name_component_echo_shift}_tau{tau_list[tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}.pdf"
                    )
                plot_corr_echo_shift_diff_to_mfg_path = os.path.join(
                    plot_patient_folder_path, 
                    f"correlation_{plot_name_component_gy}_to_{plot_name_component_echo_shift}_tau{tau_list[tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}.pdf"
                    )
                plot_corr_mfg_to_log_magnitude_path = os.path.join(
                    plot_patient_folder_path, 
                    f"correlation_{plot_name_component_log_sig_rat}_to_{plot_name_component_gy}_tau{tau_list[tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}.pdf"
                    )
            else:
                plot_corr_echo_shift_diff_to_log_magnitude_path = None
                plot_corr_echo_shift_diff_to_mfg_path = None
                plot_corr_mfg_to_log_magnitude_path = None

            bbox_kwargs = dict(boxstyle="square", facecolor="w")
            text_kwargs = dict(x=0.05, y=0.95, s=f"Pearson Correlation: {pearson_corr_echo_shift_diff_to_log_magnitude:.2f}", 
                            fontsize=12, verticalalignment='top', bbox=bbox_kwargs)
            scatter_plot_2d_with_histograms(masked_echo_shift_difference, masked_log_magnitude_ratio, x_bins=50, y_bins=50,
                                            title=f"Primary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\; ms$",
                                            x_lim=tau_eff_difference_range, y_lim=log_mag_ratio_range,
                                            x_color_lim=tau_eff_difference_range, y_color_lim=log_mag_ratio_range,
                                            x_cmap="RdBu", y_cmap="PuOr",
                                            x_label=f"{tau_sub_label} $\\left[ {tau_sub_unit} \\right]$",
                                            y_label=f"{log_mag_ratio_label}",
                                            text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs, 
                                            x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                            figsize=figsize_scatter_plot, labels_fontsize=fontsize_scatter_plot_labels,
                                            cmap="jet", is_show=False,
                                            save_path=plot_corr_echo_shift_diff_to_log_magnitude_path)

            text_kwargs = dict(x=0.05, y=0.95, s=f"Pearson Correlation: {pearson_corr_echo_shift_diff_to_mfg:.2f}", 
                            fontsize=12, verticalalignment='top', bbox=bbox_kwargs) 
            scatter_plot_2d_with_histograms(masked_echo_shift_difference, masked_mfg_phase_nonnan, x_bins=50, y_bins=50,
                                            title=f"Primary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\; ms$",
                                            x_lim=tau_eff_difference_range, y_lim=gy_range,
                                            x_color_lim=tau_eff_difference_range, y_color_lim=gy_range,
                                            x_cmap="RdBu", y_cmap="PiYG",
                                            x_label=f"{tau_sub_label} $\\left[ {tau_sub_unit} \\right]$",
                                            y_label=f"{ped_mfg_label} $\\left[ {ped_mfg_unit.lstrip('$').rstrip('$')} \\right]$",
                                            text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs,
                                            x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                            figsize=figsize_scatter_plot, labels_fontsize=fontsize_scatter_plot_labels,
                                            cmap="jet", is_show=False,
                                            save_path=plot_corr_echo_shift_diff_to_mfg_path)


            text_kwargs = dict(x=0.05, y=0.95, s=f"Pearson Correlation: {pearson_corr_mfg_to_log_magnitude:.2f}", 
                            fontsize=12, verticalalignment='top', bbox=bbox_kwargs)   
            scatter_plot_2d_with_histograms(masked_mfg_phase_nonnan, masked_log_magnitude_ratio, x_bins=50, y_bins=50,
                                            title=f"Primary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\; ms$",
                                            x_lim=gy_range, y_lim=log_mag_ratio_range,
                                            x_color_lim=gy_range, y_color_lim=log_mag_ratio_range,
                                            x_cmap="PiYG", y_cmap="PuOr",
                                            x_label=f"{ped_mfg_label} $\\left[ {ped_mfg_unit.lstrip('$').rstrip('$')} \\right]$",
                                            y_label=f"{log_mag_ratio_label}",
                                            text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs, 
                                            x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                            figsize=figsize_scatter_plot,  labels_fontsize=fontsize_scatter_plot_labels,
                                            cmap="jet", is_show=False,
                                            save_path=plot_corr_mfg_to_log_magnitude_path)


            """plot_slice_from_3d_image(tau_echo_shift_difference, display_slice_index, plot_brain_mask_arr,
                                    title=f"{tau_sub_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\;ms$",
                                    val_range=tau_eff_difference_range, cmap="RdBu",
                                    is_show_colorbar=True, colorbar_label=tau_sub_unit,
                                    figsize=figsize_slice_plot, is_show=False)
            plot_slice_from_3d_image(tau_log_magnitude_ratio, display_slice_index, plot_brain_mask_arr,
                                    title=f"{log_mag_ratio_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\;ms$",
                                    val_range=log_mag_ratio_range, cmap="PuOr",
                                    is_show_colorbar=True, colorbar_label=log_mag_ratio_unit, 
                                    figsize=figsize_slice_plot, is_show=False)
            plot_slice_from_3d_image(mfg_phase_arr, display_slice_index, plot_brain_mask_arr,
                                    title=f"{ped_mfg_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\;ms$",
                                    val_range=gy_range, cmap="PiYG",
                                    is_show_colorbar=True, colorbar_label=ped_mfg_unit, 
                                    figsize=figsize_slice_plot, is_show=(not is_save_plots))"""


def main():
    #########
    # Setup #
    #########
    
    base_folder_path = "SPECIFY" # Have to define proper location of data for your case
    fase_subfolder_name = "FASE" # In which subfolder is FASE data saved?
    mfg_subfolder_name = "mfg_from_b0" # In which subfolder is MFG data saved?

    tau_list = [0, 6, 12, 18, 24, 30, 36, 42, 48]

    # Define factor that will be multiplied to the read MFG file. In case definition of direction is inconsistent with labels
    # Order: Read, Phase, Slice
    mfg_multiplication_factor_list = [1., 1., -1.] # In our case, the script somehow calculates "IS" instead of "SI", so have to modify slice factor

    # phase-encoding direction (ped) number 1 is the standard direction and 2 the reversed direction (with subscript rev)
    # In our data set, labels of phase encoding direction are reversed. So need to have true label opposite of file label

    patient_settings_dict_list = [ # Have to add your own settings for the plot
        dict(patient_name="SPECIFY", parallel_imaging_file_label="_P2", ped_1_true="PA", ped_2_true="AP",
             ped_1_file_label="_AP", ped_2_file_label="_PA", tau_list=tau_list,
             mfg_multiplication_factor_list=mfg_multiplication_factor_list, display_slice_index=12),
    ]

    is_use_topup_image = True
    evaluate_tau_index_list = None #[-1, 0] # Only those tau indices are used in evaluation. If None, use all

    # Settings for plot
    is_skip_plots = False
    is_save_plots = True

    if is_skip_plots:
        plot_settings_dict = None
    else:
        # Define plot settings
        log_mag_ratio_range = (-0.7, 0.7) # Unitless
        tau_eff_difference_range = (-10, 10) #milliseconds
        gy_range = (-45, 45) # micro T / m
        b0_range = (-60, 60) # Hz

        is_show_only_brain_on_plotted_slices = True # If False, plot whole slice of each parameter, if true, show only values in brain

        scatter_plot_grid_kwargs = dict(which="major", axis="both")
        stairs_kwargs = dict(color="k") # Parameters for the rim above the histograms. If None given, don't draw it

        figsize_slice_plot = (4, 3)
        figsize_scatter_plot = (4.416, 4.416) # (5, 5)
        fontsize_scatter_plot_labels = 11

        tau_sub_label = f"$\\tau_{{eff}}-\\tau_{{eff,rev}}$"
        tau_sub_unit = "ms"
        log_mag_ratio_label = f"$ln \\left(S / S_{{rev}}\\right)$" # f"$ln \left(\\frac{{ S }}{{ S_{{rev}} }} \\right)$"
        log_mag_ratio_unit = ""
        ped_mfg_label = f"$G_{{y}}$"
        ped_mfg_unit = "$\\mu T / m$" # "$\\frac{\\mu T}{m}$"

        plot_folder_path = "SPECIFY"
        plot_name_component_log_sig_rat = "log_sig_rat"
        plot_name_component_echo_shift = "echo_shift_sub"
        plot_name_component_gy = "ped_mfg"
        plot_name_component_read_mfg = "read_mfg"
        plot_name_component_phase_mfg = "phase_mfg"
        plot_name_component_slice_mfg = "slice_mfg"
        plot_name_component_b0 = "b0"

        plot_settings_dict = dict(
            log_mag_ratio_range=log_mag_ratio_range, tau_eff_difference_range=tau_eff_difference_range,
            gy_range=gy_range, b0_range=b0_range,
            is_show_only_brain_on_plotted_slices=is_show_only_brain_on_plotted_slices,
            scatter_plot_grid_kwargs=scatter_plot_grid_kwargs, stairs_kwargs=stairs_kwargs,
            figsize_slice_plot=figsize_slice_plot, figsize_scatter_plot=figsize_scatter_plot,
            fontsize_scatter_plot_labels=fontsize_scatter_plot_labels,
            tau_sub_label=tau_sub_label, tau_sub_unit=tau_sub_unit, log_mag_ratio_label=log_mag_ratio_label,
            log_mag_ratio_unit=log_mag_ratio_unit, ped_mfg_label=ped_mfg_label, ped_mfg_unit=ped_mfg_unit,
            )

        if is_save_plots:
            plot_settings_dict["plot_folder_path"] = plot_folder_path
            plot_settings_dict["plot_name_component_log_sig_rat"] = plot_name_component_log_sig_rat
            plot_settings_dict["plot_name_component_echo_shift"] = plot_name_component_echo_shift
            plot_settings_dict["plot_name_component_gy"] = plot_name_component_gy
            plot_settings_dict["plot_name_component_read_mfg"] = plot_name_component_read_mfg
            plot_settings_dict["plot_name_component_phase_mfg"] = plot_name_component_phase_mfg
            plot_settings_dict["plot_name_component_slice_mfg"] = plot_name_component_slice_mfg
            plot_settings_dict["plot_name_component_b0"] = plot_name_component_b0

    for patient_settings_dict in patient_settings_dict_list:
        # Modify patient dict so it aligns with the input format of the process_patient function.
        input_settings_dict = patient_settings_dict.copy()

        patient_name = patient_settings_dict["patient_name"]
        
        input_folder = os.path.join(base_folder_path, patient_name, fase_subfolder_name)
        mfg_input_folder = os.path.join(base_folder_path, patient_name, mfg_subfolder_name)

        input_settings_dict["input_folder"] = input_folder
        input_settings_dict["mfg_input_folder"] = mfg_input_folder
        input_settings_dict["evaluate_tau_index_list"] = evaluate_tau_index_list
        input_settings_dict["is_use_topup_image"] = is_use_topup_image

        process_patient(input_settings_dict=input_settings_dict, plot_settings_dict=plot_settings_dict)

    print("Done")


if __name__ == "__main__":
    main()