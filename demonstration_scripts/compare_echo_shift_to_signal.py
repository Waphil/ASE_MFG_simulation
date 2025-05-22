# This script plots the slice visualizations and scatter plots of MFGs, logarithmic signal ratios and echo shifts
# for acquired data in nifti format. Used to create figures 7 & 8 in paper.


import os
import numpy as np
import scipy as sp
import pandas as pd
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import sys
sys.path.append('./')
from utilities.constants_and_helpers import convert_mfg_from_Hz_per_mm_to_microtesla_per_m
from visualization.basic_plot import basic_multiline_plot
from visualization.scatter_plot import scatter_plot_2d_kde_colors, scatter_plot_2d_with_histograms
from visualization.image_show_plot import plot_slice_from_3d_image

SLOPE_CONFIDENCE_LEVEL = 0.95

lin_func = lambda x, a, b: a * x + b
lin_func_nointercept = lambda x, a: a * x


# Two-sided inverse Students t-distribution, used to calculate confidence interval of linear regression
# Taken from here: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
tinv = lambda p, df: abs(sp.stats.t.ppf(p / 2, df))


def analyze_slope(x_arr, y_arr, confidence_level=0.95, is_fit_intercept=True):
    """
    Analyze the slope between two variables using a conventional linear least squares fit.
    :param x_arr: ndarray shape (n, ): x-values of data points
    :param y_arr: ndarray shape (n, ): y-values of data points
    :param confidence_level: (optional) float: confidence level for returned confidence interval (default 0.95)
    :param is_fit_intercept: (optional) bool: Whether to fit an intercept or not. If False, intercept is 0.
    :return: tuple of floats: slope, intercept, pearson correlation, pvalue,
    upper confidence interval boundary for slope, lower confidence interval boundary for slope,
    upper confidence interval boundary for intercept, lower confidence interval boundary for intercept
    """
    if is_fit_intercept:
        linregress_result = sp.stats.linregress(x_arr, y_arr)

        slope = linregress_result.slope
        intercept = linregress_result.intercept
        rvalue = linregress_result.rvalue
        pvalue = linregress_result.pvalue

        slope_stderr = linregress_result.stderr
        intercept_stderr = linregress_result.intercept_stderr
    else:
        popt, pcov = sp.optimize.curve_fit(lin_func_nointercept, x_arr, y_arr)

        slope = popt[0]
        intercept = 0.
        rvalue = np.nan
        pvalue = np.nan

        slope_stderr = np.sqrt(pcov[0, 0])
        intercept_stderr = np.nan

    if confidence_level is None:
        return slope, intercept, rvalue, pvalue
    else:
        t_factor = tinv(1 - confidence_level, x_arr.size - 2)  # Number of standard deviations in confidence interval

        slope_confidence_lower = slope - t_factor * slope_stderr
        slope_confidence_upper = slope + t_factor * slope_stderr

        intercept_confidence_lower = intercept - t_factor * intercept_stderr
        intercept_confidence_upper = intercept + t_factor * intercept_stderr

        return (slope, intercept, rvalue, pvalue, slope_confidence_lower, slope_confidence_upper,
                intercept_confidence_lower, intercept_confidence_upper)


def analyze_slope_robust(x_arr, y_arr, confidence_level=0.95):
    """
    Analyze the slope between two variables using a robust Theil-Sen linear fit.
    :param x_arr: ndarray shape (n, ): x-values of data points
    :param y_arr: ndarray shape (n, ): y-values of data points
    :param confidence_level: (optional) float: confidence level for returned confidence interval (default 0.95)
    :return: tuple of floats: slope, intercept, pearson correlation, pvalue,
    upper confidence interval boundary for slope, lower confidence interval boundary for slope,
    upper confidence interval boundary for intercept, lower confidence interval boundary for intercept
    """
    linregress_result = sp.stats.theilslopes(y_arr, x_arr, alpha=confidence_level, method="separate")

    slope = linregress_result.slope
    intercept = linregress_result.intercept
    rvalue = np.nan
    pvalue = np.nan

    if confidence_level is None:
        return slope, intercept, rvalue, pvalue
    else:
        slope_confidence_lower = linregress_result.low_slope
        slope_confidence_upper = linregress_result.high_slope

        intercept_confidence_lower = np.nan
        intercept_confidence_upper = np.nan

        return (slope, intercept, rvalue, pvalue, slope_confidence_lower, slope_confidence_upper,
                intercept_confidence_lower, intercept_confidence_upper)


def process_patient(input_settings_dict, plot_settings_dict=None, output_settings_dict=None):
    """
    For a given patient, load images and other parameter maps, create plots of one slice for each parameter, calculate
    slope between parameters, create scatter plots between parameters and save slope results as csv.
    :param input_settings_dict: dict: defines the input directory and things such as the phase encoding direction,
    parallel imaging factor, taus, and whether the topup corrected data should be used.
    :param plot_settings_dict: (optional) dict: defines the settings for the plots such as units and labels of variables
    and the output folder path. If None is given, skip plots.
    :param output_settings_dict: (optional) dict: Define folder where slope results should be saved as csv. If None is
    given, skip saving slope results.
    """
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

    if output_settings_dict is not None:
        output_csv_folder = output_settings_dict.get("output_csv_folder", None)

        if not os.path.exists(output_csv_folder):
            os.mkdir(output_csv_folder)

        # Define paths for output csvs
        echoshift_name = "echoshift"
        mfg_name = "mfg"
        logsigrat_name = "logsigrat"
        log_magnitude_to_echo_shift_diff_csv_path = os.path.join(
            output_csv_folder,
            f"{logsigrat_name}_to_{echoshift_name}_slope_results_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}.csv"
        )
        echo_shift_diff_to_mfg_csv_path = os.path.join(
            output_csv_folder,
            f"{echoshift_name}_to_{mfg_name}_slope_results_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}.csv"
        )
        log_magnitude_to_mfg_csv_path = os.path.join(
            output_csv_folder,
            f"{logsigrat_name}_to_{mfg_name}_slope_results_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}.csv"
        )
    else:
        output_csv_folder = None

    print("-------------------------------------------------------------------")
    print(f"Evaluating patient from folder: {input_folder}")

    #############
    # Load data #
    #############

    if not is_use_topup_image:
        magnitude_image_file_path_1 = os.path.join(input_folder,
                                                   f"FASE{ped_1_file_label}{parallel_imaging_file_label}_raw.nii.gz")
        magnitude_image_file_path_2 = os.path.join(input_folder,
                                                   f"FASE{ped_2_file_label}{parallel_imaging_file_label}_raw.nii.gz")
        echo_shift_time_file_path_1 = os.path.join(input_folder,
                                                   f"FASE{ped_1_file_label}{parallel_imaging_file_label}_raw_echoshift.nii.gz")
        echo_shift_time_file_path_2 = os.path.join(input_folder,
                                                   f"FASE{ped_2_file_label}{parallel_imaging_file_label}_raw_echoshift.nii.gz")

        log_magnitude_ratio_output_file_path = os.path.join(input_folder,
                                                            f"FASE_{ped_1_file_label}{parallel_imaging_file_label}_raw_phaseencode_log_magnitude_ratio.nii.gz")
    else:
        magnitude_image_file_path_1 = os.path.join(input_folder,
                                                   f"FASE{ped_1_file_label}{parallel_imaging_file_label}.nii.gz")
        magnitude_image_file_path_2 = os.path.join(input_folder,
                                                   f"FASE{ped_2_file_label}{parallel_imaging_file_label}.nii.gz")
        echo_shift_time_file_path_1 = os.path.join(input_folder,
                                                   f"FASE{ped_1_file_label}{parallel_imaging_file_label}_echoshift.nii.gz")
        echo_shift_time_file_path_2 = os.path.join(input_folder,
                                                   f"FASE{ped_2_file_label}{parallel_imaging_file_label}_echoshift.nii.gz")

        log_magnitude_ratio_output_file_path = os.path.join(input_folder,
                                                            f"FASE{ped_1_file_label}{parallel_imaging_file_label}_phaseencode_log_magnitude_ratio.nii.gz")

    brain_mask_path = os.path.join(input_folder,
                                   f"FASE{ped_1_file_label}{parallel_imaging_file_label}_0_aftertopup_brain_mask_ero.nii.gz")

    b0_path = os.path.join(mfg_input_folder,
                           f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}.nii.gz")

    # Load in MFGs from B0 map
    if ped_1_true in ["AP", "PA"]:
        mfg_read_path = os.path.join(mfg_input_folder,
                                     f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_0.nii.gz")
        mfg_phase_path = os.path.join(mfg_input_folder,
                                      f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_1.nii.gz")
        mfg_slice_path = os.path.join(mfg_input_folder,
                                      f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_2.nii.gz")
    else:
        mfg_read_path = os.path.join(mfg_input_folder,
                                     f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_1.nii.gz")
        mfg_phase_path = os.path.join(mfg_input_folder,
                                      f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_0.nii.gz")
        mfg_slice_path = os.path.join(mfg_input_folder,
                                      f"b0map_unwrapped_reg{ped_1_file_label}{parallel_imaging_file_label}_gradient_along_axis_2.nii.gz")

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
    mfg_read_arr = convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg_read_image.get_fdata()) * \
                   mfg_multiplication_factor_list[0]
    mfg_phase_arr = convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg_phase_image.get_fdata()) * \
                    mfg_multiplication_factor_list[1]
    mfg_slice_arr = convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg_slice_image.get_fdata()) * \
                    mfg_multiplication_factor_list[2]

    log_magnitude_ratio = np.log(magnitude_arr_1 / magnitude_arr_2)
    echo_shift_difference = echo_shift_arr_1 - echo_shift_arr_2

    # Test: do sub-selection in brain mask. Maybe that helps us identify problems
    # Threshold for mfg
    # mfg_phase_threshold = 28.78 #microT/m, corresponds to gamma Gy/v of 0.1 in regular case.
    # brain_mask_arr = brain_mask_arr & (np.abs(mfg_phase_arr) < mfg_phase_threshold)
    mfg_slice_threshold = 30  # microT/m
    brain_mask_arr = brain_mask_arr & (np.abs(mfg_slice_arr) < mfg_slice_threshold)

    masked_mfg_read_arr = mfg_read_arr[brain_mask_arr]
    masked_mfg_phase_arr = mfg_phase_arr[brain_mask_arr]
    masked_mfg_slice_arr = mfg_slice_arr[brain_mask_arr]

    # t0_to_tau_slope, t0_to_tau_intercept = analyze_linear_relationship(echo_shift_difference, np.array(tau_list), mask=brain_mask_arr)
    # masked_t0_to_tau_slope = t0_to_tau_slope[brain_mask_arr]

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
        # If plots are not skipped, perform the visualization plots of the variables on one slice per image here.

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

            plot_slice_from_3d_image(echo_shift_difference[:, :, :, evaluate_tau_index], display_slice_index,
                                     plot_brain_mask_arr,
                                     title=f"{tau_sub_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[evaluate_tau_index]} \\;ms$",
                                     val_range=tau_eff_difference_range, cmap="RdBu",
                                     is_show_colorbar=True, colorbar_label=tau_sub_unit, figsize=figsize_slice_plot,
                                     is_show=False,
                                     save_path=plot_echo_shift_diff_slice_path)
            plot_slice_from_3d_image(log_magnitude_ratio[:, :, :, evaluate_tau_index], display_slice_index,
                                     plot_brain_mask_arr,
                                     title=f"{log_mag_ratio_label}\nPrimary PED = {ped_1_true}, $\\tau = {tau_list[evaluate_tau_index]} \\;ms$",
                                     val_range=log_mag_ratio_range, cmap="PuOr",
                                     is_show_colorbar=True, colorbar_label=log_mag_ratio_unit,
                                     figsize=figsize_slice_plot, is_show=False,
                                     save_path=plot_logsigrat_slice_path)
        # show_slice(t0_to_tau_slope[:,:,display_slice_index], (-0.3, 0.3), "PRGn", f"slope of delta t0 to tau")

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

    # Define variable names and variable lists used to save results as csv
    csv_title_string_list = ["tau", "slope", "intercept", "rvalue", "pvalue", "slope_confidence_lower",
                             "slope_confidence_upper",
                             "intercept_confidence_lower", "intercept_confidence_upper"]
    log_magnitude_to_echo_shift_diff_linregress_list = []
    echo_shift_diff_to_mfg_linregress_list = []
    log_magnitude_to_mfg_linregress_list = []

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

        # Calculate pearson correlation between variables.
        pearson_result_log_magnitude_to_echo_shift_diff = sp.stats.pearsonr(
            masked_echo_shift_difference, masked_log_magnitude_ratio
        )
        pearson_result_echo_shift_diff_to_mfg = sp.stats.pearsonr(
            masked_mfg_phase_nonnan, masked_echo_shift_difference
        )
        pearson_result_log_magnitude_to_mfg = sp.stats.pearsonr(
            masked_mfg_phase_nonnan, masked_log_magnitude_ratio
        )

        pearson_corr_log_magnitude_to_echo_shift_diff = pearson_result_log_magnitude_to_echo_shift_diff.statistic
        pearson_corr_ci_log_magnitude_to_echo_shift_diff = pearson_result_log_magnitude_to_echo_shift_diff.confidence_interval(
            confidence_level=.95)
        pearson_corr_cilow_log_magnitude_to_echo_shift_diff = pearson_corr_ci_log_magnitude_to_echo_shift_diff.low
        pearson_corr_cihigh_log_magnitude_to_echo_shift_diff = pearson_corr_ci_log_magnitude_to_echo_shift_diff.high
        pearson_corr_echo_shift_diff_to_mfg = pearson_result_echo_shift_diff_to_mfg.statistic
        pearson_corr_ci_echo_shift_diff_to_mfg = pearson_result_echo_shift_diff_to_mfg.confidence_interval(
            confidence_level=.95)
        pearson_corr_cilow_echo_shift_diff_to_mfg = pearson_corr_ci_echo_shift_diff_to_mfg.low
        pearson_corr_cihigh_echo_shift_diff_to_mfg = pearson_corr_ci_echo_shift_diff_to_mfg.high
        pearson_corr_log_magnitude_to_mfg = pearson_result_log_magnitude_to_mfg.statistic
        pearson_corr_ci_log_magnitude_to_mfg = pearson_result_log_magnitude_to_mfg.confidence_interval(
            confidence_level=.95)
        pearson_corr_cilow_log_magnitude_to_to_mfg = pearson_corr_ci_log_magnitude_to_mfg.low
        pearson_corr_cihigh_log_magnitude_to_mfg = pearson_corr_ci_log_magnitude_to_mfg.high

        print(
            f"Pearson correlation for log_mag_rat to echo_shift_diff [95% CI]: {pearson_corr_log_magnitude_to_echo_shift_diff:.4f} "
            f"[{pearson_corr_cilow_log_magnitude_to_echo_shift_diff:.4f}, {pearson_corr_cihigh_log_magnitude_to_echo_shift_diff:.4f}]")
        print(f"Pearson correlation for mfg to echo_shift_diff [95% CI]: {pearson_corr_echo_shift_diff_to_mfg:.4f} "
              f"[{pearson_corr_cilow_echo_shift_diff_to_mfg:.4f}, {pearson_corr_cihigh_echo_shift_diff_to_mfg:.4f}]")
        print(f"Pearson correlation for log_mag_rat to mfg [95% CI]: {pearson_corr_log_magnitude_to_mfg:.4f} "
              f"[{pearson_corr_cilow_log_magnitude_to_to_mfg:.4f}, {pearson_corr_cihigh_log_magnitude_to_mfg:.4f}]")

        # Perform the robust linear fits for all three variable combinations.
        slope_result_log_magnitude_to_echo_shift_diff = analyze_slope_robust(masked_echo_shift_difference,
                                                                             masked_log_magnitude_ratio,
                                                                             confidence_level=SLOPE_CONFIDENCE_LEVEL)
        slope_result_echo_shift_diff_to_mfg = analyze_slope_robust(masked_mfg_phase_nonnan,
                                                                   masked_echo_shift_difference,
                                                                   confidence_level=SLOPE_CONFIDENCE_LEVEL)
        slope_result_log_magnitude_to_mfg = analyze_slope_robust(masked_mfg_phase_nonnan, masked_log_magnitude_ratio,
                                                                 confidence_level=SLOPE_CONFIDENCE_LEVEL)

        # Append results of linear slope analysis to the lists
        log_magnitude_to_echo_shift_diff_linregress_list.append(
            [tau_list[tau_index]] + list(slope_result_log_magnitude_to_echo_shift_diff))
        echo_shift_diff_to_mfg_linregress_list.append([tau_list[tau_index]] + list(slope_result_echo_shift_diff_to_mfg))
        log_magnitude_to_mfg_linregress_list.append([tau_list[tau_index]] + list(slope_result_log_magnitude_to_mfg))

        if not is_skip_plots:
            # If plots are not skipped, perform the scatter plots here.
            if is_save_plots:
                plot_corr_log_magnitude_to_echo_shift_diff_path = os.path.join(
                    plot_patient_folder_path,
                    f"correlation_{plot_name_component_log_sig_rat}_to_{plot_name_component_echo_shift}_tau{tau_list[tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}.pdf"
                )
                plot_corr_echo_shift_diff_to_mfg_path = os.path.join(
                    plot_patient_folder_path,
                    f"correlation_{plot_name_component_echo_shift}_to_{plot_name_component_gy}_tau{tau_list[tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}.pdf"
                )
                plot_corr_log_magnitude_to_mfg_path = os.path.join(
                    plot_patient_folder_path,
                    f"correlation_{plot_name_component_log_sig_rat}_to_{plot_name_component_gy}_tau{tau_list[tau_index]}_{plot_name_component_true_ped}{plot_name_component_parallel_imaging}.pdf"
                )
            else:
                plot_corr_log_magnitude_to_echo_shift_diff_path = None
                plot_corr_echo_shift_diff_to_mfg_path = None
                plot_corr_log_magnitude_to_mfg_path = None

            # Show slopes on plot. Also show slope as text.
            slope_x_arr = np.linspace(np.amin(masked_echo_shift_difference), np.amax(masked_echo_shift_difference))
            slope_y_arr = slope_result_log_magnitude_to_echo_shift_diff[0] * slope_x_arr + \
                          slope_result_log_magnitude_to_echo_shift_diff[1]
            bbox_kwargs = dict(boxstyle="square", facecolor="w")
            text_kwargs = dict(x=0.05, y=0.95,
                               s=f"slope: ${slope_result_log_magnitude_to_echo_shift_diff[0]:.3f} \\left[ 1/ms \\right]$",
                               fontsize=12, verticalalignment='top', bbox=bbox_kwargs)
            ax, _, _ = scatter_plot_2d_with_histograms(masked_echo_shift_difference, masked_log_magnitude_ratio,
                                                       x_bins=50, y_bins=50,
                                                       title=f"Primary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\; ms$",
                                                       x_lim=tau_eff_difference_range, y_lim=log_mag_ratio_range,
                                                       x_color_lim=tau_eff_difference_range,
                                                       y_color_lim=log_mag_ratio_range,
                                                       x_cmap="RdBu", y_cmap="PuOr",
                                                       x_label=f"{tau_sub_label} $\\left[ {tau_sub_unit} \\right]$",
                                                       y_label=f"{log_mag_ratio_label}",
                                                       text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs,
                                                       x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                                       figsize=figsize_scatter_plot,
                                                       labels_fontsize=fontsize_scatter_plot_labels,
                                                       cmap="jet", is_show=False,
                                                       save_path=None)  # plot_corr_log_magnitude_to_echo_shift_diff_path
            ax.plot(slope_x_arr, slope_y_arr, color="k", marker="")
            plt.savefig(plot_corr_log_magnitude_to_echo_shift_diff_path)
            plt.close()

            slope_x_arr = np.linspace(np.amin(masked_mfg_phase_nonnan), np.amax(masked_mfg_phase_nonnan))
            slope_y_arr = slope_result_echo_shift_diff_to_mfg[0] * slope_x_arr + slope_result_echo_shift_diff_to_mfg[1]
            # text_kwargs = dict(x=0.05, y=0.95, s=f"Pearson Correlation: {pearson_corr_echo_shift_diff_to_mfg:.2f}",
            # text_kwargs = dict(x=0.05, y=0.95, s=f"$\\rho={pearson_corr_echo_shift_diff_to_mfg:.2f}$, $a={slope_result_echo_shift_diff_to_mfg[0]:.4f}$",
            #                fontsize=12, verticalalignment='top', bbox=bbox_kwargs)
            text_kwargs = dict(x=0.05, y=0.95,
                               s=f"slope: ${slope_result_echo_shift_diff_to_mfg[0]:.3f} \\left[ ms \\; m/\\mu T \\right]$",
                               fontsize=12, verticalalignment='top', bbox=bbox_kwargs)
            ax, _, _ = scatter_plot_2d_with_histograms(masked_mfg_phase_nonnan, masked_echo_shift_difference, x_bins=50,
                                                       y_bins=50,
                                                       title=f"Primary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\; ms$",
                                                       x_lim=gy_range, y_lim=tau_eff_difference_range,
                                                       x_color_lim=gy_range, y_color_lim=tau_eff_difference_range,
                                                       x_cmap="PiYG", y_cmap="RdBu",
                                                       x_label=f"{ped_mfg_label} $\\left[ {ped_mfg_unit.lstrip('$').rstrip('$')} \\right]$",
                                                       y_label=f"{tau_sub_label} $\\left[ {tau_sub_unit} \\right]$",
                                                       text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs,
                                                       x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                                       figsize=figsize_scatter_plot,
                                                       labels_fontsize=fontsize_scatter_plot_labels,
                                                       cmap="jet", is_show=False,
                                                       save_path=None)  # plot_corr_echo_shift_diff_to_mfg_path
            ax.plot(slope_x_arr, slope_y_arr, color="k", marker="")
            plt.savefig(plot_corr_echo_shift_diff_to_mfg_path)
            plt.close()

            slope_x_arr = np.linspace(np.amin(masked_mfg_phase_nonnan), np.amax(masked_mfg_phase_nonnan))
            slope_y_arr = slope_result_log_magnitude_to_mfg[0] * slope_x_arr + slope_result_log_magnitude_to_mfg[1]
            # text_kwargs = dict(x=0.05, y=0.95, s=f"Pearson Correlation: {pearson_corr_log_magnitude_to_mfg:.2f}",
            # text_kwargs = dict(x=0.05, y=0.95, s=f"$\\rho={pearson_corr_log_magnitude_to_mfg:.2f}$, $a={slope_result_log_magnitude_to_mfg[0]:.4f}$",
            #                fontsize=12, verticalalignment='top', bbox=bbox_kwargs)
            text_kwargs = dict(x=0.05, y=0.95,
                               s=f"slope: ${slope_result_log_magnitude_to_mfg[0]:.4f} \\left[ m/\\mu T \\right]$",
                               fontsize=12, verticalalignment='top', bbox=bbox_kwargs)
            ax, _, _ = scatter_plot_2d_with_histograms(masked_mfg_phase_nonnan, masked_log_magnitude_ratio, x_bins=50,
                                                       y_bins=50,
                                                       title=f"Primary PED = {ped_1_true}, $\\tau = {tau_list[tau_index]} \\; ms$",
                                                       x_lim=gy_range, y_lim=log_mag_ratio_range,
                                                       x_color_lim=gy_range, y_color_lim=log_mag_ratio_range,
                                                       x_cmap="PiYG", y_cmap="PuOr",
                                                       x_label=f"{ped_mfg_label} $\\left[ {ped_mfg_unit.lstrip('$').rstrip('$')} \\right]$",
                                                       y_label=f"{log_mag_ratio_label}",
                                                       text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs,
                                                       x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                                       figsize=figsize_scatter_plot,
                                                       labels_fontsize=fontsize_scatter_plot_labels,
                                                       cmap="jet", is_show=False,
                                                       save_path=None)  # plot_corr_log_magnitude_to_mfg_path
            ax.plot(slope_x_arr, slope_y_arr, color="k", marker="")
            plt.savefig(plot_corr_log_magnitude_to_mfg_path)
            plt.close()

    if output_csv_folder is not None:
        # Convert the slope fit results to pandas data frames and write them as csvs.
        log_magnitude_to_echo_shift_diff_df = pd.DataFrame(log_magnitude_to_echo_shift_diff_linregress_list,
                                                           columns=csv_title_string_list)
        echo_shift_diff_to_mfg_df = pd.DataFrame(echo_shift_diff_to_mfg_linregress_list, columns=csv_title_string_list)
        log_magnitude_to_mfg_df = pd.DataFrame(log_magnitude_to_mfg_linregress_list, columns=csv_title_string_list)

        log_magnitude_to_echo_shift_diff_df.to_csv(log_magnitude_to_echo_shift_diff_csv_path, ";", index=False)
        echo_shift_diff_to_mfg_df.to_csv(echo_shift_diff_to_mfg_csv_path, ";", index=False)
        log_magnitude_to_mfg_df.to_csv(log_magnitude_to_mfg_csv_path, ";", index=False)


def plot_relationship_slopes_for_taus(base_folder_path, patient_settings_dict_list, csv_subfolder_name):
    """
    Visualize the dependence of the fitted slopes between the variables on tau for each investigated image series.
    :param base_folder_path: string: path to where the data is located.
    :param patient_settings_dict_list: list of dicts: Each entry contains information on the patient settings.
    :param csv_subfolder_name: float: Name of the subfolder where csv files are saved.
    """
    parallel_imaging_linestyle_dict = {"_P2": "-", "_P4": "--"}
    patient_marker_dict = {"volunteer01": "o", "volunteer01_session2": "o", "volunteer02": "x"}
    patient_color_dict = {"volunteer01": np.array([90. / 255., 196. / 255., 223. / 255.]),
                          "volunteer01_session2": np.array([90. / 255., 196. / 255., 223. / 255.]),
                          # "volunteer02" : np.array([147./255., 198./255., 63./255.])}
                          "volunteer02": np.array([19. / 255., 51. / 255., 155. / 255.])}
    ped_color_dict = {"AP": "b", "PA": "b", "RL": "r", "LR": "r"}
    ped_marker_dict = {"AP": "o", "PA": "o", "RL": "x", "LR": "x"}
    patient_label_dict = {"volunteer01": "Volunteer 1", "volunteer01_session2": None, "volunteer02": "Volunteer 2"}
    ped_label_dict = {"AP": "PA-AP", "PA": "PA-AP", "RL": "RL-LR", "LR": "RL-LR"}
    parallel_imaging_label_dict = {"_P2": "P=2", "_P4": "P=4"}

    figsize = (4.416, 3.4)
    grid_kwargs = dict(which="major", axis="both")
    fontsize_plot_labels = 11

    # Do plot of slope data for each patient
    fig_log_magnitude_to_echo_shift_diff, ax_log_magnitude_to_echo_shift_diff = plt.subplots(1, 1, layout="constrained",
                                                                                             figsize=figsize)
    fig_echo_shift_diff_to_mfg, ax_echo_shift_diff_to_mfg = plt.subplots(1, 1, layout="constrained", figsize=figsize)
    fig_log_magnitude_to_mfg, ax_log_magnitude_to_mfg = plt.subplots(1, 1, layout="constrained", figsize=figsize)

    for patient_settings_dict in patient_settings_dict_list:
        patient_name = patient_settings_dict["patient_name"]
        parallel_imaging_file_label = patient_settings_dict["parallel_imaging_file_label"]
        ped_1_true = patient_settings_dict["ped_1_true"]
        ped_2_true = patient_settings_dict["ped_2_true"]

        output_csv_folder = os.path.join(base_folder_path, patient_name, csv_subfolder_name)

        # Define paths for output csvs
        echoshift_name = "echoshift"
        mfg_name = "mfg"
        logsigrat_name = "logsigrat"
        log_magnitude_to_echo_shift_diff_csv_path = os.path.join(
            output_csv_folder,
            f"{logsigrat_name}_to_{echoshift_name}_slope_results_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}.csv"
        )
        echo_shift_diff_to_mfg_csv_path = os.path.join(
            output_csv_folder,
            f"{echoshift_name}_to_{mfg_name}_slope_results_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}.csv"
        )
        log_magnitude_to_mfg_csv_path = os.path.join(
            output_csv_folder,
            f"{logsigrat_name}_to_{mfg_name}_slope_results_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}.csv"
        )

        log_magnitude_to_echo_shift_diff_df = pd.read_csv(log_magnitude_to_echo_shift_diff_csv_path, sep=";")
        echo_shift_diff_to_mfg_df = pd.read_csv(echo_shift_diff_to_mfg_csv_path, sep=";")
        log_magnitude_to_mfg_df = pd.read_csv(log_magnitude_to_mfg_csv_path, sep=";", )

        x_variable = "tau"
        y_variable = "slope"  # "rvalue" #"slope"
        y_uncertainty_variable = "slope_confidence"  # None #"slope_confidence"
        if y_variable == "slope":
            y_unit_string_log_magnitude_to_echo_shift_diff = "$\\left[ 1/ms \\right]$"
            y_unit_string_echo_shift_diff_to_mfg = "$\\left[ ms \\; m/\\mu T \\right]$"
            y_unit_string_log_magnitude_to_mfg = "$\\left[ m/\\mu T \\right]$"
        elif y_variable == "intercept":
            y_unit_string_log_magnitude_to_echo_shift_diff = ""
            y_unit_string_echo_shift_diff_to_mfg = "$\\left[ ms \\right]$"
            y_unit_string_log_magnitude_to_mfg = ""
        else:
            y_unit_string_log_magnitude_to_echo_shift_diff = ""
            y_unit_string_echo_shift_diff_to_mfg = ""
            y_unit_string_log_magnitude_to_mfg = ""

        tau_arr = log_magnitude_to_echo_shift_diff_df[x_variable].values
        echo_log_magnitude_to_echo_shift_diff_variable_arr = log_magnitude_to_echo_shift_diff_df[y_variable].values
        echo_shift_diff_to_mfg_variable_arr = echo_shift_diff_to_mfg_df[y_variable].values
        log_magnitude_to_mfg_variable_arr = log_magnitude_to_mfg_df[y_variable].values
        if not y_uncertainty_variable is None:
            y_uncertainty_variable_lower = f"{y_uncertainty_variable}_lower"
            y_uncertainty_variable_upper = f"{y_uncertainty_variable}_upper"
            if y_uncertainty_variable in log_magnitude_to_echo_shift_diff_df:
                print(
                    f"Warning: confidence {y_uncertainty_variable} is saved here with old format (i.e. just one value rather than lower and upper.)")
                log_magnitude_to_echo_shift_diff_variable_err_arr = log_magnitude_to_echo_shift_diff_df[
                    y_uncertainty_variable].values
            else:
                y_err_lower = echo_log_magnitude_to_echo_shift_diff_variable_arr - log_magnitude_to_echo_shift_diff_df[
                    y_uncertainty_variable_lower].values
                y_err_upper = log_magnitude_to_echo_shift_diff_df[
                                  y_uncertainty_variable_upper].values - echo_log_magnitude_to_echo_shift_diff_variable_arr
                log_magnitude_to_echo_shift_diff_variable_err_arr = np.stack((y_err_lower, y_err_upper), axis=0)
            if y_uncertainty_variable in echo_shift_diff_to_mfg_df:
                print(
                    f"Warning: confidence {y_uncertainty_variable} is saved here with old format (i.e. just one value rather than lower and upper.)")
                echo_shift_diff_to_mfg_variable_err_arr = echo_shift_diff_to_mfg_df[y_uncertainty_variable].values
            else:
                y_err_lower = echo_shift_diff_to_mfg_variable_arr - echo_shift_diff_to_mfg_df[
                    y_uncertainty_variable_lower].values
                y_err_upper = echo_shift_diff_to_mfg_df[
                                  y_uncertainty_variable_upper].values - echo_shift_diff_to_mfg_variable_arr
                echo_shift_diff_to_mfg_variable_err_arr = np.stack((y_err_lower, y_err_upper), axis=0)
            if y_uncertainty_variable in log_magnitude_to_mfg_df:
                print(
                    f"Warning: confidence {y_uncertainty_variable} is saved here with old format (i.e. just one value rather than lower and upper.)")
                log_magnitude_to_mfg_variable_err_arr = log_magnitude_to_mfg_df[y_uncertainty_variable].values
            else:
                y_err_lower = log_magnitude_to_mfg_variable_arr - log_magnitude_to_mfg_df[
                    y_uncertainty_variable_lower].values
                y_err_upper = log_magnitude_to_mfg_df[
                                  y_uncertainty_variable_upper].values - log_magnitude_to_mfg_variable_arr
                log_magnitude_to_mfg_variable_err_arr = np.stack((y_err_lower, y_err_upper), axis=0)
        else:
            log_magnitude_to_echo_shift_diff_variable_err_arr = None
            echo_shift_diff_to_mfg_variable_err_arr = None
            log_magnitude_to_mfg_variable_err_arr = None

        # Define labels and colors of values
        x_label = "$\\tau$ $\\left[ ms \\right]$"
        y_label_log_magnitude_to_echo_shift_diff = f"{y_variable} $ln \\left(S / S_{{rev}}\\right)$ to $\\tau_{{eff}}-\\tau_{{eff,rev}}$ {y_unit_string_log_magnitude_to_echo_shift_diff}"
        y_label_echo_shift_diff_to_mfg = f"{y_variable} $\\tau_{{eff}}-\\tau_{{eff,rev}}$ to $G_{{y}}$ {y_unit_string_echo_shift_diff_to_mfg}"
        y_label_log_magnitude_to_mfg = f"{y_variable} $ln \\left(S / S_{{rev}}\\right)$ to $G_{{y}}$ {y_unit_string_log_magnitude_to_mfg}"

        label = f"{patient_name.split('_')[0]}_{ped_1_true}_{ped_2_true}{parallel_imaging_file_label}"
        linestyle = parallel_imaging_linestyle_dict[parallel_imaging_file_label]
        color = patient_color_dict[patient_name]
        marker = ped_marker_dict[ped_1_true]

        if log_magnitude_to_echo_shift_diff_variable_err_arr is not None:
            log_magnitude_to_echo_shift_diff_variable_err_arr = np.expand_dims(
                log_magnitude_to_echo_shift_diff_variable_err_arr, axis=0)
            # log_magnitude_to_echo_shift_diff_variable_err_arr = log_magnitude_to_echo_shift_diff_variable_err_arr.reshape((1, -1))
        if echo_shift_diff_to_mfg_variable_err_arr is not None:
            echo_shift_diff_to_mfg_variable_err_arr = np.expand_dims(echo_shift_diff_to_mfg_variable_err_arr, axis=0)
            # echo_shift_diff_to_mfg_variable_err_arr = echo_shift_diff_to_mfg_variable_err_arr.reshape((1, -1))
        if log_magnitude_to_mfg_variable_err_arr is not None:
            log_magnitude_to_mfg_variable_err_arr = np.expand_dims(log_magnitude_to_mfg_variable_err_arr, axis=0)
            # log_magnitude_to_mfg_variable_err_arr = log_magnitude_to_mfg_variable_err_arr.reshape((1, -1))

        basic_multiline_plot(tau_arr, echo_log_magnitude_to_echo_shift_diff_variable_arr.reshape((1, -1)), [label],
                             x_err_arr=None, y_err_arr=log_magnitude_to_echo_shift_diff_variable_err_arr,
                             ax=ax_log_magnitude_to_echo_shift_diff, figsize=None, colors=[color],
                             linestyles=[linestyle], markers=[marker], alphas=None,
                             title=None, x_label=x_label, y_label=y_label_log_magnitude_to_echo_shift_diff,
                             grid_kwargs=None, ticklabel_kwargs=None,
                             is_use_scalar_formatter=False, x_tick_major_spacing=None, y_tick_major_spacing=None,
                             x_tick_minor_spacing=None, y_tick_minor_spacing=None,
                             x_scale=None, y_scale=None, x_lim=None, y_lim=None, labels_fontsize=fontsize_plot_labels,
                             legend_title=None, legend_kwargs=None, is_show=False)
        basic_multiline_plot(tau_arr, echo_shift_diff_to_mfg_variable_arr.reshape((1, -1)), [label],
                             x_err_arr=None, y_err_arr=echo_shift_diff_to_mfg_variable_err_arr,
                             ax=ax_echo_shift_diff_to_mfg, figsize=None, colors=[color], linestyles=[linestyle],
                             markers=[marker], alphas=None,
                             title=None, x_label=x_label, y_label=y_label_echo_shift_diff_to_mfg, grid_kwargs=None,
                             ticklabel_kwargs=None,
                             is_use_scalar_formatter=False, x_tick_major_spacing=10, y_tick_major_spacing=None,
                             x_tick_minor_spacing=None, y_tick_minor_spacing=None,
                             x_scale=None, y_scale=None, x_lim=None, y_lim=None, labels_fontsize=fontsize_plot_labels,
                             legend_title=None, legend_kwargs=None, is_show=False)
        basic_multiline_plot(tau_arr, log_magnitude_to_mfg_variable_arr.reshape((1, -1)), [label],
                             x_err_arr=None, y_err_arr=log_magnitude_to_mfg_variable_err_arr,
                             ax=ax_log_magnitude_to_mfg, figsize=None, colors=[color], linestyles=[linestyle],
                             markers=[marker], alphas=None,
                             title=None, x_label=x_label, y_label=y_label_log_magnitude_to_mfg, grid_kwargs=None,
                             ticklabel_kwargs=None,
                             is_use_scalar_formatter=False, x_tick_major_spacing=10, y_tick_major_spacing=None,
                             x_tick_minor_spacing=None, y_tick_minor_spacing=None,
                             x_scale=None, y_scale=None, x_lim=None, y_lim=None, labels_fontsize=fontsize_plot_labels,
                             legend_title=None, legend_kwargs=None, is_show=False)

        print("ok")

    unique_patient_name_list = []
    unique_ped_1_true_list = []
    unique_parallel_imaging_file_label_list = []
    for patient_settings_dict in patient_settings_dict_list:
        patient_name = patient_settings_dict["patient_name"]
        ped_1_true = patient_settings_dict["ped_1_true"]
        parallel_imaging_file_label = patient_settings_dict["parallel_imaging_file_label"]
        if patient_name not in unique_patient_name_list:
            unique_patient_name_list.append(patient_name)
        if ped_1_true not in unique_ped_1_true_list:
            unique_ped_1_true_list.append(ped_1_true)
        if parallel_imaging_file_label not in unique_parallel_imaging_file_label_list:
            unique_parallel_imaging_file_label_list.append(parallel_imaging_file_label)

    # Define custom legend to make it pretty

    legend_handles_patients = [
        mpatches.Patch(color=patient_color_dict[patient_name],
                       label=patient_label_dict[patient_name])
        for patient_name in unique_patient_name_list if patient_label_dict[patient_name] is not None
    ]
    legend_handles_peds = [
        Line2D([0], [0], color="gray", marker=ped_marker_dict[ped_1_true], linestyle="",
               label=ped_label_dict[ped_1_true]) for ped_1_true in unique_ped_1_true_list
    ]
    legend_handles_parallel_imaging = [
        Line2D([0], [0], color="gray", marker="",
               linestyle=parallel_imaging_linestyle_dict[parallel_imaging_file_label],
               label=parallel_imaging_label_dict[parallel_imaging_file_label]) for parallel_imaging_file_label in
        unique_parallel_imaging_file_label_list
    ]

    legend_handles = legend_handles_patients + legend_handles_peds + legend_handles_parallel_imaging

    ax_log_magnitude_to_echo_shift_diff.legend(handles=legend_handles, framealpha=1.)
    ax_echo_shift_diff_to_mfg.legend(handles=legend_handles, framealpha=1.)
    ax_log_magnitude_to_mfg.legend(handles=legend_handles, framealpha=1.)

    # Add grid

    ax_log_magnitude_to_echo_shift_diff.grid(**grid_kwargs)
    ax_echo_shift_diff_to_mfg.grid(**grid_kwargs)
    ax_log_magnitude_to_mfg.grid(**grid_kwargs)
    plt.show()

def main():
    #########
    # Setup #
    #########
    
    base_folder_path = "SPECIFY" # Have to define proper location of data for your case
    fase_subfolder_name = "FASE" # In which subfolder is FASE data saved?
    mfg_subfolder_name = "mfg_from_b0" # In which subfolder is MFG data saved?
    csv_subfolder_name = "slope_csv" # In which subfolder are slope results saved as csv?

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
    is_skip_data_analysis = True # If you don't want to repeat the linear regressions etc and just plot the slopes based
    # on the csvs. Particularly useful for robust fits, as Theil-Sen is not fast.
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

    if not is_skip_data_analysis:
        # Run script for each patient
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

            output_csv_folder = os.path.join(base_folder_path, patient_name, csv_subfolder_name)
            output_settings_dict = dict(output_csv_folder=output_csv_folder)

            process_patient(input_settings_dict=input_settings_dict, plot_settings_dict=plot_settings_dict,
                            output_settings_dict=output_settings_dict)

    plot_relationship_slopes_for_taus(base_folder_path, patient_settings_dict_list, csv_subfolder_name)

    print("Done")


if __name__ == "__main__":
    main()