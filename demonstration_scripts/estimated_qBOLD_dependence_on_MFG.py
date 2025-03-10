# This scripts plots the simulated streamlined qBOLD parameters from simulated ASE signal curves subject to a
# phase-encoding diretion MFG. Used to create Figure 6 of the paper.

import numpy as np
import time
import matplotlib.pyplot as plt
import sys
sys.path.append("./")
from visualization.basic_plot import basic_multiline_plot
from signal_simulation.epi_ase_signal_simulation import epi_ase_signal_with_mfgs
from signal_simulation.apparent_oef_calculation import calculate_oef_r2prime_dbv_from_signal_curve
from signal_simulation.qbold_static_dephasing_function import (
    calculate_static_dephasing_regime_signal_decay_with_precalculated_integral,
    calculate_random_cylinder_integral_fast)
from utilities.constants_and_helpers import DEFAULT_GAMMA, calculate_oef_conversion_constant

def main():
    n_steps_precalc = int(1e4)

    # Look at range of gamma Gy / v
    n_gammaGytov = 1001
    gammaGytov_min = -.3#5*200/143.9
    gammaGytov_max = .3#4*200/143.9
    gammaGytov_arr = np.linspace(gammaGytov_min, gammaGytov_max, n_gammaGytov)

    voxel_size_arr = np.array([3.5, 3.5, 2.9])  # mm # Important for magnetic field gradient concerns
    phase_encoding_axis = 1
    read_axis = 0

    te = 0.
    t2 = 100  # ms

    tau = np.linspace(-50, 50, num=1001)
    tail_start_tau = 45. #20. # 45.

    n_phase_0 = 11
    phase_0_min = -0.898 # rad/mm
    phase_0_max = 0.898 # rad/mm
    phase_0_init = 0. # rad/mm
    phase_0_arr = np.linspace(phase_0_min, phase_0_max, num=n_phase_0, endpoint=True)

    # Tissue parameters
    t2_min = 10
    t2_max = 500
    t2_arr = np.logspace(np.log(t2_min), np.log(t2_max), num=7, endpoint=True, base=np.e)

    true_dbv = 0.03
    true_oef = 0.4

    ######################
    # Calculation starts #
    ######################

    true_r2prime = true_dbv * true_oef / calculate_oef_conversion_constant() * 1000 # Have to convert to Hz
    kspace_velocity = 0.077 #1. # Doesn't matter, as we modulate the Gy accordingly # Wrong! It actually matters because of signal dropout
    G_phase_arr = gammaGytov_arr*kspace_velocity/(DEFAULT_GAMMA * 1e-12)

    print(f"Precalculating qBOLD Integral with {n_steps_precalc} steps.")
    start_time = time.time()
    precalc_x_arr = np.linspace(-30, 30, num=n_steps_precalc)
    precalc_integral_arr = calculate_random_cylinder_integral_fast(precalc_x_arr)
    end_time = time.time()
    print(f"Precalculation of Integrals took: {end_time - start_time}")

    static_dephasing_f = calculate_static_dephasing_regime_signal_decay_with_precalculated_integral

    f_kwargs_base = {
        # "dbv" : ,
        # "oef" : ,
        "precalc_x_arr": precalc_x_arr,
        "precalc_integral_arr": precalc_integral_arr,
    }

    def signal_intensities_shiftedecho(dbv, oef, t2, tau, te, kspace_velocity, G_read=0., G_phase=0., G_slice=0.,
                                       phase_0=0.):
        mfg_vec = np.zeros((3,))
        mfg_vec[read_axis] = G_read
        mfg_vec[phase_encoding_axis] = G_phase
        mfg_vec[-1] = G_slice
        f_kwargs = f_kwargs_base.copy()
        f_kwargs["dbv"] = dbv
        f_kwargs["oef"] = oef
        signal_values = epi_ase_signal_with_mfgs(tau, static_dephasing_f, 1 / t2, kspace_velocity, s0=1., te=te,
                                                 mfg_vec=mfg_vec, voxel_size_vec=voxel_size_arr,
                                                 phase_axis=phase_encoding_axis, read_axis=read_axis,
                                                 phase_0=phase_0, f_kwargs=f_kwargs, is_simulate_mfg_dephasing=True)
        return signal_values

    def signal_intensities_shiftedecho_nodephase(dbv, oef, t2, tau, te, kspace_velocity, G_read=0., G_phase=0.,
                                                 G_slice=0., phase_0=0.):
        mfg_vec = np.zeros((3,))
        mfg_vec[read_axis] = G_read
        mfg_vec[phase_encoding_axis] = G_phase
        mfg_vec[-1] = G_slice
        f_kwargs = f_kwargs_base.copy()
        f_kwargs["dbv"] = dbv
        f_kwargs["oef"] = oef
        signal_values = epi_ase_signal_with_mfgs(tau, static_dephasing_f, 1 / t2, kspace_velocity, s0=1., te=te,
                                                 mfg_vec=mfg_vec, voxel_size_vec=voxel_size_arr,
                                                 phase_axis=phase_encoding_axis, read_axis=read_axis,
                                                 phase_0=phase_0, f_kwargs=f_kwargs, is_simulate_mfg_dephasing=False)
        return signal_values

    def signal_intensities_nomfg(dbv, oef, t2, tau, te, kspace_velocity, G_read=0., G_phase=0., G_slice=0.,
                                 phase_0=0.):
        mfg_vec = np.zeros((3,))
        mfg_vec[read_axis] = 0.
        mfg_vec[phase_encoding_axis] = 0.
        mfg_vec[-1] = 0.
        f_kwargs = f_kwargs_base.copy()
        f_kwargs["dbv"] = dbv
        f_kwargs["oef"] = oef
        signal_values = epi_ase_signal_with_mfgs(tau, static_dephasing_f, 1 / t2, kspace_velocity, s0=1., te=te,
                                                 mfg_vec=mfg_vec, voxel_size_vec=voxel_size_arr,
                                                 phase_axis=phase_encoding_axis, read_axis=read_axis,
                                                 phase_0=phase_0, f_kwargs=f_kwargs, is_simulate_mfg_dephasing=True)
        return signal_values

    dbv_result_arr = np.zeros((phase_0_arr.size, gammaGytov_arr.size))
    oef_result_arr = np.zeros((phase_0_arr.size, gammaGytov_arr.size))
    r2prime_result_arr = np.zeros((phase_0_arr.size, gammaGytov_arr.size))

    def helper_func_mfg(G_phase, phase_0=0.):
        signal_arr = signal_intensities_shiftedecho(true_dbv, true_oef, t2, tau, te,
                                                    kspace_velocity, G_read=0., G_phase=G_phase, G_slice=0.,
                                                    phase_0=phase_0)
        measured_oef, measured_r2prime, measured_dbv = calculate_oef_r2prime_dbv_from_signal_curve(tau, signal_arr,
                                                                                                   tail_start_tau)
        return measured_oef, measured_r2prime, measured_dbv

    for i, phase_0 in enumerate(phase_0_arr):
        print(f"phase_0 ({i+1}/{phase_0_arr.size}): {phase_0}")
        # with multiprocessing.Pool() as pool:
        # values = pool.map(helper_func, t_acq_arr)
        # values = pool.starmap(helper_func, [(t_acq, t2) for t_acq in t_acq_arr])
        oef_values, r2prime_values, dbv_values = zip(
            *[helper_func_mfg(G_phase, phase_0) for G_phase in G_phase_arr])

        dbv_result_arr[i] = np.array(dbv_values)
        oef_result_arr[i] = np.array(oef_values)
        r2prime_result_arr[i] = np.array(r2prime_values)

    ########
    # Plot #
    ########
    figsize = (4.416, 3.4) #(4.8, 3.6)#(5.4, 3.9)#(5.2, 3.9)

    colors = plt.cm.jet(np.linspace(0, 1, phase_0_arr.size))
    colors_list = list(colors) + ["k"]

    true_val_list = [true_dbv, true_oef, true_r2prime]
    var_name_list = ["DBV", "OEF", "R_2'"]
    var_unit_list = [None, None, "Hz"]
    y_lim_list = [None, (-0.05, 2.55), None]
    y_tick_spacing_list = [0.02, 0.5, 2.]
    data_arr_list = [dbv_result_arr, oef_result_arr, r2prime_result_arr]

    for true_val, var_name, var_unit, y_lim, y_tick_spacing, data_arr in zip(
            true_val_list, var_name_list, var_unit_list, y_lim_list, y_tick_spacing_list, data_arr_list
    ):
        new_data_arr = np.vstack([data_arr, np.ones(G_phase_arr.size) * true_val])

        linestyles = ["-" for data in new_data_arr]
        linestyles[-1] = "--"

        alphas = [1. for data in new_data_arr]
        alphas[-1] = 0.4

        if phase_0_arr.size == 1:
            label_list = [f"${var_name}$ from\nstreamlined qBOLD\nwith simulated signal"] + [f"True ${var_name}$"]
            legend_title = None
        else:
            label_list = [fr"${phase_0:.2f}$" for phase_0 in phase_0_arr] + [f"True\nvalue"]#[f"True ${var_name}$"]
            legend_title = r"$k_{0,y} \left[\frac{rad}{mm}\right]$"

        #title = f"${var_name}$ as Function of $\\frac{{\\gamma G_y}}{{v}}$.\nSimulation Settings: DBV: {true_dbv}, OEF: {true_oef},\n$\\tau_{{min}}$: ${tau[0]} \\; ms$, $\\tau_{{max}}$: ${tau[-1]} \\; ms$"
        #title = f"Simulated ${var_name}$ Obtained Using Streamlined qBOLD"
        #title = f"Streamlined qBOLD: Simulated ${var_name}$"
        title = None

        legend_kwargs = dict(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, fancybox=True, shadow=True)

        x_label = r"$\gamma G_y / v$"#r"$\frac{\gamma G_y}{v}$"
        y_label = f"Apparent ${var_name}$"
        label_fontsize = 11 #None #12

        if not var_unit is None:
            y_label = f"{y_label} $\\left[ {var_unit} \\right]$"

        grid_kwargs = dict(which="major", axis="both")

        basic_multiline_plot(gammaGytov_arr, new_data_arr, label_list, figsize=figsize, colors=colors_list,
                             linestyles=linestyles, alphas=alphas,
                             title=title, x_label=x_label, y_label=y_label, grid_kwargs=grid_kwargs,
                             y_lim=y_lim, x_tick_major_spacing=0.1, y_tick_major_spacing=y_tick_spacing,
                             labels_fontsize=label_fontsize,
                             x_scale=None, y_scale=None, legend_title=legend_title, legend_kwargs=legend_kwargs)

if __name__ == "__main__":
    main()