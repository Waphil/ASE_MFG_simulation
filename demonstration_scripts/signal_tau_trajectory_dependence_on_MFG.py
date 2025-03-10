# This scripts plots the simulated ASE signal curves subject to a phase-encoding diretion MFG. Used to create
# Figure 5 of the paper.

import numpy as np
from signal_simulation.epi_ase_signal_simulation import epi_ase_signal_with_mfgs
from signal_simulation.qbold_static_dephasing_function import (
    calculate_random_cylinder_integral_fast,
    calculate_static_dephasing_regime_signal_decay_with_precalculated_integral
)
from visualization.basic_plot import basic_multiline_plot
import matplotlib.pyplot as plt

def main():
    #######################
    # Simulation Settings #
    #######################

    n_steps_precalc = int(1e4)

    voxel_size_arr = np.array([3.5, 3.5, 2.9])  # mm # Important for magnetic field gradient concerns
    phase_encoding_axis = 1
    read_axis = 0

    te = 74.
    tau = np.linspace(-50, 50, num=10000)

    true_dbv = 0.03 #0.
    true_oef = 0.4

    #What if we simulate t2 = t2' ?
    #true_r2prime = true_dbv * true_oef / calculate_oef_conversion_constant()
    #t2 = 1 / true_r2prime

    n_G_phase_simulations = 11
    G_phase_min = -50  # -200 # micro T/ m
    G_phase_max = 50  # 200 # micro T/ m
    G_phase_arr = np.linspace(G_phase_min, G_phase_max, num=n_G_phase_simulations, endpoint=True)

    n_phase_0_simulations = 11
    phase_0_min = -np.pi/voxel_size_arr[phase_encoding_axis]#-2*np.pi/voxel_size_arr[phase_encoding_axis] # rad/mm
    phase_0_max = np.pi/voxel_size_arr[phase_encoding_axis]#2*np.pi/voxel_size_arr[phase_encoding_axis] # rad/mm
    phase_0_arr = np.linspace(phase_0_min, phase_0_max, num=n_phase_0_simulations, endpoint=True)

    # Define which scenarios should be simulated for G phase
    G_phase_scenario_list = [
        # (t2[ms], epi_readout_time [ms], g_phase_arr[microT/m], phase_0_for_G_phase_plot[rad/mm])
        (100, 23.3, G_phase_arr, 0.),
        (np.inf, 23.3, G_phase_arr, 0.),
        (100, 23.3/2, G_phase_arr, 0.),
        (100, 23.3, G_phase_arr, 0.5),
        (100, 23.3, G_phase_arr, -0.5),
    ]
    # Define which scenarios should be simulated for phase_0
    phase_0_scenario_list = [
        # (t2[ms], epi_readout_time [ms], phase_0_arr[rad/mm], G_phase_for_phase_0_plot[microT/m])
        (100, 23.3, phase_0_arr, 0.),
        (np.inf, 23.3, phase_0_arr, 0.),
        (100, 23.3/2, phase_0_arr, 0.),
        (100, 23.3, phase_0_arr, 50),
        (100, 23.3, phase_0_arr, -50),
    ]

    #######################
    # Prepare Simulations #
    #######################

    print(f"Precalculating qBOLD Integral with {n_steps_precalc} steps.")
    import time
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

    def signal_intensities_shiftedecho(dbv, oef, t2, tau, te, epi_readout_time, G_read=0., G_phase=0., G_slice=0.,
                                       phase_0=0.):
        kspace_velocity = np.inf if epi_readout_time == 0 else 2 * np.pi / (
                    voxel_size_arr[phase_encoding_axis] * epi_readout_time)
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


    ###################
    # Run Simulations #
    ###################

    G_phase_signal_arr_list = [] # Store the simulated signal arrays for each scenario here
    for t2, epi_readout_time, g_phase_arr, phase_0_for_G_phase_plot in G_phase_scenario_list:
        k_space_velocity = 2 * np.pi / (epi_readout_time * voxel_size_arr[phase_encoding_axis]) if not epi_readout_time == 0 else np.inf

        G_phase_signal_arr = np.zeros((G_phase_arr.size, tau.size))
        for i, G_phase in enumerate(G_phase_arr):
            G_phase_signal_arr[i, :] = signal_intensities_shiftedecho(true_dbv, true_oef, t2, tau, te,
                                                              epi_readout_time, G_read=0., G_phase=G_phase,
                                                              G_slice=0., phase_0=phase_0_for_G_phase_plot)
        G_phase_signal_arr_list.append(G_phase_signal_arr)

    phase_0_signal_arr_list = [] # Store the simulated signal arrays for each scenario here
    for t2, epi_readout_time, phase_0_arr, G_phase_for_phase_0_plot in phase_0_scenario_list:
        k_space_velocity = 2 * np.pi / (epi_readout_time * voxel_size_arr[phase_encoding_axis]) if not epi_readout_time == 0 else np.inf

        phase_0_signal_arr = np.zeros((phase_0_arr.size, tau.size))
        for i, phase_0 in enumerate(phase_0_arr):
            phase_0_signal_arr[i, :] = signal_intensities_shiftedecho(true_dbv, true_oef, t2, tau, te,
                                                              epi_readout_time, G_read=0., G_phase=G_phase_for_phase_0_plot,
                                                              G_slice=0., phase_0=phase_0)
        phase_0_signal_arr_list.append(phase_0_signal_arr)

    ########
    # Plot #
    ########

    for (t2, epi_readout_time, g_phase_arr, phase_0_for_G_phase_plot), G_phase_signal_arr in zip(G_phase_scenario_list,
                                                                                                 G_phase_signal_arr_list):
        k_space_velocity = 2 * np.pi / (epi_readout_time * voxel_size_arr[phase_encoding_axis]) if not epi_readout_time == 0 else np.inf

        colors = plt.cm.jet(np.linspace(0, 1, G_phase_arr.size))

        colors_list = list(colors) + ["k"]
        linestyles = ["-" for data in G_phase_signal_arr]

        alphas = [1. for data in G_phase_signal_arr]

        label_list = [f"{G_phase:.0f}" for G_phase in G_phase_arr]

        """title = (f"Signal ($\\tau$) curve for different $G_{{y}}$ values"
                 f"\nSimulation Settings: $DBV: {true_dbv}$, $OEF: {true_oef}$, $T_2: {t2:.0f}\\:ms$,"
                 f"\n$v: {k_space_velocity:.3f} \\frac{{rad}}{{ms \\: mm}}$, "
                 f"$\\Delta y: {voxel_size_arr[phase_encoding_axis]:.1f}mm$, "
                 f"$k_{{0,y}}: {phase_0_for_G_phase_plot:.1f} \\frac{{rad}}{{mm}}$")"""
        print(f"---------------------------------------------------------------------------------------"
              f"\nSimulation Settings: $DBV: {true_dbv}$, $OEF: {true_oef}$, $T_2: {t2:.0f}\\:ms$,"
              f"\n$v: {k_space_velocity:.3f} \\frac{{rad}}{{ms \\: mm}}$, "
              f"$\\Delta y: {voxel_size_arr[phase_encoding_axis]:.1f}mm$, "
              f"$k_{{0,y}}: {phase_0_for_G_phase_plot:.1f} \\frac{{rad}}{{mm}}$")

        figsize = (4.416, 3.128) #(4.8, 3.4)#(5.4, 3.9)#(5.2, 3.9)

        title = None
        legend_title = r"$G_y \left[\frac{\mu T}{m}\right]$"

        legend_kwargs = dict(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, fancybox=True, shadow=True)

        x_label = r"$\tau \; \left[ms\right]$"
        y_label = r"Signal $\left[A.U.\right]$"
        label_fontsize = 11 #None #12

        grid_kwargs = dict(which="both", axis="both")
        ticklabel_kwargs = dict(axis="both", scilimits=(-3, 3))
        x_tick_major_spacing = 20.
        y_tick_major_spacing = 0.02
        y_tick_minor_spacing = 100.

        basic_multiline_plot(tau, G_phase_signal_arr, label_list, figsize=figsize, colors=colors_list,
                             linestyles=linestyles, alphas=alphas,
                             title=title, x_label=x_label, y_label=y_label, grid_kwargs=grid_kwargs,
                             ticklabel_kwargs=ticklabel_kwargs, is_use_scalar_formatter=True,
                             labels_fontsize=label_fontsize, x_tick_major_spacing=x_tick_major_spacing,
                             y_tick_major_spacing=y_tick_major_spacing, y_tick_minor_spacing=y_tick_minor_spacing,
                             y_scale="log", x_scale=None, legend_title=legend_title, legend_kwargs=legend_kwargs)


    for (t2, epi_readout_time, phase_0_arr, G_phase_for_phase_0_plot), phase_0_signal_arr in zip(phase_0_scenario_list,
                                                                                                 phase_0_signal_arr_list):
        k_space_velocity = 2 * np.pi / (epi_readout_time * voxel_size_arr[phase_encoding_axis]) if not epi_readout_time == 0 else np.inf

        rss_signal_arr = np.sqrt(np.nansum(phase_0_signal_arr*phase_0_signal_arr, axis=0))
        #rmss_signal_arr = np.sqrt(np.nanmean(phase_0_signal_arr*phase_0_signal_arr, axis=0))
        rms_signal_arr = np.sqrt(np.nansum(phase_0_signal_arr*phase_0_signal_arr, axis=0)/phase_0_signal_arr.shape[0])
        new_phase_0_signal_arr = np.vstack([phase_0_signal_arr, rms_signal_arr])

        colors = plt.cm.jet(np.linspace(0, 1, phase_0_arr.size))

        colors_list = list(colors) + ["k"]
        linestyles = ["-" for data in phase_0_signal_arr] + ["--"]

        alphas = [1. for data in phase_0_signal_arr] + [0.8]

        label_list = [fr"${phase_0:.2f}$" for phase_0 in phase_0_arr] + [None]#["Root\nMean\nSquare"]

        """title = (f"Signal ($\\tau$) curve for different $k_{{0,y}}$ values"
                 f"\nSimulation Settings: $DBV: {true_dbv}$, $OEF: {true_oef}$, $T_2: {t2:.0f}\\:ms$,"
                 f"\n$v: {k_space_velocity:.3f} \\frac{{rad}}{{ms \\: mm}}$, "
                 f"$\\Delta y: {voxel_size_arr[phase_encoding_axis]:.1f}mm$, "
                 f"$G_{{y}}: {G_phase_for_phase_0_plot:.0f} \\frac{{\\mu T}}{{m}}$")"""
        print(f"---------------------------------------------------------------------------------------"
              f"\nSimulation Settings: $DBV: {true_dbv}$, $OEF: {true_oef}$, $T_2: {t2:.0f}\\:ms$,"
              f"\n$v: {k_space_velocity:.3f} \\frac{{rad}}{{ms \\: mm}}$, "
              f"$\\Delta y: {voxel_size_arr[phase_encoding_axis]:.1f}mm$, "
              f"$G_{{y}}: {G_phase_for_phase_0_plot:.0f} \\frac{{\\mu T}}{{m}}$")

        figsize = (4.6, 3.128) #(5.0, 3.4)  # (5.4, 3.9)#(5.2, 3.9)

        title = None
        legend_title = r"$k_{0,y} \left[\frac{rad}{mm}\right]$"

        legend_kwargs = dict(loc='center left', bbox_to_anchor=(1.02, 0.5), ncol=1, fancybox=True, shadow=True)

        x_label = r"$\tau \; \left[ms\right]$"
        y_label = r"Signal $\left[A.U.\right]$"
        label_fontsize = 11 #None#12

        grid_kwargs = dict(which="both", axis="both")
        ticklabel_kwargs = dict(axis="both", scilimits=(-3, 3))
        x_tick_major_spacing = 20.
        y_tick_major_spacing = 0.04
        y_tick_minor_spacing = 100.

        basic_multiline_plot(tau, new_phase_0_signal_arr, label_list, figsize=figsize, colors=colors_list,
                             linestyles=linestyles, alphas=alphas,
                             title=title, x_label=x_label, y_label=y_label, grid_kwargs=grid_kwargs,
                             ticklabel_kwargs=ticklabel_kwargs, is_use_scalar_formatter=True,
                             labels_fontsize=label_fontsize, x_tick_major_spacing=x_tick_major_spacing,
                             y_tick_major_spacing=y_tick_major_spacing, y_tick_minor_spacing=y_tick_minor_spacing,
                             x_scale=None, y_scale="log", legend_title=legend_title, legend_kwargs=legend_kwargs)

if __name__ == "__main__":
    main()