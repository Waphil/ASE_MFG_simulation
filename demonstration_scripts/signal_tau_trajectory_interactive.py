# This scripts plots the simulated ASE signal curves subject to a phase-encoding diretion MFG and lets the user change
# parameters on the fly to see the effects

import numpy as np
import sys
sys.path.append("./")
from signal_simulation.epi_ase_signal_simulation import epi_ase_signal_with_mfgs
from signal_simulation.qbold_static_dephasing_function import \
    (calculate_static_dephasing_regime_signal_decay_with_precalculated_integral,
     calculate_random_cylinder_integral_fast)
from visualization.slider_plot import create_slider_plot, create_slider_dict

def main():

    t2_min = 1.
    t2_max = 3000.
    t2_init = 100.  # 3000 # 70 # ms

    t2prime_min = 0.1
    t2prime_max = 500.
    t2prime_init = 120.  # 0.1 #120 # ms

    epi_readout_time_min = 0.
    epi_readout_time_max = 100.
    epi_readout_time_init = 23.3  # 20.

    te_min = 0.
    te_max = 200.
    te_init = 74.

    dbv_min = 0.
    dbv_max = 0.2
    dbv_init = 0.03

    oef_min = 0.
    oef_max = 1.0
    oef_init = 0.4

    G_read_min = -200
    G_read_max = 200
    G_read_init = 0.  # micro T/ m

    G_phase_min = -200
    G_phase_max = 200
    G_phase_init = 0.  # micro T/ m

    G_slice_min = -200
    G_slice_max = 200
    G_slice_init = 0.  # micro T/ m

    voxel_size_arr = np.array([3.5, 3.5, 2.9])  # mm # Important for magnetic field gradient concerns
    phase_encoding_axis = 1
    read_axis = 0

    phase_0_min = -np.pi / voxel_size_arr[phase_encoding_axis]
    phase_0_max = np.pi / voxel_size_arr[phase_encoding_axis]
    phase_0_init = 0.  # rad/mm

    tau_arr = np.linspace(-50, 50, num=1000)

    precalc_x_arr = np.linspace(-30, 30, num=1000)
    precalc_integral_arr = calculate_random_cylinder_integral_fast(precalc_x_arr)

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

    def signal_intensities_shiftedecho_nodephase(dbv, oef, t2, tau, te, epi_readout_time, G_read=0., G_phase=0.,
                                                 G_slice=0., phase_0=0.):
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
                                                 phase_0=phase_0, f_kwargs=f_kwargs, is_simulate_mfg_dephasing=False)
        return signal_values

    def signal_intensities_nomfg(dbv, oef, t2, tau, te, epi_readout_time, G_read=0., G_phase=0., G_slice=0.,
                                 phase_0=0.):
        kspace_velocity = np.inf if epi_readout_time == 0 else 2 * np.pi / (
                    voxel_size_arr[phase_encoding_axis] * epi_readout_time)
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
                                                 phase_0=0., f_kwargs=f_kwargs, is_simulate_mfg_dephasing=False)
        return signal_values

    slider_dict_list = [
        create_slider_dict("DBV", dbv_min, dbv_max, dbv_init, "dbv"),
        create_slider_dict("OEF", oef_min, oef_max, oef_init, "oef"),
        create_slider_dict("T2", t2_min, t2_max, t2_init, "t2"),
        create_slider_dict("EPI readout time [ms]", epi_readout_time_min, epi_readout_time_max, epi_readout_time_init,
                           "epi_readout_time"),
        create_slider_dict("TE", te_min, te_max, te_init, "te"),
        create_slider_dict("MFG read [microT/m]", G_read_min, G_read_max, G_read_init, "G_read"),
        create_slider_dict("MFG phase [microT/m]", G_phase_min, G_phase_max, G_phase_init, "G_phase"),
        create_slider_dict("MFG slice [microT/m]", G_slice_min, G_slice_max, G_slice_init, "G_slice"),
        create_slider_dict("Phase 0 [rad/mm]", phase_0_min, phase_0_max, phase_0_init, "phase_0"),
    ]
    plot_ax_kwargs = {
        "xlabel": "tau [ms]",
        "ylabel": "Signal intensity [A.U.]",
        "yscale": "log",
    }
    plot_ax_grid_kwargs = {
        "which": "major",
        "axis": "both",
    }
    curve_list = [
        signal_intensities_nomfg,
        signal_intensities_shiftedecho_nodephase,
        signal_intensities_shiftedecho,
    ]
    lineplot_kwargs_list = [
        {"color": "k", "linestyle": "--", "alpha": 0.8, "label": "No MFG"},
        {"color": "b", "linestyle": "-", "alpha": 1.0, "label": "echo shift but no dephase"},
        {"color": "r", "linestyle": "-", "alpha": 1.0, "label": "echo shift and dephase"},
    ]
    create_slider_plot(curve_list, tau_arr, "tau", slider_dict_list,
                       plot_ax_kwargs=plot_ax_kwargs, plot_ax_grid_kwargs=plot_ax_grid_kwargs,
                       lineplot_kwargs_list=lineplot_kwargs_list)

if __name__ == "__main__":
    main()