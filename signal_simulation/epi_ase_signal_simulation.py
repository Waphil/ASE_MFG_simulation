# Simulate the signal trajectory of the EPI ASE sequence, including MFGs. But not including PSF considerations.

import numpy as np
import sys
sys.path.append("./")
from signal_simulation.echo_shift_simulation import calculate_tau_eff
from signal_simulation.mfg_dephasing_simulation import simulate_mfg_dephasing_one_direction, simulate_mfg_signal_dropout_one_direction

def epi_ase_signal_with_mfgs(tau, f, r2, kspace_velocity, s0=1., te=0., mfg_vec=None, voxel_size_vec=None, phase_axis=0, read_axis=1,
                             phase_0=0., f_kwargs=None, is_simulate_mfg_dephasing=True):
    if mfg_vec is not None: # If an mfg vector is defined, calculate the tau_eff based on the mfg in phase encoding direction
        tau_eff = calculate_tau_eff(tau, kspace_velocity, G_phase=mfg_vec[phase_axis], phase_0=phase_0)
    else:
        tau_eff = tau
    if is_simulate_mfg_dephasing and mfg_vec is not None and voxel_size_vec is not None:
        # If macroscopic dephasing should be simulated and the appropriate quantities are given, the dephasing factor is calculated.
        # Dephasing should only be done in slice-direction with this formalism, as the echo shift accurately models in-plane dephasing.
        # Determine slice axis
        slice_axis_list = [i for i in range(mfg_vec.size) if (i != read_axis) & (i != phase_axis)]
        if len(slice_axis_list) > 1:
            raise ValueError(f"Provided too many dimensions ({mfg_vec.size}), definition of slice axis is unclear.")
        else:
            slice_axis = slice_axis_list[0]
        dephasing_factor = simulate_mfg_dephasing_one_direction(gradient=mfg_vec[slice_axis], time=tau_eff,
                                                                voxel_size=voxel_size_vec[slice_axis], initial_phase_gradient=None)
        # Could add initial phase gradient in slice-direction
        # Remove negative and zero dephasing, as we would see no signal. On a log scale, these values could also cause problems
        dephasing_factor[dephasing_factor <= 0] = np.nan 
    else:
        dephasing_factor = np.ones_like(tau_eff)
    # For voxels where the echo has been shifted outside of the sampled k-space plane (in read and phase direction), set the the signal to nan
    if mfg_vec is not None and voxel_size_vec is not None:
        dephasing_factor *= simulate_mfg_signal_dropout_one_direction(mfg_vec[read_axis], tau_eff, voxel_size_vec[read_axis], phase_0=0.)
        dephasing_factor *= simulate_mfg_signal_dropout_one_direction(mfg_vec[phase_axis], tau_eff, voxel_size_vec[phase_axis], phase_0=phase_0)
    return s0 * np.exp(-r2 * te) * modified_f_with_echo_shift(tau, f, r2, tau_eff=tau_eff, f_kwargs=f_kwargs) * dephasing_factor

def modified_f_with_echo_shift(tau, f, r2, tau_eff=None, f_kwargs=None):
    """Simulate the signal behaviour around the spin echo (f function) taking into account the echo shift impact on r2 decay.

    Args:
        tau (float): array or scalar: the nominal time difference from the spin echo when the signal is read out.
        f (function): function that takes as input a time t and additional arguments f_kwargs, which describes the signal behavior around
        the spin echo. For example in MRI an f(t) = exp(-R2' abs(t)) is often assumed.
        r2 (float): rate of irreversible transversal signal decay in 1/ms.
        tau_eff (float, optional): array or scalar: the actual time difference from the spin echo when center of kspace is read out.
        Should have same shape as tau. Defaults to tau.
        f_kwargs (dict, optional): keyword arguments for f. Depend on the input of f. Defaults to None.

    Returns:
        array or scalar: the relative signal observed at the different nominal tau values, taking into account the echo shift. same shape as tau.
    """
    if tau_eff is None:
        tau_eff = tau

    return np.exp(-r2 * (tau_eff - tau)) * f(t=tau_eff, **f_kwargs)
