import numpy as np
import sys
sys.path.append('./')
from utilities.constants_and_helpers import calculate_gamma_G_product

def simulate_mfg_dephasing_one_direction(gradient, time, voxel_size, initial_phase_gradient=None):
    """Calculates the sinc-shaped dephasing from differential larmor frequency due to a magnetic field gradient 
    within a homogenous voxel with sharp boundaries.

    Args:
        gradient (float): array or scalar: Magnetic field gradient strength in one direction, in microtesla per meter
        time (float): array or scalar: time point after spin echo (in ASE) or after excitation (in GRE) when dephasing is evaluated, in ms.
        voxel_size (float): Size of voxel in direction of gradient, in mm.

    Returns:
        float: Multiplicative intensity modification factor due to dephasing. Unitless. Should be multiplied with signal.
    """
    if initial_phase_gradient is None:
        initial_phase_gradient = np.zeros_like(gradient * time)
    return np.sinc((calculate_gamma_G_product(gradient) * time + initial_phase_gradient) * voxel_size / (2 * np.pi))

def simulate_mfg_dephasing_n_directions(gradient_arr, time, voxel_size_arr, initial_phase_gradient_arr=None): #THIS IS WRONG!  
    """Calculates the sinc-shaped dephasing from differential larmor frequency due to magnetic field gradients 
    within a homogenous voxel with sharp boundaries in n directions.

    Args:
        gradient_arr (float array): Magnetic field gradient strength in each direction, in microtesla per meter.
        Should have shape (n, ...) where n is the number of directions, and the ... can show an arbitrary further shape
        time (float): array or scalar: time point after spin echo (in ASE) or after excitation (in GRE) when dephasing is evaluated, in ms.
        voxel_size_arr (float array): Sizes of voxel in each direction of gradients, in mm. Should have shape (n,).

    Returns:
        float: Multiplicative intensity modification factor due to dephasing. Unitless. Should be multiplied with signal.
    """
    if initial_phase_gradient_arr is None:
        initial_phase_gradient_arr = np.zeros_like(gradient_arr)
    dephasing_factor_list = [simulate_mfg_dephasing_one_direction(gradient, time, voxel_size, initial_phase_gradient)
                             for gradient, voxel_size, initial_phase_gradient
                             in zip(gradient_arr, voxel_size_arr, initial_phase_gradient_arr)]
    return np.prod(dephasing_factor_list, axis=0) #TODO: test

def simulate_mfg_signal_dropout_one_direction(gradient, time, voxel_size, phase_0=0.):
    """Calculate if the signal drops out of k-space for the given situation.

    Args:
        gradient (float): array or scalar: Magnetic field gradient strength in one direction, in microtesla per meter
        time (float): time point after spin echo (in ASE) or after excitation (in GRE) when dephasing is evaluated, in ms.
        voxel_size (float): Size of voxel in direction of gradient, in mm.
        phase_0 (float, optional): Phase of the signal at spin echo in rad/mm. Defaults to 0.

    Returns:
        array in shape of product between gradient and time: Elements are 1 if there is no dropout and nan if there is dropout
    """
    dropout_condition = np.array(np.abs(calculate_gamma_G_product(gradient) * time + phase_0) > np.pi/voxel_size)
    if dropout_condition.size == 1:
        scale_value = np.nan if dropout_condition else 1.
    else:
        scale_value = np.ones_like(dropout_condition, dtype=np.float32)
        scale_value[dropout_condition] = np.nan
    return scale_value
