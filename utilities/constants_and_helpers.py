import numpy as np

# Values defined like in Cherukara et al, 2019. doi: 10.1016/j.neuroimage.2019.116106
DEFAULT_B0 = 3.0 # T
DEFAULT_DELTA_CHI = 0.264*1e-6 # Unitless
DEFAULT_HCT = 0.40 # Unitless
DEFAULT_GAMMA = 2*np.pi*42.576*1e6 # In rad Hz / T

def calculate_oef_conversion_constant(b0=DEFAULT_B0, delta_chi=DEFAULT_DELTA_CHI, Hct=DEFAULT_HCT, gamma=DEFAULT_GAMMA):
    """Calculates the conversion factor between OEF and R2prime/DBV. Is recurring in many places of the qBOLD model.

    Args:
        b0 (float, optional): Magnetic field strength in T. Defaults to DEFAULT_B0.
        delta_chi (float, optional): Susceptibility Difference between tissue and deoxyhaemoglobin contained in blood vessels. 
        Defaults to DEFAULT_DELTA_CHI.
        Hct (float, optional): Fractional Hematocrit. Defaults to DEFAULT_HCT.
        gamma (float, optional): Gyromagnetic ratio of protons in rad * Hz / T. Defaults to DEFAULT_GAMMA. 

    Returns:
        float: Conversion factor between OEF and R2prime/DBV, in unit ms
    """
    return 1/(4./3.*np.pi*gamma*b0*delta_chi*Hct) * 1000 # Convert to miliseconds

def calculate_tc(oef, b0=DEFAULT_B0, delta_chi=DEFAULT_DELTA_CHI, Hct=DEFAULT_HCT, gamma=DEFAULT_GAMMA): 
    """Calculate the caracteristic time used in the qBOLD model. Also expressed as 1/(delta omega)

    Args:
        oef (float): Oxygen Extraction Fraction
        b0 (float, optional): Magnetic field strength in T. Defaults to DEFAULT_B0.
        delta_chi (float, optional): Susceptibility Difference between tissue and deoxyhaemoglobin contained in blood vessels. 
        Defaults to DEFAULT_DELTA_CHI.
        Hct (float, optional): Fractional Hematocrit. Defaults to DEFAULT_HCT.
        gamma (float, optional): Gyromagnetic ratio of protons in rad * Hz / T. Defaults to DEFAULT_GAMMA. 

    Returns:
        float: characteristic time as function of OEF in units ms
    """
    return calculate_oef_conversion_constant(b0, delta_chi, Hct, gamma)/oef

def calculate_epi_kspace_velocity(n_phase_encoding_steps, echo_spacing, parallel_imaging_factor, k_max):
    """Calculates the k-space velocity for an EPI readout with the given parameters.

    Args:
        n_phase_encoding_steps (integer): Number of phase encoding steps to be acquired. Shold not include effect of parallel imaging or partial fourier (that is handled separetely).
        echo_spacing (float): Time between readout of subsequent phase-encoding lines during the EPI readout. Should not include effect of parallel imaging.
        parallel_imaging_factor (float): Parallel imaging acceleration factor.
        k_max (float): The largest k-space value that is written (generally pi/voxel_size)

    Returns:
        float: k-space velocity value. Unit is [unit of k_max] / [unit of echo_spacing]
    """
    if echo_spacing == 0:
        k_space_velocity = np.inf # Use this to simulate infinite k-space velocity
    else:
        k_space_velocity = 2 * k_max * parallel_imaging_factor / (n_phase_encoding_steps * echo_spacing)
    return k_space_velocity

def calculate_epi_readout_duration(n_phase_encoding_steps, echo_spacing, parallel_imaging_factor):
    """Calculates the time duration of the EPI readout. Warning: Partial Fourier not yet supported

    Args:
        n_phase_encoding_steps (integer): Number of phase encoding steps to be acquired. Shold not include effect of parallel imaging or partial fourier (that is handled separetely).
        echo_spacing (float): Time between readout of subsequent phase-encoding lines during the EPI readout. Should not include effect of parallel imaging.
        parallel_imaging_factor (float): Parallel imaging acceleration factor.

    Returns:
        float: Total readout time of EPI readout.
    """
    return 2./calculate_epi_kspace_velocity(n_phase_encoding_steps, echo_spacing, parallel_imaging_factor, 1.)


def calculate_gamma_G_product(G_phase):
    """Calculate product of gamma and G with the correct units.

    Args:
        G_phase (float): Magnetic field gradient strength in microtesla / m

    Returns:
        float: Product of gamma and G_phase in units of rad / (mm ms), which is the same unit as k-space velocity
    """
    product = DEFAULT_GAMMA * G_phase * 1e-12
    return product

def convert_mfg_from_Hz_per_mm_to_microtesla_per_m(mfg):
    """Converts a magnetic field gradient value from Hz per mm to microtesla per m

    Args:
        mfg (float): Magnetic field gradient strength in Hz / mm

    Returns:
        float: Magnetic field gradient strength in microtesla / m
    """
    converted_mfg = mfg * 2 * np.pi / DEFAULT_GAMMA * 1e9
    return converted_mfg

def calculate_voxel_size_arr_from_affine(affine_mat):
    """Calculate voxel spacing in each image direction from the affine matrix.

    Args:
        affine_mat (ndarray): numeric array of shape (4, 4), representing the affine matrix of a nifti.

    Returns:
        ndarray: numeric array of shape (3,), representing the size of a voxel in all three spatial dimensions.
    """
    return np.linalg.norm(affine_mat[:-1,:-1], axis=0)