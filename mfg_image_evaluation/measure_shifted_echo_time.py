# These scripts implement the simulated signal dropout method to estimate the shifted echo location in k-space.

import numpy as np
import scipy as sp
import sys
sys.path.append('./')
from utilities.fourier_image_processing import fourier_transform_complex_image_arr, inverse_fourier_transform_complex_kspace_arr

def calculate_echo_shift_time(complex_image_arr, read_encode_axis, phase_encode_axis, epi_acquisition_time, gaussian_sigma=None):
    # Calculate at each voxel which k-space index shows the largest signal dip
    signal_dip_index_arr = calculate_signal_dip_index(complex_image_arr, read_encode_axis, phase_encode_axis, gaussian_sigma=gaussian_sigma)

    # Calculate the temporal spacing between phase-encoding k-space lines (usually echo spacing / parallel imaging factor)
    n_phase_encode_steps = complex_image_arr.shape[phase_encode_axis]
    delta_time = epi_acquisition_time/n_phase_encode_steps

    center_kspace_index = np.ceil(n_phase_encode_steps/2)

    # Calculate the time shift from the intended echo time where the k-space center seems to be at each voxel.
    return (signal_dip_index_arr - center_kspace_index) * delta_time


def calculate_signal_dip_index(complex_image_arr, read_encode_axis, phase_encode_axis, gaussian_sigma=None):

    # For testing, TODO: clean up
    import matplotlib.pyplot as plt

    # Fourier transform the complex image into a complex k-space image with the same shape.
    # Apply just a 2D fourier transform along the read and phase encode axes.
    complex_kspace_arr = fourier_transform_complex_image_arr(complex_image_arr, read_encode_axis, phase_encode_axis)

    # Create a new array that will store the signal magnitude of every voxel for every left-out k-space line
    n_phase_encode_steps = complex_image_arr.shape[phase_encode_axis]
    shape_signal_shift_image_arr = (n_phase_encode_steps,) + complex_image_arr.shape

    signal_shift_image_arr = np.zeros(shape_signal_shift_image_arr, dtype=np.float64)

    # For every phase encoding step index, leave out one k-space line, then calculate the signal magnitude at every voxel.
    for signal_shift_index in range(n_phase_encode_steps):
        leave_one_out_complex_kspace_arr = complex_kspace_arr.copy()
        # Set one k-space line in the phase encoding direction to 0
        #leave_one_out_complex_kspace_arr = 
        #np.put_along_axis(leave_one_out_complex_kspace_arr, indices=signal_shift_index, values=0., axis=phase_encode_axis) # Does not work. #TODO: fix it
        if phase_encode_axis == 0:
            leave_one_out_complex_kspace_arr[signal_shift_index] = 0.
        elif phase_encode_axis == 1:
            leave_one_out_complex_kspace_arr[:,signal_shift_index] = 0.
        elif phase_encode_axis == 2:
            leave_one_out_complex_kspace_arr[:,:,signal_shift_index] = 0.
        elif phase_encode_axis == 3:
            leave_one_out_complex_kspace_arr[:,:,:,signal_shift_index] = 0.

        # Transform the modified k-space array back into an image
        leave_one_out_complex_image_arr = inverse_fourier_transform_complex_kspace_arr(leave_one_out_complex_kspace_arr, 
                                                                                       read_encode_axis, phase_encode_axis)
        
        # Write the signal magnitude to the appropriate index of the shift array
        signal_shift_image_arr[signal_shift_index] = np.absolute(leave_one_out_complex_image_arr)

    # Apply gaussian filter to signal shifts along the omitted axis in order to improve the consistency
    if not gaussian_sigma is None: # If no sigma is given, omit this step
        signal_shift_image_arr = sp.ndimage.gaussian_filter1d(signal_shift_image_arr, sigma=gaussian_sigma, axis=0,
                                                              mode="nearest")
    
    # Calculate at which phase encoding index the signal has dipped the lowest by encluding that k-space line.
    signal_dip_index_arr = np.argmin(signal_shift_image_arr, axis=0)

    # Uncomment if you want to debug: Test if the signal dip looks fine
    #if np.size(signal_shift_image_arr.shape) >= 5:
    #    plt.figure()
    #    for i in range(signal_shift_image_arr.shape[-1]):
    #        #plt.plot(np.arange(n_phase_encode_steps) - np.ceil(n_phase_encode_steps/2), signal_shift_image_arr[:, 31, 27, 8, i],
    #        #         label=f"Tau index {i}")
    #        plt.plot(np.arange(n_phase_encode_steps) - np.ceil(n_phase_encode_steps/2), signal_shift_image_arr[:, 26, 45, 19, i],
    #                 label=f"Tau index {i}")
    #    plt.xlabel("omitted k space location")
    #    plt.ylabel("Signal magnitude")
    #    plt.legend()
    #    plt.show()

    return signal_dip_index_arr
