import scipy as sp

def fourier_transform_complex_image_arr(complex_image_arr, read_encode_axis, phase_encode_axis):
    complex_kspace_arr = sp.fft.fftshift(sp.fft.fftn(complex_image_arr, axes=(read_encode_axis, phase_encode_axis)),
                                         axes=(read_encode_axis, phase_encode_axis))

    return complex_kspace_arr

def inverse_fourier_transform_complex_kspace_arr(complex_kspace_arr, read_encode_axis, phase_encode_axis):
    complex_image_arr = sp.fft.ifftn(sp.fft.ifftshift(complex_kspace_arr, axes=(read_encode_axis, phase_encode_axis)),
                                     axes=(read_encode_axis, phase_encode_axis))

    return complex_image_arr


