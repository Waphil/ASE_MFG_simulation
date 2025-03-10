import os
import numpy as np
import nibabel as nib
import sys
sys.path.append('./')
from utilities.constants_and_helpers import calculate_voxel_size_arr_from_affine

def calculate_central_differences_in_image_axes(image_arr, voxel_size_arr=None):
    """Calculates the central differences of the image along every axis. 
    If voxel_size_arr is specified, the gradient is scaled by that distance in each axis, otherwise a spacing of 1 between voxels is assumed.

    Args:
        image_arr (ndarray): numeric, of any shape. contains image information (e.g. B0 field)
        voxel_size_arr (ndarray, optional): numeric, should be of shape (n,) where n is the number of dimensions of image_arr. 
        Defaults to None, which corresponds to a distance of 1 between voxels.

    Returns:
        tuple with n elements, where n is the number of dimensions of image_arr: Every element is an array of the same shape as
        image_arr, corresponding to the gradient of the image along the respective axis
    """
    if voxel_size_arr is None:
        voxel_size_arr = np.ones(image_arr.ndim)
    #diff_arr_list = [np.gradient(image_arr, voxel_size_arr[i_ax], axis=i_ax) for i_ax in range(image_arr.ndim)] # Old version, clumsier
    diff_arr_list = np.gradient(image_arr, *voxel_size_arr) # New version, elegant usage of the function

    return tuple(diff_arr_list)

def calculate_mfgs(nifti):
    """Calculates macroscopic magnetic field gradients along each axis for given input B0 field nifti using central differences.
    Returns a nifti representing the gradient for every axis.

    Args:
        nifti (nifti): nifti loaded by nibabel representing a B0 field, dimensions can be arbitrary.

    Returns:
        tuple of niftis: number of elements is the number of dimensions of the image and each element is a nifti representing the 
        gradient in that axis.
    """
    image_arr = nifti.get_fdata()
    affine = nifti.affine
    voxel_size_arr = calculate_voxel_size_arr_from_affine(affine)
    mfg_arr_tup = calculate_central_differences_in_image_axes(image_arr, voxel_size_arr)

    mfg_nifti_tup = tuple(nib.Nifti1Image(mfg_arr, affine, nifti.header) for mfg_arr in mfg_arr_tup)
    return mfg_nifti_tup

def load_nifti_and_calculate_mfgs_and_save(in_path, out_path_list=None):
    """Loads the B0 nifti from the specified path, calculates the magnetic field gradients along each image axis and writes those 
    informations into nifti files. If no output paths are specified, they are automatically created in the input folder. In that case,
    the nifti_ending appropriate for the input path may have to be specified.

    Args:
        in_path (string): file path to the input nifti file.
        out_path_list (list of strings, optional): file path to the output nifti files representing the gradients along each axis. 
        Defaults to None, meaning the pahts are automatically generated in the input directory. 
    """
    nifti = nib.load(in_path)

    mfg_nifti_tup = calculate_mfgs(nifti)

    # If no list of output paths are given, create paths automatically.
    if out_path_list is None:
        file_name_start_index = in_path.rfind(os.sep) + 1
        file_ending_start_index = in_path.find(".")
        
        in_file_name = in_path[file_name_start_index:]
        in_file_name_no_ending = in_path[file_name_start_index:file_ending_start_index]
        in_folder_name = in_path[:file_name_start_index]
        file_ending = in_path[file_ending_start_index:]

        out_file_name_list = [f"{in_file_name_no_ending}_gradient_along_axis_{i}{file_ending}" for i in range(len(mfg_nifti_tup))]
        out_path_list = [f"{in_folder_name}{out_file_name}" for out_file_name in out_file_name_list]

    for out_path, mfg_nifti in zip(out_path_list, mfg_nifti_tup):
        nib.save(mfg_nifti, out_path)
