import os
import json
import nibabel as nib

IN_PLANE_PARALLEL_REDUCTION_FACTOR_STR = "ParallelReductionFactorInPlane"
EFFECTIVE_ECHO_SPACING_STR = "EffectiveEchoSpacing"
ACTUAL_ECHO_SPACING_STR = "ActualEchoSpacing"
PIXEL_BANDWIDTH_PHASE_STR = "PixelBandwidthPhase"
TOTAL_READOUT_TIME_STR = "TotalReadoutTime"
PHASE_ENCODING_DIRECTION_DICOM_STR = "InPlanePhaseEncodingDirectionDICOM"
PHASE_ENCODING_POLARITY_GE_STR = "PhaseEncodingPolarityGE"
PHASE_ENCODING_DIRECTION_STR = "PhaseEncodingDirection"
VOXEL_SIZE_X_STR = "VoxelSizeX"
VOXEL_SIZE_Y_STR = "VoxelSizeY"
MATRIX_SIZE_X_STR = "MatrixSizeX"
MATRIX_SIZE_Y_STR = "MatrixSizeY"
SLICE_NUMBER_STR = "SliceNumber"

PHASE_ENCODING_DIRECTION_COLUMN_STR = "COL"
PHASE_ENCODING_DIRECTION_ROW_STR = "ROW"
PHASE_ENCODING_POLARITY_FLIPPED_STR = "Flipped"
PHASE_ENCODING_POLARITY_UNFLIPPED_STR = "Unflipped"

PHASE_ENCODING_DIRECTION_AP_STR = "AP"
PHASE_ENCODING_DIRECTION_RL_STR = "RL"

TAU_STRING = "tau"
JSON_ENDING_STR = ".json"
FASE_FILENAME_STR = "-FASE"
NIFTI_ENDING_STR = ".nii.gz"

VOXEL_SIZE_ROUND_DIGITS = 3
NIFTI_SLICE_COORDINATE_INDEX = 2

def create_phase_encoding_direction_string(phase_encoding_direction_dicom, phase_encoding_polarity_ge):
    """Creates a string consisting of two letters that indicates the phase encoding direction (e.g. AP, PA, LR, RL).
    Note: It appears that the file names can be wrong about the actually chosen direction.
    The convention adopted with this function has been validated by visual inspection of the images.
    Above the sinuses, there should be a positive B0 inhomogeneity (due to the paramagnetic nature). 
    This should lead to a distortion in the phase encoding direction, as the phase there is "too large".

    Args:
        phase_encoding_direction_dicom (string): dicom string indicating phase encoding direction (e.g. "COL", "ROW")
        phase_encoding_polarity_ge (string): dicom string indicating phase encoding polarity (e.g. "Flipped", "Unflipped")
    """
    if phase_encoding_direction_dicom == PHASE_ENCODING_DIRECTION_COLUMN_STR:
        ped_str = PHASE_ENCODING_DIRECTION_AP_STR
    elif phase_encoding_direction_dicom == PHASE_ENCODING_DIRECTION_ROW_STR:
        ped_str = PHASE_ENCODING_DIRECTION_RL_STR

    if phase_encoding_polarity_ge == PHASE_ENCODING_POLARITY_FLIPPED_STR:
        ped_str = ped_str[::-1]
    return ped_str

def load_fase_imageparam_json(path, is_read_voxel_size_from_corresponding_nifti=False, is_read_tau_from_name=False):
    """Loads a json file as a python dictionary and adds certain elements that contain useful information about the MR acquisition.
    If chosen, try to identify the tau value from the file name.

    Args:
        path (string): path to the json file.
        is_read_voxel_size_from_corresponding_nifti (bool, optional): If True, a nifti file with the same name in the same folder will
        be loaded to get the voxel size. That data is written to the dict as well. Defaults to False.
        is_read_tau_from_name (bool, optional): If True, try to identify tau from the file name. Defaults to False.
    """
    json_dict = load_imageparam_json(path, is_read_voxel_size_from_corresponding_nifti)

    # If want to read tau, split the file name at the _. Then convert the first text after an element containing "-FASE" to a float.
    # If that works, the tau has been identified
    if is_read_tau_from_name:
        tau_half = None

        file_name = path.split(os.sep)[-1]
        split_file_name_list = file_name.rstrip(JSON_ENDING_STR).split("_")
        for i, split_file_name_str in enumerate(split_file_name_list):
            if FASE_FILENAME_STR in split_file_name_str:
                try:
                    tau_half = float(split_file_name_list[i+1])
                except ValueError:
                    pass
        
        if not tau_half is None:
            json_dict[TAU_STRING] = 2 * tau_half
        else:
            # Assume that we have a reverse phase encoding file here, which will have tau 0.
            json_dict[TAU_STRING] = 0.

    return json_dict

def get_pe_polarity(json_dict, path):
    phase_encoding_polarity_ge = json_dict.get(PHASE_ENCODING_POLARITY_GE_STR)
    if phase_encoding_polarity_ge is None:
        file_name = path.split(os.sep)[-1]
        if "PA" in file_name:
            phase_encoding_polarity_ge = PHASE_ENCODING_POLARITY_UNFLIPPED_STR
        else:
            phase_encoding_polarity_ge = PHASE_ENCODING_POLARITY_FLIPPED_STR
    return phase_encoding_polarity_ge

def load_imageparam_json(path, is_read_voxel_size_from_corresponding_nifti=False):
    """Loads a json file as a python dictionary and adds certain elements that contain useful information about the MR acquisition.

    Args:
        path (string): path to the json file.
        is_read_voxel_size_from_corresponding_nifti (bool, optional): If True, a nifti file with the same name in the same folder will
        be loaded to get the voxel size. That data is written to the dict as well. Defaults to False.

    Returns:
        dictionary: A dictionary that contains the sequence parameter names as keys and the corresponding parameters as values.
    """
    json_dict = load_json(path)
    
    # Calculate actual echo spacing and add it to the dict. If no parallel imaging factor is found, use 1.
    in_plane_p_factor = json_dict.get(IN_PLANE_PARALLEL_REDUCTION_FACTOR_STR, 1)
    effective_echo_spacing = json_dict[EFFECTIVE_ECHO_SPACING_STR]
    json_dict[ACTUAL_ECHO_SPACING_STR] = in_plane_p_factor * effective_echo_spacing

    # Calculate the pixel bandwidth in phase encoding direction and add it to the dict
    total_readout_time = json_dict[TOTAL_READOUT_TIME_STR]
    json_dict[PIXEL_BANDWIDTH_PHASE_STR] = 1 / total_readout_time

    # Calculate the human readable phase encoding direction string and add it to the dict
    phase_encoding_direction_dicom = json_dict[PHASE_ENCODING_DIRECTION_DICOM_STR]
    phase_encoding_polarity_ge = get_pe_polarity(json_dict, path)
    phase_encoding_direction = create_phase_encoding_direction_string(phase_encoding_direction_dicom, phase_encoding_polarity_ge)
    json_dict[PHASE_ENCODING_DIRECTION_STR] = phase_encoding_direction

    if is_read_voxel_size_from_corresponding_nifti:
        # Define nifti path
        nifti_path = f"{path.rstrip(JSON_ENDING_STR)}{NIFTI_ENDING_STR}"

        # Load nifti with nibabel
        nifti = nib.load(nifti_path)

        # Read voxel and matrix size parameters from nifti
        voxel_size_x = round(abs(nifti.affine[0, 0]), VOXEL_SIZE_ROUND_DIGITS)
        voxel_size_y = round(abs(nifti.affine[1, 1]), VOXEL_SIZE_ROUND_DIGITS)
        matrix_size_x = nifti.shape[0]
        matrix_size_y = nifti.shape[1]
        slice_number = nifti.shape[NIFTI_SLICE_COORDINATE_INDEX]
        
        json_dict[VOXEL_SIZE_X_STR] = voxel_size_x
        json_dict[VOXEL_SIZE_Y_STR] = voxel_size_y
        json_dict[MATRIX_SIZE_X_STR] = matrix_size_x
        json_dict[MATRIX_SIZE_Y_STR] = matrix_size_y
        json_dict[SLICE_NUMBER_STR] = slice_number


    return json_dict


def load_json(path):
    """Loads a json file and returns the contents as a python dictionary.

    Args:
        path (string): path to the json file.

    Returns:
        dictionary: A dictionary that contains the json parameter names as keys and the corresponding entries as values.
    """
    with open(path) as f:
        json_dict = json.load(f)
    
    return json_dict
