"""
Script can be called through the terminal. Creates an appropriate acqparam file to be used for topup.
Originally written by Philipp.
"""
import argparse
import sys
sys.path.append("./")
from utilities.load_json import (load_fase_imageparam_json, PHASE_ENCODING_DIRECTION_AP_STR, PHASE_ENCODING_DIRECTION_RL_STR,
                                 PHASE_ENCODING_DIRECTION_STR, TOTAL_READOUT_TIME_STR)

def create_acq_string(phase_encoding_direction, total_readout_time):
    """Create the acquisition string in the correct format for topup.
    The first three arguments relate to the phase encoding direction, the last argument is the total readout time.
    The phase encoding direction convention adopted here was validated by looking at the estimated B0 output of topup.
    Furthermore, the convention is the same as shown here: https://ftp.nmr.mgh.harvard.edu/pub/dist/freesurfer/tutorial_packages/centos6/fsl_507/doc/wiki/topup(2f)Faq.html

    Args:
        phase_encoding_direction (string): The phase encoding direction in human readable form, i.e. one of these: AP, PA, RL, LR
        total_readout_time (string): Total readout time of the EPI readout (usually defined as echo spacing * (phase encoding steps - 1)).
        This factor should be declared as if there was no partial fourier.

    Returns:
        string: String of four numbers, the first three define the phase encoding direction, the last one the total readout time.
        The values are separated by spaces so they are in the right format for topup.
    """
    num_list = [0, 0, 0, total_readout_time]
    if phase_encoding_direction == PHASE_ENCODING_DIRECTION_AP_STR:
        num_list[1] = -1
    if phase_encoding_direction == PHASE_ENCODING_DIRECTION_AP_STR[::-1]:
        num_list[1] = 1
    if phase_encoding_direction == PHASE_ENCODING_DIRECTION_RL_STR:
        num_list[0] = 1
    if phase_encoding_direction == PHASE_ENCODING_DIRECTION_RL_STR[::-1]:
        num_list[0] = -1
    num_list = [str(num) for num in num_list]
    acq_string = " ".join(num_list)
    return acq_string

def main(path1, outpath, path2=None):
    jsond_dict1 = load_fase_imageparam_json(path1)
    if path2 is not None:
        jsond_dict2 = load_fase_imageparam_json(path2)
    else:
        jsond_dict2 = jsond_dict1.copy()
        jsond_dict2[PHASE_ENCODING_DIRECTION_STR] = jsond_dict1[PHASE_ENCODING_DIRECTION_STR][::-1]
    
    acq_string_1 = create_acq_string(jsond_dict1[PHASE_ENCODING_DIRECTION_STR], jsond_dict1[TOTAL_READOUT_TIME_STR])
    acq_string_2 = create_acq_string(jsond_dict2[PHASE_ENCODING_DIRECTION_STR], jsond_dict2[TOTAL_READOUT_TIME_STR])

    combined_acq_string = "\n".join([acq_string_1, acq_string_2])

    with open(outpath, "w") as f:
        f.write(combined_acq_string)

if __name__ == '__main__':
    # If input 2 is not given, automatically choose the inverse phase encoding direction of input 1.
    parser = argparse.ArgumentParser(description="create")
    parser.add_argument('--paramsin1', type=str, required=True, help="Image parameter json input path, first image")
    parser.add_argument('--paramsin2', type=str, required=False, help="Image parameter json input path, image with reversed phase encoding direction")
    parser.add_argument('--paramsout', type=str, required=True, help="Acquisition parameter txt output path")
    args = parser.parse_args()
    main(args.paramsin1, args.paramsout, args.paramsin2)