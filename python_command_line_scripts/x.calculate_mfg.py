"""
Script can be called through the terminal. Calculates the gradients in all directions for an image (intended for B0 map).
Originally written by Philipp.
"""
import argparse
import sys
sys.path.append("./")
from mfg_image_evaluation.calculate_mfg_from_b0 import load_nifti_and_calculate_mfgs_and_save


def main(in_path, out_path_list=None):
    load_nifti_and_calculate_mfgs_and_save(in_path, out_path_list)


if __name__ == '__main__':
    # If input 2 is not given, automatically choose the inverse phase encoding direction of input 1.
    parser = argparse.ArgumentParser(description="create")
    parser.add_argument('--inpath', type=str, required=True, help="Path to input image")
    parser.add_argument('--outpaths', type=str, required=False, help="Path list to output gradient images.")
    args = parser.parse_args()
    main(args.inpath, args.outpaths)