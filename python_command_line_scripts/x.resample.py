"""
Script can be called through the terminal. Enables resampling of images using SimpleITK.
Originally written by Fatmeh annotated by Philipp.
"""

import SimpleITK as sitk
import numpy as np
import argparse
import os

 

def main(path, reference, resultpath, mask=None, maskout=None):
    #Load image
    img = sitk.ReadImage(path)
    ref = sitk.ReadImage(reference)

    #resample image and optionally mask
    resample = sitk.ResampleImageFilter()
    resample.SetOutputDirection(ref.GetDirection())
    resample.SetOutputOrigin(ref.GetOrigin())
    resample.SetInterpolator(sitk.sitkLinear)
    resample.SetTransform(sitk.Transform(3, sitk.sitkIdentity))

    resample.SetSize(ref.GetSize())
    resample.SetOutputSpacing(ref.GetSpacing())

    resample.SetNumberOfThreads(1)

    resampled_img = resample.Execute(img)
    sitk.WriteImage(resampled_img, resultpath)

    if mask is not None:
        mask_image = sitk.ReadImage(mask)
        resampled_mask = resample.Execute(mask_image)
        sitk.WriteImage(resampled_mask, maskout)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Resample an image to 1x1x1 voxels")
    parser.add_argument('--img', type=str, required=True, help="Image path")
    parser.add_argument('--ref', type=str, required=True, help="Image path to resample to")
    parser.add_argument('--mask', type=str, required=False, help="Mask path")
    parser.add_argument('--imgout', type=str, required=True, help="Image output path")
    parser.add_argument('--maskout', type=str, required=False, help="Mask output path")
    args = parser.parse_args()
    main(args.img, args.ref, args.imgout)