import numpy as np
import matplotlib.pyplot as plt


def show_slice(image_arr, intensity_range, cmap, title=None, scale_description=None):
    plt.figure()
    if not title is None:
        plt.title(title)
    plt.imshow(image_arr.T[::-1], interpolation=None, vmin=intensity_range[0], vmax=intensity_range[-1], cmap=cmap)
    plt.xticks([])
    plt.yticks([])
    cbar = plt.colorbar()
    if not scale_description is None:
        cbar.set_label(scale_description)

def plot_slice_from_3d_image(image_3d_arr, slice_index, mask_3d_arr=None, ax=None, figsize=None,
                             title=None, cmap=None, val_range=None, is_show_colorbar=False, colorbar_label=None,
                             is_transform_nifti_coords=True, is_show=True, save_path=None):
    # Allow ax to be passed, which means users can plot this content in a subplot somewhere.
    # If no ax is given, create it.
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        is_created_new_axis = True
    else:
        is_created_new_axis = False
    if cmap is None:
        cmap = "jet" # Default cmap is jet

    if title is not None:
        ax.set_title(title)

    # Select only the chosen slice from the image
    display_image_arr = image_3d_arr.copy()[:, :, slice_index]
    # If a mask is provided, set all values outside of mask to NaN, so they don't appear on image
    if mask_3d_arr is not None:
        display_image_arr[~(mask_3d_arr[:, :, slice_index])] = np.nan

    # Have to reformat array if it came from nifti, so it looks correct
    if is_transform_nifti_coords:
        display_image_arr = display_image_arr.T[::-1]

    # If no val_range is given, use min and max of array to be displayed.
    if val_range is None:
        val_range = [display_image_arr.nanmin(), display_image_arr.nanmax()]

    ax.imshow(display_image_arr, interpolation=None, vmin=np.amin(val_range), vmax=np.amax(val_range), cmap=cmap)

    if is_show_colorbar:
        colorbar = plt.colorbar(ax.images[0], ax=ax)
        if colorbar_label is not None:
            colorbar.set_label(colorbar_label)

    # Make sure there are no ticks, because we want to just show a clean image
    ax.set_xticks([])
    ax.set_yticks([])

    if is_created_new_axis:
        plt.tight_layout()

    if is_show:
        plt.show()

    if save_path is not None:
        plt.savefig(save_path)
        plt.close()
