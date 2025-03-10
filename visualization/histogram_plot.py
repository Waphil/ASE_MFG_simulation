import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm


def compare_histograms(dataframe, value_param, x_axis_differentiator_param=None, 
                       y_axis_differentiator_param=None, hue_differentiator_param=None, bins=None):
    if x_axis_differentiator_param is None:
        unique_x_axes_descriptions = None
        n_x_axes = 1
    else:
        unique_x_axes_descriptions = np.unique(dataframe.loc[:,x_axis_differentiator_param])
        n_x_axes = len(unique_x_axes_descriptions)
    if y_axis_differentiator_param is None:
        unique_y_axes_descriptions = None
        n_y_axes = 1
    else:
        unique_y_axes_descriptions = np.unique(dataframe.loc[:,y_axis_differentiator_param])
        n_y_axes = len(unique_y_axes_descriptions)
    if hue_differentiator_param is None:
        unique_hue_descriptions = None
        n_hues = 1
    else:
        unique_hue_descriptions = np.unique(dataframe.loc[:,hue_differentiator_param])
        n_hues = len(unique_hue_descriptions)
    if bins is None:
        bins = 100
    fig, axes = plt.subplots(n_y_axes, n_x_axes, sharex=True, sharey=True, squeeze=False)

    for i_x in range(n_x_axes):
        # Select only the part of the data frame that corresponds to the currently selected x axis parameter. Repeat with y and hue.
        if not unique_x_axes_descriptions is None:
            x_dataframe = dataframe.loc[dataframe[x_axis_differentiator_param] == unique_x_axes_descriptions[i_x]]
        else:
            x_dataframe = dataframe
        for i_y in range(n_y_axes):
            if not unique_y_axes_descriptions is None:
                y_dataframe = x_dataframe.loc[x_dataframe[y_axis_differentiator_param] == unique_y_axes_descriptions[i_y]]
            else:
                y_dataframe = x_dataframe
            for i_hue in range(n_hues):
                if not unique_hue_descriptions is None:
                    hue_label = unique_hue_descriptions[i_hue]
                    hue_dataframe = y_dataframe.loc[y_dataframe[hue_differentiator_param] == unique_hue_descriptions[i_hue]]
                else:
                    hue_label = None
                    hue_dataframe = y_dataframe
                # Select the values corresponding to the chosen subset of the data frame
                value_arr = hue_dataframe.loc[:,value_param]
                axes[i_y, i_x].hist(value_arr, label=hue_label, alpha=0.5, bins=bins)
                axes[i_y, i_x].legend()
    
    
    # Set axes labels:
    if not unique_x_axes_descriptions is None:
        for i_x in range(n_x_axes):
            axes[0, i_x].set_title(unique_x_axes_descriptions[i_x])
    if not unique_y_axes_descriptions is None:
        for i_y in range(n_y_axes):
            axes[i_y, 0].set_ylabel(unique_y_axes_descriptions[i_y])
    
    plt.show()

def plot_histogram(val_arr, bins, is_along_y_axis=False, val_lim=None, cmap=None, val_color_lim=None, ax=None,
                   figsize=None, title=None, edgecolor=None, stairs_kwargs=None, x_ticks=None, y_ticks=None, is_show=True,
                   save_path=None):
    # Allow ax to be passed, which means users can plot this content in a subplot somewhere.
    # If no ax is given, create it.
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        is_created_new_axis = True
    else:
        is_created_new_axis = False
    if val_lim is None:
        val_lim = (np.nanmin(val_arr), np.nanmax(val_arr))
    if val_color_lim is None:
        val_color_lim = (np.nanmin(val_arr), np.nanmax(val_arr))

    if title is not None:
        ax.set_title(title)
    
    norm = Normalize(np.amin(val_color_lim), np.amax(val_color_lim))

    if is_along_y_axis:
        orientation = "horizontal"
    else:
        orientation = "vertical"

    # If bins is a single value (i.e. number of bins), make sure that this number of bins exists within the limits.
    if np.size(bins) == 1: 
        bins = np.linspace(np.amin(val_lim), np.amax(val_lim), num=(bins+1), endpoint=True)

    N, plotted_bins, patches = ax.hist(val_arr, bins, orientation=orientation, edgecolor=edgecolor)

    if stairs_kwargs is not None:
        ax.stairs(N, plotted_bins, orientation=orientation, **stairs_kwargs)

    if x_ticks is not None:
        ax.set_xticks(x_ticks)
    if y_ticks is not None:
        ax.set_yticks(y_ticks)

    # get correct center locations of bins (for colors)
    plotted_bin_centers = (plotted_bins[:-1] + plotted_bins[1:])/2. 

    # Set color for each patch
    for patch, val in zip(patches, plotted_bin_centers):
        color = cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(val)
        patch.set_facecolor(color)

    if is_created_new_axis:
        plt.tight_layout()
    
    if is_show:
        plt.show()

    if save_path is not None:
        plt.savefig(save_path)
        plt.close()

def main():
    val_arr = np.arange(500) 
    bins = 20
    plot_histogram(val_arr, bins, val_lim=None, cmap=None, val_color_lim=(100, 300), ax=None, figsize=None, title=None, is_show=True)

if __name__ == "__main__":
    main()