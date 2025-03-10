import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
import sys
sys.path.append('./')
from visualization.histogram_plot import plot_histogram

def makeColours(val_arr, cmap=None, val_range=None):
    if cmap is None:
        cmap = "jet" # Default cmap is jet

    # If no val_range is given, use min and max of array.
    if val_range is None:
        val_range = [val_arr.min(), val_arr.max()]

    norm = Normalize(vmin=np.amin(val_range), vmax=np.amax(val_range))

    # Create colors for each data point based on the value
    color_arr = cm.ScalarMappable(norm=norm, cmap=cmap).to_rgba(val_arr)

    return color_arr

def scatter_plot_2d_kde_colors(arr_1, arr_2, cmap=None, val_range=None, **kwargs):
    # Code on KDE colors from here: https://stackoverflow.com/questions/19064772/visualization-of-scatter-plots-with-overlapping-points-in-matplotlib

    kde_2d = sp.stats.gaussian_kde(np.vstack([arr_1, arr_2]))
    val_arr = kde_2d.evaluate(np.vstack([arr_1, arr_2]))

    scatter_plot_2d_cmap_colors(arr_1, arr_2, val_arr, cmap=cmap, val_range=val_range, **kwargs)

def scatter_plot_2d_cmap_colors(arr_1, arr_2, val_arr, cmap=None, val_range=None, **kwargs):
    colors = makeColours(val_arr, cmap=cmap, val_range=val_range)

    scatter_plot_2d(arr_1, arr_2, colors=colors, **kwargs)

def scatter_plot_2d(arr_1, arr_2, colors=True, ax=None, figsize=None,
                    title=None, alpha=None, x_label=None, y_label=None, x_lim=None, y_lim=None,
                    labels_fontsize=None,
                    grid_kwargs=None, is_show=True, save_path=None):
    # Allow ax to be passed, which means users can plot this content in a subplot somewhere.
    # If no ax is given, create it.
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        is_created_new_axis = True
    else:
        is_created_new_axis = False

    if alpha is None:
        alpha = 1.

    if title is not None:
        ax.set_title(title)

    ax.scatter(arr_1, arr_2, color=colors, alpha=alpha)

    if x_label is not None:
        ax.set_xlabel(x_label, fontsize=labels_fontsize)
    if y_label is not None:
        ax.set_ylabel(y_label, fontsize=labels_fontsize)
    if grid_kwargs is not None:
        ax.grid(**grid_kwargs)
    if x_lim is not None:
        ax.set_xlim(x_lim)
    if y_lim is not None:
        ax.set_ylim(y_lim)

    if is_created_new_axis:
        plt.tight_layout()

    if is_show:
        plt.show()

    if save_path is not None:
        plt.savefig(save_path)
        plt.close()


def scatter_plot_2d_with_histograms(x_arr, y_arr, x_bins, y_bins, figsize=None, title=None,
                                    x_lim=None, y_lim=None, x_color_lim=None, y_color_lim=None,
                                    x_cmap=None, y_cmap=None, x_edgecolor=None, y_edgecolor=None,
                                    x_stairs_kwargs=None, y_stairs_kwargs=None,
                                    is_show=True, text_kwargs=None, save_path=None, **scatter_kwargs):
    # Plot procedure inspired from here: https://matplotlib.org/stable/gallery/lines_bars_and_markers/scatter_hist.html
    fig = plt.figure(figsize=figsize)
    # Add a gridspec with two rows and two columns and a ratio of 1 to 4 between
    # the size of the marginal Axes and the main Axes in both directions.
    # Also adjust the subplot parameters for a square plot.
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)
    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

    # Plot histograms
    plot_histogram(x_arr, x_bins, is_along_y_axis=False, val_lim=x_lim, cmap=x_cmap, val_color_lim=x_color_lim,
                   ax=ax_histx, edgecolor=x_edgecolor, stairs_kwargs=x_stairs_kwargs, is_show=False)
    plot_histogram(y_arr, y_bins, is_along_y_axis=True, val_lim=y_lim, cmap=y_cmap, val_color_lim=y_color_lim,
                   ax=ax_histy, edgecolor=y_edgecolor, stairs_kwargs=y_stairs_kwargs, is_show=False)

    # Plot scatter plot
    scatter_plot_2d_kde_colors(x_arr, y_arr, ax=ax, is_show=False, **scatter_kwargs)

    # Hide ticks from histogram to make it cleaner
    plt.setp(ax_histx.get_xticklabels(), visible=False)
    plt.setp(ax_histx.get_yticklabels(), visible=False)
    plt.setp(ax_histy.get_xticklabels(), visible=False)
    plt.setp(ax_histy.get_yticklabels(), visible=False)

    # Add text, if given
    if text_kwargs is not None:
        ax.text(transform=ax.transAxes, **text_kwargs)

    # Set limits
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)

    # Add title, if given
    if title is not None:
        plt.suptitle(title)

    gs.tight_layout(fig)

    if is_show:
        plt.show()

    if save_path is not None:
        plt.savefig(save_path)
        plt.close()

def main():
    import numpy.random as random
    import scipy as sp

    x_arr = random.uniform(0., 100., 1000)
    y_arr = random.normal(0., 30., 1000) + x_arr

    pearson_corr = sp.stats.pearsonr(x_arr, y_arr).statistic

    bbox_kwargs = dict(boxstyle="square", facecolor="w")
    text_kwargs = dict(x=0.05, y=0.95, s=f"Pearson Correlation: {pearson_corr:.2f}", fontsize=14,
                       verticalalignment='top', bbox=bbox_kwargs)
    scatter_plot_grid_kwargs = dict(which="major", axis="both")
    stairs_kwargs = dict(color="k")

    scatter_plot_2d_with_histograms(x_arr, y_arr, x_bins=30, y_bins=30, figsize=None, title="Test Plot",
                                    x_lim=None, y_lim=None, x_color_lim=(30, 50), y_color_lim=None,
                                    x_cmap="RdBu", y_cmap="PuOr", #x_edgecolor="black", y_edgecolor="black",
                                    x_label="$x$", y_label="$f(x)$",
                                    text_kwargs=text_kwargs, grid_kwargs=scatter_plot_grid_kwargs,
                                    x_stairs_kwargs=stairs_kwargs, y_stairs_kwargs=stairs_kwargs,
                                    is_show=True, cmap="jet")

if __name__ == "__main__":
    main()
