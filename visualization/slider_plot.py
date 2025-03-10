import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider

SLIDER_DICT_LABEL_LABEL = "label"
SLIDER_DICT_MIN_LABEL = "min_val"
SLIDER_DICT_MAX_LABEL = "max_val"
SLIDER_DICT_INITIAL_LABEL = "initial_val"
SLIDER_DICT_FUNCLABEL_LABEL = "func_label"
SLIDER_DICT_KWARGS_LABEL = "kwargs"

def create_slider_dict(label, min_val, max_val, initial_val, func_label, **kwargs):
    slider_dict = {
        SLIDER_DICT_LABEL_LABEL : label, # The label that should be shown on screen
        SLIDER_DICT_MIN_LABEL : min_val, # Minimum selectable value
        SLIDER_DICT_MAX_LABEL : max_val, # Maximum selectable value
        SLIDER_DICT_INITIAL_LABEL : initial_val, # Initial value
        SLIDER_DICT_FUNCLABEL_LABEL : func_label, # The label that corresponds to the variable name in the plotted function definition
        SLIDER_DICT_KWARGS_LABEL : kwargs, # Keyword arguments to define the format of the slider
    }
    return slider_dict

def calculate_function_y_values(func, x_values, x_values_func_label, param_value_list, param_value_func_label_list):
    input_dict = {
        param_value_func_label : param_value for param_value_func_label, param_value in zip(param_value_func_label_list, param_value_list)
    }
    input_dict[x_values_func_label] = x_values
    y_values = func(**input_dict)
    return y_values

import time

def calculate_function_y_values_using_sliders(func, x_values, x_values_func_label, slider_list, slider_dict_list):
    #start_time = time.time()
    param_value_list = [slider.val for slider in slider_list]
    param_value_func_label_list = [slider_dict.get(SLIDER_DICT_FUNCLABEL_LABEL) for slider_dict in slider_dict_list]
    input_dict = {
        param_value_func_label : param_value for param_value_func_label, param_value in zip(param_value_func_label_list, param_value_list)
    }
    input_dict[x_values_func_label] = x_values
    #end_time = time.time()
    #print(f"Updating of Plots before function evaluation took: {end_time-start_time}")
    #start_time = time.time()
    y_values = func(**input_dict)
    #end_time = time.time()
    #print(f"Updating of Plots function evaluation took: {end_time-start_time}")
    return y_values

def plot_function_with_slider_values(ax, func, x_values, x_values_func_label, slider_list, slider_dict_list, lineplot_kwargs=None):
    y_values = calculate_function_y_values_using_sliders(func, x_values, x_values_func_label, slider_list, slider_dict_list)
    line = ax.plot(x_values, y_values, **lineplot_kwargs)[0]
    return line

def create_slider_plot(func_list, x_values, x_values_func_label, slider_dict_list, 
                       plot_ax_kwargs=None, plot_ax_grid_kwargs=None, lineplot_kwargs_list=None):
    if plot_ax_kwargs is None:
        plot_ax_kwargs = dict()
    if plot_ax_grid_kwargs is None:
        plot_ax_grid_kwargs = dict()
    if lineplot_kwargs_list is None:
        lineplot_kwargs_list = [dict() for func in func_list]

    n_sliders = len(slider_dict_list)

    # Define plot
    fig = plt.figure(layout="constrained")
    # Set up gridspec which will be filled by plot and sliders
    gridspec = fig.add_gridspec(n_sliders, 2)
    # Set up plot axis
    plot_ax = fig.add_subplot(gridspec[:, 0])
    # Set up slider axes
    slider_ax_list = [fig.add_subplot(gridspec[i, 1]) for i in range(n_sliders)]

    # Apply all given axis parameters to the plot axis
    plot_ax.set(**plot_ax_kwargs)
    plot_ax.grid(**plot_ax_grid_kwargs)

    # Define the sliders
    slider_list = []
    for slider_dict, slider_ax in zip(slider_dict_list, slider_ax_list):
        slider = Slider(
            ax=slider_ax,
            label=slider_dict.get(SLIDER_DICT_LABEL_LABEL),
            valmin=slider_dict.get(SLIDER_DICT_MIN_LABEL),
            valmax=slider_dict.get(SLIDER_DICT_MAX_LABEL),
            valinit=slider_dict.get(SLIDER_DICT_INITIAL_LABEL),
            **(slider_dict.get(SLIDER_DICT_KWARGS_LABEL))
        )
        slider_list.append(slider)
    
    # Do plot with initial values
    line_list = [
        plot_function_with_slider_values(plot_ax, func, x_values, x_values_func_label, slider_list, slider_dict_list, 
                                         lineplot_kwargs=lineplot_kwargs) for func, lineplot_kwargs in zip(func_list, lineplot_kwargs_list)
    ]
    
    # Define update function if any sliders are changed
    def update(val):
        y_max_global = 0.
        y_min_global = np.inf
        for func, line in zip(func_list, line_list):
            y_values = calculate_function_y_values_using_sliders(func, x_values, x_values_func_label, slider_list, slider_dict_list)
            line.set_ydata(y_values)
            y_max_global = np.maximum(np.nanmax(y_values), y_max_global)
            y_min_global = np.minimum(np.nanmin(y_values), y_min_global)
        y_span = y_max_global-y_min_global
        plot_ax.set_ylim((y_min_global-0.1*y_span, y_max_global+0.1*y_span))
        fig.canvas.draw_idle()

    # Connect update function to slider change events.
    for slider in slider_list:
        slider.on_changed(update)
    
    plt.show()

