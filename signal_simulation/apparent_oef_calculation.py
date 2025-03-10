import numpy as np
from scipy.optimize import curve_fit
import sys
sys.path.append("./")
from utilities.constants_and_helpers import calculate_oef_conversion_constant

def r2prime_func(t, a, r2prime):
    return a * np.exp(-r2prime*t)

def parameter_guess_r2prime(t_arr, s_arr):
    """Provide an initial guess for the least squares fit of the signal decay at the tail of the tau-axis.
    Assume a mono-exponential decay. s(t)=a*exp(-r2prime*t)

    Args:
        t_arr (float array): Time points at which the signal decay tail is evaluated.
        s_arr (float array): Signals obtained at those time points in the signal decay tail. Must have same shape as t_arr.

    Returns:
        tuple (float, float): Guesses for the parameters a and r2prime of monoexponential decay of form a*exp(-r2prime*t).
    """
    #highest_nonnan_index = np.arange() # TODO: Think about implementing robustness to NaN here
    r2prime = (np.log(s_arr[-1])-np.log(s_arr[0]))/(t_arr[0]-t_arr[-1])
    a = s_arr[0]/np.exp(-r2prime*t_arr[0])
    return (a, r2prime)

def fit_r2prime_on_curve_tail(tau_arr, signal_arr, tail_start_tau):
    """ Perform a nonlinear least squares fit of the signal decay at the tail of the tau-axis.
    Assume a mono-exponential decay. s(t)=a*exp(-r2prime*t). Yields a and r2prime.

    Args:
        tau_arr (float arr): Time shift points at which the signal is provided.
        signal_arr (float arr): Signals obtained at those time shift points. Must have same shape as tau_arr.
        tail_start_tau (float): Time beyond which we assume to be in the mono-exponential tail. Is used to estimate r2prime.
        Needs to be in same unit as tau_arr

    Returns:
        tuple (float, float): Fits for the parameters a and r2prime of monoexponential decay of form a*exp(-r2prime*t) when 
        looking at high-tau tail of signal. Unit: 1/ms
    """
    tail_tau_arr = tau_arr[tau_arr >= tail_start_tau]
    tail_signal_arr = signal_arr[tau_arr >= tail_start_tau]

    p0 = parameter_guess_r2prime(tail_tau_arr, tail_signal_arr)

    tail_nonnan_tau_arr = tail_tau_arr[~np.isnan(tail_signal_arr)]
    tail_nonnan_signal_arr = tail_signal_arr[~np.isnan(tail_signal_arr)]

    if tail_nonnan_signal_arr.size < 10:
        return np.array([np.nan, np.nan])

    popt, pcov = curve_fit(r2prime_func, tail_nonnan_tau_arr, tail_nonnan_signal_arr, p0)

    return popt

def calculate_dbv(tau_arr, signal_arr, extrapolated_s0):
    """Use monoexponential signal extrapolated to tau=0 and measured signal to calculate DBV.

    Args:
        tau_arr (float arr): Time shift points at which the signal is provided.
        signal_arr (float arr): Signals obtained at those time shift points. Must have same shape as tau_arr.
        extrapolated_s0 (float): Signal extrapolated to tau=0 based on the monoexponential decay at high tau.

    Returns:
        float: fractional deoxygenated blood volume
    """
    zero_index = np.argmin(np.absolute(tau_arr))
    s0 = signal_arr[zero_index]
    return np.log(extrapolated_s0)-np.log(s0)

def calculate_oef_r2prime_dbv_from_signal_curve(tau_arr, signal_arr, tail_start_tau):
    """Define monoexponential signal tail of signal(tau) curve. Obtain R2' from slope at tail. Extrapolate slope to
    tau=0 and determine DBV from difference to observed intersection

    Args:
        tau_arr (_type_): _description_
        signal_arr (_type_): _description_
        tail_start_tau (_type_): _description_

    Returns:
        _type_: _description_
    """
    extrapolated_s0, r2prime = fit_r2prime_on_curve_tail(tau_arr, signal_arr, tail_start_tau)
    dbv = calculate_dbv(tau_arr, signal_arr, extrapolated_s0)
    oef = calculate_oef_conversion_constant()*r2prime/dbv # calculate_oef_conversion_constant needs r2prime in 1/ms
    r2prime_hz = r2prime*1000. # Convert r2prime from unit 1/ms to unit Hz
    return oef, r2prime_hz, dbv



