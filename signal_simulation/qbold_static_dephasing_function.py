import warnings
import numpy as np
from scipy import integrate # Numerical integration
from scipy.special import j0 # Bessel Function of order 0
import multiprocessing # Library to speed up calculation
import sys
sys.path.append("./")
from utilities.constants_and_helpers import calculate_tc

def calculate_random_cylinder_integral(x):
    """Calculates the integral over u for the qBOLD signal corresponding to randomly oriented,
    infinitely long cylinders with uniform magnetic susceptiblity as introduced by 
    Yablonskiy & Haacke (1994).

    Args:
        x (float): x=t/tc, time shift of the signal compared to the spin echo location divided by 
        the characteristic time tc. (sometimes expressed instead as t * delta omega)

    Returns:
        float: Integral value, unitless
    """
    func = lambda u: ((2.+u) * np.sqrt(1 - u)) / (3. * u * u) * (1 - j0(1.5 * u * x))
    integral = integrate.quad(func, 0, 1)[0]
    return integral

def calculate_random_cylinder_integral_fast(x_arr):
    """Use multiprocessing to calculate the random cylinder integral quickly for an array of x.
    This may be slower than simply iterating over x for low sizes of x_arr (< 1000 elements)

    Args:
        x_arr (ndarray, float): A 1 dimensional array containing values x for which the random cylinder integral should be calculated. 
        x=t/tc, time shift of the signal compared to the spin echo location divided by 
        the characteristic time tc. (sometimes expressed instead as t * delta omega)
    """
    if x_arr.size == 1:
        values = calculate_random_cylinder_integral(x_arr)
    else:
        with multiprocessing.Pool() as pool:
            values = pool.map(calculate_random_cylinder_integral, x_arr)
    values = np.array(values)

    return values

def calculate_static_dephasing_regime_signal_decay(t, dbv, oef):
    """Calculates signal trajectory around spin echo for static dephasing regime.

    Args:
        t (float): ndarray with 1 dimension or scalar, the time from the spin echo (in ms) where the signal is calculated
        dbv (float): deoxygenated blood volume fraction. Unitless
        oef (float): oxygen extraction fraction. Unitless

    Returns:
        ndarray with 1 dimension or scalar: the relavtive signal intensity expected at time t from spin echo. same shape as t.
    """
    tc = calculate_tc(oef=oef)
    t = np.array(t)
    return np.exp(-dbv*calculate_random_cylinder_integral_fast(t/tc))

def calculate_static_dephasing_regime_signal_decay_with_precalculated_integral(t, dbv, oef, precalc_x_arr, precalc_integral_arr):
    """Calculates signal trajectory around spin echo for static dephasing regime. Uses precaluclated values of the integral
    to save computation time.

    Args:
        t (float): ndarray with 1 dimension or scalar, the time from the spin echo (in ms) where the signal is calculated
        dbv (float): deoxygenated blood volume fraction. Unitless
        oef (float): oxygen extraction fraction. Unitless
        precalc_x_arr (ndarray): 1 dimensional. values of t/tc for which the static dephasing integral had been precalculated.
        precalc_integral_arr (ndarray): 1 dimensional. precalculated dephasing integral values. Must have same shape as precalc_x_arr

    Returns:
        ndarray with 1 dimension or scalar: the relavtive signal intensity expected at time t from spin echo. same shape as t.
    """
    tc = calculate_tc(oef=oef)
    t = np.array(t)
    x = t/tc
    if np.any(x < np.amin(precalc_x_arr)) or np.any(x > np.amax(precalc_x_arr)): # Check if value can be interpolated from precalculated range
        warnings.warn("The precalculated integral range does not cover the desired t values! Default to full calculation of integral")
        return calculate_static_dephasing_regime_signal_decay(t, dbv, oef)
    return np.exp(-dbv*np.interp(x, precalc_x_arr, precalc_integral_arr))

def generate_static_dephasing_regime_signal_decay_function(dbv, oef):
    f = lambda tau_arr: calculate_static_dephasing_regime_signal_decay(dbv, oef, tau_arr)
    return f

def main():
    import time

    #x_arr = np.linspace(-10, 10, num=10001, endpoint=True)
    tau_arr = np.linspace(-10, 10, num=10001, endpoint=True)

    starttime = time.time()
    #values_arr = [calculate_random_cylinder_integral(x) for x in x_arr]
    #values_arr = calculate_random_cylinder_integral_fast(x_arr)
    values_arr = calculate_static_dephasing_regime_signal_decay(tau_arr, 0.03, 0.4)
    #value = calculate_static_dephasing_regime_signal_decay(10, 0.03, 0.4)
    endtime = time.time()
    print(f"Calculation for {values_arr.size} points took {endtime-starttime} s")

if __name__ == "__main__":
    main()