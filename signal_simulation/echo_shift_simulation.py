import numpy as np
from utilities.constants_and_helpers import calculate_gamma_G_product

def calculate_tau_eff(tau, kspace_velocity, G_phase=0., phase_0=0):
    """Calculate the actual time (from the spin echo) when the data corresponding to the center of k-space is read out.

    Args:
        tau (float): nominal ASE shift tau in ms
        kspace_velocity (float): k-space velocity of EPI train in phase direction, in rad / (mm ms)
        G_phase (float, optional): Strength of magnetic field gradient occurring at the investigated location (in micro T/m). Defaults to 0.
        phase_0 (float, optional): Phase of the signal at spin echo in rad/mm. Defaults to 0.

    Returns:
        float: The actual time tau_eff after the spin echo when the center of k-space is written
    """
    return (tau - phase_0/kspace_velocity)/(1 + calculate_gamma_G_product(G_phase)/kspace_velocity)

def calculate_tau_intercept_left(k_max, kspace_velocity, G_phase=0., phase_0=0):
    return -(-k_max*(1 + calculate_gamma_G_product(G_phase) / kspace_velocity) + phase_0)/calculate_gamma_G_product(G_phase)

def calculate_tau_intercept_right(k_max, kspace_velocity, G_phase=0., phase_0=0):
    return -(k_max*(1 + calculate_gamma_G_product(G_phase) / kspace_velocity) + phase_0)/calculate_gamma_G_product(G_phase)

def calculate_tau_min(k_max, kspace_velocity, G_phase=0., phase_0=0):
    if G_phase == 0:
        return -np.inf
    if G_phase < 0:
        return calculate_tau_intercept_left(k_max, kspace_velocity, G_phase, phase_0)
    return calculate_tau_intercept_right(k_max, kspace_velocity, G_phase, phase_0)

def calculate_tau_max(k_max, kspace_velocity, G_phase=0., phase_0=0):
    if G_phase == 0:
        return np.inf
    if G_phase < 0:
        return calculate_tau_intercept_right(k_max, kspace_velocity, G_phase, phase_0)
    return calculate_tau_intercept_left(k_max, kspace_velocity, G_phase, phase_0)