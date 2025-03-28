# Introduction

__ASE_MFG_simulation, Version 0.1.0__

This code simulates the signal intensity of asymmetric spin echo (ASE) Magnetic Resonance (MR) images with echo-planar
imaging (EPI) readout.

The effects of macroscopic magnetic field gradients (MFGs) are modelled as a shift in the time point where the center 
of k-space is read out. This results in a mix of R<sub>2</sub> and R<sub>2</sub>' decay of the signal, depending on the local MFGs. The figure
below shows simulated signal trajectories as a function of the offset $\tau$ for different MFGs in phase-encoding direction
G<sub>y</sub>.

![s_tau_phase_mfg_simulation_notitle_largefont](https://github.com/user-attachments/assets/6f409074-c418-4501-882e-d748b14f37b5)

Deriving quantitative Blood Oxygen Level Dependent contrast (qBOLD) measurements from ASE images is subject to bias if 
the effect of MFGs is not appropriately corrected.

The paper describing the theory behind this approach, and the results in a volunteer cohort, has been submitted for 
publication. This page will be updated with a link to it, if the publication is successful.

______________

# Credit:

This work was carried out in a collaboration between University Hospital Zurich and Erasmus Medical Center Rotterdam.

If you use this code, please cite the above mentioned paper.

In case of questions, contact Philipp Wallimann 
(philipp.wallimann@usz.ch).

______________

# Dependencies:

__For Simulations__

_Python libraries_: numpy, scipy, matplotlib

__For Image Processing__

_Python libraries_: numpy, scipy, matplotlib, nibabel

_Other_: FSL (the FMRIB Software Library)

______________

# Info:

Code clean-up is still ongoing.

Imaging data acquired as part of this work is not public, so all scripts modifying images are not directly functional.
The corresponding code can be adjusted to your own needs. This notably includes all Preprocessing shell scripts in the 
folder [preprocessing_shell_scripts](preprocessing_shell_scripts).

Simulations are functional based on this code alone. Results can be seen when running 
[signal_tau_trajectory_dependence_on_MFG.py](demonstration_scripts/signal_tau_trajectory_dependence_on_MFG.py),
[estimated_qBOLD_dependence_on_MFG.py](demonstration_scripts/estimated_qBOLD_dependence_on_MFG.py),
or 
[signal_tau_trajectory_interactive.py](demonstration_scripts/signal_tau_trajectory_interactive.py).
