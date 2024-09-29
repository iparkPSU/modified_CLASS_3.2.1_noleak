
from classy import Class
import numpy as np
import itertools


# Define parameter ranges
Omega_Lambda2 = np.arange(0.000016, 0.000018, 0.000001)
Omega_Lambda3 = np.arange(-0.015, -0.013, 0.001)
h_range = np.arange(0.644, 0.752, 0.002)
omega_b_range = np.arange(0.0218, 0.0224, 0.0002)
omega_cdm_range = np.arange(0.118, 0.122, 0.002)
A_s_range = np.arange(2.09e-09, 2.13e-09, 0.01e-09)
n_s_range = np.arange(0.956, 0.968, 0.002)
tau_reio_range = np.arange(0.054, 0.057, 0.0002)

# Compute Cartesian product of parameter ranges
param_sets = itertools.product(Omega_Lambda2, Omega_Lambda3, h_range, omega_b_range, omega_cdm_range,
                                A_s_range, n_s_range, tau_reio_range)




    
    # Iterate over parameter sets
    for i, param_set in enumerate(param_sets):
        # Create an instance of the CLASS wrapper
        cosmo = Class()

        # Set the parameters to the cosmological code
        cosmo.set({
            'output': 'tCl lCl',
            'l_max_scalars': 2508,
            'lensing': 'yes',
            'h': param_set[0],
            'Omega_Lambda2': param_set[1],
            'Omega_Lambda3': param_set[2],
            'omega_b': param_set[3],
            'omega_cdm': param_set[4],
            'A_s': param_set[5],
            'n_s': param_set[6],
            'tau_reio' : param_set[7],
        })

        # Run the code
        cosmo.compute()

        # Access the lensed Cls up to l=2508
        cls = cosmo.lensed_cl(2508)

        # Extract data for computing distance
        ell = cls['ell'][2:]  # Start from ell=2
        tt = cls['tt'][2:]  # Start from ell=2
        ttres = tt * 7.42563e12 * ell * (ell + 1) / (2 * np.pi)  # Compute ttres


        # Clean CLASS
        cosmo.struct_cleanup()
        cosmo.empty()

        # Delete CLASS instance to clear memory
        del cosmo

        # Clear distance_param_sets list to release memory
        distance_param_sets = []



