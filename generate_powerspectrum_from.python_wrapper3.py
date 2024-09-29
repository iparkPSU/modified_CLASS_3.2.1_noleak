from classy import Class
import numpy as np
import itertools
import os

# Define parameter ranges
Omega_Lambda2 = np.arange(0.00000017, 0.00000018, 0.0000001)
Omega_Lambda3 = np.arange(-0.013, -0.012, 0.001)
h_range = np.arange(0.67810, 0.752, 0.1)
omega_b_range = np.arange(0.02238280, 0.0224, 0.01)
omega_cdm_range = np.arange(0.1201075, 0.122, 0.01)
A_s_range = np.arange(2.100549e-09, 2.13e-09, 0.1e-09)
n_s_range = np.arange(0.9660499, 0.968, 0.01)
tau_reio_range = np.arange(0.05430842, 0.057, 0.01)

# Directory to save output files
output_directory = '/home/namu_2nd_samsung/Documents/ONGOING/class_public/output'

# Ensure the output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

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
        'Omega_Lambda2': param_set[0],  # Corrected indices
        'Omega_Lambda3': param_set[1],
        'h': param_set[2],
        'omega_b': param_set[3],
        'omega_cdm': param_set[4],
        'A_s': param_set[5],
        'n_s': param_set[6],
        'tau_reio': param_set[7],
    })

    # Run the code
    cosmo.compute()

    # Access the lensed Cls up to l=2508
    cls = cosmo.lensed_cl(2508)

    # Extract data for computing distance
    ell = cls['ell'][2:]  # Start from ell=2
    tt = cls['tt'][2:]    # Start from ell=2
    ttres = tt * ell * (ell + 1) / (2 * np.pi)  # Compute ttres

    # Generate file name based on parameter values
    file_name = f"{param_set[0]:.6f}_{param_set[1]:.6f}_{param_set[2]:.6f}_{param_set[3]:.6f}_{param_set[4]:.6f}_{param_set[5]:.10e}_{param_set[6]:.6f}_{param_set[7]:.6f}.txt"
    file_path = os.path.join(output_directory, file_name)

    # Open the file and write the parameter values in the first line
    with open(file_path, 'w') as f:
        # Write header: the 8 parameters in the first line
        f.write(f"{param_set[0]} {param_set[1]} {param_set[2]} {param_set[3]} {param_set[4]} {param_set[5]} {param_set[6]} {param_set[7]}\n")
        
        # Write data: ell and ttres in two columns
        for l, tr in zip(ell, ttres):
            f.write(f"{l} {tr}\n")

    # Clean CLASS
    cosmo.struct_cleanup()
    cosmo.empty()

    # Delete CLASS instance to clear memory
    del cosmo

    # Optional: Print progress every 100 parameter sets
    if i % 100 == 0:
        print(f"Processed {i} parameter sets")

print("All parameter sets processed and results saved.")




