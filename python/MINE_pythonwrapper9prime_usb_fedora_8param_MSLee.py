
from classy import Class
import numpy as np
import itertools


# Define parameter ranges

Omega_Lambda2 = np.arange(0.00000016, 0.00000018, 0.00000001)
Omega_Lambda3 = np.arange(-0.015, -0.013, 0.001)
h_range = np.arange(0.65, 0.74, 0.01)
omega_b_range = np.arange(0.0218, 0.0224, 0.0002)
omega_cdm_range = np.arange(0.118, 0.122, 0.002)
A_s_range = np.arange(2.09e-09, 2.13e-09, 0.01e-09)
n_s_range = np.arange(0.956, 0.968, 0.002)
tau_reio_range = np.arange(0.054, 0.057, 0.001)

# Compute Cartesian product of parameter ranges
param_sets = itertools.product(Omega_Lambda2, Omega_Lambda3, h_range, omega_b_range, omega_cdm_range,
                                A_s_range, n_s_range, tau_reio_range)

# Load data from the second file
second_file_data = np.loadtxt("./COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01_headerremoved.txt")
second_file_ttres = second_file_data[:, 1]  # Extract second column from the second file

# Write distances and parameter sets to file
output_file = "./resultoutput/distances.txt"

with open(output_file, 'w') as file:
    file.write("Distance, Omega_Lambda2, Omega_Lambda3, h, omega_b, omega_cdm, A_s, n_s, tau_reio\n")

    param_sets, param_sets_copy = itertools.tee(param_sets)
    total_loops = sum(1 for _ in param_sets_copy)
    print(f"Total loops: {total_loops}")
    
    # Iterate over parameter sets
    for i, param_set in enumerate(param_sets):
        # Create an instance of the CLASS wrapper
        cosmo = Class()

        # Set the parameters to the cosmological code
        
        cosmo.set({
            'output': 'tCl lCl',
            'l_max_scalars': 2508,
            'lensing': 'yes',            
            'Omega_Lambda2': param_set[0],
            'Omega_Lambda3': param_set[1],
            'h': param_set[2],
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

        # Compute Euclidean distance
        euclidean_distance = np.linalg.norm(ttres - second_file_ttres)

        # Write distance and parameter set to file
        file.write(f"{euclidean_distance}, {param_set[0]}, {param_set[1]}, {param_set[2]}, {param_set[3]}, {param_set[4]}, {param_set[5]}, {param_set[6]}, {param_set[7]}\n")
        if i % 100 == 0:
            print(f"[{i}] Distance: {euclidean_distance}")

        # Clean CLASS
        cosmo.struct_cleanup()
        cosmo.empty()

        # Delete CLASS instance to clear memory
        del cosmo

        # Clear distance_param_sets list to release memory
        distance_param_sets = []




# Read distances from the file
with open(output_file, 'r') as file:
    lines = file.readlines()

# Remove header line
lines = lines[1:]

# Sort lines based on distance
lines.sort(key=lambda x: float(x.split(',')[0]))

# Select the first 100 rows
top_100_rows = lines[:100]

# Write the selected rows to first100min_distances.txt
#output_top_100_file = "/home/namu_2nd_samsung/Documents/ONGOING/class_public/python/resultoutput/first100min_distances.txt"
output_top_100_file = "./resultoutput/first100min_distances.txt"
with open(output_top_100_file, 'w') as file:
    # Write header line
    file.write("Distance, Omega_Lambda2, Omega_Lambda3, h, omega_b, omega_cdm, A_s, n_s, tau_reio\n")
    # Write selected rows
    for line in top_100_rows:
        file.write(line)

print("Done")
