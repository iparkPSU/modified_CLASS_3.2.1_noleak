import numpy as np
import random
from classy import Class
import matplotlib.pyplot as plt

# Define parameter ranges
h_range = np.arange(62, 75.5, 1)
omega_b_range = np.arange(0.02200, 0.02280, 0.0001)
omega_cdm_range = np.arange(0.1100, 0.1440, 0.01)

# Generate random sets of parameters
random_params = random.sample(list(zip(h_range, omega_b_range, omega_cdm_range)), 3)

# Create a single plot for all parameter sets
plt.figure(figsize=(10, 6))

# Iterate over random parameter sets
for i, params in enumerate(random_params):
    # Create an instance of the CLASS wrapper
    cosmo = Class()

    # Set the parameters to the cosmological code
    cosmo.set({
        'output': 'tCl lCl',
        'l_max_scalars': 2000,
        'lensing': 'yes',
        'A_s': 2.3e-9,
        'n_s': 0.9624,
        'h': params[0],
        'omega_b': params[1],
        'omega_cdm': params[2]
    })

    # Run the whole code
    cosmo.compute()

    # Access the lensed Cls up to l=2000
    cls = cosmo.lensed_cl(2000)

    # Extract data for plotting
    ell = cls['ell'][10:]  # Start from ell=2
    tt = cls['tt'][10:]  # Start from ell=2

    # Plot the results with a different color for each parameter set
    plt.plot(ell, tt * ell * (ell + 1) / (2 * 3.14159), label=f'Set {i+1}')

    # Clean CLASS
    cosmo.struct_cleanup()
    cosmo.empty()

# Add title, labels, legend, and grid to the plot
plt.title('CMB Power Spectrum for Random Parameter Sets')
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$\ell(\ell+1)C_{\ell}/(2\pi)$ [$\mu K^2$]')
plt.legend()
plt.grid(True)
plt.show()


#I think I shoudl consider directly (wihout worrying about putting the results into a database) computing the distances between these results and Planck one.