import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Creating a dataframe to store values
idx = ["Ac", "Ca", "Ce", "Es", "Ir", "Fe", "Rh", "Sr", "Th", "Yb"]

input_values = pd.DataFrame(index=idx)
input_values['12-6 LJ sigma'] = [3.843, 4.025, 3.734, 4.133, 2.785, 2.590, 2.757, 4.379, 3.683, 3.924]
input_values['12-6 LJ epsilon_0'] = [6.51, 3.36, 6.38, 2.88, 9.20, 6.00, 7.84, 3.40, 8.47, 2.71]
input_values['12-6 LJ A'] = (67550000, 60750000, 46870000, 71540000, 2003000, 546700, 1512000, 169000000, 52760000, 38160000)
input_values['12-6 LJ B'] = (41940, 28570, 34590, 28710, 8586, 3620, 6886, 47950, 42280, 20340)

input_values['9-6 LJ sigma'] = [3.907, 4.088, 3.797, 4.195, 2.836, 2.645, 2.807, 4.445, 3.746, 4.001]
input_values['9-6 LJ epsilon_0'] = [5.40, 2.80, 5.28, 2.39, 7.48, 4.79, 6.38, 2.85, 7.01, 2.24]

output_values = pd.DataFrame(index=idx)

output_values['Experimental'] = [5.311, 5.5884, 5.1610, 5.7500, 3.6394, 3.8392, 3.8032, 6.0840, 5.0842, 5.4848]

# Conversion factor from kcal/mol to eV
kcal_to_eV = 0.04336

# Define the Lennard-Jones potential function
def lj_potential_9_6(r1, epsilon1, sigma1):
    r6 = (sigma1 / r1) ** 6
    r9 = (sigma1 / r1) ** 9
    return epsilon1 * ((2 * r9) - (3 * r6))

# Generate FCC lattice positions
def generate_fcc_positions_9_6(lattice_constant_9_6, repetitions_9_6):
    positions_9_6 = []
    basis_9_6 = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    for x in range(repetitions_9_6):
        for y in range(repetitions_9_6):
            for z in range(repetitions_9_6):
                for b in basis_9_6:
                    pos_9_6 = (np.array([x, y, z]) + b) * lattice_constant_9_6
                    positions_9_6.append(pos_9_6)
    return np.array(positions_9_6)

# Calculate the total potential energy per atom
def total_potential_energy_per_atom_9_6(lattice_constant_9_6, epsilon1, sigma1, repetitions_9_6=3, cutoff1=2.5):
    positions_9_6 = generate_fcc_positions_9_6(lattice_constant_9_6, repetitions_9_6)
    num_atoms_9_6 = len(positions_9_6)
    total_energy_9_6 = 0.0
    cutoff_distance_9_6 = cutoff1 * sigma1  # Use scalar sigma for simplicity

    for i in range(num_atoms_9_6):
        for j in range(i + 1, num_atoms_9_6):
            r_vec_9_6 = positions_9_6[j] - positions_9_6[i]
            r_vec_9_6 -= np.round(r_vec_9_6 / (lattice_constant_9_6 * repetitions_9_6)) * (lattice_constant_9_6 * repetitions_9_6)
            r1 = np.linalg.norm(r_vec_9_6)
            if r1 < cutoff_distance_9_6:
                total_energy_9_6 += lj_potential_9_6(r1, epsilon1, sigma1)

    return total_energy_9_6 / num_atoms_9_6

# Optimize the lattice constant for each metal
for metal_9_6 in idx:
    epsilon1 = input_values.loc[metal_9_6, '9-6 LJ epsilon_0'] * kcal_to_eV
    sigma1 = input_values.loc[metal_9_6, '9-6 LJ sigma']

    def objective_9_6(lattice_constant_9_6):
        return total_potential_energy_per_atom_9_6(lattice_constant_9_6, epsilon1, sigma1)

    # Initial guess for the lattice constant (in Å)
    initial_guess_9_6 = 5.0
    result_9_6 = minimize(objective_9_6, initial_guess_9_6, method='Nelder-Mead')

    # Store the optimal lattice constant
    optimal_lattice_constant_9_6 = result_9_6.x[0]
    output_values.loc[metal_9_6, '9-6'] = optimal_lattice_constant_9_6
    print(f'{metal_9_6}: Optimal lattice constant: {optimal_lattice_constant_9_6:.3f} Å')

    # Display the initial table of values
print(input_values)
# Display the final output
print(5*output_values)