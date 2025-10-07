import numpy as np
from scipy.optimize import minimize

# Lennard-Jones potential parameters for Argon
epsilon = 0.282  # eV
sigma = 3.843  # Å

# Define the Lennard-Jones potential function
def lj_potential(r, epsilon, sigma):
    r6 = (sigma / r)**6
    r12 = r6 * r6
    return epsilon * ((r12) - (2*r6))

# Generate FCC lattice positions
def generate_fcc_positions(lattice_constant, repetitions):
    positions = []
    basis = np.array([[0, 0, 0],
                      [0.5, 0.5, 0],
                      [0.5, 0, 0.5],
                      [0, 0.5, 0.5]])
    for x in range(repetitions):
        for y in range(repetitions):
            for z in range(repetitions):
                for b in basis:
                    pos = (np.array([x, y, z]) + b) * lattice_constant
                    positions.append(pos)
    return np.array(positions)

# Calculate the total potential energy per atom
def total_potential_energy_per_atom(lattice_constant, repetitions=3, cutoff=2.5):
    positions = generate_fcc_positions(lattice_constant, repetitions)
    num_atoms = len(positions)
    total_energy = 0.0
    cutoff_distance = cutoff * sigma

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            r_vec = positions[j] - positions[i]
            r_vec -= np.round(r_vec / (lattice_constant * repetitions)) * (lattice_constant * repetitions)
            r = np.linalg.norm(r_vec)
            if r < cutoff_distance:
                total_energy += lj_potential(r, epsilon, sigma)

    return total_energy / num_atoms

# Objective function for minimization
def objective(lattice_constant):
    return total_potential_energy_per_atom(lattice_constant)

# Initial guess for the lattice constant (in Å)
initial_guess = 5

# Perform the minimization
result = minimize(objective, initial_guess, method='Nelder-Mead')

# Optimal lattice constant
optimal_lattice_constant = result.x[0]
print(f'Optimal lattice constant: {optimal_lattice_constant:.3f} Å')