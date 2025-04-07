import numpy as np
import matplotlib.pyplot as plt
import sys

try:
    filetype = sys.argv[1]
except IndexError:
    print("Enter the name of the json file to parse")

# List of specific filenames
files = [f"{filetype}_{i}.dat" for i in range(26)] + [f"{filetype}_target.dat"]

for file in files:
    try:
        # Load the data, skipping the first line (description)
        data = np.loadtxt(file, skiprows=1)

        # Assume the first column is positions and the second column is potentials
        positions = data[:, 0]
        potentials = data[:, 1]

        # Plot the data
        plt.plot(positions, potentials, label=str(file))
        plt.title(f"{str(filetype)} vs Positions")
        plt.ylabel(f"{str(filetype)}")

    except OSError:
        print(f"File {file} not found. Skipping.")

if filetype == 'lj':
    x = np.linspace(np.min(positions), np.max(positions), 500)

    lj_potential = 1/x**12 - 1/x**6
    lj_potential -= lj_potential[-1]  # Normalize to the first value
    plt.plot(x, lj_potential, '-x', label="LJ potential")
    
# Add labels, legend, and title
plt.xlabel("Position")
plt.legend()
plt.show()
