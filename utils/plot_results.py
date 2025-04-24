import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

try:
    filetype = sys.argv[1]
except IndexError:
    print("Enter the name of the file to parse")


isForce = False
column_index = 1
plot_type = filetype
if filetype == 'forces':
    filetype = 'pot'
    column_index = 2

if filetype != 'err':
    # Automatically find all files matching the pattern
    files = glob.glob(f"{filetype}_*.dat")
    for file in files:
        try:
            # Load the data, skipping the first line (description)
            data = np.loadtxt(file, skiprows=1)

            # Assume the first column is positions and the second column is potentials
            positions = data[:, 0]
            y_data = data[:, column_index]

            # Plot the data
            plt.plot(positions, y_data, label=str(file))
            plt.title(f"{str(plot_type)} vs Positions")
            plt.ylabel(f"{str(filetype)}")


        except OSError:
            print(f"File {file} not found. Skipping.")


if filetype == 'err':
    try:
        # Load the data, skipping the first line (description)
        data = np.loadtxt(f"{filetype}.dat", skiprows=1)

        # Assume the first column is the data to plot
        values = data[:, 1]

        # Plot the data against the line index
        plt.scatter(range(len(values)), values, label=f"{filetype}.dat")
        plt.title(f"{str(filetype)} vs Line Index")
        plt.ylabel(f"{str(filetype)}")

    except OSError:
        print(f"File {filetype}.dat not found. Skipping.")

if filetype == 'lj':
    x = np.linspace(np.min(positions), np.max(positions), 500)

    lj_potential = 1/x**12 - 1/x**6
    lj_potential -= lj_potential[-1]  # Normalize to the first value
    plt.plot(x, lj_potential, '-x', label="LJ potential")
    
# Add labels, legend, and title
plt.xlabel("Position")
plt.grid(linestyle='--')
plt.savefig(f"{filetype}.svg", format="svg")
plt.legend()
plt.show()
