import numpy as np
import os
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Calculate minimum distances from configurations")
parser.add_argument(
    "--directory", "-d",
    type=str,
    required=True,
    help="The required home directory prefix where 'inputs' and 'outputs' will be created"
)
args = parser.parse_args()

Path = os.path.join(args.directory, 'outputs')

iterations_prefix = 'iteration_'

files = [f for f in os.listdir(Path) if f.startswith(iterations_prefix)]

fig, axs = plt.subplots(1, 3, figsize=(18, 6))  # Create a single figure with 3 horizontal subplots

files = sorted(files, key=lambda x: float(x.split('_')[1].split('.')[0]) if x.split('_')[1].split('.')[0].lstrip('-').isdigit() else float('inf'))
print("Sorted files:", files)
for i, file in enumerate(files):
    radii, gr, pot, forces = np.loadtxt(os.path.join(Path, file), unpack=True)

    axs[0].plot(radii, gr, label=f'{i-1}')
    axs[1].plot(radii, pot, label=f'{i-1}')

# Configure the first subplot for g(r)
axs[0].set_title('g(r)')
axs[0].set_xlabel(r'$r/\sigma$')
axs[0].set_ylabel('g(r)')
axs[0].legend()

print("Plotted g(r) and Potential for all iterations.")

# Configure the second subplot for Potential
axs[1].set_title('Potential')
axs[1].set_xlabel(r'$r/\sigma$')
axs[1].set_ylabel(r'$\beta u(r)$')
axs[1].legend()

# Load convergence data
iteration, _, err, iteration_diff, _, _, _ = np.loadtxt(os.path.join(Path, 'convergence_data.dat'), unpack=True)

# Configure the third subplot for Convergence
axs[2].semilogy(iteration, err, label='Error')
axs[2].semilogy(iteration, iteration_diff, label='Iteration Diff')
axs[2].set_title('Convergence')
axs[2].set_xlabel('Iteration')
axs[2].set_ylabel('Log Scale')
axs[2].legend()

plt.tight_layout()
plt.savefig(os.path.join(Path, '00Convergence_plots.pdf'), dpi=300)  
plt.show()
