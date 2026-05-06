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

iterations_prefix = 'gr_'

files = [f for f in os.listdir(Path) if f.startswith(iterations_prefix)]

fig, axs = plt.subplots(1, 1, figsize=(6, 6))  # Create a single figure with 3 horizontal subplots

#files = sorted(files, key=lambda x: float(x.split('_')[1].split('.')[0]) if x.split('_')[1].split('.')[0].lstrip('-').isdigit() else float('inf'))
print("files:", files)
for i, file in enumerate(files):
    radii, gr= np.loadtxt(os.path.join(Path, file), unpack=True)

    axs.plot(radii, gr, label=f'{i}')

# Configure the first subplot for g(r)
axs.set_title('g(r)')
axs.set_xlabel(r'$r/\sigma$')
axs.set_ylabel('g(r)')
axs .legend()

print("Plotted g(r) and Potential for all iterations.")
plt.tight_layout()
plt.show()