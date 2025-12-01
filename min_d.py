import numpy as np
import os
import sys
from scipy.spatial.distance import pdist
from tqdm import tqdm
from multiprocessing import Pool
import glob
from numba import njit
import matplotlib.pyplot as plt

def find_config_files(directory='./configs/', pattern='*.dat'):
    """Find all LJ configuration files in the directory"""
    files = glob.glob(os.path.join(directory, pattern))
    # Extract numbers from filenames and sort numerically
    files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    return files


@njit
def minDistance_from_positions(positions, box_size, dim=2):
    N = len(positions)
    dist_nd_sq = np.zeros(N * (N - 1) // 2)  # to match the result of pdist
    idx = 0

    for i in range(N - 1):
        for j in range(i + 1, N):
            dist_sq = 0.0
            for d in range(dim):
                dist_1d = positions[i, d] - positions[j, d]
                dist_1d -= box_size * np.round(dist_1d / box_size)  # Apply periodic boundary conditions
                dist_sq += dist_1d ** 2
            dist_nd_sq[idx] = dist_sq
            idx += 1

    dist_nd = np.sqrt(dist_nd_sq)
    return np.min(dist_nd)

def minDistance_from_file(file_path, dim=2):
    box_size = np.loadtxt(file_path, comments=None, max_rows=1, usecols=-1)
    positions = np.loadtxt(file_path, skiprows=1)

    min_distance = minDistance_from_positions(positions, box_size, dim)
    file_number = int(file_path.split('_')[-1].split('.')[0])
    return file_number, min_distance

def main():
    inputs_path = "./inputs/"
    output_path = "./outputs/"
    os.makedirs(output_path, exist_ok=True)
    files = find_config_files(inputs_path+"configs/")
    if not files:
        print("No config files found in ./configs/")
        return

    print(f"Found {len(files)} files to process")

    # Process in parallel with progress bar
    with Pool() as pool:
        results = list(tqdm(pool.imap(minDistance_from_file, files), 
                        total=len(files),
                        desc="Processing RDFs"))
    
    # Extract file numbers and minimum distances
    file_numbers, min_distances = zip(*results)

    plt.figure(figsize=(10, 6))
    plt.subplot(1, 2, 1)
    plt.scatter(file_numbers, min_distances, marker='o')
    plt.xlabel('File Number')
    plt.ylabel(r'$d_{\text{min}} / \sigma$')
    plt.title(f'Minimum Distance = {np.min(min_distances):.4f}' + r' $\sigma$')
    plt.grid(True)

    # Set y-axis limits
    dr = 0.002
    y_min, y_max = min(min_distances) - dr, max(min_distances) +dr
    plt.ylim(y_min, y_max)

    # Create histogram
    plt.subplot(1, 2, 2)
    bins = np.arange(y_min, y_max + 0.002, 0.002)
    plt.hist(min_distances, bins=bins, alpha=0.7, orientation='horizontal')
    plt.ylabel(r'$d_{\text{min}} / \sigma$')
    plt.xlabel('Frequency')
    plt.title('Histogram of Minimum Distances')
    plt.grid(True)

    # Set the same y-axis limits
    plt.ylim(y_min, y_max)

    # Save and show
    plt.tight_layout()
    plt.savefig(output_path+"min_distances.pdf", dpi=300)
    plt.show()

    sorted_file_numbers = [file_number for file_number, _ in sorted(results, key=lambda x: x[1])]
    np.savetxt(inputs_path+"ordered_wt.dat", sorted_file_numbers, 
           header='wt ordered by min dist', fmt='%d')
        
    return

if __name__ == "__main__":
    main()