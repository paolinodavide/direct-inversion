import numpy as np
import os
import sys
from scipy.spatial.distance import pdist
from tqdm import tqdm
from multiprocessing import Pool
import glob
import matplotlib.pyplot as plt

def find_lj_config_files(directory='./configs/', pattern='lj_*.dat'):
    """Find all LJ configuration files in the directory"""
    files = glob.glob(os.path.join(directory, pattern))
    # Extract numbers from filenames and sort numerically
    files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    return files

def calculate_fileNb_minDistance(file_path, dim=2) :
    box_size = np.loadtxt(file_path, comments=None, max_rows=1, usecols=-1)
    positions = np.loadtxt(file_path, skiprows=1)

    N = len(positions)
    
    dist_nd_sq = np.zeros(N * (N - 1) // 2)  # to match the result of pdist
    for d in range(dim):
        pos_1d = positions[:, d][:, np.newaxis]  # shape (N, 1)
        dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
        dist_1d[dist_1d > box_size * 0.5] -= box_size
        dist_nd_sq += dist_1d ** 2  # d^2 = dx^2 + dy^2 + dz^2
    dist_nd = np.sqrt(dist_nd_sq)

    file_number = int(file_path.split('_')[-1].split('.')[0])
    return file_number, np.min(dist_nd)


def main():
    files = find_lj_config_files()
    if not files:
        print("No LJ config files found in ./configs/")
        return

    print(f"Found {len(files)} files to process")

    # Process in parallel with progress bar
    with Pool() as pool:
        results = list(tqdm(pool.imap(calculate_fileNb_minDistance, files), 
                        total=len(files),
                        desc="Processing RDFs"))
    
    # Extract file numbers and minimum distances
    file_numbers, min_distances = zip(*results)

    # Plot the results

    plt.figure(figsize=(10, 6))
    plt.plot(file_numbers, min_distances, marker='o', linestyle='-', color='b')
    plt.xlabel('File Number')
    plt.ylabel('Minimum Distance')
    plt.title('Minimum Distance vs File Number')
    plt.grid(True)
    plt.savefig('min_distance_plot.png')  # Save the plot as a PNG file
    plt.show()

    sorted_file_numbers = [file_number for file_number, _ in sorted(results, key=lambda x: x[1])]
    np.savetxt('ordered_wt.dat', sorted_file_numbers, 
           header='wt ordered by min dist', fmt='%d')
        
    return

if __name__ == "__main__":
    main()