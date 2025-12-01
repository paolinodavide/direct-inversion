import os
import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
import glob
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
from numba import njit
import sys

def find_lj_config_files(directory='./configs/', pattern='*.dat'):
    """Find all LJ configuration files in the directory"""
    files = glob.glob(os.path.join(directory, pattern))
    # Extract numbers from filenames and sort numerically
    files.sort(key=lambda x: int(x.split('_')[-1].split('.')[0]))
    return files

def rdf_from_file(filename, dr, r_max=10):
    """Process a single file to calculate RDF"""
    try:
        with open(filename, 'r') as f:
            # Read first line to get box size
            first_line = f.readline().strip()
            box_size = float(first_line.split('l_box = ')[1])

            # Read remaining lines for coordinates
            data = []
            for line in f:
                if line.strip():  # skip empty lines
                    x, y = map(float, line.strip().split())
                    data.append([x, y])
        positions = np.array(data)

        # Calculate RDF
        r, g_r = calculate_rdf_numba(positions, box_size, dr, r_max)

        # Save results
        file_num = os.path.basename(filename).split('_')[1].split('.')[0]
        output_folder = 'outputs/rdfs/'
        os.makedirs(output_folder, exist_ok=True)
        
        output_filename = f'{output_folder}/g_r_lj_{file_num}.dat'
        np.savetxt(output_filename, np.column_stack((r, g_r)), 
                  header='# r g(r)')

        return g_r, r

    except Exception as e:
        print(f"Error processing {filename}: {str(e)}")
        return None

def calculate_rdf_scipy(positions, box_size, dr=0.02, r_max=10, dim=2):
    """Optimized RDF calculation"""
    N = len(positions)
    rho = N / (box_size ** 2)
    if r_max > box_size / 2:
        r_max = box_size / 2
    bins = np.arange(0, r_max + dr, dr)
    r = bins[:-1] + dr/2
    g_r = np.zeros_like(r)

    # Vectorized distance calculation with PBC
    dist_nd_sq = np.zeros(N * (N - 1) // 2)
    for d in range(dim):
        pos_1d = positions[:, d][:, np.newaxis]
        dist_1d = pdist(pos_1d)
        dist_1d[dist_1d > box_size * 0.5] -= box_size
        dist_nd_sq += dist_1d ** 2
    
    dist_nd = np.sqrt(dist_nd_sq)
    valid_dist = dist_nd[dist_nd < r_max]
    
    # Vectorized bin counting
    bin_indices = (valid_dist // dr).astype(int)
    unique_bins, counts = np.unique(bin_indices, return_counts=True)
    g_r[unique_bins] = counts * 2  # count i-j and j-i

    # Normalization
    ring_areas = np.pi * (bins[1:]**2 - bins[:-1]**2)
    g_r /= (N * rho * ring_areas)

    return r, g_r

@njit
def calculate_rdf_numba(positions, box_size, dr=0.02, r_max=10):
    """Numba-optimized RDF calculation"""
    N, dim = positions.shape
    rho = N / (box_size ** dim)
    
    if r_max > box_size / 2:
        r_max = box_size / 2
    
    bins = np.arange(0, r_max + dr, dr)
    r = bins[:-1] + dr/2
    g_r = np.zeros_like(r)
    
    # Calculate all pairwise distances with PBC
    n_pairs = N * (N - 1) // 2
    dist_nd = np.zeros(n_pairs)
    
    for i in range(N):
        for j in range(i + 1, N):
            idx = i * (2 * N - i - 1) // 2 + j - i - 1
            dist_sq = 0.0
            for d in range(dim):
                dx = positions[i, d] - positions[j, d]
                dx = dx - box_size * np.round(dx / box_size)
                dist_sq += dx * dx
            dist_nd[idx] = np.sqrt(dist_sq)
    
    # Bin the distances
    valid_dist = dist_nd[dist_nd <= r_max + dr]
    bin_indices = (valid_dist / dr).astype(np.int64)
    
    # Count occurrences in each bin
    counts = np.zeros(len(r), dtype=np.int64)
    for idx in bin_indices:
        if idx < len(counts):
            counts[idx] += 1
    
    # Apply symmetry factor (count each pair only once)
    g_r = counts * 2  # multiply by 2 to count i-j and j-i
    
    # Normalization
    if dim == 2:
        ring_areas = np.pi * (bins[1:]**2 - bins[:-1]**2)
    else:  # dim == 3
        ring_areas = (4/3) * np.pi * (bins[1:]**3 - bins[:-1]**3)
    
    g_r = g_r / (N * rho * ring_areas)
    return r, g_r

def process_file(args):
    file, dr = args
    return rdf_from_file(file, dr) 

def main():
    if len(sys.argv) > 1:
        dr = float(sys.argv[1])
    else:
        print("Usage: python gr_histo_davide.py <dr> ")
        sys.exit(1)
    
    # Find all files to process
    inputs_path = "./inputs/"
    output_path = "./outputs/"
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(output_path+'/rdfs/', exist_ok=True)
    files = find_lj_config_files(inputs_path+"configs/")
    if not files:
        print("No LJ config files found in ./configs/")
        return
    
    print(f"Found {len(files)} files to process")
    
    # Process in parallel with progress bar
    with Pool() as pool:
        results = list(tqdm(pool.imap(process_file, [(file, dr) for file in files]), 
                                  total=len(files),
                                  desc="Processing RDFs"))

    
    # Combine results (if needed)
    valid_results = [res for res in results if res is not None]
    if valid_results:
        g_total = np.mean([g for g, _ in valid_results], axis=0)
        var_g = np.std([g for g, _ in valid_results], axis=0)
        r = valid_results[0][1]
        
        # Save combined RDF

        np.savetxt(output_path+'/rdfs/g_r_h_avg.dat', 
                  np.column_stack((r, g_total, var_g)),
                  header='# r g(r) var_g(r)')
        print("All files processed successfully")
    else:
        print("No valid results were generated")
    


    plt.figure(figsize=(8, 6))
    plt.scatter(r, g_total, s=10, label='g(r)')
    plt.xlabel(r'$r/\sigma$', fontsize=14)
    plt.ylabel('g(r)', fontsize=14)
    plt.title('Radial Distribution Function', fontsize=16)
    plt.tight_layout()
    plt.savefig(output_path + '/rdfs/00gr_histo.pdf', dpi=300)
    plt.show()

if __name__ == '__main__':
    main()