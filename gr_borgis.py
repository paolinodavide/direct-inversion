import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
from numba import njit
from tqdm import tqdm
from cProfile import Profile
import time
import pstats
import shutil

import gr_iteration as it

global NUM_THREADS
NUM_THREADS = 12

def initialize_potential(pot_length, r_bin, x_low):
    if True:
        u = it.get_pot('lj_full', pot_length, r_bin, x_low, 1.0)
        x = it.get_x(u, x_low, r_bin, pot_length)
    return u, x


def gBorgis_numpy(file, box_size, binned_force, bins):
    """ This function computes the Borgis g(r) using numpy for a given configuration file.
    It's suboptimal for large systems, but it's easy to follow and understand."""

    config = np.load(file)
    num_particles, dimensions = config.shape

    # Compute pairwise displacements with Minimum Image Convention (MIC)
    displacements = config[:, np.newaxis, :] - config[np.newaxis, :, :]
    displacements -= box_size * np.round(displacements / box_size)

    # Extract upper triangular indices for unique pairs
    i_indices, j_indices = np.triu_indices(num_particles, k=1)
    pairwise_vectors = displacements[i_indices, j_indices]
    
    squared_distances = np.sum(pairwise_vectors**2, axis=1)
    pairwise_distances = np.sqrt(squared_distances) # Memory allocation can be impoved

    # Bin indices for distances
    bin_indices = np.searchsorted(bins, pairwise_distances) - 1

    # Compute force magnitudes using precomputed binned forces
    force_magnitudes = binned_force[bin_indices]
    #np.nan_to_num(force_magnitudes, copy=False, nan=0.0, posinf=0.0, neginf=0.0)

    # Calculate pairwise forces
    pairwise_forces = (force_magnitudes[:, np.newaxis] * pairwise_vectors) / pairwise_distances[:, np.newaxis]

    # Accumulate total forces on each particle
    total_forces = np.zeros_like(config)
    np.add.at(total_forces, i_indices, pairwise_forces)
    np.add.at(total_forces, j_indices, -pairwise_forces)

    # Compute force differences for each pair
    force_differences = total_forces[i_indices] - total_forces[j_indices]

    # Calculate dot product and delta values
    dot_products = np.einsum("ij,ij->i", force_differences, pairwise_vectors)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        delta_values = np.divide(dot_products, squared_distances, out=np.zeros_like(dot_products), where=squared_distances > 0)

    # Bin and accumulate delta values
    valid_bins = (bin_indices >= 0) & (bin_indices < len(bins))
    gr_contributions = np.bincount(bin_indices[valid_bins], weights=delta_values[valid_bins], minlength=len(bins))


    return np.cumsum(gr_contributions)

@njit
def periodic_boundary_conditions(dx, box_length, pythonRound=True):
    """Fastest PBC implementation using floor division."""
    if pythonRound is True:
        return dx - box_length * np.round(dx / box_length)
    else:
        half_box = 0.5 * box_length
        if dx >= 0:
            return dx - box_length * np.floor((dx + half_box) / box_length)
        else:
            return dx - box_length * np.ceil((dx - half_box) / box_length)

@njit
def grBorgis_notNorm(particle_positions, box_length, min_radius, r_bin, num_bins, force_div_r, rlow, r_cut, method='out'):
    """Compute the Borgis g(r) using Numba for a given configuration file. The normalization is not included."""
    num_particles, dimensions = particle_positions.shape
    net_forces = np.zeros((num_particles, dimensions), dtype=np.float64)
    borgis_contributions = np.zeros(num_bins, dtype=np.float64)
    binlow = int(rlow / r_bin)

    # 1) Compute net forces
    for i in range(num_particles):
        for j in range(i+1, num_particles):
            if i != j:
                dx_ij = particle_positions[i, 0] - particle_positions[j, 0]
                dx_ij = periodic_boundary_conditions(dx_ij, box_length)
                dy_ij = particle_positions[i, 1] - particle_positions[j, 1]
                dy_ij = periodic_boundary_conditions(dy_ij, box_length)

                d_ij_squared = dx_ij * dx_ij + dy_ij * dy_ij
                if d_ij_squared == 0.0 or d_ij_squared > r_cut * r_cut:  # Skip self-interaction and cut-off
                    continue
                d_ij = np.sqrt(d_ij_squared)
                binIdx_ij = int(d_ij/ r_bin)

                if binIdx_ij > binlow: #Already below the cut-off
                    alpha = d_ij / r_bin - binIdx_ij
                    force_magnitude = alpha * force_div_r[binIdx_ij - binlow] + (1 - alpha) * force_div_r[binIdx_ij - binlow + 1]
                elif binIdx_ij <=  binlow:
                    alpha = binlow - binIdx_ij + 1
                    force_magnitude = force_div_r[0] + alpha * (force_div_r[0] - force_div_r[1])
                else:
                    force_magnitude = 0.0

                net_forces[i, 0] += dx_ij * force_magnitude
                net_forces[i, 1] += dy_ij * force_magnitude
                net_forces[j, 0] -= dx_ij * force_magnitude  # Newton's 3rd law
                net_forces[j, 1] -= dy_ij * force_magnitude
   
    # 2) Compute Borgis Delta contributions
    for i in range(num_particles):
        for j in range(i + 1, num_particles):
            dx_ij = particle_positions[i, 0] - particle_positions[j, 0]
            dx_ij = periodic_boundary_conditions(dx_ij, box_length)
            dy_ij = particle_positions[i, 1] - particle_positions[j, 1]
            dy_ij = periodic_boundary_conditions(dy_ij, box_length)
            d_ij_squared = dx_ij * dx_ij + dy_ij * dy_ij
            
            binIndex_ij = int( np.sqrt(d_ij_squared) / r_bin) 

            # Delta_ij = ((Fi − Fj) ⋅ rij) / r²
            force_diff_x = net_forces[i, 0] - net_forces[j, 0]
            force_diff_y = net_forces[i, 1] - net_forces[j, 1]
            Delta_ij = (force_diff_x * dx_ij + force_diff_y * dy_ij) / d_ij_squared
        
            if binIndex_ij >= num_bins:
                binIndex_ij = num_bins - 1 if method != 'in' else num_bins 
                #binIndex_ij = int(r_cut / r_bin) - 1 #Older and wrong version
            
            borgis_contributions[binIndex_ij] += Delta_ij 
    
    # 3) Cumulative sum of contributions
    if method == 'in':
        g0 = np.cumsum(borgis_contributions)
        return g0
    elif method == 'out':
        gInf = np.cumsum(borgis_contributions[::-1])[::-1]
        return gInf
    elif method == 'both':
        return borgis_contributions
    else:
        raise ValueError("Invalid method. Choose 'in', 'out', or 'both'.")

def cmpt_borgis_from_file(file_path, box_length, min_radius, bin_width, num_bins, force_bins, rlow, r_cut, method):

    particle_positions = np.load(file_path)
    gr = grBorgis_notNorm(particle_positions, box_length, min_radius, bin_width, num_bins, force_bins, rlow, r_cut, method)
    return gr


def gBorgis_parallel(box_length,force_div_r, x_low, max_radius, bin_width, x_cut, method='out', min_radius=0,
    configs_path="../configs_npy/lj", ordered_indices_file="../ordered_wt.dat", max_config_nb=-1):
    """    
    Compute the Borgis radial distribution function (g(r)) in parallel over multiple configuration files.
    Notice that the function is not normalized. The normalization should be done outside this function."""

    # Calculate the number of bins
    num_bins = int((max_radius - min_radius) / bin_width)

    # Convert force_div_r to a numpy array
    x_current = np.asarray(force_div_r)

    # Load ordered indices from file
    ordered_indices = np.loadtxt(ordered_indices_file, dtype=int)
    ordered_indices = ordered_indices[:max_config_nb]

    # Generate list of configuration files
    config_files = [
        f"{configs_path}_{i}.npy"
        for i in ordered_indices
        if os.path.exists(f"{configs_path}_{i}.npy")
    ]
    print(f"Evaluating Borgis g(r) for {len(config_files)} files...")

    # Warm up compilation with the first configuration file
    _ = cmpt_borgis_from_file(config_files[0], box_length, min_radius, bin_width, num_bins, x_current, x_low, x_cut, method)

    # Parallel computation using multiprocessing
    with Pool(NUM_THREADS) as pool:
        borgis_results = list(
            pool.imap(partial(cmpt_borgis_from_file,
                    box_length=box_length,
                    min_radius=min_radius,
                    bin_width=bin_width,
                    num_bins=num_bins,
                    force_bins=x_current,
                    rlow=x_low,
                    r_cut=x_cut,
                    method=method),
                config_files))

    # Compute and return the mean g(r)
    return np.mean(borgis_results, axis=0)


def gBorgis_parallel_variance(box_length,force_div_r, x_low, max_radius, bin_width, x_cut, prefactor,
    configs_path="../configs_npy/lj", ordered_indices_file="../ordered_wt.dat", max_config_nb=-1):
    """    
    Compute the Borgis radial distribution function (g(r)) in parallel over multiple configuration files.
    Notice that the function REQUIRES the normalization factor to be passed as an argument."""

    min_radius = 0.00
    # Calculate the number of bins
    num_bins = int((max_radius - min_radius) / bin_width)

    # Convert force_div_r to a numpy array
    x_current = np.asarray(force_div_r)

    # Load ordered indices from file
    ordered_indices = np.loadtxt(ordered_indices_file, dtype=int)
    ordered_indices = ordered_indices[:max_config_nb]

    # Generate list of configuration files
    config_files = [
        f"{configs_path}_{i}.npy"
        for i in ordered_indices
        if os.path.exists(f"{configs_path}_{i}.npy")]

    print(f"Evaluating Borgis g(r) for {len(config_files)} files...")

    # Warm up compilation with the first configuration file
    _ = cmpt_borgis_from_file(config_files[0], box_length, min_radius, bin_width, num_bins, x_current, x_low, x_cut, 'in')

    # Parallel computation using multiprocessing
    with Pool(NUM_THREADS) as pool:
        borgis_results = list(
            pool.imap(partial(cmpt_borgis_from_file,
                    box_length=box_length,
                    min_radius=min_radius,
                    bin_width=bin_width,
                    num_bins=num_bins,
                    force_bins=x_current,
                    rlow=x_low,
                    r_cut=x_cut,
                    method='both'),
                config_files))
    
    gr0 = []
    grInf = []
    for result in borgis_results:
        gr0.append(np.cumsum(result)/prefactor)
        grInf.append(1 - np.cumsum(result[::-1])[::-1]/prefactor)
    gr0 = np.array(gr0)
    grInf = np.array(grInf)


    Delta = grInf - gr0
    meanDelta = np.mean(Delta, axis=0)
    varDelta = np.var(Delta, axis=0)

    gr0_mean = np.mean(gr0, axis=0)
    cov0 = np.mean(gr0 * Delta, axis=0) - (meanDelta * gr0_mean)

    lambda0 =  - cov0 / varDelta
    grOpt = (1- lambda0) * gr0 + ( lambda0 )* grInf

    # plt.plot(np.var(grOpt, axis=0), label='Variance of g(r)')
    # plt.xlabel('Distance (r)')
    # plt.ylabel('Variance')
    # plt.title('Variance of g(r)')
    # plt.legend()
    # plt.grid()
    # plt.show()
    
    return np.mean(grOpt, axis=0)

def main():
    # Example usage
    box_length = 60.0
    min_radius = 0.00
    bin_width = 0.002
    max_radius = 10.0
    x_cut = 2.5
    method = 'out'
    x_low = 0.99
    pot_length = int((x_cut - x_low) / bin_width)

    ordered_indices_file = "../ordered_wt.dat"
    configs_path = "../configs_npy/lj"
    max_config_nb = -1


    # Load the binned force data (example)
    potential = it.get_pot('lj_full', pot_length, bin_width, x_low, 1.0)
    force_div_r = it.get_x(potential, x_low, bin_width, pot_length)


    # Convert force_div_r to a numpy array
    x_current = np.asarray(force_div_r)

    # Load ordered indices from file
    ordered_indices = np.loadtxt(ordered_indices_file, dtype=int)
    ordered_indices = ordered_indices[:max_config_nb]

    # Generate list of configuration files
    config_files = [
        f"{configs_path}_{i}.npy"
        for i in ordered_indices
        if os.path.exists(f"{configs_path}_{i}.npy")
    ]
    

    num_bins = int((max_radius - min_radius) / bin_width)
    
    config = np.load(config_files[0])
    N, dimensions = config.shape
    prefactor = 2 *np.pi * N**2 / box_length**2 /4

    print("Warm up compilation\n")
    _ = cmpt_borgis_from_file(config_files[0], box_length, min_radius, bin_width, num_bins, x_current, x_low, x_cut, 'in')


    start_time = time.time()
    gr_parallel = gBorgis_parallel(box_length, x_current, x_low, max_radius, bin_width, x_cut, method, min_radius=min_radius,
    configs_path=configs_path, ordered_indices_file=ordered_indices_file, max_config_nb=max_config_nb)

    if method == 'out':
        gr_parallel = 1 - gr_parallel/prefactor
    else:
        gr_parallel = gr_parallel/prefactor

    end_time = time.time()
    print(f"Time taken for parallel computation: {end_time - start_time:.2f} seconds")
    
    start_time = time.time()
    gr_parallel_opt = gBorgis_parallel_variance(box_length, x_current, x_low, max_radius, bin_width, x_cut, prefactor,
    configs_path=configs_path, ordered_indices_file=ordered_indices_file, max_config_nb=max_config_nb)
    end_time = time.time()
    print(f"Time taken for parallel computation with variance: {end_time - start_time:.2f} seconds")

    #Plot the results
    plt.figure(figsize=(8, 6))
    plt.plot(np.arange(min_radius, max_radius, bin_width), gr_parallel, label=f'g(r) from Borgis {method}', color='blue')
    plt.plot(np.arange(min_radius, max_radius, bin_width), gr_parallel_opt, label='g(r) from Borgis with variance', color='red')
    plt.xlabel('Distance (r)')
    plt.ylabel('g(r)')
    plt.title('Borgis g(r) with and without variance')
    plt.legend()
    plt.grid()
    plt.show()
    return

if __name__ == "__main__":
    main()
