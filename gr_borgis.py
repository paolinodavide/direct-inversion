import os
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial
from numba import njit
from cProfile import Profile
import pstats

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


@njit(fastmath=True, cache=True) #Clear the cache if you want to recompile
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
                dx_ij -= box_length * np.round(dx_ij / box_length)
                dy_ij = particle_positions[i, 1] - particle_positions[j, 1]
                dy_ij -= box_length * np.round(dy_ij / box_length)

                d_ij_squared = dx_ij * dx_ij + dy_ij * dy_ij
                if d_ij_squared == 0.0 or d_ij_squared >= r_cut * r_cut:  # Skip self-interaction and cut-off
                    continue
                d_ij = np.sqrt(d_ij_squared)
                binIdx_ij = int(d_ij/ r_bin)

                if binIdx_ij > binlow:
                    alpha = d_ij / r_bin - binIdx_ij
                    force_magnitude = alpha * force_div_r[binIdx_ij - binlow] + (1 - alpha) * force_div_r[binIdx_ij - binlow + 1]
                elif binIdx_ij <= binlow:
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
            dx_ij -= box_length * np.round(dx_ij / box_length)
            dy_ij = particle_positions[i, 1] - particle_positions[j, 1]
            dy_ij -= box_length * np.round(dy_ij / box_length)

            d_ij_squared = dx_ij * dx_ij + dy_ij * dy_ij
            if d_ij_squared == 0.0:
                continue
            d_ij = np.sqrt(d_ij_squared)

            binIndex_ij = int((d_ij - min_radius) / r_bin)
            if binIndex_ij < 0 or binIndex_ij >= num_bins:
                continue
            # Delta_ij = ((Fi − Fj) ⋅ rij) / r²
            force_diff_x = net_forces[i, 0] - net_forces[j, 0]
            force_diff_y = net_forces[i, 1] - net_forces[j, 1]
            dot_product = force_diff_x * dx_ij + force_diff_y * dy_ij
            borgis_contributions[binIndex_ij] += dot_product / d_ij_squared

    # 3) Cumulative sum of contributions
    if method == 'in':
        borgis_contributions = np.cumsum(borgis_contributions)
    elif method == 'out':
        borgis_contributions = np.cumsum(borgis_contributions[::-1])[::-1]
    return borgis_contributions

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
    if max_config_nb > 0:
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
    # Add this after function definition
    #_ = grBorgis_notNorm(np.random.rand(10,2), 1.0, 0.1, 0.01, 100, np.zeros(100), 0.5, 2.5, 'out')

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

def main():
    # Example usage
    box_length = 60.0
    min_radius = 0.00
    bin_width = 0.002
    max_radius = 10.0
    x_cut = 2.5
    method = 'out'
    rlow = 0.93
    pot_length = int((x_cut - rlow) / bin_width)

    # Load the binned force data (example)
    potential = it.get_pot('lj_full', pot_length, bin_width, rlow, 1.0)
    force_div_r = it.get_x(potential, rlow, bin_width, pot_length)

    profiler = Profile()
    profiler.enable()
    
    # Compute gBorgis
    gr_result = gBorgis_parallel(box_length, force_div_r, rlow, max_radius, bin_width, x_cut, method=method)

    profiler.disable()
    stats = pstats.Stats(profiler)
    stats.sort_stats(pstats.SortKey.CUMULATIVE)
    with open("profiling.log", "w") as log_file:
        stats.stream = log_file
        stats.print_stats()
    
    

if __name__ == "__main__":
    main()
