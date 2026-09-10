using LoopVectorization
using StaticArrays
using Statistics
using ThreadsX
using JSON
using DelimitedFiles
using LinearAlgebra
using Plots

abstract type IntegrationMethod end
struct In <: IntegrationMethod end
struct Out <: IntegrationMethod end
struct Both <: IntegrationMethod end

"""
Precomputed pair data for a single particle pair.
Stores the minimum geometric information needed for both force and Borgis evaluation.
"""
struct PairData{D}
    i::Int32
    j::Int32
    rVec::SVector{D, Float64}
    r::Float64              # = sqrt(dot(rVec, rVec))
end

"""
Cached pair data for one configuration snapshot.
The pairs vector is partitioned: pairs[1:n_force_pairs] have r ≤ r_cut,
pairs[n_force_pairs+1:end] have r_cut < r ≤ max_dist.
"""
struct SnapshotCache{D}
    n_particles::Int
    pairs::Vector{PairData{D}}
    n_force_pairs::Int
end

"""
Apply minum image convention to a separation vector given the box length and its inverse.
"""
@inline function wrap_pbc_distances(separation::Float64, box_length::Float64, inv_box_length::Float64)::Float64
    return separation - box_length * round(separation * inv_box_length)
end

"""
    Compute the minimum image distance vectore between two positions under PBC.
"""
@inline function pbc_distance(pos_i::SVector{D, Float64}, pos_j::SVector{D, Float64}, box_sizes::SVector{D, Float64}) where D
    d = pos_i - pos_j

    # "Map" compiles to a single CPU instruction loop. 
    image = map((x, L) -> x - round(x / L) * L, d, box_sizes)

    #d -= image # PBCs
    return image, dot(image, image)
end

struct CellList{D}
    n_cells_1d::Int
    cell_width::Float64
    cell_starts::Vector{Int}
    particle_indices::Vector{Int}
    cell_neighbors::Vector{Vector{Int}}
end

function get_forward_neighbor_offsets(D::Int)
    offsets = SVector{D, Int}[]
    if D == 2
        for dx in -1:1
            for dy in -1:1
                if dx > 0 || (dx == 0 && dy > 0)
                    push!(offsets, SVector(dx, dy))
                end
            end
        end
    elseif D == 3
        for dx in -1:1
            for dy in -1:1
                for dz in -1:1
                    if dx > 0 || (dx == 0 && dy > 0) || (dx == 0 && dy == 0 && dz > 0)
                        push!(offsets, SVector(dx, dy, dz))
                    end
                end
            end
        end
    end
    return offsets
end

@inline function get_cell_index(coords::SVector{D, Int}, n_cells_1d::Int)::Int where D
    cid = 1
    mult = 1
    for d in 1:D
        # Wrap coordinates periodically
        wrapped_c = mod(coords[d], n_cells_1d)
        cid += wrapped_c * mult
        mult *= n_cells_1d
    end
    return cid
end

function cell_index_to_coords(cid::Int, n_cells_1d::Int, D::Int)
    coords = Vector{Int}(undef, D)
    temp = cid - 1
    for d in 1:D
        coords[d] = temp % n_cells_1d
        temp = div(temp, n_cells_1d)
    end
    return SVector{D, Int}(coords)
end

function compute_cell_neighbors(n_cells_1d::Int, D::Int)
    total_cells = n_cells_1d^D
    neighbor_offsets = get_forward_neighbor_offsets(D)
    
    cell_neighbors = [Int[] for _ in 1:total_cells]
    
    for cid in 1:total_cells
        coords = cell_index_to_coords(cid, n_cells_1d, D)
        
        for offset in neighbor_offsets
            neighbor_coords = coords + offset
            neighbor_cid = get_cell_index(neighbor_coords, n_cells_1d)
            push!(cell_neighbors[cid], neighbor_cid)
        end
    end
    return cell_neighbors
end

function build_cell_list(positions::Vector{SVector{D, Float64}}, box_length::Float64, r_cut::Float64) where D
    n_cells_1d = max(1, floor(Int, box_length / r_cut))
    cell_width = box_length / n_cells_1d
    total_cells = n_cells_1d^D
    
    cell_counts = zeros(Int, total_cells)
    
    for pos in positions
        c_coords = map(x -> clamp(floor(Int, mod(x, box_length) / cell_width), 0, n_cells_1d - 1), pos)
        cid = 1
        mult = 1
        for d in 1:D
            cid += c_coords[d] * mult
            mult *= n_cells_1d
        end
        cell_counts[cid] += 1
    end
    
    cell_starts = ones(Int, total_cells + 1)
    for cid in 1:total_cells
        cell_starts[cid + 1] = cell_starts[cid] + cell_counts[cid]
    end
    
    particle_indices = Vector{Int}(undef, length(positions))
    current_offsets = copy(cell_starts)
    
    for (i, pos) in enumerate(positions)
        c_coords = map(x -> clamp(floor(Int, mod(x, box_length) / cell_width), 0, n_cells_1d - 1), pos)
        cid = 1
        mult = 1
        for d in 1:D
            cid += c_coords[d] * mult
            mult *= n_cells_1d
        end
        
        offset = current_offsets[cid]
        particle_indices[offset] = i
        current_offsets[cid] += 1
    end
    
    cell_neighbors = compute_cell_neighbors(n_cells_1d, D)
    
    return CellList{D}(n_cells_1d, cell_width, cell_starts, particle_indices, cell_neighbors)
end


""" 
Compute g(r) using the Force prescription of Borgis et al.
Method specifyes the integration direction:
- "in": integrate from 0 outward
- "out": integrate from infty inward
- "both": return raw integrand
"""
function grForce_notNorm_svectorized(particle_positions::Matrix{Float64},
    box_length::Float64,
    bin_width::Float64,
    num_bins_gr::Int,
    force_over_r::Vector{Float64},
    r_min_interaction::Float64,
    r_cutoff_interaction::Float64,
    method::IntegrationMethod;
    core_strength::Int=13)
    
    num_particles, num_dimensions = size(particle_positions)
    if num_dimensions == 2
        return grForce_notNorm_svectorized_impl(Val(2), particle_positions, box_length, bin_width, num_bins_gr, force_over_r, r_min_interaction, r_cutoff_interaction, method; core_strength=core_strength)
    elseif num_dimensions == 3
        return grForce_notNorm_svectorized_impl(Val(3), particle_positions, box_length, bin_width, num_bins_gr, force_over_r, r_min_interaction, r_cutoff_interaction, method; core_strength=core_strength)
    else
        throw(ArgumentError("Unsupported physical dimension: $num_dimensions. Only 2D and 3D are supported."))
    end
end

function grForce_notNorm_svectorized_impl(::Val{dim},
    particle_positions::Matrix{Float64},
    box_length::Float64,
    bin_width::Float64,
    num_bins_gr::Int,
    force_over_r::Vector{Float64},
    r_min_interaction::Float64,
    r_cutoff_interaction::Float64,
    method::IntegrationMethod;
    core_strength::Int=13) where dim
    
    num_particles = size(particle_positions, 1)
    inv_bin_width = 1.0 / bin_width

    # Convert positions to static vectors of correct physical dimension
    positions = [SVector{dim, Float64}(particle_positions[i, :]) for i in 1:num_particles]
    total_forces = zeros(SVector{dim, Float64}, num_particles)

    # 1. Build cell list for force calculation (using r_cutoff_interaction as cutoff)
    cell_list_force = build_cell_list(positions, box_length, r_cutoff_interaction)

    # 2. Evaluate total forces using the cell list
    evaluate_total_forces!(total_forces, positions, cell_list_force, box_length, force_over_r, r_min_interaction, r_cutoff_interaction, bin_width; core_strength=core_strength)

    # 3. Build cell list for Borgis contributions (using max_distance as cutoff)
    max_dist = num_bins_gr * bin_width
    cell_list_borgis = build_cell_list(positions, box_length, max_dist)

    # 4. Compute Borgis contributions using the cell list (truncated integration at max_distance)
    borgis_contributions = zeros(Float64, num_bins_gr)
    compute_borgis_contributions!(borgis_contributions, positions, total_forces, cell_list_borgis, box_length, inv_bin_width, num_bins_gr)

    return integrate_borgis_contributions(borgis_contributions, method)
end

"""
Compute contributions B_ij in the Borgis g(r) formula, which involves the relative forces and positions of particle pairs.
"""
@inline function compute_borgis_contributions!(
    borgis_contributions::Vector{Float64}, 
    positions::Vector{SVector{D, Float64}},     
    total_forces::Vector{SVector{D, Float64}},
    cell_list::CellList{D},
    box_length::Float64, 
    inv_bin_width::Float64, 
    num_bins_gr::Int
    ) where D
    #fill!(borgis_contributions, 0.0)
    box_sizes = fill(box_length, SVector{D, Float64})
    total_cells = cell_list.n_cells_1d^D
    max_dist2 = (num_bins_gr * (1.0 / inv_bin_width))^2
    
    @inbounds for cid in 1:total_cells
        start_i = cell_list.cell_starts[cid]
        end_i = cell_list.cell_starts[cid + 1] - 1
        
        # Pairs within the same cell
        for idx_i in start_i:end_i
            i = cell_list.particle_indices[idx_i]
            pos_i = positions[i]
            force_i = total_forces[i]
            
            for idx_j in (idx_i + 1):end_i
                j = cell_list.particle_indices[idx_j]
                pos_j = positions[j]
                force_j = total_forces[j]
                
                rVec_ij, r2_ij = pbc_distance(pos_i, pos_j, box_sizes)
                if max_dist2 < r2_ij
                    continue
                end
                
                r_ij = sqrt(r2_ij)
                radial_bin_index = floor(Int, r_ij * inv_bin_width)
                target_bin = radial_bin_index + 1
                if target_bin > num_bins_gr
                    continue # Discard pairs beyond max_distance
                end
                
                force_diff = force_i - force_j
                borgis_delta = borgis_delta_calculation(force_diff, rVec_ij, r2_ij, r_ij)
                borgis_contributions[target_bin] += borgis_delta
            end
        end
        
        # Pairs between cell cid and neighboring cells
        for neighbor_cid in cell_list.cell_neighbors[cid]
            start_j = cell_list.cell_starts[neighbor_cid]
            end_j = cell_list.cell_starts[neighbor_cid + 1] - 1
            
            for idx_i in start_i:end_i
                i = cell_list.particle_indices[idx_i]
                pos_i = positions[i]
                force_i = total_forces[i]
                
                for idx_j in start_j:end_j
                    j = cell_list.particle_indices[idx_j]
                    pos_j = positions[j]
                    force_j = total_forces[j]
                    
                    rVec_ij, r2_ij = pbc_distance(pos_i, pos_j, box_sizes)
                    if max_dist2 < r2_ij
                        continue
                    end
                    
                    r_ij = sqrt(r2_ij)
                    radial_bin_index = floor(Int, r_ij * inv_bin_width)
                    target_bin = radial_bin_index + 1
                    if target_bin > num_bins_gr
                        continue # Discard pairs beyond max_distance
                    end
                    
                    force_diff = force_i - force_j
                    borgis_delta = borgis_delta_calculation(force_diff, rVec_ij, r2_ij, r_ij)
                    borgis_contributions[target_bin] += borgis_delta
                end
            end
        end
    end
    return nothing
end

@inline function integrate_borgis_contributions(contributions::Vector{Float64}, ::In)
    cumsum!(contributions, contributions)
end
@inline function integrate_borgis_contributions(contributions::Vector{Float64}, ::Out)
    reverse!(contributions)
    cumsum!(contributions, contributions)
    reverse!(contributions)
end
@inline function integrate_borgis_contributions(contributions::Vector{Float64}, ::Both)
    return contributions
end

# 2D Specialization
@inline function borgis_delta_calculation(force_diff::SVector{2, Float64}, rVec_ij::SVector{2, Float64}, r2_ij::Float64, r_ij::Float64)::Float64  
    return dot(force_diff, rVec_ij) / r2_ij
end
@inline function borgis_delta_calculation(force_diff::SVector{3, Float64}, rVec_ij::SVector{3, Float64}, r2_ij::Float64, r_ij::Float64)::Float64
    return dot(force_diff, rVec_ij) / (r2_ij * r_ij)
end
@inline function borgis_delta_calculation(force_diff::SVector{D, Float64}, rVec_ij::SVector{D, Float64}, r2_ij::Float64, r_ij::Float64)::Float64 where D
    return dot(force_diff, rVec_ij) / (r_ij^D)
end

# ══════════════════════════════════════════════════════════════════════════════
# In-memory pair caching: precompute geometry once, evaluate per iteration
# ══════════════════════════════════════════════════════════════════════════════

"""
    precompute_snapshot_pairs(Val(D), positions_matrix, box_length, bin_width, num_bins_gr, r_cut)

Build a partitioned pair list for one snapshot. Uses a single cell list at 
max_dist = num_bins_gr * bin_width (the larger cutoff). Pairs with r ≤ r_cut 
are stored first, followed by pairs with r_cut < r ≤ max_dist.
"""
function precompute_snapshot_pairs(::Val{D},
        positions_matrix::Matrix{Float64},
        box_length::Float64, bin_width::Float64,
        num_bins_gr::Int, r_cut::Float64) where D

    N = size(positions_matrix, 1)
    positions = [SVector{D,Float64}(positions_matrix[i,:]) for i in 1:N]
    max_dist = num_bins_gr * bin_width

    # One cell list at the larger cutoff
    cell_list = build_cell_list(positions, box_length, max_dist)
    box_sizes = fill(box_length, SVector{D, Float64})
    max_dist2 = max_dist^2
    r_cut2 = r_cut^2
    total_cells = cell_list.n_cells_1d^D

    force_pairs = PairData{D}[]
    far_pairs   = PairData{D}[]

    @inbounds for cid in 1:total_cells
        start_i = cell_list.cell_starts[cid]
        end_i   = cell_list.cell_starts[cid + 1] - 1

        # Same-cell pairs
        for idx_i in start_i:end_i
            i = cell_list.particle_indices[idx_i]
            for idx_j in (idx_i + 1):end_i
                j = cell_list.particle_indices[idx_j]
                rVec, r2 = pbc_distance(positions[i], positions[j], box_sizes)
                r2 > max_dist2 && continue
                r = sqrt(r2)
                p = PairData{D}(Int32(i), Int32(j), rVec, r)
                if r2 <= r_cut2
                    push!(force_pairs, p)
                else
                    push!(far_pairs, p)
                end
            end
        end

        # Neighbor-cell pairs
        for neighbor_cid in cell_list.cell_neighbors[cid]
            start_j = cell_list.cell_starts[neighbor_cid]
            end_j   = cell_list.cell_starts[neighbor_cid + 1] - 1
            for idx_i in start_i:end_i
                i = cell_list.particle_indices[idx_i]
                for idx_j in start_j:end_j
                    j = cell_list.particle_indices[idx_j]
                    rVec, r2 = pbc_distance(positions[i], positions[j], box_sizes)
                    r2 > max_dist2 && continue
                    r = sqrt(r2)
                    p = PairData{D}(Int32(i), Int32(j), rVec, r)
                    if r2 <= r_cut2
                        push!(force_pairs, p)
                    else
                        push!(far_pairs, p)
                    end
                end
            end
        end
    end

    n_force = length(force_pairs)
    pairs = vcat(force_pairs, far_pairs)

    return SnapshotCache{D}(N, pairs, n_force)
end

"""
    precompute_all_snapshots(directory, box_length, bin_width, num_bins_gr, r_cut)

Parallel precomputation of pair lists for all binary config files in a directory.
"""
function precompute_all_snapshots(
        directory::String, box_length::Float64, bin_width::Float64,
        num_bins_gr::Int, r_cut::Float64)

    file_paths = String[]
    for (root, _, files) in walkdir(directory)
        for file in files
            endswith(file, ".bin") && push!(file_paths, joinpath(root, file))
        end
    end
    isempty(file_paths) && throw(ArgumentError("No binary files found in directory: $directory"))

    # Infer dimension from first file
    sample = read_particle_positions_binary(file_paths[1])
    dim = size(sample, 2)

    caches = ThreadsX.map(file_paths) do fp
        pos = read_particle_positions_binary(fp)
        precompute_snapshot_pairs(Val(dim), pos, box_length, bin_width, num_bins_gr, r_cut)
    end
    return caches
end

"""
    evaluate_gr_from_cache(cache, force_over_r, num_bins_gr, r_min, r_cut, bin_width, method; core_strength)

Evaluate g(r) for a single snapshot using precomputed pair data.
Force loop iterates only over force pairs (1:n_force_pairs).
Borgis loop iterates over all pairs.
"""
function evaluate_gr_from_cache(
        cache::SnapshotCache{D},
        force_over_r::Vector{Float64},
        num_bins_gr::Int,
        r_min::Float64,
        r_cut::Float64,
        bin_width::Float64,
        method::IntegrationMethod;
        core_strength::Int=13) where D

    inv_bin_width = 1.0 / bin_width

    # ── Step 1: forces — only force pairs ──
    total_forces = zeros(SVector{D,Float64}, cache.n_particles)

    @inbounds for idx in 1:cache.n_force_pairs
        p = cache.pairs[idx]

        f = if p.r < r_min
            force_magnitude_below_rmin(p.r, r_min, force_over_r, inv_bin_width;
                                       core_strength=core_strength)
        else
            force_magnitude_between_bins(p.r, r_min, r_cut, force_over_r, inv_bin_width)
        end

        total_forces[p.i] += f * p.rVec
        total_forces[p.j] -= f * p.rVec
    end

    # ── Step 2: Borgis contributions — ALL pairs ──
    borgis = zeros(Float64, num_bins_gr)

    @inbounds for p in cache.pairs
        bin = floor(Int, p.r * inv_bin_width) + 1
        bin > num_bins_gr && continue

        Δf = total_forces[p.i] - total_forces[p.j]
        borgis[bin] += borgis_delta_calculation(Δf, p.rVec, p.r * p.r, p.r)
    end

    return integrate_borgis_contributions(borgis, method)
end

"""
    evaluate_gr_from_caches(caches, force_over_r, num_bins_gr, r_min, r_cut, bin_width, box_length, method; core_strength)

Parallel evaluation and averaging of g(r) across all cached snapshots.
Replaces `gr_force_from_dir_parallel_binary` inside the iteration loop.
"""
function evaluate_gr_from_caches(
        caches::Vector{<:SnapshotCache},
        force_over_r::Vector{Float64},
        num_bins_gr::Int,
        r_min::Float64, r_cut::Float64,
        bin_width::Float64,
        box_length::Float64,
        method::IntegrationMethod;
        core_strength::Int=13)

    results = ThreadsX.map(caches) do cache
        evaluate_gr_from_cache(cache, force_over_r, num_bins_gr,
                               r_min, r_cut, bin_width, method;
                               core_strength=core_strength)
    end

    if method != Both()
        file_count = length(results)
        avg = sum(results) ./ file_count
        var = sum(x -> x .^ 2, results) ./ file_count .- avg .^ 2
        return avg, var
    else
        # Both() mixing logic: blend inner and outer estimators
        N = caches[1].n_particles
        prefactor = compute_prefactor(N, box_length, length(caches[1].pairs[1].rVec))
        inv_prefactor = 1.0 / prefactor

        num_results = length(results)
        num_bins = length(results[1])
        grOpt_sum = zeros(Float64, num_bins)
        half_bins = div(num_bins, 2)
        temp_cumsum = Vector{Float64}(undef, num_bins)

        for result in results
            cumsum!(temp_cumsum, result)
            total_sum = temp_cumsum[end]
            for i in 1:num_bins
                val_gr0 = temp_cumsum[i] * inv_prefactor
                current_rev_cumsum = total_sum - (i > 1 ? temp_cumsum[i-1] : 0)
                val_grInf = 1.0 - (current_rev_cumsum * inv_prefactor)
                if i <= half_bins
                    grOpt_sum[i] += val_gr0
                else
                    grOpt_sum[i] += val_grInf
                end
            end
        end

        grOpt_estimator = (grOpt_sum ./ num_results) .* prefactor
        lambda0 = vcat(zeros(half_bins), ones(num_bins - half_bins))
        return grOpt_estimator, lambda0
    end
end

"""
Compute total force acting on each particle.
"""
@inline function evaluate_total_forces!(total_forces::Vector{SVector{D, Float64}},
    positions::Vector{SVector{D, Float64}},
    cell_list::CellList{D},
    box_length::Float64, force_over_r::Vector{Float64},
    r_min_interaction::Float64, 
    r_cutoff_interaction::Float64, 
    bin_width::Float64;
    core_strength::Int=0) where D

    num_particles = length(positions)
    box_sizes = fill(box_length, SVector{D, Float64})
    inv_bin_width = 1.0 / bin_width
    r_cut2 = r_cutoff_interaction^2
    total_cells = cell_list.n_cells_1d^D

    fill!(total_forces, zero(SVector{D, Float64}))

    @inbounds for cid in 1:total_cells
        start_i = cell_list.cell_starts[cid]
        end_i = cell_list.cell_starts[cid + 1] - 1
        
        # Pairs within the same cell
        for idx_i in start_i:end_i
            i = cell_list.particle_indices[idx_i]
            pos_i = positions[i]
            
            for idx_j in (idx_i + 1):end_i
                j = cell_list.particle_indices[idx_j]
                pos_j = positions[j]
                
                rVec_ij, r2_ij = pbc_distance(pos_i, pos_j, box_sizes)
                if r_cut2 < r2_ij
                    continue
                end
                
                r_ij = sqrt(r2_ij)
                f_magnitude = if r_ij < r_min_interaction
                    force_magnitude_below_rmin(r_ij, r_min_interaction, force_over_r, inv_bin_width; core_strength=core_strength)
                else
                    force_magnitude_between_bins(r_ij, r_min_interaction, r_cutoff_interaction, force_over_r, inv_bin_width)
                end
                
                total_forces[i] += f_magnitude * rVec_ij
                total_forces[j] -= f_magnitude * rVec_ij
            end
        end
        
        # Pairs between cell cid and its neighboring cells
        for neighbor_cid in cell_list.cell_neighbors[cid]
            start_j = cell_list.cell_starts[neighbor_cid]
            end_j = cell_list.cell_starts[neighbor_cid + 1] - 1
            
            for idx_i in start_i:end_i
                i = cell_list.particle_indices[idx_i]
                pos_i = positions[i]
                
                for idx_j in start_j:end_j
                    j = cell_list.particle_indices[idx_j]
                    pos_j = positions[j]
                    
                    rVec_ij, r2_ij = pbc_distance(pos_i, pos_j, box_sizes)
                    if r_cut2 < r2_ij
                        continue
                    end
                    
                    r_ij = sqrt(r2_ij)
                    f_magnitude = if r_ij < r_min_interaction
                        force_magnitude_below_rmin(r_ij, r_min_interaction, force_over_r, inv_bin_width; core_strength=core_strength)
                    else
                        force_magnitude_between_bins(r_ij, r_min_interaction, r_cutoff_interaction, force_over_r, inv_bin_width)
                    end
                    
                    total_forces[i] += f_magnitude * rVec_ij
                    total_forces[j] -= f_magnitude * rVec_ij
                end
            end
        end
    end
end

@inline function force_magnitude_below_rmin(r_ij::Float64, r_min::Float64, force_over_r::Vector{Float64}, inv_bin_width::Float64; core_strength::Int=0)::Float64
    if core_strength == 0
        return 0.0
    elseif core_strength == 1
        bin_width = 1.0 / inv_bin_width
        a = r_min / r_ij
        return a * force_over_r[1] + a*(1-a)*inv_bin_width * r_min* ((force_over_r[3] - force_over_r[2]))
    else
        return force_over_r[1] * (r_min/r_ij)^core_strength
    end
end

@inline function force_magnitude_between_bins(r_ij::Float64, r_min::Float64, r_max::Float64, force_over_r::Vector{Float64}, inv_bin_width::Float64)::Float64
    radial_bin_index = floor(Int, r_ij * inv_bin_width)
    force_table_index = clamp(radial_bin_index - floor(Int, r_min * inv_bin_width) + 1, 1, length(force_over_r)-1)
    interpolation_weight = r_ij * inv_bin_width - radial_bin_index
    return (1.0 - interpolation_weight) * force_over_r[force_table_index] + interpolation_weight * force_over_r[force_table_index+1]
end

function lennard_jones_force_div_r(r, epsilon=1.0, sigma=1.0)
    """
    Compute Lennard-Jones force divided by r: F(r)/r = 24ϵ/σ² * (2(σ/r)^14 - (σ/r)^8)
    """
    if r == 0.0
        return 0.0
    end
    sr = sigma / r
    sr8 = sr^8
    sr14 = sr^14
    # F(r)/r = 24ε[2σ¹²/r¹⁴ - σ⁶/r⁸]
    return (24.0 * epsilon / (sigma^2)) * (2.0 * sr14 - sr8)
end

function compute_force_div_r_bins(r_bin, num_bins, rlow, epsilon, sigma, r_cut)
    """
    Precompute force_div_r for each bin center
    """
    force_div_r = zeros(Float64, num_bins)

    for bin_idx in 1:num_bins
        r = rlow + (bin_idx - 0.5) * r_bin  # bin center
        if r <= r_cut && r > 0.0
            force_div_r[bin_idx] = lennard_jones_force_div_r(r, epsilon, sigma)
        else
            force_div_r[bin_idx] = 0.0
        end
    end

    return force_div_r
end

function read_particle_positions(filename)
    """
    Read particle positions from file.
    Expected format: one particle per line, x y [z] coordinates
    """
    data = readdlm(filename, comments=true)

    num_particles = size(data, 1)
    dimensions = size(data, 2)

    particle_positions = zeros(Float64, num_particles, dimensions)
    for i in 1:num_particles
        for d in 1:dimensions
            particle_positions[i, d] = data[i, d]
        end
    end

    return particle_positions
end


function gr_force_from_dir(directory::String, box_length::Float64, r_bin::Float64, num_bins::Int, force_div_r::Vector{Float64}, rlow::Float64, r_cut::Float64, method::String="out")
    """
    Compute Borgis g(r) from all configuration files in a directory
    """
    borgis_gr_total = zeros(Float64, num_bins)
    file_count = 0

    start_time = time()
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if endswith(file, ".dat")  # Assuming .dat files contain particle positions
                filepath = joinpath(root, file)
                println("Processing file: $filepath")

                particle_positions = read_particle_positions(filepath)

                borgis_gr_unnormalized = grForce_notNorm_svectorized(particle_positions, box_length, r_bin, num_bins, force_div_r, rlow, r_cut, method)

                borgis_gr_total .+= borgis_gr_unnormalized
                file_count += 1
            end
        end
    end
    println("Processed $file_count files in $(time() - start_time) seconds.")

    if file_count > 0
        borgis_gr_average = borgis_gr_total ./ file_count
        return borgis_gr_average, borgis_gr_total .^ 2 ./ file_count .- borgis_gr_average .^ 2
    else
        throw(ArgumentError("No valid configuration files found in directory: $directory"))
    end
end

function gr_force_from_dir_parallel(directory::String, box_length::Float64, r_bin::Float64, num_bins::Int, force_div_r::Vector{Float64}, rlow::Float64, r_cut::Float64, method::String="out")
    """
    Compute Borgis g(r) from all configuration files in a directory - Parallel version
    """

    # Collect all file paths first
    file_paths = String[]
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if endswith(file, ".dat")
                push!(file_paths, joinpath(root, file))
            end
        end
    end

    if isempty(file_paths)
        throw(ArgumentError("No valid configuration files found in directory: $directory"))
    end

    # Process files in parallel
    results = ThreadsX.map(file_paths) do filepath
        particle_positions = read_particle_positions(filepath)
        grForce_notNorm_svectorized(particle_positions, box_length, r_bin, num_bins, force_div_r, rlow, r_cut, method)
    end

    # Combine results
    borgis_gr_total = sum(results)
    file_count = length(file_paths)
    borgis_gr_average = borgis_gr_total ./ file_count

    # Calculate variance
    squared_sum = sum(x -> x .^ 2, results)
    variance = squared_sum ./ file_count .- borgis_gr_average .^ 2

    (borgis_gr_average, variance)

    return borgis_gr_average, variance
end

function read_particle_positions_binary(filename)
    """
    Read particle positions from binary format - 5-10x faster
    """
    open(filename) do io
        num_particles = read(io, Int32)
        dimensions = read(io, Int32)
        
        # Read all data at once into a pre-allocated array
        particle_positions = Matrix{Float64}(undef, num_particles, dimensions)
        read!(io, particle_positions)
        
        return particle_positions
    end
end

function convert_to_binary(ascii_dir::String, binary_dir::String)
    """
    Convert all .dat files to binary format
    """
    mkpath(binary_dir)
    
    for (root, dirs, files) in walkdir(ascii_dir)
        for file in files
            if endswith(file, ".dat")
                ascii_path = joinpath(root, file)
                binary_path = joinpath(binary_dir, replace(file, ".dat" => ".bin"))
                
                # Read ASCII
                data = readdlm(ascii_path, comments=true)
                num_particles, dimensions = size(data)
                
                # Write binary
                open(binary_path, "w") do io
                    write(io, Int32(num_particles))    # Header: particle count
                    write(io, Int32(dimensions))       # Header: dimensions
                    write(io, Float64.(data))          # Binary data
                end
                
                println("Converted: $ascii_path → $binary_path")
            end
        end
    end
end


function gr_force_from_dir_parallel_binary(directory::String, box_length::Float64, r_bin::Float64, 
    num_bins::Int, force_div_r::Vector{Float64}, 
    rlow::Float64, r_cut::Float64, method::IntegrationMethod; core_strength::Int=13)
    """
    Optimized version using binary files
    """

    # Collect all binary file paths
    file_paths = String[]
    for (root, dirs, files) in walkdir(directory)
        for file in files
            if endswith(file, ".bin")  # Now looking for binary files
            push!(file_paths, joinpath(root, file))
            end
        end
    end

    if isempty(file_paths)
    throw(ArgumentError("No binary files found in directory: $directory"))
    end

    results = ThreadsX.map(file_paths) do filepath
        particle_positions = read_particle_positions_binary(filepath)
        res = grForce_notNorm_svectorized(particle_positions, box_length, r_bin, num_bins, force_div_r, rlow, r_cut, method; core_strength=core_strength)
        return res
    end

    # Combine results (same as before)
    if method != Both()
        borgis_gr_total = sum(results)
        file_count = length(file_paths)
        borgis_gr_average = borgis_gr_total ./ file_count

        squared_sum = sum(x -> x .^ 2, results)
        variance = squared_sum ./ file_count .- borgis_gr_average .^ 2
        return borgis_gr_average, variance
    elseif method == Both()
        # 1. Setup metadata
        sample_pos = read_particle_positions_binary(file_paths[1])
        N, dim = size(sample_pos)
        prefactor = compute_prefactor(N, box_length, dim)
        inv_prefactor = 1.0 / prefactor

        # 2. Pre-allocate matrices to avoid push! and reduce(hcat)
        num_results = length(results)
        num_bins = length(results[1])
        
        # We can compute grOpt directly to save memory
        grOpt_sum = zeros(Float64, num_bins)
        half_bins = div(num_bins, 2)
        
        # Pre-allocate a reusable buffer for cumsum calculations
        temp_cumsum = Vector{Float64}(undef, num_bins)

        for result in results
            # Calculate gr0 (forward cumsum)
            cumsum!(temp_cumsum, result)
            
            # Calculate grInf (backward cumsum logic)
            # Instead of multiple reverses, we use the fact that:
            # rev_cumsum[end] is the total sum.
            total_sum = temp_cumsum[end]
            
            for i in 1:num_bins
                # Logic: gr0 is temp_cumsum[i] * inv_prefactor
                # Logic: grInf is 1 - (total_sum - (i > 1 ? temp_cumsum[i-1] : 0)) * inv_prefactor
                
                val_gr0 = temp_cumsum[i] * inv_prefactor
                
                # Efficiently compute the reverse cumsum value without creating new arrays
                current_rev_cumsum = total_sum - (i > 1 ? temp_cumsum[i-1] : 0)
                val_grInf = 1.0 - (current_rev_cumsum * inv_prefactor)
                
                # Apply lambda0 logic: use gr0 for first half, grInf for second half
                if i <= half_bins
                    grOpt_sum[i] += val_gr0
                else
                    grOpt_sum[i] += val_grInf
                end
            end
        end

        # Average the sum and scale back by prefactor
        grOpt_estimator = (grOpt_sum ./ num_results) .* prefactor
        
        # Generate lambda0 for return
        lambda0 = vcat(zeros(half_bins), ones(num_bins - half_bins))
        
        return grOpt_estimator, lambda0
    end
end

# Map the string to the type
function get_method_type(method::String)
    if method == "in"    return In()
    elseif method == "out"  return Out()
    elseif method == "both" return Both()
    else
        throw(ArgumentError("Unknown method: $method"))
    end
end


function compute_prefactor(N_particles, box_length, dimension)
    if dimension == 2
        return (2 * π * N_particles^2) / (box_length^dimension)
    elseif dimension == 3
        return (4 * π * N_particles^2) / (box_length^dimension)
    else
        throw(ArgumentError("Unsupported dimension: $dimension. Only 2D and 3D are supported."))
    end
end

function main()
    dummy_params = Dict(
        "directory" => "inputs",
        "box_length" => 10.0,
        "r_bin" => 0.01,
        "num_bins" => 1000,
        "r_low" => 0.1,
        "r_cut" => 5.0,
        "epsilon" => 1.0,
        "sigma" => 1.0,
        "method" => "out"
    )

    params = JSON.parsefile(joinpath(dummy_params["directory"], "params.json"))

    N = params["N_particles"]
    L_box = params["L_box"]    
    r_bin = params["bin_width"]
    r_low = params["r_low"]
    r_cut = params["r_high"]
    epsilon = dummy_params["epsilon"]
    sigma = dummy_params["sigma"]
    method = params["method_force_formula"]
    T = params["Temperature"]
    dimension = params["dimensions"]

    num_bins = params["qdim_max"]

    config_dir_binary = joinpath(dummy_params["directory"], "configs_bin")


    # Calculation

    force_div_r = compute_force_div_r_bins(r_bin, num_bins, r_low, epsilon, sigma, r_cut)

    borgis_gr_average, variance = gr_force_from_dir_parallel_binary(
        config_dir_binary, L_box, r_bin, num_bins, force_div_r, r_low, r_cut, method
    )

    # Normalization
    prefactor = compute_prefactor(N, L_box, dimension)  
    if method == "out"
        gr_borgis = 1 .- (borgis_gr_average ./ prefactor)
    else
        gr_borgis = borgis_gr_average ./ prefactor
    end

    radii = [(i - 0.5) * r_bin for i in 1:num_bins]

    # Plot results
    plot(radii, gr_borgis, xlabel="r", ylabel="g(r)", title="Borgis g(r) - Force Method", legend=false)
    savefig(joinpath(dummy_params["directory"], "gr_borgis_force.png"))
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
