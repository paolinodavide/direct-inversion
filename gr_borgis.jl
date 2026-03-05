using LoopVectorization
using StaticArrays
using Revise
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


@inline function wrap_pbc_distances(separation::Float64, box_length::Float64, inv_box_length::Float64)::Float64
    """Apply minimum image convention for periodic boundary conditions."""
    return separation - box_length * round(separation * inv_box_length)
end

@inline function pbc_distance(pos_i::SVector{D, Float64}, pos_j::SVector{D, Float64}, box_size::Float64) where D
    """Compute the minimum image distance vector between two positions under PBC."""
    distance = pos_i - pos_j
    distance -= box_size * round.(distance / box_size)
    return distance, dot(distance, distance)
end


""" 
Compute g(r) using the Force prescription of Borgis et al.
Method specifyes the integration direction:
- "in": integrate from 0 outward
- "out": integrate from infty inward
- "both": return raw integrand
"""
function grForce_notNorm_svectorized(particle_positions::Array{Float64, D},
    box_length::Float64,
    bin_width::Float64,
    num_bins_gr::Int,
    force_over_r::Vector{Float64},
    r_min_interaction::Float64,
    r_cutoff_interaction::Float64,
    method::IntegrationMethod;
    core_strength::Int=13) where D
    
    num_particles, num_dimensions = size(particle_positions)
    inv_bin_width = 1.0 / bin_width

    positions = [SVector{D, Float64}(particle_positions'[:, i]) for i in 1:num_particles]
    total_forces = zeros(SVector{D, Float64}, num_particles)
    #println("Converted pos to SVectors")

    evaluate_total_forces!(total_forces, positions, box_length, force_over_r, r_min_interaction, r_cutoff_interaction, bin_width; core_strength=core_strength)

    borgis_contributions = compute_borgis_contributions(positions, total_forces, box_length, inv_bin_width, num_bins_gr)

    return integrate_borgis_contributions(borgis_contributions, method)
end

@inline function compute_borgis_contributions(positions::Vector{SVector{D, Float64}}, total_forces::Vector{SVector{D, Float64}}, box_length::Float64, inv_bin_width::Float64, num_bins_gr::Int) where D
    num_particles = length(positions)
    borgis_contributions = zeros(Float64, num_bins_gr)
    @inbounds for i in 1:num_particles-1
        pos_i = positions[i]
        force_i = total_forces[i]

        for j in i+1:num_particles
            pos_j = positions[j]
            force_j = total_forces[j]

            rVec_ij, r2_ij= pbc_distance(pos_i, pos_j, box_length)
            if iszero(r2_ij)
                continue
            end

            r_ij = sqrt(r2_ij)
            radial_bin_index = floor(Int, r_ij * inv_bin_width)
            target_bin = clamp(radial_bin_index + 1, 1, num_bins_gr)

            force_diff = force_i - force_j
            borgis_delta = borgis_delta_calculation(force_diff, rVec_ij, r2_ij, r_ij)

            borgis_contributions[target_bin] += borgis_delta
        end
    end
    return borgis_contributions
end

# "in" method
@inline function integrate_borgis_contributions(contributions::Vector{Float64}, ::In)
    cumsum!(contributions, contributions)
end

# "out" method
@inline function integrate_borgis_contributions(contributions::Vector{Float64}, ::Out)
    reverse!(contributions)
    cumsum!(contributions, contributions)
    reverse!(contributions)
end

# "both" method (identity)
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

@inline function evaluate_total_forces!(total_forces::Vector{SVector{D, Float64}},
    positions::Vector{SVector{D, Float64}},
    box_length::Float64, force_over_r::Vector{Float64},
    r_min_interaction::Float64, 
    r_cutoff_interaction::Float64, 
    bin_width::Float64;
    core_strength::Int=0) where D

    num_particles = length(positions)
    inv_bin_width = 1.0 / bin_width

    @inbounds for i in 1:num_particles-1
        pos_i = positions[i]

        for j in i+1:num_particles
            pos_j = positions[j]

            rVec_ij, r2_ij= pbc_distance(pos_i, pos_j, box_length)
            if r2_ij == 0.0 || r2_ij > r_cutoff_interaction^2
                continue
            end

            r_ij = sqrt(r2_ij)


            #f_magnitude::Float64 = 0.0
            f_magnitude = if r_ij < r_min_interaction
                force_magnitude_below_rmin(r_ij, r_min_interaction, force_over_r; core_strength=core_strength)
            else
                force_magnitude_between_bins(r_ij, r_min_interaction, r_cutoff_interaction, force_over_r, inv_bin_width)
            end

            total_forces[i] += f_magnitude * rVec_ij
            total_forces[j] -= f_magnitude * rVec_ij
        end
    end
end

@inline function force_magnitude_below_rmin(r_ij::Float64, r_min::Float64, force_over_r::Vector{Float64}; core_strength::Int=0)::Float64
    if core_strength == 0
        return 0.0
    elseif core_strength == 1
        a = r_min / r_ij
        return a * force_over_r[1] + a*(1-a)*inv_bin_width * (force_over_r[2] - force_over_r[1] + bin_width * force_over_r[1]/r_min)
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


function grForce_notNorm(particle_positions::Matrix{Float64},
    box_length::Float64,
    bin_width::Float64,
    num_bins::Int,
    force_over_r::Vector{Float64},
    r_min::Float64,
    r_cutoff::Float64,
    method::String="out")
    num_particles, num_dimensions = size(particle_positions)
    total_forces = zeros(Float64, num_particles, num_dimensions)
    borgis_contributions = zeros(Float64, num_bins)

    # Precompute cutoff and binning parameters
    cutoff_squared = r_cutoff^2
    inv_bin_width = 1.0 / bin_width
    min_bin_skip = Int(floor(r_min * inv_bin_width))
    inv_box_length = 1.0 / box_length

    # Transpose positions for cache-efficient column-major access
    positions_transposed = particle_positions'

    # Phase 1: Compute net forces on all particles
    @inbounds for particle_i in 1:num_particles
        position_i = @view positions_transposed[:, particle_i]

        for particle_j in (particle_i+1):num_particles
            position_j = @view positions_transposed[:, particle_j]

            # Compute minimum image separation with periodic boundary conditions
            separation_squared = 0.0
            @simd for dim in 1:num_dimensions
                displacement = position_i[dim] - position_j[dim]
                displacement = wrap_pbc_distances(displacement, box_length, inv_box_length)
                separation_squared += displacement * displacement
            end

            # Skip particles outside cutoff or at same position
            (separation_squared == 0.0 || separation_squared > cutoff_squared) && continue

            interparticle_distance = sqrt(separation_squared)

            # Interpolate force_over_r table
            force_magnitude = if interparticle_distance < r_min
                force_over_r[1] + (interparticle_distance - r_min) * inv_bin_width * (force_over_r[2] - force_over_r[1])
            elseif interparticle_distance < r_cutoff
                radial_bin_index = Int(floor(interparticle_distance * inv_bin_width))
                force_table_index = radial_bin_index - min_bin_skip + 1
                interpolation_weight = interparticle_distance * inv_bin_width - radial_bin_index
                (1.0 - interpolation_weight) * force_over_r[force_table_index] + interpolation_weight * force_over_r[force_table_index+1]
            else
                0.0
            end


            # Apply forces to both particles (Newton's 3rd law)
            @simd for dim in 1:num_dimensions
                displacement = position_i[dim] - position_j[dim]
                displacement = wrap_pbc_distances(displacement, box_length, inv_box_length)
                total_forces[particle_i, dim] += force_magnitude * displacement
                total_forces[particle_j, dim] -= force_magnitude * displacement
            end
        end
    end

    # Phase 2: Compute Borgis integrand contributions
    @inbounds for particle_i in 1:num_particles
        position_i = @view positions_transposed[:, particle_i]
        force_on_i = @view total_forces[particle_i, :]

        for particle_j in (particle_i+1):num_particles
            position_j = @view positions_transposed[:, particle_j]
            force_on_j = @view total_forces[particle_j, :]

            # Compute separation vector and relative force
            separation_squared = 0.0
            force_dot_displacement = 0.0

            @simd for dim in 1:num_dimensions
                displacement = position_i[dim] - position_j[dim]
                displacement = wrap_pbc_distances(displacement, box_length, inv_box_length)
                separation_squared += displacement * displacement

                # Compute relative force difference
                relative_force = force_on_i[dim] - force_on_j[dim]
                force_dot_displacement += relative_force * displacement
            end

            # Skip overlapping particles
            (separation_squared == 0.0) && continue

            interparticle_distance = sqrt(separation_squared)
            radial_bin_index = Int(floor(interparticle_distance * inv_bin_width))

            # Determine target bin with bounds checking
            target_bin = clamp(radial_bin_index + 1, 1, num_bins)

            # Compute Borgis delta function: ∇·F_irreducible component
            if num_dimensions == 2
                borgis_delta = force_dot_displacement / separation_squared
            elseif num_dimensions == 3
                borgis_delta = force_dot_displacement / (separation_squared * interparticle_distance)
            else
                borgis_delta = force_dot_displacement / interparticle_distance^num_dimensions
            end
            borgis_contributions[target_bin] += borgis_delta
        end
    end

    # Phase 3: Compute cumulative integral based on method
    if method == "in"
        # Cumulative integral from r_min outward
        cumsum!(borgis_contributions, borgis_contributions)
        return borgis_contributions
    elseif method == "out"
        # Cumulative integral from r_cutoff inward  
        reverse!(borgis_contributions)
        cumsum!(borgis_contributions, borgis_contributions)
        reverse!(borgis_contributions)
        return borgis_contributions
    elseif method == "both"
        # Return the raw integrand
        return borgis_contributions
    else
        throw(ArgumentError("Invalid integration method. Choose 'in', 'out', or 'both'."))
    end
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
