using LoopVectorization
using StaticArrays
using Statistics
using ThreadsX
using JSON
using DelimitedFiles
using LinearAlgebra
using Plots

# ---------------------------------------------------------------------------
# Integration method dispatch types
# ---------------------------------------------------------------------------

abstract type IntegrationMethod end
struct In   <: IntegrationMethod end
struct Out  <: IntegrationMethod end
struct Both <: IntegrationMethod end

function get_method_type(method::String)::IntegrationMethod
    method == "in"   && return In()
    method == "out"  && return Out()
    method == "both" && return Both()
    throw(ArgumentError("Unknown integration method: \"$method\". Use \"in\", \"out\", or \"both\"."))
end
get_method_type(m::IntegrationMethod) = m   # passthrough for already-typed values


# ---------------------------------------------------------------------------
# Periodic boundary conditions
# ---------------------------------------------------------------------------

@inline function pbc_distance(pos_i::SVector{D,Float64}, pos_j::SVector{D,Float64},
                               box_length::Float64) where {D}
    Δ = pos_i - pos_j
    Δ -= box_length .* round.(Δ ./ box_length)
    return Δ, dot(Δ, Δ)
end


# ---------------------------------------------------------------------------
# Borgis g(r) integrand (dimension-generic)
# ---------------------------------------------------------------------------

@inline borgis_weight(force_diff::SVector{2,Float64}, r_vec::SVector{2,Float64}, r2, r) =
    dot(force_diff, r_vec) / r2

@inline borgis_weight(force_diff::SVector{3,Float64}, r_vec::SVector{3,Float64}, r2, r) =
    dot(force_diff, r_vec) / (r2 * r)

@inline borgis_weight(force_diff::SVector{D,Float64}, r_vec::SVector{D,Float64}, r2, r) where {D} =
    dot(force_diff, r_vec) / r^D


# ---------------------------------------------------------------------------
# Integration accumulation
# ---------------------------------------------------------------------------

@inline function integrate_borgis_contributions(contrib::Vector{Float64}, ::In)
    cumsum!(contrib, contrib)
end

@inline function integrate_borgis_contributions(contrib::Vector{Float64}, ::Out)
    reverse!(contrib)
    cumsum!(contrib, contrib)
    reverse!(contrib)
end

@inline integrate_borgis_contributions(contrib::Vector{Float64}, ::Both) = contrib


# ---------------------------------------------------------------------------
# Force table helpers
# ---------------------------------------------------------------------------

@inline function force_below_rmin(r_ij::Float64, force_over_r::Vector{Float64};
                                   core_strength::Int=0)::Float64
    core_strength == 0 && return 0.0
    return force_over_r[1] * (first(force_over_r) == 0.0 ? 0.0 : 1.0)   # guarded fallback
end

# The general repulsive-core model: F/r ∝ (r_min/r)^core_strength
@inline function force_below_rmin_core(r_ij::Float64, r_min::Float64,
                                        force_over_r::Vector{Float64},
                                        core_strength::Int)::Float64
    core_strength == 0 && return 0.0
    return force_over_r[1] * (r_min / r_ij)^core_strength
end

@inline function force_between_bins(r_ij::Float64, r_min::Float64,
                                     force_over_r::Vector{Float64},
                                     inv_bin_width::Float64)::Float64
    bin_idx   = floor(Int, r_ij * inv_bin_width)
    skip      = floor(Int, r_min * inv_bin_width)
    tbl_idx   = clamp(bin_idx - skip + 1, 1, length(force_over_r) - 1)
    weight    = r_ij * inv_bin_width - bin_idx
    return muladd(weight, force_over_r[tbl_idx + 1], (1.0 - weight) * force_over_r[tbl_idx])
end


# ---------------------------------------------------------------------------
# Phase 1 – net forces
# ---------------------------------------------------------------------------

function evaluate_total_forces!(total_forces::Vector{SVector{D,Float64}},
                                  positions::Vector{SVector{D,Float64}},
                                  box_length::Float64,
                                  force_over_r::Vector{Float64},
                                  r_min::Float64,
                                  r_cutoff::Float64,
                                  bin_width::Float64;
                                  core_strength::Int=0) where {D}
    N             = length(positions)
    inv_bin_width = 1.0 / bin_width
    r_cut2        = r_cutoff^2

    @inbounds for i in 1:N-1
        pos_i = positions[i]
        for j in i+1:N
            r_vec, r2 = pbc_distance(pos_i, positions[j], box_length)
            (r2 == 0.0 || r2 > r_cut2) && continue

            r = sqrt(r2)
            f = if r < r_min
                force_below_rmin_core(r, r_min, force_over_r, core_strength)
            else
                force_between_bins(r, r_min, force_over_r, inv_bin_width)
            end

            total_forces[i] += f * r_vec
            total_forces[j] -= f * r_vec
        end
    end
end


# ---------------------------------------------------------------------------
# Phase 2 – Borgis integrand bins
# ---------------------------------------------------------------------------

function compute_borgis_contributions(positions::Vector{SVector{D,Float64}},
                                        total_forces::Vector{SVector{D,Float64}},
                                        box_length::Float64,
                                        inv_bin_width::Float64,
                                        num_bins::Int) where {D}
    N       = length(positions)
    contrib = zeros(Float64, num_bins)

    @inbounds for i in 1:N-1
        pos_i = positions[i]
        f_i   = total_forces[i]
        for j in i+1:N
            r_vec, r2 = pbc_distance(pos_i, positions[j], box_length)
            iszero(r2) && continue

            r       = sqrt(r2)
            bin     = clamp(floor(Int, r * inv_bin_width) + 1, 1, num_bins)
            Δf      = f_i - total_forces[j]
            contrib[bin] += borgis_weight(Δf, r_vec, r2, r)
        end
    end
    return contrib
end


# ---------------------------------------------------------------------------
# Main single-snapshot entry point
# ---------------------------------------------------------------------------

function grForce_notNorm_svectorized(particle_positions::Matrix{Float64},
                                      box_length::Float64,
                                      bin_width::Float64,
                                      num_bins::Int,
                                      force_over_r::Vector{Float64},
                                      r_min::Float64,
                                      r_cutoff::Float64,
                                      method::IntegrationMethod;
                                      core_strength::Int=13)
    N, D          = size(particle_positions)
    inv_bin_width = 1.0 / bin_width

    positions     = [SVector{D,Float64}(particle_positions'[:, i]) for i in 1:N]
    total_forces  = zeros(SVector{D,Float64}, N)

    evaluate_total_forces!(total_forces, positions, box_length, force_over_r,
                            r_min, r_cutoff, bin_width; core_strength=core_strength)

    contrib = compute_borgis_contributions(positions, total_forces, box_length,
                                            inv_bin_width, num_bins)

    return integrate_borgis_contributions(contrib, method)
end


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

function read_particle_positions_binary(filename::String)::Matrix{Float64}
    open(filename) do io
        np  = read(io, Int32)
        dim = read(io, Int32)
        buf = Matrix{Float64}(undef, np, dim)
        read!(io, buf)
        return buf
    end
end

function convert_to_binary(ascii_dir::String, binary_dir::String)
    mkpath(binary_dir)
    for (root, _, files) in walkdir(ascii_dir)
        for file in files
            endswith(file, ".dat") || continue
            src  = joinpath(root, file)
            dst  = joinpath(binary_dir, replace(file, ".dat" => ".bin"))
            data = readdlm(src, comments=true)
            np, dim = size(data)
            open(dst, "w") do io
                write(io, Int32(np))
                write(io, Int32(dim))
                write(io, Float64.(data))
            end
            println("Converted: $src → $dst")
        end
    end
end


# ---------------------------------------------------------------------------
# Prefactor
# ---------------------------------------------------------------------------

function compute_prefactor(N_particles::Int, box_length::Float64, dimension::Int)::Float64
    dimension == 2 && return (2π * N_particles^2) / box_length^2
    dimension == 3 && return (4π * N_particles^2) / box_length^3
    throw(ArgumentError("Unsupported dimension $dimension. Only 2 and 3 are supported."))
end


# ---------------------------------------------------------------------------
# grOpt accumulation for Both() method
# ---------------------------------------------------------------------------

function accumulate_grOpt!(grOpt_sum::Vector{Float64},
                             result::Vector{Float64},
                             inv_prefactor::Float64,
                             half_bins::Int)
    nb         = length(result)
    temp       = cumsum(result)
    total_sum  = temp[end]

    @inbounds for i in 1:nb
        if i <= half_bins
            grOpt_sum[i] += temp[i] * inv_prefactor
        else
            rev_val       = total_sum - (i > 1 ? temp[i-1] : 0.0)
            grOpt_sum[i] += 1.0 - rev_val * inv_prefactor
        end
    end
end


# ---------------------------------------------------------------------------
# Primary public API – parallel binary directory processing
# ---------------------------------------------------------------------------

function gr_force_from_dir_parallel_binary(directory::String,
                                             box_length::Float64,
                                             r_bin::Float64,
                                             num_bins::Int,
                                             force_div_r::Vector{Float64},
                                             rlow::Float64,
                                             r_cut::Float64,
                                             method;          # String or IntegrationMethod
                                             core_strength::Int=13)
    method_type = get_method_type(method)

    # Collect binary files
    file_paths = [joinpath(root, f)
                  for (root, _, files) in walkdir(directory)
                  for f in files if endswith(f, ".bin")]

    isempty(file_paths) &&
        throw(ArgumentError("No .bin files found in: $directory"))

    # Parallel per-snapshot computation
    results = ThreadsX.map(file_paths) do path
        pos = read_particle_positions_binary(path)
        grForce_notNorm_svectorized(pos, box_length, r_bin, num_bins, force_div_r,
                                     rlow, r_cut, method_type; core_strength=core_strength)
    end

    nf = length(results)

    if !(method_type isa Both)
        total    = sum(results)
        average  = total ./ nf
        variance = sum(x -> x .^ 2, results) ./ nf .- average .^ 2
        return average, variance

    else  # Both() – return grOpt estimator + lambda0 mask
        sample_pos      = read_particle_positions_binary(file_paths[1])
        N, dim          = size(sample_pos)
        prefactor       = compute_prefactor(N, box_length, dim)
        inv_prefactor   = 1.0 / prefactor
        half_bins       = div(num_bins, 2)

        grOpt_sum = zeros(Float64, num_bins)
        for r in results
            accumulate_grOpt!(grOpt_sum, r, inv_prefactor, half_bins)
        end

        grOpt_estimator = (grOpt_sum ./ nf) .* prefactor
        lambda0         = vcat(zeros(half_bins), ones(num_bins - half_bins))
        return grOpt_estimator, lambda0
    end
end


# ---------------------------------------------------------------------------
# Lennard-Jones force table builder
# ---------------------------------------------------------------------------

function lennard_jones_force_div_r(r::Float64, epsilon=1.0, sigma=1.0)::Float64
    r == 0.0 && return 0.0
    sr  = sigma / r
    sr8 = sr^8; sr14 = sr^14
    return (24.0 * epsilon / sigma^2) * (2.0 * sr14 - sr8)
end

function compute_force_div_r_bins(r_bin::Float64, num_bins::Int, rlow::Float64,
                                    epsilon::Float64, sigma::Float64, r_cut::Float64)
    [begin
        r = rlow + (i - 0.5) * r_bin
        (r > 0.0 && r <= r_cut) ? lennard_jones_force_div_r(r, epsilon, sigma) : 0.0
     end for i in 1:num_bins]
end


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

function main()
    input_dir = "inputs"
    params    = JSON.parsefile(joinpath(input_dir, "params.json"))

    N         = params["N_particles"]
    L_box     = params["L_box"]
    r_bin     = params["bin_width"]
    r_low     = params["r_low"]
    r_cut     = params["r_high"]
    method    = params["method_force_formula"]
    dimension = params["dimensions"]
    num_bins  = params["qdim_max"]

    epsilon   = 1.0
    sigma     = 1.0
    config_dir_binary = joinpath(input_dir, "configs_bin")

    force_div_r = compute_force_div_r_bins(r_bin, num_bins, r_low, epsilon, sigma, r_cut)

    borgis_avg, variance = gr_force_from_dir_parallel_binary(
        config_dir_binary, L_box, r_bin, num_bins, force_div_r, r_low, r_cut, method
    )

    prefactor = compute_prefactor(N, L_box, dimension)
    gr_borgis = method == "out" ? 1.0 .- (borgis_avg ./ prefactor) :
                                  borgis_avg ./ prefactor

    radii = [(i - 0.5) * r_bin for i in 1:num_bins]
    plot(radii, gr_borgis, xlabel="r", ylabel="g(r)",
         title="Borgis g(r) - Force Method", legend=false)
    savefig(joinpath(input_dir, "gr_borgis_force.png"))
end


if abspath(PROGRAM_FILE) == @__FILE__ 
    main()
end