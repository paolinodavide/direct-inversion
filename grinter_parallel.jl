#forceIBI/grinter_parallel.jl
using JSON
using DelimitedFiles
using LinearAlgebra
using ArgParse

include("utils.jl")
include("gr_borgis.jl")


"""
Run the inversion loop for the iterative Boltzmann inversion method.
This function parses input parameters, loads target data, initializes potentials, and runs the inversion loop.
"""
function main()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--directory", "-d"
            help = "Directory containing input folder with params.json and target g(r) file. Default is current directory."
            arg_type = String
            default = "."
    end
    parsed_args = parse_args(s)
    HomeDir = parsed_args["directory"]
    
    # Load parameters
    params = JSON.parsefile(joinpath(HomeDir, "inputs/params.json"))
    
    # Extract parameters with clear grouping
    N_particles = params["N_particles"]::Int
    L_box = Float64(params["L_box"])
    dimensions = params["dimensions"]::Int
    T = Float64(params["Temperature"])

    # Radial parameters
    r_low = params["r_low"]::Float64
    r_high = params["r_high"]::Float64
    bin_width = params["bin_width"]::Float64
    binlow = floor(Int, r_low / bin_width) + 1
    binhigh = floor(Int, r_high / bin_width) + 1
    max_distance = min(L_box / 2.0, 10.0)
    num_bins_gr = floor(Int, max_distance / bin_width)

    # File paths
    path_target = joinpath(HomeDir, "inputs")
    config_dir = joinpath(path_target, params["config_dir"]::String)
    target_file = params["target_gr_file"]::String
    initial_pot = params["init_pot_type"]::String
    target_pot = params["target_pot_type"]::String
    number_config = params["n_inversion_snapshots"]::Int
    wt_file_path = joinpath(path_target, params["wt_file"]::String)
    
    # Optimization parameters
    max_iter = params["max_iter"]::Int
    learning_rate = Float64(params["learning_rate"])
    target_tol = Float64(params["target_precision"])
    iteration_tol = Float64(params["iteration_precision"])
    method_force_formula = params["method_force_formula"]::String
    method_type = get_method_type(method_force_formula)
    core_strength = Int(params["core_strength"]) #Default is Integer but can be Float. Adjust as needed.
    shift_gr = params["shift_gr"]::Bool

    # Create output directory
    mkpath(joinpath(HomeDir, "outputs"))

    # Load target data
    data = readdlm(joinpath(path_target, target_file), comments=true)
    r_values, full_target_gr = data[:, 1], data[:, 2]
    binlow = findfirst(r -> r >= r_low, r_values)
    binhigh = findlast(r -> r <= r_high, r_values)

    # Then define r_range based on what you actually extracted
    r_range = r_values[binlow:binhigh]
    r_low = r_range[1]
    r_high = r_range[end]
    if length(r_range) != length(r_values[binlow:binhigh])
        @warn "Mismatch in r_range length and target g(r) range length"
        println("binlow: $binlow, binhigh: $binhigh")
        binhigh = binhigh - 1
        println("Adjusted binhigh to $binhigh")
    end
    save_gr_data(r_values, full_target_gr, joinpath(HomeDir, "outputs/gr_target.dat"))

    gr_target = @view full_target_gr[binlow:binhigh]
    gr_current = copy(gr_target)
    
    # Initialize potentials
    βu_target = get_potential_from_name(target_pot, T, r_range)
    f_target = f_over_r_from_potential(βu_target, r_low, bin_width)
    
    βu_current = if initial_pot == "mean_force"
        - log.(abs.(gr_target))
    else
        get_potential_from_name(initial_pot, T, r_range)
    end
    βu_current .-= βu_current[end]
    f_current = f_over_r_from_potential(βu_current, r_low, bin_width)

    # Handle binary configuration files
    config_dir = ensure_binary_configs(config_dir, wt_file_path, number_config)

    
    # Precompute constants
    prefactor = compute_prefactor(N_particles, L_box, dimensions)
    delta_target = gr_target[1]
    
    # Initialize tracking
    start_time = time()
    convergence_file = joinpath(HomeDir, "outputs/convergence_data.dat")
    initialize_convergence_file(convergence_file)
    
    # Main optimization loop
    for iteration in 0:max_iter
        gr_old = copy(gr_current)
        
        # βu_t → gr_t
        gr_notNorm, _ = gr_force_from_dir_parallel_binary(
            config_dir, L_box, bin_width, num_bins_gr, 
            f_current, r_low, r_high, method_type;
            core_strength=core_strength
        )
        if any(isnan.(gr_notNorm))
            @error "Invalid g(r) values detected. Check force calculations and consider increasing r_low."
            exit(1)
        end

        # Normalize and update gr
        if method_force_formula == "out"
            gr_normalized = 1.0 .- gr_notNorm ./ prefactor
        else
            gr_normalized = gr_notNorm ./ prefactor
        end

        # Rescale gr
        gr_current = @view gr_normalized[binlow:binhigh]
        g_min = minimum(gr_current)
        g_end = gr_current[end]
        save_iteration_data(iteration, r_range, gr_current, βu_current, f_current, HomeDir)

        # Check convergence
        error, iteration_diff, potential_increase = compute_convergence_metrics(gr_current, gr_target, gr_old)
        
        @info "Iteration $iteration" error=error iteration_diff=iteration_diff potential_increase=potential_increase g_min=g_min delta=(gr_current[1] - delta_target)
        
        # Save iteration data      
        append_convergence_data(convergence_file, iteration, time() - start_time, error, iteration_diff, potential_increase, g_end, g_min)
        
        if error < target_tol || iteration_diff < iteration_tol || iteration == max_iter
            @info "Convergence achieved" iteration=iteration
            save_gr_data(r_values, gr_normalized, joinpath(HomeDir, "outputs/gr_final.dat"))
            
            break
        end
        save_gr_data(r_values, gr_normalized, joinpath(HomeDir, "outputs/gr_$(iteration).dat"))

        # u_t → u_t+1
        update_potential!(βu_current, gr_current, gr_target, learning_rate, shift_gr)
        f_current = f_over_r_from_potential(βu_current, r_low, bin_width)
    end
    
    # Save final target data
    save_target_data(r_range, gr_target, βu_target, f_target, joinpath(HomeDir, "outputs/iteration_-1.dat"))
    cp(joinpath(HomeDir, "inputs/params.json"), joinpath(HomeDir, "outputs/00params.json"), force=true)
    @info "Optimization completed in $(time() - start_time) seconds"
end

"""
Updates the potential based on the current and target g(r) values using the Schommer update rule.
"""
function update_potential!(
    βu_t::AbstractVector{Float64}, 
    gr_t::AbstractVector{Float64}, 
    gr_tgt::AbstractVector{Float64}, 
    learning_rate::Float64, 
    correct_offset::Bool=false; 
    small_number::Float64=1e-10)
    # Ensure all vectors have the same dimensions
    if !(length(βu_t) == length(gr_t) == length(gr_tgt))
        throw(DimensionMismatch("Vectors βu_t, gr_t, and gr_tgt must have the same length."))
    end

    min_index = findmin(gr_t)[2]
    g_min = gr_t[min_index]
    Delta = 0.0
    if correct_offset
        Delta = g_min - gr_tgt[min_index] 
    end

    @. gr_t = gr_t - Delta 
    if minimum(gr_t) < 0 || any(isnan.(gr_t))
        @error "Negative g(r) or NaN detected. Retry with higher r_low"
        exit(1)
    end

    @. βu_t = βu_t + learning_rate * log(gr_t / gr_tgt + small_number)
    βu_t .-= βu_t[end]  
end

"""
Computes force/r from the potential using finite differences.
"""
function f_over_r_from_potential(potential::Vector{Float64}, r_low::Float64, bin_width::Float64)::Vector{Float64}
    pot_length = length(potential)
    f_over_r = similar(potential)

    @inbounds for i in 2:(pot_length - 1)
        r = r_low + (i - 1) * bin_width
        f_over_r[i] = - (potential[i + 1] - potential[i - 1]) / (2.0 * bin_width * r)
    end

    # 1. Handle the start (i = 1)
    r_start = r_low
    # Force = -dV/dr, then divide by r
    f_over_r[1] = -(-3.0 * potential[1] + 4.0 * potential[2] - potential[3]) / (2.0 * bin_width * r_start)

    # 2. Handle the end (i = pot_length)
    r_end = r_low + (pot_length - 1) * bin_width
    f_over_r[end] = -(3.0 * potential[end] - 4.0 * potential[end-1] + potential[end-2]) / (2.0 * bin_width * r_end)

    return f_over_r
end

"""
Computes convergence metrics: error, iteration difference, and potential increase.
"""
function compute_convergence_metrics(gr_current, gr_target, gr_old)
    error = sum((gr_current .- gr_target).^2) / length(gr_target)
    iteration_diff = sum((gr_current .- gr_old).^2) / length(gr_current)
    potential_increase = sum((log.(abs.(gr_current ./ gr_target))).^2) / length(gr_current)
    return error, iteration_diff, potential_increase
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end