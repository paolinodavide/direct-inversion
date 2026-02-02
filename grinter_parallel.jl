using JSON
using DelimitedFiles
using LinearAlgebra
using Plots

include("utils.jl")
include("gr_borgis.jl")

function main()
    # Load parameters
    params = JSON.parsefile("inputs/params.json")
    
    # Extract parameters with clear grouping
    N_particles = params["N_particles"]::Int
    L_box = params["L_box"]::Float64
    dimensions = params["dimensions"]::Int
    T = params["Temperature"]::Float64

    # Radial parameters
    r_low = params["r_low"]::Float64
    r_high = params["r_high"]::Float64
    bin_width = params["bin_width"]::Float64
    binlow = floor(Int, r_low / bin_width) + 1
    binhigh = floor(Int, r_high / bin_width) + 1
    max_distance = min(L_box / 2.0, 10.0)
    num_bins_gr = floor(Int, max_distance / bin_width)

    # File paths
    path_target = "inputs"
    config_dir = joinpath(path_target, params["config_dir"]::String)
    target_file = params["target_gr_file"]::Stringlength
    initial_pot = params["init_pot_type"]::String
    target_pot = params["target_pot_type"]::String
    number_config = params["n_inversion_snapshots"]::Int
    wt_file_path = joinpath(path_target, params["wt_file"]::String)
    
    # Optimization parameters
    max_iter = params["max_iter"]::Int
    learning_rate = params["learning_rate"]::Float64
    target_tol = params["target_precision"]::Float64
    iteration_tol = params["iteration_precision"]::Float64
    method_force_formula = params["method_force_formula"]::String
    core_strength = params["core_strength"]::Int
    shift_gr = params["shift_gr"]::Bool

    # Create output directory
    mkpath("outputs")

    # Load target data
    data = readdlm(joinpath(path_target, target_file), comments=true)
    r_values, full_target_gr = data[:, 1], data[:, 2]
    r_range = r_low:bin_width:r_high
    #print length r range and r_Values[binlow:binhigh]
    if length(r_range) != length(r_values[binlow:binhigh])
        binhigh = binlow + length(r_range) - 1
    end
    r_range
    save_gr_data(r_values, full_target_gr, "gr_target.dat")

    gr_target = @view full_target_gr[binlow:binhigh]
    gr_current = copy(gr_target)
    
    # Initialize potentials
    βu_target = get_potential_from_name(target_pot, T, r_low, r_high, bin_width)
    f_target = f_over_r_from_potential(βu_target, r_low, bin_width)
    
    βu_current = if initial_pot == "mean_force"
        - log.(abs.(gr_target))
    else
        get_potential_from_name(initial_pot, T, r_low, r_high, bin_width)
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
    convergence_file = "outputs/convergence_data.dat"
    initialize_convergence_file(convergence_file)
    
    # Main optimization loop
    for iteration in 0:max_iter
        gr_old = copy(gr_current)
        
        # βu_t → gr_t
        gr_notNorm, _ = gr_force_from_dir_parallel_binary(
            config_dir, L_box, bin_width, num_bins_gr, 
            f_current, r_low, r_high, method_force_formula; 
            core_strength=core_strength
        )
        
        # Normalize and update gr
        if method_force_formula == "out"
            gr_normalized = 1.0 .- gr_notNorm ./ prefactor
        else
            gr_normalized = gr_notNorm ./ prefactor
        end

        # Rescale gr
        gr_current = @view gr_normalized[binlow:binhigh]
        g_min = minimum(gr_current)
        save_iteration_data(iteration, r_range, gr_current, βu_current, f_current)

        # Check convergence
        error, iteration_diff = compute_convergence_metrics(gr_current, gr_target, gr_old)
        
        @info "Iteration $iteration" error=error iteration_diff=iteration_diff g_min=g_min delta=(gr_current[1] - delta_target)
        
        # Save iteration data      
        append_convergence_data(convergence_file, iteration, time() - start_time, error, iteration_diff, delta_target, g_min)
        
        if error < target_tol || iteration_diff < iteration_tol
            @info "Convergence achieved" iteration=iteration

            save_gr_data(r_values, gr_normalized, "gr_final.dat")
            break
        end

        # u_t → u_t+1
        update_potential!(βu_current, gr_current, gr_target, learning_rate, shift_gr)
        f_current = f_over_r_from_potential(βu_current, r_low, bin_width)
    end
    
    # Save final target data
    save_target_data(r_range, gr_target, βu_target, f_target)
    cp("inputs/params.json", "outputs/00params.json", force=true)
    @info "Optimization completed in $(time() - start_time) seconds"
end


function update_potential!(βu_t, gr_t, gr_tgt, learning_rate, correct_offset::Bool=false; small_number::Float64=1e-10)
    min_index = findmin(gr_t)[2]
    g_min = gr_t[min_index]
    Delta = 0.0
    correct_offset && (Delta = g_min - gr_tgt[min_index])

    @. gr_t = gr_t - Delta
    if minimum(gr_t) < 0
        @error "Negative g(r). Retry with higher r_low"
        exit(1)
    end

    @. βu_t = βu_t + learning_rate * log(gr_t / gr_tgt )
    βu_t .-= βu_t[end]  
end

function compute_convergence_metrics(gr_current, gr_target, gr_old)
    error = sum((gr_current .- gr_target).^2) / length(gr_target)
    iteration_diff = sum((gr_current .- gr_old).^2) / length(gr_current)
    return error, iteration_diff
end

function f_over_r_from_potential(potential::Vector{Float64}, r_low::Float64, bin_width::Float64)::Vector{Float64}
    """Computes force/r from the potential."""
    pot_length = length(potential)
    f_over_r = similar(potential)

    @inbounds for i in 2:(pot_length - 1)
        r = r_low + (i - 1) * bin_width
        f_over_r[i] = - (potential[i + 1] - potential[i - 1]) / (2.0 * bin_width * r)
    end

    # Use linear extrapolation for boundary values
    f_over_r[1] = 2.0 * f_over_r[2] - f_over_r[3]
    f_over_r[end] = 2.0 * f_over_r[end - 1] - f_over_r[end - 2]

    return f_over_r
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end