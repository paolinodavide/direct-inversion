using LinearAlgebra

"""
Module containing functions related to the potential part of the iteration program.
"""

function lj_x(x::Float64, T::Float64)::Float64
    """Force derived from the Lennard-Jones potential/x."""
    return 2.0 / T * (12.0 * x^(-14) - 6.0 * x^(-8))
end

function lj_full(x::Float64, T::Float64)::Float64
    """Lennard-Jones potential."""
    return 4.0 / T * (x^(-12) - x^(-6))
end

function lj_rep(x::Float64, T::Float64)::Float64
    """Repulsive part of the Lennard-Jones potential."""
    return 4.0 / T * (x^(-12))
end

function lj_att(x::Float64, T::Float64)::Float64
    """Attractive part of the Lennard-Jones potential."""
    return 4.0 / T * (-x^(-6))
end

function lj_test(x::Float64, T::Float64)::Float64
    """Test potential."""
    return 4.0 / T * (x^(-12.00001) - x^(-6))
end

function shoulder_potential(x::Float64, T::Float64)::Float64
    """Shoulder potential."""
    pot = 1/x^14 - 0.5*tanh(10 * (x - 2.5))  
    return pot / T
end

function wca_potential(x::Float64, T::Float64)::Float64
    """Weeks-Chandler-Andersen potential."""
    sigma = 1.0
    epsilon = 1.0 

    if x < 2.0^(1/6) * sigma
        return 4.0 * epsilon / T * ( (sigma / x)^12 - (sigma / x)^6 ) + epsilon / T
    else
        return 0.0
    end
end

function get_potential_from_name(name::String, T::Float64, r_low::Float64, r_high::Float64, bin_width)::Vector{Float64}
    """Choose the potential to initiate the iteration loop. Returns the potential as an array."""
    radii = collect(r_low:bin_width:r_high)
    pot_list = zeros(Float64, length(radii))

    potential_function = Dict(
        "lj_full" => lj_full,
        "lj_rep" => lj_rep,
        "lj_att" => lj_att,
        "lj_test" => lj_test,
        "wca" => lj_full,
        "r3" => (x, T) -> 1.0 / T * x^(-3), 
        "sh" => shoulder_potential
    )

    if name == "zero"
        return pot_list
    elseif haskey(potential_function, name)
        func = potential_function[name]
        for (i, x) in enumerate(radii)
            pot_list[i] = func(x, T)
        end
        pot_list .-= pot_list[end]  # Offset adjustment
        return pot_list
    else
        println("Unknown potential name: ", name)
        return pot_list
    end
end

### IO Functions ###

function initialize_convergence_file(filename)
    open(filename, "w") do io
        write(io, "# iteration\ttime\terror\titeration_difference\tpotential_increase\tdelta_tgt\tphi\n")
    end
end

function append_convergence_data(filename, iteration, elapsed_time, error, iteration_diff, potential_increase, delta_target, phi) 
    open(filename, "a") do io
        writedlm(io, [[iteration elapsed_time error iteration_diff potential_increase delta_target phi]])
    end
end

function save_iteration_data(iteration, r_range, gr_current, βu_current, f_current, home_dir::String)
    output_data = hcat(collect(r_range), gr_current, βu_current, f_current)
    open(joinpath(home_dir, "outputs/iteration_$iteration.dat"), "w") do io
        write(io, "# r\tgr\tu\tf/r\n")
        writedlm(io, output_data)
    end
end

function save_target_data(r_range, gr_target, βu_target, f_target, filename::String = "outputs/iteration_-1.dat")
    output_data = hcat(collect(r_range), gr_target, βu_target, f_target)
    open(filename, "w") do io
        write(io, "# r\tgr\tu\tf/r\n")
        writedlm(io, output_data)
    end
end

function save_gr_data(r_range, gr_data, filename)
    output_data = hcat(collect(r_range), gr_data)
    open(filename, "w") do io
        write(io, "# r\tgr\n")
        writedlm(io, output_data)
    end
end

# Helper functions for better organization
function ensure_binary_configs(config_dir, wt_file_path, number_config)
    if !any(endswith.(readdir(config_dir), ".bin"))
        @info "Converting configuration files to binary format"
        bin_dir = config_dir * "_bin"
        convert_to_binary(config_dir, bin_dir)
        config_dir = bin_dir
    end
    
    # Filter and copy selected weight files
    inversion_dir = joinpath(config_dir, "../inversion_bin")
    mkpath(inversion_dir)
    
    wt_files = readdir(config_dir)
    order_wt = readdlm(wt_file_path, comments=true)[:, 1]
    selected_wt = Set(order_wt[1:number_config])
    
    for file in wt_files
        if endswith(file, ".bin")
            file_weight = parse(Float64, split(basename(file), ['_', '.'])[2])
            if file_weight in selected_wt
                dest_file = joinpath(inversion_dir, file)
                !isfile(dest_file) && cp(joinpath(config_dir, file), dest_file)
            end
        end
    end
    
    return inversion_dir
end