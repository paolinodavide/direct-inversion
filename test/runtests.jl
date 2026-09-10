# test/runtests.jl
using JSON
using DelimitedFiles
using LinearAlgebra
using StaticArrays
using Random

include(joinpath(@__DIR__, "..", "gr_borgis.jl"))

const TEST_DIR = joinpath(@__DIR__, "test_env")
const BASELINE_DIR = joinpath(@__DIR__, "baseline")

function setup_test_environment()
    rm(TEST_DIR, force=true, recursive=true)
    mkpath(TEST_DIR)
    mkpath(joinpath(TEST_DIR, "inputs"))
    mkpath(joinpath(TEST_DIR, "inputs", "configs"))
    
    # 1. Parse the first 2 snapshots from examples/lj_92_2.dat
    examples_file = joinpath(@__DIR__, "..", "examples", "lj_92_2.dat")
    if !isfile(examples_file)
        error("examples/lj_92_2.dat not found!")
    end
    
    box_length = 30.0
    num_snapshots_to_extract = 2
    snapshots_extracted = 0
    extracted_timesteps = String[]
    
    open(examples_file) do io
        lines = readlines(io)
        i = 1
        while i <= length(lines)
            if startswith(lines[i], "ITEM: TIMESTEP")
                timestep = strip(lines[i+1])
                push!(extracted_timesteps, timestep)
                
                # Find number of atoms
                while !startswith(lines[i], "ITEM: NUMBER OF ATOMS")
                    i += 1
                end
                N = parse(Int, lines[i+1])
                
                # Find atoms start
                while !startswith(lines[i], "ITEM: ATOMS")
                    i += 1
                end
                
                # Write config file config_<timestep>.dat
                filepath = joinpath(TEST_DIR, "inputs", "configs", "config_$(timestep).dat")
                open(filepath, "w") do out_io
                    write(out_io, "# x \t y \t l_box = $(box_length)\n")
                    for _ in 1:N
                        i += 1
                        parts = split(strip(lines[i]))
                        # xs and ys are at index 3 and 4 (scaled coordinates)
                        x = parse(Float64, parts[3]) * box_length
                        y = parse(Float64, parts[4]) * box_length
                        write(out_io, "$(x)\t$(y)\n")
                    end
                end
                
                snapshots_extracted += 1
                if snapshots_extracted >= num_snapshots_to_extract
                    break
                end
            end
            i += 1
        end
    end
    
    if snapshots_extracted < num_snapshots_to_extract
        error("Failed to extract $num_snapshots_to_extract snapshots from $examples_file")
    end
    
    N_particles = 841

    # 2. Generate ordered_wt.dat containing the extracted timesteps
    open(joinpath(TEST_DIR, "inputs", "ordered_wt.dat"), "w") do io
        write(io, "# wt ordered by min dist\n")
        for ts in extracted_timesteps
            write(io, "$(ts)\n")
        end
    end

    # 3. Generate a dummy target g(r) file
    bin_width = 0.01
    r_values = collect(bin_width:bin_width:10.0)
    target_gr = [r < 0.8 ? 0.0 : 1.0 for r in r_values]
    
    open(joinpath(TEST_DIR, "inputs", "gr_weighted.dat"), "w") do io
        write(io, "# r\tg(r)\n")
        writedlm(io, hcat(r_values, target_gr))
    end

    # 4. Generate params.json
    params = Dict(
        "N_particles" => N_particles,
        "n_inversion_snapshots" => num_snapshots_to_extract,
        "L_box" => box_length,
        "dimensions" => 2,
        "bin_width" => bin_width,
        "config_dir" => "configs",
        "wt_file" => "ordered_wt.dat",
        "r_high" => 5.0,
        "x_min" => 0.8,
        "r_low" => 0.8,
        "r_max" => 10.0,
        "target_gr_file" => "gr_weighted.dat",
        "target_precision" => 1e-15,
        "iteration_precision" => 1e-15,
        "output_file" => "gr_final.dat",
        "max_iter" => 3,
        "method_force_formula" => "out",
        "Temperature" => 2.0,
        "init_pot_type" => "wca",
        "target_pot_type" => "lj_full",
        "learning_rate" => 0.1,
        "core_strength" => 14,
        "shift_gr" => true
    )
    
    open(joinpath(TEST_DIR, "inputs", "params.json"), "w") do io
        JSON.print(io, params, 2)
    end
end

function run_inversion()
    # Execute the inversion loop as a subprocess to ensure clean state
    cmd = `julia --project=$(joinpath(@__DIR__, "..")) $(joinpath(@__DIR__, "..", "grinter_parallel.jl")) -d $(TEST_DIR)`
    run(cmd)
end

function generate_baseline()
    println("Setting up test environment...")
    setup_test_environment()
    
    println("Running inversion to generate baseline data...")
    run_inversion()
    
    println("Saving outputs to baseline directory...")
    rm(BASELINE_DIR, force=true, recursive=true)
    mkpath(BASELINE_DIR)
    
    # Copy output files we want to preserve
    outputs_to_copy = ["gr_target.dat", "convergence_data.dat"]
    for file in readdir(joinpath(TEST_DIR, "outputs"))
        if endswith(file, ".dat")
            cp(joinpath(TEST_DIR, "outputs", file), joinpath(BASELINE_DIR, file))
        end
    end
    println("Baseline successfully generated in $(BASELINE_DIR)!")
end

function compare_with_baseline()
    println("Setting up test environment...")
    setup_test_environment()
    
    println("Running inversion to verify against baseline...")
    run_inversion()
    
    println("Comparing outputs...")
    
    baseline_files = filter(f -> endswith(f, ".dat"), readdir(BASELINE_DIR))
    
    mismatch_found = false
    for file in baseline_files
        baseline_path = joinpath(BASELINE_DIR, file)
        current_path = joinpath(TEST_DIR, "outputs", file)
        
        if !isfile(current_path)
            println("❌ Error: Output file $(file) was not generated!")
            mismatch_found = true
            continue
        end
        
        # Read the numerical values and compare them
        baseline_data = readdlm(baseline_path, comments=true)
        current_data = readdlm(current_path, comments=true)
        
        if file == "convergence_data.dat"
            # Ignore the elapsed time column (second column) since it varies between runs
            baseline_data = baseline_data[:, [1, 3, 4, 5, 6, 7]]
            current_data = current_data[:, [1, 3, 4, 5, 6, 7]]
        end
        
        if size(baseline_data) != size(current_data)
            println("❌ Error: Size mismatch in $(file)!")
            mismatch_found = true
            continue
        end
        
        # Compare within a tight numerical tolerance
        max_diff = maximum(abs.(baseline_data .- current_data))
        if max_diff > 1e-12
            println("❌ Error: Numerical difference in $(file) exceeds tolerance! Max diff = $(max_diff)")
            mismatch_found = true
        else
            println("✅ File $(file) matches baseline exactly (max diff = $(max_diff))")
        end
    end
    
    # 2. Run Borgis vs Histogram verification unit test
    if !test_borgis_vs_histogram()
        mismatch_found = true
    end
    
    if mismatch_found
        println("❌ Regression tests FAILED!")
        exit(1)
    else
        println("🎉 All regression tests PASSED!")
    end
end

# Simple Lennard-Jones force divided by r helper
function lennard_jones_force_div_r(r, epsilon=1.0, sigma=1.0)
    if r == 0.0
        return 0.0
    end
    sr = sigma / r
    sr8 = sr^8
    sr14 = sr^14
    return (24.0 * epsilon / (sigma^2)) * (2.0 * sr14 - sr8)
end

function test_borgis_vs_histogram()
    println("\nRunning Borgis vs Histogram verification unit test...")
    
    # 1. Parse the first snapshot from examples/lj_92_2.dat (thermodynamic equilibrium configuration)
    examples_file = joinpath(@__DIR__, "..", "examples", "lj_92_2.dat")
    if !isfile(examples_file)
        println("❌ Error: examples/lj_92_2.dat not found!")
        return false
    end
    
    box_length = 30.0
    positions = SVector{2, Float64}[]
    
    open(examples_file) do io
        lines = readlines(io)
        atoms_start = findfirst(l -> startswith(l, "ITEM: ATOMS"), lines)
        if atoms_start === nothing
            error("Invalid LAMMPS dump format in examples/lj_92_2.dat")
        end
        
        for i in (atoms_start + 1):length(lines)
            line = strip(lines[i])
            if isempty(line) || startswith(line, "ITEM:")
                break
            end
            parts = split(line)
            # xs and ys are at index 3 and 4
            x = parse(Float64, parts[3]) * box_length
            y = parse(Float64, parts[4]) * box_length
            push!(positions, SVector(x, y))
        end
    end
    
    N = length(positions)
    dim = 2
    println("  Loaded $N equilibrium particles from examples/lj_92_2.dat")
    
    # 2. Compute Histogram g(r)
    bin_width = 0.02
    max_distance = 5.0
    num_bins = floor(Int, max_distance / bin_width)
    counts = zeros(Int, num_bins)
    box_sizes = fill(box_length, SVector{dim, Float64})
    
    for i in 1:N-1
        pos_i = positions[i]
        for j in i+1:N
            pos_j = positions[j]
            rVec, r2 = pbc_distance(pos_i, pos_j, box_sizes)
            r = sqrt(r2)
            if r < max_distance
                bin_idx = floor(Int, r / bin_width) + 1
                if bin_idx <= num_bins
                    counts[bin_idx] += 1
                end
            end
        end
    end
    
    rho = N / (box_length^dim)
    r_vals = [(i - 0.5) * bin_width for i in 1:num_bins]
    gr_histogram = zeros(num_bins)
    for i in 1:num_bins
        r_inner = (i - 1) * bin_width
        r_outer = i * bin_width
        shell_volume = π * (r_outer^2 - r_inner^2)
        gr_histogram[i] = (counts[i] * 2.0) / (N * rho * shell_volume)
    end
    
    # 3. Compute Borgis g(r)
    # Precompute forces for LJ
    forces = zeros(SVector{dim, Float64}, N)
    r_cut = 3.5
    T = 2.0 # Temperature
    for i in 1:N-1
        pos_i = positions[i]
        for j in i+1:N
            pos_j = positions[j]
            rVec, r2 = pbc_distance(pos_i, pos_j, box_sizes)
            r = sqrt(r2)
            if r < r_cut && r > 0.0
                f_mag = lennard_jones_force_div_r(r, 1.0, 1.0) / T
                forces[i] += f_mag * rVec
                forces[j] -= f_mag * rVec
            end
        end
    end
    
    # Include gr_borgis.jl functions indirectly
    # Already included at top of file
    
    borgis_contributions = zeros(Float64, num_bins)
    max_dist = num_bins * bin_width
    compute_borgis_contributions!(borgis_contributions, positions, forces, box_length, 1.0/bin_width, num_bins)
    integrate_borgis_contributions(borgis_contributions, Out())
    prefactor = compute_prefactor(N, box_length, dim)
    gr_borgis = 1.0 .- borgis_contributions ./ prefactor
    
    # 4. Compare the two curves for r in [0.85, 4.5]
    valid_indices = findall(r -> r > 0.85 && r < 4.5, r_vals)
    hist_sub = gr_histogram[valid_indices]
    borgis_sub = gr_borgis[valid_indices]
    
    mae = sum(abs.(hist_sub .- borgis_sub)) / length(valid_indices)
    mean_hist = sum(hist_sub) / length(valid_indices)
    mean_borgis = sum(borgis_sub) / length(valid_indices)
    
    cov = sum((hist_sub .- mean_hist) .* (borgis_sub .- mean_borgis))
    var_hist = sum((hist_sub .- mean_hist).^2)
    var_borgis = sum((borgis_sub .- mean_borgis).^2)
    correlation = cov / sqrt(var_hist * var_borgis)
    
    println("  Pearson Correlation:      $(round(correlation, digits=4))")
    println("  Mean Absolute Error (MAE): $(round(mae, digits=4))")
    
    if correlation > 0.85 && mae < 0.20
        println("✅ Borgis vs Histogram verification PASSED!")
        return true
    else
        println("❌ Borgis vs Histogram verification FAILED! Correlation: $correlation, MAE: $mae")
        return false
    end
end

function main()
    if "--generate-baseline" in ARGS
        generate_baseline()
    else
        if !isdir(BASELINE_DIR)
            println("No baseline found. Generating baseline first...")
            generate_baseline()
        else
            compare_with_baseline()
        end
    end
end

main()
