import json
import os

# Create the inputs directory if it doesn't exist
os.makedirs("inputs", exist_ok=True)

# Define the dummy parameters
params = {
    "N_particles": 842,
    "n_inversion_snapshots": 500,
    "L_box": 30.0,
    "dimensions": 2,
    "bin_width": 0.002,
    "config_dir": "configs_bin",
    "wt_file": "ordered_wt.dat",
    "r_high": 2.5,
    "x_min": 0.881,
    "r_low": 0.86,
    "target_gr_file": "gr_weighted.dat",
    "target_precision": 1e-6,
    "iteration_precision": 1e-9,
    "output_file": "gr_final.dat",
    "max_iter": 500,
    "method_force_formula": "out",
    "Temperature": 2.0,
    "init_pot_type": "mean_force",
    "target_pot_type": "lj_full",
    "learning_rate": 0.2,
    "core_strength": 14,
    "shift_gr": True
}

# Write to JSON file
if os.path.exists("inputs/params.json"):
    print("File inputs/params.json already exists. Stopping.")
else:
    with open("inputs/params.json", "w") as f:
        json.dump(params, f, indent=2)
    print("Created inputs/params.json")