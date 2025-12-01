import os
import json

# Parameters to be written to the JSON file
params = {
    "N_particles": 1000,
    "n_inversion_snapshots": 125,
    "L_box": 60.0,
    "dimensions": 2,
    "bin_width": 0.002,
    "config_dir": "configs_bin",
    "wt_file": "ordered_wt.dat",
    "r_high": 2.5,
    "x_min": 0.5859,
    "r_low": 0.65,
    "target_gr_file": "gr_weighted.dat",
    "target_precision": 0.0,
    "iteration_precision": 1e-9,
    "output_file": "gr_final.dat",
    "max_iter": 50,
    "method_force_formula": "out",
    "Temperature": 1.0,
    "init_pot_type": "mean_force",
    "target_pot_type": "lj_full",
    "learning_rate": 0.2
}

# Directory and file paths
inputs_dir = "./inputs"
params_file = os.path.join(inputs_dir, "params.json")

# Create the inputs directory if it does not exist
os.makedirs(inputs_dir, exist_ok=True)

# Write the parameters to the JSON file
with open(params_file, "w") as f:
    json.dump(params, f, indent=4)

print(f"Initialized params.json in {inputs_dir}")