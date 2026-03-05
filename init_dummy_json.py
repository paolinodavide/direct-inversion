import json
import os
import argparse

parser = argparse.ArgumentParser(description="Initialize dummy params.json for forceIBI")
parser.add_argument(
    "--directory", "-d",
    type=str, 
    required=True,  # This makes the flag mandatory
    help="The required home directory prefix where 'inputs' will be created"
)
args = parser.parse_args()

# Create the inputs directory if it doesn't exist
os.makedirs(os.path.join(args.directory, "inputs"), exist_ok=True)

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
if os.path.exists(os.path.join(args.directory, "inputs/params.json")):
    print("File inputs/params.json already exists. Stopping.")
else:
    with open(os.path.join(args.directory, "inputs/params.json"), "w") as f:
        json.dump(params, f, indent=2)
    print(f"Created {os.path.join(args.directory, 'inputs/params.json')} with dummy parameters.")