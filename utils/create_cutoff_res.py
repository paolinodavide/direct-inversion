import numpy as np
import os
import json

def main():
    prefix = "pot"
    try:
        largest_file = max(
            (f for f in os.listdir('.') if f.startswith(prefix) and f.endswith('.dat') and f != f"{prefix}_target.dat"),
            key=lambda x: int(x[len(prefix) + 1:-4]),
            default=None
        )
        if largest_file:
            r, pot, force = np.loadtxt(largest_file, unpack=True)
            #print(f"Loaded data from {largest_file}")
        else:
            print("No matching .dat files found.")
    except ValueError:
        print("Error processing file names.")
    
    target_file = f"{prefix}_target.dat"
    try:
        r, target_pot, target_force = np.loadtxt(target_file, unpack=True)
        #print(f"Loaded target data from {target_file}")
    except OSError:
        print(f"File {target_file} not found. Skipping.")
    
    L2norm_pot = np.linalg.norm(pot - target_pot)
    L2norm_force = np.linalg.norm(force - target_force)
    SNR_pot =  10 * np.log10(np.linalg.norm(target_pot) / L2norm_pot)
    SNR_force = 10 * np.log10(np.linalg.norm(target_force) / L2norm_force)
    
    # Load parameters from config.json
    config_file = "config.json"
    try:
        with open(config_file, 'r') as f:
            params = json.load(f)
        #print(f"Loaded parameters from {config_file}")
    except (OSError, json.JSONDecodeError) as e:
        print(f"Error loading {config_file}: {e}")
        return

    output_file = f"cutoff_results.dat"
    try:
        if os.path.exists(output_file):
            mode = 'a'  # Append mode if file exists
        else:
            mode = 'w'  # Write mode if file does not exist
        
        with open(output_file, mode) as f:
            if mode == 'w':
                f.write(f"# rcut\tn_it\tL2_pot\tSNR_pot\tL2_force\tSNR_force\ttarget_file\tinit_pot\n")
            largest_nb = int(largest_file[len(prefix) + 1:-4]) if largest_file else 0
            f.write(f"{(params['x_cut']):.3f}\t{int(largest_nb)}\t{L2norm_pot:.4e}\t{SNR_pot:.2f}\t{L2norm_force:.4e}\t{SNR_force:.2f}\t # {params['target_file']}\t{params['init_pot']}\n")
        print(f"Summary updated in {output_file}")
    except KeyError as e:
        print(f"Error: Missing key in parameters: {e}")
    except OSError as e:
        print(f"Error writing to {output_file}: {e}")

if __name__ == "__main__":
    main()