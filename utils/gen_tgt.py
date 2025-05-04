import numpy as np
import json
import matplotlib.pyplot as plt
from scipy.interpolate import make_smoothing_spline
from scipy.signal import savgol_filter

# Constants - named in UPPER_CASE for better identification
INPUT_FILENAME = 'rdfs/g_r_h_avg.dat'
OUTPUT_HISTO_JSON = './src/histo.json'
OUTPUT_SMOOTH_JSON = './src/smooth.json'

# Window parameters - grouped together for clarity
WINDOW1_MIN = 0.80
WINDOW1_MAX = 1.15

WINDOW2_MIN = 1.02
WINDOW2_MAX = 5

OVERLAP_CENTER = 1.124
R_MAX = 10
DR = 0.002

p1 = 0.999999
p2 = 0.9999

def load_target_data(filename) -> dict:
    """Load target data from file and return as dictionary."""

    target_data = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                target_data[parts[0]] = parts[1]
    return target_data

def create_window_arrays():
    """Create position arrays for interpolation windows."""
    
    window1_positions = np.arange(WINDOW1_MIN, WINDOW1_MAX + DR, DR)
    window2_positions = np.arange(WINDOW2_MIN, WINDOW2_MAX + DR, DR)
    return window1_positions, window2_positions

def filter_data_by_window(positions, values, window_min, window_max):
    """Filter data points within specified window."""

    filtered_pos = []
    filtered_val = []
    for pos, val in zip(positions, values):
        pos_float = float(pos)
        if window_min <= pos_float <= window_max:
            filtered_pos.append(pos_float)
            filtered_val.append(float(val))
    return filtered_pos, filtered_val

def merge_windows(pos1, val1, pos2, val2, overlap_center):
    """Merge data from two windows with overlap handling."""

    merged_pos = []
    merged_val = []
    
    # Add points from first window up to overlap center
    for p, v in zip(pos1, val1):
        if p < overlap_center:
            merged_pos.append(p)
            merged_val.append(v)
        elif p == overlap_center:
            # Find matching point in second window
            try:
                j = pos2.index(overlap_center)
                merged_pos.append(p)
                merged_val.append(0.5 * v + 0.5 * val2[j])
            except ValueError:
                continue
    
    # Add points from second window after overlap center
    for p, v in zip(pos2, val2):
        if p > overlap_center:
            merged_pos.append(p)
            merged_val.append(v)
    
    return merged_pos, merged_val

def create_smoothed_gr(merged_val, r_max, dr) -> dict:
    """Create final smoothed dictionary for output."""

    all_positions = np.arange(0., r_max, dr)
    index_threshold = int(WINDOW1_MIN / dr)
    result = {}
    
    for i, pos in enumerate(all_positions):
        if pos < WINDOW1_MIN:
            result[pos] = 0
        elif WINDOW1_MIN <= pos < WINDOW2_MAX:
            result[pos] = merged_val[i - index_threshold]
        else:
            result[pos] = 1
    return result

def main():
    # Read the averaged g(r) and r values from the input file
    r_values = []
    g_values = []
    with open(INPUT_FILENAME, 'r') as input_file:
        for line in input_file:
            if line.strip() and not line.startswith('#'):
                r_val, g_r_val = map(float, line.split())
                r_values.append(r_val)
                g_values.append(g_r_val)
    r_values = np.array(r_values)
    g_values = np.array(g_values)

    # Step 3: Create interpolation windows
    pos_window1, pos_window2 = create_window_arrays()

    # Step 4: Filter data for each window
    pos_red1, g_red1 = filter_data_by_window(r_values, g_values, WINDOW1_MIN, WINDOW1_MAX)
    pos_red2, g_red2 = filter_data_by_window(r_values, g_values, WINDOW2_MIN, WINDOW2_MAX)

    spline1 = make_smoothing_spline(pos_red1, g_red1, lam=(1-p1))
    spline2 = make_smoothing_spline(pos_red2, g_red2, lam=(1-p2))
    scipy_pos, scipy_val = merge_windows(pos_window1, spline1(pos_window1), pos_window2, spline2(pos_window2), OVERLAP_CENTER)


    win1_spline = spline1(pos_window1)
    win2_spline = spline2(pos_window2)


    # Plot the RDF
    plt.figure(figsize=(10, 6))

    plt.plot(scipy_pos, scipy_val, color='red', label='Global Spline', linewidth=2)
    plt.plot(pos_window1, win1_spline, color='green', label=f'p1 = {p1}')
    plt.plot(pos_window2, win2_spline, color='purple', label=f'p2 = {p2}')

    # Plot the points used for interpolation (pos_red1 and pos_red2)
    plt.scatter(pos_red1, g_red1, color='blue', s=10)
    plt.scatter(pos_red2, g_red2, color='blue', s=10)


    # Add vertical lines for reference
    plt.axvline(WINDOW1_MAX, color='red', linestyle='--', label=f'{WINDOW1_MAX}', alpha=0.5)
    plt.axvline(WINDOW2_MIN, color='blue', linestyle='--', label=f'{WINDOW2_MIN}', alpha=0.5)
    plt.axvline(OVERLAP_CENTER, color='green', linestyle='--', label=f'{OVERLAP_CENTER}', alpha=0.5)

    # Add labels, legend, and grid
    plt.xlabel("r")
    plt.ylabel("g(r)")
    plt.title("Scipy Method vs Interpolation Points")
    plt.legend()
    plt.grid()

    plt.show()

    final_gr = create_smoothed_gr(scipy_val, R_MAX, DR)
    # Write scipy_pos and scipy_val to the output file
    with open('gs_target.dat', 'w') as output_file:
        output_file.write("# r\tg(r)\n")
        for r, g_r in final_gr.items():
            output_file.write(f"{r:.6f}\t{g_r:.6f}\n")

if __name__ == '__main__':
    main()