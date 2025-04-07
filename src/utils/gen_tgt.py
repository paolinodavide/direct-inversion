import numpy as np
import json
import sys
import matplotlib.pyplot as plt
import matplotlib.mathtext
from scipy.interpolate import make_smoothing_spline

# Constants - named in UPPER_CASE for better identification
INPUT_FILENAME = '../rdfs_0.002/g_r_h_avg.dat'
OUTPUT_HISTO_JSON = 'histo.json'
OUTPUT_SMOOTH_JSON = 'smooth.json'

# Window parameters - grouped together for clarity
WINDOW1_MIN = 0.80
WINDOW1_MAX = 1.084
WINDOW2_MIN = 1.0
WINDOW2_MAX = 4.0

OVERLAP_CENTER = 1.024
R_MAX = 10
DR = 0.002

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


################### Figure settings ######################################
# # good latex font
# from matplotlib import rc
# rc('text', usetex = True)
# #style
# plt.style.use('figure')
# # put figure cleverly
# plt.rcParams.update({'figure.autolayout': True})


def main():
    ################# Save or not save #######################################
    try:
        sys.argv[1]
    except IndexError:
        sys.argv.append(False)

    if (sys.argv[1] == 'S'):
        save = True
    else:
        save = False


    # Step 1: Load and save target data
    target_data = load_target_data(INPUT_FILENAME)
    with open(OUTPUT_HISTO_JSON, 'w') as f:
        json.dump(target_data, f)
    
    # Step 2: Prepare data for interpolation
    with open(OUTPUT_HISTO_JSON) as f:
        target = json.load(f)
    
    positions = list(target.keys())
    g_values = list(target.values())
    
    # Step 3: Create interpolation windows
    pos_window1, pos_window2 = create_window_arrays()
    
    # Step 4: Filter data for each window
    pos_red1, g_red1 = filter_data_by_window(positions, g_values, WINDOW1_MIN, WINDOW1_MAX)
    pos_red2, g_red2 = filter_data_by_window(positions, g_values, WINDOW2_MIN, WINDOW2_MAX)
    

    spline1 = make_smoothing_spline(pos_red1, g_red1, lam=(1-0.999999))
    spline2 = make_smoothing_spline(pos_red2, g_red2, lam=(1-0.99995))
    
    # Step 6: Merge windows
    merged_pos, merged_val = merge_windows(pos_window1, spline1(pos_window1), pos_window2, spline2(pos_window2), OVERLAP_CENTER)

    # Step 7: Add short and long range behavior
    smoothed_dict = create_smoothed_gr(merged_val, R_MAX, DR)
    
    # Step 8: Save results
    with open(OUTPUT_SMOOTH_JSON, 'w') as f:
        json.dump(smoothed_dict, f)

    ###################### Figure ###############################################
    #
    fig, ax = plt.subplots()
    #
    ax.plot(pos_red1, g_red1, color = 'blue', marker = 'o', markersize = '5', \
            linestyle = 'none', fillstyle = 'none', label = r'$g_{histo}$')
    ax.plot(pos_red2, g_red2, color = 'blue', marker = 'o', markersize = '5', \
            linestyle = 'none', fillstyle = 'none')
    #ax.plot(xi1, g_test1, color = 'red', marker = 'o', markersize = '5', linestyle = 'none', label = r'$g_{smooth}$')
    #ax.plot(xi2, g_test2, color = 'red', marker = 'o', markersize = '5', linestyle = 'none')
    ax.plot(merged_pos, merged_val, color = 'red', label = r'$g_{smooth}$')
    # Add vertical lines for overlap center and window limits
    ax.axvline(x=OVERLAP_CENTER, linestyle='--', linewidth=0.5, label='Overlap Center')
    ax.axvline(x=WINDOW1_MAX, linestyle='--', linewidth=0.5, label='Window 1 Upper Limit')
    ax.axvline(x=WINDOW2_MIN, linestyle='--', linewidth=0.5, label='Window 2 Lower Limit')

    ax.set_xlabel(r'$r$')
    ax.set_ylabel(r'$g(r)$')
    ax.set_title(r'$g(r)$')
    ax.grid(True)
    ax.legend()

    plt.subplots_adjust()

    # see or save
    if (save == True):
        plt.savefig('TestCsaps2.pdf', transparent = False)
    else:
        plt.show()

if __name__ == '__main__':
    main()