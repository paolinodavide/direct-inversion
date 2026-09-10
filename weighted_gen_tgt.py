import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_smoothing_spline
import sys
import argparse
import os

def get_weights(r, gr):
    dr = r[1] - r[0]
    mask = gr > 0
    g_masked, r_masked = gr[mask], r[mask]

    max_g = np.max(g_masked)
    crossings = np.where(np.diff(np.sign(g_masked - 1)))[0]

    # Filter crossings to ensure at least 0.1 distance between them
    filtered_crossings = [crossings[0]] if len(crossings) > 0 else []
    for i in range(1, len(crossings)):
        if r_masked[crossings[i]] - r_masked[filtered_crossings[-1]] >= 0.1:
            filtered_crossings.append(crossings[i])
    crossings = np.array(filtered_crossings)

    if len(crossings) >= 1:
        print(f"First crossing with gMax/2 = {max_g / 2:.2f} at r = {r_masked[crossings[0]]}")
    else:
        print(f"No crossing of g_masked with {max_g / 2:.2f}.")

    threshold = max_g / 10
    oscillation_start_index = next(
        (i for i in range(len(g_masked)) if np.all(np.abs(g_masked[i:] - 1) <= threshold)), None
    )
    oscillation_start_index = round(len(g_masked) * 0.75)
    if oscillation_start_index is not None:
        print(f"g(r) starts oscillating with amplitude = {threshold:.2f} at r = {r_masked[oscillation_start_index]}")
    else:
        print("g_masked does not oscillate around 1 within the given range.")


    weights = np.ones_like(r)
    if len(crossings) >= 1:
        weights[r < r_masked[crossings[0]]] *= 10**2
    if oscillation_start_index is not None:
        weights[r > r_masked[oscillation_start_index]] *= 10**1
        weights[ (r_masked[crossings[0]] < r) & (r < r_masked[oscillation_start_index]) ] *= 10**0

    return weights

def main():
    parser = argparse.ArgumentParser(description="Generate weighted target g(r) for forceIBI")
    parser.add_argument(
        "--directory", "-d",
        type=str,
        required=True,
        help="The required home directory prefix where 'inputs' and 'outputs' will be created"
    )
    parser.add_argument(
        "--dr",
        type=float,
        default=0.002,
        help="New bin width for output g(r) grid (default: 0.002)"
    )
    args = parser.parse_args()

    rdf_path = os.path.join(args.directory, "outputs/rdfs/")
    r, g_r, var_g = np.loadtxt(os.path.join(rdf_path, 'g_r_h_avg.dat'), unpack=True)
    plt.plot(r, g_r, 'o', label='Histo RDF', markersize=3)

    try:
        r_old, gr_tgt_old = np.loadtxt('./gs_target.dat', unpack=True)
        plt.plot(r_old, gr_tgt_old, 'o', label='Target RDF', markersize=3)
    except OSError:
        print("File 'gs_target.dat' not found.")

    try:
        low_exp, high_exp = map(int, sys.argv[1:4])
        print(f"Using command line arguments: low_exp={low_exp}, high_exp={high_exp}")
        weights = np.ones_like(r)
        weights[r < 1.25] = 10**low_exp
        weights[r > 5] = 10**high_exp
    except (ValueError, IndexError):
        weights = get_weights(r, g_r)

    spline = make_smoothing_spline(r, g_r, weights)

    # 1. Obtain old dr
    dr_old = r[1] - r[0]
    # 2. Calculate old bin edges to evaluate min_r and max_r
    min_r = r[0] - dr_old / 2.0
    max_r = r[-1] + dr_old / 2.0
    # 3. Produce new bin centers with new dr which become new radii
    dr_new = args.dr
    new_radii = np.arange(min_r + dr_new / 2.0, max_r, dr_new)

    g_smoothed = spline(new_radii)
    if (g_min := np.min(g_smoothed)) < -1e-5:
        print(f"\nWarning: g_smoothed has a minimum value of {g_min}, below the threshold.")

    plt.plot(new_radii, g_smoothed, label=f'Spline g({new_radii[-1]:.1f}) = {g_smoothed[-1]:.3f}')
    plt.xlabel(r'$r/\sigma$')
    plt.ylabel('g(r)')
    plt.legend(loc='lower right', frameon=False)
    plt.tick_params(direction="in", top=True, right=True)
    plt.xlim(0, max_r)
    plt.tight_layout()
    plt.savefig(os.path.join(rdf_path, '01_spline_gr.pdf'), dpi=300)
    plt.savefig(os.path.join(rdf_path, '01_spline_gr.png'))
    plt.show()

    output_file = os.path.join(args.directory, "inputs/gr_weighted.dat")
    np.savetxt(output_file, np.column_stack((new_radii, g_smoothed)), delimiter='\t', header="# r\tg(r)", comments='')
    print(f"Weighted spline saved to {output_file}.")

if __name__ == "__main__":
    main()
