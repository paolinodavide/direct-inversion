import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_smoothing_spline
import sys
from scipy.signal import find_peaks

def get_weights(r, gr):
    dr = r[1] - r[0]
    # Mask to filter out non-positive values of g(r)
    mask = gr > 0
    g_masked = gr[mask]
    r_masked = r[mask]

    # Find indices where g_masked crosses 1
    crossings = np.where(np.diff(np.sign(g_masked - 1)))[0]

    # Filter crossings to ensure at least 0.1 distance between them
    filtered_crossings = [crossings[0]] if len(crossings) > 0 else []
    for i in range(1, len(crossings)):
        if r_masked[crossings[i]] - r_masked[filtered_crossings[-1]] >= 0.1:
            filtered_crossings.append(crossings[i])
    crossings = np.array(filtered_crossings)

    # Check if there is a second crossing
    if len(crossings) > 1:
        second_crossing_index = crossings[1]
        second_crossing_r = r_masked[second_crossing_index]
        print(f"The second crossing of g_masked with 1 occurs at r = {second_crossing_r}")
    else:
        print("There is no second crossing of g_masked with 1.")

    # Find the moment where g_masked oscillates around 1 with an amplitude of 0.02
    oscillation_start_index = None
    for i in range(len(g_masked)):
        if np.all(np.abs(g_masked[i:] - 1) <= 10 * dr):  # Check the next 10 points
            oscillation_start_index = i
            break

    if oscillation_start_index is not None:
        oscillation_start_r = r_masked[oscillation_start_index]
        print(f"g_masked starts oscillating around 1 with an amplitude of 10*dr at r = {oscillation_start_r}")
    else:
        print("g_masked does not oscillate around 1 with an amplitude of 10*dr within the given range.")

    weigths = np.ones_like(r)
    if len(crossings) > 1:
        weigths[r < r_masked[crossings[1]]] *= 10**3
    if oscillation_start_index is not None:
        weigths[r > oscillation_start_r] *= 10**-2

    # plt.figure()
    # plt.plot(r, weigths, label='Weights')
    # plt.xlabel('r')
    # plt.ylabel('Weight')
    # plt.legend()
    # plt.grid('--')
    # plt.show()
    precision = dr
    print(f"Lambda set to {precision:e}")
    return weigths, dr


def main():
    r, g_r = np.loadtxt('./rdfs/g_r_h_avg.dat', unpack=True)
    plt.plot(r, g_r, 'o', label='Noisy RDF', markersize=3)

    # Load target RDF if available
    try:
        r, gr_tgt_old = np.loadtxt('./gs_target.dat', unpack=True)
        plt.plot(r, gr_tgt_old, 'o', label='Target RDF', markersize=3)
    except OSError:
        print("File 'gs_target.dat' not found. ")
        
    # Import three integers from sys.argv or set preset values
    try:
        low_exp, high_exp, lambda_exp = map(int, sys.argv[1:4])
        print(f"Using command line arguments: low_exp={low_exp}, high_exp={high_exp}, lambda_exp={lambda_exp}")
        weights = np.ones_like(r)
        weights[r < 1.25] = 10**low_exp
        weights[r > 5] = 10**high_exp
        precision = 10**lambda_exp
    except (ValueError, IndexError):
        weights, precision = get_weights(r, g_r)
        
        
    # Weighted spline
    spline = make_smoothing_spline(r, g_r, w=weights, lam=precision)  # Adjust 's' as needed

    # Plot
    plt.plot(r, spline(r), '--', label='Weighted Spline')
    plt.xlabel('r')
    plt.ylabel('g(r)')
    plt.legend()
    plt.grid('--')
    plt.show()

    output_file = 'gr_weighted.dat'
    np.savetxt(output_file, np.column_stack((r, spline(r))), delimiter='\t', header="# r\tg(r)", comments='')
    print(f"Weighted spline saved to {output_file}.")

    return

if __name__ == "__main__":
    main()
