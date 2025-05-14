import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_smoothing_spline
import sys

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
    except (ValueError, IndexError):
        low_exp, high_exp, lambda_exp = 3, -2, -3
    print(f"Using low_exp={low_exp}, high_exp={high_exp}, lambda_exp={lambda_exp}")
    
    # Weighted spline
    weights = np.ones_like(r)
    weights[r < 1.25] = 10**low_exp
    weights[r > 5] = 10**high_exp
    spline = make_smoothing_spline(r, g_r, w=weights, lam=10**lambda_exp)  # Adjust 's' as needed

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
