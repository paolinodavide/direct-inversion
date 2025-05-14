import numpy as np
import matplotlib.pyplot as plt

def main():
    r, gr_converged_1 = np.loadtxt("gr_50.dat", unpack=True)
    r, gr_target_1 = np.loadtxt("gr_target.dat", unpack=True)

    r, gr_converged_2 = np.loadtxt("gr_50_weight.dat", unpack=True)
    r, gr_target_2 = np.loadtxt("gr_target_weight.dat", unpack=True)

    l2_norm = np.linalg.norm(gr_converged_1 - gr_converged_2) / len(gr_converged_1)
    print("=== L2 Norms ===")
    print(f"L2 norm between the two converged datasets:\t {l2_norm:.6e}")
    
    l2_norm_target = np.linalg.norm(gr_target_1 - gr_target_2) / len(gr_target_1)
    print(f"L2 norm between the two target datasets:\t {l2_norm_target:.6e}")

    l2_norm_converged_target_1 = np.linalg.norm(gr_converged_1 - gr_target_1) / len(gr_converged_1)
    print(f"L2 norm between gr_converged_1 and gr_target_1:\t {l2_norm_converged_target_1:.6e}")

    l2_norm_converged_target_2 = np.linalg.norm(gr_converged_2 - gr_target_2) / len(gr_converged_2)
    print(f"L2 norm between gr_converged_2 and gr_target_2:\t {l2_norm_converged_target_2:.6e}")
    
    r, gr_h_avg = np.loadtxt("../rdfs/g_r_h_avg.dat", unpack=True)

    l2_norm_h_avg_target_1 = np.linalg.norm(gr_h_avg - gr_target_1) / len(gr_h_avg)
    print(f"L2 norm between gr_h_avg and gr_target_1:\t {l2_norm_h_avg_target_1:.6e}")

    l2_norm_h_avg_target_2 = np.linalg.norm(gr_h_avg - gr_target_2) / len(gr_h_avg)
    print(f"L2 norm between gr_h_avg and gr_target_2:\t {l2_norm_h_avg_target_2:.6e}")
    print("================\n")
    return

if __name__ == "__main__":
    main()