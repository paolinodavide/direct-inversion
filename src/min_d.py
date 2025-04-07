import numpy as np
import sys
from multiprocessing import Pool
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import os

def min_distance_bruteforce(file_number, dim=2):
         # Read positions from file
    filename = f'./ljs/lj_{file_number}000.dat'
    with open(filename, 'r') as f:
        # Read first line to get box size
        first_line = f.readline().strip()
        box_size = float(first_line.split('l_box = ')[1])

        # Read remaining lines for coordinates
        data = []
        for line in f:
            if line.strip():  # skip empty lines
                x, y = map(float, line.strip().split())
                data.append([x, y])
    positions = np.array(data)

    box_size = 60
    min_distance = box_size
    N = len(positions)
    for i in range(N):
        for j in range(i+1, N):
            dx = positions[j,0] - positions[i,0]
            dy = positions[j, 1] - positions[i,1]

            # Apply periodic boundary conditions
            dx -= box_size * round(dx / box_size)
            dy -= box_size * round(dy / box_size)

            distance = np.sqrt(dx**2 + dy**2)
            if distance < min_distance:
                min_distance = distance
    print(f"Minimum distance for file {file_number}: {min_distance}")
    return min_distance

def min_distance(file_number, dim=2):
    # Read positions from file
    filename = f'./ljs/lj_{file_number}000.dat'
    with open(filename, 'r') as f:
        # Read first line to get box size
        first_line = f.readline().strip()
        box_size = float(first_line.split('l_box = ')[1])

        # Read remaining lines for coordinates
        data = []
        for line in f:
            if line.strip():  # skip empty lines
                x, y = map(float, line.strip().split())
                data.append([x, y])
    positions = np.array(data)


    N = len(positions)
    box_size = 60

    dist_nd_sq = np.zeros(N * (N - 1) // 2)  # to match the result of pdist
    for d in range(dim):
        pos_1d = positions[:, d][:, np.newaxis]  # shape (N, 1)
        dist_1d = pdist(pos_1d)  # shape (N * (N - 1) // 2, )
        dist_1d[dist_1d > box_size * 0.5] -= box_size
        dist_nd_sq += dist_1d ** 2  # d^2 = dx^2 + dy^2 + dz^2
    dist_nd = np.sqrt(dist_nd_sq)

    min_distance = np.min(dist_nd)
    print(f"Minimum distance for file {file_number}: {min_distance}")
    return min_distance


if __name__ == "__main__":
    # Give the name of the simulation file as an argument
    try:
        sys.argv[1]
    except IndexError:
        print("Error: Missing #cores & step")
        sys.exit(1)


    with Pool(processes=int(sys.argv[1])) as pool:
        results = pool.map(min_distance, range(5, 805, int(sys.argv[2])))
    
    minimum_distance = np.min(results)
    print(f"\nMinimum distances for files: {minimum_distance}")
    # Plot the results
    steps = range(5, 805, int(sys.argv[2]))
    plt.plot(steps, results, marker='o')
    plt.xlabel('File Number')
    plt.ylabel('Minimum Distance')
    plt.title('Minimum Distance vs File Number')
    plt.grid()
    plt.show()
