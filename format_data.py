# A script to transform the LAMMPS output file in a file sequence compatible
# with the postprocess codes: double column list of particles' positions
# with one file for each dump time
# Call this file from a .sh script to generate the appropriate directory

import os
import sys
import argparse
import numpy as np

# Set up argument parsing to take -i (input) and -d (directory)
parser = argparse.ArgumentParser(description="Parse LAMMPS dump and generate configs.")
parser.add_argument("-i", "--input", required=True, help="Path to the LAMMPS simulation file")
parser.add_argument("-d", "--directory", required=True, help="Output directory (<YOUR_DIR>)")
args = parser.parse_args()

simu_file = args.input
out_dir = os.path.join(args.directory, "inputs")  # Base output directory for all generated files

# Define the dynamic paths based on your <YOUR_DIR>
# Note: If these need to sit inside an "inputs" folder to match your README, 
# you can change this to os.path.join(out_dir, "inputs", "configs")
configs_dir = os.path.join(out_dir, "configs")
configs_npy_dir = os.path.join(out_dir, "configs_npy") # Or configs_bin
times_file_path = os.path.join(out_dir, "gr_wt.dat")

# Ensure the base directory exists
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Common part of the file name (using basename to ignore directory paths if provided)
prefix = os.path.basename(simu_file).split('_')[0] + '_'
suffix = '.dat'

List_column = []
with open(simu_file, "r") as f:
    List_data = f.read()
List_lines = List_data.split('\n')
for line in List_lines:
    if (line.split() != []):
        List_column.append(line.split())

index = 0
simu_time = ''
copy = False
read_lbox = True

# This block is only if the times have to be saved in a file
save_times = open(times_file_path, "w")
save_times.write("# Waiting times \n")

while (index < len(List_column)):
    if copy:
        while (List_column[index][0] != 'ITEM:'):
            x = float(List_column[index][2])*l_box
            y = float(List_column[index][3])*l_box
            save_results.write(str(x) + '\t' + str(y) + '\n')
            index += 1
            try:
                (List_column[index][0] != 'ITEM:')
            except IndexError:
                # The program has reached the end of the file
                break
        save_results.close()
        copy = False
        try:
            (List_column[index][0] != 'ITEM:')
        except IndexError:
            # The program has reached the end of the file
            break
    if (List_column[index][0] == 'ITEM:'):
        if (List_column[index][1] == 'TIMESTEP'):
            simu_time = List_column[index + 1][0]
            filename = prefix + simu_time + suffix

            # This block is only if the times have to be saved in a file
            save_times.write(simu_time + '\n')

            index += 1
        elif (List_column[index][1] == 'BOX') and read_lbox:
            index +=1
            l_box = float(List_column[index][1])
            read_lbox = False
            index += 3
        elif (List_column[index][1] == 'ATOMS'):
            copy = True
            if not os.path.exists(configs_dir):
                os.makedirs(configs_dir)
            save_results = open(os.path.join(configs_dir, filename), "w")
            save_results.write('# x \t y \t l_box = ' + str(l_box) + '\n')
            index += 1
        else:
            index += 1
    else:
        index += 1

save_times.close()
print(f"Configs files created in {configs_dir}")

# Convert the .dat files to .npy files
# Create the output directory if it doesn't exist
if not os.path.exists(configs_npy_dir):
    os.makedirs(configs_npy_dir)

# Iterate through all .dat files in the configs directory and process them
for file_name in filter(lambda f: f.endswith(".dat"), os.listdir(configs_dir)):
    file_path = os.path.join(configs_dir, file_name)
    
    # Load data directly into a numpy array, skipping header lines
    data_array = np.loadtxt(file_path, comments="#", usecols=(0, 1))
    
    # Save the numpy array to a .npy file
    npy_file_path = os.path.join(configs_npy_dir, file_name.replace(".dat", ".npy"))
    np.save(npy_file_path, data_array)

print(f"Conversion to .npy files completed. Saved in {configs_npy_dir}")