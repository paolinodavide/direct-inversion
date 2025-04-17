# A script to transform the LAMMPS output file in a file sequence compatible
# with the postprocess codes: double column list of particles' positions
# with one file for each dump time
# Call this file from a .sh script to generate the appropriate directory

import sys
import os

# Give the name of the simulation file as an argument
try:
    sys.argv[1]
except IndexError:
    print("Error: Missing file name")
    sys.exit(1)

simu_file = str(sys.argv[1])

# Common part of the file name
prefix = 'lj_'
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

#This block is only if the times have to be saved in a file
save_times = open("gr_wt.dat", "w")
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

            #This block is only if the times have to be saved in a file
            save_times.write(simu_time + '\n')

            index += 1
        elif (List_column[index][1] == 'BOX') and read_lbox:
            index +=1
            l_box = float(List_column[index][1])
            read_lbox = False
            index += 3
        elif (List_column[index][1] == 'ATOMS'):
            copy = True
            if not os.path.exists("./configs"):
                os.makedirs("./configs")
            save_results = open(f"./configs/{filename}", "w")
            save_results.write('# x \t y \t l_box = ' + str(l_box) + '\n')
            index += 1
        else:
            index += 1
    else:
        index += 1

save_times.close()
