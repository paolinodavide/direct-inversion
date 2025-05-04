import sys
import json

# Convert the output json files into properly named .dat files
# Needs the name of the json file to parse

try:
    filename = str(sys.argv[1])
except IndexError:
    print("Enter the name of the json file to parse")

str_test = ''
for i in range(-5, 0):
    str_test += filename[i]
if (str_test != '.json'):
    raise ValueError('Wrong extension in the file name')

dict_res = {}
with open(filename, "r") as f:
    dict_res = json.load(f)
list_keys = list(key for key in dict_res.keys())

suffix = filename.replace('.json', '')

if (suffix == 'err'):
    data_file = suffix + ".dat"
    save_data = open(data_file, 'w')
    save_data.write('# Delta(g)^2 \t Delta(u)^2 \t delta_tgt \t delta_pot' + \
            '\t alpha_pot \t Delta(iter)^2 \n')
    for i in range(len(dict_res['gr'])):
        save_data.write(str(dict_res['gr'][i]) + "\t" + str(dict_res['pot'][i]) + \
                "\t" + str(dict_res['delta_tgt'][i]) + "\t" + \
                str(dict_res['delta_pot'][i]) + "\t" + str(dict_res['alpha_pot'][i]) + \
                "\t" + str(dict_res['iteration'][i]) + "\n")
    save_data.close()

for key in list_keys:
    if (suffix == 'positions'):
        data_file = suffix + "_" + str(key) + ".dat"
        save_data = open(data_file, 'w')
        if (key == 'wt'):
            save_data.write("# waiting time list \n")
            for i in range(len(dict_res[key])):
                save_data.write(str(dict_res[key][i]) + "\n")
            save_data.close()
        else:
            save_data.write("# x \t y \n")
            for i in range(len(dict_res[key])):
                save_data.write(str(dict_res[key][i][0]) + "\t" + \
                        str(dict_res[key][i][1]) + "\n")
            save_data.close()
    elif (suffix == 'pot'):
        data_file = suffix + "_" + str(key) + ".dat"
        save_data = open(data_file, 'w')
        save_data.write("# r \t u(r) \t r*force(r) \n")
        for j in range(len(dict_res[key])):
            save_data.write(str(dict_res[key][j][0]) + "\t" + str(dict_res[key][j][1]) + \
                    "\t" + str(dict_res[key][j][2]) + "\n")
        save_data.close()
    elif (suffix == 'gr'):
        i_best = dict_res['best']
        if (key != 'best'):
            data_file = suffix + "_" + str(key) + ".dat"
            save_data = open(data_file, 'w')
            save_data.write("# r \t g(r) \t best index: " + str(i_best) + "\n")
            for j in range(len(dict_res[key])):
                save_data.write(str(dict_res[key][j][0]) + "\t" + \
                        str(dict_res[key][j][1]) + "\n")
        save_data.close()
    elif (suffix == 'err'):
        pass
    else:
        print("Error: this file is not appropriate")

