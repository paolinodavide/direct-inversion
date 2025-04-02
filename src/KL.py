import sys
import json
import numpy as np

# Computes the Kullback-Leibler divergence and information

try:
    filename = str(sys.argv[1])
except IndexError:
    print("Enter the name of the json file to parse.")

str_test = ''
for i in range(-5, 0):
    str_test += filename[i]
if (str_test != '.json'):
    raise ValueError('Wrong extension in the filename')

dict_res = {}
with open(filename, "r") as f:
    dict_res = json.load(f)
list_keys = list(key for key in dict_res.keys())

suffix = filename.replace('.json', '')

with open('config.json', 'r') as f:
    dict_conf = json.load(f)
r_bin = dict_conf['r_bin']
x_low = dict_conf['x_low']
x_cut = dict_conf['x_cut']
n_low = int(x_low/r_bin)
n_cut = int(x_cut/r_bin)

save_data = open('KL.dat', "w")
save_data.write('# I_KL \t J_KL \n')
for key in list_keys:
    if ((key == 'target') or (key == 'best')):
        continue
    else:
        data_file = 'gr_' + str(key) + '.dat'
        IKL = 0
        JKL = 0
        for j in range(n_low, n_cut):
            IKL += dict_res[key][j][1]*np.log(dict_res[key][j][1]/dict_res['target'][j][1])*r_bin
            JKL += (dict_res[key][j][1] - dict_res['target'][j][1])*np.log( \
                    dict_res[key][j][1]/ dict_res['target'][j][1])*r_bin
        save_data.write(str(IKL) + '\t' + str(JKL) + '\n')
save_data.close()
