import json
import sys
import numpy as np

with open("config.json") as f:
    params = json.load(f)

wt_file = params['wt_file']
#wt_file = '../../Morse_T1/CONF_Mo/gr_wt.dat'
prefix_file = params['prefix_file']
n_part = params['n_a_part'] + params['n_b_part']
max_nb_config = params['n_max_wt']

# Initialise waiting times
cor_wt = []
with open(wt_file, 'r') as f:
    wt_data = f.read()
wt_lines = wt_data.split("\n")

# EDITTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
for i, line in enumerate(wt_lines): #
    if i > max_nb_config + 1: 
        break #
    if ((line.split() != []) and (line.split()[0][0] != '#')):
        cor_wt.append(int(line.split()[0]))
dict_pos = {}
n_correl = len(cor_wt)

#array_wt = np.zeros(n_correl)
List_pos = np.zeros((n_part, 2))
for i in range(0, n_correl):
#    array_wt[i] = cor_wt[i]
    wt = cor_wt[i]
    filename = str(prefix_file) + str(int(wt)) + '.dat'
    with open(filename, "r") as f:
        List_data = f.read()
    List_lines = List_data.split("\n")
    index = 0
    for j in range(0, len(List_lines)):
        if ((List_lines[j].split() != []) and (List_lines[j].split()[0] != '#')):
            List_pos[index][0] = float(List_lines[j].split()[0])
            List_pos[index][1] = float(List_lines[j].split()[1])
            index += 1
    dict_pos[wt] = List_pos.tolist()
#dict_pos['wt'] = array_wt
dict_pos['wt'] = cor_wt

json.dump(dict_pos, open('positions.json', 'w'))
