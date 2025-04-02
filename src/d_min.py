import sys
import numpy as np
from gr_pair import *
import json

cor_wt = []
with open('test_wt.dat', 'r') as f:
    List_wt_data = f.read()
List_wt_lines = List_wt_data.split('\n')
for line in List_wt_lines:
    if (line.split() != []) and (line.split()[0] != '#'):
        cor_wt.append(float(line.split()[0]))
n_correl_wt_new = len(cor_wt)

with open("config.json") as f:
    params = json.load(f)
l_box = params['l_box']
prefix_file = params['prefix_file']

d_min = l_box

with open("positions.json") as f:
    dict_pos=json.load(f)

print('Computing minimal distance ...')

for i in range(0, n_correl_wt_new):
    index_time = int(cor_wt[i])
    print(prefix_file, index_time)
    List_pos = dict_pos[str(index_time)]
    for x1 in List_pos:
        for x2 in List_pos:
            if (x1 != x2):
                dx = float(x2[0]) - float(x1[0])
                dy = float(x2[1]) - float(x1[1])
                if dist(dx, dy) < d_min:
                    d_min = dist(dx, dy)
                    print(d_min)

print('d_min = ', d_min)
