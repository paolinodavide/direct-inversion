import sys
import numpy as np
from gr_pair import *
import json

with open("config.json") as f:
    params = json.load(f)
l_box = params['l_box']
wt_file = params['wt_file']

cor_wt = []
with open(wt_file, 'r') as f:
    List_wt_data = f.read()
List_wt_lines = List_wt_data.split('\n')
for line in List_wt_lines:
    if (line.split() != []) and (line.split()[0] != '#'):
        cor_wt.append(float(line.split()[0]))
n_correl_wt = len(cor_wt)

with open("positions.json") as f:
    dict_pos=json.load(f)

dict_tmp = {}


print('Computing minimal distances ...')

for i in range(0, n_correl_wt):
    index_time = int(cor_wt[i])
    print('d_min ?', index_time)
    d_min = l_box
    List_pos = dict_pos[str(index_time)]
    for x1 in List_pos:
        for x2 in List_pos:
            if (x1 != x2):
                dx = float(x2[0]) - float(x1[0])
                dy = float(x2[1]) - float(x1[1])
                if dist(dx, dy) < d_min:
                    d_min = dist(dx, dy)
                    print(d_min)
    dict_tmp[index_time] = d_min

dict_min = dict(sorted(dict_tmp.items(), key=lambda item: item[1]))

json.dump(dict_min, open('dmin.json', 'w+'))
