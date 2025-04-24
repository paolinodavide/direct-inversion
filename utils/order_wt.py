import numpy as np 
import json 

with open("config.json") as f:
    params = json.load(f)

wt_file = params['wt_file']

with open(wt_file, 'r') as f:
    List_wt_data = f.read()

List_wt_lines = List_wt_data.split('\n')

header = List_wt_lines[0]
data_lines = List_wt_lines[1:]
np.random.shuffle(data_lines)
List_wt_lines = [header] + data_lines

output_file = '../randomized_wt.dat'

with open(output_file, 'w') as f:
    f.write('\n'.join(List_wt_lines))