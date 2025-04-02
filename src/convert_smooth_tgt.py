import numpy as np
import sys
import json

with open("smooth.json") as f:
    dict_smooth = json.load(f)

save_tgt = open("gs_target.dat", "w+")
save_tgt.write("# r \t g(r) \n")
for key in dict_smooth.keys():
    save_tgt.write(str(key) + "\t" + str(dict_smooth[key]) + "\n")
save_tgt.close()
