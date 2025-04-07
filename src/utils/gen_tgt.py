import numpy as np
import json
import sys
from csaps import csaps
import matplotlib.pyplot as plt
import matplotlib.mathtext

################# Save or not save #######################################
try:
    sys.argv[1]
except IndexError:
    sys.argv.append(False)

if (sys.argv[1] == 'S'):
    save = True
else:
    save = False

################### Figure settings ######################################
# # good latex font
# from matplotlib import rc
# rc('text', usetex = True)
# #style
# plt.style.use('figure')
# # put figure cleverly
# plt.rcParams.update({'figure.autolayout': True})

################# Program ################################################
input_filename = './rdfs_0.002/g_r_h_avg.dat'

# Create a json file with the g(r) target
dict_tgt = {}
with open(input_filename, 'r') as f:
    list_data = f.read()
list_lines = list_data.split("\n")
for line in list_lines:
    if (line.split() != []) and (line.split()[0] != '#'):
        dict_tgt[line.split()[0]] = line.split()[1]
json.dump(dict_tgt, open('histo.json', 'w+'))

with open('histo.json') as f:
    target = json.load(f)
positions = list(target.keys())
g_tgt = list(target.values())

# New version: cut in two
xi1 = []
for i in range(0,143):
    xi1.append(0.80 + i*0.002)
xi2 = []
for i in range(0,1501):
    xi2.append(1.000 + i*0.002)

pos_red1 = []
g_red1 = []
pos_red2 = []
g_red2 = []
for i in range(0, len(positions)):
    if float(positions[i])>= 0.80 and float(positions[i]) <= 1.084:
        pos_red1.append(float(positions[i]))
        g_red1.append(float(g_tgt[i]))
    if float(positions[i])>= 1.000 and float(positions[i]) <= 4.0:
        pos_red2.append(float(positions[i]))
        g_red2.append(float(g_tgt[i]))

g_test1 = csaps(pos_red1, g_red1, xi1, smooth = 0.999999)
g_test2 = csaps(pos_red2, g_red2, xi2, smooth = 0.99995)

###### Redefine the g(r) function in the range 0.8-4.0

pos_red3 = []
g_red3 = []
for i in range(0, len(pos_red1)):
    if float(pos_red1[i])<= 1.024:
        pos_red3.append(float(pos_red1[i]))
        g_red3.append(float(g_test1[i]))
    elif float(pos_red1[i] == 1.024):
        pos_red3.append(float(pos_red1[i]))
        j = 0
        for k in range(0, len(pos_red2)):
            if (pos_red2[k] == 1.024):
                j = k
                print(j)
        g_red3.append(0.5*float(g_test1[i]) + 0.5*float(g_test2[j]))


for i in range(0, len(pos_red2)):
    if float(pos_red2[i]) > 1.024:
        pos_red3.append(float(pos_red2[i]))
        g_red3.append(float(g_test2[i]))

##pos = np.linspace(0., 10., 5000)
pos = []
for i in range(0, 5001):
    pos.append(i*0.002)


dict_test = {}
for i in range(0, len(pos)):
    if pos[i]>= 0.80 and pos[i] < 4.0: ## modified
        dict_test[pos[i]] = g_red3[i-400] #-400 ???????
    elif pos[i] < 0.8:
        dict_test[pos[i]] = 0
    else:
        dict_test[pos[i]] = 1

json.dump(dict_test, open('smooth.json', 'w+'))

###################### Figure ###############################################
#
fig, ax = plt.subplots()
#
ax.plot(pos_red1, g_red1, color = 'blue', marker = 'o', markersize = '5', \
        linestyle = 'none', fillstyle = 'none', label = r'$g_{histo}$')
ax.plot(pos_red2, g_red2, color = 'blue', marker = 'o', markersize = '5', \
        linestyle = 'none', fillstyle = 'none')
#ax.plot(xi1, g_test1, color = 'red', marker = 'o', markersize = '5', linestyle = 'none', label = r'$g_{smooth}$')
#ax.plot(xi2, g_test2, color = 'red', marker = 'o', markersize = '5', linestyle = 'none')
ax.plot(pos_red3, g_red3, color = 'red', label = r'$g_{smooth}$')
ax.set_xlabel(r'$r$')
ax.set_xlim(0.80, 4.0)

ax.legend()

plt.subplots_adjust()

# see or save
if (save == True):
    plt.savefig('TestCsaps2.pdf', transparent = False)
else:
    plt.show()
