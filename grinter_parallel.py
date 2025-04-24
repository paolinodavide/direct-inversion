import json
import sys
import numpy as np
from gr_pair import *
from gr_iteration import *
import xgboost
from gr_pair_parallel import *

# Load parameters
with open("config.json") as f:
    params = json.load(f)

n_part = params['n_a_part'] + params['n_b_part']
pi = np.pi
l_box = params['l_box']
rho = n_part/l_box**2
r_bin = params['r_bin']
qdim = params['qdim_max']
x_cut = params['x_cut']
x_min = params['x_min']
x_low = params['x_low']
x_max = 10
binmin = int(x_min/r_bin)
binlow = int(x_low/r_bin)
bincut = int(x_cut/r_bin)
binmax = int(10/r_bin)
target_file = params['target_file']
fileout = params['output_file']
max_iter = params['max_iter']
method = params['method']
T = params['Temperature']
init_pot = params['init_pot']
delta_reg = params['delta_reg']
INIT = params['INIT']
GR = params['GR']

# store the data
dict_pot = {}
dict_gr = {}
dict_err = {}

# Initialise waiting times
with open("positions.json") as f:
    dict_pos = json.load(f)
cor_wt = dict_pos['wt']
n_correl_wt = len(cor_wt)

# Initialise the potential
pot_length = bincut - binlow
u_tmp = get_pot('lj_full', pot_length, r_bin, x_low, T)
x_tmp = init_x(pot_length, r_bin, x_low, T)
u_tgt = np.zeros(pot_length)
for j in range(0, pot_length):
    u_tgt[j] = u_tmp[j]

dict_pot['target'] = []
for j in range(0, len(u_tmp)):
    dict_pot['target'].append([x_low + j*r_bin, u_tmp[j], x_tmp[j]]) 
# Initialize g(r)
prefactor = 2.*pi*n_part*rho*n_correl_wt/4.
g_tgt = np.zeros(qdim)
if (GR == 'FORCE'):
    gr_tmp = gen_grBorgis_cpp(n_part, qdim, l_box, r_bin, x_max,  \
            x_cut, binmin, binlow, x_tmp)
    for wt in cor_wt:
        for k in range(0, qdim):
            gr_tmp['tot'][k] += gr_tmp[wt][k]

    for i in range(0, binmax):
        if (method == 'in'):
            g_tgt[i] = 1. - g_tmp['tot'][i]/prefactor

elif (GR == 'HISTO'):
    size = 0
    with open(target_file, 'r') as f:
        list_target_data = f.read()
    list_target_lines = list_target_data.split("\n")
    for line in list_target_lines:
        if (line.split() != []) and (line.split()[0] != '#'):
            size += 1
    list_target = np.zeros((size, 2))
    i = 0
    for line in list_target_lines:
        if (line.split() != []) and (line.split()[0] != '#'):
            list_target[i] = [float(line.split("\t")[0]), float(line.split("\t")[1])]
            i += 1
    for j in range(0, qdim):
        g_tgt[j] = list_target[j][1]

# Prepare the error computation
g_error_tgt = np.zeros(pot_length)
for i in range(binlow, bincut):
    g_error_tgt[i - binlow] = g_tgt[i]
delta_tgt = g_tgt[binlow]

dict_gr['target'] = []
for j in range(0, qdim):
    dict_gr['target'].append([j*r_bin, g_tgt[j]])

# Initialise the potential
#if (INIT == 'True'):
if (INIT):
    u_tmp = get_pot(init_pot, pot_length, r_bin, x_low, T)
    x_tmp = get_x(u_tmp, x_low, r_bin, pot_length)
else:
    # Initialize u on previous data
    data_u = np.genfromtxt(f'{str(init_pot)}.dat', delimiter = '\t', usecols = (0,1,2))
    for i in range(0, pot_length):
        u_tmp[i] = data_u[i][1]
        x_tmp[i] = data_u[i][2]

dict_pot['0'] = []
for j in range(0, len(u_tmp)):
    dict_pot['0'].append([x_low + j*r_bin, u_tmp[j], x_tmp[j]])

precision = params['target_precision']
iter_index = 1
final_index = 1
Dg2 = 1.

dict_err['gr'] = []
dict_err['pot'] = []
dict_err['delta_tgt'] = []
dict_err['delta_pot'] = []
dict_err['alpha_pot'] = []

#while (precision > eps):
while True:

    if (iter_index == max_iter + 1):
        print('Error: Target precision is not reached')
        gr_final = gr_res
        break

    print('Iteration #: ', iter_index, '...')
    gr_tmp = gen_grBorgis_cpp(n_part, qdim, l_box, r_bin, x_max, \
            x_cut, binmin, binlow, x_tmp)
    #for wt in cor_wt:
    #    for k in range(0, qdim):
    #        gr_tmp['tot'][k] += gr_tmp[wt][k]
    delta_pot = 1 - gr_tmp[binlow]/prefactor

    gr_res =np.zeros(qdim)
    for j in range(0, binmax):
        if (method == 'in'):
            gr_res[j] = 1. - gr_tmp[j]/prefactor
            gr_res[j] = abs( ((1. - delta_tgt)*gr_res[j] + (delta_tgt - delta_pot) \
                    )/(1. - delta_pot) )
        elif (method == 'out'):
            gr_res[j] = gr_tmp[j]/prefactor
            gr_res[j] = gr_res[j]/(delta_pot)

    g_error_tmp = np.zeros(pot_length)
    for j in range(binlow, bincut):
        g_error_tmp[j - binlow] = gr_res[j]
    error1 = get_error(g_error_tgt, g_error_tmp, x_low, x_cut, r_bin)
    error2 = get_error(u_tgt, u_tmp, x_low, x_cut, r_bin)

    dict_err['gr'].append(error1)
    dict_err['pot'].append(error2)
    dict_err['delta_tgt'].append(delta_tgt)
    dict_err['delta_pot'].append(delta_pot)

    print('precision: ', error1, 'Error pot:', error2)
    print('delta_tgt:', delta_tgt, 'delta_pot:', delta_pot)
    print('Delta alpha_pot : ', (1. - delta_tgt)/(1 - delta_pot) - 1.)
    if (error1 < Dg2):
        final_index = iter_index
        gr_final = gr_res
        print('best index = ', final_index)
        Dg2 = error1

    dict_gr[str(iter_index)] = []
    for j in range(0, qdim):
        dict_gr[str(iter_index)].append([j*r_bin, gr_res[j]])

    alpha_pot = (1. - delta_tgt)/(1 - delta_pot)
    u_tmp = update_pot(pot_length, delta_reg, T, g_tgt, gr_res, u_tmp, r_bin, \
            x_low, alpha_pot)
    x_tmp = get_x(u_tmp, x_low, r_bin, pot_length)

    dict_pot[str(iter_index)] = []
    for j in range(0, len(u_tmp)):
        dict_pot[str(iter_index)].append([x_low + j*r_bin, u_tmp[j], x_tmp[j]])

    dict_err['alpha_pot'].append(alpha_pot)

    if (error1 < precision):
        print('Target precision reached at iteration ', iter_index)
        print('Relative error: ', precision)
        gr_final = gr_res
        break
    
    iter_index += 1

#pool.close()

print('Best index: ', final_index)
dict_gr['best'] = final_index

print('... Saving results ...')
json.dump(dict_gr, open('gr.json', 'w+'))
json.dump(dict_pot, open('pot.json', 'w+'))
json.dump(dict_err, open('err.json', 'w+'))
