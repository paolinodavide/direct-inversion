import json
import sys
#import ipdb
import numpy as np
from gr_function import *
from time import time

# A script to compute g(r) from a set of positions
# For now, the script only works in 2d
dimension = 2

# imports parameters from configuration file
# see definition in make_config.py
with open("config.json") as f:
    params = json.load(f)

# Initialisation of parameters
# number of particles
n_part = params['n_a_part'] + params['n_b_part']
pi = np.pi
# simulation box' size
l_box = params['l_box']
# particle's density
rho = n_part/l_box**2
log_lin = params['log_lin']
n_totstep = params['n_totstep']
n_correl_wt = params['n_correl_wt']
prefix_file = params['prefix_file']
r_bin = params['r_bin']
qdim = params['qdim']
qdim_max = params['qdim_max']
x_cut = params['x_cut']

# Generate the waiting time list and store it
#n_correl_wt_new, cor_wt = gen_wt(log_lin, n_correl_wt, n_totstep)
#write_wt = open("gr_wt.dat", "w") # "a" for append, does not overwrite, "w" for write
#write_wt.write("# Waiting times, n_correl_wt_new = " + str(n_correl_wt_new) + "\n")
#for j in cor_wt:
#    write_wt.write(str(j) + '\n')
#write_wt.close()

# In this version of the code, we read the waiting times from the simulations
cor_wt = []
with open('list_wt.dat', 'r') as f:
    List_wt_data = f.read()
List_wt_lines = List_wt_data.split('\n')
for line in List_wt_lines:
    if (line.split() != []) and (line.split()[0] != '#'):
        cor_wt.append(float(line.split()[0]))
n_correl_wt_new = len(cor_wt)
# the initial configuration is not thermalised, we remove it
#n_correl_wt_new -= 1

# For test purpose, we just test on one file
#cor_wt = [6000]
#n_correl_wt_new = 1

# Pair correlation function from Borgis' formula
# both with the inner and the outer definitions
g_hist = np.zeros(qdim)
g_dum = np.zeros(qdim)
g_inter = np.zeros(qdim)

# Computation of g(r)
# set the prefactor: n_correl_wt_ because average over time, 2 and 4pi from the
# formula, and 4 because N_a N_b
prefactor = n_correl_wt_new*pi*4.*n_part*rho/4.
prefactor2 = pi*4.*n_part*rho/4.
for i in range(0, n_correl_wt_new):
    # define the time in the name of the data file
    index_time = int(cor_wt[i])
    print(prefix_file, index_time)
    List_pos = gen_pos(prefix_file, index_time)

    # generate the correlation function
    print('analysis')
    g_old = g_hist
    g_hist = gen_grHisto(List_pos, n_part, qdim_max, l_box, r_bin, g_old)

    g_inter = gen_grHisto(List_pos, n_part, qdim_max, l_box, r_bin, g_dum)
    write_inter = open("gr_h_slj_" + str(i) + ".dat", "w+")
    write_inter.write('# r \t g_hist(r) \n')
    for j in range(0, qdim + 1):
        if (float(j)*r_bin <= 10):
            write_inter.write(str(float(j)*r_bin) + '\t' + str(g_inter[j]/prefactor2) + '\n')
    write_inter.close()

# write the results in output file
write_results = open("gr_h_" + "slj" + ".dat", "w")
write_results.write('# r \t g_hist(r) \n')
for j in range(0, qdim + 1):
    if (float(j)*r_bin <= 10):
        write_results.write(str(float(j)*r_bin) + '\t' + str(g_hist[j]/prefactor) + '\n')
write_results.close()

