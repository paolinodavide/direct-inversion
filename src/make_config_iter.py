import json

# A script to generate the configuration file
# for the computation of  g(r)

params = {}
# number of particles of type A
params['n_a_part'] = 2402
# number of particles of type B
params['n_b_part'] = 0
# fixes the length of the waiting time axis
params['n_totstep'] = 800000
# initial length of list of waiting times
params['n_correl_wt'] = 1000
# log time axis or linear axis (False -> linear)
params['log_lin'] = False
# size of the simulation box
params['l_box'] = 60
# size of the histogram bin for g(r) (?)
params['r_bin'] = 0.002
# prefix for the location of datafiles
params['prefix_file'] = '/home/davide/OneDrive/POLITO/Magistrale/Internship/W2_MD/ljs/lj_'
# output file for waitingtimes
params['wt_file'] = 'list_wt.dat'
# cutoff in the lj potential
params['x_cut'] = 2.5
# size of the r axis in the g(r) plot
params['qdim_max'] = 5000
params['qdim'] = 60000
# minimum length for the potential
params['x_min'] = 0.878
params['x_low'] = 0.934
# file for target g(r)
params['target_file'] = 'gs_target.dat'
# target precision
params['target_precision'] = 0
params['output_file'] = 'gr_final.dat'
params['max_iter'] = 15
params['method'] = 'in'
params['Temperature'] = 1
params['init_pot'] = 'pot_15.dat'
params['delta_reg'] = 0.05
# Decide if begins from an old config, or from scratch
# Needs to store it as a string, json does not like booleans
params['INIT'] = False
params['GR'] = 'HISTO'

# Generate configuration file
json.dump(params, open('config.json', 'w'))
