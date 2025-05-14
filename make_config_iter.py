import json

# A script to generate the configuration file
# for the computation of  g(r)

params = {}
# number of particles of type A
params['n_a_part'] = 2025
# number of particles of type B
params['n_b_part'] = 0
# fixes the length of the waiting time axis
params['n_totstep'] = 500000
# initial length of list of waiting times
params['n_correl_wt'] = 1000
# maximum number of waiting times
params['n_max_wt'] = 125
# log time axis or linear axis (False -> linear)
params['log_lin'] = False
# size of the simulation box
params['l_box'] = 60
# size of the histogram bin for g(r) 
params['r_bin'] = 0.002
# prefix for the location of datafiles
params['prefix_file'] = '../configs_npy/wca' #Recall to change this
# output file for waitingtimes
params['wt_file'] = '../ordered_wt.dat'
# cutoff in the lj potential
params['x_cut'] = 2.5
# size of the r axis in the g(r) plot
params['qdim_max'] = 5_000
params['qdim'] = 60000
# minimum length for the potential
params['x_min'] = 0.881
params['x_low'] = 0.93
# file for target g(r)
params['target_file'] = '../gs_target.dat'
# target precision
params['target_precision'] = 0
params['output_file'] = 'gr_final.dat'
params['max_iter'] = 50
params['method'] = 'out'
params['Temperature'] = 1
params['init_pot'] = 'lj_full'
params['target_pot'] = 'lj_full'
params['delta_reg'] = 0.2
# Decide if begins from an old config, or from scratch
# Needs to store it as a string, json does not like booleans
params['INIT'] = True
params['GR'] = 'HISTO'

json.dump(params, open('config.json', 'w'))
