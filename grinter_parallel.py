import json
import numpy as np
import gr_pair as pair
import gr_iteration as it
import borgis as grNb
import time 

def load_config():
    with open("config.json") as f:
        return json.load(f)

def load_positions():
    with open("positions.json") as f:
        return json.load(f)['wt']

def initialize_potential(params, pot_length, r_bin, x_low):
    if params['INIT']:
        u = it.get_pot(params['init_pot'], pot_length, r_bin, x_low, params['Temperature'])
        x = it.get_x(u, x_low, r_bin, pot_length)
    else:
        data = np.genfromtxt(f"{params['init_pot']}.dat", delimiter='\t', usecols=(0, 1, 2))
        u = data[:, 1]
        x = data[:, 2]
    return u, x

def load_target_gr(params, qdim):
    g_target = np.zeros(qdim)
    if params['GR'] == 'HISTO':
        with open(params['target_file'], 'r') as f:
            lines = [line for line in f if line.strip() and not line.startswith('#')]
        data = np.array([[float(v) for v in line.strip().split()] for line in lines])
        g_target[:len(data)] = data[:, 1]
    return g_target

def compute_prefactor(n_part, rho, n_correl_wt, temperature=1.0):
    """ No need anymore to divide by n_correl_wt, grNb returns the average 
    The 4 is needed as the LJ potential in it.get_pot is divided by 4.
    Temperature is incuded in the potential."""
    return 2 * np.pi * rho * n_part / 4

def build_dict(entries):
    return [[x, y, z] for x, y, z in zip(*entries)]

def compute_error_metrics(g_target, g_current, u_target, u_current, x_low, x_cut, r_bin):
    e_gr = pair.get_error(g_target, g_current, x_low, x_cut, r_bin)
    e_pot = pair.get_error(u_target, u_current, x_low, x_cut, r_bin)
    return e_gr, e_pot

def main():
    params = load_config()
    cor_wt = load_positions()

    n_part = params['n_a_part'] + params['n_b_part']
    l_box = params['l_box']
    rho = n_part / l_box**2
    r_bin = params['r_bin']
    qdim = params['qdim_max']
    x_min, x_low, x_cut = params['x_min'], params['x_low'], params['x_cut']
    binmin, binlow, bincut, binmax = map(lambda x: int(x / r_bin), [x_min, x_low, x_cut, 10])
    pot_length = bincut - binlow

    g_target_full = load_target_gr(params, qdim)
    g_target = g_target_full[binlow:bincut]
    delta_tgt = g_target_full[binlow]

    u_target = it.get_pot('lj_full', pot_length, r_bin, x_low, params['Temperature'])
    x_target = it.init_x(pot_length, r_bin, x_low, params['Temperature'])

    dict_gr = {'target': [[j * r_bin, g_target_full[j]] for j in range(qdim)]}
    dict_pot = {'target': [[x_low + j*r_bin, u_target[j], x_target[j]] for j in range(pot_length)]}
    dict_err = {'gr': [], 'pot': [], 'delta_tgt': [], 'delta_pot': [], 'alpha_pot': []}

    u_current, x_current = initialize_potential(params, pot_length, r_bin, x_low)
    dict_pot['0'] = [[x_low + j*r_bin, u_current[j], x_current[j]] for j in range(pot_length)]

    prefactor = compute_prefactor(n_part, rho, len(cor_wt))
    precision = params['target_precision']
    max_iter = params['max_iter']
    method = params['method']
    delta_reg = params['delta_reg']
    T = params['Temperature']

    best_error = 1.
    best_index = 1
    iteration = 1

    startTime = time.time()
    while True:
        if iteration > max_iter:
            print("Error: Target precision not reached.")

            break
        if iteration > 1:
            g_previous = g_current

        print(f"Iteration #{iteration} ...")


        x_current = np.asarray(x_current, dtype=np.float64)

        gr_raw = grNb.gBorgis_parallel(l_box, x_current, x_low, 10.0, r_bin, x_cut , method='out', min_radius=0, 
                configs_path=params['prefix_file'], ordered_indices_file=params['wt_file'], max_config_nb=params['n_max_wt'])


        gr_res = np.zeros(qdim)
        for j in range(binmax):
            if method == 'in':
                delta_pot = gr_raw[binlow] / prefactor
                temp =  gr_raw[j] / prefactor 
                gr_res[j] = abs(((1. - delta_tgt) * temp + (delta_tgt - delta_pot)) / (1. - delta_pot))
            elif method == 'out':
                temp =  1 - gr_raw[j] / prefactor #edit eh 1-
                delta_pot = 1 - gr_raw[binlow] / prefactor
                gr_res[j] = abs(((1. - delta_tgt) * temp + (delta_tgt - delta_pot)) / (1. - delta_pot))

        g_current = gr_res[binlow:bincut]

        err_gr = pair.get_error(g_target, g_current, x_low, x_cut, r_bin)
        err_pot = pair.get_error(u_target, u_current, x_low, x_cut, r_bin)
        err_iteration = pair.get_error(g_current, g_previous, x_low, x_cut, r_bin) if iteration > 1 else 0
        
        dict_err['iteration'] = dict_err.get('iteration', [])
        dict_err['iteration'].append(err_iteration)
        dict_err['gr'].append(err_gr)
        dict_err['pot'].append(err_pot)
        dict_err['delta_tgt'].append(delta_tgt)
        dict_err['delta_pot'].append(delta_pot)

        print(f"Precision: {err_gr:.4e}, Error pot: {err_pot:.4e}")
        print(f"Delta tgt: {delta_tgt:.4e}, Delta pot: {delta_pot:.4e}")
        print(f"Delta alpha_pot: {(1. - delta_tgt)/(1 - delta_pot) - 1.:.4e}")
        print(f"Error iteration: {err_iteration:.4e}")

        if err_gr < best_error:
            best_error = err_gr
            best_index = iteration
            gr_final = gr_res
            print("Best index =", best_index)

        dict_gr[str(iteration)] = [[j * r_bin, gr_res[j]] for j in range(qdim)]

        alpha = (1. - delta_tgt) / (1 - delta_pot)
        dict_err['alpha_pot'].append(alpha)

        u_current = it.update_pot(pot_length, delta_reg, T, g_target_full, gr_res, u_current, r_bin, x_low, alpha)
        x_current = it.get_x(u_current, x_low, r_bin, pot_length)

        dict_pot[str(iteration)] = [[x_low + j*r_bin, u_current[j], x_current[j]] for j in range(pot_length)]

        if err_gr < precision:
            print(f"Target precision reached at iteration {iteration}")
            break

        iteration += 1

    endTime = time.time()
    print(f"Total time: {(endTime - startTime) / 60:.2f} minutes")
    print("Best index:", best_index)
    dict_gr['best'] = best_index

    print("... Saving results ...")
    json.dump(dict_gr, open('gr.json', 'w+'))
    json.dump(dict_pot, open('pot.json', 'w+'))
    json.dump(dict_err, open('err.json', 'w+'))

if __name__ == "__main__":
    main()
