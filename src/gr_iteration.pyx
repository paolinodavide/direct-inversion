import numpy as np
import cython

""" Module containing functions related to the potential part of the iteration
program."""

cdef double lj_x(double x, double T):
	""" Force derived from the Lennard-Jones potential/x."""

	return 1./T*(12.*x**(-14) - 6.*x**(-8) )

cdef double lj_force(double x, double T):
	""" Lennard-Jones potential."""

	return 1./T*( x**(-12) - x**(-6) )

cdef double lj_rep(double x, double T):
	""" Repulsive part of the Lennard-Jones potential."""

	return 1./T*( x**(-12) )

cdef double lj_att(double x, double T):
	""" Attractive part of the Lennard-Jones potential."""

	return 1./T*( - x**(-6) )

cdef double lj_test(double x, double T):
	""" Test potential."""

	return 1./T*( x**(-12.00001) - x**(-6) )

cpdef double[:] get_pot(str name, int pot_length, double r_bin, double x_min, double T):
	""" Choose the potential to initiate the iteration loop. 
	Returns the potential as an array. """

	cdef double[:] pot_list = np.zeros(pot_length)
	cdef int i
	cdef double x, offset

	if (name == 'zero'):
		return pot_list
	elif (name == 'lj_full'):
		for i in range(0, pot_length):
			x = x_min + i*r_bin
			pot_list[i] = lj_force(x, T)
		offset = pot_list[pot_length - 1]
		for i in range(0, pot_length):
			pot_list[i] -= offset
		return pot_list
	elif (name == 'lj_rep'):
		for i in range(0, pot_length):
			x = x_min + i*r_bin
			pot_list[i] = lj_rep(x, T)
		offset = pot_list[pot_length - 1]
		for i in range(0, pot_length):
			pot_list[i] -= offset
		return pot_list
	elif (name == 'lj_att'):
		for i in range(0, pot_length):
			x = x_min + i*r_bin
			pot_list[i] = lj_att(x, T)
		offset = pot_list[pot_length - 1]
		for i in range(0, pot_length):
			pot_list[i] -= offset
		return pot_list
	elif (name == 'lj_test'):
		for i in range(0, pot_length):
			x = x_min + i*r_bin
			pot_list[i] = lj_test(x, T)
		offset = pot_list[pot_length - 1]
		for i in range(0, pot_length):
			pot_list[i] -= offset
		return pot_list

cpdef double[:] init_x(int pot_length, double r_bin, double x_min, double T):
	""" Choose the force/r to initiate the iteration loop. 
	Returns the force/r as an array. """

	cdef double[:] x_list = np.zeros(pot_length)
	cdef int i
	cdef double x

	for i in range(0, pot_length):
		x = x_min + i*r_bin
		x_list[i] = lj_x(x, T)
	return x_list


cpdef double[:] get_x(double[:] potential, double x_min, double r_bin, int pot_length):
	""" Computes force/r from the potential. """

	cdef int i
	cdef double[:] list_x = np.zeros(pot_length)

	for i in range(1, pot_length - 1):
		list_x[i] = - (potential[i + 1] - potential[i - 1])/( \
			2.*r_bin*(x_min + i*r_bin))
	
	list_x[0] = 2.*list_x[1] - list_x[2]
	list_x[pot_length - 1] = 2.*list_x[pot_length - 2] - list_x[pot_length - 3]

	return list_x

cpdef double[:] update_pot(int pot_length, double delta_reg, double T, \
		double[:] g_target, double[:] g_tmp, double[:] u_tmp, double r_bin, \
		double x_min, double alpha_pot):
	""" Update the potential from the g(r) difference. delta_reg can be used
	to regulate the iteration loop by forcing smaller steps."""

	cdef int i, binmin
	cdef double[:] pot_new = np.zeros(pot_length)
	# cutoff effect on noise in g(r)
	cdef double delta_noise = 0.0000001, offset

	for i in range(0, pot_length):
		pot_new[i] = u_tmp[i]

	binmin = int(x_min/r_bin)

	for i in range(0, pot_length):
		# convert the position to get the right index
		pot_new[i] += delta_reg*np.log( abs(g_tmp[i + binmin] + \
			delta_noise)/abs(g_target[i + binmin] + delta_noise) )
		#pot_new[i] /= (1. - delta_pot)
		pot_new[i] *= alpha_pot

	offset = pot_new[pot_length - 1]
	for i in range(0, pot_length):
		pot_new[i] -= offset

	return pot_new
