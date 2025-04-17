import numpy as np
import cython
from gr_iteration import *
""" Module containing all the functions for the computation of g(r) """

cdef double make_periodic(double d_pos, double l_box):
	""" Returns a position corrected to account for the periodic boundary
	conditions in a square box of size l_box """

	if (d_pos >= l_box/2.):
		return d_pos - l_box
	elif (d_pos <= -l_box/2.):
		return d_pos + l_box
	else:
		return d_pos

cpdef double dist(double d_x,double d_y):
	""" returns the distance between two points from the coordinate differences"""

	return np.sqrt(d_x**2 + d_y**2)

cdef double[:,:] comput_frc(double[:,:] List_pos, int n_part, float x_cut, 
		double l_box, double r_bin, int binmin, int binlow, double[:] rforce):
	""" Compute the list of forces acting on the particles from their positions 
	and a force profile given in the form of a list. Careful: the program is written in 
	such a way that the force profile is computed on the appropriate grid."""

	cdef double[:,:] List_frc = np.zeros((n_part,2))
	cdef int i, j, binloc
	cdef double rsq, d_x, d_y, xalpha

	for i in range(0, n_part):
		for j in range(0, n_part):
			if (i != j):
				d_x = float(List_pos[i][0]) - float(List_pos[j][0])
				d_y = float(List_pos[i][1]) - float(List_pos[j][1])
				d_x = make_periodic(d_x, l_box)
				d_y = make_periodic(d_y, l_box)
				rsq = dist(d_x, d_y)
				binloc = int(rsq/r_bin)

				if (rsq <= x_cut) and (binloc > binlow):
					xalpha = rsq/r_bin - binloc
					try:
						f_temp = xalpha*rforce[binloc - binlow] \
						+ (1 - xalpha)*rforce[binloc - binlow + 1]
					except IndexError:
						f_temp = rforce[binloc - binlow]
					List_frc[i][0] += d_x*f_temp
					List_frc[i][1] += d_y*f_temp
				elif (binloc <= binlow):
					xalpha = binlow - binloc + 1
					f_temp = rforce[0] + xalpha*(rforce[0] - rforce[1])
					List_frc[i][0] += d_x*f_temp
					List_frc[i][1] += d_y*f_temp

	return List_frc

cpdef double[:] gen_grBorgis(str method, double[:,:] List_pos, int n_part, int qdim, \
		double l_box, double r_bin, double x_max, double x_cut, int binmin, \
		int binlow, double[:] g_old, double[:] rforce):
	""" Compute the g(r) from the Borgis' formula, from the inner and the outer
	integration, up to a prefactor which is set in the main program with the
	appropriate parameters. qdim_max is a cutoff for the particle's counting.
	List_pos and List_frc are the lists of positions and forces for each particle"""

	cdef int i,j, k, binloc
	cdef double d_x, d_y, rsq, Delta
	cdef double[:] g_new = np.zeros(qdim)
	cdef double[:,:] List_frc = np.zeros((n_part, 2))

	for i in range(0, qdim):
		g_new[i] = g_old[i]

	List_frc = comput_frc(List_pos, n_part, x_cut, l_box, r_bin, binmin, binlow, rforce)

	for i in range(0, n_part):
		for j in range(i+1, n_part):

			d_x = float(List_pos[i][0]) - float(List_pos[j][0])
			d_y = float(List_pos[i][1]) - float(List_pos[j][1])
			d_x = make_periodic(d_x, l_box)
			d_y = make_periodic(d_y, l_box)
			rsq = dist(d_x, d_y)
			
			binloc = int(rsq/r_bin)
			#if (binloc >= binlow):
			# incremental value of the integrand in Borgis' formula
			Delta = (List_frc[i][0] - List_frc[j][0])*d_x \
			+ (List_frc[i][1] - List_frc[j][1])*d_y
			Delta /= rsq**2
			#else:
			#	Delta = 0

			if (binloc < x_max/r_bin):
				if (method == 'in'):
					#for k in range(0, binloc + 1):
#					print(i,j,binloc,round(Delta,4))
					for k in range(0, binloc + 1):
						g_new[k] += Delta
				elif (method == 'out'):
					for k in range(binloc + 1, int(x_max/r_bin)):
						try:
							g_new[k] += Delta
						except IndexError:
							print(k)
			else:
				if (method == 'in'):
					for k in range(0, int(x_cut/r_bin) + 1):
						g_new[k] += Delta

	return g_new

cpdef double get_error(double[:] array1, double[:] array2, double x_min, double x_max,
		double r_bin):
	""" Computes the total relative error between two arrays. """

	cdef int l = len(array1), i, j, k
	cdef double Delta = 0

	j = int(x_min/r_bin)
	k = int(x_max/r_bin)
	l = k - j

	for i in range(0, l):
		Delta += (array1[i] - array2[i])**2
	Delta /= l

	return Delta
