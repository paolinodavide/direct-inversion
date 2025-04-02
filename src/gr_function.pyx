import numpy as np
from time import time
import cython
""" Module containing all the functions for the computation of g(r) """

cpdef gen_wt(bint log_lin, int n_correl_wt, int n_totstep, int tw_min = 1):
	""" Generates the waiting time list, cor_wt 
	It returns the new value of n_correl_wt, n_correl_wt_new as well
	The argument tw_min is optional and used only """

	cdef int k
	cdef int n_correl_wt_new
	# The two following numbers are only for the log case
	# i_max runs through the list. i_min is only incremented if 
	# list[i_min + 1] != list[i_min]
	cdef int i_min = 1
	cdef int i_max = 1
	cdef double[:] cor_wt = np.zeros(n_correl_wt)

	if log_lin:
		for k in range(1, n_correl_wt):
			cor_wt[k] = (n_totstep/tw_min)**(k/n_correl_wt)
			# erases the values that appear several times due to numerical
			# error. n_correl_wt_new is the bew value of n_correl_wt
			# Is it still needed now that python instead of FORTRAN which deals with
			# non integer ?

			# the first item of the list is trivial
			while (i_max <= n_correl_wt):
				if ( cor_wt[i_max] != cor_wt[i_min] ):
					i_min += 1
					cor_wt[i_min] = cor_wt[i_max]
					i_max += 1
				else:
					i_max += 1

		# Compute the new correlation time
		if (n_correl_wt == 0):
			n_correl_wt_new = 0
		else:
			n_correl_wt_new = i_min
	else:
	# Default case: linear time axis
		n_correl_wt_new = n_correl_wt
		for k in range(1, n_correl_wt):
			cor_wt[k] = k*n_totstep/n_correl_wt

	return n_correl_wt_new, cor_wt

cpdef list gen_pos(str prefix_file, int index_time):
	""" Generates the list of particle's positions from the appropriate data file
	Returns a list of positions List_pos """

	cdef str filename
	cdef list List_pos
	cdef str List_data
	cdef list List_lines

	# build the appropriate file name
	#filename = str(prefix_file) + str(index_time) for Ludovic's files
	filename = str(prefix_file) + str(index_time) + '.dat'
	List_pos = []
	with open(filename, "r") as f:
		List_data = f.read()
	List_lines = List_data.split('\n')
	for line in List_lines:
		if (line.split() != []) and (line.split()[0] != '#'):
			List_pos.append(line.split())

	return List_pos

cdef double lj_force(double x):
	""" Force from the Lennard-Jones potential """

	return 12.*(x)**(-14.) - 6.*(x)**(-8.)

cdef double make_periodic(double d_pos, double l_box):
	""" Returns a position corrected to account for the periodic boundary
	conditions in a square box of size l_box """

	if (d_pos >= l_box/2.):
		return d_pos - l_box
	elif (d_pos <= -l_box/2.):
		return d_pos + l_box
	else:
		return d_pos

cdef double dist(double d_x,double d_y):
	""" returns the distance between two points from the coordinate differences"""

	return np.sqrt(d_x**2 + d_y**2)

cpdef list gen_force(list List_pos, int n_part, float x_cut, double l_box):
	""" Compute the list of forces acting on the particles from their positions """

	cdef list List_frc = []
	cdef int i,j
	cdef double x, y, f_x, f_y, rsq, d_x, d_y

	for i in range(0, n_part):
		x = float(List_pos[i][0])
		y = float(List_pos[i][1])
		f_x = 0
		f_y = 0

		for j in range(0, n_part):
			if (i != j):
				d_x = x - float(List_pos[j][0])
				d_y = y - float(List_pos[j][1])
				d_x = make_periodic(d_x, l_box)
				d_y = make_periodic(d_y, l_box)
				rsq = dist(d_x, d_y)

				if (rsq <= x_cut):
					f_temp = lj_force(rsq)
					f_x = f_x + d_x*f_temp
					f_y = f_y + d_y*f_temp
		List_frc.append([f_x, f_y])

	return List_frc

cpdef gen_grBorgis(list List_pos, list List_frc, int n_part, int qdim_max, double l_box, float r_bin, double[:] g_in_old, double[:] g_out_old):
	""" Compute the g(r) from the Borgis' formula, from the inner and the outer
	integration, up to a prefactor which is set in the main program with the
	appropriate parameters. qdim_max is a cutoff for the particle's counting.
	List_pos and List_frc are the lists of positions and forces for each particle"""

#	cdef double[:] g_in_new = g_in_old
#	cdef double[:] g_out_new = g_out_old
	cdef int i,j, k, binloc
	cdef int qdim = len(g_in_old)
	cdef double x, y, d_x, d_y, rsq, Delta
	cdef double[:] g_in_new = np.zeros(qdim)
	cdef double[:] g_out_new = np.zeros(qdim)

	for i in range(0,qdim):
		g_in_new[i] = g_in_old[i]
		g_out_new[i] = g_out_old[i]

	for i in range(0, n_part):
		x = float(List_pos[i][0])
		y = float(List_pos[i][1])

		for j in range(i+1, n_part):

			d_x = x - float(List_pos[j][0])
			d_y = y - float(List_pos[j][1])
			d_x = make_periodic(d_x, l_box)
			d_y = make_periodic(d_y, l_box)
			rsq = dist(d_x, d_y)
			
			# incremental value of the integrand in Borgis' formula
			Delta = (List_frc[i][0] - List_frc[j][0])*d_x \
			+ (List_frc[i][1] - List_frc[j][1])*d_y
			Delta /= rsq**2

			binloc = int(rsq/r_bin)

			if (binloc <= qdim_max):
				for k in range(0, binloc):
					g_in_new[k] += Delta
				for k in range(binloc, qdim_max):
					g_out_new[k] += Delta
			else:
				for k in range(0, qdim_max):
					g_in_new[k] += Delta

	return g_in_new, g_out_new

cpdef double[:] gen_grHisto(list List_pos, int n_part, int qdim_max, double l_box, float r_bin, double[:] g_old):
	""" Compute the g(r) from the histograms, 
	up to a prefactor which is set in the main program with the
	appropriate parameters. qdim_max is a cutoff for particle's counting.
	List_pos is the list of positions for each particle"""

	cdef double[:] g_new = g_old
	cdef int i,j, k, binloc
	cdef double x, y, d_x, d_y, rsq, Delta

	for i in range(0, n_part):
		x = float(List_pos[i][0])
		y = float(List_pos[i][1])

		for j in range(i+1, n_part):

			d_x = x - float(List_pos[j][0])
			d_y = y - float(List_pos[j][1])
			d_x = make_periodic(d_x, l_box)
			d_y = make_periodic(d_y, l_box)
			rsq = dist(d_x, d_y)

			# distance in units of r_bin
			binloc = int(rsq/r_bin)

			Delta = 1/rsq/r_bin

			if (binloc <= qdim_max):
				g_new[binloc] += Delta

	return g_new

cpdef double[:] patch_gr_end(double[:] g_old, int qdim_max, double prefactor, double r_bin):
	""" To be applied to the out function of the Borgis formula.
	Corrects the end value of g(r).
	Carefull the data is not normalized !"""

	cdef double[:] g_new = g_old
	cdef int n = qdim_max, l_patch = int(n/500)
	cdef int i
	cdef double delta_1 = 0., alpha = 5.
	cdef double beta

	# Computes the error on the end value. Averages on 0.5% of the points
	for i in range(n - l_patch, n):
		delta_1 += g_old[i]/prefactor - 1
	delta_1 /= l_patch

	for i in range(1,n):
		# Careful with i=0
		# Too violent in the low r region
		# g_new[i] = g_old[i]/(1 + delta_1)
		beta = - delta_1/(1. + delta_1)
		g_new[i] = g_old[i]*(1. + beta*np.exp(- alpha/(i*r_bin)))

	return g_new

cpdef double[:] patch_gr_inout(double[:] g_in_old, double[:] g_out_old, int qdim_max, double prefactor, double r_bin):
	""" To be applied to the out function of the Borgis formula.
	Corrects the end value of g(r) by interpolating the in and out formulas"""

	cdef double[:] g_new = np.zeros(qdim_max)
	cdef int i
	cdef double alpha = 3.

	for i in range(1, qdim_max):
		# Careful with i=0
		# Too violent in the low r region
		# g_new[i] = g_old[i]/(1 + delta_1)
		g_new[i] = (1 - g_in_old[i]/prefactor)*(np.exp( - (alpha/(i*r_bin))**4 )) \
			+ g_out_old[i]/prefactor*(1 - np.exp( - (alpha/(i*r_bin))**4 ))

	return g_new

cpdef double get_delta_0(list list_gr, int l):

	""" Returns the offset at the beginning of the curve of gr computed with the
	Borgis formula with inside integration. l is the length over which
	the value is averaged. """

	cdef int n = 0
	cdef double delta_0 = 0

	while (n < l):
		delta_0 += list_gr[n]
		n += 1
	delta_0 /= n

	return delta_0

cpdef double get_delta_1(list list_gr, int l):

	""" Returns the offset at the end of the curve of gr computed with the
	Borgis formula with outside integration. l is the length over which
	the value is averaged. """

	cdef int n = 0
	cdef double delta_1 = 0
	cdef int length = len(list_gr)

	while (n < l):
		#Carefull, last value is 0
		delta_1 += list_gr[length - n - 2]
		n += 1
	delta_1 /= n

	return delta_1
