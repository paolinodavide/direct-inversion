#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <omp.h>
#include <fstream>
#include <jsoncpp/json/json.h>
#include <vector>
#include <map>
#include <cmath>

using namespace std;

#define THREAD_NUM 8

// For compilation see the associated .sh file

double dist_cpp(double dx, double dy) {

	return sqrt( pow(dx,2) + pow(dy,2));
}

double make_periodic_cpp(double d_pos, double l_box) {
	
	if (d_pos >= l_box/2.){
		return d_pos - l_box;
	} else if (d_pos <= -l_box/2.) {
		return d_pos + l_box;
	} else {
		return d_pos;
	}

}

vector<vector<double>> comput_frc_cpp (vector<vector<double>> const List_pos,
		int const n_part, double const x_cut, double const l_box,
		double const r_bin, int const binmin, int const binlow,
		//pybind11::array_t<double> const rforce) {
		vector <double> const rforce_vector) {
	 
	vector<vector<double>> List_frc;

	//double *rforce_pointer = (double*) rforce.request().ptr;

	int i, j, binloc;
	double  rsq, d_x, d_y, xalpha, f_temp;

	for (i=0; i < n_part; i++) {
		List_frc.push_back({0., 0.});
	}

	for (i=0; i < n_part; i++) {
		for (j=0; j < n_part; j++) {
			if (i != j) {
				d_x = List_pos[i][0] - List_pos[j][0];
				d_y = List_pos[i][1] - List_pos[j][1];
				d_x = make_periodic_cpp(d_x, l_box);
				d_y = make_periodic_cpp(d_y, l_box);
				rsq = dist_cpp(d_x, d_y);
				binloc = (int) (rsq/r_bin);

				if ((rsq <= x_cut) and (binloc > binlow)) {
					xalpha = rsq/r_bin - binloc;
					//f_temp = xalpha*rforce_pointer[binloc - binlow]
					//	+ (1 - xalpha)*rforce_pointer[binloc - binlow + 1];
					f_temp = xalpha*rforce_vector[binloc - binlow]
						+ (1 - xalpha)*rforce_vector[binloc - binlow + 1];
					List_frc[i][0] = List_frc[i][0] + d_x*f_temp;
					List_frc[i][1] = List_frc[i][1] + d_y*f_temp;
				} else if (binloc <= binlow){
					xalpha = binlow - binloc + 1;
					//f_temp = rforce_pointer[0] +
					//	xalpha*(rforce_pointer[0] - rforce_pointer[1]);
					f_temp = rforce_vector[0] +
						xalpha*(rforce_vector[0] - rforce_vector[1]);
					List_frc[i][0] = List_frc[i][0] + d_x*f_temp;
					List_frc[i][1] = List_frc[i][1] + d_y*f_temp;
				}
			}
		}
	}

	//delete rforce_pointer;

	return List_frc;
}

vector<double> gen_gtmp (int const n_part, int const qdim,
	       	vector<vector<double>> const List_pos,
		vector<vector<double>> const List_frc,
	       	double const x_max, double const r_bin,
		double const x_cut, double const l_box) {

	int i, j, k, binloc;
	double rsq, Delta, d_x, d_y;
	vector<double> g_tmp;

	for (k=0; k < qdim; k++) {
		g_tmp.push_back(0.);
	}

	for (i=0; i < n_part; i++) {
		for (j=i+1; j < n_part; j++) {

			d_x = List_pos[i][0]
				- List_pos[j][0];
			d_y = List_pos[i][1] 
				- List_pos[j][1];
			d_x = make_periodic_cpp(d_x, l_box);
			d_y = make_periodic_cpp(d_y, l_box);
			rsq = dist_cpp(d_x, d_y);

			binloc = (int) (rsq/r_bin);
			Delta = (List_frc[i][0] 
					- List_frc[j][0])*d_x
				+ (List_frc[i][1] 
						- List_frc[j][1])*d_y;
			Delta /= pow(rsq, 2);

			if (binloc < (int) (x_max/r_bin)) {
				for (k=0; k < binloc + 1; k++){
					g_tmp[k] += Delta;
				}
			} else {
				for (k=0; k< (int) (x_cut/r_bin); k++) {
					g_tmp[k] += Delta;
				}
			}
		}
	}
	return g_tmp;
}

vector <double> gen_force (pybind11::array_t<double> const rforce, int const pot_length) {

	double *rforce_pointer = (double*) rforce.request().ptr;
	vector <double> force_vector;

	for (int k=0; k < pot_length; k++) {
		force_vector.push_back(rforce_pointer[k]);
	}

	return force_vector;
}

pybind11::array_t<double> gen_grBorgis_cpp(int const n_part, int const qdim,
		double const l_box, double const r_bin, double const x_max,
		double const x_cut, int const binmin, int const binlow,
		pybind11::array_t<double> const rforce) {

	pybind11::array_t<double> g_new({qdim});
	double *g_pointer = (double*) g_new.request().ptr;
	for (int k=0; k < qdim; k++) {
		g_pointer[k] = 0.;
	}

	int i, k, wt, N;
	map <int, vector<vector<double>>> List_frc;
	map <int, vector<vector<double>>> List_pos;
	map <int, vector<double>> g_tmp;

	int bincut = (int) (x_cut/r_bin);
	int pot_length = bincut - binlow;
	vector <double> rforce_vector;
	rforce_vector = gen_force(rforce, pot_length);

	//Initialise g(r)
	//for (k=0; k < qdim; k++) {
	//	g_new.push_back(0.);
	//}

	//Read positions
	Json::Value positions;
	std::ifstream positions_file("positions.json", std::ifstream::binary);
	positions_file >> positions;

	vector<int> List_wt;
	int number_wt = positions["wt"].size();
	for (i=0; i < number_wt; i++){
		List_wt.push_back( positions["wt"][i].asInt() );
	}
	for (N=0; N < number_wt; N++){
		wt = List_wt[N];
		for (i=0; i < n_part; i++){
			List_pos[wt].push_back({
					positions[to_string(wt)][i][0].asDouble(),
					positions[to_string(wt)][i][1].asDouble()
					});
		}
	}

	for (N = 0; N < number_wt; N++){
		wt = List_wt[N];
		//cout << wt << endl;
		for (k=0; k < qdim; k++) {
			g_tmp[wt].push_back(0.);
		}
	}

	#pragma omp parallel for num_threads(THREAD_NUM)
	for (N = 0; N < number_wt; N++){

		//List_frc[List_wt[N]] = comput_frc_cpp(List_pos[List_wt[N]], n_part,
		//		x_cut, l_box, r_bin, binmin, binlow, rforce);

		List_frc[List_wt[N]] = comput_frc_cpp(List_pos[List_wt[N]], n_part,
				x_cut, l_box, r_bin, binmin, binlow, rforce_vector);

		g_tmp[List_wt[N]] = gen_gtmp(n_part, qdim, List_pos[List_wt[N]],
				List_frc[List_wt[N]], x_max, r_bin, x_cut, l_box);
	}

	#pragma omp critical
	for (N=0; N < number_wt; N++) {
		wt = List_wt[N];
		for (k=0; k < qdim; k++) {
			g_pointer[k] += g_tmp[wt][k];
		}
	}

	//delete g_pointer;

	return g_new;
}

PYBIND11_MODULE(gr_pair_parallel, m) {
	m.doc() = "PyBind11 library.";

	m.def("gen_grBorgis_cpp", &gen_grBorgis_cpp,
			"");
}
