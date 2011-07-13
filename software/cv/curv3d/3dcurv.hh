/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/

#ifndef CURV3D_HEADER_DEF
#define CURV3D_HEADER_DEF

#include <vector>
using std::vector;

#include "mesh3d.hh"


struct Curv3D_Options
{
	double lambda;
	double gamma;
	int xsize,ysize,zsize;

	bool opposing_constraints;
	bool use_volume;

	Curv3D_Options()
	{
		lambda = 1;
		gamma = 0;
		xsize = 16;
		ysize = 16;
		zsize = 16;

		opposing_constraints = true;
		use_volume = true;
	}
};

void curvature_3d_lp(Mesh3D& mesh, const Curv3D_Options& options, 
	const vector<double>& cell_data,
	const vector<int>& face_fixed,
	const vector<int>& cell_fixed);

void parse_lp_solution(vector<int>& labels, const vector<double>& cost, const double* solution);

void curvature_3d_qpbo(Mesh3D& mesh, const Curv3D_Options& options, 
	const vector<double>& cell_data,
	const vector<int>& cell_fixed);


#endif