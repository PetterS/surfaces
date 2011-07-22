/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/
//
// Uses C++0x features, compile with g++/gcc option -std=c++0x
// or with Microsoft Visual Studio 2010 (preferred)
//
//

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <new>
#include <cstdlib>
#include <ctime>

#include "Petter-Color.h"

#include "3dcurv.hh"
#include "lp_constraints.hh"

#include <coin/ClpSimplex.hpp>

void curvature_3d_lp(Mesh3D& mesh, const Curv3D_Options& options, 
	const vector<double>& cell_data,
	const vector<int>& face_fixed,
	const vector<int>& cell_fixed)
{
	using namespace std;
	using namespace Petter;

	//
	// Linear programming
	//

	const auto opposing_constraints = options.opposing_constraints;
	const auto gamma = options.gamma;
	const auto lambda = options.lambda;
	const auto use_volume = options.use_volume;
	const auto xsize = options.xsize;
	const auto ysize = options.ysize;
	const auto zsize = options.zsize;

	if (options.gamma > 0) {
		statusTry("Generating face pairs...");
		mesh.generate_face_pair_list();
		statusOK();
	}

	//
	//LP matrices
	//
	size_t nEntries = 50000000; //Ususally big enough
	size_t nConstraints = mesh.nedges(); 
	size_t nVars = 2*mesh.nfaces();

	if (opposing_constraints) {
		nConstraints += mesh.nfaces();
	}

	if (gamma > 0) {
		nVars += 2*mesh.npairs();
		nConstraints +=  2 * 4*mesh.nfaces(); //Assumes every face has 4 edges
	}

	size_t cells_var_offset = nVars;
	if (use_volume) {
		nVars += mesh.ncells();
		nConstraints += mesh.nfaces();
	}

	cout << setw(9) << nVars << " variables and " << endl 
		<< setw(9) << nConstraints << " constraints." << endl;

	//Description of sparse matrix
	vector<int> rows;
	vector<int> cols;
	vector<double> values;
	rows.reserve(nEntries);
	cols.reserve(nEntries);
	values.reserve(nEntries);

	//Other LP parameters
	vector<double> rhs_lower(nConstraints, 0.0);
	vector<double> rhs_upper(nConstraints, 0.0);
	vector<double> var_lb(nVars, 0.0); // \  -
	vector<double> var_ub(nVars, 1.0); // - variables between 0 and 1
	vector<double> cost(nVars, 0.0);

	vector<int> labels(nVars, 0);

	//Adds a value to the sparse matrix
	auto add_element = [&rows,&cols,&values](size_t row, size_t col, double value) 
	{
		rows.push_back(int(row));
		cols.push_back(int(col));
		values.push_back(value);
	};
	//Changes the right-hand side of the constraints
	auto change_rhs = [&rhs_lower, &rhs_upper](size_t row, double lower, double upper) 
	{
		rhs_lower[row] = lower;
		rhs_upper[row] = upper;
	};

	//
	// Fixed faces
	//
	for (auto f=0; f < 2*mesh.nfaces() && f < face_fixed.size(); ++f) {
		if (face_fixed[f] == 0) {
			var_ub[f] = 0;
		}
		else if (face_fixed[f] == 1) {
			var_lb[f] = 1;
			labels[f] = 1; //For visualization
		}
	}

	//
	// Fixed cells and data term
	//
	if (use_volume) {
		//Fixed cells
		for (auto c=0; c<mesh.ncells() && c<cell_fixed.size(); ++c) {
			if (cell_fixed[c] == 0) {
				var_ub[cells_var_offset + c] = 0;
			}
			else if (cell_fixed[c] == 1) {
				var_lb[cells_var_offset + c] = 1;
			}
		}
		//Data term
		for (auto c=0; c<mesh.ncells() && c<cell_data.size(); ++c) {
			cost[cells_var_offset + c] = cell_data[c];
		}
	}

	//
	// Area term
	//
	statusTry("Area term");
	for (size_t fi=0;fi<mesh.nfaces();++fi) {
		double area = mesh.face_area(fi);
		double weight = lambda*area;
		cost[2*fi]   = weight;
		cost[2*fi+1] = weight;
	}
	statusOK();

	//
	// Curvature term
	//
	if (gamma > 0) {
		statusTry("Curvature term");
		for (size_t ip=0;ip<mesh.npairs();++ip) {    
			double weight = gamma * mesh.curvature_weight(ip);
			cost[2*mesh.nfaces() + 2*ip]   = weight;
			cost[2*mesh.nfaces() + 2*ip+1] = weight;
		}
		statusOK();
	}

	statusTry("Drawing data term...");
	mesh.draw_faces("x3d/fixed_faces.x3d",&labels[0]);
  mesh.draw_faces("x3d/fixed_faces.xhtml",&labels[0], true);
	statusOK();


	//
	// Linear constraints
	//
	size_t offset;
	edge_consistency_constraints(add_element,mesh,offset);

	if (opposing_constraints) {
		opposing_faces_constraints(add_element,change_rhs,mesh,offset);
	}

	if (gamma > 0) {
		surface_continuation_constraints(add_element,mesh,offset);
	}

	if (use_volume) {
		volume_face_constraints(add_element,mesh,offset, cells_var_offset);
	}


	//
	// Solve linear program 
	//
	statusTry("Loading problem...");
	CoinPackedMatrix coinMatrix(false,&rows[0],&cols[0],&values[0], CoinBigIndex(values.size()) );
	ClpSimplex lpSolver;
	lpSolver.loadProblem (coinMatrix, &var_lb[0], &var_ub[0], &cost[0], &rhs_lower[0], &rhs_upper[0]);

	lpSolver.setLogLevel(0);
	statusOK();

	statusTry("Solving LP...");
	int error = lpSolver.dual();
	if (error != 0) {
		throw runtime_error("Clp failed");
	}
	statusOK();

	cerr << "-------------------------------------" << endl;

	//
	// Obtain segmentation
	//
	parse_lp_solution(labels,cost,lpSolver.primalColumnSolution());

	int active_faces = 0;
	for (uint i=0; i < 2*mesh.nfaces(); i++) { 
		if (labels[i])
			active_faces++;
	}
	cerr << "Active faces : " << active_faces << endl;

	if (gamma > 0) {
		int active_face_pairs = 0;
		for (size_t i=2*mesh.nfaces(); i < 2*mesh.nfaces() + 2*mesh.npairs(); i++) { 
			if (labels[i])
				active_face_pairs++;
		}
		cerr << "Active face pairs : " << active_face_pairs << endl;
	}

	if (use_volume) {
		int active_cells = 0;
		for (size_t i=cells_var_offset; i < cells_var_offset + mesh.ncells(); ++i) {
			if (labels[i]) {
				active_cells++;
			}
		}
		cerr << "Active cells : " << active_cells << endl;

		statusTry("Drawing result (cells)...");
		mesh.draw("x3d/cells.x3d",&labels[cells_var_offset]);
    mesh.draw("x3d/cells.xhtml",&labels[cells_var_offset], true);
		statusOK();
	}

	statusTry("Drawing result (faces)...");
	mesh.draw_faces("x3d/faces.x3d",&labels[0]);
  mesh.draw_faces("x3d/faces.xhtml",&labels[0], true);
	statusOK();

	//
	// Verify solution for feasibility
	//
	if (use_volume) {
		statusTry("Verifying face variables...");
		for (size_t f=0; f<mesh.nfaces(); ++f) {
			if (labels[2*f] == 1 || labels[2*f+1]) {
				//Cells of this face
				const vector<uint>& cells = mesh.adjacent_cells(f);
				if (cells.size() != 2) {
					//Face at the boundary
					continue;
				}
				if (labels[cells_var_offset + cells[0]] == labels[cells_var_offset + cells[1]]) {
					throw runtime_error("face/cell infeasibility");
				}
			}
		}
		statusOK();
	}

	if (gamma > 0) {
		statusTry("Verifying face pair variables...");
		for (size_t p=0; p<mesh.npairs(); ++p) {
			if (labels[2*mesh.nfaces() + 2*p]   == 1 ||
				labels[2*mesh.nfaces() + 2*p+1] == 1) {
					//Face pair active, check that corresponding faces are also active
					uint f1 = mesh.pair(p).first_face_idx_;
					uint f2 = mesh.pair(p).second_face_idx_;
					if (labels[2*f1]==0 && labels[2*f1+1]==0) {
						throw runtime_error("face pair/face infeasibility");
					}
					if (labels[2*f2]==0 && labels[2*f2+1]==0) {
						throw runtime_error("face pair/face infeasibility");
					}
			}
		}
		statusOK();
	}
}


void parse_lp_solution(vector<int>& labels, const vector<double>& cost, const double* solution)
{
	using namespace std;
	using namespace Petter;

	//
	// Obtain segmentation
	//
	double energy = 0; //Total energy of the thresholded LP

	const double int_thresh = 0.00005;
	double largest_zero = 0;
	double smallest_one = 1;
	int fractional = 0;
	int ones = 0;

	//Function to update the above statistics
	auto check_var = [&](double var) -> int
	{
		if (var*(1-var) > int_thresh) fractional++;

		if (var < 0.5) {
			if (var > largest_zero) {
				largest_zero = var;
			}
			return 0;
		}
		else {
			if (var < smallest_one) {
				smallest_one = var;  
			}
			ones++;
			return 1;
		}
	};



	for (uint i=0; i < cost.size(); i++) { 
		labels[i] = check_var(solution[i]);
		energy += cost[i] * labels[i];
	}

	cerr << ones << " ones and " << (cost.size()-ones) << " zeros out of " << cost.size() << endl;

	if (fractional == 0) {
		cerr << GREEN << "Global optimum" << NORMAL << endl;
	}
	else {
		cerr << YELLOW << "Possibly sub-optimal solution" << NORMAL << endl;
	}
	cerr << "energy: " << energy << std::endl;
	cerr << "fractional: " << fractional << " of " << cost.size();
	cerr << std::endl;

	if (fractional > 0) {
		cerr << "    largest zero: " << setprecision(16) << largest_zero << " log err : " << log10(abs(largest_zero)) << std::endl;
		cerr << "    smallest one: " << setprecision(16) << smallest_one << " log err : " << log10(abs(smallest_one-1.0)) << std::endl;
	}
}
