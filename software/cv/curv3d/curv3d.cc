/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/
//
// Uses C++0x features, compile with g++/gcc option -std=c++0x
// or with Microsoft Visual Studio 2010 (preferred)
//
//

//Get access to standard functions
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <new>
#include <cstdlib>
#include <ctime>
#include <map>
#include <set>

//Mesh3D class 
#include "mesh3d.hh"
//Functions to solve 3D problems
#include "3dcurv.hh"
//Write in color to console
#include "Petter-Color.h"

#include <armadillo>


//Simple routine for conversion of strings
//used for the command line
template <typename T>
T convert_string(const std::string s) 
{
	std::istringstream is(s);
	T result;
	is >> result;
	if (!is) {
		throw std::runtime_error("conversion of \"" + s + "\" failed"); 
	}
	return result;
}

//Test routine for mean curvature
void test_curvature(std::map<std::string,std::string>& cmd_line)
{
	using namespace std;
	using namespace Petter;
	using namespace arma;

	cout << WHITE << "Test of mean curvature." << NORMAL << endl;
	cout << "-------------------------------------" << endl;

	int n = convert_string<int>(cmd_line["-res"]);
	double r = convert_string<double>(cmd_line["-r"]);
	double h = convert_string<double>(cmd_line["-h"]);

	cout << "r = " << r << endl;
	cout << "h = " << h << endl;

	statusTry("Creating theta,phi...");

	statusOK();


	mat x(n,n);
	mat y(n,n);
	mat z(n,n);
	double true_energy;
	double true_area;

	if (cmd_line.find("-cylinder") != cmd_line.end()) {
		true_area = 2.0*math::pi()*r*h;
		double k1 = 0.0;
		double k2 = 1/r;
		double H = 0.5*(k1+k2);
		double H2 = H*H;
		double K = k1*k2;
		true_energy = true_area * (k1*k1 + k2*k2);
		statusTry("Creating cylinder...");
		mat phi = linspace<mat>(0,2*math::pi(), n);
		mat zz   = linspace<mat>(0,h, n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; ++j) {
				x(i,j) = r*cos(phi(i));
				y(i,j) = r*sin(phi(i));
				z(i,j) = zz(j);
			}
		}
	}
	else if (cmd_line.find("-catenoid") != cmd_line.end()) {
		true_area = 0;
		true_energy = 0;

		statusTry("Creating catenoid...");
		mat phi = linspace<mat>(0,2*math::pi(), n);
		mat zz   = linspace<mat>(-h/2,h/2, n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; ++j) {
				double r = cosh(zz(j));
				x(i,j) = r*cos(phi(i));
				y(i,j) = r*sin(phi(i));
				z(i,j) = zz(j);
			}
		}


		true_area = 2*math::pi()* (h/2 + sinh(h)/2);
		true_area = 0;
		int N = 20000;
		zz = linspace<mat>(-h/2,h/2, N);
		mat area(N, 1);
		mat energy(N, 1);
		for (int i=0;i<N;++i) {
			double u = zz(i);
			double r = cosh(u);
			double rp = sinh(u);

			double k1 = -1/(cosh(u)*cosh(u));
			double k2 = 1/(cosh(u)*cosh(u));
			double H = 0.5*(k1+k2);
			double H2 = H*H;
			double K = k1*k2;

			double c =  k1*k1+k2*k2;

			area(i) = 2*math::pi() * r * sqrt(1+rp*rp);
			energy(i) = area(i) * c;
		}
		true_area = area(0) + area(N-1);
		true_energy = energy(0) + energy(N-1);
		for (int i=1;i<N-1;++i) {
			true_area   += 2*area(i);
			true_energy += 2*energy(i);
		}
		true_area   *= h/(2*N);
		true_energy *= h/(2*N);
	}
	else if (cmd_line.find("-sphere") != cmd_line.end()) {
		true_area = 4.0*math::pi()*r*r;
		double k1 = 1/r;
		double k2 = 1/r;
		double H = 0.5*(k1+k2);
		double H2 = H*H;
		double K = k1*k2;
		true_energy = true_area * (k1*k1 + k2*k2);
		statusTry("Creating sphere...");
		mat theta = linspace<mat>(0,math::pi(), n);
		mat phi   = linspace<mat>(0,2*math::pi(), n);
		for (int i=0; i<n; i++) {
			for (int j=0; j<n; ++j) {
				x(i,j) = r*sin(theta(i))*cos(phi(j));
				y(i,j) = r*sin(theta(i))*sin(phi(j));
				z(i,j) = r*cos(theta(i));
			}
		}
	}
	else {
		throw runtime_error("No shape specified");
	}
	statusOK();

	Mesh3D mesh;
	statusTry("Adding points to mesh...");
	Mat<uint> p(n,n);
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; ++j) {
			p(i,j) = mesh.find_or_add_point(x(i,j),y(i,j),z(i,j));
			if (i>0 && j>0) {

				//Triangles
				/*
				set<uint> point_set;
				vector<uint> points;
				points.push_back(p(i,j));
				point_set.insert(p(i,j));
				uint p1 = p(i-1,j);
				if (point_set.find(p1) == point_set.end()) {
				points.push_back(p1);
				point_set.insert(p1);
				}
				p1 = p(i-1,j-1);
				if (point_set.find(p1) == point_set.end()) {
				points.push_back(p1);
				point_set.insert(p1);
				}

				vector<uint> edges;
				for (int ind=1;ind<points.size();++ind) {
				uint l1 = mesh.find_or_add_edge(points[ind],points[ind-1]);
				edges.push_back(l1);
				}
				uint l1 =  mesh.find_or_add_edge(points[points.size()-1],points[0]);
				edges.push_back(l1);
				if (edges.size() == 3) {
				mesh.find_or_add_face(edges);
				}

				point_set.clear();
				points.clear();
				points.push_back(p(i,j));
				point_set.insert(p(i,j));
				p1 = p(i,j-1);
				if (point_set.find(p1) == point_set.end()) {
				points.push_back(p1);
				point_set.insert(p1);
				}
				p1 = p(i-1,j-1);
				if (point_set.find(p1) == point_set.end()) {
				points.push_back(p1);
				point_set.insert(p1);
				}

				edges.clear();
				for (int ind=1;ind<points.size();++ind) {
				uint l1 = mesh.find_or_add_edge(points[ind],points[ind-1]);
				edges.push_back(l1);
				}
				l1 =  mesh.find_or_add_edge(points[points.size()-1],points[0]);
				edges.push_back(l1);
				if (edges.size() == 3) {
				mesh.find_or_add_face(edges);
				}
				*/


				//Squares

				set<uint> point_set;
				vector<uint> points;
				points.push_back(p(i,j));
				point_set.insert(p(i,j));
				uint p1 = p(i-1,j);
				if (point_set.find(p1) == point_set.end()) {
					points.push_back(p1);
					point_set.insert(p1);
				}
				p1 = p(i-1,j-1);
				if (point_set.find(p1) == point_set.end()) {
					points.push_back(p1);
					point_set.insert(p1);
				}
				p1 = p(i,j-1);
				if (point_set.find(p1) == point_set.end()) {
					points.push_back(p1);
					point_set.insert(p1);
				}

				vector<uint> edges;
				for (int ind=1;ind<points.size();++ind) {
					uint l1 = mesh.find_or_add_edge(points[ind],points[ind-1]);
					edges.push_back(l1);
				}
				uint l1 =  mesh.find_or_add_edge(points[points.size()-1],points[0]);
				edges.push_back(l1);
				if (edges.size() == 3 || edges.size() == 4) {
					mesh.find_or_add_face(edges);
				}
				//*/


			}
		}
	}
	statusOK();

	double area = 0;
	for (size_t i=0;i<mesh.nfaces();++i) {
		area += mesh.face_area(i);
	}

	cout << "True area : " << true_area << endl;
	cout << "Computed area : " << area << endl;


	statusTry("Building face pairs...");
	mesh.generate_face_pair_list();
	statusOK();
	cout << mesh.npairs() << " pairs." << endl;

	statusTry("Calculating curvature...");
	double energy = 0;
	for (size_t i=0;i<mesh.npairs();++i) {
		energy += mesh.curvature_weight(i);
	}
	statusOK();

	cout << "Integral of curvature should be " << true_energy << endl;
	cout << "Energy : " << energy << endl;
	cout << "Factor : " << energy/true_energy << endl;

	statusTry("Saving mesh...");
	vector<int> face_labels(2*mesh.nfaces(), 1);
	for (int i=0;i<mesh.nfaces();++i) {
		vec3 n = mesh.face_normal(i);
		Mesh3DPoint p = mesh.face_center(i);
		vec3 v;
		v[0] = p.x_;
		v[1] = p.y_;
		v[2] = p.z_;
		if (dot(n,v) > 0) {
			face_labels[2*i] = 1;
		}
		else {
			face_labels[2*i+1] = 1;
		}
	}
	mesh.draw_faces("out/sphere.x3d", &face_labels[0]);
	statusOK();
}

int main_program(int argc, char** argv) {
	using namespace std;
	using namespace Petter;
	using namespace arma;

	if ( (argc == 2 && string(argv[1])=="-h") || argc==1 ) {

		cout << "Armadillo version: " << arma_version::as_string() << endl;
		if (arma_config::debug) {
			cerr << "Armadillo is in debug mode." << endl;
		}
		else {
			cerr << "Armadillo is not in debug mode." << endl;
		}

		cerr << "USAGE: " << argv[0] << std::endl
			<< "  -size <int>      :  size of mesh" << endl
			<< "  -lambda <double> :  area weight" << endl
			<< "  -gamma <double>  :  curvature weight" << endl
			<< "  -cells           :  use volume variables " << endl
			<< "  -cross           :  use cross data term (instead of the face circles)" << endl
			<< "  -cubic           :  use cubes instead of tetrahedrons " << endl
			<< "  -subdivide       :  subdivides cells further " << endl
			<< "  -qpbo            :  use QPBO instead of linear programming " << endl
			<< endl
			;
		return 0;
	}

	map<string,string> cmd_line;

	//Default parameters
	cmd_line["-gamma"] = "1";
	cmd_line["-lambda"] = "0";
	cmd_line["-res"] = "30";
	cmd_line["-r"] = "10";
	cmd_line["-h"] = "10";

	//Read command line into map
	for (size_t i=1;i<argc;++i) {
		string cmd = argv[i];
		if (i<argc-1 && argv[i+1][0]!='-') {
			cmd_line[cmd] = argv[i+1];
			++i;
		}
		else {
			cmd_line[cmd] = "";
		}
	}

	//If -test is in the command line, run some tests and exit
	if (cmd_line.find("-test") != cmd_line.end()) {
		test_curvature(cmd_line);
		return 0;
	}


	cout << WHITE << "3D mean curvature." << NORMAL << endl;
	cout << "-------------------------------------" << endl;


	int size = convert_string<int>(cmd_line["-size"]);
	Mesh3D mesh;
	int xsize = size;
	int ysize = size;
	int zsize = size;
	if (cmd_line.find("-cross")==cmd_line.end()) {
		zsize = size/2 - 5;
	}

	//
	// Options
	//
	Curv3D_Options options;
	options.xsize = xsize;
	options.ysize = ysize;
	options.zsize = zsize;
	options.gamma = convert_string<double>(cmd_line["-gamma"]);
	options.lambda = convert_string<double>(cmd_line["-lambda"]);
	options.use_volume = cmd_line.find("-cells") != cmd_line.end();

	cout << "lambda=" << options.lambda << " gamma=" << options.gamma << endl;

	cout << "Mesh size : " << xsize << "x" << ysize << "x" << zsize << endl;
	if (cmd_line.find("-cubic") != cmd_line.end()) {
		statusTry("Building cubic mesh...");
		generate_3D_cubic_mesh(xsize,ysize,zsize, mesh);
	}
	else {
		statusTry("Building mesh...");
		generate_3D_mesh(xsize,ysize,zsize, mesh);
	}
	statusOK();

	cout << setw(9) << mesh.npoint() << " points, " << endl << setw(9) << mesh.nedges() << " edges, " 
		<< endl << setw(9) << mesh.nfaces() << " faces and " << endl << setw(9) << mesh.ncells() 
		<< " cells" << endl;

	if (cmd_line.find("-points") != cmd_line.end()) {
		statusTry("Subdividing mesh cells...");
		mesh.subdivide_cells();
		statusOK();
		cout << setw(9) << mesh.npoint() << " points, " << endl << setw(9) << mesh.nedges() << " edges, " 
			<< endl << setw(9) << mesh.nfaces() << " faces and " << endl << setw(9) << mesh.ncells() 
			<< " cells" << endl;
	}

	if (size < 3) {
		statusTry("Drawing mesh...");
		mesh.draw("3Dmesh.x3d");
		statusOK();

		//Too small mesh to do anything with, but the image
		//is useful as illustration
		return 0;
	}


	// Fixed variables
	vector<int> face_fixed(2*mesh.nfaces(), -1);
	vector<int> cell_fixed(mesh.ncells(), -1);
	// Cell data term
	vector<double> cell_data(mesh.ncells(), 0.0);

	//
	// Fixed data term
	//
	statusTry("Fixed labels");
	if (cmd_line.find("-cross") == cmd_line.end()) {
		//
		// We fix face variables
		//
		double xc = double(xsize)/2 + 0.5;
		double yc = double(ysize)/2 + 0.5;
		double zc = double(zsize)/2 + 0.5;
		double marg = 1;

		for (size_t fi=0;fi<mesh.nfaces();++fi) {
			auto center = mesh.face_center(fi);
			double x = center.x_;
			double y = center.y_;
			double z = center.z_;
			auto normal = mesh.face_normal(fi);

			bool in_cylinder = (x-xc)*(x-xc) + (y-yc)*(y-yc) <= size*size/6;

			if ( abs(z) < 1e-6 /* z == 0 */                  &&
				in_cylinder)
			{
				if (normal[2] < 0 /* normal pointing outwards */ ) {
					face_fixed[2*fi] = 1;
				}
				else {
					face_fixed[2*fi+1] = 1;
				}
			}

			else if ( abs(z-zsize) < 1e-6 /* z == zsize */        &&
				in_cylinder )
			{
				if (normal[2] > 0 /* normal pointing outwards */ ) {
					face_fixed[2*fi] = 1;
				}
				else {
					face_fixed[2*fi+1] = 1;
				}
			}

		}
	}
	else {
		//
		// We fix cell variables
		//
		if (!options.use_volume) {
			throw runtime_error("Cross requires volume variables");
		}

		double xc = double(xsize)/2 + 0.5;
		double yc = double(ysize)/2 + 0.5;
		double zc = double(zsize)/2 + 0.5;
		double marg = 1;

		for (size_t ic=0;ic<mesh.ncells();++ic) {
			auto center = mesh.cell_center(ic);
			double x = center.x_;
			double y = center.y_;
			double z = center.z_;
			if ( x>marg && x<xsize-marg && y>marg && y<ysize-marg && z>marg && z<zsize-marg &&
				// (abs(x-xc)<0.5 || abs(y-yc)<0.5 || abs(z-zc)<0.5) 
				( (abs(x-xc)<0.5 || abs(y-yc)<0.5) &&
				(abs(x-xc)<0.5 || abs(z-zc)<0.5) &&
				(abs(z-zc)<0.5 || abs(y-yc)<0.5) )
				) {
					//
					// Center cross
					//
					cell_fixed[ic] = 1;
			}
			else {
				cell_data[ic] = 0.001;
			}
		}

		//
		//Background at the edges
		//
		for (size_t fi=0;fi<mesh.nfaces();++fi) {
			const std::vector<uint>& cells = mesh.adjacent_cells(fi);
			if (cells.size() == 1) {
				//This face is at an edge
				cell_fixed[cells[0]] = 0;
			}
		}
	}
	statusOK();



	//
	// Solve problem
	//
	cerr << endl;
	if (cmd_line.find("-qpbo") == cmd_line.end() ) {
		// Linear programming
		cerr << WHITE << "Using linear programming" << NORMAL << endl;
		cerr << "-------------------------------------" << endl;
		curvature_3d_lp(mesh,options,cell_data,face_fixed,cell_fixed);
	}
	else {
		// QPBO
		cerr << WHITE << "Using QPBO" << NORMAL << endl;
		cerr << "-------------------------------------" << endl;
		curvature_3d_qpbo(mesh,options,cell_data,cell_fixed);
	}

	return 0;
}


int main(int argc, char** argv) 
{
	using namespace std;
	using namespace Petter;
	try {
		return main_program(argc,argv);
	}
	catch (runtime_error& e) {
		statusFailed();
		cerr << "Run-time error : " << e.what() << endl;
	}
	catch (bad_alloc&) {
		statusFailed();
		cerr << RED << "Out of memory." << endl;
	}
	catch (exception& e) {
		cerr << RED << "Exception : " << e.what() << NORMAL << endl;
	}
	catch (...) {
		cerr << RED << "Unknown error" << NORMAL << endl;
	}
	return 1;
}
