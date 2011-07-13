/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/
//
// Uses C++0x features, compile with g++ option -std=c++0x
// or with Microsoft Visual Studio 2010 (preferred)
//
//

#include <iostream>
#include <stdexcept>
#include <new>
#include <cstdlib>
#include <ctime>
#include <typeinfo>

#include "Petter-Color.h"

#include "3dcurv.hh"

#include "QPBO.h" 
#include "HOCR.h"

void err_function(char * err)
{
	throw std::runtime_error(err);
}

void curvature_3d_qpbo(Mesh3D& mesh, const Curv3D_Options& options, 
	const vector<double>& cell_data,
	const vector<int>& cell_fixed)
{
	using namespace std;
	using namespace Petter;

	typedef long long real;

	auto& tinfo = typeid(real);
	cout << "QPBO type " << tinfo.name()  << endl;

	if (!options.use_volume) {
		throw runtime_error("QPBO requires cell variables");
	}

	auto xsize = options.xsize;
	auto ysize = options.ysize;
	auto zsize = options.zsize;

	int required_nodes = int(mesh.ncells());
	int required_edges = 1000;

	statusTry("Generating face pairs...");
	mesh.generate_face_pair_list();
	statusOK();

	statusTry("Creating QPBO object...");
	QPBO<real> qpbo(required_nodes,required_edges, err_function);
	statusOK();
	statusTry("Creating HOCR object...");
	HOCR<real,4,QPBO<real> > hocr(qpbo);
	statusOK();
	statusTry("Adding nodes...");
	hocr.AddNode(int(mesh.ncells()));
	statusOK();

	//
	// Data term
	//
	statusTry("Data term");

	for (size_t ic=0;ic<mesh.ncells();++ic) {
		if (cell_fixed.at(ic) == 0) {
			hocr.AddUnaryTerm(int(ic),0,real(1e9));
		}
		else if (cell_fixed.at(ic) == 1) {
			hocr.AddUnaryTerm(int(ic), real(1e9),0);
		}
		else {
			hocr.AddUnaryTerm(int(ic), 0, real(cell_data.at(ic)) );
		}
	}


	statusOK();


	//
	// Area regularizing term
	//
	statusTry("Area terms...");
	int nterms = 0;
	for (size_t fi=0;fi<mesh.nfaces();++fi) {
		real cost = real( options.lambda*mesh.face_area(fi) );

		const std::vector<uint>& cells = mesh.adjacent_cells(fi);
		if (cells.size() == 2) {
			hocr.AddPairwiseTerm(cells[0],cells[1],0,cost,cost,0);
			nterms++;
		}
	}
	statusOK();
	cout << "Number of area terms : " << nterms << endl;

	//
	// Curvature regularizing term
	//
	int threecliques = 0;
	int fourcliques  = 0;

	if (options.gamma > 0) {
		statusTry("Curvature terms...");
		for (size_t j=0; j<mesh.npairs(); ++j) {
			uint first  = mesh.pair(j).first_face_idx_;
			uint second = mesh.pair(j).second_face_idx_;
			uint e  = mesh.pair(j).common_edge_idx_;

			if (mesh.adjacent_cells(first).size() != 2 ||
				mesh.adjacent_cells(second).size() != 2) continue;

			//Four adjacent cells
			uint f1 = mesh.adjacent_cells(first)[0];
			uint f2 = mesh.adjacent_cells(first)[1];
			uint s1 = mesh.adjacent_cells(second)[0];
			uint s2 = mesh.adjacent_cells(second)[1];

			//Face normals
			auto fn = mesh.face_normal(first);
			auto sn = mesh.face_normal(second);
			//Face centers
			Mesh3DPoint fc = mesh.face_center(first);
			Mesh3DPoint sc = mesh.face_center(second);

			//Make them point toward the same side
			//Are vectors (close to) parallel?
			if ( norm(cross(fn,sn), 2) < 1e-4 ) {
				if ( dot(fn,sn) < 0 ) {
					//Swap orientation
					fn = -fn;
				}
			}
			else {
				//Interior point
				Mesh3DPoint interior = fc;
				interior += sc;
				interior /= 2;
				//Vectors to interior
				auto intf = interior - fc;
				auto ints = interior - sc;

				//Is the normal pointing the "correct way"?
				if ( dot(intf,fn) < 0) {
					fn = -fn;
				}
				if ( dot(ints,sn) < 0) {
					sn = -sn;
				}
			}


			//Cell centers
			Mesh3DPoint f1c = mesh.cell_center(f1);
			Mesh3DPoint f2c = mesh.cell_center(f2);
			Mesh3DPoint s1c = mesh.cell_center(s1);
			Mesh3DPoint s2c = mesh.cell_center(s2);

			//Vectors from the face centers to the cell centers
			auto f1n = f1c - fc;
			auto f2n = f2c - fc;
			auto s1n = s1c - sc;
			auto s2n = s2c - sc;

			//f1 on the same side as s1?
			//The projections on the normals have different sign
			if ( dot(f1n,fn) * dot(s1n,sn) < 0) {
				swap(f1,f2);
			}


			real weight = real( options.gamma * mesh.curvature_weight(j) );

			if (abs(weight) <= 1e-9) {
				continue; //No reason to add this 
			}

			if (f1 == s1) {
				//Three-clique
				threecliques++;
				//000
				//001
				//010
				//011 -- cost
				//100 -- cost
				//101
				//110
				//111
				real vals[8] = {0,0,0,weight, weight,0,0,0};
				int vars[3] = {f1, f2, s2};
				hocr.AddHigherTerm(3, vars, vals);
			}
			else if (f2 == s2) {
				//Three-clique
				threecliques++;
				//011 -- cost
				//100 -- cost
				real vals[8] = {0,0,0,weight, weight,0,0,0};
				int vars[3] = {f2, f1, s1};
				hocr.AddHigherTerm(3, vars, vals);
			}
			else if (f1 == s2 || f2 == s1) {
				//The swap above should have prevented this
				throw runtime_error("f1 == s2 || f2 == s1");
			}
			else {
				//Four-clique
				fourcliques++;

				//0101 -- cost
				//1010 -- cost
				real vals[16] = {0,0,0,0, 0,weight,0,0, 0,0,weight,0, 0,0,0,0};
				int vars[4] = {f1, f2, s1, s2};
				hocr.AddHigherTerm(4,vars,vals);
			}
		}
		statusOK();
		cout << threecliques << " three-cliques and " << fourcliques << " four-cliques."  << endl;
	}



	int nedges = 0;
	QPBO<real>::EdgeId e;
	for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e))
	{
		nedges++;
	}
	cerr << "Number of QPBO nodes : " << qpbo.GetNodeNum() << endl;
	cerr << "Number of QPBO edges : " << nedges << endl;


	statusTry("Merging parallel edges....");
	qpbo.MergeParallelEdges();
	nedges = 0;
	for (e=qpbo.GetNextEdgeId(-1); e>=0; e=qpbo.GetNextEdgeId(e)) {
		nedges++;
	}
	statusOK();
	cerr << "Number of QPBO edges : " << nedges << endl;

	//
	// Run optimizer
	//
	statusTry("Running QPBO...");
	qpbo.Solve();
	qpbo.ComputeWeakPersistencies();
	int unlabelled=0;
	int* labels = new int[mesh.ncells()];
	for (int i=0; i<mesh.ncells(); ++i) {
		labels[i] = qpbo.GetLabel(i);
		if (qpbo.GetLabel(i) < 0)
			unlabelled++;
	}
	statusOK();

	cerr << "Unlabelled regions    : " << unlabelled << " (" << 100*double(unlabelled)/double(mesh.ncells()) << "%)" << endl;

	if (unlabelled > 0) {
		int *mapping = new int[qpbo.GetNodeNum()];
		int *tmp_mapping = new int[qpbo.GetNodeNum()];
		for (int i = 0; i < qpbo.GetNodeNum(); i++) {
			mapping[i] = i * 2;
			tmp_mapping[i] = i * 2;
		}

		// Probe options
		QPBO<real>::ProbeOptions options;
		options.C = real(1e14);
		options.dilation = 1;
		options.weak_persistencies = 1;
		//options.iters = 1; //Only v.1.1


		statusTry("Running QPBO-P...");
		qpbo.Probe(mapping, options);
		qpbo.ComputeWeakPersistencies();
		/*for (int iter=1;iter<=5;++iter) {
		qpbo.Probe(tmp_mapping, options);
		qpbo.ComputeWeakPersistencies();
		qpbo.MergeMappings(nvars,mapping,tmp_mapping);
		}*/
		statusOK();


		// Read out entire labelling again (as weak persistencies may have changed)
		unlabelled = 0;
		for (int nodeCount = 0; nodeCount < mesh.ncells(); nodeCount++) {
			labels[nodeCount] = (int)qpbo.GetLabel(mapping[nodeCount]/2);
			if (labels[nodeCount] >= 0) {
				labels[nodeCount] = (labels[nodeCount] + mapping[nodeCount]) % 2;
			}
			else {
				unlabelled++;
			}
		}

		if (0) {
			statusTry("Running QPBO-I...");
			for (int iter=1;iter<=30;++iter) {
				qpbo.Improve();
			}
			statusOK();

			// Read out entire labelling again (as weak persistencies may have changed)
			for (int nodeCount = 0; nodeCount < mesh.ncells(); nodeCount++) {
				labels[nodeCount] = (int)qpbo.GetLabel(mapping[nodeCount]/2);
				if (labels[nodeCount] >= 0) {
					labels[nodeCount] = (labels[nodeCount] + mapping[nodeCount]) % 2;
				}
			}
		}

		delete[] mapping;
		delete[] tmp_mapping;
	}

	cerr << "Unlabelled regions    : " << unlabelled << " (" << 100*double(unlabelled)/double(mesh.ncells()) << "%)" << endl;

	statusTry("Drawing mesh...");
	mesh.draw("test.x3d",labels);
	statusOK();

	delete[] labels;
}
