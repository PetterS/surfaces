/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/
//
// Uses C++0x features, compile with g++/gcc option -std=c++0x
// or with Microsoft Visual Studio 2010 (preferred)
//
//


//////////////////////////////////
// Functions to add constraints //
//////////////////////////////////
// 
// Variables are assumed to be ordered as follows
//
//   Order | Offset           | Number of vars       | Description
//   ------|------------------|----------------------|--------------
//     1   | 0                | 2*mesh.nfaces()      | Variables for every face of the mesh.
//     2   | 2*mesh.nfaces()  | 2*face_pairs.size()  | Variables for every pair of faces. (if present)
//     3   | cells_var_offset | mesh.ncells()        | Variables for every volume element
//
// This means that cells_var_offset is equal to either 2*mesh.nfaces() (no face pairs)
// or equal to 2*mesh.nfaces() + 2*face_pairs.size() (with variables for face pairs)
//
//
// The functions use their argument add_element to add entries to the linear programming
// constraint matrix. They will start adding rows at offset and update the offset afterwards.
// Thus, the variable offset will always contain the row at which the next constraint should
// be added.
//


#ifndef LP_CONSTRAINTS_HH
#define LP_CONSTRAINTS_HH

#include <vector>
#include <stdexcept>
using std::vector;

#include "mesh3d.hh"

#include "Petter-Color.h"

//
// Edge consistency constraints
//
template<typename Adder>
void edge_consistency_constraints(const Adder& add_element, const Mesh3D& mesh, size_t& offset)
{
	Petter::statusTry("Edge consistency constraints...");
	for (uint e=0; e < mesh.nedges(); ++e) {
		//First point of the edge
		uint p1 = mesh.edge(e).from_idx_;
		uint p2 = mesh.edge(e).to_idx_;

		//Faces adjacent to this edge
		const vector<uint>& faces = mesh.adjacent_faces(e);
		for (size_t fi=0; fi<faces.size(); ++fi) {
			uint f = faces[fi];   

			//
			//How does this edge agree with the normal?
			//
			int incidence = 0;
			auto normal = mesh.face_normal(f);
			Mesh3DPoint center = mesh.face_center(f);
			auto v1 = mesh.point(p1) - center;
			auto v2 = mesh.point(p2) - center;
			auto c = cross(v1,v2);
			if ( dot(c,normal) < 0) {
				incidence = -1;
			}
			else {
				incidence = 1;
			}

			//
			// Add to constraint matrix 
			//
			if (incidence == 1) {
				add_element(e, 2*f,    1);
				add_element(e, 2*f+1, -1);
			}
			else {
				add_element(e, 2*f,   -1);
				add_element(e, 2*f+1,  1);
			}

		}
	}
	offset += mesh.nedges();
	Petter::statusOK();
}

//
// Opposing faces constraints
//
template<typename Adder, typename RHSChanger>
void opposing_faces_constraints(const Adder& add_element, const RHSChanger& change_rhs, const Mesh3D& mesh, size_t& offset)
{
	Petter::statusTry("Opposing faces constraints...");
	for (uint f=0;f<mesh.nfaces();++f) {
		add_element(offset+f, 2*f, 1);
		add_element(offset+f, 2*f+1, 1);
		change_rhs(offset+f, 0, 1);
	}
	offset += mesh.nfaces();
	Petter::statusOK();
}

//
// Surface continuation constraints
//
template<typename Adder>
void surface_continuation_constraints(const Adder& add_element, const Mesh3D& mesh, size_t& offset)
{
	Petter::statusTry("Surface continuation constraints...");

	size_t con = 0;
	for (uint f=0; f < mesh.nfaces(); ++f) {
		//Edges of this face
		const vector<uint>& edges = mesh.face(f).edge_idx_;

		for (size_t j=0; j < edges.size(); ++j) {
			//Index of this edge
			uint e = edges[j]; 

			add_element(offset+con  , 2*f  ,  1);
			add_element(offset+con+1, 2*f+1,  1);

			//The pairs adjacent to this 
			//vector<uint> pairs = edge_adjacent_pairs[e];
			const vector<uint>& pairs = mesh.adjacent_pairs(e);

			for (size_t jj=0; jj<pairs.size(); ++jj) {
				//Index of this pair
				uint p = pairs[jj];
				const Mesh3DFacePair& pair = mesh.pair(p);

				//Sanity check
				if (e != pair.common_edge_idx_) {
					throw std::runtime_error("e != pair.common_edge_idx_");
				}

				if (f == pair.first_face_idx_) {
					//Normal of first face determines orientation of face pair
					add_element(offset+con  , 2*mesh.nfaces() + 2*p  ,  -1);
					add_element(offset+con+1, 2*mesh.nfaces() + 2*p+1,  -1);
				}
				else if (f == pair.second_face_idx_) {
					//Do the normals of the face pair match?
					int match = mesh.match( pair.first_face_idx_, pair.second_face_idx_, e);
					if (match == 1) {
						add_element(offset+con  , 2*mesh.nfaces() + 2*p  ,  -1);
						add_element(offset+con+1, 2*mesh.nfaces() + 2*p+1,  -1);
					}
					else {
						add_element(offset+con  , 2*mesh.nfaces() + 2*p+1,  -1);
						add_element(offset+con+1, 2*mesh.nfaces() + 2*p  ,  -1);
					}
				}
				else {
					//This edge pair has nothing to do with f
				}
			}

			con+=2;
		}
	}
	offset += con;
	Petter::statusOK();

	std::cerr << con << " continuation constraints. " << std::endl;
}


template<typename Adder>
void volume_face_constraints(const Adder& add_element, const Mesh3D& mesh, size_t& offset, size_t cells_var_offset)
{
	Petter::statusTry("Volume/face consistency constraints...");
	for (uint f=0; f < mesh.nfaces(); ++f) {
		//Cells of this face
		const vector<uint>& cells = mesh.adjacent_cells(f);
		if (cells.size() != 2) {
			//Face at the boundary
			continue;
		}

		add_element(offset+f  , 2*f  ,  1);
		add_element(offset+f  , 2*f+1, -1);

		uint c1 = cells[0];
		uint c2 = cells[1];
		auto normal = mesh.face_normal(f);
		auto center1 = mesh.cell_center(c1);
		auto center2 = mesh.cell_center(c2);
		auto face_center = mesh.face_center(f);
		auto vc1 = center1 - face_center; //From face to the center of v1
		auto vc2 = center2 - face_center;
		if ( dot(vc1,normal) > 0) {
			if ( dot(vc2,normal) >= 0) {
				throw std::runtime_error("Cells on same side");
			}

			//The face normal points towards c1
			add_element(offset+f   , cells_var_offset + c1,   1);
			add_element(offset+f   , cells_var_offset + c2,  -1);

		}
		else {
			if ( dot(vc2,normal) <= 0) {
				throw std::runtime_error("Cells on same side");
			}

			//The face normal points away from c1
			add_element(offset+f   , cells_var_offset + c1,  -1);
			add_element(offset+f   , cells_var_offset + c2,   1);
		}
	}
	Petter::statusOK();

	offset += mesh.nfaces();
}

#endif
