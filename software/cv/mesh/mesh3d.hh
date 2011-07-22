/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/

#ifndef MESH_3D_HH
#define MESH_3D_HH

//Type used a lot
typedef unsigned int uint;

#include <iostream>
#include <vector>
#include <cstddef>
using std::size_t;

#include <armadillo>

//
// Holds a point in the mesh
//
struct Mesh3DPoint 
{
	Mesh3DPoint(double x, double y, double z);
	Mesh3DPoint();

	bool operator==(const Mesh3DPoint& cmp) const;
	void operator+=(const Mesh3DPoint& rhs);
	void operator/=(double rhs);
	void operator-=(const Mesh3DPoint& rhs);
	void operator*=(double rhs);
	void operator^=(double rhs);
	double dist(Mesh3DPoint) const;

	double x_;
	double y_;
	double z_;
};

//Subtracting two points gives a vector from p2 to p1
arma::vec3 operator-(const Mesh3DPoint& p1, const Mesh3DPoint& p2);

//
// Holds an edge in the mesh, only makes sense together with its mesh,
// since it holds indices
//
struct Mesh3DEdge {

	Mesh3DEdge(uint from_idx, uint to_idx); 

	uint from_idx_;
	uint to_idx_;
};

//
// Holds a face in the mesh, only makes sense together with its mesh,
// since it holds indices
//
struct Mesh3DFace {
	std::vector<uint> edge_idx_;

	bool operator==(const std::vector<uint>& edge_idx) const;
};

//
// Holds a pair of faces in the mesh, only makes sense together with its mesh,
// since it holds indices
//
struct Mesh3DFacePair {
	uint first_face_idx_;
	uint second_face_idx_;
	uint common_edge_idx_;
};

//
// Holds a cell (volume element) in the mesh, only makes sense together with its mesh,
// since it holds indices
//
struct Mesh3DCell {
	std::vector<uint> face_idx_;
};

//
// Main class. Holds a complete 3D mesh.
//
class Mesh3D {
public:

	Mesh3D();

	//
	// Points
	//
	const Mesh3DPoint& point(size_t idx) const;
	size_t npoint() const {return point_.size();}
	void add_point(const Mesh3DPoint& new_point);
	uint find_point(Mesh3DPoint point) const;
	uint find_or_add_point(Mesh3DPoint point);
	uint find_or_add_point(double x, double y, double z);

	//
	// Edges
	//
	const Mesh3DEdge& edge(size_t idx) const;
	size_t nedges() const {return edge_.size();} 
	void add_edge(uint from_idx, uint to_idx);
	uint find_or_add_edge(uint to, uint from);

	Mesh3DPoint edge_center(size_t idx) const;
	double edge_length(size_t edge_idx) const;

	//
	// Faces
	//
	const Mesh3DFace& face(size_t idx) const;
	size_t nfaces() const {return face_.size();} 
	void add_face(const std::vector<uint>& edge_indices);
	uint find_or_add_face(const std::vector<uint>& edge_indices);

	Mesh3DPoint face_center(size_t idx) const;
	double face_area(size_t face_idx) const;
	double face_perim(size_t face_idx) const;
	arma::vec3 face_normal(size_t face_idx) const;

	//
	// Pairs of faces
	//
	const Mesh3DFacePair& pair(size_t ip) const;
	size_t npairs() const;
	//Must be called before any face pair function are called.
	//If the mesh is modified, it probably needs to be called again.
	void generate_face_pair_list(); 

	//Calculates the curvature weight for a face pair
	double curvature_weight(size_t pair_index) const;

	//
	// Cells (volume elements)
	//
	const Mesh3DCell& cell(size_t idx) const;
	size_t ncells() const {return cell_.size();}
	void add_cell(const std::vector<uint>& face_indices);

	Mesh3DPoint cell_center(size_t idx) const;

	//
	// Adjacency functions
	//
	//Returns all cells touching the specified face
	// (usually 2, but equal to 1 at the borders)
	const std::vector<uint>& adjacent_cells(size_t fi) const;
	//Returns all faces touching the specified edge
	const std::vector<uint>& adjacent_faces(size_t ei) const;
	//Returns all edges touching the specified point
	const std::vector<uint>& adjacent_edges(size_t pi) const;

	//Returns the pairs containing the specified edge
	const std::vector<uint>& adjacent_pairs(size_t edge_index) const;

	//
	// Matches the normals of two faces with
	// respect to a common edge
	//
	//  returns +1 or -1
	//
	int match(size_t f1, size_t f2, size_t e) const;

	//
	// Drawing functions
	//
	void draw(std::string filename, int* cell_labels=0, bool html_file=false) const;
	void draw_faces(std::string filename, int* face_labels=0, bool html_file=false) const;

	//
	// Refinement functions
	//
	void subdivide_cells();

protected:

	void X3D_start(std::ofstream& out, bool html_file) const;
	void X3D_end(std::ofstream& out, bool html_file) const;

	void get_polygon_points(size_t face_idx, std::vector<uint>& point_indices) const;
	void get_polygon_points(size_t face_idx, std::vector<Mesh3DPoint>& points) const;

	std::vector<Mesh3DPoint> point_;
	std::vector<Mesh3DEdge> edge_;
	std::vector<Mesh3DFace> face_;
	std::vector<Mesh3DCell> cell_;

	std::vector<std::vector<uint> > point_adjacent_edges_;
	std::vector<std::vector<uint> > edge_adjacent_faces_;
	std::vector<std::vector<uint> > face_adjacent_cells_;

	std::vector<Mesh3DFacePair>  face_pairs_;
	std::vector<std::vector<uint> > edge_adjacent_pairs_;


	double min_x_;
	double min_y_;
	double min_z_;
	double max_x_;
	double max_y_;
	double max_z_;
};



//
// Functions to generate 3D meshes
//
void generate_3D_mesh(size_t xdim, size_t ydim, size_t zdim, Mesh3D& mesh);
void generate_3D_cubic_mesh(size_t xdim, size_t ydim, size_t zdim, Mesh3D& mesh);

#endif