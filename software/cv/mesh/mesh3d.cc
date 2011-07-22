/*** Written by Petter Strandmark as an employee of Lund University, Sweden, September 2010 ***/
//
// Uses C++0x features, compile with g++/gcc option -std=c++0x
// or with Microsoft Visual Studio 2010 (preferred)
//
//


#include "mesh3d.hh"

#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <map>
#include <stdexcept>
#include <limits>
using namespace std;

#define assert(expr) if(!(expr)){std::string err = "Assertion failed: "; err+=(#expr); throw std::runtime_error(err.c_str());}

//
// Mesh3DPoint
//
bool Mesh3DPoint::operator==(const Mesh3DPoint& cmp) const
{  
	return ( (std::abs(x_-cmp.x_) + std::abs(y_-cmp.y_) + std::abs(z_-cmp.z_)) <= 1e-6);
}

Mesh3DPoint::Mesh3DPoint() : x_(0), y_(0), z_(0) {}
Mesh3DPoint::Mesh3DPoint(double x, double y, double z) : x_(x), y_(y), z_(z) {}
void Mesh3DPoint::operator+=(const Mesh3DPoint& rhs) {x_+=rhs.x_; y_+=rhs.y_; z_+=rhs.z_; }
void Mesh3DPoint::operator/=(double rhs) {x_/=rhs; y_/=rhs; z_/=rhs; }
void Mesh3DPoint::operator-=(const Mesh3DPoint& rhs) {x_-=rhs.x_; y_-=rhs.y_; z_-=rhs.z_; }
void Mesh3DPoint::operator*=(double rhs) {x_*=rhs; y_*=rhs; z_*=rhs; }
void Mesh3DPoint::operator^=(double rhs) {x_=pow(x_,rhs); y_=pow(y_,rhs); z_=pow(z_,rhs); }
double Mesh3DPoint::dist(Mesh3DPoint to) const
{
	to -= *this;
	to ^= 2;
	return sqrt(to.x_ + to.y_ +  to.z_);
}

arma::vec3 operator-(const Mesh3DPoint& p1, const Mesh3DPoint& p2)
{
	arma::vec3 v;
	v[0] = p1.x_ - p2.x_;
	v[1] = p1.y_ - p2.y_;
	v[2] = p1.z_ - p2.z_;
	return v;
}


//
// Mesh3DEdge
//
Mesh3DEdge::Mesh3DEdge(uint from_idx, uint to_idx) : from_idx_(from_idx), to_idx_(to_idx) {}

//
// Mesh3DFace
//
bool Mesh3DFace::operator==(const std::vector<uint>& edge_idx) const
{
	std::set<uint> edges1(edge_idx_.begin(), edge_idx_.end()),
		edges2(edge_idx.begin(), edge_idx.end());
	return edges1 == edges2;
}



//
// Mesh3D
//

Mesh3D::Mesh3D() 
{
}

const Mesh3DCell& Mesh3D::cell(size_t idx) const {return cell_[idx];}
const Mesh3DFace& Mesh3D::face(size_t idx) const {return face_[idx];}
const Mesh3DEdge& Mesh3D::edge(size_t idx) const {return edge_[idx];}
const Mesh3DPoint& Mesh3D::point(size_t idx) const {return point_[idx];}

Mesh3DPoint Mesh3D::cell_center(size_t idx) const
{
	Mesh3DPoint p;
	const vector<uint>& face_idx = cell_[idx].face_idx_;
	for (int i=0;i<face_idx.size();++i) {
		p += face_center(face_idx[i]);
	}
	p /= double(face_idx.size());

	return p;
}

Mesh3DPoint Mesh3D::face_center(size_t idx) const
{
	Mesh3DPoint p;
	const vector<uint>& edge_idx = face_[idx].edge_idx_;
	for (int i=0;i<edge_idx.size();++i) {
		p += edge_center(edge_idx[i]);
	}
	p /= double(edge_idx.size());

	return p;
}

Mesh3DPoint Mesh3D::edge_center(size_t idx) const
{
	Mesh3DPoint p;
	p += point_[edge_[idx].from_idx_];
	p += point_[edge_[idx].to_idx_];
	p/=2;
	return p;
}

const std::vector<uint>& Mesh3D::adjacent_cells(size_t fi) const
{
	return face_adjacent_cells_[fi];
}
const std::vector<uint>& Mesh3D::adjacent_faces(size_t ei) const
{
	return edge_adjacent_faces_[ei];
}
const std::vector<uint>& Mesh3D::adjacent_edges(size_t pi) const
{
	return point_adjacent_edges_[pi];
}

void Mesh3D::add_point(const Mesh3DPoint& new_point) 
{
	point_.push_back(new_point);
	point_adjacent_edges_.push_back(std::vector<uint>());

	max_x_ = std::max(max_x_,new_point.x_);
	max_y_ = std::max(max_x_,new_point.y_); 
	max_z_ = std::max(max_z_,new_point.z_); 
	min_x_ = std::min(min_x_,new_point.x_);
	min_y_ = std::min(min_y_,new_point.y_);
	min_z_ = std::min(min_z_,new_point.z_); 
}
uint Mesh3D::find_point(Mesh3DPoint point) const
{
	for (uint i=0; i < point_.size(); i++) {
		if (point_[i] == point)
			return i;
	}

	return std::numeric_limits<uint>::max();
}
uint Mesh3D::find_or_add_point(Mesh3DPoint point) 
{
	uint ind = find_point(point);
	if (ind < std::numeric_limits<uint>::max()) {
		return ind;
	}
	else {
		add_point(point);
		return uint( point_.size()-1 );
	}
}
uint Mesh3D::find_or_add_point(double x, double y, double z)
{	
	return find_or_add_point(Mesh3DPoint(x,y,z));
}

uint Mesh3D::find_or_add_edge(uint to, uint from)
{
	std::vector<uint>::const_iterator it;

	//
	// Check whether this edge already exists
	//
	for (it = point_adjacent_edges_[to].begin(); it != point_adjacent_edges_[to].end(); it++) {
		uint i = *it;
		const Mesh3DEdge& edge = edge_[i];

		if (edge.from_idx_ == from && edge.to_idx_ == to) 
			return i;
		else if (edge.from_idx_ == to && edge.to_idx_ == from)
			return i;
	}

	//If not, add it
	add_edge(to, from);
	return uint( edge_.size() - 1 );
}

void Mesh3D::add_edge(uint from_idx, uint to_idx)
{
	assert(from_idx < point_.size());
	assert(to_idx < point_.size());

	edge_.push_back(Mesh3DEdge(from_idx,to_idx));
	edge_adjacent_faces_.push_back(std::vector<uint>());

	uint idx = uint( edge_.size()-1 );
	point_adjacent_edges_[from_idx].push_back(idx);
	point_adjacent_edges_[to_idx].push_back(idx);
}

void Mesh3D::add_face(const std::vector<uint>& edge_indices)
{
	face_.push_back(Mesh3DFace());
	face_.back().edge_idx_ = edge_indices;
	face_adjacent_cells_.push_back(std::vector<uint>());

	uint idx = uint( face_.size()-1 );
	std::set<uint> all_points;
	std::set<uint> point_indices;
	for (size_t i=0; i < edge_indices.size(); i++) {
		uint edge_idx = edge_indices[i];
		assert(edge_idx < edge_.size());
		edge_adjacent_faces_[edge_idx].push_back(idx);    

		uint next_edge_idx = edge_indices[(i+1) % edge_indices.size()]; 

		all_points.insert(edge_[edge_idx].to_idx_);
		all_points.insert(edge_[edge_idx].from_idx_);

		point_indices.clear();
		point_indices.insert(edge_[edge_idx].to_idx_);
		point_indices.insert(edge_[edge_idx].from_idx_);
		point_indices.insert(edge_[next_edge_idx].to_idx_);
		point_indices.insert(edge_[next_edge_idx].from_idx_);

		if (point_indices.size() != 3) {
			throw runtime_error("Not a proper face. Exiting...");
		}
	}

	//This code may be removed if LAPACK is not available since it
	//just checks that the face is actually a face
	arma::mat points(3,all_points.size());
	size_t c = 0;
	auto first_point = point_[*all_points.begin()];
	for (auto p=all_points.begin(); p != all_points.end(); ++p) {
		points.col(c++) = first_point - point_[*p];
	}
	int rank = arma::rank(points);
	if (rank != 2) {
		std::stringstream sout;
		sout << "Points do not constitute a plane face, rank(points)==" << rank;
		throw runtime_error(sout.str().c_str());
	}
}

uint Mesh3D::find_or_add_face(const std::vector<uint>& edge_indices)
{
	if (edge_indices.size() == 0) {
		return UINT_MAX;
	}

	std::vector<uint>::const_iterator it;
	uint edge = edge_indices[0];

	//
	// Check if this face already exists
	//
	for (it = edge_adjacent_faces_[edge].begin(); it != edge_adjacent_faces_[edge].end(); it++) {
		uint i = *it;
		const Mesh3DFace& face = face_[i];

		if (face == edge_indices) 
			return i;
	}

	//If not, add it
	add_face(edge_indices);
	return uint( face_.size()-1 );
}

void Mesh3D::add_cell(const std::vector<uint>& face_indices)
{
	cell_.push_back(Mesh3DCell());
	cell_.back().face_idx_ = face_indices;
	uint idx = uint( cell_.size()-1 );

	for (size_t i=0; i < face_indices.size(); i++) {
		face_adjacent_cells_[face_indices[i]].push_back(idx);
	}
}

void Mesh3D::generate_face_pair_list() 
{
	face_pairs_.clear();
	edge_adjacent_pairs_.clear();
	edge_adjacent_pairs_.reserve(edge_.size());

	for (uint i=0; i < edge_.size(); i++) {
		edge_adjacent_pairs_.push_back( vector<uint>() );

		uint nFaces= uint( edge_adjacent_faces_[i].size() );

		for (uint j=0; j < nFaces; j++) {
			for (uint j2 = j+1; j2 < nFaces; j2++) {
				Mesh3DFacePair pair;
				pair.first_face_idx_  = edge_adjacent_faces_[i][j];
				pair.second_face_idx_ = edge_adjacent_faces_[i][j2];
				pair.common_edge_idx_ = i;
				face_pairs_.push_back(pair);
				edge_adjacent_pairs_[i].push_back( uint(face_pairs_.size()-1) );
			}
		}
	}
}

size_t Mesh3D::npairs() const 
{
	if (face_pairs_.size() == 0) {
		throw runtime_error("No face pairs generated");
	}
	return face_pairs_.size();
}
const Mesh3DFacePair& Mesh3D::pair(size_t ip) const 
{
	if (face_pairs_.size() == 0) {
		throw runtime_error("No face pairs generated");
	}
	return face_pairs_.at(ip);
}
const std::vector<uint>& Mesh3D::adjacent_pairs(size_t edge_index) const
{
	if (face_pairs_.size() == 0) {
		throw runtime_error("No face pairs generated");
	}
	return edge_adjacent_pairs_.at(edge_index);
}


double triangle_area(const Mesh3DPoint& p1, const Mesh3DPoint& p2, const Mesh3DPoint& p3)
{
	arma::vec3 d1;
	arma::vec3 d2;

	d1[0] = p2.x_ - p1.x_;
	d1[1] = p2.y_ - p1.y_;
	d1[2] = p2.z_ - p1.z_;

	d2[0] = p3.x_ - p1.x_;
	d2[1] = p3.y_ - p1.y_;
	d2[2] = p3.z_ - p1.z_;

	return 0.5 * norm( cross(d1,d2), 2);
}
void Mesh3D::get_polygon_points(size_t face_idx, std::vector<uint>& point_indices) const {

	point_indices.clear();

	std::map<uint,uint> point_count;

	size_t nEdges = face_[face_idx].edge_idx_.size();

	for (uint j=0; j < nEdges; j++) {
		uint edge_idx = face_[face_idx].edge_idx_[j];
		uint next_edge_idx = face_[face_idx].edge_idx_[(j+1) % nEdges];

		point_count.clear();
		point_count[edge_[edge_idx].to_idx_]++;
		point_count[edge_[edge_idx].from_idx_]++;
		point_count[edge_[next_edge_idx].to_idx_]++;
		point_count[edge_[next_edge_idx].from_idx_]++;

		assert(point_count.size() == 3);

		for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
			if (it->second == 2) {
				point_indices.push_back(it->first);
				break;
			}
		}
	}
}

void Mesh3D::get_polygon_points(size_t face_idx, std::vector<Mesh3DPoint>& points) const 
{
	std::vector<uint> points_idx;
	get_polygon_points(face_idx,points_idx);
	points.clear();
	points.reserve(points_idx.size());
	for (uint i=0;i<points_idx.size(); ++i) {
		points.push_back( point_[points_idx[i]] );
	}
}


double Mesh3D::edge_length(size_t edge_idx) const
{
	const Mesh3DPoint& p1 = point_[ edge_[edge_idx].from_idx_];
	const Mesh3DPoint& p2 = point_[ edge_[edge_idx].to_idx_];
	return p1.dist(p2);
}

double Mesh3D::face_area(size_t face_idx) const 
{

	assert(face_[face_idx].edge_idx_.size() > 2);

	std::vector<uint> points;
	get_polygon_points(face_idx,points);

	double sum = 0.0;
	for (uint j=2; j < points.size(); j++) {
		sum += triangle_area(point_[points[0]],point_[points[j-1]],point_[points[j]]);
	}

	return sum;
}

double Mesh3D::face_perim(size_t face_idx) const
{
	auto& edges = face_[face_idx].edge_idx_;

	double perim = 0.0;
	for (auto e = edges.begin(); e != edges.end(); ++e) {
		perim += edge_length(*e);
	}
	return perim;
}

double Mesh3D::curvature_weight(size_t pair_index) const
{
	const Mesh3DFacePair& pair = face_pairs_.at(pair_index);

	size_t first  = pair.first_face_idx_;
	size_t second = pair.second_face_idx_;
	size_t e  = pair.common_edge_idx_;

	//Face normals
	auto fn = face_normal(first);
	auto sn = face_normal(second);

	if (match(first,second,e) < 0) {
		sn = -sn;
	}

	//Calculate angle between vectors
	double cosangle = (dot(fn,sn))/norm(fn,2)/norm(sn,2);
	double angle = 0;
	//This is to prevent NaN, because cosangle
	//might be equal to 1 + O(eps), which would
	//result in NaN
	if (cosangle < 1) { 
		angle = acos(cosangle);
	}

	//double weight = gamma*angle*angle;
	double el = edge_length(e);
	double twosinhalf = 2.0*sin(angle/2);

	//4.0 for squares? 3.0 for triangles
	double weight = (3.0*el*el)/(2.0*(face_area(first) + face_area(second))) * twosinhalf*twosinhalf;

	//Check for NaN
	if (weight != weight) {
#define print(expr) cout << #expr << " = " << (expr) << endl;
		print(weight);
		print(angle);
		print(cosangle);
		print(dot(fn,sn));
		throw runtime_error("NaN");
	}

	//return weight;


	//Areas of squares
	double fa = face_area(first);
	double sa = face_area(second);

	//Perimiters
	double fp = face_perim(first);
	double sp = face_perim(second);

	//Face centers 
	auto fc = face_center(first);
	auto sc = face_center(second);

	//Arc lengths
	auto& p1 = point_[edge_[e].from_idx_];
	auto& p2 = point_[edge_[e].to_idx_];
	auto fv1 = p1 - fc;
	auto fv2 = p2 - fc;
	double farc = acos(dot(fv1,fv2)/norm(fv1,2)/norm(fv2,2));
	auto sv1 = p1 - sc;
	auto sv2 = p2 - sc;
	double sarc = acos(dot(sv1,sv2)/norm(sv1,2)/norm(sv2,2));

	//Number of edges
	double fne = face_[first].edge_idx_.size();
	double sne = face_[second].edge_idx_.size();

	//Heights of triangles
	double fh = 2.0*fa/el;
	if (fne == 4)
		fh = fa/el;
	double sh = 2.0*sa/el;
	if (sne == 4) 
		sh = sa/el;

	//Radius of curvature
	double r = (fh+sh)/(2*angle);

	if (angle > 1e-5) {
		//cout << "el = " << el << " fh = " << fh << " sh = " << sh << " angle = " << angle << " r = " << r << endl;
	}

	double curv = 1/(r*r);

	weight = fa  * farc / (2*arma::math::pi())  * curv  
		+ sa  * sarc / (2*arma::math::pi())  * curv;


	return weight;
}

arma::vec3 Mesh3D::face_normal(size_t face_idx) const
{
	std::vector<uint> point_indices;
	get_polygon_points(face_idx, point_indices);
	if (point_indices.size() < 3) {
		throw runtime_error("Not enough points in face");
	}

	arma::vec3 u;
	u[0] = point_[point_indices[0]].x_ - point_[point_indices[1]].x_;
	u[1] = point_[point_indices[0]].y_ - point_[point_indices[1]].y_;
	u[2] = point_[point_indices[0]].z_ - point_[point_indices[1]].z_;
	arma::vec3 v;
	v[0] = point_[point_indices[0]].x_ - point_[point_indices[2]].x_;
	v[1] = point_[point_indices[0]].y_ - point_[point_indices[2]].y_;
	v[2] = point_[point_indices[0]].z_ - point_[point_indices[2]].z_;

	arma::vec3 n = cross(u,v);
	n = n / norm(n,2);
	return n;
}

//
// Matches the normals of two faces with
// respect to an edge
//
int Mesh3D::match(size_t first, size_t second, size_t e) const
{
	//Face normals
	auto fn = face_normal(first);
	auto sn = face_normal(second);
	//Face centers
	auto fc = face_center(first);
	auto sc = face_center(second);

	//Make them point toward the same side
	//Are vectors (close to) parallel?
	if ( norm( cross(fn,sn), 2) < 1e-4 ) {
		if ( dot(fn,sn) < 0 ) {
			//Normals point in opposite direction
			return -1;
		}
		else {
			//Normals point in same directions
			return 1;
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
		if ( dot(intf,fn)*dot(ints,sn) < 0) {
			//Normals point in opposite direction
			return -1;
		}
		else {
			//Normals point in same directions
			return 1;
		}
	}
}

void Mesh3D::X3D_start(ofstream& out, bool html_file) const
{
  if (html_file) {
    out << "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.0 Strict//EN\" \"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd\">" << endl;
    out << "<html xmlns=\"http://www.w3.org/1999/xhtml\">" << endl;
    out << "<head>" << endl;
    out << "  <meta http-equiv=\"Content-Type\"/>" << endl;
    out << "  <title>Viewer Test</title>" << endl;
    out << "  <link rel=\"stylesheet\" type=\"text/css\" href=\"x3dom/x3dom.css\"/>" << endl;
    out << "  <script type=\"text/javascript\" src=\"x3dom/x3dom.js\"> </script>" << endl;
    out << "</head>" << endl;
    out << "<body>" << endl;
    out << "<h1>Program output</h1>" << endl;
    out << "<X3D xmlns=\"http://www.web3d.org/specifications/x3d-namespace\" showStat=\"false\" showLog=\"false\" x=\"0px\" y=\"0px\" width=\"500px\" height=\"500px\">" << endl;
  }
  else {
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;
	  out << "<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\" \"http://www.web3d.org/specifications/x3d-3.0.dtd\">" << endl;
	  out << "<X3D>" << endl;
  }
	out << "<Scene>" << endl;
	out << "	<Background skyColor='1 1 1' />" << endl;
	out << "	<Viewpoint position=\'0 0 ";
	out << 1.4*max_z_;
	out << "\' />" << endl;
	out << "	<Group>" << endl;
	out << "		<DirectionalLight direction=\'0 -1 0\' intensity=\'0.4\'/> " << endl;
	out << "		<DirectionalLight direction=\'-1 0 0\' intensity=\'0.4\'/> " << endl;
	out << "		<DirectionalLight direction=\'0 0 -1\' intensity=\'0.4\'/> " << endl;
	out << "		<DirectionalLight direction=\'0 1 0\' intensity=\'0.4\'/> " << endl;
	out << "		<DirectionalLight direction=\'1 0 0\' intensity=\'0.4\'/> " << endl;
	out << "		<DirectionalLight direction=\'0 0 1\' intensity=\'0.4\'/> " << endl;
}

void Mesh3D::X3D_end(ofstream& out, bool html_file) const
{
	out << "	</Group>" << endl;
	out << "  </Scene>" << endl;
	out << "</X3D>" << endl;
  if (html_file) {
    out << "</body>" << endl;
    out << "</html>" << endl;
  }
}

void Mesh3D::draw(std::string filename, int* cell_labels, bool html_file) const
{
	ofstream out(filename.c_str());
	X3D_start(out, html_file);


	if (edge_.size() >= 3000 || edge_.size()==0) {
		//
		// Draw points
		//
		out << "	  <Shape>" << endl;
		out << "		<Appearance>" << endl;
		out << "		  <Material emissiveColor='1 0 0'/>" << endl;
		out << "		</Appearance>" << endl;
		out << "		<PointSet>" << endl;
		out << "		  <Coordinate point=\'";
		size_t ip;
		for (ip=0; ip<point_.size()-1; ++ip) {
			out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_ << ", ";
		}
		if (ip < point_.size()) {
			out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_;
		}
		out << "\'/>" << endl;
		out << "		</PointSet>" << endl;
		out << "	  </Shape>" << endl;

	} else {
		//
		// Draw edges
		//
		for (size_t ie=0; ie < edge_.size(); ++ie) {
			out << "	<Shape>" << endl;
			out << "		<Appearance>" << endl;
			out << "			<Material diffuseColor=\'0 0 0\' emissiveColor=\'1 0.5 0\'/>" << endl;
			out << "		</Appearance>" << endl;
			out << "		<IndexedLineSet coordIndex=\'0 1\'>" << endl;
			out << "			<Coordinate point='";
			Mesh3DPoint from = point_[edge_[ie].from_idx_];
			Mesh3DPoint to = point_[edge_[ie].to_idx_];
			out << from.x_ << " " << from.y_ << " " << from.z_ << ", ";
			out << to.x_ << " " << to.y_ << " " << to.z_ << ", ";
			out << "\'/>" << endl;
			out << "		</IndexedLineSet>" << endl;
			out << "	</Shape>" << endl;
		}
	}

	//
	// Draw faces
	//

	if (cell_labels) {
    out << "<Shape>" << endl;
	  out << "	<Appearance>" << endl;
	  out << "		<Material diffuseColor=\'0 0.5 1\'/>" << endl;
	  out << "	</Appearance>" << endl;

    out << "	<IndexedFaceSet coordIndex=\'";
		for (uint i=0; i < face_.size(); i++) {
			const vector<uint>& adjacent = adjacent_cells(i);
			if (adjacent.size()==2 && cell_labels[adjacent[0]] != cell_labels[adjacent[1]]) {

        //
        // Get all points
        //
        size_t nEdges = face_[i].edge_idx_.size();
        std::vector<uint> points;
	      for (uint k=0; k < nEdges; k++) {

		      std::map<uint,uint> point_count;

		      uint edge_idx = face_[i].edge_idx_[k];
		      uint next_edge_idx = face_[i].edge_idx_[(k+1) % nEdges];

		      point_count[edge_[edge_idx].to_idx_]++;
		      point_count[edge_[edge_idx].from_idx_]++;
		      point_count[edge_[next_edge_idx].to_idx_]++;
		      point_count[edge_[next_edge_idx].from_idx_]++;

		      if (point_count.size() != 3) {
			      throw runtime_error("point_count.size() != 3");
		      }
		      uint end_point = std::numeric_limits<uint>::max();
		      uint start_point = std::numeric_limits<uint>::max();

		      for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
			      if (it->second == 2) {
				      end_point = it->first;
				      break;
			      }
		      }

		      start_point = (edge_[edge_idx].to_idx_ == end_point) ? 
			      edge_[edge_idx].from_idx_ : edge_[edge_idx].to_idx_;

		      points.push_back(end_point);
	      }
        // We have all points

		    for (size_t i=0; i<points.size(); ++i) {
			    out << points[i] << " ";
		    }
        out << "-1, ";
		    for (int i=int(points.size()-1); i>=0; --i) {
			    out << points[i] << " ";
		    }
        out << "-1, ";
			}
		}
    out << "\'>" << endl;

    out << "    <Coordinate point=\'";
    size_t ip;
    for (ip=0; ip<point_.size()-1; ++ip) {
		  out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_ << ", ";
	  }
	  if (ip < point_.size()) {
		  out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_;
	  }
    out << "\'/>" << endl;

    out << "  </IndexedFaceSet>" << endl;
	  out << "</Shape>"<<endl;
	}

	X3D_end(out, html_file);
}


void Mesh3D::draw_faces(std::string filename, int* face_labels, bool html_file) const
{
	ofstream out(filename.c_str());
	X3D_start(out, html_file);


	if (edge_.size() >= 3000) {
		//
		// Draw points
		//
		out << "	  <Shape>" << endl;
		out << "		<Appearance>" << endl;
		out << "		  <Material emissiveColor='1 0 0'/>" << endl;
		out << "		</Appearance>" << endl;
		out << "		<PointSet>" << endl;
		out << "		  <Coordinate point=\'";
		size_t ip;
		for (ip=0; ip<point_.size()-1; ++ip) {
			out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_ << ", ";
		}
		if (ip < point_.size()) {
			out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_;
		}
		out << "\'/>" << endl;
		out << "		</PointSet>" << endl;
		out << "	  </Shape>" << endl;

	} else {
		//
		// Draw edges
		//
		for (size_t ie=0; ie < edge_.size(); ++ie) {
			out << "	<Shape>" << endl;
			out << "		<Appearance>" << endl;
			out << "			<Material diffuseColor=\'0 0 0\' emissiveColor=\'1 0.5 0\'/>" << endl;
			out << "		</Appearance>" << endl;
			out << "		<IndexedLineSet coordIndex=\'0 1\'>" << endl;
			out << "			<Coordinate point='";
			Mesh3DPoint from = point_[edge_[ie].from_idx_];
			Mesh3DPoint to = point_[edge_[ie].to_idx_];
			out << from.x_ << " " << from.y_ << " " << from.z_ << ", ";
			out << to.x_ << " " << to.y_ << " " << to.z_ << ", ";
			out << "\'/>" << endl;
			out << "		</IndexedLineSet>" << endl;
			out << "	</Shape>" << endl;
		}
	}

	//
	// Draw faces
	//
	if (face_labels) {

    out << "<Shape>" << endl;
	  out << "	<Appearance>" << endl;
	  out << "		<Material diffuseColor=\'0 0.5 1\'/>" << endl;
	  out << "	</Appearance>" << endl;

    out << "	<IndexedFaceSet coordIndex=\'";
		for (uint i=0; i < face_.size(); i++) {
			const vector<uint>& adjacent = adjacent_cells(i);
			if (face_labels[2*i] || face_labels[2*i+1]) {
        auto& f = face_[i];

        //
        // Get all points
        //
        size_t nEdges = face_[i].edge_idx_.size();
        std::vector<uint> points;
	      for (uint k=0; k < nEdges; k++) {

		      std::map<uint,uint> point_count;

		      uint edge_idx = face_[i].edge_idx_[k];
		      uint next_edge_idx = face_[i].edge_idx_[(k+1) % nEdges];

		      point_count[edge_[edge_idx].to_idx_]++;
		      point_count[edge_[edge_idx].from_idx_]++;
		      point_count[edge_[next_edge_idx].to_idx_]++;
		      point_count[edge_[next_edge_idx].from_idx_]++;

		      if (point_count.size() != 3) {
			      throw runtime_error("point_count.size() != 3");
		      }
		      uint end_point = std::numeric_limits<uint>::max();
		      uint start_point = std::numeric_limits<uint>::max();

		      for (std::map<uint,uint>::iterator it = point_count.begin(); it != point_count.end(); it++) {
			      if (it->second == 2) {
				      end_point = it->first;
				      break;
			      }
		      }

		      start_point = (edge_[edge_idx].to_idx_ == end_point) ? 
			      edge_[edge_idx].from_idx_ : edge_[edge_idx].to_idx_;

		      points.push_back(end_point);
	      }
        // We have all points



        if (face_labels[2*i]) {
		      for (size_t i=0; i<points.size(); ++i) {
			      out << points[i] << " ";
		      }
          out << "-1, ";
        }
        else {
		      for (int i=int(points.size()-1); i>=0; --i) {
			      out << points[i] << " ";
		      }
          out << "-1, ";
        }

			}
		}
    out << "\'>" << endl;

    out << "    <Coordinate point=\'";
    size_t ip;
    for (ip=0; ip<point_.size()-1; ++ip) {
		  out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_ << ", ";
	  }
	  if (ip < point_.size()) {
		  out << point_[ip].x_ << " " << point_[ip].y_ << " " << point_[ip].z_;
	  }
    out << "\'/>" << endl;

    out << "  </IndexedFaceSet>" << endl;
	  out << "</Shape>"<<endl;
	}

	X3D_end(out, html_file);
}




void Mesh3D::subdivide_cells()
{
	throw std::runtime_error("subdivide_cells not implemented"); 
}


void generate_3D_cubic_mesh(size_t xdim, size_t ydim, size_t zdim, Mesh3D& mesh)
{
	for (int ix=0;ix < xdim; ++ix) {
		for (int iy=0;iy < ydim; ++iy) {
			for (int iz=0;iz < zdim; ++iz) {
				double x=ix,y=iy,z=iz;
				uint p1,p2,p3,p4,p5,p6,p7,p8;
				p1 = mesh.find_or_add_point(x  ,y  ,z);
				p2 = mesh.find_or_add_point(x+1,y  ,z);
				p3 = mesh.find_or_add_point(x  ,y+1,z);
				p4 = mesh.find_or_add_point(x+1,y+1,z);
				p5 = mesh.find_or_add_point(x  ,y  ,z+1);
				p6 = mesh.find_or_add_point(x+1,y  ,z+1);
				p7 = mesh.find_or_add_point(x  ,y+1,z+1);
				p8 = mesh.find_or_add_point(x+1,y+1,z+1);

				//Cubic edges
				uint e1 = mesh.find_or_add_edge(p1,p2);
				uint e2 = mesh.find_or_add_edge(p1,p3);
				uint e3 = mesh.find_or_add_edge(p2,p4);
				uint e4 = mesh.find_or_add_edge(p3,p4);

				uint e5 = mesh.find_or_add_edge(p1,p5);
				uint e6 = mesh.find_or_add_edge(p2,p6);
				uint e7 = mesh.find_or_add_edge(p4,p8);
				uint e8 = mesh.find_or_add_edge(p3,p7);

				uint e9  = mesh.find_or_add_edge(p5,p6);
				uint e10 = mesh.find_or_add_edge(p5,p7);
				uint e11 = mesh.find_or_add_edge(p7,p8);
				uint e12 = mesh.find_or_add_edge(p6,p8);

				std::vector<uint> edges,faces;

				//Cube tetraeder
				edges.push_back(e1);
				edges.push_back(e2);
				edges.push_back(e4);
				edges.push_back(e3);
				uint f1 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e1);
				edges.push_back(e6);
				edges.push_back(e9);
				edges.push_back(e5);
				uint f2 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e7);
				edges.push_back(e12);
				edges.push_back(e6);
				edges.push_back(e3);
				uint f3 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e4);
				edges.push_back(e8);
				edges.push_back(e11);
				edges.push_back(e7);
				uint f4 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e5);
				edges.push_back(e10);
				edges.push_back(e8);
				edges.push_back(e2);
				uint f5 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e12);
				edges.push_back(e11);
				edges.push_back(e10);
				edges.push_back(e9);
				uint f6 = mesh.find_or_add_face(edges);
				faces.clear();
				faces.push_back(f1);
				faces.push_back(f2);
				faces.push_back(f3);
				faces.push_back(f4);
				faces.push_back(f5);
				faces.push_back(f6);
				mesh.add_cell(faces);
			}
		}
	}
}

void generate_3D_mesh(size_t xdim, size_t ydim, size_t zdim, Mesh3D& mesh)
{
	for (int ix=0;ix < xdim; ++ix) {
		for (int iy=0;iy < ydim; ++iy) {
			for (int iz=0;iz < zdim; ++iz) {
				//
				// To make the edges and faces fit, every other
				// voxel needs to be mirrored
				//
				int xflip=1,yflip=1,zflip=1;
				if (ix%2==1) xflip=-1;
				if (iy%2==1) yflip=-1;
				if (iz%2==1) zflip=-1;
				bool flip = xflip*yflip*zflip==-1;

				double x=ix,y=iy,z=iz;

				uint p1,p2,p3,p4,p5,p6,p7,p8;
				if (!flip) {
					p1 = mesh.find_or_add_point(x  ,y  ,z);
					p2 = mesh.find_or_add_point(x+1,y  ,z);
					p3 = mesh.find_or_add_point(x  ,y+1,z);
					p4 = mesh.find_or_add_point(x+1,y+1,z);
					p5 = mesh.find_or_add_point(x  ,y  ,z+1);
					p6 = mesh.find_or_add_point(x+1,y  ,z+1);
					p7 = mesh.find_or_add_point(x  ,y+1,z+1);
					p8 = mesh.find_or_add_point(x+1,y+1,z+1);
				}
				else {
					//Mirrored version
					p2 = mesh.find_or_add_point(x  ,y  ,z);
					p1 = mesh.find_or_add_point(x+1,y  ,z);
					p4 = mesh.find_or_add_point(x  ,y+1,z);
					p3 = mesh.find_or_add_point(x+1,y+1,z);
					p6 = mesh.find_or_add_point(x  ,y  ,z+1);
					p5 = mesh.find_or_add_point(x+1,y  ,z+1);
					p8 = mesh.find_or_add_point(x  ,y+1,z+1);
					p7 = mesh.find_or_add_point(x+1,y+1,z+1);
				}

				//Cubic edges
				uint e1 = mesh.find_or_add_edge(p1,p2);
				uint e2 = mesh.find_or_add_edge(p1,p3);
				uint e3 = mesh.find_or_add_edge(p2,p4);
				uint e4 = mesh.find_or_add_edge(p3,p4);

				uint e5 = mesh.find_or_add_edge(p1,p5);
				uint e6 = mesh.find_or_add_edge(p2,p6);
				uint e7 = mesh.find_or_add_edge(p4,p8);
				uint e8 = mesh.find_or_add_edge(p3,p7);

				uint e9  = mesh.find_or_add_edge(p5,p6);
				uint e10 = mesh.find_or_add_edge(p5,p7);
				uint e11 = mesh.find_or_add_edge(p7,p8);
				uint e12 = mesh.find_or_add_edge(p6,p8);

				//Tetrahedonal edges
				uint e13 = mesh.find_or_add_edge(p2,p5);
				uint e14 = mesh.find_or_add_edge(p5,p3);
				uint e15 = mesh.find_or_add_edge(p3,p2);

				uint e16 = mesh.find_or_add_edge(p5,p8);
				uint e17 = mesh.find_or_add_edge(p8,p3);
				uint e18 = mesh.find_or_add_edge(p3,p5); //same

				uint e19 = mesh.find_or_add_edge(p2,p5); //same
				uint e20 = mesh.find_or_add_edge(p5,p8); //same
				uint e21 = mesh.find_or_add_edge(p8,p2);

				uint e22 = mesh.find_or_add_edge(p2,p3); //same
				uint e23 = mesh.find_or_add_edge(p3,p8); //same
				uint e24 = mesh.find_or_add_edge(p8,p2); //same

				std::vector<uint> edges,faces;

				//Middle tetraeder
				edges.push_back(e13);
				edges.push_back(e14);
				edges.push_back(e15);
				uint f1 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e20);
				edges.push_back(e13);
				edges.push_back(e24);
				uint f2 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e14);
				edges.push_back(e20);
				edges.push_back(e23);
				uint f3 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e23);
				edges.push_back(e21);
				edges.push_back(e15);
				uint f4 = mesh.find_or_add_face(edges);
				faces.clear();
				faces.push_back(f1);
				faces.push_back(f2);
				faces.push_back(f3);
				faces.push_back(f4);
				mesh.add_cell(faces);

				//p1 tetraeder
				edges.clear();
				edges.push_back(e1);
				edges.push_back(e13);
				edges.push_back(e5);
				uint f5 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e5);
				edges.push_back(e14);
				edges.push_back(e2);
				uint f6 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e2);
				edges.push_back(e15);
				edges.push_back(e1);
				uint f7 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e14);
				edges.push_back(e13);
				edges.push_back(e15);
				uint f8 = mesh.find_or_add_face(edges);
				faces.clear();
				faces.push_back(f5);
				faces.push_back(f6);
				faces.push_back(f7);
				faces.push_back(f8);
				mesh.add_cell(faces);

				//p4 tetraeder
				edges.clear();
				edges.push_back(e3);
				edges.push_back(e15);
				edges.push_back(e4);
				uint f9 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e7);
				edges.push_back(e24);
				edges.push_back(e3);
				uint f10 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e24);
				edges.push_back(e23);
				edges.push_back(e15);
				uint f11 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e4);
				edges.push_back(e23);
				edges.push_back(e7);
				uint f12 = mesh.find_or_add_face(edges);
				faces.clear();
				faces.push_back(f9);
				faces.push_back(f10);
				faces.push_back(f11);
				faces.push_back(f12);
				mesh.add_cell(faces);

				//p6 tetraeder
				edges.clear();
				edges.push_back(e24);
				edges.push_back(e12);
				edges.push_back(e6);
				uint f13 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e12);
				edges.push_back(e16);
				edges.push_back(e9);
				uint f14 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e6);
				edges.push_back(e9);
				edges.push_back(e13);
				uint f15 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e16);
				edges.push_back(e24);
				edges.push_back(e13);
				uint f16 = mesh.find_or_add_face(edges);
				faces.clear();
				faces.push_back(f13);
				faces.push_back(f14);
				faces.push_back(f15);
				faces.push_back(f16);
				mesh.add_cell(faces);


				//p7 tetraeder
				edges.clear();
				edges.push_back(e11);
				edges.push_back(e10);
				edges.push_back(e16);
				uint f17 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e23);
				edges.push_back(e8);
				edges.push_back(e11);
				uint f18 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e14);
				edges.push_back(e10);
				edges.push_back(e8);
				uint f19 = mesh.find_or_add_face(edges);
				edges.clear();
				edges.push_back(e16);
				edges.push_back(e14);
				edges.push_back(e23);
				uint f20 = mesh.find_or_add_face(edges);
				faces.clear();
				faces.push_back(f17);
				faces.push_back(f18);
				faces.push_back(f19);
				faces.push_back(f20);
				mesh.add_cell(faces);

			}
		}
	}

}

