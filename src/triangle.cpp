#include <triangle.hpp>

using namespace rt;

triangle_t::triangle_t(material_t* _mat, Eigen::Vector3d _verts[3])
{ 
	for(int i=0;i<3;i++)
	{
		verts[i] = _verts[i];
	}
	mat = _mat;
}

triangle_t::~triangle_t() { }

bool triangle_t::intersect(hit_t& result, const ray_t& _ray) const
{

	Eigen::Vector3d O = _ray.origin;
	Eigen::Vector3d D = _ray.direction;


	  Eigen::Vector3d e1, e2;  //Edge1, Edge2
	  Eigen::Vector3d P, Q, T;
	  double det, inv_det, u, v;
	  double t;

	  //Find vectors for two edges sharing V0
	  e1 = verts[1] - verts[0];
	  e2 = verts[2] - verts[0];
	  //Begin calculating determinant - also used to calculate u parameter
	  P = D.cross(e2);
	  //if determinant is near zero, ray lies in plane of triangle or ray is parallel to plane of triangle
	  det = e1.dot(P);
	  //NOT CULLING
	  if(det > -EPSILON && det < EPSILON) return false;
	  inv_det = 1.0 / det;

	  //calculate distance from V0 to ray origin
	  T = O - verts[0];

	  //Calculate u parameter and test bound
	  u = T.dot(P) * inv_det;
	  //The intersection lies outside of the triangle
	  if(u < 0.0 || u > 1.0) return false;

	  //Prepare to test v parameter
	  Q = T.cross(e1);

	  //Calculate V parameter and test bound
	  v = D.dot(Q) * inv_det;
	  //The intersection lies outside of the triangle
	  if(v < 0.0 || u + v  > 1.0) return false;

	  t = e2.dot(Q) * inv_det;

	  if(t > EPSILON) { //ray intersection
		result = hit_t(this,t);
	  }
	  else
	  {
	  	return false;
	  }
	  
	  return  t >= _ray.mint && t <= _ray.maxt;
}

Eigen::Vector3d triangle_t::get_normal(Eigen::Vector3d& _p) const
{
	return ((verts[1] - verts[0]).cross(verts[2] - verts[0])).normalized();
}

material_t* triangle_t::get_material(void) const
{
	return mat;
}

void triangle_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Object Properties: -------------------------------"<<std::endl;
	stream<<"Type: triangle"<<std::endl;
	stream<<"Vertex 0: "<<verts[0].format(CommaInitFmt)<<std::endl;
	stream<<"Vertex 1: "<<verts[1].format(CommaInitFmt)<<std::endl;
	stream<<"Vertex 2: "<<verts[2].format(CommaInitFmt)<<std::endl<<std::endl;
}