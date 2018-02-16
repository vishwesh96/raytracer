#include <torus.hpp>
#include <poly34.h>
#include <cfloat>

using namespace rt;
ray_t transform_to_local(ray_t world_ray , Eigen::Vector3d center, Eigen::Vector3d axis);

torus_t::torus_t(material_t* _mat):center(0.0,0.0,0.0),R(1.0), r(0.0), axis(0.0,1.0,0.0), mat(_mat) { }
torus_t::torus_t(material_t* _mat, Eigen::Vector3d _c, double _R, double _r, Eigen::Vector3d _axis)
{
	center = _c;
	R = _R;
	r = _r;
	axis = _axis.normalized();
	mat = _mat;
}

torus_t::~torus_t() { }

bool torus_t::intersect(hit_t& result, const ray_t& _ray) const
{

	ray_t local_ray = transform_to_local(_ray,center,axis);

	Eigen::Vector3d e = local_ray.origin;
	Eigen::Vector3d d = local_ray.direction;
	
	double iR = R*R - r*r;
	double oR = 4*R*R;
  	double l = e.squaredNorm() + iR;
  	double e_dot_d = e.dot(d);
	double c4 = pow(d.squaredNorm(),2);
	double c3 = 4*d.squaredNorm()*e_dot_d;
	double c2 = 2*(2*pow((e_dot_d),2) + d.squaredNorm()*l) - oR*(pow(d[0],2)+pow(d[1],2));
	double c1 = 4*e_dot_d*l - 2*oR*(e[0]*d[0]+e[1]*d[1]);
	double c0 = pow(l,2) - oR*(pow(e[0],2)+pow(e[1],2));

	double roots[4];
	int num_roots = SolveP4(roots,c3/c4,c2/c4,c1/c4,c0/c4);

	if(num_roots == 0)
	{
		return false;		
	} 
	bool intersect = false;
	double min_t = DBL_MAX;
 	for(int i=0;i<num_roots;i++)
	{
		double t = roots[i];
		if(t > EPSILON && t < min_t)
		{
			min_t = t;
			intersect = true;
		}
	}
	if(intersect)
	{
		result = hit_t(this,min_t);
	}
	else
	{
		return false;
	}
	return  min_t >= _ray.mint && min_t <= _ray.maxt;
}

Eigen::Vector3d torus_t::get_normal(Eigen::Vector3d& _p) const
{

	Eigen::Vector3d p = _p - center;
	double cos_theta = axis.dot(p)/(p.norm());
	double sin_theta = sqrt(1-cos_theta*cos_theta);

	Eigen::Vector3d normal = (p-R*(p.normalized() - cos_theta*axis)/sin_theta)/r;
	normal.normalize();
	return normal;
}

material_t* torus_t::get_material(void) const
{
	return mat;
}

void torus_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Object Properties: -------------------------------"<<std::endl;
	stream<<"Type: torus"<<std::endl;
	stream<<"Center: "<<center.format(CommaInitFmt)<<std::endl;
	stream<<"Outer radius: "<<R<<std::endl;
	stream<<"Inner radius: "<<r<<std::endl;
	stream<<"Axis: "<<axis.format(CommaInitFmt)<<std::endl<<std::endl;
}

ray_t transform_to_local(ray_t world_ray , Eigen::Vector3d center, Eigen::Vector3d axis) 
{
	ray_t local_ray(world_ray);
	Eigen::Vector3d z(0.0,0.0,1.0);
	Eigen::Vector3d v = axis.cross(z);
	double s = v.norm();
	double c = z.dot(axis);

	Eigen::Matrix3d V;
	V << 0.0, -v[2], v[1],
		 v[2], 0.0, -v[0],
		 -v[1], v[0], 0.0;
	Eigen::Matrix3d R = Matrix3d::Identity() + V + V*V*(1-c)/(s*s);

	local_ray.origin = (R * (world_ray.origin - center));
	local_ray.direction = (R * world_ray.direction).normalized();
	return local_ray;    
}