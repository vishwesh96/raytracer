#include <sphere.hpp>

using namespace rt;

sphere_t::sphere_t(material_t* _mat):center(0.0,0.0,0.0),radius(1.0),mat(_mat) { }
sphere_t::sphere_t(material_t* _mat, Eigen::Vector3f _c, float _r): center(_c), radius(_r), mat(_mat) { }

sphere_t::~sphere_t() { }

bool sphere_t::intersect(hit_t& result, const ray_t& _ray) const
{
	Vector3f r2c = center - _ray.origin;
	const float b = r2c.dot(_ray.direction);
	float d = b*b - r2c.dot(r2c) + radius*radius;

	if (d < 0)
		return false;
	else
		d = sqrt(d);

	float t;

	t=b-d;

	if (!is_zero(t) && t > _ray.mint)
	{
		result = hit_t(this, t);
	}
	else 
	{
		t = b+d;
		if (!is_zero(t) && t > _ray.mint) 
			result = hit_t(this,t);
		else return false;
	}

	return t >= _ray.mint && t <= _ray.maxt;

}

Eigen::Vector3f sphere_t::get_normal(Eigen::Vector3f& _p) const
{
	Eigen::Vector3f normal = _p - center;
	normal.normalize();

	return normal;
}

material_t* sphere_t::get_material(void) const
{
	return mat;
}

void sphere_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Object Properties: -------------------------------"<<std::endl;
	stream<<"Type: Sphere"<<std::endl;
	stream<<"Center: "<<center.format(CommaInitFmt)<<std::endl;
	stream<<"Radius: "<<radius<<std::endl<<std::endl;
}