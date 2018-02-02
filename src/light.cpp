#include <iostream>
#include <light.hpp>

using namespace rt;

color_t vector3f_to_colour_t(const Vector3f);
color_t multiply_color_t_vector3f(color_t col, Vector3f vec); 

light_t::light_t() { }
light_t::~light_t() { }


point_light_t::point_light_t(const Vector3f& _pos, const Vector3f& _col, const float _ka): pos(_pos), col(_col), ka(_ka) 
{ }

point_light_t::~point_light_t()
{ }

color_t point_light_t::direct(const Vector3f& hitpt, const Vector3f& normal, const material_t* mat, const scene_t* scn) const
{
	color_t kd = mat -> get_diffuse();
	color_t ks = mat -> get_specular();

	Vector3f incident  = (pos - hitpt).normalized();
	Vector3f reflected =  (2 * incident.dot(normal) * normal - incident).normalized();
	Vector3f eye = scn -> cam -> get_eye();
	Vector3f view = (eye - hitpt).normalized();

	ray_t shadow_ray = ray_t(hitpt + BIAS * normal, incident);
	bool in_shadow = false;

	std::vector<object_t*> objects = scn -> objs;
	for(unsigned int  i=0;i<objects.size();i++) {
		hit_t result;
		in_shadow = in_shadow || objects[i]->intersect(result,shadow_ray);
	}

	color_t ambient = vector3f_to_colour_t(ka * col);
	color_t diffuse(0.0);
	color_t specular(0.0);

	if(!in_shadow) {
		diffuse = multiply_color_t_vector3f(kd, col * incident.dot(normal)).clamp();
		specular = multiply_color_t_vector3f(ks,col * pow(std::max(double(reflected.dot(view)),0.0),mat -> get_shininess()));
	}
	color_t total = ambient + diffuse + specular;
	return total;
}
		

void point_light_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Light Properties: -------------------------------"<<std::endl;
	stream<<"Type: Point Light"<<std::endl;
	stream<<"Position: "<<pos.format(CommaInitFmt)<<std::endl;
	stream<<"Color: "<<col.format(CommaInitFmt)<<std::endl;
	stream<<"Ambient Coefficient: "<<ka<<std::endl<<std::endl;
}

color_t vector3f_to_colour_t(const Vector3f vec) {
	return color_t(vec[0],vec[1],vec[2]);
}

color_t multiply_color_t_vector3f(color_t col, Vector3f vec) {
	return color_t(col.r() * vec[0], col.g() * vec[1], col.b() * vec[2]);
}