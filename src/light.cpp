#include <iostream>
#include <light.hpp>

using namespace rt;

color_t Vector3d_to_colour_t(const Vector3d);
color_t multiply_color_t_Vector3d(color_t col, Vector3d vec); 

light_t::light_t() { }
light_t::~light_t() { }


point_light_t::point_light_t(const Vector3d& _pos, const Vector3d& _col, const double _ka): pos(_pos), col(_col), ka(_ka) 
{ }

point_light_t::~point_light_t()
{ }

color_t point_light_t::direct(const Vector3d& hitpt, const Vector3d& normal, const material_t* mat, const scene_t* scn) const
{

	//printf("%f\n", (hitpt- Vector3d(0.0,100003.60,35.0)).norm() );
	color_t kd = mat -> get_diffuse();
	color_t ks = mat -> get_specular();


	Vector3d incident  = (pos - hitpt).normalized();
	Vector3d reflected =  (2 * incident.dot(normal) * normal - incident).normalized();
	Vector3d eye = scn -> cam -> get_eye();
	Vector3d view = (eye - hitpt).normalized();

	ray_t shadow_ray = ray_t(hitpt + BIAS * normal, incident);
	shadow_ray.maxt = (pos - hitpt).norm();

	bool in_shadow = false;
	std::vector<object_t*> objects = scn -> objs;
	for(unsigned int  i=0;i<objects.size();i++) {
		hit_t result;
		in_shadow = in_shadow || objects[i]->intersect(result,shadow_ray);
	}

	color_t ambient = multiply_color_t_Vector3d(kd, Vector3d_to_colour_t(ka * col));

	color_t diffuse(0.0);
	color_t specular(0.0);

	if(!in_shadow) {
		diffuse = multiply_color_t_Vector3d(kd, col * incident.dot(normal)).clamp();
		specular = multiply_color_t_Vector3d(ks, col * pow(std::max(double(reflected.dot(view)),0.0),mat -> get_shininess()));
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

area_light_t::area_light_t(const Vector3d& _center, const double& _radius, const Vector3d& _normal, const int& _num_samples, const double _ka)
{
	center = _center;
	radius = _radius;
	normal = _normal.normalized();
	num_samples = _num_samples;
	col = _col;
	ka = _ka;
}

color_t area_light_t::direct(const Vector3d& _hitpt, const Vector3d& _normal, const material_t* _mat, const scene_t* _scn) const
{
	std::vector<Eigen::Vector3d> sample_points = sample(num_samples);

	color_t color(0.0);

	for(int i =0 ; i < num_samples ; i++)
	{
		point_light_t * point_light = new point_light_t(sample_points[i], col, ka);
		color += point_light->direct(_hitpt, _normal, _mat, _scn);
		delete point_light;
	}

	color /= num_samples;
	return color;
}

std::vector<Eigen::Vector3d> area_light_t::sample(int _num_samples){
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	std::vector<Eigen::Vector3d> samples;
	for(int i=0;i<_num_samples;i++){
		double theta = 2 * M_PI * distribution(generator);
		double r = sqrt(distribution(generator)) * radius;
		Vector3d e_x = center.cross(normal).normalized();
		Vector3d e_y = normal.cross(e_x).normalized();	
		samples.push_back(center + r*cos(theta)*e_x + r*sin(theta)*e_y);
	}
	return samples;
}

void area_light_t::print(std::ostream &stream) const
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Light Properties: -------------------------------"<<std::endl;
	stream<<"Type: Area Light"<<std::endl;
	stream<<"Center: "<<center.format(CommaInitFmt)<<std::endl;
	stream<<"Radius: "<<radius<<std::endl;
	stream<<"Normal: "<<normal.format(CommaInitFmt)<<std::endl;
	stream<<"Color: "<<col.format(CommaInitFmt)<<std::endl;
	stream<<"Ambient Coefficient: "<<ka<<std::endl<<std::endl;
}

color_t Vector3d_to_colour_t(const Vector3d vec) {
	return color_t(vec[0],vec[1],vec[2]);
}

color_t multiply_color_t_Vector3d(color_t col, Vector3d vec) {
	return color_t(col.r()*vec[0], col.g()*vec[1], col.b()*vec[2]);
}