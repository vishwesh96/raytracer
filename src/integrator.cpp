#include <iostream>

#include <integrator.hpp>

using namespace rt;
Eigen::Vector3d uniform_sample_hemisphere(double u, double v);
Eigen::Vector3d cosine_sample_hemisphere(double u, double v);
Eigen::Vector3d specular_sample_hemisphere(double u, double v, int n);
Eigen::Vector3d transform_to_world(Eigen::Vector3d& local_direction, Eigen::Vector3d normal);
double color_norm(color_t col); 

color_t whitted_integrator_t::radiance(const scene_t* _scn, ray_t& _ray, int d) const
{
	if(d >= this->depth) 
		return _scn->img->get_bgcolor();
	int depth = d+1;
	bool object_intersection=false;
	std::vector<object_t*>::const_iterator oit;
	hit_t hit, minhit;
	Eigen::Vector3d hitpt, normal;

	for (oit=_scn->objs.begin(); oit!=_scn->objs.end(); oit++)
	{
		if ((*oit)->intersect(hit, _ray))
		{
		  _ray.maxt = hit.second;
		  minhit = hit;
		  
		  hitpt = _ray.origin+_ray.maxt*_ray.direction;
		  normal = (*oit)->get_normal(hitpt);
		  object_intersection=true;
		}
		
	}
	
	color_t d_col(0.0);
	if(object_intersection)
	{
		simplemat_t * mat = (simplemat_t *)minhit.first->get_material();
		std::list<light_t*>::const_iterator lit;
		for(lit=_scn->lits.begin(); lit!=_scn->lits.end(); lit++)
		{
			d_col += (*lit)->direct(hitpt, normal, mat , _scn);
		}

		Eigen::Vector3d incident = _ray.direction.normalized();
		double incident_dot_normal = incident.dot(normal);
		
		Eigen::Vector3d transmitted;
		ray_t transmitted_ray;
		double nr=1.0;

		if(incident_dot_normal > 0.0){		
			if(mat->get_is_transmit()){
				nr = mat->get_eta();
			}
			normal = -normal;
			incident_dot_normal = -incident_dot_normal;
		} else {
			if(mat->get_is_transmit()){
				nr = 1.0/mat->get_eta();
			}
		}

		Eigen::Vector3d reflected =  (incident - 2 * incident_dot_normal * normal).normalized();
		ray_t reflected_ray(hitpt + BIAS * normal,reflected);


		color_t reflect_radiance = color_t(0.0);

		if (mat->get_is_reflect() || mat->get_is_transmit()) reflect_radiance = radiance(_scn,reflected_ray,depth);

		if(mat->get_is_reflect())	// if reflecting surface
			d_col += mat->get_reflect() * reflect_radiance;

		if(mat->get_is_transmit()){ 
			double D = 1 - (nr * nr *(1 - incident_dot_normal*incident_dot_normal));

			if(D >= 0.0){		// Transmission
				transmitted = (nr * incident - (nr * incident_dot_normal + sqrt(D)) * normal).normalized();
				transmitted_ray = ray_t(hitpt - BIAS * normal ,transmitted);
				d_col += mat->get_transmit() * radiance(_scn,transmitted_ray,depth);
			} else {			// TIR
				if(mat->get_is_reflect())	// if reflecting remove the previous one added and add 1.0
					d_col += (color_t(1.0) - mat->get_reflect()) * reflect_radiance;
				else 
					d_col +=  reflect_radiance;
			}
		}
	}
	else d_col = _scn->img->get_bgcolor();

	return d_col;
}

color_t monte_carlo_integrator_t::radiance(const scene_t* _scn, ray_t& _ray, int d) const
{
	if(d >= this->depth) 
		return _scn->img->get_bgcolor();
	int depth = d+1;
	bool object_intersection=false;
	std::vector<object_t*>::const_iterator oit;
	hit_t hit, minhit;
	Eigen::Vector3d hitpt, normal;
	std::default_random_engine generator(std::random_device{}());
	std::uniform_real_distribution<double> distribution(0.0,1.0);

	for (oit=_scn->objs.begin(); oit!=_scn->objs.end(); oit++)
	{
		if ((*oit)->intersect(hit, _ray))
		{
		  _ray.maxt = hit.second;
		  minhit = hit;
		  
		  hitpt = _ray.origin+_ray.maxt*_ray.direction;
		  normal = (*oit)->get_normal(hitpt);
		  object_intersection=true;
		}
		
	}

	// Area light intersection, return area light color
	bool light_intersection=false;
	light_hit_t light_hit;
	for(std::list<light_t*>::const_iterator lit = _scn->lits.begin(); lit!=_scn->lits.end(); lit++)
	{
		if((*lit)->intersect(light_hit, _ray))
		{
			_ray.maxt = light_hit.second;
			light_intersection=true;

		}
	}

	if(light_intersection)
	{
		Eigen::Vector3d col = light_hit.first->get_color();
		return color_t(col[0],col[1],col[2]);
	}

	color_t indirect_illumination(0.0);
	color_t direct_illumination(0.0);

	if(object_intersection)
	{
		simplemat_t * mat = (simplemat_t *)minhit.first->get_material();
		color_t kd = mat -> get_diffuse();
		color_t ks = mat -> get_specular();
		int n = mat -> get_shininess();
		Eigen::Vector3d incident = _ray.direction.normalized();
		double incident_dot_normal = incident.dot(normal);
		Eigen::Vector3d reflected =  (incident - 2 * incident_dot_normal * normal).normalized();
		Eigen::Vector3d transmitted;
		ray_t transmitted_ray;
		double nr=1.0;
		bool tir = false;

		if(incident_dot_normal > 0.0){		
			if(mat->get_is_transmit()){
				nr = mat->get_eta();
			}
			normal = -normal;
			incident_dot_normal = -incident_dot_normal;
		} else {
			if(mat->get_is_transmit()){
				nr = 1.0/mat->get_eta();
			}
		}

		double D = 1 - (nr * nr *(1 - incident_dot_normal*incident_dot_normal));

		if(D >= 0.0){		// Transmission
			transmitted = (nr * incident - (nr * incident_dot_normal + sqrt(D)) * normal).normalized();
		} else {			// TIR
			tir = true;
		}
	
		// std::list<light_t*>::const_iterator lit;
		// for(lit=_scn->lits.begin(); lit!=_scn->lits.end(); lit++)
		// {
		// 	direct_illumination += (*lit)->direct(hitpt, normal, mat , _scn);
		// }
		
		if(mat->get_is_reflect() || tir) 
		{
			ray_t reflected_ray(hitpt + BIAS*normal, reflected);
			return radiance(_scn,reflected_ray,depth);
		}

		if(mat->get_is_transmit()) 
		{
			transmitted_ray = ray_t(hitpt - BIAS * normal ,transmitted);
			return radiance(_scn,transmitted_ray,depth);

		}
		
		double t = distribution(generator);
		double p = distribution(generator);
		double q = distribution(generator);
		Eigen::Vector3d local_direction;
		Eigen::Vector3d sample_brdf_direction;
		Eigen::Vector3d sample_btdf_direction;
		int diffuse = 0.0, specular = 0.0;
		if( t < color_norm(kd)/(color_norm(kd) + color_norm(ks)))
		{
			diffuse = 1.0;
			local_direction = cosine_sample_hemisphere(p,q);
			sample_brdf_direction = transform_to_world(local_direction, normal);
			sample_btdf_direction = transform_to_world(local_direction, -normal);
		} 
		else
		{
			specular = 1.0;
			local_direction = specular_sample_hemisphere(p,q,n);
			sample_brdf_direction = transform_to_world(local_direction, reflected);
			sample_btdf_direction = transform_to_world(local_direction, transmitted);
		}

		double cos_theta_brdf = sample_brdf_direction.dot(normal);
		double cos_theta_btdf = sample_btdf_direction.dot(-normal);
		color_t reflectance = mat->get_reflect();
		color_t transmittance = mat->get_transmit();
		ray_t sampled_brdf_ray(hitpt + BIAS*normal, sample_brdf_direction);
		ray_t sampled_btdf_ray(hitpt + BIAS*normal, sample_btdf_direction);

		color_t brdf_illumination = color_t(0.0);
		color_t btdf_illumination = color_t(0.0);

		if(reflectance.r() >= EPSILON || reflectance.g() >= EPSILON || reflectance.b() >= EPSILON)
		{
			brdf_illumination =  radiance(_scn, sampled_brdf_ray, depth);
		}
		if(transmittance.r() >= EPSILON || transmittance.g() >= EPSILON || transmittance.b() >= EPSILON)
		{
			btdf_illumination = radiance(_scn, sampled_btdf_ray, depth);
		}
		
		indirect_illumination += reflectance * brdf_illumination * (diffuse * kd * M_PI + specular * ks * M_PI * (n+2)/(n+1) * cos_theta_brdf) 
							 + transmittance * btdf_illumination * (diffuse * kd * M_PI + specular * ks * M_PI * (n+2)/(n+1) * cos_theta_btdf);

	}
	else return _scn->img->get_bgcolor();

	return indirect_illumination;

}

Eigen::Vector3d uniform_sample_hemisphere(double u, double v)
{
    double r = sqrt(1.0f - u * u);
    double phi = 2 * M_PI * v;
    return Eigen::Vector3d(cos(phi) * r, sin(phi) * r, u).normalized();
}

Eigen::Vector3d cosine_sample_hemisphere(double u, double v)
{

	double r = sqrt(u);
    double theta = 2 * M_PI * v;
    return Eigen::Vector3d(r * cos(theta), r * sin(theta), sqrt(std::max(0.0, 1 - u))).normalized();
}

Eigen::Vector3d specular_sample_hemisphere(double u, double v, int n)
{

	double r = sqrt(1 - pow(u, 2.0/(n+1)));
    double theta = 2 * M_PI * v;
    return Eigen::Vector3d(r * cos(theta), r * sin(theta), sqrt(pow(u, 1.0/(n+1)))).normalized();
}

Eigen::Vector3d transform_to_world(Eigen::Vector3d& local_direction, Eigen::Vector3d normal) 
{
	Eigen::Vector3d z(0.0,0.0,1.0);
	Eigen::Vector3d v = z.cross(normal);
	double s = v.norm();
	double c = z.dot(normal);

	Eigen::Matrix3d V;
	V << 0.0,-v[2],v[1],v[2],0.0,-v[0],-v[1],v[0],0;
	Eigen::Matrix3d R = Matrix3d::Identity() + V + V*V*(1-c)/(s*s);

	return (R * local_direction).normalized();
}

double color_norm(color_t col)
{
	return sqrt(pow(col.r(),2) + pow(col.g(),2) + pow(col.b(),2));
}