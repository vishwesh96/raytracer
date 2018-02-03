#include <iostream>

#include <integrator.hpp>

using namespace rt;

color_t whitted_integrator_t::radiance(const scene_t* _scn, ray_t& _ray, int& d) const
{
	const float BIAS = 0.5f;
	if(d >= this->depth) 
		return _scn->img->get_bgcolor();
	int depth = d+1;
	bool found_intersection=false;
	std::vector<object_t*>::const_iterator oit;
	hit_t hit, minhit;
	Eigen::Vector3f hitpt, normal;

	for (oit=_scn->objs.begin(); oit!=_scn->objs.end(); oit++)
	{
		if ((*oit)->intersect(hit, _ray))
		{
		  _ray.maxt = hit.second;
		  minhit = hit;
		  
		  hitpt = _ray.origin+_ray.maxt*_ray.direction;
		  normal = (*oit)->get_normal(hitpt);

		  found_intersection=true;
		}
		
	}
	
	color_t d_col(0.0);
	if(found_intersection)
	{
		simplemat_t * mat = (simplemat_t *)minhit.first->get_material();
		std::list<light_t*>::const_iterator lit;
		for(lit=_scn->lits.begin(); lit!=_scn->lits.end(); lit++)
		{
			// d_col += multiply_color_t_color_t(mat->get_diffuse(),(*lit)->direct(hitpt, normal, mat , _scn));
			d_col += (*lit)->direct(hitpt, normal, mat , _scn);
		}

		Vector3f incident = _ray.direction.normalized();
		float incident_dot_normal = incident.dot(normal);
		
		Vector3f transmitted;
		ray_t transmitted_ray;
		float nr=1.0;

		// if(incident_dot_normal < 0.0){		//entering the object
		// 	nr = 1.0/mat->get_eta();
		// } else {							//exiting the object
		// 	nr = mat->get_eta();
		// 	normal = -normal;
		// 	incident_dot_normal = -incident_dot_normal;
		// }

		if(incident_dot_normal > 0.0){		//entering the object
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

		Vector3f reflected =  (incident - 2 * incident_dot_normal * normal).normalized();
		ray_t reflected_ray(hitpt + BIAS * normal,reflected);

		if(mat->get_is_reflect())
			d_col += mat->get_reflect() * radiance(_scn,reflected_ray,depth);

		if(mat->get_is_transmit()){
			float D = 1 - (nr * nr *(1 - incident_dot_normal*incident_dot_normal));
			if(D >= 0.0){				//normal refraction
				transmitted = (nr * incident - (nr * incident_dot_normal + sqrt(D)) * normal).normalized();
				transmitted_ray = ray_t(hitpt - BIAS * normal ,transmitted);
				d_col += mat->get_transmit() * radiance(_scn,transmitted_ray,depth);
			} else {					//TIR
				if(mat->get_is_reflect())
					d_col += (color_t(1.0) - mat->get_reflect()) * radiance(_scn,reflected_ray,depth);
				else 
					d_col +=  radiance(_scn,reflected_ray,depth);
			}
		}
	}
	else d_col = _scn->img->get_bgcolor();

	return d_col;
}
