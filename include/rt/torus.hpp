#pragma once

#include <material.hpp>
#include <object.hpp>
#include <ray.hpp>
#include <utils.hpp>


namespace rt
{
	/**
	 * \brief The torus object class.
	 **/
	class torus_t : public object_t
	{
	private:
		/// Torus center
		Eigen::Vector3d center;

		/// Torus outer radius
		double R;

		/// Torus inner radius
		double r;

		//// Torus axis direction 
		Eigen::Vector3d axis;

		/// torus material
		material_t* mat;

	public:
		/// Constructor
		torus_t(material_t* _mat);
		/// Constructor
		torus_t(material_t* _mat, Eigen::Vector3d _c, double _R, double _r, Eigen::Vector3d axis);
		/// Destuctor
		virtual ~torus_t();

		/// Returns the mandatory object name
		std::string get_name(void) const { return std::string("torus"); }

		/**
		* Returns true if the _ray hits this object. The hit information is returned in result. 
		* This is not valid if there is no intersection and the function returns false.
		**/
		bool intersect(hit_t& result, const ray_t& _ray) const;

		/// Returns the normal to the surface at point _p.
		Eigen::Vector3d get_normal(Eigen::Vector3d& _p) const;

		/// Returns the material for the torus.
		material_t* get_material(void) const;

		/// Prints information about the torus. to stream.
		virtual void print(std::ostream &stream) const;
	};
}
