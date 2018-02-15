#include <resolver.hpp>

#include <rt.hpp>


using namespace rt;

void rt::render(const scene_t* scn)
{
	int num_samples = scn->img->get_num_samples_per_pixel();
	unsigned int w=scn->img->get_width();
	unsigned int h=scn->img->get_height();

	for (unsigned int i=0; i<w; i++)
	{
		for (unsigned int j=0; j<h; j++)
		{
			ray_t ray;
			int d=0;
			std::vector<Eigen::Vector2d> psamples=scn->img->sample_pixel(i,j);
			color_t color(0.0);
			for(int i=0;i<num_samples*num_samples;i++){
				color_t col = scn->cam->sample_ray(ray, psamples[i]);
				color += col * scn->intg->radiance(scn, ray, d);
			}
			color /= num_samples*num_samples;
			scn->img->set_pixel(i,j,color);
		}

	}
}

int main(int argc, char **argv)
{
	if (argc != 2) 
		{
        	std::cerr << "Syntax: " << argv[0] << " <scene.xml>" << std::endl;
        	return -1;
		}

	filesystem::path path(argv[1]);

	try
	{
		if (path.extension() == "xml") 
		{
			std::string scene_filename(argv[1]);
  			rt::scene_t scn(scene_filename);

  			rt::render(&scn);

			std::string img_filename = scene_filename;
    		size_t lastdot = img_filename.find_last_of(".");
    		
    		if (lastdot != std::string::npos) img_filename.erase(lastdot, std::string::npos);

   			img_filename += ".ppm";
  	
  			scn.img->write(img_filename);
  		}
  		else
  		{
  			std::cerr<<"Error: Unknown file type."<<std::endl;
  			return -1;
  		}
  } 
  catch (const std::exception &e)
  {
  	std::cerr<<"Error: "<<e.what()<<std::endl;
  	return -1;
  }
  	
  return 0;
}
