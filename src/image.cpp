/*
    This file is part of rt.

    rt is a simple ray tracer meant to be used for teaching ray tracing.

    Copyright (c) 2018 by Parag Chaudhuri

	Some parts of rt are derived from Nori by Wenzel Jacob.
	https://github.com/wjakob/nori/

    rt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    rt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <image.hpp>

using namespace rt;

image_t::image_t(int _w, int _h, int _samples, color_t _bgc):width(_w),height(_h),samples_per_pixel(_samples),bgcolor(_bgc)
{
	aspect = double(width)/double(height);
	data = new char[width*height*3]; 
}

image_t::~image_t()
{ 
	delete data;
}

int image_t::get_width(void) const {return width; }
int image_t::get_height(void) const {return height; }
int image_t::get_num_samples_per_pixel(void) const {return samples_per_pixel;}
double image_t::get_aspect(void) const {return aspect; }

color_t image_t::get_bgcolor(void) const {return bgcolor; }

std::vector<Eigen::Vector2d> image_t::sample_pixel(unsigned int _x, unsigned int _y) const
{
	//grid
	// std::vector<Eigen::Vector2d> samples;
	// float start_x = float(_x)/width + 1.0/(2*samples_per_pixel*width);
	// float start_y = float(_y)/height + 1.0/(2*samples_per_pixel*height);
	// for(int i=0;i<samples_per_pixel;i++){
	// 	for(int j=0;j<samples_per_pixel;j++){
	// 		float center_x = start_x + float(i)/(samples_per_pixel*width);
	// 		float center_y = start_y + float(j)/(samples_per_pixel*height);
	// 		samples.push_back(Eigen::Vector2d(center_x, center_y));
	// 	}
	// }
	// return samples;

	//jitter
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution_x(-1.0/(2*samples_per_pixel*width),1.0/(2*samples_per_pixel*width));
	std::uniform_real_distribution<double> distribution_y(-1.0/(2*samples_per_pixel*height),1.0/(2*samples_per_pixel*height));

	std::vector<Eigen::Vector2d> samples;
	float start_x = float(_x)/width + 1.0/(2*samples_per_pixel*width);
	float start_y = float(_y)/height + 1.0/(2*samples_per_pixel*height);
	for(int i=0;i<samples_per_pixel;i++){
		for(int j=0;j<samples_per_pixel;j++){
			float center_x = start_x + float(i)/(samples_per_pixel*width);
			float center_y = start_y + float(j)/(samples_per_pixel*height);
			samples.push_back(Eigen::Vector2d(center_x + distribution_x(generator), center_y + distribution_y(generator)));
		}
	}
	return samples;
}

color_t image_t::get_pixel(unsigned int _x, unsigned int _y) const
{
	int pos=(_y)*width*3+(_x)*3;
	double r=double(data[pos])/255.0;
	double g=double(data[pos+1])/255.0;
	double b=double(data[pos+2])/255.0;
	return color_t(r,g,b);
}

void image_t::set_pixel(unsigned int _x, unsigned int _y, color_t _col)
{
	int pos=(_y)*width*3+(_x)*3;
	char r = to_char(_col.r());
	char g = to_char(_col.g());
	char b = to_char(_col.b());
	data[pos]=r; data[pos+1]=g; data[pos+2]=b;
}

void image_t::write(std::string filename)
{
	std::ofstream out(filename.c_str(), std::ios::binary|std::ios::out);
	out<<"P6"<<std::endl<<width<<" "<<height<<" "<<255<<std::endl;
	out.write((const char*)data,width*height*3);
	out.close();
}

void image_t::print(std::ostream &stream)
{
	Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "[ ", " ]");
	
	stream<<"Image Properties: -------------------------------"<<std::endl;
	stream<<"BG Color: "<<bgcolor.format(CommaInitFmt)<<std::endl;
	stream<<"Width: "<<width<<std::endl;
	stream<<"Height:"<<height<<std::endl;
	stream<<"aspect:"<<aspect<<std::endl<<std::endl;
}