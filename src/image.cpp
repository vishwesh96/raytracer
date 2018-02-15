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

image_t::image_t(int _w, int _h, color_t _bgc):width(_w),height(_h),bgcolor(_bgc)
{
	aspect = float(width)/float(height);
	data = new char[width*height*3]; 
}

image_t::~image_t()
{ 
	delete data;
}

int image_t::get_width(void) const {return width; }
int image_t::get_height(void) const {return height; }
float image_t::get_aspect(void) const {return aspect; }

color_t image_t::get_bgcolor(void) const {return bgcolor; }

std::vector<Eigen::Vector2d> image_t::sample_pixel(unsigned int _x, unsigned int _y, int num_samples) const
{
	//grid
	// std::vector<Eigen::Vector2d> samples;
	// float start_x = float(_x)/width + 1.0/(2*num_samples*width);
	// float start_y = float(_y)/height + 1.0/(2*num_samples*height);
	// for(int i=0;i<num_samples;i++){
	// 	for(int j=0;j<num_samples;j++){
	// 		float center_x = start_x + float(i)/(num_samples*width);
	// 		float center_y = start_y + float(j)/(num_samples*height);
	// 		samples.push_back(Eigen::Vector2d(center_x, center_y));
	// 	}
	// }
	// return samples;
	//jitter
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution_x(-1.0/(2*num_samples*width),1.0/(2*num_samples*width));
	std::uniform_real_distribution<double> distribution_y(-1.0/(2*num_samples*height),1.0/(2*num_samples*height));

	std::vector<Eigen::Vector2d> samples;
	float start_x = float(_x)/width + 1.0/(2*num_samples*width);
	float start_y = float(_y)/height + 1.0/(2*num_samples*height);
	for(int i=0;i<num_samples;i++){
		for(int j=0;j<num_samples;j++){
			float center_x = start_x + float(i)/(num_samples*width);
			float center_y = start_y + float(j)/(num_samples*height);
			samples.push_back(Eigen::Vector2d(center_x + distribution_x(generator), center_y + distribution_y(generator)));
		}
	}
	return samples;
}

color_t image_t::get_pixel(unsigned int _x, unsigned int _y) const
{
	int pos=(_y)*width*3+(_x)*3;
	float r=float(data[pos])/255.0f;
	float g=float(data[pos+1])/255.0f;
	float b=float(data[pos+2])/255.0f;
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