#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <utility>

#include "structures.h"

class Planets_class {
private:
	const planets planet_marks;
	std::vector<planet_data> Moon;
	std::vector<planet_data> Sun;
	const size_t MAX_ELEMENTS;
	double first_julian_day_Moon = 0.;
	double first_julian_day_Sun = 0.;

public:
	Planets_class(planets pl) : planet_marks(pl), MAX_ELEMENTS(270'000) {}

	void Fill_container(const std::string& str, std::vector<planet_data>& v) {
		std::string file_name = str;
		std::ifstream file(file_name);
		if (file.is_open()) {
			std::cout << "file " << file_name << " is opened" << std::endl;
			
			std::string buffer;
			while (!file.eof()) {
				file >> buffer;
				file >> buffer;

				planet_data tmp{0.,0.,0.,0., 0., 0., 0.};
				file >> tmp.date;

				file >> buffer;
				file >> buffer;
				file >> tmp.X;

				file >> buffer;
				file >> buffer;
				file >> tmp.Y;

				file >> buffer;
				file >> buffer;
				file >> tmp.Z;

				file >> buffer;
				file >> buffer;
				file >> tmp.R;

				file >> buffer;
				file >> buffer;
				file >> tmp.VX;

				file >> buffer;
				file >> buffer;
				file >> tmp.VY;

				file >> buffer;
				file >> buffer;
				file >> tmp.VZ;

				v.push_back(move(tmp));				
			}
 		}
		else {
			std::cout << "file " << file_name << " is not opened" << std::endl;
		}
	}

	void Fill_Moon(const std::string str) {
		Moon.reserve(MAX_ELEMENTS);
		Fill_container(str, Moon);
		first_julian_day_Moon = Moon[0].date;
	}

	void Fill_Sun(const std::string str) {
		Sun.reserve(MAX_ELEMENTS);
		Fill_container(str, Sun);
		first_julian_day_Moon = Sun[0].date;
	}

	My_vector Get_position(const double current_date) {


	}

	const std::vector<planet_data>& Get_vector_Moon() const {
		return Moon;
	}

	const std::vector<planet_data>& Get_vector_Sun() const {
		return Sun;
	}
};