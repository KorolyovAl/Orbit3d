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

				planet_data tmp{0.,0.,0.,0.};
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
	}

	void Fill_Sun(const std::string str) {
		Sun.reserve(MAX_ELEMENTS);
		Fill_container(str, Sun);
	}

	void Get_position() {

	}
};