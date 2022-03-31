#pragma once
#include <iostream>
#include <vector>
#include <fstream>

#include "structures.h"

class Planets_class {
private:
	const planets planet_marks;
	std::vector<My_vector> Moon;
	std::vector<My_vector> Sun;

public:
	Planets_class(planets pl) : planet_marks(pl) {}

	void Fill_container(const std::string& str, std::vector<My_vector>& v) {
		std::string file_name = str;
		std::ifstream file(file_name);
		if (file.is_open()) {
			std::cout << "file " << file_name << " is opened" << std::endl;
 		}
		else {
			std::cout << "file " << file_name << " is not opened" << std::endl;
		}
	}

	void Fill_Moon(const std::string str) {
		Fill_container(str, Moon);
	}

	void Fill_Sun(const std::string str) {
		Fill_container(str, Sun);
	}

	void Get_position() {

	}
};