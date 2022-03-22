#pragma once
#include "structures.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cctype>

void reader (struct elements_of_orbit& _orb, struct satellite& _sat) {
	std::string buffer;
	std::string line;

	std::ifstream file("initial_data.txt");
		if (file.is_open()) {
			while (!file.eof()) { // 0 - is not the end of file, 1 - is it.
				file >> buffer;
				line += buffer + " ";
				if (line == "Parameters of orbit: " || line == "Parameters of satellite: ") {
					line.clear();
				}
				else if (line == "inclination ") {
					double tmp = 0;
					file >> tmp;
					_orb.inclination = tmp;
					line.clear();
				}
				else if (line == "orbit eccentricity ") {
					double tmp = 0;
					file >> tmp;
					_orb.orb_eccentricity = tmp;
					line.clear();
				}
				else if (line == "apoapsis ") {
					double tmp = 0;
					file >> tmp;
					_orb.apoapsis = tmp;
					line.clear();
				}
				else if (line == "periapsis ") {
					double tmp = 0;
					file >> tmp;
					_orb.periapsis = tmp;
					line.clear();
				}
				else if (line == "argument of periapsis ") {
					double tmp = 0;
					file >> tmp;
					_orb.periapsis_argument = tmp;
					line.clear();
				}
				else if (line == "true anomaly ") {
					double tmp = 0;
					file >> tmp;
					_orb.true_anomaly = tmp;
					line.clear();
				}
				else if (line == "ascending long ") {
					double tmp = 0;
					file >> tmp;
					_orb.ascending_long = tmp;
					line.clear();
				}
				else if (line == "orbital major axis ") {
					double tmp = 0;
					file >> tmp;
					_orb.orb_major_ax = tmp;
					line.clear();
				}
				else if (line == "orbital minor axis ") {
					double tmp = 0;
					file >> tmp;
					_orb.orb_minor_ax = tmp;
					line.clear();
				}
				else if (line == "mass ") {
					double tmp = 0;
					file >> tmp;
					_sat.mass = tmp;
					line.clear();
				}
				else if (line == "area ") {
					double tmp = 0;
					file >> tmp;
					_sat.S = tmp;
					line.clear();
				}
				else if (line == "size ") {
					double tmp = 0;
					file >> tmp;
					_sat.d = tmp;
					line.clear();
				}
				else if (line == "Cx ") {
					double tmp = 0;
					file >> tmp;
					_sat.Cx = tmp;
					line.clear();
				}
			}
		}
		else std::cerr << "File of initial data is not opened" << std::endl;

}


