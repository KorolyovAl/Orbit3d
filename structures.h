#pragma once
#include <vector>
#include <iostream>

class My_vector {
public:
	My_vector() {
		x = 0;
		y = 0;
		z = 0;
	}
	My_vector(double x_out, double y_out, double z_out) {
		x = x_out;
		y = y_out;
		z = z_out;
	}
	std::vector<double> GetCoord() const {
		return { GetX(), GetY(), GetZ() };
	}
	double GetX() const {
		return x;
	}
	double GetY() const {
		return y;
	}
	double GetZ() const {
		return z;
	}
	double Length() const {
		return sqrt(x * x + y * y + z * z);
	}
	My_vector Normalize() const {
		double l = Length();
		return { x / l, y / l, z / l };
	}
	void Reverse() {
		x *= -1;
		y *= -1;
		z *= -1;
	}
	My_vector Get_sqr() const {
		return { x * x, y * y, z * z };
	}
	void Print() const {
		std::cout << "x: " << x << " y: " << y << " z: " << z << std::endl;
	}
private:
	double x;
	double y;
	double z;
};

//---------KEPLER_ELEMENTS_OF_ORBIT---------
struct elements_of_orbit {
	double orb_eccentricity = 0.;
	double apoapsis = 0.; // km, from surface
	double periapsis = 0.;
	double periapsis_argument = 0.;
	double true_anomaly = 0.;
	double inclination = 0.;
	double ascending_long = 0.;
	double orb_major_ax = 0.;
	double orb_minor_ax = 0.;
};
//-----------------------------------

struct satellite {
	double mass = 0.;
	double S = 0.;// pi*d*d/4 for sphere, d*d for cubesat
	double Cx = 0.;
	double d = 0.;
	double thrust = 1.9215; // newton
	bool thrust_on = false;
	int thrust_way = 1;
};

struct atmosphere {
	double f107;   // F10.7 solar activity index
	int day;    // day number sins 1,jan
	double d; // number density [1/m^3]
	double ro;  // density [kg/m^3]
	double tk;  // temperture [K]
	double pres;  // pressure [Pa]
	double mass;  // molecular mass [ce]
	double gamma; // cp/cv ratio
	double sv; // sound velosity [m/c]
	double lambda; // mean free path [m]

	atmosphere() {
		f107 = 0.;
		day = 0.;
		d = 0.;
		ro = 0.;
		tk = 0.;
		pres = 0.;
		mass = 0.;
		gamma = 0.;
		sv = 0.;
		lambda = 0.;
	}
};

struct planets {
	bool general_mark;
	std::pair<bool, std::string> Moon;
	std::pair<bool, std::string> Sun;
};

struct planet_data {
	double date;
	double X;
	double Y;
	double Z;
	double R;
	double VX;
	double VY;
	double VZ;
};

struct date_struct {
	int day;
	int month;
	int year;
};