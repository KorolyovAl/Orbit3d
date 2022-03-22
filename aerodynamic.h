#pragma once
#include "constants.h"
#include "structures.h"

const double an = 0.8; // coefficient of accommodation
const double at = 0.8; // coefficient of accommodation
double theta = 90.; // angle of attack, degree
const double Tw = 300.; // Wall temperature
const double k = 1.38E-23; // Boltzmann constant, J/K
const double atom_mass = 4.81E-26; // Air, kg
int case_aero = 0;

void Aerodynamic_check_case() {
	if (theta == 0. || theta == 90. || theta == 180. || theta == 270.) {
		case_aero = 0;
	}
	else {
		case_aero = 1;
	}
	theta = theta * (PI / 180.); // radian
}

double Calculate_Cx (const double& theta, const double& k, const double& atom_mass, const atmosphere& atmo, const double& Tw, const std::pair<My_vector, My_vector>& U) {
	// Алгоритм: Динамика разреженного газа, Коган М. Н. (на стр 345)
	double sin_theta = sin(theta);
	double sin2_theta = (1. - cos(2. * theta)) / 2.;
	double cos2_theta = (1. + cos(2. * theta)) / 2.;
	double h = atom_mass / (2. * k * atmo.tk); // atom_mass/(2*k*T)
	double S = U.second.Length() * sqrt(h); // скоростное соотношение
	double Stet = S * sin_theta; // S*sin(theta)
	double X1 = 2 * (2 - an) * sin2_theta / (S*sqrt(PI));
	double X2 = 2 * at * cos2_theta / (S*sqrt(PI));
	double X3 = an / S * sqrt(PI*Tw/atmo.tk) * sin2_theta;
	double X4 = (2 - an) * (sin2_theta + 1/(S*S) + at * cos2_theta);
	double X5 = erf(Stet);

	double Cx = ((X1 + X2) * exp(-1 * Stet * Stet) + X3 + 2 * sin_theta * X4 * X5);
	return Cx;
}
double Get_Cx(const atmosphere& atmo, const std::pair<My_vector, My_vector>& U) {
	double Cx = 0.;
	if (case_aero == 0) {
		Cx = Calculate_Cx((90. * PI / 180.), k, atom_mass, atmo, Tw, U);
	}
	else if (case_aero == 1) {
		double theta_1, theta_2;
		double Cx_1, Cx_2;
		theta_1 = (90. * PI / 180.) - theta;
		theta_2 = (0. * PI / 180.) + theta;
		Cx_1 = Calculate_Cx(theta_1, k, atom_mass, atmo, Tw, U);
		Cx_2 = Calculate_Cx(theta_2, k, atom_mass, atmo, Tw, U);

		Cx = Cx_1 + Cx_2;
	}
	return Cx;
}
