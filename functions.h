#pragma once
#include "structures.h"
#include "constants.h"
#include "myvector_operators.h"
#include "planets.h"

My_vector Calculate_coordinates_periapsis(const double& latitude, const double& longitude, const double& hight, const double& e) {
	double L = longitude * PI / 180;
	double B = latitude * PI / 180;
	double H = hight * 1'000.; // to meters

	double N = major_ax / sqrt(1 - e * e * sin(B) * sin(B));
	double x = (N + H) * cos(B) * cos(L);
	double y = (N + H) * cos(B) * sin(L);
	double z = (N + H - N * e * e) * sin(B);
	return { x, y, z };
}

std::vector<double> Calculate_elliptic_coordinates(const My_vector& coord) {
	/*
	Алгоритм взят отсюда: специальные методы прикладной геодезии, А. О. Куприянов, А. С. Корчагин, Д. А. Морозов. (на стр 6)
	Algorithm: special methods of applied geodesy, A. O. Kupriyanov, A. S. Korchagin, D. A. Morozov. (page 6)
	*/
	double x = coord.GetX();
	double y = coord.GetY();
	double z = coord.GetZ();

	double latitude;
	double longitude;
	double hight;

	double Q = sqrt(x * x + y * y); //the radius of the parallel
	double B1 = atan(z / (Q * (1 - eccentricity * eccentricity)));

	double N;
	double B_next;
	double B_old = B1;
	double dB;

	for (int i = 0; i < 5; i++) {
		N = major_ax / sqrt(1 - eccentricity * eccentricity * sin(B_old) * sin(B_old));
		hight = Q / cos(B_old) - N;
		B_next = atan(z / (Q * (1 - N * eccentricity * eccentricity / (N + hight))));
		dB = B_next * 180.0 / PI - B_old * 180.0 / PI;
		B_old = B_next;
		if (abs(dB) < 1e-8) break;
	}

	longitude = asin(y / Q) * 180.0 / PI;
	latitude = B_next * 180.0 / PI;

	//	cout << " L = " << longitude << " B = " << latitude << " H, km = " << hight / 1000 << endl;

	if (x < 0. && y > 0.) {
		longitude = 180. - longitude;
	}
	else if (x < 0. && y < 0.) {
		longitude = -180. - longitude;
	}


	return { longitude, latitude, hight };
}

My_vector Get_g(const My_vector& coord) {
	My_vector g = coord.Normalize();
	g.Reverse();
	double g_number = GE / (coord.Length() * coord.Length());
	return { g.GetX() * g_number, g.GetY() * g_number, g.GetZ() * g_number };
}

std::pair<My_vector, My_vector> Calculate_state_vector(struct elements_of_orbit& orb) {
	/*
	Алгоритм взят отсюда: В. И. Крылов, основы теории движения исз (часть первая: невозмущённое движение), (на стр 28)
	Algorithm: V. I. Krylov, Fundamentals of the theory of motion of the artificial earth satellites (part one: undisturbed motion). (page 28)
	*/
	double per_arg = orb.periapsis_argument * PI / 180;
	double T = orb.true_anomaly * PI / 180;
	double i = orb.inclination * PI / 180;
	double L = orb.ascending_long * PI / 180;

	double u = T + per_arg; // latitude argument
	double p = orb.orb_major_ax * (1 - orb.orb_eccentricity * orb.orb_eccentricity); // fokal parameter
	double r = p / (1 + orb.orb_eccentricity * cos(T)); // radius vector

	double x = r * (cos(u) * cos(L) - sin(u) * sin(L) * cos(i));
	double y = r * (sin(L) * cos(u) + sin(u) * cos(L) * cos(i));
	double z = r * (sin(u) * sin(i));

	double Vr = sqrt(GE / p) * orb.orb_eccentricity * sin(T); // radial velocity
	double V1 = sqrt(GE / p) * (1 + orb.orb_eccentricity * cos(T));

	double Vx = Vr * (cos(u) * cos(L) - sin(u) * sin(L) * cos(i)) + V1 * (-cos(L) * sin(u) - sin(L) * cos(u) * cos(i));
	double Vy = Vr * (sin(L) * cos(u) + sin(u) * cos(L) * cos(i)) + V1 * (-sin(L) * sin(u) + cos(L) * cos(u) * cos(i));
	double Vz = Vr * (sin(u) * sin(i)) + V1 * (cos(u) * sin(i));

	return { {x, y, z}, {Vx, Vy, Vz} };
}

My_vector Get_a(const My_vector& velocity, const satellite& sat, const atmosphere& atmo) {
	double a1 = (sat.Cx * atmo.ro * sat.S) / (2. * sat.mass);
	My_vector vel = (velocity.Length() * velocity.Length()) * velocity.Normalize();
	My_vector a = vel * a1 * -1;
	return a;
}

My_vector Get_thrust_a(const My_vector& coord, const My_vector& velocity, const double& thrust, const int& way) {
	My_vector acc;
	if (way == 1) { // против движения
		acc = velocity.Normalize();
		acc.Reverse();
		return acc * thrust;
	}
	else if (way == 2) { // по движению
		acc = velocity.Normalize();
		return acc * thrust;
	}
	else if (way == 3) { // направо от движения (к центру)
		acc = coord.Normalize();
		acc.Reverse();
		return acc * thrust;
	}
	else if (way == 4) { // налево от движения (от центра)
		acc = coord.Normalize();
		return acc * thrust;
	}
	else if (way == 5) { // вверх

	}
	else if (way == 6) { // вниз

	}

	return {};
}

	My_vector Get_planets_a(double current_date, std::vector<planet_data>& v, const My_vector& coord) {
		
		return {};
	}

std::pair<My_vector, My_vector> Function(const std::pair<My_vector, My_vector>& state, double dt, const satellite& sat, const atmosphere& atmo, 
	const My_vector& atm_vel, const planets& planets_struct) {
	My_vector coord = state.first;
	My_vector velocity = state.second; // orbital velocity
	My_vector atmo_velocity = velocity - atm_vel;
	My_vector g = Get_g(coord); // acceleration of free fall
	My_vector a = Get_a(atmo_velocity, sat, atmo); // acceleration of drag of atmosphere
	//My_vector a = {0.,0.,0.}; // acceleration of drag of without atmosphere

	if (sat.thrust_on == true) {
		a = a + Get_thrust_a(coord, velocity, (sat.thrust / sat.mass), sat.thrust_way); // plus thrust acceleration
	}

	if (planets_struct.general_mark == true) { // plus planets influence
		My_vector local_a = { 0., 0., 0. };
		if (planets_struct.Moon.first == true) {

		}
		if (planets_struct.Sun.first == true) {

		}
		a = a + local_a;
	}

	My_vector new_coord = velocity;
	My_vector new_velocity = g + a;

	return { new_coord, new_velocity };
}

bool Check_turns(const double& start, const double& angle, const char& direction, const bool half_turn) {
	if (direction == '+') {
		if (half_turn == false) {
			if (angle > start) return false;
			else return true;
		}
		else {
			if (angle > start) return true;
			else return false;
		}
	}
	else {
		if (half_turn == false) {
			if (angle < start) return false;
			else return true;
		}
		else {
			if (angle > start) return true;
			else return false;
		}
	}
}

std::vector<double> Calculate_keplers_elements(const My_vector& coord_, const My_vector& velocity_) {
	/*
	 * Справочное руководство по небесной механике и астродинамике. Абалакин В. К., Аксенов Е. П., Гребеников Е. А., Демин В. Г., Рябов Ю. А.
	 * Издание 2-е, дополненное и переработанное. (стр 283)
	 * Algorithm: Reference Guide to Celestial Mechanics and Astrodynamics. Abalakin V. K., Aksenov E. P., Grebenikov E. A., Demin V. G., Ryabov Yu. A.
	 * 2nd edition, expanded and revised. (page 283)
	 */

	My_vector vec_mod = coord_ / velocity_;
	My_vector k = (vec_mod) / vec_mod.Length();
	double i = acos(k.GetZ()) * 180. / PI;
	double r = coord_.Length();
	double a = 1 / (2 / r - velocity_.Length() * velocity_.Length() / GE);

	return { i, a };
}

double Calculate_line_velocity_rotation(const double H_, const double latitude_) {
	double angle = latitude_ * PI / 180.;
	double hight = 0.;
	double v1 = (major_ax * minor_ax) / sqrt(minor_ax * minor_ax + major_ax * major_ax * tan(angle) * tan(angle));
	double v2 = (minor_ax * minor_ax * hight) / sqrt(pow(minor_ax, 4) + pow(major_ax, 4) * tan(angle) * tan(angle));
	double v = (v1 + v2) * Angular_speed_Earth;
	return v;
}

My_vector Calculate_vector_velocity_rotation(const My_vector& coord_) {
	const My_vector omega{ 0, 0, Angular_speed_Earth };
	My_vector v = omega / coord_;
	//	v.Reverse();
	return v;
}

