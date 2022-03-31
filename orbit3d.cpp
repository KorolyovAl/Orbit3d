#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <cmath>

#include "constants.h"
#include "structures.h"
#include "myvector_operators.h"
#include "msis.h"
#include "reader.h"
#include "aerodynamic.h"
#include "functions.h"
#include "planets.h"

using namespace std::chrono;

int main()
{
	elements_of_orbit orb;
	satellite sat;
	atmosphere atmo;
	planets planet_marks{
		{true},							// general mark
		{true, "Moon_geocenter.txt"},	// Moon
		{true, "Sun_Earth.txt"}			// Sun
	};
	Planets_class objPlanets(planet_marks);
	reader(orb, sat);
	Aerodynamic_check_case();

	if (planet_marks.general_mark == true) {
		if (planet_marks.Moon.first == true) {
			objPlanets.Fill_Moon(planet_marks.Moon.second);
		}
		if (planet_marks.Sun.first == true) {
			objPlanets.Fill_Sun(planet_marks.Sun.second);
		}
	}
	
	sat.S = sat.d * sat.d * 3; // pi*d*d/4 for sphere, d*d for cubesat
	sat.thrust_on = false;

	ios_base::sync_with_stdio(false);
	My_vector periapsis_coord = Calculate_coordinates_periapsis(orb.inclination, orb.periapsis_argument, orb.periapsis, eccentricity);
	My_vector apoapsis_coord = Calculate_coordinates_periapsis(orb.inclination, orb.periapsis_argument - 180., orb.apoapsis, eccentricity);
	orb.orb_eccentricity = (apoapsis_coord.Length() - periapsis_coord.Length())
		/ (apoapsis_coord.Length() + periapsis_coord.Length()); // (a-p)/(a+p)

	auto start = steady_clock::now();

	orb.orb_major_ax = periapsis_coord.Length() / (1 - orb.orb_eccentricity);
	orb.orb_minor_ax = orb.orb_major_ax * sqrt(1 - orb.orb_eccentricity * orb.orb_eccentricity);
	pair<My_vector, My_vector> state = Calculate_state_vector(orb); // vector of state (coordinates and velocity)

	My_vector atm_vel; // velocity of atmosphere
	My_vector velocity;
	My_vector coordinates;

	double time(0.), timestep(1.), step(0.);
	pair<My_vector, My_vector> U, Uold;
	U = Uold = state;
	coordinates = U.first;
	velocity = U.second;
	pair<My_vector, My_vector> k1, k2, k3, k4;
	ofstream output("results/test_with_optimization.txt");
	//ofstream output ("incl_65_1-1-1-0-400.txt");

	__int64 day_count = 0; // count of days
	int day = 1; // initial day of year
	//int day_hour = 19, day_minute = 28, day_second = 34; // time of start PS-1, UTC
	int day_hour = 1, day_minute = 0, day_second = 0;
	int day_time = day_second + day_minute * 60 + day_hour * 3600; // time of current day in seconds
	int hours = 0; // count of hours
	bool half_turn = false;
	double turns = 0.;
	char direction = '+'; // direction of flight, '+' - from greenwich to 90 degree of east longitude and etc. '-' - from greenwich to 90 degree of west longitude and etc
	double initial_angle = 0.;
	double B, L, H; // B - latitude, L - longitude, H - height
	double B0, L0, H0; // initial: B - latitude, L - longitude, H - height
	double orb_period = 2 * PI * sqrt(orb.orb_major_ax * orb.orb_major_ax * orb.orb_major_ax / GE) / 60.; // orbital period, minutes

	vector<double> Elipt = Calculate_elliptic_coordinates(coordinates);
	L = L0 = Elipt[0]; B = B0 = Elipt[1]; H = H0 = Elipt[2] / 1000.;


	cout << "time: " << time << " hours: " << hours << " period: " << orb_period << endl;
	cout << "L: " << Elipt[0] << " B: " << Elipt[1] << " H : " << Elipt[2] / 1000. << endl;
	output << "day: " << day << "; F10.7: " << 200 << "; mass: " << sat.mass << "; d: " << sat.d << "; H: " << H
		<< "; apoapsis: " << orb.apoapsis << "; periapsis: " << orb.periapsis << "; initial time (UTC): hour = " << day_hour << " minute = " << day_minute << endl;

	while (H > 10.) {

		atmosMSIS00(
			H,    // altitude [km]
			B,   // latitude [deg]
			L,  // longitude [deg]
			90.,   // F10.7 solar activity index
			day,    // day number sins 1,jan
			day_time, // seconds in day (UT)
			&atmo.d, // number density [1/m^3]
			&atmo.ro,  // density [kg/m^3]
			&atmo.tk,  // Temperature [K]
			&atmo.pres,  // pressure [Pa]
			&atmo.mass,  // molecular mass [ce]
			&atmo.gamma, // cp/cv ratio
			&atmo.sv, // sound velocity [m/c]
			&atmo.lambda // mean free path [m]
		);

		if (H <= 400.) {
			timestep = 0.1;
		}

		/*		if (time >= 662.65 * 60. * 60) {
					sat.thrust_on = false;
				}
		*/
		sat.Cx = Get_Cx(atmo, U);

		atm_vel = Calculate_vector_velocity_rotation(coordinates);

		k1 = Function(U, 0., sat, atmo, atm_vel) * timestep;
		k2 = Function(U + 0.5 * k1, timestep / 2., sat, atmo, atm_vel) * timestep;
		k3 = Function(U + 0.5 * k2, timestep / 2., sat, atmo, atm_vel) * timestep;
		k4 = Function(U + k3, timestep, sat, atmo, atm_vel) * timestep;

		U = U + 1. / 6. * (k1 + 2 * k2 + 2 * k3 + k4);
		coordinates = U.first;
		velocity = U.second;

		sat.thrust_on = false;

		My_vector a = Get_a(velocity, sat, atmo);
		vector<double> Elipt = Calculate_elliptic_coordinates(coordinates);
		L = Elipt[0]; B = Elipt[1]; H = Elipt[2] / 1000.;

		int delta_t = 60; // step in time for print
		if (step == 1 * 60. / timestep) { // 86'400'000 - one day with step = 0.001 |||| print every 60 seconds
			day_count = time / 86'400;
			day_time += delta_t;
			if (day_time >= 86'400) {
				day_time = 0;
				day++;
			}
			if (day > 365) {
				day = 1;
			}
			hours = (time - day_count * 86'400) / 3'600;

			//------------------------------------------BLOK_OF_COUNT_OF_TURNS-----------------------------------------------------
			if (orb.inclination < 80.) {
				initial_angle = L0;
			}
			else {
				initial_angle = B0;
			}
			if (Check_turns(initial_angle, Elipt[0], direction, half_turn)) {
				turns += 0.5;
				if (turns - int(turns) != 0.) {
					half_turn = true;
				}
				else {
					half_turn = false;
				}
			}
			//------------------------------------------BLOK_OF_COUNT_OF_TURNS-----------------------------------------------------

			vector<double> Keplers = Calculate_keplers_elements(coordinates, velocity);
			orb_period = 2 * PI * sqrt(Keplers[1] * Keplers[1] * Keplers[1] / GE) / 60.;

			std::cout << setfill(' ') << fixed << setprecision(0);
			std::cout << "day: " << setw(3) << day_count << " hours: " << setw(2) << hours << " year day: " << day << std::endl;
			std::cout << setprecision(2) << "L: " << setw(7) << L << " B: " << setw(7) << B << " H : " << setw(7) << H;
			std::cout << setprecision(4) << " i = " << Keplers[0] << scientific << " acceleration = " << setw(7) << a.GetX() << " " 
				<< setw(7) << a.GetY() << " dens = " << atmo.ro
				<< fixed << setprecision(2) << " period: " << orb_period << " Cx: " << sat.Cx << std::endl;

			output << setfill(' ') << fixed << setprecision(0);
			output << "time: " << time << " day: " << setw(3) << day_count << " hours: " << setw(2) << hours
				<< setprecision(2) << " L: " << setw(7) << Elipt[0] << " B: " << setw(7) << Elipt[1] << " H: " << setw(7) << Elipt[2] / 1000. << " N: " << turns
				//		<< " N: " << turns << " inclination: " << Keplers[0] << " dens = " << atmo.ro << " temperature = " << atmo.tk
				<< setprecision(1) << " velocity: " << setw(7) << velocity.GetX() << " " << setw(7) << velocity.GetY()
				//		<< setprecision(1) << " position: x = " << setw(7) << U.first.x << " y =  " << setw(7) << U.first.y << " z = " << setw(7) << U.first.z
				<< setprecision(2) << " acceleration: " << scientific << setw(7) << a.GetX() << " " << setw(7) << a.GetX()
				<< fixed << setprecision(2) << " period: " << setw(4) << orb_period
				<< setprecision(4) << " inclination: " << setw(4) << Keplers[0]
				<< setprecision(2) << " Cx: " << setw(4) << sat.Cx
				<< scientific << " dens = " << atmo.ro << '\n';
			velocity.Print();
			step = 0;
		}

		step += 1;
		time += timestep;
		Uold = U;

	}

	cout << fixed << "common time: " << time << endl;
	output << fixed << "common time : " << time;

	auto end = steady_clock::now();
	auto elapsed = duration_cast<milliseconds>(end - start);
	cerr << elapsed.count() << " ms";

	return 0;
}