#pragma once
#include <cmath>

const double PI = 3.14159265;

//--------CONSTANTS_OF_EARTH---------
const double major_ax = 6'378'136.5; //Semimajor axis, meters
const double minor_ax = 6'356'751.758; //Semiminor axis, meters
const double eccentricity = sqrt(major_ax * major_ax - minor_ax * minor_ax) / major_ax;
const double eccentricity_sec = sqrt(major_ax * major_ax - minor_ax * minor_ax) / minor_ax;
//const double eccentricity = 0.;
//const double minor_ax = 6'378'136.5; //Semiminor axis, meters
//const double eccentricity_sec = 0.;
const double flattening = 1 / 298.2564151;
const double gravitational_const = 6.674184e-11; // gravitational constant, m^3 * kg^(-1) * s^(-2)
const double Earth_mass = 5.9722e24; // kg
//const double GE =  gravitational_const*Earth_mass; // geocentr gravitational constant, m^3 * s^(-2), from DE200
const double GE = 3.986004391e14; // geocentr gravitational constant, m^3 * s^(-2), from DE200
const double Angular_speed_Earth = 7.2921158553e-5; // second^-1
const double Period = 86164.090530833; // seconds
//---------------------------------