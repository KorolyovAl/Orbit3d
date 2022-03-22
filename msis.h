#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cmath>

using namespace std;

void atmosMSIS00(
 float alt,    // altitude [km]
 float glat,   // latitude [deg]
 float glong,  // longitude [deg]
 float f107,   // F10.7 solar activity index
 int   day,    // day number sins 1,jan
 int   sec,    // seconds in day (UT)
 double *ndens, // number density [1/m^3]
 double *dens,  // density [kg/m^3]
 double *temp,  // temperture [K]
 double *pres,  // pressure [Pa]
 double *mass,  // molecular mass [ce]
 double *gamma, // cp/cv ratio
 double *sound, // sound velosity [m/c]
 double *lambda // mean free path [m]
 );
