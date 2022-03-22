/* msis.cc
 * wrapper to msis 90 standart atmospfere
 * 15.03.13 A.V. Kashkovsky
 */

// http://www.brodo.de/space/nrlmsise/index.html
extern "C" {
#include "nrlmsise-00.h"
};
#define _USE_MATH_DEFINES
#include "msis.h"


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
 )
{
 struct nrlmsise_output output;
 struct nrlmsise_input input;
 struct nrlmsise_flags flags;
 struct ap_array aph;
 int i;
 int j;
 /* input values */
 for (i=0;i<7;i++) aph.a[i]=100;
 for(i=0;i<24;i++) flags.switches[i]=1;

 input.year=0;
 input.alt = alt;
 input.f107 = f107;
 input.g_lat = glat;
 input.g_long = glong;
 input.doy = day; 
 input.f107A = input.f107;

 input.sec = sec;
 input.lst = input.g_long/15.0 + input.sec/3600.0;
 if(input.lst < 0.0) input.lst += 24.0;
 if (input.lst > 24)
 {
	 input.lst -= 24.0;
 }
 input.ap=4;

// call of atmosphere
 if (alt > 500.) {
	 gtd7d(&input, &flags, &output);
 }
 else {
	  gtd7(&input, &flags, &output);
 }

  float mmass[]={
     4.0026,  // He
    15.9994,  // O
    28.01348, // N2
    31.9988,  // O2
    39.948,   // Ar
     1.00794, // H
    14.00674  // N
    };
  double diam[]={
    2.33e-10,  // He
    2.96E-10,  // O
    4.11e-10,  // N2
    4.01e-10,  // O2
    3.595e-10, // Ar
    2.33e-10,  // H
    2.96E-10   // N
    };
  double alpha[]={
    0.16, // He
    0.25, // O
    0.24, // N2
    0.27, // O2
    0.24, // Ar
    0.25, // H
    0.25  // N
    };

// d[5] - mass density, all other - number density of component
  double nd[7];
  *dens = output.d[5];
  for( i=0;i<5;i++) nd[i] = output.d[i]; 
  nd[5] = output.d[6];
  nd[6] = output.d[7];
  *temp = output.t[1];

  double Boltz;
  Boltz = 1.3806e-23;
  double dof[7];

// mean molecular mass
  *ndens=0.;
  *mass=0.;
  for( i=0;i<7;i++)
  {
   *ndens+=nd[i];
   *mass +=  mmass[i]*nd[i];
  }
  *mass /= *ndens;

/* mean free path calcculated by VHS parameters
Bird94,p. 96

$$\lambda=\sum_{s=1}^{N_{sp}}{n_s \over n}\Big[
\sum_{r=1}^{N_{sp}}\pi (d_{ref}^{rs})^2 n_r \Big({T_{ref}^{rs}\over T_t}
\Big)^{\omega_{rs}-0.5}\sqrt{1+{m_s \over m_r}}\Big]^{-1}.$$

\noindent where $T_{ref}^{rs}=T_{ref}^{r}+T_{ref}^{s}$ and
$d_{ref}^{rs}=d_{ref}^{r}+d_{ref}^{s}$,
and $\omega^{rs}=(\omega^{r}+\omega^{s})/2$ is the temperature-viscosity exponent
(if no other value is prescribed).

$$ \alpha = \omega-0.5$$
*/
  double s, dref, Tref;
  *lambda=0;
  Tref = 273./(*temp);
  for( i=0;i<7;i++)
  { s = 0;
    for(j=0;j<7;j++)
    { dref = diam[i]+diam[j];
      s+= M_PI * dref * dref *nd[j] *pow( Tref, alpha[j] ) *sqrt( 1.+mmass[i]/mmass[j]);
    }
    (*lambda)+=nd[i]/s;
  }
  (*lambda)/=*ndens;
// pressure
  *pres = *ndens*Boltz*(*temp);

// calculate degree of freedom and gamma
  for( i=0;i<7;i++) dof[i]=3.;
  dof[2] += 2. + 2.*3371./(*temp)/(exp(3371./(*temp))-1.); //N2
  dof[3] += 2. + 2.*2256./(*temp)/(exp(2256./(*temp))-1.); //O2
  *gamma=0.;
  for( i=0;i<7;i++) *gamma+=nd[i]/(*ndens)*(dof[i]+2.)/dof[i];
// sound velosity
  *sound = sqrt( *gamma * 8314.510/(*mass) * (*temp) );

}
