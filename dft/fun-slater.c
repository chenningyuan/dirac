/* fun-Slater.c:
   implementation of Slater functional and its derivatives 
   (c), Pawel Salek, pawsa@theochem.kth.se, aug 2001
   NOTE:
   this file may seem unnecessarily complex but the structure really pays off
   when implementing multiple functionals depending on different parameters.
*/

/* strictly conform to XOPEN ANSI C standard */

#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int dirac_read(const char* conf_line);
static real dirac_energy(DftDensProp* dp);
static void dirac_first (FirstDrv *ds, real factor,      DftDensProp* dp);
static void dirac_second(SecondFuncDrv *ds, real factor, DftDensProp* dp);
static void dirac_third (ThirdFuncDrv *ds, real factor,  DftDensProp* dp);

Functional SlaterFunctional = {
  "Slater",       /* name */
  fun_false,   /* not gga-corrected */
  dirac_read, 
  NULL,
  dirac_energy, 
  dirac_first,
  dirac_second,
  dirac_third
};

/* IMPLEMENTATION PART */
static int
dirac_read(const char* conf_line)
{
  dft_set_hf_weight(0.0);
  return 1;
}

static real
dirac_energy(DftDensProp* dp)
{
  const real PREF= -3.0/4.0*pow(6/M_PI, 1.0/3.0);
  return PREF*(pow(dp->rhoa,4.0/3.0) + pow(dp->rhob,4.0/3.0));
}

static void
dirac_first(FirstDrv *ds, real factor, DftDensProp* dp)
{
  ds->df1000 += -pow(6.0/M_PI*dp->rhoa, 1.0/3.0)*factor;
  ds->df0100 += -pow(6.0/M_PI*dp->rhob, 1.0/3.0)*factor;
}
static void
dirac_second(SecondFuncDrv *ds, real factor, DftDensProp* dp)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
  ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
  ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
  ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
}

/* dirac_third:
   Slater functional derivatives.
*/
static void
dirac_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
  ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
  ds->df3000 +=  PREF*pow(dp->rhoa, -5.0/3.0)*2.0/9.0*factor;
  ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
  ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
  ds->df0300 +=  PREF*pow(dp->rhob, -5.0/3.0)*2.0/9.0*factor;
}
