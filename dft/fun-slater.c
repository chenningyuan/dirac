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
#include "dcbham.h"

/* INTERFACE PART */
static int dirac_read(const char* conf_line);
static real dirac_energy(real rhoa, real rhob, real grada, real gradb);
static void dirac_first(FirstDrv *ds, real factor, 
                        real rhoa, real rhob, real ngrada, real ngradb);
static void dirac_second(SecondFuncDrv *ds, real factor,
                         real rhoa, real rhob, real ngrada, real ngradb);
static void dirac_third(ThirdFuncDrv *ds, real factor,
                        real rhoa, real rhob, real ngrada, real ngradb);

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
#ifndef NO_BACKWARD_COMP
  dcrham_.hfxfac = 0;
#endif
 return 1;
}

static real
dirac_energy(real rhoa, real rhob, real grada, real gradb)
{
  const real PREF= -3.0/4.0*pow(6/M_PI, 1.0/3.0);
  return PREF*(pow(rhoa,4.0/3.0) + pow(rhob,4.0/3.0));
}

static void
dirac_first(FirstDrv *ds, real factor, 
            real rhoa, real rhob, real ngrada, real ngradb)
{
  ds->df1000 += -pow(6.0/M_PI*rhoa, 1.0/3.0)*factor;
}
static void
dirac_second(SecondFuncDrv *ds, real factor,
             real rhoa, real rhob, real ngrada, real ngradb)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  ds->df1000 += -PREF*pow(rhoa,  1.0/3.0)*factor;
  ds->df2000 += -PREF*pow(rhoa, -2.0/3.0)/3*factor;
}

/* dirac_third:
   Slater functional derivatives.
*/
static void
dirac_third(ThirdFuncDrv *ds, real factor,
            real rhoa, real rhob, real ngrada, real ngradb)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  ds->df1000 += -PREF*pow(rhoa,  1.0/3.0)*factor;
  ds->df2000 += -PREF*pow(rhoa, -2.0/3.0)/3*factor;
  ds->df3000 +=  PREF*pow(rhoa, -5.0/3.0)*2.0/9.0*factor;
}
