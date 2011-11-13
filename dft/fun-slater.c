/*
       Copyright (c) 2004 by the authors of Dirac (see below).
       All Rights Reserved.
    
       The source code in this file is part of
       Dirac, a relativistic ab initio electronic structure program,
       release DIRAC04.0 (2004),
       written by H. J. Aa. Jensen, T. Saue, and L. Visscher
       with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
       O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
       P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
       P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
       (http://dirac.chem.sdu.dk).
    
       For a suitable BibTEX entry, see:
          http://dirac.chem.sdu.dk/doc/reference.shtml
    
       This source code is provided under a written licence and may be
       used, copied, transmitted, or stored only in accord with that
       written licence.
    
       In particular, no part of the source code or compiled modules may
       be distributed outside the research group of the licence holder.
       This means also that persons (e.g. post-docs) leaving the research
       group of the licence holder may not take any part of Dirac,
       including modified files, with him/her, unless that person has
       obtained his/her own licence.
    
       For questions concerning this copyright write to:
          dirac-admin@dirac.chem.sdu.dk
    
       For information on how to get a licence see:
          http://dirac.chem.sdu.dk
*/
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
static void dirac_first (FirstFuncDrv *ds, real factor,      DftDensProp* dp);
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
dirac_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
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
