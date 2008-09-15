/*
!dirac_copyright_start
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
!dirac_copyright_end
*/
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-Slater.c:
   implementation of Slater functional and its derivatives 
   (c), Pawel Salek, pawsa@theochem.kth.se, aug 2001
*/

/* strictly conform to XOPEN ANSI C standard */
#if !defined(SYS_DEC)
/* XOPEN compliance is missing on old Tru64 4.0E Alphas and pow() prototype
 * is not specified. */
#define _XOPEN_SOURCE          500
#define _XOPEN_SOURCE_EXTENDED 1
#endif
#include <math.h>
#include <stdio.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int slater_isgga(void) { return 0; }
static int slater_read(const char* conf_line);
static real slater_energy(const FunDensProp* dp);
static void slater_first(FunFirstFuncDrv *ds,   real factor, const FunDensProp* dp);
static void slater_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp);
static void slater_third(FunThirdFuncDrv *ds,   real factor, const FunDensProp* dp);

Functional SlaterFunctional = {
  "Slater",       /* name */
  slater_isgga,   /* gga-corrected */
  slater_read, 
  NULL,
  slater_energy, 
  slater_first,
  slater_second,
  slater_third
};

/* IMPLEMENTATION PART */
static int
slater_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

/* SLATER_THRESHOLD Only to avoid numerical problems due to raising 0
 * to a fractional power. */
static const real SLATER_THRESHOLD = 1e-20;
static real
slater_energy(const FunDensProp* dp)
{
  real ea = 0.0, eb = 0.0; 
  const real PREF= -3.0/4.0*pow(6/M_PI, 1.0/3.0);
  if (dp->rhoa >SLATER_THRESHOLD)
      ea= PREF*pow(dp->rhoa,4.0/3.0);
  if (dp->rhob >SLATER_THRESHOLD)
      eb= PREF*pow(dp->rhob,4.0/3.0);   
  return ea+eb;  
}

static void
slater_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  if (dp->rhoa>SLATER_THRESHOLD)
     ds->df1000 += -pow(6.0/M_PI*dp->rhoa, 1.0/3.0)*factor;
  if (dp->rhob>SLATER_THRESHOLD)
     ds->df0100 += -pow(6.0/M_PI*dp->rhob, 1.0/3.0)*factor;
}
static void
slater_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>SLATER_THRESHOLD) {
    ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
  }
  if (dp->rhob>SLATER_THRESHOLD) {
    ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
  }
}

/* slater_third:
   Slater functional derivatives.
*/
static void
slater_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  const real PREF = pow(6.0/M_PI, 1.0/3.0);
  if (dp->rhoa>SLATER_THRESHOLD) {
    ds->df1000 += -PREF*pow(dp->rhoa,  1.0/3.0)*factor;
    ds->df2000 += -PREF*pow(dp->rhoa, -2.0/3.0)/3*factor;
    ds->df3000 +=  PREF*pow(dp->rhoa, -5.0/3.0)*2.0/9.0*factor;
  }
  if (dp->rhob>SLATER_THRESHOLD) {
    ds->df0100 += -PREF*pow(dp->rhob,  1.0/3.0)*factor;
    ds->df0200 += -PREF*pow(dp->rhob, -2.0/3.0)/3*factor;
    ds->df0300 +=  PREF*pow(dp->rhob, -5.0/3.0)*2.0/9.0*factor;
  }
}
