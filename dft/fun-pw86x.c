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
/* fun-pw86x.c:
   The PW86 exchange functional and its derivative.
   (c) Olav Fossgaard, olav@chem.uit.no, may 2002
 
   Reference: Phys. Rev. B 33. 8800 (1986)
*/

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
static int pw86x_isgga(void) { return 1; }
static int pw86x_read(const char* conf_line);
static real pw86x_energy(DftDensProp* dp);
static void pw86x_first (FirstFuncDrv *ds,      real factor, DftDensProp* dp);
static void pw86x_second(SecondFuncDrv *ds, real factor, DftDensProp* dp);
static void pw86x_third (ThirdFuncDrv *ds,  real factor, DftDensProp* dp);

Functional PW86xFunctional = {
  "PW86x",       /* name */
  pw86x_isgga,   /* gga-corrected */
  pw86x_read, 
  NULL,
  pw86x_energy, 
  pw86x_first,
  NULL,
  NULL
};

/* IMPLEMENTATION PART */
static int
pw86x_read(const char* conf_line)
{
  dft_set_hf_weight(0);
  return 1;
}

static real
pw86x_energy(DftDensProp* dp)
{
/* Use density functional form. In case of spin polarization,
   this function will have to be called twice with arguments
   rho=2rhoa and rho=2rhob, respectively. The total energy is then
   half the sum of the returned values.
*/
  const real a = 1.0;
  const real b = 1.296;
  const real c = 14.0;
  const real d = 0.20;
/* Closed shell (See eq. (25) in reference) */
  real rho = dp->rhoa+dp->rhob, grad = dp->grada+dp->gradb;

  const real Ax = -pow(3.0/M_PI,1.0/3.0)*3.0/4.0;
  const real kf = pow(3.0*pow(M_PI,2.0)*rho,1.0/3.0);
  real s = grad/(2.0*kf*rho);
  real F = pow(a+b*pow(s,2.0)+c*pow(s,4.0)+d*pow(s,6.0),1.0/15.0); 
  return Ax*pow(rho,4.0/3.0)*F;
}

static void
pw86x_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
{
/* The energy expression is the integral int(Ax*rho**(4/3)*F). We first 
   calculate d(F)/d(rho) and d(F)/d(grad_rho) and differentiate the
   product in the last step.
*/
  const real a = 1.0;
  const real b = 1.296;
  const real c = 14.0;
  const real d = 0.20;
/* Closed shell (See eq. (25) in reference) */
  real rho = dp->rhoa+dp->rhob, grad = dp->grada+dp->gradb;

  const real Ax= -pow(3.0/M_PI,1.0/3.0)*3.0/4.0;
  const real kf= pow(3.0*M_PI*M_PI*rho,1.0/3.0);
  real  s = grad/(2.0*kf*rho);
  real  F = pow(a+b*pow(s,2.0)+c*pow(s,4.0)+d*pow(s,6.0),1.0/15.0);

  real F1 = 1.0/15.0*pow(a+b*pow(s,2.0)+c*pow(s,4.0)+d*pow(s,6.0),-14.0/15.0)
            *(2.0*b*s+4.0*c*pow(s,3.0)+6.0*d*pow(s,5.0)); /* dF/ds */

  real s1 = -4.0*s/(3.0*rho); /* ds/d(rho) */
  real s2 = 1.0/(2.0*kf*rho); /* d(s)/d(grad) */
  real G1 = F1*s1; /* dF/d(rho) */ 
  real G2 = F1*s2; /* dF/d(grad) */
   
  ds->df1000 += Ax*((4.0/3.0)*pow(rho,1.0/3.0)*F + pow(rho,4.0/3.0)*G1 )*factor;
  ds->df0010 += Ax*pow(rho,4.0/3.0)*G2*factor;
}
