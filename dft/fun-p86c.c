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
/* fun-p86x.c:
   The P86 correlation functional and its derivative. This is 
   a correction term to the LDA correlation.
   (c) Olav Fossgaard, olav@chem.uit.no, may 2002

   Reference: Phys. Rev. B 33, 8822 (1986)
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
static int p86c_isgga(void) { return 1; }
static int p86c_read(const char* conf_line);
static real p86c_energy(DftDensProp* dp);
static void p86c_first (FirstFuncDrv *ds,      real factor, DftDensProp* dp);
static void p86c_second(SecondFuncDrv *ds, real factor, DftDensProp* dp);
static void p86c_third (ThirdFuncDrv *ds,  real factor, DftDensProp* dp);

Functional P86cFunctional = {
  "P86c",       /* name */
  p86c_isgga,   /* gga-corrected */
  p86c_read, 
  NULL,
  p86c_energy, 
  p86c_first,
  NULL,
  NULL
};

/* IMPLEMENTATION PART */
static int
p86c_read(const char* conf_line)
{
 dft_set_hf_weight(0);
 return 1;
}

static const real   c0 = 1.667e-3;
static const real  cn1 = 2.568e-3;
static const real  cn2 = 2.3266e-2;
static const real  cn3 = 7.389e-6;
static const real  cd1 = 1.0;
static const real  cd2 = 8.723;
static const real  cd3 = 0.472;

static real
p86c_energy(DftDensProp* dp)
{
  const real  cd4 = 1.0e4*cn3;
  const real cinf = c0+cn1/cd1;
  const real phi1 = pow(9.0*M_PI,1.0/6.0); 
  const real phi2 = 0.11;

  real rhoa = dp->rhoa, rhob = dp->rhob; 
  real grada = dp->grada, gradb = dp->gradb; 
/* Form the spin polarization zeta from rhoa and rhob. */
  real zeta = (rhoa-rhob)/(rhoa+rhob);
  real dzeta = pow(2.0,1.0/3.0)
              *sqrt( pow((1+zeta)/2.0,5.0/3.0) + pow((1-zeta)/2.0,5.0/3.0) );

/* We use total density rho and spin polarity zeta, so we write */
  real   rho = rhoa+rhob, grad = grada+gradb;

  real    rs = pow(3.0/(4.0*M_PI),1.0/3.0)*pow(rho,-1.0/3.0);  
  real crhon = cn1+cn2*rs+cn3*pow(rs,2.0);
  real crhod = cd1+cd2*rs+cd3*pow(rs,2.0)+cd4*pow(rs,3.0);
  real  crho = c0+crhon/crhod;
  real   phi = phi1*phi2*(cinf/crho)*grad*pow(rho,-7.0/6.0);

/* Get the LDA energy and add the correction term */  
  real  vwne = VWNFunctional.func(dp);
  
  return vwne + pow(dzeta,-1.0)*exp(-phi)*crho*pow(grad,2.0)*pow(rho,-4.0/3.0);
}

static void
p86c_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
{
  const real  cd4 = 1.0e4*cn3;
  const real cinf = c0+cn1/cd1;
  const real phi1 = pow(9.0*M_PI,1.0/6.0); 
  const real phi2 = 0.11;
  real rhoa = dp->rhoa, rhob = dp->rhob; 
  real grada = dp->grada, gradb = dp->gradb; 

/* Form the spin polarization zeta from rhoa and rhob. */
  real  zeta = (rhoa-rhob)/(rhoa+rhob);
  real dzeta = pow(2.0,1.0/3.0)
  *sqrt( pow((1+zeta)/2.0,5.0/3.0) + pow((1-zeta)/2.0,5.0/3.0) );

/* We use total density rho and spin polarity zeta, so we write */
  real   rho = rhoa+rhob, grad = grada+gradb;

  real    rs = pow(3.0/(4.0*M_PI),1.0/3.0)*pow(rho,-1.0/3.0);  
  real crhon = cn1+cn2*rs+cn3*pow(rs,2.0);
  real crhod = cd1+cd2*rs+cd3*pow(rs,2.0)+cd4*pow(rs,3.0);
  real  crho = c0+crhon/crhod;
  real   phi = phi1*phi2*(cinf/crho)*grad*pow(rho,-7.0/6.0);

/* rho derivative */
/* Derive C(rho) and phi */
  real drs = (-1.0/3.0)*pow(3.0/(4.0*M_PI),1.0/3.0)*pow(rho,-4.0/3.0); /* d(rs)/d(rho) */  

  real dcrhon = cn2 + 2.0*cn3*rs;
  real dcrhod = cd2 + 2.0*cd3*rs + 3.0*cd4*pow(rs,2.0);
  real dcrho  = ( (dcrhon*crhod - crhon*dcrhod)/pow(crhod,2.0) )*drs; /* dC(rho)/d(rho) */
  
  real  dphir = -phi1*phi2*cinf*grad
               *(pow(rho,7.0/6.0)*dcrho + (7.0/6.0)*pow(rho,1.0/6.0)*crho)
               /pow(pow(rho,7.0/6.0)*crho,2.0); /* d(phi)/d(rho) */
  
/* The derivative is a constant times the derivative of the product 
   abc = exp(-phi)*crho*pow(rho,-4.0/3.0)  */
  real dprod3 = (-exp(-phi)*dphir*crho+exp(-phi)*dcrho)*pow(rho,-4.0/3.0) /* (ab)'c */
               +exp(-phi)*crho*(-4.0/3.0)*pow(rho,-7.0/3.0); /* (ab)c' */
  
/* grad derivative */
/* A constant times the derivative of the product exp(-phi)*pow(grad,2.0) */
  real dphig = phi1*phi2*(cinf/crho)*pow(rho,-7.0/6.0); /* d(phi)/d(grad) */
  real dprod2 =  -exp(-phi)*dphig*pow(grad,2.0) + exp(-phi)*2.0*grad;

/* Get the LDA derivatives and add the correction */  
  VWNFunctional.first(ds, factor, dp);
  
  ds->df1000 += factor*pow(dzeta,-1.0)*pow(grad,2.0)*dprod3;
  ds->df0010 += factor*pow(dzeta,-1.0)*crho*pow(rho,-4.0/3.0)*dprod2;
}
