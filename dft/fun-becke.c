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
/* fun-becke.c:
   implementation of Becke(88) functional and its derivatives.
   or exactly: Becke GGA correction to the functional. (total Becke(88)
   energy is E_LDA+E_BCK).
   (c) Pawel Salek, pawsa@theochem.kth.se, aug 2001

   NOTE:
   this file may seem unnecessarily complex but the structure does pay off
   when implementing multiple functionals depending on different parameters.
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
#include "dcbham.h"

/* INTERFACE PART */
static int  becke_read(const char* conf_line);
static real becke_energy(DftDensProp* dp);
static void becke_first(FirstFuncDrv *ds, real factor, DftDensProp* dp);
static void becke_second(SecondFuncDrv *ds, real factor, DftDensProp* dp);
static void becke_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp);

static real b88corr_energy(DftDensProp* dp);
static void b88corr_first(FirstFuncDrv *ds, real factor, DftDensProp* dp);
static void b88corr_second(SecondFuncDrv *ds, real factor,DftDensProp* dp);
static void b88corr_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp);

/* B88CorrFunctional is a GGA correction to Becke88 functional.
   Becke(88) functional consists of Dirac rho-dependent part and GGA
   correction.
   NOTE: Dirac is so simple that one might contemplate inserting its
   expressions in-line in Becke functional.
*/
Functional BeckeFunctional = {
  "Becke",      /* name */
  fun_true,  /* gga-corrected */
  becke_read,   /* set bloody common blocks */
  NULL,         /* reporter */
  becke_energy, 
  becke_first,
  becke_second,
  becke_third
};

Functional B88CorrFunctional = {
  "B88Corr",      /* name */
  fun_true,  /* gga-corrected */
  becke_read,   /* set bloody common blocks */
  NULL,         /* reporter */
  b88corr_energy, 
  b88corr_first,
  b88corr_second,
  b88corr_third
};

/* IMPLEMENTATION PART */

static int
becke_read(const char* conf_line)
{
    dft_set_hf_weight(0.0);
  return 1;
}

static real
becke_energy(DftDensProp* dp)
{
  return SlaterFunctional.func(dp) + B88CorrFunctional.func(dp);
}

static void
becke_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
{
  SlaterFunctional.first (ds,factor,dp);
  B88CorrFunctional.first(ds,factor,dp);
}
static void
becke_second(SecondFuncDrv *ds, real factor, DftDensProp* dp)
{
  SlaterFunctional.second (ds,factor,dp);
  B88CorrFunctional.second(ds,factor,dp);
}
static void
becke_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp)
{
  SlaterFunctional.third (ds,factor, dp);
  B88CorrFunctional.third(ds,factor, dp);
}

static const real BECKE_THRESHOLD = 1e-14;
static const real BETA = 0.0042;
/* b88corr_energy:
   note that in reality E_BCK = E_BCK,alpha + E_BCK,beta
   i.e the is linear in alpha and beta densities.
*/
static real
b88corr_energy(DftDensProp* dp)
{
   real ea,eb;
   if (dp->rhob<BECKE_THRESHOLD)
     eb = 0.0;
   else {
     real xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
     real rb = pow(dp->rhob,4.0/3.0);
     real denomb = 1.0 +6.0*xb*BETA*asinh(xb);
     eb = rb*xb*xb/denomb; 
   } 
   if (dp->rhoa<BECKE_THRESHOLD) 
     ea=0;
   else {
       real xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
       real ra = pow(dp->rhoa,4.0/3.0);
       real denoma = 1.0 +6.0*BETA*xa*asinh(xa);
       ea = ra*xa*xa/denoma;
   }
   return -BETA*(ea+eb);
}

static void
b88corr_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
{
    real alpha, aa2, asha, sq1aa2;
    real alphb, ab2, ashb, sq1ab2;
    real denoma, denoma2;
    real denomb, denomb2; 
    real alphaa10, alphaa01;
    real alphab10, alphab01;
    real ffa, ff1a;
    real ffb, ff1b;

    if (dp->rhoa >BECKE_THRESHOLD) {
      alpha = dp->grada*pow(dp->rhoa,-4.0/3.0);
      aa2 = alpha*alpha;
      asha = asinh(alpha);
      alphaa10 = -4.0/3.0*alpha/dp->rhoa;   /* dalpha/drhoa   */
      alphaa01 = pow(dp->rhoa,-4.0/3.0);    /* dalpha/dgrho  */
      sq1aa2 = sqrt(1 + aa2);
      denoma= 1 + 6*alpha*BETA*asha;
      denoma2 = denoma*denoma;
      ffa  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1a = BETA*(6*aa2*BETA - sq1aa2)/(sq1aa2*denoma2);
      ds->df1000 += factor*dp->grada* ff1a *alphaa10;
      ds->df0010 += factor*(ffa + dp->grada*ff1a*alphaa01);
    }
    if (dp->rhob >BECKE_THRESHOLD) {
        alphb = dp->gradb*pow(dp->rhob,-4.0/3.0);   
        ab2 = alphb*alphb;
        ashb = asinh(alphb);
        alphab10 = -4.0/3.0*alphb/dp->rhob; 
        alphab01 = pow(dp->rhob,-4.0/3.0);
        sq1ab2 = sqrt(1 + ab2);
        denomb= 1 + 6*alphb*BETA*ashb;      
        denomb2 = denomb*denomb;  
        ffb  = -alphb*BETA/(1+6*alphb*BETA*ashb);   
        ff1b = BETA*(6*ab2*BETA - sq1ab2)/(sq1ab2*denomb2); 
        ds->df0100 += factor*dp->gradb* ff1b *alphab10;
        ds->df0001 += factor*(ffb + dp->gradb*ff1b*alphab01); 
    } 
}

static void
b88corr_second(SecondFuncDrv *ds, real factor, DftDensProp* dp)
{
    const real BETA2 = BETA*BETA;
    real alpha, a2, a3, a4, asha, asha2, sq1a2;
    real denom, denom2, denom3;
    real alpha10, alpha20, alpha01, alpha11, alpha10_2;
    real ff, ff1, ff2;
    
    if (dp->rhoa >BECKE_THRESHOLD) {    
      alpha = dp->grada*pow(dp->rhoa,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhoa;   /* dalpha/drhoa   */
      alpha20 = -7.0/3.0*alpha10/dp->rhoa; /* d²alpha/drhoa² */
      alpha01 = pow(dp->rhoa,-4.0/3.0);    /* dalpha/dgrho  */
      alpha11 = -4.0/3.0*alpha01/dp->rhoa; /* d²alpha/(drho dgrho) */
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
	     2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ds->df1000 += factor*dp->grada* ff1 *alpha10;
      ds->df0010 += factor*(ff + dp->grada*ff1*alpha01);
      ds->df1010 += factor*(ff1*alpha10 + 
                    dp->grada*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df2000 += factor*dp->grada*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0020 += factor*(2*ff1*alpha01 + 
		    dp->grada*(ff2*alpha01*alpha01));
  }
   /* note: reuse of variables for beta part, maybe a bit .. */
  if (dp->rhob >BECKE_THRESHOLD) {    
      alpha = dp->gradb*pow(dp->rhob,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhob;   
      alpha20 = -7.0/3.0*alpha10/dp->rhob; 
      alpha01 = pow(dp->rhob,-4.0/3.0);    
      alpha11 = -4.0/3.0*alpha01/dp->rhob; 
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
	     2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ds->df0100 += factor*dp->gradb* ff1 *alpha10;
      ds->df0001 += factor*(ff + dp->gradb*ff1*alpha01);
      ds->df0101 += factor*(ff1*alpha10 + 
                    dp->gradb*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df0200 += factor*dp->gradb*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0002 += factor*(2*ff1*alpha01 + 
		    dp->gradb*(ff2*alpha01*alpha01));
    }  
}

/* b88corr_third:
   Becke functional derivatives.
   Input: rho, rhogrd.
   Output:
   df1000 - d/drho F
   df1010 - d/dgrd F
   df2000 - d^2/drho^2 F
   df1010 - d^2/drho d/dgrd F
   df0020 - d^2/dgrd^2 F
   df3000 - d^3/drho^3 F

   NOTES: Instead of passing 6+ numbers, a pointer to a structure
   ggaSecDrv could be passed instead. The numbers are closely related
   to each other.

   The Becke functional is expressed as a function of grho and alpha:
   F_B = grho*f
 
   The functional derivatives are expressed through the partial derivatives.
   For example (g==grho):
   df/drho   = g df/dalpha * dalpha/drho
   d²f/drho² = g(d²f/dalpha² * (dalpha/drho)² + df/dalpha * d²alpha/drho²)
   d³f/drho³ = g(d³f/dalpha³ * (dalpha/drho)³ 
             + 3d²f/dalpha²*(dalpha/drho)*d²alpha/drho²
             + df/dalpha * d³alpha/drho³)

   The derivatives with respect to grho==g are:
   df/dgrho  = f + g df/dalpha * dalpha/dgrho
   d²f/dgrho²= 2*df/dalpha*dalpha/dg + g(d²f/dalpha²(dalpha/dg)²)

   NOTE that d²alpha/dg² = 0, and some terms are missing above.
   d³f/dgrho³= 3*d²f/dalpha²*(dalpha/dg)² (+0)
             + (+0)+ g(d³f/dalpha³(dalpha/dg)³ (+0))

   NOTE: it adds to ds, not sets it.
   This routine has unrestricted interface but works only for
   restricted case.
*/
static void
b88corr_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp)
{
    const real BETA2 = BETA*BETA;
    real alpha, a2, a3, a4, a5, asha, asha2, sq1a2;
    real denom, denom2, denom3, denom4;
    real alpha10, alpha20, alpha30, alpha01, alpha11, alpha21, alpha10_2;
    real ff, ff1, ff2, ff3;

    if (dp->rhoa > BECKE_THRESHOLD) {
      alpha = dp->grada*pow(dp->rhoa,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      a5 = alpha*a4;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhoa;   /* dalpha/drhoa   */
      alpha20 = -7.0/3.0*alpha10/dp->rhoa; /* d²alpha/drhoa² */
      alpha30 = -10./3.0*alpha20/dp->rhoa; /* d³alpha/drho³ */
      alpha01 = pow(dp->rhoa,-4.0/3.0);    /* dalpha/dgrho  */
      alpha11 = -4.0/3.0*alpha01/dp->rhoa; /* d²alpha/(drho dgrho) */
      alpha21 = alpha20/dp->grada;         /* d³alpha/(drho² dgrho) */
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;
      denom4 = denom3*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
            2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ff3 = (6*BETA2*(6 + 5*a2 + 2*a4 
	    - 12*alpha*(6 + 16*a2 + 7*a4)*BETA*asha + 
	    36*sq1a2*BETA*(-a2*(3 + 2*a2) + 6*a5*BETA*asha - 
	    pow(1 + a2,2.0)*asha2) + 
	    36*a4*BETA2*(6*(1 + a2) + (-1 + 2*a2)*asha2)))/
            (pow(1 + a2,2.5)*denom4);

      ds->df1000 += factor*dp->grada* ff1 *alpha10;
      ds->df0010 += factor*(ff + dp->grada*ff1*alpha01);
      ds->df1010 += 
	factor*(ff1*alpha10 + dp->grada*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df2000 += factor*dp->grada*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0020 += factor*(2*ff1*alpha01 + dp->grada*(ff2*alpha01*alpha01));
    

      ds->df2010 += factor*(ff2*alpha10_2 + ff1 * alpha20 +
                          dp->grada*(ff3*alpha10_2*alpha01 + 
                                 2*ff2*alpha11*alpha10 +
                                 ff2*alpha20*alpha01 + ff1*alpha21));

      ds->df1020 += factor*(2*ff2*alpha10*alpha01 + 2*ff1*alpha11 +
                          dp->grada*(ff3*alpha01*alpha01*alpha10 + 
                                 2*ff2*alpha11*alpha01));
    
      ds->df3000 += factor*dp->grada*(ff3 * alpha10_2*alpha10 +
                                3*ff2*alpha10*alpha20 +
                                ff1*alpha30);

      ds->df0030 += factor*(3*ff2*alpha01*alpha01 +
                          dp->grada*ff3*pow(alpha01,3.0));
    }
 if (dp->rhob > BECKE_THRESHOLD) {
      alpha = dp->gradb*pow(dp->rhob,-4.0/3.0);
      a2 = alpha*alpha;
      a3 = alpha*a2;
      a4 = alpha*a3;
      a5 = alpha*a4;
      asha = asinh(alpha);
      asha2= asha*asha;
      alpha10 = -4.0/3.0*alpha/dp->rhob;   /* dalpha/drhoa   */
      alpha20 = -7.0/3.0*alpha10/dp->rhob; /* d²alpha/drhoa² */
      alpha30 = -10./3.0*alpha20/dp->rhob; /* d³alpha/drho³ */
      alpha01 = pow(dp->rhob,-4.0/3.0);    /* dalpha/dgrho  */
      alpha11 = -4.0/3.0*alpha01/dp->rhob; /* d²alpha/(drho dgrho) */
      alpha21 = alpha20/dp->gradb;         /* d³alpha/(drho² dgrho) */
      alpha10_2 = alpha10*alpha10;
      sq1a2 = sqrt(1 + a2);
      denom= 1 + 6*alpha*BETA*asha;
      denom2 = denom*denom;
      denom3 = denom2*denom;
      denom4 = denom3*denom;

      ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
      ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

      ff2 = (6*BETA2*(4*alpha + a3*(3 - 12*sq1a2*BETA) + 
            2*(pow(1 + a2,1.5) - 3*a4*BETA)*asha))/
	    (pow(1 + a2,1.5)*denom3);
    
      ff3 = (6*BETA2*(6 + 5*a2 + 2*a4 
	    - 12*alpha*(6 + 16*a2 + 7*a4)*BETA*asha + 
	    36*sq1a2*BETA*(-a2*(3 + 2*a2) + 6*a5*BETA*asha - 
	    pow(1 + a2,2.0)*asha2) + 
	    36*a4*BETA2*(6*(1 + a2) + (-1 + 2*a2)*asha2)))/
            (pow(1 + a2,2.5)*denom4);

      ds->df0100 += factor*dp->gradb* ff1 *alpha10;
      ds->df0001 += factor*(ff + dp->gradb*ff1*alpha01);
      ds->df0101 += 
	factor*(ff1*alpha10 + dp->gradb*(ff2*alpha10*alpha01+ff1*alpha11));
      ds->df0200 += factor*dp->gradb*(ff2* alpha10_2 + ff1 * alpha20);
      ds->df0002 += factor*(2*ff1*alpha01 + dp->gradb*(ff2*alpha01*alpha01));
    

      ds->df0201 += factor*(ff2*alpha10_2 + ff1 * alpha20 +
                          dp->gradb*(ff3*alpha10_2*alpha01 + 
                                 2*ff2*alpha11*alpha10 +
                                 ff2*alpha20*alpha01 + ff1*alpha21));

      ds->df0102 += factor*(2*ff2*alpha10*alpha01 + 2*ff1*alpha11 +
                          dp->gradb*(ff3*alpha01*alpha01*alpha10 + 
                                 2*ff2*alpha11*alpha01));
    
      ds->df0300 += factor*dp->gradb*(ff3 * alpha10_2*alpha10 +
                                3*ff2*alpha10*alpha20 +
                                ff1*alpha30);

      ds->df0003 += factor*(3*ff2*alpha01*alpha01 +
                          dp->gradb*ff3*pow(alpha01,3.0));
    }
}
