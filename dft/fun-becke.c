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
static real becke_energy(real rhoa, real rhob, real grada, real gradb);
static void becke_first(FirstDrv *ds, real factor, 
                        real rhoa, real rhob, real ngrada, real ngradb);
static void becke_second(SecondFuncDrv *ds, real factor,
                         real rhoa, real rhob, real ngrada, real ngradb);
static void becke_third(ThirdFuncDrv *ds, real factor,
                        real rhoa, real rhob, real ngrada, real ngradb);

static real b88corr_energy(real rhoa, real rhob, real grada, real gradb);
static void b88corr_first(FirstDrv *ds, real factor, 
                          real rhoa, real rhob, real ngrada, real ngradb);
static void b88corr_second(SecondFuncDrv *ds, real factor,
                           real rhoa, real rhob, real ngrada, real ngradb);
static void b88corr_third(ThirdFuncDrv *ds, real factor,
                          real rhoa, real rhob, real ngrada, real ngradb);

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
#ifndef NO_BACKWARD_COMP
  dcrham_.hfxfac = 0;
#endif
  return 1;
}

static real
becke_energy(real rhoa, real rhob, real grada, real gradb)
{
  return SlaterFunctional.func(rhoa,rhob,grada,gradb) +
    B88CorrFunctional.func(rhoa,rhob,grada,gradb);
}

static void
becke_first(FirstDrv *ds, real factor, 
            real rhoa, real rhob, real grada, real gradb)
{
  SlaterFunctional.first (ds,factor,rhoa,rhob,grada,gradb);
  B88CorrFunctional.first(ds,factor,rhoa,rhob,grada,gradb);
}
static void
becke_second(SecondFuncDrv *ds, real factor,
             real rhoa, real rhob, real grada, real gradb)
{
  SlaterFunctional.second (ds,factor,rhoa,rhob,grada,gradb);
  B88CorrFunctional.second(ds,factor,rhoa,rhob,grada,gradb);
}
static void
becke_third(ThirdFuncDrv *ds, real factor,
            real rhoa, real rhob, real grada, real gradb)
{
  SlaterFunctional.third (ds,factor,rhoa,rhob,grada,gradb);
  B88CorrFunctional.third(ds,factor,rhoa,rhob,grada,gradb);
}

/* b88corr_energy:
   note that in reality E_BCK = E_BCK,alpha + E_BCK,beta
   i.e the is linear in alpha and beta densities.
*/
static real
b88corr_energy(real rhoa, real rhob, real grada, real gradb)
{
    const real BETA = 0.0042;
    real alpha = (grada+gradb)*pow(0.5*(rhoa+rhob),-1.0/3.0)/(rhoa+rhob);
    real asha = asinh(alpha);
    real ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
    return (grada+gradb)*ff;
}

static void
b88corr_first(FirstDrv *ds, real factor, 
              real rhoa, real rhob, real grada, real gradb)
{
    const real BETA = 0.0042;
    real alpha, a2, asha, sq1a2;
    real denom, denom2;
    real alpha10, alpha01;
    real ff, ff1;

    if (rhoa < 1e-10)
	return;

    /* no luck, just compute it.
     * f{1,2,3} are derivatives of ebck with respect to alpha.
     */
    alpha = grada*pow(rhoa,-4.0/3.0);
    a2 = alpha*alpha;
    asha = asinh(alpha);
    alpha10 = -4.0/3.0*alpha/rhoa;   /* dalpha/drhoa   */
    alpha01 = pow(rhoa,-4.0/3.0);    /* dalpha/dgrho  */
    sq1a2 = sqrt(1 + a2);
    denom= 1 + 6*alpha*BETA*asha;
    denom2 = denom*denom;

    ff  = -alpha*BETA/(1+6*alpha*BETA*asha);
    ff1 = BETA*(6*a2*BETA - sq1a2)/(sq1a2*denom2);

    ds->df1000 += factor*grada* ff1 *alpha10;
    ds->df0010 += factor*(ff + grada*ff1*alpha01);
}

static void
b88corr_second(SecondFuncDrv *ds, real factor,
               real rhoa, real rhob, real grada, real gradb)
{
    const real BETA = 0.0042, BETA2 = BETA*BETA;
    real alpha, a2, a3, a4, asha, sq1a2;
    real denom, denom2, denom3;
    real alpha10, alpha20, alpha01, alpha11, alpha10_2;
    real ff, ff1, ff2;

    if (rhoa < 1e-10)
	return;

    /* no luck, just compute it.
     * f{1,2,3} are derivatives of ebck with respect to alpha.
     */
    alpha = grada*pow(rhoa,-4.0/3.0);
    a2 = alpha*alpha;
    a3 = alpha*a2;
    a4 = alpha*a3;
    asha = asinh(alpha);
    alpha10 = -4.0/3.0*alpha/rhoa;   /* dalpha/drhoa   */
    alpha20 = -7.0/3.0*alpha10/rhoa; /* d²alpha/drhoa² */
    alpha01 = pow(rhoa,-4.0/3.0);    /* dalpha/dgrho  */
    alpha11 = -4.0/3.0*alpha01/rhoa; /* d²alpha/(drho dgrho) */
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
    
    ds->df1000 += factor*grada* ff1 *alpha10;
    ds->df0010 += factor*(ff + grada*ff1*alpha01);
    ds->df1010 += 
	factor*(ff1*alpha10 + grada*(ff2*alpha10*alpha01+ff1*alpha11));
    ds->df2000 += factor*grada*(ff2* alpha10_2 + ff1 * alpha20);
    ds->df0020 += factor*(2*ff1*alpha01 + grada*(ff2*alpha01*alpha01));
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
b88corr_third(ThirdFuncDrv *ds, real factor,
              real rhoa, real rhob, real grada, real gradb)
{
    const real BETA = 0.0042, BETA2 = BETA*BETA;
    real alpha, a2, a3, a4, a5, asha, asha2, sq1a2;
    real denom, denom2, denom3, denom4;
    real alpha10, alpha20, alpha30, alpha01, alpha11, alpha21, alpha10_2;
    real ff, ff1, ff2, ff3;

    if (rhoa < 1e-10)
	return;

    /* no luck, just compute it.
     * f{1,2,3} are derivatives of ebck with respect to alpha.
     */
    alpha = grada*pow(rhoa,-4.0/3.0);
    a2 = alpha*alpha;
    a3 = alpha*a2;
    a4 = alpha*a3;
    a5 = alpha*a4;
    asha = asinh(alpha);
    asha2= asha*asha;
    alpha10 = -4.0/3.0*alpha/rhoa;   /* dalpha/drhoa   */
    alpha20 = -7.0/3.0*alpha10/rhoa; /* d²alpha/drhoa² */
    alpha30 = -10./3.0*alpha20/rhoa; /* d³alpha/drho³ */
    alpha01 = pow(rhoa,-4.0/3.0);    /* dalpha/dgrho  */
    alpha11 = -4.0/3.0*alpha01/rhoa; /* d²alpha/(drho dgrho) */
    alpha21 = alpha20/grada;         /* d³alpha/(drho² dgrho) */
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
				   pow(1 + a2,2)*asha2) + 
		    36*a4*BETA2*(6*(1 + a2) + (-1 + 2*a2)*asha2)))/
      (pow(1 + a2,2.5)*denom4);

    ds->df1000 += factor*grada* ff1 *alpha10;
    ds->df0010 += factor*(ff + grada*ff1*alpha01);
    ds->df1010 += 
	factor*(ff1*alpha10 + grada*(ff2*alpha10*alpha01+ff1*alpha11));
    ds->df2000 += factor*grada*(ff2* alpha10_2 + ff1 * alpha20);
    ds->df0020 += factor*(2*ff1*alpha01 + grada*(ff2*alpha01*alpha01));
    

    ds->df2010 += factor*(ff2*alpha10_2 + ff1 * alpha20 +
                          grada*(ff3*alpha10_2*alpha01 + 
                                 2*ff2*alpha11*alpha10 +
                                 ff2*alpha20*alpha01 + ff1*alpha21));

    ds->df1020 += factor*(2*ff2*alpha10*alpha01 + 2*ff1*alpha11 +
                          grada*(ff3*alpha01*alpha01*alpha10 + 
                                 2*ff2*alpha11*alpha01));
    
    ds->df3000 += factor*grada*(ff3 * alpha10_2*alpha10 +
                                3*ff2*alpha10*alpha20 +
                                ff1*alpha30);

    ds->df0030 += factor*(3*ff2*alpha01*alpha01 +
                          grada*ff3*pow(alpha01,3.0));
}
