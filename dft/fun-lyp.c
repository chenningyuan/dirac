/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* fun-lyp.c:
   implementation of LYP functional and its derivatives 
   (c) Pawel Salek, pawsa@theochem.kth.se, oct 2001
*/
#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int lyp_read(const char* conf_line);
static real lyp_energy(real rhoa, real rhob, real grada, real gradb);
static void lyp_first(FirstDrv *ds, real factor, 
                      real rhoa, real rhob, real ngrada, real ngradb);
static void lyp_second(SecondFuncDrv *ds, real factor,
                       real rhoa, real rhob, real ngrada, real ngradb);
static void lyp_third(ThirdFuncDrv *ds, real factor,
                      real rhoa, real rhob, real ngrada, real ngradb);

Functional LYPFunctional = {
  "LYP",      /* name */
  fun_true,  /* gga-corrected */
  lyp_read,   /* no extra input expected, just set the common block */
  NULL, /* reporter */
  lyp_energy, 
  lyp_first,
  lyp_second,
  lyp_third
};

/* IMPLEMENTATION PART */
static int
lyp_read(const char* conf_line)
{
  /* COMPATIBILITY BLOCK: to be removed when we chop off the dalton common 
     blocks; ATM it is useful for printing the status. */
  /*  dftcom_.hfxfac = 0; 
  dfterg_.wdftx  = 0;
  dfterg_.wdftc  = 0;
  dfterg_.wdftb  = 0;
  dfterg_.wdftl  = 1;*/
  return 1;
}

/**
   The LYP formulas are based on Miehlich et al. article (CPL 157,
   p. 200, 1989).
   The implementation works also for unrestricted case (which is more
   important than you think).
*/
static real
lyp_energy(real rhoa, real rhob, real grada, real gradb)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, rhom13, denom, omega, delta, ret, ngrad2, ngrada2, ngradb2;
    real rho2, t1, t2, t3, t4, t5, t6;

    rho = rhoa+rhob;
    rho2 = rho*rho;
    ngrad2  = (grada+gradb)*(grada+gradb);
    ngrada2 = grada*grada;
    ngradb2 = gradb*gradb;
    rhom13 = pow(rho,-1.0/3.0);
    denom = 1+D*rhom13;
    omega = exp(-C*rhom13)/denom*pow(rho,-11.0/3.0);
    delta = rhom13*(C + D/denom);
 
    t1 =  pow(2.0,11.0/3.0)*CF*(pow(rhoa,8.0/3.0) +pow(rhob,8.0/3.0));
    t2 =  (47.0 - 7.0*delta)*ngrad2/18.0;
    t3 = -(2.5 -delta/18.0)*(ngrada2+ngradb2);
    t4 =  (11.0-delta)/9.0*(rhoa*ngrada2 + rhob*ngradb2)/rho;
    t5 = -2.0/3.0*rho2*ngrad2;
    t6 = ((2.0/3.0*rho2-rhoa*rhoa)*ngradb2 +
	(2.0/3.0*rho2-rhob*rhob)*ngrada2);
    ret = -A*(4*rhoa*rhob/(denom*rho)
	      +B*omega*(rhoa*rhob*(t1+t2+t3+t4)+t5+t6)); 

    return ret;
}


/* lyp_first:
   the derivatives of LYP functional wrt alpha density and gradient, suitable
   for unrestricted calculations.
   lyprho : dF/drhoa
   lypgrad: dF/dngrada
   See Doc/dft/functionals.tex for details.
*/
static void
lyp_first(FirstDrv *ds, real factor, 
          real rhoa, real rhob, real grada, real gradb)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, expcr, grad, om, dl, f1, f2, rho2, rhoa2, rhob2;
    real rho13, drho3_2, drho13, grad2, grada2, gradb2, om_1, dl_1;
    real f0_1000, f1_1000, f1_0010, f2_1000, f2_0010, sA, sA10;

    rho = rhoa + rhob;
    rho2 = rho*rho;
    rhoa2 = rhoa*rhoa;
    rhob2 = rhob*rhob;
    grad = grada+gradb;
    grad2 = grad*grad;
    rho13  = pow(rho,1.0/3.0);
    drho13 = D+rho13;
    drho3_2= drho13*drho13;
    expcr  = exp(-C/rho13);
    grada2 = grada*grada;
    gradb2 = gradb*gradb;
    sA   = (rhoa*grada2+rhob*gradb2)/rho;
    sA10 = (grada-gradb)*grad*rhob/rho2;

    om = expcr*pow(rho,-11.0/3.0)/(1.0+D/rho13); 
    om_1 = (-11*rho*rho13+rho*(C-10*D)+rho/rho13*C*D)/
	(3*pow(rho,16.0/3.0)*drho3_2)*expcr;

    dl = C/rho13 + D/(rho13*(1+D/rho13));
    dl_1 = (-C*drho3_2-D*rho/rho13)/(3*rho*rho13*drho3_2);

    f0_1000 = 4*rhob*(D*rhoa+3*rhob*drho13)*rho13/(3*rho*rho*drho3_2);

    f1 = pow(2.0, 11.0/3.0)*CF*(pow(rhoa,8.0/3.0)+pow(rhob,8.0/3.0))
	+ (47-7*dl)*grad2/18.0 +(dl-45)*(grada2+gradb2)/18.0
	  + (11-dl)*sA/9.0;
    f1_1000 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhoa,5.0/3.0)
    	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
    	+ (11-dl)*sA10/9.0;
    f1_0010 = (47-7*dl)*grad/9.0 +
	(dl-45+(22-2*dl)*rhoa/rho)*grada/9.0;
    
    f2 = -2.0/3.0*rho2*grad2
	+(2.0/3.0*rho2-rhoa2)*gradb2+(2.0/3.0*rho2-rhob2)*grada2;
    f2_1000 = -gradb*(6*rhoa*gradb+8*rho*grada)/3.0;
    f2_0010 = -4.0/3.0*rho2*grad+(2.0/3.0*rho2-rhob2)*2*grada;

    /* the final section: begin */
    ds->df1000 += factor*(-A*f0_1000-A*B*om_1*(rhoa*rhob*f1+f2)
                          -A*B*om*(rhob*f1+rhoa*rhob*f1_1000+f2_1000)); 

    ds->df0010 += factor*(-A*B*om*(rhoa*rhob*f1_0010+f2_0010));
    /* the final section: end */                          
}

/* lyp_second:
   second order derivatives of LYP functional with respect to 
   alpha density and gradient.
*/
static void
lyp_second(SecondFuncDrv *ds, real factor,
           real rhoa, real rhob, real grada, real gradb)
{
    const real A = 0.04918, B = 0.132, C = 0.2533, D = 0.349;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    
    real rho, expcr, grad, dl, f1, f2, rho2, rhoa2, rhob2,rhoab;
    real rho13, drho13, drho3_2, drho3_3, grad2, grada2, gradb2, om_1, dl_1;
    real f0_1000, f0_1100, f1_1000, f1_0100, f1_0010, f1_0001, f1_1100;
    real f1_1010, f1_1001, f1_0020, f1_0011; 
    real f2_1000, f2_0100, f2_0010, f2_0001;
    real om_2, dl_2, f0_2000, f1_2000;
    real f2_1010, f2_1001, f2_2000, f2_1100, f2_0020, f2_0011;
    real sA, sA1000, sA0100, sA2000, sA1100;
    real rff, rff_1000, rff_0100, rff_2000, rff_1100, rff_0010, rff_1010;
    real om, omx, omx1, omx2, de, de1;

    rho = rhoa + rhob;
    rho2 = rho*rho;
    rhoa2 = rhoa*rhoa;
    rhob2 = rhob*rhob;
    rhoab = rhoa*rhob;
    grad = grada+gradb;
    grad2 = grad*grad;
    rho13  = pow(rho,1.0/3.0);
    drho13 = D+rho13;
    drho3_2= drho13*drho13;
    drho3_3= drho13*drho3_2;
    expcr  = exp(-C/rho13);
    de     = C/(3*rho*rho13); /* de = expcr'/expcr */
    de1    = -4.0*C/(9*rho2*rho13);
    grada2 = grada*grada;
    gradb2 = gradb*gradb;
    sA     = (rhoa*grada2+rhob*gradb2)/rho;
    sA1000 = (grada-gradb)*grad*rhob/rho2;
    sA0100 = (gradb-grada)*grad*rhoa/rho2;
    sA2000 = -2*sA1000/rho;
    sA1100 = (rhoa-rhob)*(grada-gradb)*grad/(rho2*rho);

    dl = C/rho13 + D/(rho13*(1+D/rho13));

    omx  = pow(rho,-11.0/3.0)/(1.0+D/rho13); 
    omx1 = (-11*rho13-10*D)/(3*rho2*rho2*rho13*drho3_2);
    omx2 = 2*(77*rho/rho13+141*D*rho13+65*D*D)/
	(9*pow(rho, 16.0/3.0)*drho3_3);

    om      = expcr*omx;
    om_1 = expcr*(omx1 + omx*de);
    om_2 = expcr*(omx2 + 2*omx1*de + omx*(de1+de*de));
    om_2 = expcr*(omx2 + 2*omx1*de + omx*(de1+de*de));

    dl_1 = (-C*drho3_2-D*rho/rho13)/(3*rho*rho13*drho3_2);

    dl_2 = (4*rho*(C+D) + 2*D*(6*rho13*C*D +2*C*D*D + rho/rho13*(6*C+D)))/
	(9*rho2*rho13*drho3_3);

    f0_1000 = 4*rhob*(D*rhoa+3*rhob*drho13)*rho13/(3*rho2*drho3_2);
    f0_2000 = -8*rhob*(rhoa*D*(D+2*rho13)+3*rhob*drho13*(2*D+3*rho13))/
	(9*rho2*rho/rho13*drho3_3);

    f0_1100 = 4*(18*rhoa*rhob*rho/rho13 + 
                  rho13*(3*rhoa2+32*rhoab+3*rhob2)*D+
                  (3*rhoa2+16*rhoab+3*rhob2)*D*D)/
        (9*rho2*rho/rho13*drho3_3);

    f1 = pow(2.0, 11.0/3.0)*CF*(pow(rhoa,8.0/3.0)+pow(rhob,8.0/3.0))
	+ (47-7*dl)*grad2/18.0 +(dl-45)*(grada2+gradb2)/18.0
	  + (11-dl)*sA/9.0;
    f1_1000 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhoa,5.0/3.0)
	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
	+ (11-dl)*sA1000/9.0;
    f1_0100 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhob,5.0/3.0)
	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
	+ (11-dl)*sA0100/9.0;

    f1_0010 = (47-7*dl)*grad/9.0 +(dl-45)*grada/9 
	+(11-dl)*rhoa/rho*grada/4.5;

    f1_0001 = (47-7*dl)*grad/9.0 +(dl-45)*gradb/9 
	+(11-dl)*rhob/rho*gradb/4.5;


    f1_2000 = pow(2.0,11.0/3.0)*CF*40.0/9.0*pow(rhoa,2.0/3.0)
	- 2*sA1000*dl_1/9.0 
	+ (grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
	+ (11-dl)*sA2000/9.0;

    f1_1100 = -2*sA0100*dl_1/18.0+(grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
        -dl_1*sA1000/9.0+(11-dl)*sA1100/9.0;

    f1_1010 = (grada*(1-2*rhoa/rho)-7.0*grad)*dl_1/9.0
    	+(11-dl)*rhob/rho2*grada/4.5;
    
    f1_1001 = -(7*grada+6*gradb)*dl_1/9.0 -
        rhob*gradb*(11-dl+rho*dl_1)/(4.5*rho2);

    f1_0020 = (2-6*dl)/9.0 +(11-dl)*rhoa/rho/4.5;
    f1_0011 = (47-7*dl)/9.0;

    f2 = -2.0/3.0*rho2*grad2
	+(2.0/3.0*rho2-rhoa2)*gradb2+(2.0/3.0*rho2-rhob2)*grada2;
    f2_1000 = -gradb*(6*rhoa*gradb+8*rho*grada)/3.0;
    f2_0100 = -grada*(6*rhob*grada+8*rho*gradb)/3.0;
    f2_0010 = -4.0/3.0*rho2*grad+(2.0/3.0*rho2-rhob2)*2*grada;
    f2_0001 = -4.0/3.0*rho2*grad+(2.0/3.0*rho2-rhoa2)*2*gradb;
    f2_2000 = -gradb*(6*gradb+8*grada)/3.0;
    f2_1100 = -8*gradb*grada/3.0;
    f2_1010 = -8.0/3.0*rho*gradb;
    f2_1001 = -8.0/3.0*rho*grada - 4*rhoa*gradb;
    f2_0020 = -2*rhob2;
    f2_0011 = -4.0/3.0*rho2;

    rff      = rhoa*rhob*f1+f2; 
    rff_1000 = rhob*f1+rhoa*rhob*f1_1000+f2_1000;
    rff_0100 = rhoa*f1+rhoa*rhob*f1_0100+f2_0100;
    rff_2000 = 2*rhob*f1_1000+rhoa*rhob*f1_2000+f2_2000; 
    rff_1100 = f1+rhob*f1_0100+rhoa*f1_1000+rhoa*rhob*f1_1100+f2_1100;
    rff_0010 = rhoa*rhob*f1_0010+f2_0010;
    rff_1010 = rhob*f1_0010+rhoa*rhob*f1_1010+f2_1010;

    /* the final section: begin */
    ds->df1000 += factor*(-A*f0_1000 -A*B*(om_1*rff+om*rff_1000));

    ds->df0010 += factor*(-A*B*om*rff_0010);

    ds->df2000 += factor*(-A*f0_2000
                          -A*B*(om_2*rff+2*om_1*rff_1000+om*rff_2000)); 

    ds->df1010 += factor*(-A*B*(om_1*rff_0010+om*rff_1010));

    ds->df0020 += factor*(-A*B*om*(rhoa*rhob*f1_0020+f2_0020));

    /* remaining, coupling derivatives */    
    ds->df1100 += factor*(-A*f0_1100
                          -A*B*(om_2*rff+om_1*rff_0100+om_1*rff_1000 + 
                                om*rff_1100));

    ds->df1001 += factor*(-A*B*(om_1*(rhoa*rhob*f1_0001+f2_0001) +
                                om*(rhob*f1_0001+rhoa*rhob*f1_1001+f2_1001)));
    ds->df0011 += factor*(-A*B*om*(rhoa*rhob*f1_0011+f2_0011));
    /* the final section: end */
}

static void
lyp_third(ThirdFuncDrv *ds, real factor,
          real rhoa, real rhob, real grada, real gradb)
{
    static const real A = 0.04918, /* B = 0.132,*/ C = 0.2533, D = 0.349;
    static const real AB = 0.04918*0.132 /* A*B */;
    const real CF = 0.3*pow(3*M_PI*M_PI,2.0/3.0);
    real rho, expcr, grad, om, dl, f1, f2, rho2, rho3, rhoa2, rhob2, rhoab;
    real rho13, drho13, drho3_2, drho3_3, grad2, grada2, gradb2, om_1, dl_1;
    real f0_1000, f0_2000, f0_1100;
    real f1_1000, f1_0100, f1_0010, f1_0001, f1_1010, f1_0020, f1_0011;
    real f1_2000, f1_1100, f1_0110, f1_1001;
    real f1_3000, f1_2100, f1_2010, f1_2001, f1_1110, f1_0120, f1_1020;
    real f2_1000, f2_0100, f2_0010, f2_0001, f2_1010, f2_1001, f2_2000;
    real f2_0110, f2_0011, f2_2010, f2_2001, f2_1100, f2_0020, f2_0120, 
        f2_1110;
    real om_2, dl_2;
    real f0_2100, f0_3000, om_3;
    real dl_3;
    real sA, sA1000, sA0100, sA2000, sA1100, sA3000;
    real rff, rff_1000, rff_0100, rff_0010, rff_0001;
    real rff_2000, rff_1100, rff_1010, rff_1001, rff_0110;
    real rff_3000, rff_2100, rff_2010, rff_2001, rff_1110;
    real omx, omx1, omx2, omx3, de, de1, de2;

    rho = rhoa + rhob;
    rho2 = rho*rho;
    rho3 = rho*rho2;
    rhoa2 = rhoa*rhoa;
    rhob2 = rhob*rhob;
    rhoab = rhoa*rhob;
    grad = grada+gradb;
    grad2 = grad*grad;
    rho13  = pow(rho,1.0/3.0);
    drho13 = D+rho13;
    drho3_2= drho13*drho13;
    drho3_3= drho13*drho3_2;
    expcr  = exp(-C/rho13);
    de     = C/(3*rho*rho13); /* de = expcr'/expcr */
    de1    = -4.0*C/(9*rho2*rho13);
    de2    = -de1*7.0/(3*rho);
    grada2 = grada*grada;
    gradb2 = gradb*gradb;
    sA     = (rhoa*grada2+rhob*gradb2)/rho;
    sA1000 = (grada-gradb)*grad*rhob/rho2;
    sA0100 = (gradb-grada)*grad*rhob/rho2;
    sA2000 = -2*sA1000/rho;
    sA1100 = (rhoa-rhob)*(grada-gradb)*grad/(rho2*rho);
    sA3000 = -3*sA2000/rho;

    om = expcr*pow(rho,-11.0/3.0)/(1.0+D/rho13); 
    dl = C/rho13 + D/(rho13*(1+D/rho13));
    
    omx  = pow(rho,-11.0/3.0)/(1.0+D/rho13); 
    omx1 = (-11*rho13-10*D)/(3*rho2*rho2*rho13*drho3_2);
    omx2 = 2*(77*rho/rho13+141*D*rho13+65*D*D)/
	(9*pow(rho, 16.0/3.0)*drho3_2*drho13);
    omx3 = -2*(1309*rho+3616*rho/rho13*D+3350*rho13*D*D+1040*D*D*D)/
	(27*pow(rho,19.0/3.0)*(drho3_2*drho3_2));

    om_1 = expcr*(omx1 + omx*de);
    om_2 = expcr*(omx2 + 2*omx1*de + omx*(de1+de*de));
    om_3 = expcr*(omx3+3*omx2*de+3*omx1*(de1+de*de)
		     +omx*(de2 + 3*de1*de+ de*de*de));

    dl_1 = (-C*drho3_2-D*rho/rho13)/(3*rho*rho13*drho3_2);

    dl_2 = (4*rho*(C+D) + 2*D*(6*rho13*C*D +2*C*D*D + rho/rho13*(6*C+D)))/
	(9*rho2*rho13*drho3_2*drho13);

    dl_3 = (-2*rho/rho13*D*D*(84*C+5*D)-4*D*(rho*(28*C+8*D)+7*C*D*D*D)
	       -28*rho13*(4*C*D*D*D+rho*(C+D)))/
	(27*rho2*rho*rho13*drho3_2*drho3_2);

    f0_1000 = 4*rhob*(D*rhoa+3*rhob*drho13)*rho13/(3*rho*rho*drho3_2);
    f0_2000 = -8*rhob*(rhoa*D*(D+2*rho13)+3*rhob*drho13*(2*D+3*rho13))/
	(9*rho2*rho/rho13*drho3_2*drho13);

    f0_1100 = 4*(18*rhoa*rhob*rho/rho13 + 
                  rho13*(3*rhoa2+32*rhoab+3*rhob2)*D+
                  (3*rhoa2+16*rhoab+3*rhob2)*D*D)/
        (9*rho2*rho/rho13*drho3_3);

    f0_2100 = 8*(-54*rhoa2*rhob-27*rhoa*rhob2+27*rhob2*rhob
                -2*rho/rho13*(3*rhoa2+65*rhoab-30*rhob2)*D 
                -rho13*(9*rhoa2+110*rhoab-45*rhob2)*D*D
                +(12*rhob2-3*rhoa2-31*rhoab)*D*D*D)/
        (27*rho2*rho2/rho13*drho3_3*drho13);
                
    f0_3000 = 8*rhob*(81*rhoa*rhob+81*rhob2+2*rho/rho13*(7*rhoa+99*rhob)*D
		      +2*rho13*(8*rhoa+81*rhob)*D*D+(5*rhoa+45*rhob)*D*D*D)/
	(27*pow(rho,11.0/3.0)*drho3_2*drho3_2);

    f1 = pow(2.0, 11.0/3.0)*CF*(pow(rhoa,8.0/3.0)+pow(rhob,8.0/3.0))
	+ (47-7*dl)*grad2/18.0 +(dl-45)*(grada2+gradb2)/18.0
	  + (11-dl)*sA/9.0;
    f1_1000 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhoa,5.0/3.0)
	+ (grada2+gradb2-7*grad*grad-2*sA)*dl_1/18.0
	+ (11-dl)*sA1000/9.0;
    f1_0100 = pow(2.0,11.0/3.0)*CF*8.0/3.0*pow(rhob,5.0/3.0)
	+ (grada2+gradb2-7*grad2-2*sA)*dl_1/18.0
	+ (11-dl)*sA0100/9.0;
    f1_0010 = (47-7*dl)*grad/9.0 +
	(dl-45+(22-2*dl)*rhoa/rho)*grada/9.0;
    f1_0001 = (47-7*dl)*grad/9.0 +
	(dl-45+(22-2*dl)*rhob/rho)*gradb/9.0;
    f1_2000 = pow(2.0,11.0/3.0)*CF*40.0/9.0*pow(rhoa,2.0/3.0)
	- sA1000*dl_1/9.0 
	+ (grada2+gradb2-7*grad*grad-2*sA)*dl_2/18.0
	- dl_1*sA1000 + (11-dl)*sA2000/9.0;

    f1_1100 = -2*sA0100*dl_1/18.0+(grada2+gradb2-7*grad2-2*sA)*dl_2/18.0
        -dl_1*sA1000/9.0+(11-dl)*sA1100/9.0;

    f1_1010 = (grada*(2.0*rhoa/rho-1.0)-7.0*grad)*dl_1/9.0
	+((11.0-dl)*rhob/rho2)*grada/4.5;
    f1_1001 = -(7*grada+6*gradb)*dl_1/9.0 -
        rhob*gradb*(11-dl+rho*dl_1)/(4.5*rho2);

    f1_0110 = (-(6*grada+7*gradb)*dl_1-2*rhoa*grada*(11-dl+rho*dl_1)/rho2)/9.;

    f1_0020 = (2-6*dl)/9.0 +(11-dl)*rhoa/rho/4.5;
    f1_0011 = (47-7*dl)/9.0;

    f1_3000 = pow(2.0,11.0/3.0)*CF*80.0/27.0*pow(rhoa,-1.0/3.0)
	- (sA2000*dl_1+sA1000*dl_2)/9.0 
	+ (grada2+gradb2-7*grad*grad-2*sA)*dl_3/18.0
	- sA*dl_2/9.0
	- sA1000*dl_2/9.0 -2*dl_1*sA2000/9.0 
	+ (11-dl)*sA3000/9.0;
    f1_2100 = -(sA1100*dl_1+sA1000*dl_2)/4.5
	+ (grada2+gradb2-7*grad2-2*sA)*dl_3/18.0 - sA0100*dl_2/9.0
	+ (11-dl)*sA3000/9.0 -dl_1*sA2000/9.0;

    f1_2010 = (grada-7.0*grad)*dl_2/9.0
	+(-dl_1*rhob/rho2-2*(11-dl)*rhob/(rho2*rho)
	  -dl_2*rhoa/rho-dl_1*rhob/rho2)*grada/4.5;

    f1_2001 = (-(7*grada+6*gradb)*dl_2+
               rhob*gradb*(4*(11-dl+rho*dl_1)-2*rho2*dl_2)/rho3)/9.0;

    f1_1110= (2*grada*(rhoa-rhob)*(11-dl+rho*dl_1)-
              rho2*(2*(4*rhoa+3*rhob)*grada+7*rho*gradb)*dl_2)/(9*rho3);

    f1_1020 = -2*dl_1/3.0 -dl_1*rhoa/rho/4.5 +(11-dl)*rhob/rho2/4.5;
    f1_0120 = -dl_1*(6.0/9.0+rhoa/(4.5*rho)) -(11-dl)*rhoa/(4.5*rho2);

    f2 = -2.0/3.0*rho2*grad2
	+(2.0/3.0*rho2-rhoa2)*gradb2+(2.0/3.0*rho2-rhob2)*grada2;
    f2_1000 = -gradb*(6*rhoa*gradb+8*rho*grada)/3.0;
    f2_0100 = -grada*(6*rhob*grada+8*rho*gradb)/3.0;
    f2_0010 = -4.0/3.0*rho2*grad+(2.0/3.0*rho2-rhob2)*2*grada;
    f2_0001 = -4.0/3.0*rho2*grad+(2.0/3.0*rho2-rhoa2)*2*gradb;
    f2_2000 = -gradb*(6*gradb+8*grada)/3.0;
    f2_1100 = -8.0/3.0*grada*gradb;
    f2_1010 = -8.0/3.0*rho*gradb;
    f2_1001 = -8.0/3.0*rho*grada-4.0*rhoa*gradb;
    f2_0110 = -8.0/3.0*rho*gradb-4.0*rhob*grada;
    f2_0011 = -4.0/3.0*rho2;
    f2_2010 = -8.0/3.0*gradb;
    f2_2001 = -8.0/3.0*grada-4*gradb;
    f2_0020 = -2*rhob2;
    f2_1110 = -8.0/3.0*gradb;
    f2_0120 = -4*rhob;

    rff      = rhoab*f1+f2;
    rff_1000 = rhob*f1+rhoab*f1_1000+f2_1000;
    rff_0100 = rhoa*f1+rhoab*f1_0100+f2_0100;
    rff_2000 = 2*rhob*f1_1000+rhoa*rhob*f1_2000+f2_2000;
    rff_1100 = f1+rhob*f1_0100+rhoa*f1_1000+rhoab*f1_1100+f2_1100;
    rff_2100 = 2*(f1_1000+rhob*f1_1100)+
        (rhoa*f1_2000+rhoab*f1_2100);
    rff_3000 = 2*rhob*f1_2000+rhob*(f1_2000+rhoa*f1_3000);
    rff_0010 = rhoab*f1_0010+f2_0010; 
    rff_0001 = rhoab*f1_0001 + f2_0001; 
    rff_1010 = rhob*f1_0010+rhoab*f1_1010+f2_1010;
    rff_1001 = rhob*f1_0001+rhoab*f1_1001+f2_1001;
    rff_0110 = rhoa*f1_0010+rhoab*f1_0110+f2_0110;
    rff_2010 = 2*rhob*f1_1010+rhoab*f1_2010+f2_2010;
    rff_2001 = 2*rhob*f1_1001+rhoab*f1_2001+f2_2001; 
    rff_1110 = f1_0010+rhob*f1_0110 +rhoa*f1_1010+rhoab*f1_1110+f2_1110;

    /* the final section: begin */
    ds->df1000 += factor*(-A*f0_1000-AB*om_1*rff-AB*om*rff_1000);

    ds->df0010 += factor*(-AB*om*(rhoa*rhob*f1_0010+f2_0010));

    ds->df2000 += factor*(-A*f0_2000
			  -AB*(om_2*rff+2*om_1*rff_1000+om*rff_2000));

    ds->df1010 += factor*(-AB*(om_1*rff_0010+om*rff_1010));

    ds->df0020 += factor*(-AB*om*(rhoa*rhob*f1_0020+f2_0020));

    ds->df3000 += factor*(-A*f0_3000
			  -AB*(om_3*rff+3*om_2*rff_1000 +
				3*om_1*rff_2000+om*rff_3000));
    ds->df2010 += factor*(-AB*(om_2*rff_0010+2*om_1*rff_1010 + om*rff_2010));
    ds->df1020 += factor*(-AB*(om_1*(rhoa*rhob*f1_0020+f2_0020)+
				om*rhob*(f1_0020+rhoa*f1_1020)));

    /* remaining derivatives of LYP */
    ds->df1100 += factor*(-A*f0_1100
                          -AB*(om_2*rff+om_1*rff_0100+om_1*rff_1000 + 
                               om*rff_1100));
    ds->df1001 += factor*(-AB*(om_1*rff_0001 + om*rff_1001));

    ds->df0011 += factor*(-AB*om*(rhoab*f1_0011+f2_0011));
    ds->df2100 += factor*(-A*f0_2100
			  -AB*(om_3*rff+om_2*(2*rff_1000+rff_0100)+
                               om_1*(rff_2000+2*rff_1100)+om*rff_2100));

    ds->df0120 += factor*(-AB*(om_1*(rhoab*f1_0020+f2_0020)+
                               om*(rhoa*f1_0020+rhoab*f1_0120+f2_0120)));

    ds->df1110 += factor*(-AB*(om_2*rff_0010+om_1*rff_0110+om_1*rff_1010 + 
                               om*rff_1110));
    ds->df2001 += factor*(-AB*(om_2*rff_0001+2*om_1*rff_1001+om*rff_2001));

    ds->df1011 += factor*(-AB*om_1*(rhoab*f1_0011+f2_0011)
                          -AB*om*(rhob*f1_0011-rhoab*7.0/9.0*dl_1
                                  -8.0/3.0*rho));

    /* the final section: end */
    /* drvs->df0030 += 0; *
     * drvs->df0021 += 0; */
}
