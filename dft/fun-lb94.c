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
/* fun-lb94.c:
   implementation of Exchange-correlation potential with correct 
   asymptotic behavior by R. van Leeuwen and E. J. Baerends;
   (The LB94 functional) and its derivatives.
   (c) B. Jansik, brano@theochem.kth.se, dec 2002
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int  lb94_isgga(void) { return 1; }
static int  lb94_read(const char* conf_line);
static real lb94_energy(DftDensProp* dens_prop);
static void lb94_first(FirstFuncDrv *ds, real factor, 
                       DftDensProp* dens_prop);
static void lb94_second(SecondFuncDrv *ds, real factor,
                        DftDensProp* dens_prop);

static void lb94_third(ThirdFuncDrv *ds, real factor,
                       DftDensProp* dens_prop);

#ifdef FOURTH_ORDER_DERIVATIVES
static void lb94_fourth(FourthFuncDrv *ds, real factor,
                        DftDensProp* dens_prop);
#endif

Functional LB94Functional = {"LB94",      /* name */
                             lb94_isgga,  /* gga-corrected */
                             lb94_read,   /* set bloody common blocks */
                             NULL,         /* reporter */
                             lb94_energy, 
                             lb94_first,
                             lb94_second,
                             lb94_third
#ifdef FOURTH_ORDER_DERIVATIVES
                             ,lb94_fourth,
#endif
};

/* IMPLEMENTATION PART */

static int
lb94_read(const char* conf_line)
{
    dft_set_hf_weight(0.0);
    return 1;
}

/* lb94_energy:
   note that in reality E_LB94 = E_LB94,alpha + E_LB94,beta
   i.e the is linear in alpha and beta densities.

   lb94 threshold is needed to avoid numerical problems on 0/0
   divisions.  The problems are small but it is better to be on the
   safe side.
*/
static const real LB94_THRESHOLD = 1e-14;
static const real BETA = 0.050;
/* BETA:
 * note that BETA parameter in LB94 functional is
 * different to BETA in Becke(88) functional.
 */

static real
lb94_energy(DftDensProp* dp)
{
    real ea,eb;
    if (dp->rhob<LB94_THRESHOLD)
        eb = 0.0;
    else {
        real xb = dp->gradb*pow(dp->rhob,-4.0/3.0);
        real rb = pow(dp->rhob,4.0/3.0);
        real denomb = 1.0 +3.0*xb*BETA*asinh(xb);
        eb = rb*xb*xb/denomb; 
    } 
    if (dp->rhoa<LB94_THRESHOLD) 
        ea=0.0;
    else {
        real xa = dp->grada*pow(dp->rhoa,-4.0/3.0);
        real ra = pow(dp->rhoa,4.0/3.0);
        real denoma = 1.0 +3.0*xa*BETA*asinh(xa);
        ea = ra*xa*xa/denoma;
    }
    return SlaterFunctional.func(dp)+VWNFunctional.func(dp)
        -BETA*(ea+eb);
}

/* LB94 functional derivatives
 * 
 * Following partitioning of LB94 has been used:
 * lb94= -BETA*rho^(4/3)*F( X(rho, grad) )
 *
 * where auxaliary functions X, F 
 * are defined as:
 *  
 * X(rho,grad) = grad/rho^(4/3)
 * F(X) = (X^2) / 1 + 3*X*ArcSinh(X)
 *
 * Derivatives were then expressed in terms of derivatives of
 * auxaliary functions.
 *
 * **********************************************************
 * Variable naming convention:
 *
 * NAME       :   Variable name
 * NAME_XY    :   X-th partial derivative with respect to rho
 * 	 	  Y-th partial derivative with respect to grad
 * 	 	  of Variable (function) NAME
 * 	 	  
 * NAME_pX    :   Xth power of NAME  (NAME^X)
 * NAME_pXfY  :   NAME^(X/Y)
 *
 * NAME_X'X'  :   X-th partial derivative with respect to X function
 * 	 	  
 * Exception  :   ASINH_X does not follow these rules
 * 		  ASINH_X = asinh(X)
 * **********************************************************
 */

static void
lb94_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
{

/* Declarations */
    real    ASINH_X;
    real    B;
    real    ComDenom;
    real    ComDenom_p2;
    real    F;
    real    F_1X;
    real    SQ1X2;
    real    SU1X2;
    real    X;
    real    X_01;
    real    X_10;
    real    X_p2;
    real    gro;
    real    ro;
    real    ro_p1f3;
    real    ro_p2;
    real    ro_p4f3;
    real    ro_p7f3;

/* Setting up Beta (In the following, BETA is called just 'B' ) */
    B = BETA;

    if (dp->rhoa >LB94_THRESHOLD) {
    
/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities 
 */
        ro = dp->rhoa;
        gro = dp->grada;


/* Powers of ro */
        ro_p1f3=pow(ro,(1.0/3.0));
        ro_p4f3=ro*ro_p1f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;


/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;

/* sqrt(1+ X*X) and its powers */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);

        ASINH_X=asinh(X);

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_01=1.0/ro_p4f3;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom -1.0) )/ ComDenom_p2;

/* Derivatives expressions */
        ds->df1000 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

        ds->df0010 += -B*factor*(ro_p4f3*F_1X*X_01);
	
    }

    /* note: reuse of variables for beta part, maybe a bit .. */
    if (dp->rhob >LB94_THRESHOLD) {
	    
/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities 
 */
        ro=dp->rhob;
        gro=dp->gradb;


/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p4f3=ro*ro_p1f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;

/* sqrt(1+ X*X) and its powers */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);

        ASINH_X=asinh(X);

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_01=1.0/ro_p4f3;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom -1.0) )/ ComDenom_p2;

/* Derivatives expressions */
        ds->df0100 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

        ds->df0001 += -B*factor*(ro_p4f3*F_1X*X_01);
    } 
    SlaterFunctional.first(ds, factor, dp);
    VWNFunctional.first(ds,   factor, dp);
}

static void
lb94_second(SecondFuncDrv *ds, real factor, DftDensProp* dp)
{

/* Declarations */
    real  ASINH_X;
    real  B;
    real  B_p2;
    real  B_p3;
    real  ComDenom;
    real  ComDenom_p2;
    real  ComDenom_p3;
    real  F;
    real  F_1X;
    real  F_2X;
    real  SQ1X2;
    real  SU1X2;
    real  SU1X2_p3f2;
    real  X;
    real  X_01;
    real  X_01_p2;
    real  X_10;
    real  X_10_p2;
    real  X_11;
    real  X_20;
    real  X_p2;
    real  X_p3;
    real  X_p4;
    real  gro;
    real  ro;
    real  ro_p10f3;
    real  ro_p1f3;
    real  ro_p2;
    real  ro_p2f3;
    real  ro_p4f3;
    real  ro_p7f3;
	

/* Setting up Beta (In the following, BETA is called just 'B' ) and its powers */
    B=BETA;
    B_p2=B*B;
    B_p3=B_p2*B;


    if (dp->rhoa >LB94_THRESHOLD) {    

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities */
        ro=dp->rhoa;
        gro=dp->grada;


/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;

/* sqrt(1+ X*X) and its powers */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;


        ASINH_X=asinh(X);

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2.0*X_p2*(SQ1X2 -9.0*B) + 3.0*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2.0+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

/* Derivatives expressions */

        ds->df1000 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

        ds->df0010 += -B*factor*(ro_p4f3*F_1X*X_01);

        ds->df2000 += -B*factor*( ((4.0/9.0)*F/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_1X*X_10 + 
                                  ro_p4f3*F_2X*X_10_p2 + ro_p4f3*F_1X*X_20 );

        ds->df0020 += -B*factor*( ro_p4f3*F_2X*X_01_p2 );

        ds->df1010 += -B*factor*( (4.0/3.0)*ro_p1f3*F_1X*X_01 + ro_p4f3*F_2X*X_01*X_10 + 
                                  ro_p4f3*F_1X*X_11 );
 
    }
  
    /* note: reuse of variables for beta part, maybe a bit .. */
    if (dp->rhob >LB94_THRESHOLD) {    

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities */	  
        ro=dp->rhob;
        gro=dp->gradb;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;

/* 1+ X*X and its powers */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;


        ASINH_X=asinh(X);

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

/* Derivatives expressions */
        ds->df0100 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

        ds->df0001 += -B*factor*(ro_p4f3*F_1X*X_01);

        ds->df0200 += -B*factor*( ((4.0/9.0)*F/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_1X*X_10 + 
                                  ro_p4f3*F_2X*X_10_p2 + ro_p4f3*F_1X*X_20 );

        ds->df0002 += -B*factor*( ro_p4f3*F_2X*X_01_p2 );

        ds->df0101 += -B*factor*( (4.0/3.0)*ro_p1f3*F_1X*X_01 + ro_p4f3*F_2X*X_01*X_10 + 
                                  ro_p4f3*F_1X*X_11 );
    }  
    SlaterFunctional.second(ds, factor, dp);
    VWNFunctional.second(ds,   factor, dp);
}

 
static void
lb94_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp) {

/* Declarations */
    real  ASINH_X;
    real  ASINH_X_p2;
    real  B;
    real  B_p2;
    real  B_p3;
    real  ComDenom;
    real  ComDenom_p2;
    real  ComDenom_p3;
    real  ComDenom_p4;
    real  F;
    real  F_1X;
    real  F_2X;
    real  F_3X;
    real  SQ1X2;
    real  SU1X2;
    real  SU1X2_p3f2;
    real  SU1X2_p5f2;
    real  X;
    real  X_01;
    real  X_01_p2;
    real  X_01_p3;
    real  X_10;
    real  X_10_p2;
    real  X_10_p3;
    real  X_11;
    real  X_20;
    real  X_21;
    real  X_30;
    real  X_p2;
    real  X_p3;
    real  X_p4;
    real  X_p5;
    real  X_p6;
    real  gro;
    real  ro;
    real  ro_p10f3;
    real  ro_p13f3;
    real  ro_p1f3;
    real  ro_p2;
    real  ro_p2f3;
    real  ro_p5f3;
    real  ro_p4f3;
    real  ro_p7f3;
    SecondFuncDrv ds_sec;

/* Setting up Beta (In the following, BETA is called just 'B' ) and its powers */
    B=BETA;
    B_p2=B*B;
    B_p3=B_p2*B;

    drv2_clear(&ds_sec);
    LB94Functional.second(&ds_sec, factor, dp);

    ds->df1000 += ds_sec.df1000;
    ds->df0010 += ds_sec.df0010;
    ds->df2000 += ds_sec.df2000;
    ds->df0020 += ds_sec.df0020;
    ds->df1010 += ds_sec.df1010;

    ds->df0100 += ds_sec.df0100;
    ds->df0001 += ds_sec.df0001;
    ds->df0200 += ds_sec.df0200;
    ds->df0002 += ds_sec.df0002;
    ds->df0101 += ds_sec.df0101;

    if (dp->rhoa >LB94_THRESHOLD) {    

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities  */
        ro = dp->rhoa;
        gro = dp->grada;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p5f3=ro*ro_p2f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;
        ro_p13f3=ro_p10f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;
        X_p5=X_p4*X;
        X_p6=X_p5*X;

/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;
        SU1X2_p5f2=SU1X2_p3f2*SU1X2;

/* asinh(X) and its powers */
        ASINH_X=asinh(X);
        ASINH_X_p2=ASINH_X*ASINH_X;

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;
        ComDenom_p4=ComDenom_p3*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_30=-(280.0/27.0)*gro/ro_p13f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;
        X_21=(28.0/9.0)/ro_p10f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;
        X_10_p3=X_10_p2*X_10;
        X_01_p3=X_01_p2*X_01;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

        F_3X=(3.0*B*(-X*(18.0 + 54.0*X_p6*B_p2 + X_p2*(26.0 - 72.0*B*SQ1X2) + 
                         X_p4*(11.0 -54.0*B*SQ1X2 + 54.0*B_p2)) + 6.0*(-SQ1X2 - 
                                                                       2.0*X_p2*SQ1X2 + 4.0*X_p6*B + X_p4*(-SQ1X2 + 7.0*B + 9.0*SQ1X2*B_p2)) 
                     *ASINH_X + 9.0*X_p5*(4.0 + X_p2)*B_p2*ASINH_X_p2)) / (SU1X2_p5f2*ComDenom_p4);

/* Derivatives expressions */
/*
  ds->df1000 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

  ds->df0010 += -B*factor*(ro_p4f3*F_1X*X_01);

  ds->df2000 += -B*factor*( ((4.0/9.0)*F/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_1X*X_10 + 
  ro_p4f3*F_2X*X_10_p2 + ro_p4f3*F_1X*X_20 );

  ds->df0020 += -B*factor*( ro_p4f3*F_2X*X_01_p2 );

  ds->df1010 += -B*factor*( (4.0/3.0)*ro_p1f3*F_1X*X_01 + ro_p4f3*F_2X*X_01*X_10 + 
  ro_p4f3*F_1X*X_11 );
*/
	
 
        ds->df3000 += -B*factor*( -((8.0/27.0)*F/ro_p5f3) + ((4.0/3.0)*F_1X*X_10/ro_p2f3) + 
                                  4.0*ro_p1f3*F_2X*X_10_p2 + ro_p4f3*F_3X*X_10_p3 + 4.0*ro_p1f3*F_1X*X_20 + 
                                  3.0*ro_p4f3*F_2X*X_10*X_20 + ro_p4f3*F_1X*X_30 );
		
        ds->df2010 += -B*factor*( ((4.0/9.0)*F_1X*X_01/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_2X*X_01*X_10 + 
                                  ro_p4f3*F_3X*X_01*X_10_p2 + (8.0/3.0)*ro_p1f3*F_1X*X_11 + 
                                  2.0*ro_p4f3*F_2X*X_10*X_11 + ro_p4f3*F_2X*X_01*X_20 + 
                                  ro_p4f3*F_1X*X_21 );

        ds->df1020 += -B*factor*( (4.0/3.0)*ro_p1f3*F_2X*X_01_p2 +  
                                  ro_p4f3*F_3X*X_01_p2*X_10 +  
                                  2*ro_p4f3*F_2X*X_01*X_11  );

        ds->df0030 += -B*factor*( ro_p4f3*F_3X*X_01_p3  );

    }
  
    if (dp->rhob >LB94_THRESHOLD) {    

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities  */
        ro = dp->rhob;
        gro = dp->gradb;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p5f3=ro*ro_p2f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;
        ro_p13f3=ro_p10f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;
        X_p5=X_p4*X;
        X_p6=X_p5*X;

/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;
        SU1X2_p5f2=SU1X2_p3f2*SU1X2;

/* asinh(X) and its powers */
        ASINH_X=asinh(X);
        ASINH_X_p2=ASINH_X*ASINH_X;

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;
        ComDenom_p4=ComDenom_p3*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_30=-(280.0/27.0)*gro/ro_p13f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;
        X_21=(28.0/9.0)/ro_p10f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;
        X_10_p3=X_10_p2*X_10;
        X_01_p3=X_01_p2*X_01;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

        F_3X=(3.0*B*(-X*(18.0 + 54.0*X_p6*B_p2 + X_p2*(26.0 - 72.0*B*SQ1X2) + 
                         X_p4*(11.0 -54.0*B*SQ1X2 + 54.0*B_p2)) + 6.0*(-SQ1X2 - 
                                                                       2.0*X_p2*SQ1X2 + 4.0*X_p6*B + X_p4*(-SQ1X2 + 7.0*B + 9.0*SQ1X2*B_p2)) 
                     *ASINH_X + 9.0*X_p5*(4.0 + X_p2)*B_p2*ASINH_X_p2)) / (SU1X2_p5f2*ComDenom_p4);

/* Derivatives expressions */
/*
  ds->df0100 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

  ds->df0001 += -B*factor*(ro_p4f3*F_1X*X_01);

  ds->df0200 += -B*factor*( ((4.0/9.0)*F/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_1X*X_10 + 
  ro_p4f3*F_2X*X_10_p2 + ro_p4f3*F_1X*X_20 );

  ds->df0002 += -B*factor*( ro_p4f3*F_2X*X_01_p2 );

  ds->df0101 += -B*factor*( (4.0/3.0)*ro_p1f3*F_1X*X_01 + ro_p4f3*F_2X*X_01*X_10 + 
  ro_p4f3*F_1X*X_11 );
 
*/
	

        ds->df0300 += -B*factor*( -((8.0/27.0)*F/ro_p5f3) + ((4.0/3.0)*F_1X*X_10/ro_p2f3) + 
                                  4.0*ro_p1f3*F_2X*X_10_p2 + ro_p4f3*F_3X*X_10_p3 + 4.0*ro_p1f3*F_1X*X_20 + 
                                  3.0*ro_p4f3*F_2X*X_10*X_20 + ro_p4f3*F_1X*X_30 );
		
        ds->df0201 += -B*factor*( ((4.0/9.0)*F_1X*X_01/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_2X*X_01*X_10 + 
                                  ro_p4f3*F_3X*X_01*X_10_p2 + (8.0/3.0)*ro_p1f3*F_1X*X_11 + 
                                  2.0*ro_p4f3*F_2X*X_10*X_11 + ro_p4f3*F_2X*X_01*X_20 + 
                                  ro_p4f3*F_1X*X_21 );

        ds->df0102 += -B*factor*( (4.0/3.0)*ro_p1f3*F_2X*X_01_p2 +  
                                  ro_p4f3*F_3X*X_01_p2*X_10 +  
                                  2*ro_p4f3*F_2X*X_01*X_11  );

        ds->df0003 += -B*factor*( ro_p4f3*F_3X*X_01_p3 );
    }
    SlaterFunctional.third(ds, factor, dp);
    VWNFunctional.third(ds,   factor, dp);

}

#ifdef FOURTH_ORDER_DERIVATIVES  
static void
lb94_fourth(FourthFuncDrv *ds, real factor, DftDensProp* dp) {

/* Declarations */
    real  ASINH_X, ASINH_X_p2;
    real  B, B_p2, B_p3;
    real  ComDenom;
    real  ComDenom_p2;
    real  ComDenom_p3;
    real  ComDenom_p4;
    real  ComDenom_p5;
    real  F;
    real  F_1X;
    real  F_2X;
    real  F_3X;
    real  F_4X;
    real  SQ1X2;
    real  SU1X2;
    real  SU1X2_p3;
    real  SU1X2_p3f2;
    real  SU1X2_p5f2;
    real  SU1X2_p7f2;
    real  X;
    real  X_01;
    real  X_01_p2;
    real  X_01_p3;
    real  X_01_p4;
    real  X_10;
    real  X_10_p2;
    real  X_10_p3;
    real  X_10_p4;
    real  X_11;
    real  X_11_p2;
    real  X_20;
    real  X_20_p2;
    real  X_21;
    real  X_30;
    real  X_31;
    real  X_40;
    real  X_p2;
    real  X_p3;
    real  X_p4;
    real  X_p5;
    real  X_p6;
    real  gro;
    real  ro;
    real  ro_p10f3;
    real  ro_p13f3;
    real  ro_p16f3;
    real  ro_p1f3;
    real  ro_p2;
    real  ro_p2f3;
    real  ro_p5f3;
    real  ro_p4f3;
    real  ro_p7f3;
    real  ro_p8f3;

/* Setting up Beta (In the following, BETA is called just 'B' ) and its powers */
    B=BETA;
    B_p2=B*B;
    B_p3=B_p2*B;

    if (dp->rhoa >LB94_THRESHOLD) {    

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities  */
        ro = dp->rhoa;
        gro = dp->grada;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p5f3=ro*ro_p2f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p8f3=ro_p7f3*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;
        ro_p13f3=ro_p10f3*ro;
        ro_p16f3=ro_p13f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;
        X_p5=X_p4*X;
        X_p6=X_p5*X;

/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;
        SU1X2_p5f2=SU1X2_p3f2*SU1X2;
        SU1X2_p3 = SU1X2*SU1X2*SU1X2;
        SU1X2_p7f2=SU1X2_p5f2*SU1X2;

/* asinh(X) and its powers */
        ASINH_X=asinh(X);
        ASINH_X_p2=ASINH_X*ASINH_X;

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;
        ComDenom_p4=ComDenom_p3*ComDenom;
        ComDenom_p5=ComDenom_p4*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_30=-(280.0/27.0)*gro/ro_p13f3;
        X_40=(3640.0/81.0)*gro/ro_p16f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;
        X_21=(28.0/9.0)/ro_p10f3;
        X_31=-(280.0/27.0)/ro_p13f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;
        X_10_p3=X_10_p2*X_10;
        X_01_p3=X_01_p2*X_01;
        X_10_p4=X_10_p2*X_10_p2;
        X_01_p4=X_01_p2*X_01_p2;
        X_11_p2=X_11*X_11;
        X_20_p2=X_20*X_20;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

        F_3X=(3.0*B*(-X*(18.0 + 54.0*X_p6*B_p2 + X_p2*(26.0 - 72.0*B*SQ1X2) + 
                         X_p4*(11.0 -54.0*B*SQ1X2 + 54.0*B_p2)) + 6.0*(-SQ1X2 - 
                                                                       2.0*X_p2*SQ1X2 + 4.0*X_p6*B + X_p4*(-SQ1X2 + 7.0*B + 9.0*SQ1X2*B_p2)) 
                     *ASINH_X + 9.0*X_p5*(4.0 + X_p2)*B_p2*ASINH_X_p2)) / (SU1X2_p5f2*ComDenom_p4);

        F_4X=(3.0*B*(X*(20.0 + 7.0*X_p2 + 2.0*X_p4) + 3.0*(48.0 + 148.0*X_p2 + 77.0*X_p4 + 22.0*X_p6)*B*ASINH_X + 
                     81.0*X_p6*B_p3*ASINH_X*(12.0*SU1X2 + (2.0*X_p2 -3.0)*ASINH_X_p2) + 
                     9.0*B_p2*(-12.0*X_p3*SU1X2*(8.0 + 5.0*X_p2) - X*(48.0 + 160.0*X_p2 + 227.0*X_p4 + 70.0*X_p6)*ASINH_X_p2) + 
                     6.0*SQ1X2*B*(X*(36.0 + X_p2*(32.0 + X_p2*(11.0 + 108.0*SU1X2*B_p2))) + 
                                  3.0*B*ASINH_X*(-2.0*X_p2*(18.0 + 58.0*X_p2 + 25.0*X_p4) + 3.0*ASINH_X * 
                                                 (X_p5*(-4.0 + 11.0*X_p2)*B -4.0*SU1X2_p3*ASINH_X))))) / 
            (SU1X2_p7f2*ComDenom_p5);

/* Derivatives expressions */
        ds->df10000 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

        ds->df00100 += -B*factor*(ro_p4f3*F_1X*X_01);

        ds->df20000 += -B*factor*( ((4.0/9.0)*F/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_1X*X_10 + 
                                   ro_p4f3*F_2X*X_10_p2 + ro_p4f3*F_1X*X_20 );

        ds->df00200 += -B*factor*( ro_p4f3*F_2X*X_01_p2 );

        ds->df10100 += -B*factor*( (4.0/3.0)*ro_p1f3*F_1X*X_01 + ro_p4f3*F_2X*X_01*X_10 + 
                                   ro_p4f3*F_1X*X_11 );
 
 
        ds->df30000 += -B*factor*( -((8.0/27.0)*F/ro_p5f3) + ((4.0/3.0)*F_1X*X_10/ro_p2f3) + 
                                   4.0*ro_p1f3*F_2X*X_10_p2 + ro_p4f3*F_3X*X_10_p3 + 4.0*ro_p1f3*F_1X*X_20 + 
                                   3.0*ro_p4f3*F_2X*X_10*X_20 + ro_p4f3*F_1X*X_30 );
		
        ds->df20100 += -B*factor*( ((4.0/9.0)*F_1X*X_01/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_2X*X_01*X_10 + 
                                   ro_p4f3*F_3X*X_01*X_10_p2 + (8.0/3.0)*ro_p1f3*F_1X*X_11 + 
                                   2.0*ro_p4f3*F_2X*X_10*X_11 + ro_p4f3*F_2X*X_01*X_20 + 
                                   ro_p4f3*F_1X*X_21 );

        ds->df10200 += -B*factor*( (4.0/3.0)*ro_p1f3*F_2X*X_01_p2 +  
                                   ro_p4f3*F_3X*X_01_p2*X_10 +  
                                   2*ro_p4f3*F_2X*X_01*X_11  );

        ds->df00300 += -B*factor*( ro_p4f3*F_3X*X_01_p3 );

        ds->df40000 += -B*factor*( 
                                   ((40.0/81.0)*F/ro_p8f3) - ((32.0/27.0)*F_1X*X_10/ro_p5f3) + 
                                   ((8.0/3.0)*F_2X*X_10_p2/ro_p2f3) + (16.0/3.0)*ro_p1f3*F_3X*X_10_p3 + 
                                   ro_p4f3*F_4X*X_10_p4 + ((8.0/3.0)*F_1X*X_20/ro_p2f3) + 
                                   16.0*ro_p1f3*F_2X*X_10*X_20 + 6.0*ro_p4f3*F_3X*X_10_p2*X_20 + 
                                   3.0*ro_p4f3*F_2X*X_20_p2 + (16.0/3.0)*ro_p1f3*F_1X*X_30 + 
                                   4.0*ro_p4f3*F_2X*X_10*X_30 + ro_p4f3*F_1X*X_40 );

        ds->df30100 += -B*factor*( 
                                   -((8.0/27.0)*F*X_01/ro_p5f3) + ((4.0/3.0)*F_2X*X_01*X_10/ro_p2f3) + 
                                   4.0*ro_p1f3*F_3X*X_01*X_10_p2 + ro_p4f3*F_4X*X_01*X_10_p3 + 
                                   ((4.0/3.0)*F_1X*X_11/ro_p2f3) + 8.0*ro_p1f3*X_10*X_11 + 
                                   3.0*ro_p4f3*F_3X*X_10_p2*X_11 + 4.0*ro_p1f3*F_2X*X_01*X_20 + 
                                   3.0*ro_p4f3*F_3X*X_01*X_10*X_20 + 3.0*ro_p4f3*F_2X*X_11*X_20 + 
                                   4.0*ro_p1f3*F_1X*X_21 + 3.0*ro_p4f3*F_2X*X_10*X_21 + 
                                   3.0*ro_p4f3*F_2X*X_10*X_21 + ro_p4f3*F_2X*X_01*X_30 + 
                                   ro_p4f3*F_1X*X_31 );

        ds->df20200 += -B*factor*( 
                                   ((4.0/9.0)*F_2X*X_01_p2/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_3X*X_01_p2*X_10 + 
                                   ro_p4f3*F_4X*X_01_p2*X_10_p2 + (16.0/3.0)*ro_p1f3*F_2X*X_01*X_11 + 
                                   4.0*ro_p4f3*F_3X*X_01*X_10*X_11 + 2.0*ro_p4f3*F_2X*X_11_p2 + 
                                   ro_p4f3*F_3X*X_01_p2*X_20 + 2.0*ro_p4f3*F_2X*X_01*X_21 );

        ds->df10300 += -B*factor*( 
                                   (4.0/3.0)*ro_p1f3*F_3X*X_01_p3 + ro_p4f3*F_4X*X_01_p3*X_10 + 
                                   3.0*ro_p4f3*F_3X*X_01_p2*X_11 );

        ds->df00400 += -B*factor*( 
                                   ro_p4f3*F_4X*X_01_p4 );
    }
  
    if (dp->rhob >LB94_THRESHOLD) {    

/* Setting up ro(density) and gro(gradient of density) 
 * This is done since the code reuses variables for alpha and beta 
 * densities  */
        ro = dp->rhob;
        gro = dp->gradb;

/* Powers of ro */
        ro_p1f3=pow(ro,1.0/3.0);
        ro_p2f3=ro_p1f3*ro_p1f3;
        ro_p4f3=ro*ro_p1f3;
        ro_p5f3=ro*ro_p2f3;
        ro_p2=ro*ro;
        ro_p7f3=ro_p2*ro_p1f3;
        ro_p8f3=ro_p7f3*ro_p1f3;
        ro_p10f3=ro_p7f3*ro;
        ro_p13f3=ro_p10f3*ro;
        ro_p16f3=ro_p13f3*ro;

/* X auxaliary function and its powers */
        X=gro/ro_p4f3;

        X_p2=X*X;
        X_p3=X_p2*X;
        X_p4=X_p3*X;
        X_p5=X_p4*X;
        X_p6=X_p5*X;

/* (1+ X*X) and its powers (often used term) */
        SU1X2=1.0 + X_p2;
        SQ1X2=sqrt(SU1X2);
        SU1X2_p3f2=SU1X2*SQ1X2;
        SU1X2_p5f2=SU1X2_p3f2*SU1X2;
        SU1X2_p3 = SU1X2*SU1X2*SU1X2;
        SU1X2_p7f2=SU1X2_p5f2*SU1X2;

/* asinh(X) and its powers */
        ASINH_X=asinh(X);
        ASINH_X_p2=ASINH_X*ASINH_X;

/* ComDenom - frequently used Common Denominator from LB94 energy expression
 * and its powers */
        ComDenom=1.0+ 3.0*B*X*ASINH_X;
        ComDenom_p2=ComDenom*ComDenom;
        ComDenom_p3=ComDenom_p2*ComDenom;
        ComDenom_p4=ComDenom_p3*ComDenom;
        ComDenom_p5=ComDenom_p4*ComDenom;

/* Derivatives of X and powers of derivatives */
        X_10=-(4.0/3.0)*gro/ro_p7f3;
        X_20=(28.0/9.0)*gro/ro_p10f3;
        X_30=-(280.0/27.0)*gro/ro_p13f3;
        X_40=(3640.0/81.0)*gro/ro_p16f3;
        X_01=1.0/ro_p4f3;
        X_11=-(4.0/3.0)/ro_p7f3;
        X_21=(28.0/9.0)/ro_p10f3;
        X_31=-(280.0/27.0)/ro_p13f3;

        X_10_p2=X_10*X_10;
        X_01_p2=X_01*X_01;
        X_10_p3=X_10_p2*X_10;
        X_01_p3=X_01_p2*X_01;
        X_10_p4=X_10_p2*X_10_p2;
        X_01_p4=X_01_p2*X_01_p2;
        X_11_p2=X_11*X_11;
        X_20_p2=X_20*X_20;

/* F(X) auxaliary function */
        F=X_p2/ComDenom;

/* F(X) derivatives with respect to X */
        F_1X=X*(2.0 - (3.0*X_p2*B/SQ1X2) + (ComDenom - 1.0) )/ ComDenom_p2;

        F_2X=(2.0*SQ1X2 + 2*X_p2*(SQ1X2 -9.0*B) + 3*X_p4*B*(-5.0 +6.0*B*SQ1X2) - 
              9.0*X_p3*(2+X_p2)*B_p2*ASINH_X ) / (SU1X2_p3f2*ComDenom_p3);

        F_3X=(3.0*B*(-X*(18.0 + 54.0*X_p6*B_p2 + X_p2*(26.0 - 72.0*B*SQ1X2) + 
                         X_p4*(11.0 -54.0*B*SQ1X2 + 54.0*B_p2)) + 6.0*(-SQ1X2 - 
                                                                       2.0*X_p2*SQ1X2 + 4.0*X_p6*B + X_p4*(-SQ1X2 + 7.0*B + 9.0*SQ1X2*B_p2)) 
                     *ASINH_X + 9.0*X_p5*(4.0 + X_p2)*B_p2*ASINH_X_p2)) / (SU1X2_p5f2*ComDenom_p4);

        F_4X=(3.0*B*(X*(20.0 + 7.0*X_p2 + 2.0*X_p4) + 3.0*(48.0 + 148.0*X_p2 + 77.0*X_p4 + 22.0*X_p6)*B*ASINH_X + 
                     81.0*X_p6*B_p3*ASINH_X*(12.0*SU1X2 + (2.0*X_p2 -3.0)*ASINH_X_p2) + 
                     9.0*B_p2*(-12.0*X_p3*SU1X2*(8.0 + 5.0*X_p2) - X*(48.0 + 160.0*X_p2 + 227.0*X_p4 + 70.0*X_p6)*ASINH_X_p2) + 
                     6.0*SQ1X2*B*(X*(36.0 + X_p2*(32.0 + X_p2*(11.0 + 108.0*SU1X2*B_p2))) + 
                                  3.0*B*ASINH_X*(-2.0*X_p2*(18.0 + 58.0*X_p2 + 25.0*X_p4) + 3.0*ASINH_X * 
                                                 (X_p5*(-4.0 + 11.0*X_p2)*B -4.0*SU1X2_p3*ASINH_X))))) / 
            (SU1X2_p7f2*ComDenom_p5);

/* Derivatives expressions */
        ds->df01000 += -B*factor*((4.0/3.0)*ro_p1f3*F + ro_p4f3*F_1X*X_10);

        ds->df00010 += -B*factor*(ro_p4f3*F_1X*X_01);

        ds->df02000 += -B*factor*( ((4.0/9.0)*F/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_1X*X_10 + 
                                   ro_p4f3*F_2X*X_10_p2 + ro_p4f3*F_1X*X_20 );

        ds->df00020 += -B*factor*( ro_p4f3*F_2X*X_01_p2 );

        ds->df01010 += -B*factor*( (4.0/3.0)*ro_p1f3*F_1X*X_01 + ro_p4f3*F_2X*X_01*X_10 + 
                                   ro_p4f3*F_1X*X_11 );
 
 
        ds->df03000 += -B*factor*( -((8.0/27.0)*F/ro_p5f3) + ((4.0/3.0)*F_1X*X_10/ro_p2f3) + 
                                   4.0*ro_p1f3*F_2X*X_10_p2 + ro_p4f3*F_3X*X_10_p3 + 4.0*ro_p1f3*F_1X*X_20 + 
                                   3.0*ro_p4f3*F_2X*X_10*X_20 + ro_p4f3*F_1X*X_30 );
		
        ds->df02010 += -B*factor*( ((4.0/9.0)*F_1X*X_01/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_2X*X_01*X_10 + 
                                   ro_p4f3*F_3X*X_01*X_10_p2 + (8.0/3.0)*ro_p1f3*F_1X*X_11 + 
                                   2.0*ro_p4f3*F_2X*X_10*X_11 + ro_p4f3*F_2X*X_01*X_20 + 
                                   ro_p4f3*F_1X*X_21 );

        ds->df01020 += -B*factor*( (4.0/3.0)*ro_p1f3*F_2X*X_01_p2 +  
                                   ro_p4f3*F_3X*X_01_p2*X_10 +  
                                   2*ro_p4f3*F_2X*X_01*X_11  );

        ds->df00030 += -B*factor*( ro_p4f3*F_3X*X_01_p3 );

        ds->df04000 += -B*factor*( 
                                   ((40.0/81.0)*F/ro_p8f3) - ((32.0/27.0)*F_1X*X_10/ro_p5f3) + 
                                   ((8.0/3.0)*F_2X*X_10_p2/ro_p2f3) + (16.0/3.0)*ro_p1f3*F_3X*X_10_p3 + 
                                   ro_p4f3*F_4X*X_10_p4 + ((8.0/3.0)*F_1X*X_20/ro_p2f3) + 
                                   16.0*ro_p1f3*F_2X*X_10*X_20 + 6.0*ro_p4f3*F_3X*X_10_p2*X_20 + 
                                   3.0*ro_p4f3*F_2X*X_20_p2 + (16.0/3.0)*ro_p1f3*F_1X*X_30 + 
                                   4.0*ro_p4f3*F_2X*X_10*X_30 + ro_p4f3*F_1X*X_40 );

        ds->df03010 += -B*factor*( 
                                   -((8.0/27.0)*F*X_01/ro_p5f3) + ((4.0/3.0)*F_2X*X_01*X_10/ro_p2f3) + 
                                   4.0*ro_p1f3*F_3X*X_01*X_10_p2 + ro_p4f3*F_4X*X_01*X_10_p3 + 
                                   ((4.0/3.0)*F_1X*X_11/ro_p2f3) + 8.0*ro_p1f3*X_10*X_11 + 
                                   3.0*ro_p4f3*F_3X*X_10_p2*X_11 + 4.0*ro_p1f3*F_2X*X_01*X_20 + 
                                   3.0*ro_p4f3*F_3X*X_01*X_10*X_20 + 3.0*ro_p4f3*F_2X*X_11*X_20 + 
                                   4.0*ro_p1f3*F_1X*X_21 + 3.0*ro_p4f3*F_2X*X_10*X_21 + 
                                   3.0*ro_p4f3*F_2X*X_10*X_21 + ro_p4f3*F_2X*X_01*X_30 + 
                                   ro_p4f3*F_1X*X_31 );

        ds->df02020 += -B*factor*( 
                                   ((4.0/9.0)*F_2X*X_01_p2/ro_p2f3) + (8.0/3.0)*ro_p1f3*F_3X*X_01_p2*X_10 + 
                                   ro_p4f3*F_4X*X_01_p2*X_10_p2 + (16.0/3.0)*ro_p1f3*F_2X*X_01*X_11 + 
                                   4.0*ro_p4f3*F_3X*X_01*X_10*X_11 + 2.0*ro_p4f3*F_2X*X_11_p2 + 
                                   ro_p4f3*F_3X*X_01_p2*X_20 + 2.0*ro_p4f3*F_2X*X_01*X_21 );

        ds->df01030 += -B*factor*( 
                                   (4.0/3.0)*ro_p1f3*F_3X*X_01_p3 + ro_p4f3*F_4X*X_01_p3*X_10 + 
                                   3.0*ro_p4f3*F_3X*X_01_p2*X_11 );

        ds->df00040 += -B*factor*( ro_p4f3*F_4X*X_01_p4 );
	
    }
    SlaterFunctional.fourth(ds, factor, dp);
    VWNFunctional.fourth(ds,   factor, dp);
}
#endif
