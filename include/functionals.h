#if defined(__CVERSION__) /* THE C VERSION OF THE HEADERS */
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* The variables, structures and functions related to computation
   of functional and their derivatives.
   (c) Pawel Salek, pawsa@theochem.kth.se. 2001.07.13

   NOTE1: the derivatives are computed with respect to the density,
   and SQUARE of the density gradient. This is a choice. It stems from
   the fact that the factors involved in the derivative vector
   distribution depend on the square of the density gradient.

   NOTE2: C version is included once per file, Fortran version -
   multiple times.
*/
#ifndef _FUNCTIONALS_H_
#define _FUNCTIONALS_H_

/* define the basic floating-point variable type used by the Fortran code */
typedef double real;

typedef struct Functional_ Functional;
extern Functional* selected_func;
typedef struct FirstDrv_  FirstDrv;
typedef struct SecondDrv_ SecondDrv;
typedef struct ThirdDrv_  ThirdDrv;

typedef struct SecondFuncDrv_ SecondFuncDrv;
typedef struct ThirdFuncDrv_  ThirdFuncDrv;

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.  
 * zeta_i = |\nabla\rho_i|²
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */
struct FirstDrv_ {
    real df1000;  /* d/drho F     */
    real df0010;  /* d/zeta F     */
};

/* SecondDrv:  matrix  of  second  order functional  derivatives  with
 * respect  to two  parameters: density  rho and  SQUARE  of the
 * density gradient zeta.  The derivatives are computed for alpha-alpha
 * or beta-beta spin-orbital block (i.e. include triplet flag).
 */
struct SecondDrv_ {
    real fR; /* d/drho  F */
    real fZ; /* d/dzeta F */
    real fRR; /* d/drho² F */
    real fRZ; /* d/(drho dzeta) F */
    real fZZ; /* d/dzeta² F */
};

/* SecondFuncDrv: this structure is used by functional derivative
 * evaluation procedures. Do not include "triplet" transformation.
 */
struct SecondFuncDrv_ {
    real df1000;  /* d/drho F          */
    real df0010;  /* d/|zeta| F        */
    real df2000;  /* d/drho_alpha^2 F  */
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df0020;  /* d/dzeta^2 F                  */
    real df0011;  /* d2/dzeta_a zeta_b F          */
};

/* ThirdDrv: matrix of third derivatives with respect to two parameters:
   density rho and SQUARE of the density gradient zeta.
*/
struct ThirdDrv_ {
    real fR;   /* d/drho  F */
    real fZ;   /* d/dzeta F */
    real fRR;  /* d/drho² F */
    real fRZ;  /* d/(drho dzeta) F */
    real fZZ;  /* d/dzeta² F */
    real fRRR; /* d/drho³ F */
    real fRRZ; /* d/(drho² dzeta) F */
    real fRZZ; /* d/(drho dzeta²) F */
    real fZZZ; /* d/dzeta³ F */
};

/* ThirdFuncDrv: matrix of third derivatives with respect to five
   parameters: density rho_alpha and SQUARE of the density gradient
   zeta.  and mu.
*/

struct ThirdFuncDrv_ {
    real df1000;  /* d/drho F          */
    real df0010;  /* d/|zeta| F        */
    real df2000;  /* d/drho_alpha^2 F  */
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df0020;  /* d/dzeta^2 F                  */
    real df0011;  /* d2/dzeta_a zeta_b F          */

    real df3000; /* d/drho_alpha^3 F              */
    real df2100; /* d/d(rho_alpha^2rho_beta) F    */
    real df2010; /* d/drhoa^2 d/dzetaa F */
    real df1110; /* d/drho^2 d/dzeta F */
    real df2001; /* d/drhoa^2 d/dzetab F */
    real df1020; /* d/drho d/dzeta^2 F */
    real df0120; /* d/drho d/dzeta^2 F */
    real df0030; /* d/dzeta^3 F        */
    real df1011;
    real df0021;
};

/* EnergyFunc: the function returning the energy for given densities
   and gradients. Note that some functionals(like LYP) depend explicitely
   on separately alpha and beta densities
*/
typedef int (*IsGGAFunc)(void);
typedef int (*ReadInputFunc)(const char* conf_string);
typedef void (*ReportFunc)(void);
typedef real (*EnergyFunc)(real rhoa, real rhob, real grada, real gradb);
typedef void (*FirstOrderFun)(FirstDrv *ds, real factor, 
                              real rhoa, real rhob, real ngrada, real ngradb);

typedef void (*SecondOrderFun)(SecondFuncDrv *ds, real factor,
                               real rhoa, real rhob, real ngrada, real ngradb);

typedef void (*ThirdOrderFun)(ThirdFuncDrv *ds, real factor,
                              real rhoa, real rhob, real ngrada, real ngradb);

struct Functional_ {
    const char* name; /* descriptive functional name (usually 5 characters) */
    IsGGAFunc       is_gga;
    ReadInputFunc   read;
    ReportFunc      report;
    EnergyFunc      func;
    FirstOrderFun   first;
    SecondOrderFun  second;
    ThirdOrderFun   third;
};

int  dft_is_gga(void);
void dftreport_(void);
void fort_print(const char* format, ...);
void dftlistfuncs_(void);
real dftenergy_(const real* rho, const real* grad);
void dftpot0_(FirstDrv *ds, const real* w, const real* rho, const real* grad);
void dftpot1_(SecondDrv*ds, const real* w, const real* rho, const real* grad,
              const int* triplet);
void dftpot2_(ThirdDrv *ds, real factor,  real rho, real grad, int isgga);

#define drv1_clear(gga) { (gga)->df1000 = (gga)->df0010 = 0; }
void drv2_clear(SecondFuncDrv* gga); /* set all components to 0 */
void drv3_clear(ThirdFuncDrv* gga);  /* set all components to 0 */
int fun_false(void); /* always returns 0 */
int fun_true(void);  /* always returns 1 */

/* The list of functionals */
/* generic functionals */
extern Functional SlaterFunctional;
extern Functional VWNFunctional;
extern Functional B88CorrFunctional;
extern Functional BeckeFunctional;
extern Functional LYPFunctional;
extern Functional PW86xFunctional;
extern Functional P86cFunctional;

/* mixed functionals */
extern Functional XAlphaFunctional;
extern Functional RhoFunctional;
extern Functional LDAFunctional;
extern Functional GGAKeyFunctional;
extern Functional BLYPFunctional;
extern Functional B3LYPFunctional;
extern Functional WAHFunctional;

/* the list of the functionals */
extern Functional* available_functionals[];

#endif /* _FUNCTIONALS_H_ */
#else /* THE FORTRAN VERSION OF THE HEADERS */

      LOGICAL  DFT_IS_GGA, DFT_INPUT
      EXTERNAL DFT_IS_GGA, DFT_INPUT
c     The derivatives with respect to rho and zeta
      INTEGER FR0, FZ0, FRR, FRZ, FZZ, FRRR, FRRZ, FRZZ, FZZZ
      INTEGER DERVS1, DERVS2, DERVS3
      PARAMETER(FR0=1,  FZ0=2,   FRR=3,  FRZ=4)
      PARAMETER(FZZ=5,  FRRR=6,  FRRZ=7, FRZZ=8, FZZZ=9)
      PARAMETER(DERVS1=3, DERVS2=5, DERVS3=9)

c     The derivatives are a vector that can be declared as
c     DIMENSION DRVS(DERVS3)
c     where 3 stands for the functional derivative order.
c     Second order derivatives can be declared as
c     DIMENSION DRVS(DERVS2)
c     etc.
c     The derivatives can be accessed as DRVS(F1000), DRVS(F0010) etc.
c
      INTEGER F1000, F0010, F2000, F1100, F1010, F1001, F0020
      INTEGER F3000, F2100, F2010, F1110, F2001
      INTEGER F1020, F0120, F0030, F1011, F0021
      INTEGER FUN_DERVS1, FUN_DERVS2, FUN_DERVS3
      PARAMETER(F1000=1, F0010=2, F2000=3, F1100=4, F1010=5, F1001=6)
      PARAMETER(F0020=7)
      PARAMETER(F3000=8, F2100=9, F2010=10,F1110=11,F2001=12)
      PARAMETER(F1020=13,F0120=14,F0030=15,F1011=16,F0021=17)
      PARAMETER(FUN_DERVS1=2, FUN_DERVS2=7, FUN_DERVS3=17)
#endif /* __CVERSION__ */
