#if defined(__CVERSION__) /* THE C VERSION OF THE HEADERS */
/*-*-mode: C; c-indentation-style: "bsd"; c-basic-offset: 4; -*-*/
/* The variables, structures and functions related to computation
   of functional and their derivatives.
   (c) Pawel Salek, pawsa@theochem.kth.se. 2001.07.13

   Extension for Embedding
   Written by Christoph Jacob, 2004

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

/* DftDensProp structure contains properties of the density that are
   needed for functional evaluation and possibly other purposes.
*/
typedef struct DftDensProp_ DftDensProp;
struct DftDensProp_ {
  /* These variables are used in "normal DFT" for the density.
     When do embedding, they are used for the TOTAL rho          */

    real rhoa,  rhob;
    real grada, gradb; /* norms of the density gradient, not squares */
    real gradab;       /* scalar product of grada and gradb */

    /* real current[3] or something may come in the future :-) */

  /* This is for the densities of the (two) subsystems when doing embedding.
     In "normal" calculations it is not used and should be NULL) */
    DftDensProp *subsystem[2];
};

typedef struct Functional_ Functional;

typedef struct FirstDrv_  FirstDrv;
typedef struct SecondDrv_ SecondDrv;
typedef struct ThirdDrv_  ThirdDrv;

typedef struct FirstFuncDrv_  FirstFuncDrv;
typedef struct SecondFuncDrv_ SecondFuncDrv;
typedef struct ThirdFuncDrv_  ThirdFuncDrv;

/* FirstDrv: matrix of first order derivatives with respect to two
 * parameters: density rho and SQUARE of the gradient of density grho.  
 * zeta_i = |\nabla\rho_i|²
 * mu     = |\nabla\rho_\alpha||\nabla\rho_\beta|
 */

/* Comment on use of FirstDrv_->fZ (by Christoph Jacob) :
 * It says here, that fZ sould be d/dzeta F (zeta = |\nabla\rho|^2)
 * (which I think would be a clean way of programming)
 *
 * However, what is actually done in function _dftpot0 (file functionals.c)
 * is something different:
 * fZ is NOT d/dzeta F, but fZ = (d/dzeta F) * | \nabla\rho |
 * When the result is used in DFTDRV (file dftdrv.F) this is corrected
 * by dividing by | \nabla\rho | again (in addition, we start messing
 * around with a factor 2)
 *
 * The code is working and I dont know where else these values are
 * used and what I have seen, the same problem applys in SecondDrv_
 * and so on, so I wont try to fix it.
 *
 * But it took me some time to figure this out and I got very confused,
 * so maybe these lines are helpful when someone else looks at the code!
 */ 

struct FirstDrv_ {
    real fR;  /* d/drho F     */
    real fZ;  /* d/zeta F     */
};

struct FirstFuncDrv_ {
    real df1000;  /* d/drho F     */
    real df0100;
    real df0010;  /* d/zeta F     */
    real df0001;
    real df00001;
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
    real df1000;  /* d/drho_alpha F               */
    real df0100;  /* d/drho_beta F                */
    real df0010;  /* d/|zeta_alpha| F             */
    real df0001;  /* d/|zeta_beta| F              */
    real df00001;
    real df2000;  /* d/drho_alpha^2 F             */
    real df0200;  /* d/drho_beta^2 F              */
    real df0020;  /* d/dzeta_alpha^2 F          */
    real df0002;  /* d/dzeta_beta^2 F             */
    real df00002;
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df0011;  /* d2/dzeta_a zeta_b F          */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df10001;
    real df0110;  /* d/drho_beta d/dzeta_alpha F  */ 
    real df0101;  /* d/drho_beta d/dzeta_beta F   */
    real df01001;
  /* df00101 and df00011 not included */ 
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
    real df1000;   /* d/drho F          */
    real df0100; 
    real df0010;   /* d/|zeta| F        */
    real df0001;
    real df00001;
     
    real df2000;  /* d/drho_alpha^2 F  */
    real df0200; 
    real df1100;  /* d/(drho_alpha drho_beta) F   */
    real df1010;  /* d/drho_alpha d/dzeta_alpha F */
    real df1001;  /* d/drho_alpha d/dzeta_beta F  */
    real df0101;
    real df0110;
    real df10001;
    real df01001;   
    real df0020;  /* d/dzeta^2 F                  */
    real df0002;
    real df0011;  /* d2/dzeta_a zeta_b F          */

    real df3000; /* d/drho_alpha^3 F              */
    real df0300; 
    real df2100; /* d/d(rho_alpha^2rho_beta) F    */
    real df1200;
    real df2010; /* d/drhoa^2 d/dzetaa F */
    real df2001;
    real df0201; 
    real df0210;     
    real df1110; /* d/drho^2 d/dzeta F */
    real df1101;
    real df1020; /* d/drho d/dzeta^2 F */
    real df1002;
    real df0102;
    real df0120; 
    real df0030; /* d/dzeta^3 F        */
    real df0003;
    real df1011;
    real df0111;
    real df0021;
    real df0012;
    real df20001;
    real df02001;
    real df11001;  
  /* most derivatives for fith variable not included *
   * this makes functionals valid only to square grad rho 
   */
};

/* helper functions for FuncDrv structures */ 
void drv1_clear(FirstFuncDrv* gga);  /* set all components to 0 */
void drv2_clear(SecondFuncDrv* gga); /* set all components to 0 */
void drv3_clear(ThirdFuncDrv* gga);  /* set all components to 0 */

/* EnergyFunc: the function returning the energy for given densities
   and gradients. Note that some functionals(like LYP) depend explicitely
   on separately alpha and beta densities
*/
typedef int (*IsGGAFunc)(void);
typedef int (*ReadInputFunc)(const char* conf_string);
typedef void (*ReportFunc)(void);
typedef real (*EnergyFunc)(DftDensProp* dens_prop);
typedef void (*FirstOrderFun)(FirstFuncDrv *ds, real factor,
                              DftDensProp* dns_prp);

typedef void (*SecondOrderFun)(SecondFuncDrv *ds, real factor,
                               DftDensProp* dens_prop);

typedef void (*ThirdOrderFun)(ThirdFuncDrv *ds, real factor,
                              DftDensProp* dens_prop);

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

/* these are the functions that will be called from Fortran */
int  dftisgga_(void);
int  dftisemb_(void);
void dftreport_(void);
void dftlistfuncs_(void);
real dftenergy_(const real* rho, const real* grad);
void dftpot0_(FirstDrv *ds, const real* w, const real* rho, const real* grad);
void dftpot1_(SecondDrv*ds, const real* w, const real* rho, const real* grad,
              const int* triplet);
void dftpot2_(ThirdDrv *ds, real factor,  real rho, real grad, int isgga);

/* helper functions */
void fort_print(const char* format, ...);
int fun_false(void); /* always returns 0 */
int fun_true(void);  /* always returns 1 */

/* The list of functionals */
/* generic functionals */
extern Functional SlaterFunctional;
extern Functional VWNFunctional;
extern Functional B88CorrFunctional;
extern Functional BeckeFunctional;
extern Functional LB94Functional;
extern Functional LYPFunctional;
extern Functional PBEFunctional;
extern Functional PW86xFunctional;
extern Functional P86cFunctional;
extern Functional PW92Functional;

/* mixed functionals */
extern Functional XAlphaFunctional;
extern Functional RhoFunctional;
extern Functional LDAFunctional;
extern Functional GGAKeyFunctional;
extern Functional BLYPFunctional;
extern Functional B3LYPFunctional;
extern Functional WAHFunctional;

/* generic functionals for embedding */
extern Functional EmbeddingFunctional;
extern Functional EmbInteractionFunctional;

/* kinetic energy functionals */
extern Functional TF_KinFunctional;     /* Thomas-Fermi */
extern Functional VW_KinFunctional;     /* TF + Von Weizaecker */

/* the list of the functionals */
extern Functional* available_functionals[];       /* XC-functionlas */
extern Functional* available_kin_functionals[];   /* kinetic energy functionals */

#endif /* _FUNCTIONALS_H_ */
#else /* THE FORTRAN VERSION OF THE HEADERS */

      LOGICAL  DFTISGGA
      EXTERNAL DFTISGGA, DFTINPUT
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
