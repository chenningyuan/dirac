/* fun-vwn.c:
   implementation of VWN functional and its derivatives 
   (c), Pawel Salek, pawsa@theochem.kth.se, sep 2001
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
#include "dcbham.h"

/* INTERFACE PART */
static int  vwn_read(const char* conf_line);
static real vwn_energy(real rhoa, real rhob, real grada, real gradb);
static void vwn_first(FirstDrv *ds, real factor, 
                        real rhoa, real rhob, real ngrada, real ngradb);
static void vwn_second(SecondFuncDrv *ds, real factor,
                         real rhoa, real rhob, real ngrada, real ngradb);
static void vwn_third(ThirdFuncDrv *ds, real factor,
                        real rhoa, real rhob, real ngrada, real ngradb);

Functional VWNFunctional = {
  "VWN",      /* name */
  fun_false,  /* not gga-corrected */
  vwn_read,   /* no extra input expected, just set the common block */
  NULL,
  vwn_energy, 
  vwn_first,
  vwn_second,
  vwn_third
};

/* IMPLEMENTATION PART */
static int
vwn_read(const char* conf_line)
{
  dcrham_.hfxfac = 0; 
  return 1;
}

static void
vwn_en_pot(real* enpot, real rho, int potentp)
{
  const real DCRS = pow(3.0/(4*M_PI),1.0/6.0),
    AI   =  0.0621814,
    BI   =  3.72744,
    CI   = 12.9352,
    X0I  = -0.10498,
    Q    = sqrt(4*CI - BI*BI),
    XF0I = X0I*X0I + BI*X0I + CI,
    YF0I = Q/(BI + 2*X0I),
    B    = X0I/XF0I,
    C    = XF0I*YF0I;
  real rho13 = pow(rho, 1.0/3.0);
  real x, xf, xfx, yf, yfx;
  /*
   *   See Table 5 in VWN paper
   *
   */
  x   = DCRS/sqrt(rho13);
  xf  = x*x + BI*x + CI;
  xfx = 2*x + BI;
  yf  = Q/xfx;
  yfx = -2*yf/xfx;
  enpot[0]  = 0.5*AI*(2*log(x) + log(xf)*(BI*B - 1)
                      +2*BI*(1.0/Q - X0I/C)*atan(yf)
                      -2*(BI*B)*log(x - X0I));
  if (potentp) {
    real xrho = -x/(6*rho);
    real ex1 = 0.5*AI*(2.0/x + xfx*(BI*B - 1)/xf
                       +2*BI*(1.0/Q-X0I/C)*yfx/(1+yf*yf)
                       -2*(BI*B)/(x-X0I));
    enpot[1] = enpot[0] + rho*ex1*xrho;
  }
}

static real
vwn_energy(real rhoa, real rhob, real grada, real gradb)
{
  real ep[2];
  vwn_en_pot(ep, rhoa+rhob, 0);
  return ep[0]*(rhoa+rhob);
}

static void
vwn_first(FirstDrv *ds, real factor, 
          real rhoa, real rhob, real ngrada, real ngradb)
{
  real ep[2];
  vwn_en_pot(ep, rhoa+rhob, 1);
  ds->df1000 += ep[1]*factor;
}

static void
vwn_second(SecondFuncDrv *ds, real factor,
           real rhoa, real rhob, real ngrada, real ngradb)
{
  static const real AI=0.0621814, BI=3.72744, CI=12.9352, X0I=-0.10498;
  const real
    Q    = sqrt(4.0*CI - BI*BI),
    XF0I = X0I*X0I + BI*X0I + CI,
    YF0I = Q/(BI + 2*X0I),
    DCRS = pow(3.0/(4*M_PI),1.0/6.0),
    B    = X0I/XF0I,
    C    = XF0I*YF0I,
    ACON = B*BI - 1.0,
    BCON = 2*ACON + 2,
    CCON = 2*BI*(1.0/Q - X0I/C);

  real rho = rhoa+rhob, rho13 = pow(rho, 1.0/3.0);
  real x    = DCRS/sqrt(rho13);
  real xf   = x*x + BI*x+CI;
  real xfx  = 2.0*x + BI;
  real yf   = Q/xfx;

  real e1   = 2*log(x) 
    + ACON*log(xf)
    - BCON*log(x - X0I)
    + CCON*atan(yf);
  real ex1  = 2.0/x 
    + ACON*xfx/xf
    - BCON/(x - X0I)
    - CCON*(2*yf/xfx)/(1 + yf*yf);
  real exx1 = 2/(x*x) 
    - ACON*(2.0/xf - xfx*xfx/(xf*xf))
    - BCON/((x - X0I)*(x - X0I))
    - CCON*8*Q*xfx/(pow(Q*Q + xfx*xfx,2.0));
  real xrho = -x/(6.0*rho);

  real drv2 = 0.5*AI*xrho*(2*ex1 - rho*(exx1*xrho + 7*ex1/(6.0*rho)))*factor;

  ds->df1000 += 0.5*AI*(e1 + rho*ex1*xrho)*factor;
  ds->df2000 += drv2;
  ds->df1100 += drv2;
}

static void
vwn_third(ThirdFuncDrv *ds, real factor,
          real rhoa, real rhob, real ngrada, real ngradb)
{
  static const real AI=0.0621814,
    BI   =  3.72744,
    CI   = 12.9352,
    X0I  = -0.10498;
  const real 
    Q    = sqrt(4*CI - BI*BI),
    XF0I = X0I*X0I + BI*X0I + CI,
    YF0I = Q/(BI + 2*X0I),
    DCRS = pow(3.0/(4*M_PI),1.0/6.0),
    B    = X0I/XF0I,
    C    = XF0I*YF0I,
    ACON = B*BI - 1.0,
    BCON = 2*ACON + 2.0,
    CCON = 2*BI*(1.0/Q - X0I/C);

  real rho   = rhoa+rhob;
  real rho13 = pow(rho,1.0/3.0);
  real x     = DCRS/sqrt(rho13);
  real xrho  = -DCRS*pow(rho,-7.0/6.0)/6.0;
  real xxrho = +DCRS*pow(rho,-13.0/6.0)*7.0/36.0;
  real xf   = x*x + BI*x+CI;
  real xfx  = 2*x + BI;
  real yf   = Q/xfx;

  real e1   = 2*log(x) 
    + ACON*log(xf)
    - BCON*log(x - X0I)
    + CCON*atan(yf);
  real ex1  = 2.0/x 
    + ACON*xfx/xf
    - BCON/(x - X0I)
    - CCON*(2*yf/xfx)/(1.0 + yf*yf);
  real exx1 = -2.0/(x*x) 
    + ACON*(2.0/xf - xfx*xfx/(xf*xf))
    + BCON/((x - X0I)*(x - X0I))
    + CCON*8*Q*xfx/((Q*Q + xfx*xfx)*(Q*Q + xfx*xfx));
  real exxx1 = 4.0/(x*x*x)
    + ACON*(2.0*xfx/(xf*xf))*(xfx*xfx/xf-3.0)
    - BCON*2.0/pow(x - X0I,3.0)
    + CCON*16.0*Q*(Q*Q-3.0*xfx*xfx)/pow(Q*Q + xfx*xfx,3.0);

  real vd0  = 0.5*AI*(e1 + rho*ex1*xrho);
  real vd1  = 0.5*AI*xrho*(2*ex1 + rho*exx1*xrho - ex1*7.0/6.0);
    
  real vd2  = 0.5*AI*xxrho*(2*ex1 + rho*exx1*xrho - 7.0/6.0*ex1)
          + 0.5*AI*xrho*(2*exx1*xrho 
                         +exx1*xrho + rho*(exxx1*xrho*xrho+exx1*xxrho)
                         -7.0/6.0*exx1*xrho);

  ds->df1000 += factor*vd0;
  ds->df2000 += factor*vd1;
  ds->df1100 += factor*vd1;
  ds->df3000 += factor*vd2;
  ds->df2100 += factor*vd2;
}

