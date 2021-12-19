/* dirac_copyright_start */
/*
 *
 *     Copyright (c) by the authors of DIRAC.
 *
 *     This program is free software; you can redistribute it and/or
 *     modify it under the terms of the GNU Lesser General Public
 *     License version 2.1 as published by the Free Software Foundation.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *     Lesser General Public License for more details.
 *
 *     If a copy of the GNU LGPL v2.1 was not distributed with this
 *     code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
 */
/* dirac_copyright_end */

/*
 *
 * LBalpha
 *
 * a minute modification of the LB94 potential for the SAOP asymptotic correction
 * so please have a look in fun-lb94.c
 *
 * the only difference are ALPHA (= 1.19) and a modified BETA (= 0.01)
 * also have a look at equation (2.1) in J. Chem. Phys. 112, 1344 (2000).
 *
 * implementation: Radovan Bast - last modification 21/04/2005
 *
 */

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/*
 *
 * interface
 *
 */

static int  lbalpha_isgga(void) { return 1; }
static int  lbalpha_read(const char* conf_line);
static real lbalpha_energy(const FunDensProp* dens_prop);
static void lbalpha_first(FunFirstFuncDrv *ds, real factor, 
                          const FunDensProp* dens_prop);
static void lbalpha_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dens_prop);
static void lbalpha_third(FunThirdFuncDrv *ds, real factor,
                          const FunDensProp* dens_prop);
#ifdef FOURTH_ORDER_DERIVATIVES
static void lbalpha_fourth(FunFourthFuncDrv *ds, real factor,
                          const FunDensProp* dens_prop);
#endif

Functional LBalphaFunctional = {"LBalpha",   /* name */
                             lbalpha_isgga,  /* gga-corrected */
                             lbalpha_read,   /* set common blocks */
                             NULL,           /* reporter */
                             lbalpha_energy, 
                             lbalpha_first,
                             lbalpha_second,
                             lbalpha_third
#ifdef FOURTH_ORDER_DERIVATIVES
                            ,lbalpha_fourth
#endif
};

/*
 *
 * implementation
 *
 */

static int
lbalpha_read(const char* conf_line)
{
  fun_set_hf_weight(0.0);
  return 1;
}

static const real LBalpha_THRESHOLD = 1e-14;
static const real ALPHA = 1.19;
static const real BETA = 0.01;

static real
lbalpha_energy(const FunDensProp* dp)
{
  return SlaterFunctional.func(dp)+VWNFunctional.func(dp);
}

static void
lbalpha_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
  real vx;

/*real rho    = dp->rhoa + dp->rhob;*/
  real rho    = dp->rhoa;
  real rho13  = pow(rho, 1.0/3.0);
/*real grad   = dp->grada + dp->gradb;*/
  real grad   = dp->grada;
  real rho43  = rho*rho13;
  real scaled_grad, sg2;
  scaled_grad = grad/(rho43>1e-13 ? rho43 : 1e-13);
  sg2         = scaled_grad*scaled_grad;

  vx = -BETA*rho13*sg2/
      (1.0+3.0*BETA*scaled_grad*asinh(scaled_grad));

  ds->df1000 += vx*factor;
  ds->df0100 += vx*factor;

  SlaterFunctional.first(ds, ALPHA*factor, dp);
  VWNFunctional.first(ds, factor, dp);
}

void
lbalphapot_(real *ds, const real* weight, const real* rho, const real* grad)
{
    FunFirstFuncDrv drvs;
    FunDensProp dp;

    /* HP does not grok C99's { *rhoa, *rhob, *grada, *gradb }; */
    dp.rhoa   = dp.rhob  = *rho *0.5;
    dp.grada  = dp.gradb = *grad *0.5;
    dp.gradab = dp.grada * dp.gradb;
    dp.subsystem[0] = NULL; dp.subsystem[1] = NULL;

    drv1_clear(&drvs);
    lbalpha_first(&drvs, *weight, &dp);

    ds[0] = drvs.df1000;
    ds[1] = drvs.df0010 + 0.5*drvs.df00001* *grad;
}

static void
lbalpha_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
/*
 *
 * assuming that LR and higher are eqivalent to ALDA like in LB94 ...
 *
 */
  SlaterFunctional.second(ds, factor, dp);
  VWNFunctional.second(ds, factor, dp);
}

static void
lbalpha_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
  SlaterFunctional.third(ds, factor, dp);
  VWNFunctional.third(ds, factor, dp);
}

#ifdef FOURTH_ORDER_DERIVATIVES  
static void
lbalpha_fourth(FourthFuncDrv *ds, real factor, const FunDensProp* dp)
{
  SlaterFunctional.fourth(ds, factor, dp);
  VWNFunctional.fourth(ds, factor, dp);
}
#endif
