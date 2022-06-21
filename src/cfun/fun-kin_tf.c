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

/* Automatically generated functional code: kin_tf
   Maxima input:
    >> PI: 3.14159265358979312;
    >> ctf: (3/10)*(3*PI^2)^(2/3);
    >> 
    >> rho(rhoa,rhob):= rhoa+rhob;
    >> 
    >> // Thomas-Fermi kinetic energy functional 
    >> tf(rhoa,rhob):= ctf*rho(rhoa,rhob)^(5/3);
    >> 
    >> K(rhoa,grada,rhob,gradb,gradab):= tf(rhoa,rhob);
*/

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"

/* INTERFACE PART */
static int kin_tf_read(const char* conf_line);
static real kin_tf_energy(const FunDensProp* dp);
static void kin_tf_first(FunFirstFuncDrv *ds, real factor, 
                         const FunDensProp* dp);
static void kin_tf_second(FunSecondFuncDrv *ds, real factor,
                          const FunDensProp* dp);
static void kin_tf_third(FunThirdFuncDrv *ds, real factor,
                         const FunDensProp* dp);

Functional TF_KinFunctional = {
  "kin_tf",
  fun_false,
  kin_tf_read,
  NULL,
  kin_tf_energy,
  kin_tf_first,
  kin_tf_second,
  kin_tf_third
};

/* IMPLEMENTATION PART */
static int
kin_tf_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}

static real
kin_tf_energy(const FunDensProp* dp)
{
    return 2.871234000188191 * pow(dp->rhob+dp->rhoa, 5.0/3.0);
}

static void
kin_tf_first(FunFirstFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t;
    
    t = 4.785390000313651 * pow(dp->rhob+dp->rhoa, 2.0/3.0);
    
    ds->df1000 += factor * t;
    ds->df0100 += factor * t;

}

static void
kin_tf_second(FunSecondFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t1, t2, t3;

    t1 = dp->rhob + dp->rhoa;
    t2 = 4.785390000313651 * pow(t1, 2.0/3.0);
    t3 = 3.190260000209101 / pow(t1, 1.0/3.0);

    ds->df1000 += factor*t2;
    ds->df0100 += factor*t2;
    ds->df2000 += factor*t3;
    ds->df1100 += factor*t3;
    ds->df0200 += factor*t3;
}

static void
kin_tf_third(FunThirdFuncDrv *ds, real factor, const FunDensProp* dp)
{
    real t1, t2, t3, t4;
 
    t1 = dp->rhob + dp->rhoa;
    t2 = 4.785390000313651 * pow(t1, 2.0/3.0);
    t3 = 3.190260000209101 / pow(t1, 1.0/3.0);
    t4 = -1.0634200000697  / pow(t1, 4.0/3.0);

    ds->df1000 += factor * t2;
    ds->df0100 += factor * t2;
    ds->df2000 += factor * t3;
    ds->df1100 += factor * t3;
    ds->df0200 += factor * t3;
    ds->df3000 += factor * t4; 
    ds->df2100 += factor * t4; 
    ds->df1200 += factor * t4; 
    ds->df0300 += factor * t4;

}
