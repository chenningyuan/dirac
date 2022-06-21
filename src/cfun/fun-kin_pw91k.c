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

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
#define LOG log
#define EXP exp
#define ABS fabs
#define ASINH asinh
#define SQRT sqrt

/* INTERFACE PART */
static int pw91k_read(const char* conf_line);
static real pw91k_energy(const FunDensProp* dp);

Functional PW91k_KinFunctional = {
  "PW91k",
  fun_true,
  pw91k_read,
  NULL,
  pw91k_energy,
  NULL,
  NULL,
  NULL
};

/* IMPLEMENTATION PART */
static int
pw91k_read(const char* conf_line)
{
    fun_set_hf_weight(0);
    return 1;
}


static real
pw91k_energy(const FunDensProp* dp)
{
    real zk;
    real rhoa  = dp->rhoa;
    real rhob  = dp->rhob;
    real grada = dp->grada;
    real gradb = dp->gradb;

    real A1    =   0.093907;
    real A2    =   0.26608;
    real A3    =   0.0809615;
    real A4    = 100.00;
    real A     =  76.320;
    real B1    =   0.57767e-4;

    real Cf    = (3.0/10.0)*pow((3.0/(M_PI*M_PI)),2.0/3.0);

    real kfa = pow(6*M_PI*M_PI*rhoa,1.0/3.0); 
    real kfb = pow(6*M_PI*M_PI*rhob,1.0/3.0);

    real ya  = ABS(grada)/(2*kfa*rhoa);
    real yb  = ABS(grada)/(2*kfb*rhob);

    real Fya = ( 1 + A1*ya*ASINH(A*ya) + ya*ya*(A2 - A3*EXP(-A4*ya*ya)) ) / ( 1 + A1*ya*ASINH(A*ya) + B1*ya*ya*ya*ya );
    real Fyb = ( 1 + A1*yb*ASINH(A*yb) + yb*yb*(A2 - A3*EXP(-A4*yb*yb)) ) / ( 1 + A1*yb*ASINH(A*yb) + B1*yb*yb*yb*yb );

    zk = pow(2,2.0/3.0) * Cf * ( pow(rhoa,5.0/3.0)*Fya + pow(rhob,5.0/3.0)*Fyb );

    return zk;
}

