/* Automatically generated functional code: pbex
   Maxima input:
    >> pi:3.14159265358979312;
    >> 
    >> // parameters for pbex 
    >> R:0.804;
    >> d:0.066725;
    >> mu:d*pi^2/3;
    >> Sa:xa/(2*(6*pi^2)^(1/3));
    >> Sb:xb/(2*(6*pi^2)^(1/3));
    >> 
    >> // functions for pbex 
    >> F(S):=1+R-R/(1+(mu*S^2)/R);
    >> Ea(n):=-3/4*(3/pi)^(1/3)*n^(4/3)*F(Sa);
    >> Eb(n):=-3/4*(3/pi)^(1/3)*n^(4/3)*F(Sb);
    >> 
    >> // kernel 
    >> K(rhoa,grada,rhob,gradb):=0.5*(Ea(2*rhoa)+Eb(2*rhob));
    >> 
    >> 
*/

// add "extern Functional pbexFunctional;" to 'functionals.h'
// add "&pbexFunctional," to 'functionals.c'

#include <math.h>
#include <stddef.h>

#define __CVERSION__

#include "functionals.h"
#include "dcbham.h"

#define LOG log
#define ABS fabs
#define ASINH asinh
#define SQRT sqrt

/* INTERFACE PART */
static int pbex_read(const char* conf_line);
//static real pbex_energy(real rhoa, real rhob, real grada, real gradb);
//static void pbex_first(FirstDrv *ds, real factor, 
//                        real rhoa, real rhob, real ngrada, real ngradb);
//static void pbex_second(SecondFuncDrv *ds, real factor,
//                        real rhoa, real rhob, real ngrada, real ngradb);
//static void pbex_third(ThirdFuncDrv *ds, real factor,
//                        real rhoa, real rhob, real ngrada, real ngradb);
static real pbex_energy(DftDensProp* dp);
static void pbex_first(FirstFuncDrv *ds, real factor, DftDensProp* dp);
static void pbex_second(SecondFuncDrv *ds, real factor,DftDensProp* dp);
static void pbex_third(ThirdFuncDrv *ds, real factor, DftDensProp* dp);


Functional pbexFunctional = {
  "pbex",
  fun_true,
  pbex_read,
  NULL,
  pbex_energy,
  pbex_first,
  pbex_second,
  pbex_third
};

/* IMPLEMENTATION PART */
static int
pbex_read(const char* conf_line)
{
  dft_set_hf_weight(0.0);
  return 1;
}


static real
pbex_energy(DftDensProp* dp)
{
    real t[2],zk;
    real rhoa,rhob,grada,gradb;

    rhoa=dp->rhoa;
    rhob=dp->rhob;
    grada=dp->grada;
    gradb=dp->gradb;

    t[1] = pow(2.0E0,3.3333333333333333333E-1);
    zk = 5.0E-1*(-1.4771175327640446419E0*t[1]*pow(rhob,1.3333333333333333333E0)*(1.804E0-8.04E-1/(4.492799480231090353E-3*pow(gradb,2.0E0)/pow(rhob,2.6666666666666666667E0)+1.0E0))-1.4771175327640446419E0*t[1]*pow(rhoa,1.3333333333333333333E0)*(1.804E0-8.04E-1/(4.492799480231090353E-3*pow(grada,2.0E0)/pow(rhoa,2.6666666666666666667E0)+1.0E0)));
    return zk;
}

static void
pbex_first(FirstFuncDrv *ds, real factor, DftDensProp* dp)
{
    real t[5];
    real dfdra, dfdga;
    real rhoa,rhob,grada,gradb;

    rhoa=dp->rhoa;
    rhob=dp->rhob;
    grada=dp->grada;
    gradb=dp->gradb;

    t[1] = pow(2.0E0,3.3333333333333333333E-1);
    t[2] = pow(grada,2.0E0);
    t[3] = 4.492799480231090353E-3*t[2]/pow(rhoa,2.6666666666666666667E0)+1.0E0;
    t[4] = 1/pow(t[3],2.0E0);
    dfdra = 5.0E-1*(1.4228426342100786056E-2*t[1]*t[2]*t[4]/pow(rhoa,2.3333333333333333333E0)-1.9694900436853928357E0*t[1]*(1.804E0-8.04E-1/t[3])*pow(rhoa,3.3333333333333333333E-1));
    dfdga = -5.3356598782877950386E-3*t[1]*t[4]*grada/pow(rhoa,1.3333333333333333333E0);
    ds->df1000 += factor*dfdra;
    ds->df0010 += factor*dfdga;
}

static void
pbex_second(SecondFuncDrv *ds, real factor,DftDensProp* dp)
{
    real t[9];
    real dfdra, dfdga;
    real d2fdraga, d2fdrara, d2fdgaga;
    real rhoa,rhob,grada,gradb;

    rhoa=dp->rhoa;
    rhob=dp->rhob;
    grada=dp->grada;
    gradb=dp->gradb;

    t[1] = pow(2.0E0,3.3333333333333333333E-1);
    t[2] = pow(grada,2.0E0);
    t[3] = 4.492799480231090353E-3*t[2]/pow(rhoa,2.6666666666666666667E0)+1.0E0;
    t[4] = 1/pow(t[3],2.0E0);
    t[5] = 1/pow(rhoa,2.3333333333333333333E0);
    t[6] = 1.804E0-8.04E-1/t[3];
    t[7] = 1/pow(rhoa,1.3333333333333333333E0);
    t[8] = 1/pow(t[3],3.0E0);
    dfdra = 5.0E-1*(1.4228426342100786056E-2*t[1]*t[2]*t[4]*t[5]-1.9694900436853928357E0*t[1]*t[6]*pow(rhoa,3.3333333333333333333E-1));
    dfdga = -5.3356598782877950386E-3*t[1]*t[4]*t[7]*grada;
    d2fdraga = 5.0E-1*(1.4228426342100786056E-2*t[1]*t[4]*t[5]*grada-2.5570186589718707735E-4*t[1]*t[8]*pow(grada,3.0E0)/pow(rhoa,5.0E0));
    d2fdrara = 5.0E-1*(3.4093582119624938689E-4*t[1]*t[8]*pow(grada,4.0E0)/pow(rhoa,6.0E0)-1.4228426342100792453E-2*t[1]*t[2]*t[4]/pow(rhoa,3.3333333333333333333E0)-6.5649668122846422854E-1*t[1]*t[6]/pow(rhoa,6.6666666666666666667E-1));
    d2fdgaga = 9.5888199711445165929E-5*t[1]*t[2]*t[8]/pow(rhoa,4.0E0)-5.3356598782877950386E-3*t[1]*t[4]*t[7];
    ds->df1000 += factor*dfdra;
    ds->df0010 += factor*dfdga;
    ds->df1010 += factor*d2fdraga;
    ds->df2000 += factor*d2fdrara;
    ds->df0020 += factor*d2fdgaga;
}

static void
pbex_third(ThirdFuncDrv *ds, real factor,DftDensProp* dp)
{
    real t[16];
    real dfdra, dfdga;
    real d2fdraga, d2fdrara, d2fdgaga;
    real d3fdraraga, d3fdragaga, d3fdrarara, d3fdgagaga;
    real rhoa,rhob,grada,gradb;

    rhoa=dp->rhoa;
    rhob=dp->rhob;
    grada=dp->grada;
    gradb=dp->gradb;

    t[1] = pow(2.0E0,3.3333333333333333333E-1);
    t[2] = pow(grada,2.0E0);
    t[3] = 4.492799480231090353E-3*t[2]/pow(rhoa,2.6666666666666666667E0)+1.0E0;
    t[4] = 1/pow(t[3],2.0E0);
    t[5] = 1/pow(rhoa,2.3333333333333333333E0);
    t[6] = 1.804E0-8.04E-1/t[3];
    t[7] = 1/pow(rhoa,1.3333333333333333333E0);
    t[8] = pow(grada,3.0E0);
    t[9] = 1/pow(t[3],3.0E0);
    t[10] = 1/pow(rhoa,5.0E0);
    t[11] = pow(grada,4.0E0);
    t[12] = 1/pow(rhoa,6.0E0);
    t[13] = 1/pow(rhoa,3.3333333333333333333E0);
    t[14] = 1/pow(rhoa,4.0E0);
    t[15] = 1/pow(t[3],4.0E0);
    dfdra = 5.0E-1*(1.4228426342100786056E-2*t[1]*t[2]*t[4]*t[5]-1.9694900436853928357E0*t[1]*t[6]*pow(rhoa,3.3333333333333333333E-1));
    dfdga= -5.3356598782877950386E-3*t[1]*t[4]*t[7]*grada;
    d2fdraga = 5.0E-1*(1.4228426342100786056E-2*t[1]*t[4]*t[5]*grada-2.5570186589718707735E-4*t[1]*t[10]*t[8]*t[9]);
    d2fdrara = 5.0E-1*(-6.5649668122846422854E-1*t[1]*t[6]/pow(rhoa,6.6666666666666666667E-1)+3.4093582119624938689E-4*t[1]*t[11]*t[12]*t[9]-1.4228426342100792453E-2*t[1]*t[13]*t[2]*t[4]);
    d2fdgaga = 9.5888199711445165929E-5*t[1]*t[14]*t[2]*t[9]-5.3356598782877950386E-3*t[1]*t[4]*t[7];
    d3fdraraga = 5.0E-1*(-9.1905376815760176701E-6*t[1]*t[15]*pow(grada,5.0E0)/pow(rhoa,8.6666666666666666667E0)-3.3199661464901841981E-2*t[1]*t[13]*t[4]*grada+1.6194451506821847512E-3*t[1]*t[12]*t[8]*t[9]);
    d3fdragaga = 5.0E-1*(6.8929032611820120636E-6*t[1]*t[11]*t[15]/pow(rhoa,7.6666666666666666667E0)-1.0228074635887482797E-3*t[1]*t[10]*t[2]*t[9]+1.4228426342100786056E-2*t[1]*t[4]*t[5]);
    d3fdrarara = 5.0E-1*(1.2254050242101353377E-5*t[1]*t[15]*pow(grada,6.0E0)/pow(rhoa,9.6666666666666666667E0)-2.3865507483737462337E-3*t[1]*t[11]*t[9]/pow(rhoa,7.0E0)+5.3751832847936318761E-2*t[1]*t[2]*t[4]/pow(rhoa,4.3333333333333333333E0)+4.376644541523094857E-1*t[1]*t[6]/pow(rhoa,1.6666666666666666667E0));
    d3fdgagaga = 2.8766459913433548026E-4*t[1]*t[14]*t[9]*grada-2.5848387229432553279E-6*t[1]*t[15]*t[8]/pow(rhoa,6.6666666666666666667E0);
    ds->df1000 += factor*dfdra;
    ds->df0010 += factor*dfdga;
    ds->df1010 += factor*d2fdraga;
    ds->df2000 += factor*d2fdrara;
    ds->df0020 += factor*d2fdgaga;
    ds->df2010 += factor*d3fdraraga;
    ds->df1020 += factor*d3fdragaga;
    ds->df3000 += factor*d3fdrarara;
    ds->df0030 += factor*d3fdgagaga;
}
