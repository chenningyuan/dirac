#include <cstdlib>

#ifndef NO_STDCXX
#include <iostream>
using namespace std;
#endif

#include "taylor.h"

#ifndef M_PI // M_PI is not standard for some reason
#define M_PI 3.14159265358979323846
#endif
#include "xc_fun.h"
#include "lda_fun.h"
#include "blyp_fun.h"
#include "pbe_fun.h"

static double xc_params[XC_PARAMS_LEN];

void xc_set_functional(const double params[XC_PARAMS_LEN])
{
  for (int i=0;i<XC_PARAMS_LEN;i++)
    xc_params[i] = params[i];
}

extern "C"
void xc_set_functional_(const double *params, int *len)
{
  assert(*len == XC_PARAMS_LEN && "C++/Fortran headers out of sync");
  xc_set_functional(params);
}

// This routine can be used independently of the number
// of variables, i.e. LDA/GGA spin polarized/unpolarized
template<class num>
static void sum_functionals(num &res,
			    const num &R,
			    const num &S,
			    const num &Z,
			    const num &X,
			    const num &Y,
			    const num &rhoa,
			    const num &rhob,
			    const num &grada2,
			    const num &gradab,
			    const num &gradb2,
			    const num &zeta)
{
  res = 0;

  if (xc_params[XC_VWN5_CORRELATION] != 0)
    res += xc_params[XC_VWN5_CORRELATION]*
      vwn5_correlation(R,S,zeta);

  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*
      slater_exchange(rhoa,rhob); 

  if (xc_params[XC_LYP_CORRELATION] != 0)
    res += xc_params[XC_LYP_CORRELATION]*
      lyp_correlation(R,rhoa,rhob,Z,grada2,gradb2);

  if (xc_params[XC_PBE_CORRELATION] != 0)
    res += xc_params[XC_PBE_CORRELATION]*
      pbe_correlation(R,S,Z,rhoa,rhob,zeta);

 if (xc_params[XC_OLD_PBE_CORRELATION] != 0)
    res += xc_params[XC_OLD_PBE_CORRELATION]*
      pbe_correlation_unpolarized(R,Z);

  if (xc_params[XC_BECKE_EXCHANGE] != 0)
    res += xc_params[XC_BECKE_EXCHANGE]*
      becke_exchange(rhoa,rhob,grada2,gradb2);

  if (xc_params[XC_PBE_EXCHANGE] != 0)
    res += xc_params[XC_PBE_EXCHANGE]*
      pbe_exchange(rhoa,rhob,grada2,gradb2);

  if (xc_params[XC_RPBE_EXCHANGE] != 0)
    res += xc_params[XC_RPBE_EXCHANGE]*
      rpbe_exchange(rhoa,rhob,grada2,gradab,gradb2);

  if (xc_params[XC_PW92_CORRELATION] != 0)
    res += xc_params[XC_PW92_CORRELATION]*
      pw92_correlation(R,zeta);

  res.deriv_facs();
}
			    

template<class num, int N>
static void gga_helper(num *result, const num d[5])
{
  taylor<num,5,N> &res = *reinterpret_cast<taylor<num,5,N> *>(result);
  const taylor<num,5,N>
#ifndef NEW_SCALING
    R(d[0],0,XC_FUNCTIONAL_SCALING),
    S(d[1],1,XC_FUNCTIONAL_SCALING),
    Z(d[2],2,XC_FUNCTIONAL_SCALING),
    X(d[3],3,XC_FUNCTIONAL_SCALING),
    Y(d[4],4,XC_FUNCTIONAL_SCALING),
#else
    R(d[0],0,d[0]),
    S(d[1],1,d[0]),
    Z(d[2],2,d[2]),
    X(d[3],3,d[2]),
    Y(d[4],4,d[2]),
#endif
    // Some functionals are easier described in terms of alpha and beta
    // quantities, prepare those here:
    rhoa = 0.5*(R+S),
    rhob = 0.5*(R-S),
    grada2 = 0.25*(Z+X+2*Y),
    gradab = 0.25*(Z-X),
    gradb2 = 0.25*(Z+X-2*Y),
    zeta = S/R;
  // Here we can in principle precalculate other common quantities
  // such as pow(rho,-1/3)

  // Add upp all active functionals
  sum_functionals(res,R,S,Z,X,Y,rhoa,rhob,grada2,gradab,gradb2,zeta);
}

template<class num, int N>
static void lda_helper(num *result, const num d[2])
{
  taylor<num,2,N> &res = *reinterpret_cast<taylor<num,2,N> *>(result);
  const taylor<num,2,N> 
    R(d[0],0,XC_FUNCTIONAL_SCALING),
    S(d[1],1,XC_FUNCTIONAL_SCALING),
    rhoa = 0.5*(R+S),
    rhob = 0.5*(R-S),
    zeta = S/R,
    nil(0);
  sum_functionals(res,R,S,nil,nil,nil,rhoa,rhob,nil,nil,nil,zeta);
}

void xc_eval_lda(double *result, int order, const double densvars[2])
{
  static void (*f[])(double *, const double *) = 
    { lda_helper<double,0>,lda_helper<double,1>,
#if XC_MAX_ORDER>=2
      lda_helper<double,2>,
#if XC_MAX_ORDER>=3
      lda_helper<double,3>,
#if XC_MAX_ORDER>=4
      lda_helper<double,4>,
#if XC_MAX_ORDER>=5
      lda_helper<double,5>,
#if XC_MAX_ORDER>=6
      lda_helper<double,6>,
#if XC_MAX_ORDER>=7
      lda_helper<double,7>,
#if XC_MAX_ORDER>=8
      lda_helper<double,8>
#if XC_MAX_ORDER>=9
#error Derivative order too large, update xc_fun.cpp
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
};
  assert(order>=0);
  assert(order<6);
  f[order](result,densvars);
}

extern "C"
void xc_eval_lda_(double *result, int *order, const double *densvars)
{
  xc_eval_lda(result,*order,densvars);
}

void xc_eval_gga(double *result, int order, const double densvars[5])
{
  static void (*f[])(double *, const double *) = 
    { gga_helper<double,0>,gga_helper<double,1>,
#if XC_MAX_ORDER>=2
      gga_helper<double,2>,
#if XC_MAX_ORDER>=3
      gga_helper<double,3>,
#if XC_MAX_ORDER>=4
      gga_helper<double,4>,
#if XC_MAX_ORDER>=5
      gga_helper<double,5>,
#if XC_MAX_ORDER>=6
      gga_helper<double,6>,
#if XC_MAX_ORDER>=7
      gga_helper<double,7>,
#if XC_MAX_ORDER>=8
      gga_helper<double,8>
#if XC_MAX_ORDER>=9
#error Derivative order too large, update xc_fun.cpp
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
};
  assert(order>=0);
#ifndef NO_STDCXX
  if (order>XC_MAX_ORDER)
    {
      cerr << "ERROR in xc_gga_eval(), order " << order << " greater than "
	   << "XC_MAX_ORDER (" << XC_MAX_ORDER << "), please reconfigure " 
	   << "libxc_fun." << endl;
      exit(EXIT_FAILURE);
    }
#else
  assert(order<=XC_MAX_ORDER);
#endif
  f[order](result,densvars);
}

extern "C"
void xc_eval_gga_(double *result, int *order, const double *densvars)
{
  xc_eval_gga(result,*order,densvars);
}

#ifdef WITH_QD

void xc_eval_gga_qd(qd_real *result, int order, const qd_real densvars[5])
{
  static void (*f[])(qd_real *, const qd_real *) = 
    { gga_helper<qd_real,0>,gga_helper<qd_real,1>,
#if XC_MAX_ORDER>=2
      gga_helper<qd_real,2>,
#if XC_MAX_ORDER>=3
      gga_helper<qd_real,3>,
#if XC_MAX_ORDER>=4
      gga_helper<qd_real,4>,
#if XC_MAX_ORDER>=5
      gga_helper<qd_real,5>,
#if XC_MAX_ORDER>=6
      gga_helper<qd_real,6>,
#if XC_MAX_ORDER>=7
      gga_helper<qd_real,7>,
#if XC_MAX_ORDER>=8
      gga_helper<qd_real,8>
#if XC_MAX_ORDER>=9
#error Derivative order too large, update xc_fun.cpp
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
};
  assert(order>=0);
#ifndef NO_STDCXX
  if (order>XC_MAX_ORDER)
    {
      cerr << "ERROR in xc_gga_eval(), order " << order << " greater than "
	   << "XC_MAX_ORDER (" << XC_MAX_ORDER << "), please reconfigure " 
	   << "libxc_fun." << endl;
      exit(EXIT_FAILURE);
    }
#else
  assert(order<=XC_MAX_ORDER);
#endif
  f[order](result,densvars);
}

#endif

// radovan: i copy-pasted the rest of this file from
//          autodiff/taylor.cpp
//          to have the taylor index available
//          in fortran unit tests
//          not sure this is how you planned it Ulf
//          and maybe it is already available under taylor
//          so feel free to undo this
static int taylor_len(int nvar, int ndeg)
  {
    int len = 1;
    for (int k=1;k<=nvar;k++)
      {
        len *= ndeg + k;
        len /= k;
      }
    return len;
  }

// Fortran interface to taylor_len
extern "C"
int taylor_len_(int *nvar, int *ndeg)
{
  return taylor_len(*nvar,*ndeg);
}

int taylor_index(int nvar, const int exponents[])
{
  int N = 0;
  for (int i=0;i<nvar;i++)
    N += exponents[i];
  int i = 0, idx = 0;
  N--;
  while (N >= 0)
    {
      idx += taylor_len(nvar-i,N);
      N -= exponents[i];
      i++;
    }
  return idx;
}

// Fortran version
extern "C"
int taylor_index_(const int *nvar, const int exponents[])
{
  return 1+taylor_index(*nvar,exponents);
}
