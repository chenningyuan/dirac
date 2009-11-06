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

template<class T> 
struct densvars
{
  // Variables for expressing functionals, these are redundant because
  // different functionals have different needs.
  T R, S, Z, X, Y;
  T rhoa, rhob, grada2, gradab, gradb2;
  T zeta; //S/R
};

//Did we already set up tables?
static int is_setup = 0;
//Which primary variables will we use?
static enum xc_mode current_mode = XC_ORDER_RS;
//Current density functional type and parameters
static enum xc_fun_type current_type = XC_GGA;
static double xc_params[XC_PARAMS_LEN];
//Dispatch table for different orders, set by setup()
static void (*helper_tab[XC_MAX_ORDER+1])(double *, const double *);

static void setup(void);

void xc_set_functional(enum xc_fun_type funtype, 
		       const double params[XC_PARAMS_LEN])
{
  if (!is_setup)
    setup();
  assert(funtype == XC_LGA or funtype == XC_GGA);
  current_type = funtype;
  for (int i=0;i<XC_PARAMS_LEN;i++)
    xc_params[i] = params[i];
  //TODO: for LDA, check that the user is not asking for
  //GGA functionals.
}

extern "C"
void xc_set_functional_(const int *funtype, const double *params, int *len)
{
  assert(*len == XC_PARAMS_LEN && "C++/Fortran headers out of sync");
  xc_set_functional((enum xc_fun_type)(*funtype),params);
}

static int fun_nvars(enum xc_fun_type funtype)
{
  switch(funtype)
    {
    case XC_LDA:
      return 2;
    case XC_GGA:
      return 5;
    default:
      assert(0 && "Invalid funtype");
      return -1;
    }
}

// This routine can be used independently of the number
// of variables, i.e. LDA/GGA spin polarized/unpolarized
template<class num>
static void sum_functionals(num &res,
			    const densvars<num> &d)
{
  res = 0;

  if (xc_params[XC_VWN5_CORRELATION] != 0)
    res += xc_params[XC_VWN5_CORRELATION]*
      vwn5_correlation(d.R,d.S,d.zeta);

  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*
      slater_exchange(d.rhoa,d.rhob); 

  if (xc_params[XC_LYP_CORRELATION] != 0)
    res += xc_params[XC_LYP_CORRELATION]*
      lyp_correlation(d.R,d.rhoa,d.rhob,d.Z,d.grada2,d.gradb2);

  if (xc_params[XC_PBE_CORRELATION] != 0)
    res += xc_params[XC_PBE_CORRELATION]*
      pbe_correlation(d.R,d.S,d.Z,d.rhoa,d.rhob,d.zeta);

 if (xc_params[XC_OLD_PBE_CORRELATION] != 0)
    res += xc_params[XC_OLD_PBE_CORRELATION]*
      pbe_correlation_unpolarized(d.R,d.Z);

  if (xc_params[XC_BECKE_EXCHANGE] != 0)
    res += xc_params[XC_BECKE_EXCHANGE]*
      becke_exchange(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_PBE_EXCHANGE] != 0)
    res += xc_params[XC_PBE_EXCHANGE]*
      pbe_exchange(d.rhoa,d.rhob,d.grada2,d.gradb2);

  if (xc_params[XC_RPBE_EXCHANGE] != 0)
    res += xc_params[XC_RPBE_EXCHANGE]*
      rpbe_exchange(d.rhoa,d.rhob,d.grada2,d.gradab,d.gradb2);

  if (xc_params[XC_PW92_CORRELATION] != 0)
    res += xc_params[XC_PW92_CORRELATION]*
      pw92_correlation(d.R,d.zeta);

  res.deriv_facs();
}
			    
template<class num, int N>
static void gga_setdens_rs(densvars<taylor<num,5,N> >&d, const num dd[5])
{
  d.R = taylor<num,5,N>(dd[0],0);
  d.S = taylor<num,5,N>(dd[1],1);
  d.Z = taylor<num,5,N>(dd[2],2);
  d.X = taylor<num,5,N>(dd[3],3);
  d.Y = taylor<num,5,N>(dd[4],4);
  d.rhoa = 0.5*(d.R+d.S);
  d.rhob = 0.5*(d.R-d.S);
  d.grada2 = 0.25*(d.Z+d.X+2*d.Y);
  d.gradab = 0.25*(d.Z-d.X);
  d.gradb2 = 0.25*(d.Z+d.X-2*d.Y);
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void gga_setdens_ab(densvars<taylor<num,5,N> >&d, const num dd[5])
{
  d.rhoa = taylor<num,5,N>(dd[0],0);
  d.rhob = taylor<num,5,N>(dd[1],1);
  d.grada2 = taylor<num,5,N>(dd[2],2);
  d.gradb2 = taylor<num,5,N>(dd[3],3);
  d.gradab = taylor<num,5,N>(dd[4],4);
  d.R = d.rhoa + d.rhob;
  d.S = d.rhoa - d.rhob;
  d.Z = d.grada2+2*d.gradab+d.gradb2; 
  d.X = d.grada2-2*d.gradab+d.gradb2; 
  d.Y = d.grada2-d.gradb2;
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void lda_setdens_rs(densvars<taylor<num,2,N> >&d, const num dd[2])
{
  d.R = taylor<num,2,N>(dd[0],0);
  d.S = taylor<num,2,N>(dd[1],1);
  d.rhoa = 0.5*(d.R+d.S);
  d.rhob = 0.5*(d.R-d.S);
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void lda_setdens_ab(densvars<taylor<num,2,N> >&d, const num dd[5])
{
  d.rhoa = taylor<num,2,N>(dd[0],0);
  d.rhob = taylor<num,2,N>(dd[1],1);
  d.R = d.rhoa + d.rhob;
  d.S = d.rhoa - d.rhob;
  d.zeta = d.S/d.R;
}

template<class num, int N>
static void fun_helper(num *result, const num d[5])
{
  if (current_type == XC_LDA)
    {
      taylor<num,2,N> &res = *reinterpret_cast<taylor<num,2,N> *>(result);
      densvars<taylor<num,2,N> > dv;
      if (current_mode == XC_ORDER_RS)
	lda_setdens_rs(dv,d);
      else
	lda_setdens_ab(dv,d);
      sum_functionals(res,dv);
    }
  else if (current_type == XC_GGA)
    {
      taylor<num,5,N> &res = *reinterpret_cast<taylor<num,5,N> *>(result);
      densvars<taylor<num,5,N> > dv;
      if (current_mode == XC_ORDER_RS)
	gga_setdens_rs(dv,d);
      else
	gga_setdens_ab(dv,d);
      sum_functionals(res,dv);
    }
  else
    {
      assert(0 && "BUG: Invalid functional type");
    }
}

void xc_eval(double *result, int order, const double densvars[2])
{
  assert(order>=0);
  assert(order<=XC_MAX_ORDER);
  helper_tab[order](result,densvars);
}

extern "C"
void xc_eval_(double *result, int *order, const double *densvars)
{
  xc_eval(result,*order,densvars);
}

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

int xc_len(enum xc_fun_type funtype, int order)
{
  return taylor_len(fun_nvars(funtype),order);
}

int xc_len_(int *funtype, int *order)
{
  return xc_len((enum xc_fun_type)*funtype,*order);
}

int xc_index(enum xc_fun_type funtype, const int exponents[])
{
  int nvar = fun_nvars(funtype);
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
int xc_index_(const int *funtype, const int exponents[])
{
  return 1+xc_index((enum xc_fun_type)*funtype,exponents);
}


template<int N>
void setup_order(void)
{
  helper_tab[N] = fun_helper<double,N>;
  setup_order<N-1>();
}

template<>
void setup_order<-1>(void) {}


void setup()
{
  setup_order<XC_MAX_ORDER>();
  is_setup = 1;
}

