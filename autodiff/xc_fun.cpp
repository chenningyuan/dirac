#include "taylor.h"
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
void xc_set_functional_(const double *params)
{
  xc_set_functional(params);
}

template<int N>
static void lda_helper(double *result, const double d[2])
{
  taylor<double,2,N> &res = *reinterpret_cast<taylor<double,2,N> *>(result);
  const taylor<double,2,N> 
    R(d[0],0,XC_FUNCTIONAL_SCALING), 
    S(d[1],1,XC_FUNCTIONAL_SCALING);
  res = 0;
  if (xc_params[XC_VWN_CORRELATION] != 0)
    res += xc_params[XC_VWN_CORRELATION]*vwn_correlation(R,S);
  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*slater_exchange(R,S);

  res.deriv_facs();
}

// Insert zeros to get five variable layout for
// lda results.
template<int N>
static taylor<double,5,N> as_gga(const taylor<double,2,N> &lda)
{
  taylor<double,5,N> gga = 0;
  int exponents[5] = {0,0,0,0,0};
  for (int i = 0; i < lda.size; i++)
    {
      lda.exponents(i,exponents);
      gga[gga.term_index(exponents)] = lda[i];
    }
  return gga;
}

template<int N>
static void gga_helper(double *result, const double d[5])
{
  taylor<double,5,N> &res = *reinterpret_cast<taylor<double,5,N> *>(result);
  const taylor<double,5,N> 
    R(d[0],0,XC_FUNCTIONAL_SCALING),
    S(d[1],1,XC_FUNCTIONAL_SCALING),
    Z(d[2],2,XC_FUNCTIONAL_SCALING),
    X(d[3],3,XC_FUNCTIONAL_SCALING),
    Y(d[4],4,XC_FUNCTIONAL_SCALING);
  // Some functionals are easier described in terms of alpha and beta
  // quantities, prepare those here:
  const taylor<double,5,N> 
    rhoa = 0.5*(R+S),
    rhob = 0.5*(R-S),
    grada2 = 0.25*(Z+X+2*Y),
    gradab = 0.25*(Z-X),
    gradb2 = 0.25*(Z+X-2*Y);
  // For lda:
  const taylor<double,2,N> 
    Rlda(d[0],0,XC_FUNCTIONAL_SCALING),
    Slda(d[1],1,XC_FUNCTIONAL_SCALING);
  // Add upp all active functionals
  res = 0;
  if (xc_params[XC_VWN_CORRELATION] != 0)
    res += xc_params[XC_VWN_CORRELATION]*
      as_gga(vwn_correlation(Rlda,Slda));

  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*
      as_gga(slater_exchange(Rlda,Slda)); 

  if (xc_params[XC_LYP_CORRELATION] != 0)
    res += xc_params[XC_LYP_CORRELATION]*
      lyp_correlation(R,rhoa,rhob,Z,grada2,gradb2);

  if (xc_params[XC_PBE_CORRELATION] != 0)
    res += xc_params[XC_PBE_CORRELATION]*
      pbe_correlation(R,Z);

  if (xc_params[XC_BECKE_EXCHANGE] != 0)
    res += xc_params[XC_BECKE_EXCHANGE]*
      becke_exchange(rhoa,rhob,grada2,gradb2);

  if (xc_params[XC_PBE_EXCHANGE] != 0)
    res += xc_params[XC_PBE_EXCHANGE]*
      pbe_exchange(rhoa,rhob,grada2,gradb2);

  res.deriv_facs();
}

void xc_eval_lda(double *result, int order, const double densvars[2])
{
  static void (*f[XC_MAX_ORDER+1])(double *, const double *) = 
    {lda_helper<0>,lda_helper<1>,lda_helper<2>,lda_helper<3>,
     lda_helper<4>,lda_helper<5>};
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
  static void (*f[XC_MAX_ORDER+1])(double *, const double *) = 
    {gga_helper<0>,gga_helper<1>,gga_helper<2>,gga_helper<3>,
     gga_helper<4>,gga_helper<5>};
  assert(order>=0);
  assert(order<6);
  f[order](result,densvars);
}

extern "C"
void xc_eval_gga_(double *result, int *order, const double *densvars)
{
  xc_eval_gga(result,*order,densvars);
}

