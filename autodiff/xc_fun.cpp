#include "xc_fun.h"
#include "auto_functionals.h"

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
static void lda_helper(double *result, const double densvars[2])
{
  taylor<2,N> &res = *reinterpret_cast<taylor<2,N> *>(result);
  res = 0;
  if (xc_params[XC_VWN_CORRELATION] != 0)
    res += xc_params[XC_VWN_CORRELATION]*vwn_correlation<N>(densvars);
  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*slater_exchange<N>(densvars);

  res.deriv_facs();
}

template<int N>
static void gga_helper(double *result, const double densvars[5])
{
  taylor<5,N> &res = *reinterpret_cast<taylor<5,N> *>(result);
  res = 0;
  if (xc_params[XC_LYP_CORRELATION] != 0)
    res += xc_params[XC_LYP_CORRELATION]*lyp_correlation<N>(densvars);
  if (xc_params[XC_PBE_CORRELATION] != 0)
    res += xc_params[XC_PBE_CORRELATION]*pbe_correlation<N>(densvars);
  if (xc_params[XC_BECKE_EXCHANGE] != 0)
    res += xc_params[XC_BECKE_EXCHANGE]*becke_exchange<N>(densvars);
  if (xc_params[XC_PBE_EXCHANGE] != 0)
    res += xc_params[XC_PBE_EXCHANGE]*pbe_exchange<N>(densvars);

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

