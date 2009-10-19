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

template<class num, int N>
static void gga_helper(num *result, const num d[5])
{
  taylor<num,5,N> &res = *reinterpret_cast<taylor<num,5,N> *>(result);
  const taylor<num,5,N> 
    R(d[0],0,XC_FUNCTIONAL_SCALING),
    S(d[1],1,XC_FUNCTIONAL_SCALING),
    Z(d[2],2,XC_FUNCTIONAL_SCALING),
    X(d[3],3,XC_FUNCTIONAL_SCALING),
    Y(d[4],4,XC_FUNCTIONAL_SCALING);
  // Some functionals are easier described in terms of alpha and beta
  // quantities, prepare those here:
  const taylor<num,5,N> 
    rhoa = 0.5*(R+S),
    rhob = 0.5*(R-S),
    grada2 = 0.25*(Z+X+2*Y),
    gradab = 0.25*(Z-X),
    gradb2 = 0.25*(Z+X-2*Y);
  // Here we can in principle precalculate other common quantities
  // such as pow(rho,-1/3)

  // Add upp all active functionals
  res = 0;
  if (xc_params[XC_VWN_CORRELATION] != 0)
    res += xc_params[XC_VWN_CORRELATION]*
      vwn_correlation(R,S);

  if (xc_params[XC_SLATER_EXCHANGE] != 0)
    res += xc_params[XC_SLATER_EXCHANGE]*
      slater_exchange(R,S); 

  if (xc_params[XC_LYP_CORRELATION] != 0)
    res += xc_params[XC_LYP_CORRELATION]*
      lyp_correlation(R,rhoa,rhob,Z,grada2,gradb2);

  if (xc_params[XC_PBE_CORRELATION] != 0)
    res += xc_params[XC_PBE_CORRELATION]*
      pbe_correlation_new(R,S,Z,rhoa,rhob);

  if (xc_params[XC_BECKE_EXCHANGE] != 0)
    res += xc_params[XC_BECKE_EXCHANGE]*
      becke_exchange(rhoa,rhob,grada2,gradb2);

  if (xc_params[XC_PBE_EXCHANGE] != 0)
    res += xc_params[XC_PBE_EXCHANGE]*
      pbe_exchange(rhoa,rhob,grada2,gradb2);

  if (xc_params[XC_RPBE_EXCHANGE] != 0)
    res += xc_params[XC_RPBE_EXCHANGE]*
      rpbe_exchange(rhoa,rhob,grada2,gradab,gradb2);

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
    {gga_helper<double,0>,gga_helper<double,1>,gga_helper<double,2>,
     gga_helper<double,3>,gga_helper<double,4>,gga_helper<double,5>};
  assert(order>=0);
  assert(order<6);
  f[order](result,densvars);
}

extern "C"
void xc_eval_gga_(double *result, int *order, const double *densvars)
{
  xc_eval_gga(result,*order,densvars);
}

#ifdef XC_FUN_NUMTEST

#include <qd/qd_real.h>
#include <iostream>
using namespace std;

template<class num, int Nvar, int Ndeg>
void printpoly(ostream &dst, const polynomial<num,Nvar,Ndeg> &p)
{
  int exps[Nvar] = {0};
  for (int i=0;i<p.size;i++)
    {
      p.exponents(i,exps);
      for (int j=0;j<Nvar;j++)
	cout << exps[j] << " ";
      cout << " " << p[i] << endl;
    }
}

template<class num, int Nvar, int Ndeg>
void print_rms_diff(const polynomial<num,Nvar,Ndeg> &p1, 
		    const polynomial<num,Nvar,Ndeg> &p2)
{
  cout << "Order \tRMS \t\tRMS_Error\tRelative_Error" << endl;
  for (int n=0;n<=Ndeg;n++)
    {
      num rms = 0.0, avg = 0.0;
      for (int i=polylen(Nvar,n-1);i<polylen(Nvar,n);i++)
	{
	  rms += pow(p1[i] - p2[i],2);
	  avg += p1[i]*p1[i] + p2[i]*p2[i];
	}
      avg = sqrt(avg/(2*polylen(Nvar-1,n)));
      rms = sqrt(rms/polylen(Nvar-1,n));
      cout << n << " \t" << avg << " \t" << rms <<  " \t" << rms/avg << endl;
    }
}

#define TEST_ORDER 6

int main(void)
{
  unsigned int cw;
  fpu_fix_start(&cw);
  taylor<double, 5, TEST_ORDER> res_d;  
  taylor<qd_real, 5, TEST_ORDER> res_qd, cmp;
  double dd[5] =   {0.2,0.1,0.1,0.0,0.1};
  qd_real dqd[5] = {0.2,0.1,0.1,0.0,0.1};
  double fun[XC_PARAMS_LEN] = {0};
  /*  fun[XC_SLATER_EXCHANGE] = 1;
      fun[XC_VWN_CORRELATION] = 1; */
  fun[XC_PBE_EXCHANGE] = 1;
  xc_set_functional(fun);
  gga_helper<double,TEST_ORDER>(res_d.c,dd);
  gga_helper<qd_real,TEST_ORDER>(res_qd.c,dqd);
  res_d.convert_to(cmp);
  print_rms_diff(cmp,res_qd);
  fpu_fix_end(&cw);
  return 0;
}

#endif
