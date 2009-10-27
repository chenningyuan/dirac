#include <iostream>
using namespace std;
#include "taylor.h"
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
void xc_set_functional_(const double *params)
{
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
    R(d[0],0,XC_FUNCTIONAL_SCALING),
    S(d[1],1,XC_FUNCTIONAL_SCALING),
    Z(d[2],2,XC_FUNCTIONAL_SCALING),
    X(d[3],3,XC_FUNCTIONAL_SCALING),
    Y(d[4],4,XC_FUNCTIONAL_SCALING),
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
  static void (*f[XC_MAX_ORDER+1])(double *, const double *) = 
    {lda_helper<double,0>,lda_helper<double,1>,
     lda_helper<double,2>,lda_helper<double,3>,
     lda_helper<double,4>,lda_helper<double,5>};
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

#ifdef WITH_QD
#include <qd/qd_real.h>
#endif
#include <iostream>
using namespace std;


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

#define TEST_ORDER 5

int main_gga(void)
{
#ifdef WITH_QD
  unsigned int cw;
  fpu_fix_start(&cw);
  taylor<qd_real, 5, TEST_ORDER> res_qd, cmp;
  qd_real dqd[5] = {1.0,0.0,1.0,0.0,0.0};
#endif
  double dd[5] =   {1.0,0.0,1.0,0.0,0.0};
  taylor<double, 5, TEST_ORDER> res_d;  
  double fun[XC_PARAMS_LEN] = {0};

  cout << "Input 5 density variables: ";
  cin >> dd[0] >> dd[1] >> dd[2] >> dd[3] >> dd[4];
#ifdef WITH_QD
  for (int i=0;i<5;i++)
    dqd[i] = dd[i];
#endif
  /*
  fun[XC_PBE_EXCHANGE] = 0;
  fun[XC_OLD_PBE_CORRELATION] = 1;
  xc_set_functional(fun);
  gga_helper<double,TEST_ORDER>(res_d.c,dd);
  cout << "Old PBE:" << endl;
  printpoly(cout,res_d);
  fun[XC_OLD_PBE_CORRELATION] = 0;
  */

  fun[XC_VWN5_CORRELATION] = 1;
  fun[XC_SLATER_EXCHANGE] = 1;
  xc_set_functional(fun);
  gga_helper<double,TEST_ORDER>(res_d.c,dd);
  //  cout << "New PBE:" << endl;
  //printpoly(cout,res_d);


#ifdef WITH_QD
  gga_helper<qd_real,TEST_ORDER>(res_qd.c,dqd);
  printpoly(cout,res_qd);
  res_d.convert_to(cmp);
  print_rms_diff(cmp,res_qd);
  fpu_fix_end(&cw);
#endif
  return 0;
}

int main(void)
{
#ifdef WITH_QD
  unsigned int cw;
  fpu_fix_start(&cw);
  taylor<qd_real, 2, TEST_ORDER> res_qd, cmp;
  qd_real dqd[5] = {1.0,0.0,1.0,0.0,0.0};
#endif
  double dd[5] =   {1.0,0.0,1.0,0.0,0.0};
  taylor<double, 2, TEST_ORDER> res_d;  
  double fun[XC_PARAMS_LEN] = {0};

  cout << "Input 5 density variables: ";
  cin >> dd[0] >> dd[1] >> dd[2] >> dd[3] >> dd[4];
#ifdef WITH_QD
  for (int i=0;i<5;i++)
    dqd[i] = dd[i];
#endif
  /*
  fun[XC_PBE_EXCHANGE] = 0;
  fun[XC_OLD_PBE_CORRELATION] = 1;
  xc_set_functional(fun);
  gga_helper<double,TEST_ORDER>(res_d.c,dd);
  cout << "Old PBE:" << endl;
  printpoly(cout,res_d);
  fun[XC_OLD_PBE_CORRELATION] = 0;
  */

  fun[XC_VWN5_CORRELATION] = 1;
  fun[XC_SLATER_EXCHANGE] = 1;
  xc_set_functional(fun);
  lda_helper<double,TEST_ORDER>(res_d.c,dd);
  //  cout << "New PBE:" << endl;
  cout << "Double:" << endl;
  printpoly(cout,res_d);


#ifdef WITH_QD
  lda_helper<qd_real,TEST_ORDER>(res_qd.c,dqd);
  printpoly(cout,res_qd);
  res_d.convert_to(cmp);
  print_rms_diff(cmp,res_qd);
  fpu_fix_end(&cw);
#endif
  return 0;
}

#endif
