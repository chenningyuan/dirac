#include <iostream>
#include <qd/qd_real.h>
#include "taylor.h"
#include "xc_fun.h"

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

#define TEST_ORDER 2

int main(void)
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

  fun[XC_SLATER_EXCHANGE] = 1;
  fun[XC_VWN5_CORRELATION] = 1;
  //fun[XC_BECKE_EXCHANGE] = 1;
  xc_set_functional(XC_RSZXY,fun);
  cout.precision(20);
  cin >> dd[0] >> dd[1] >> dd[2] >> dd[3] >> dd[4];
#ifdef WITH_QD
  for (int i=0;i<5;i++)
    dqd[i] = dd[i];
#endif
  xc_eval(res_d.c,TEST_ORDER,dd);
  cout << "XC (double):" << endl;
  printpoly(cout,res_d);
#ifdef WITH_QD
  xc_eval_qd(res_qd.c,TEST_ORDER,dqd);
  cout << "XC (quad double):" << endl;
  printpoly(cout,res_qd);
  fpu_fix_end(&cw);
#endif
  return 0;
}
