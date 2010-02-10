#include <iostream>
#ifdef WITH_WD
#include <qd/qd_real.h>
#endif
#include <cstdlib>
#include <cstdio>
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

template<class num, int Ndeg>
void print_ref(ostream &dst, num d[5], const polynomial<num,5,Ndeg> &p)
{
  cout << "rhoa= " << d[0] << " rhob= " << d[1] << " sigmaaa= " << d[2] << " sigmaab= " << d[4] << " sigmabb= " << d[3] << endl;
  cout << " zk            = " << p[0] << endl;
  if (Ndeg==0)
    return;
  cout << " vrhoa         = " << p[1] << endl;
  cout << " vrhob         = " << p[2] << endl;
  cout << " vsigmaaa      = " << p[3] << endl;
  cout << " vsigmaab      = " << p[5] << endl;
  cout << " vsigmabb      = " << p[4] << endl;
  if (Ndeg == 1)
    return;
  cout << endl;
  cout << " v2rhoa2       = " << p[6] << endl;
  cout << " v2rhoab       = " << p[7] << endl;
  cout << " v2rhob2       = " << p[11] << endl;
  cout << endl;
  cout << " v2rhoasigmaaa = " << p[8] << endl;
  cout << " v2rhoasigmaab = " << p[10] << endl;
  cout << " v2rhoasigmabb = " << p[9] << endl;
  cout << " v2rhobsigmaaa = " << p[12] << endl;
  cout << " v2rhobsigmaab = " << p[14] << endl;
  cout << " v2rhobsigmabb = " << p[13] << endl;
  cout << endl;
  cout << " v2sigmaaa2    = " << p[15] << endl;
  cout << " v2sigmaaaab   = " << p[17] << endl;
  cout << " v2sigmaaabb   = " << p[16] << endl;
  cout << " v2sigmaab2    = " << p[20] << endl;
  cout << " v2sigmaabbb   = " << p[19] << endl;
  cout << " v2sigmabb2    = " << p[18] << endl;


}

#define TEST_ORDER XC_MAX_ORDER

int main(int argc, char *argv[])
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

#ifdef WITH_SINGLE
  float fd[5] =   {1.0,0.0,1.0,0.0,0.0};
  taylor<float, 5, TEST_ORDER> res_f; 
#endif
  if (argc < 3)
    {
      cerr << "Usage: funeval FUNCTIONAL WEIGHT [FUNCTIONAL WEIGHT].." << endl;
      return -1;
    }
  for (int i=1;i+1<argc;i+=2)
    {
      int p = xc_parse_functional(argv[i]);
      if (p<0)
	{
	  cerr << "Unknown functional " << argv[i] << ", quitting." << endl;
	  return EXIT_FAILURE;
	}
      double w;
      if (sscanf(argv[i+1],"%lf",&w) != 1)
	{
	  cerr << argv[i+1] << " is not a float, quitting." << endl;
	  return EXIT_FAILURE;
	}
      cerr << "Setting " << argv[i] << " to " << w << endl;
      fun[p] = w;
    }
  xc_set_functional(XC_ABFGH,fun);
  cerr << "Input density: " << endl;
  cout.precision(20);
  cin >> dd[0] >> dd[1] >> dd[2] >> dd[3] >> dd[4];
#ifdef WITH_QD
  for (int i=0;i<2;i++)
    dqd[i] = dd[i];
#endif
#ifdef WITH_SINGLE
  for (int i=0;i<2;i++)
    fd[i] = dd[i];
#endif
  xc_eval(res_d.c,TEST_ORDER,dd);
  cout << "XC (double):" << endl;
  print_ref(cout,dd,res_d);
#ifdef WITH_SINGLE
  xc_eval_single(res_f.c,TEST_ORDER,fd);
  cout << "XC (single):" << endl;
  printpoly(cout,res_f);
#endif
#ifdef WITH_QD
  xc_eval_qd(res_qd.c,TEST_ORDER,dqd);
  cout << "XC (quad double):" << endl;
  printpoly(cout,res_qd);
  fpu_fix_end(&cw);
#endif
  return 0;
}
