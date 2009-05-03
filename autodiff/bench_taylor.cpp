#include <iostream>
#include <time.h>
#include "taylor.h"

using namespace std;

// For Boys see also Arndt, p 671
template<int Ndeg>
void boys_taylor(taylor<1,Ndeg>& res, double x0)
{
  taylor<1,Ndeg> sqrtt;
  pow_taylor(sqrtt,x0,0.5);
  res = (0.5*sqrt(M_PI))/sqrtt;
  res *= erf(sqrtt);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> boys(const taylor<Nvar,Ndeg>& t)
{
  taylor<1,Ndeg> tmp;
  boys_taylor(tmp,t[0]);
  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

double boys(const double& t)
{
  double sqrtt = sqrt(t);
  return 0.5*sqrt(M_PI)/sqrtt*erf(sqrtt);
}

#define INNER 10000
#ifndef NVAR
#define NVAR 3
#endif
#ifndef NDEG
#define NDEG 5
#endif
#define NSEC 10

int main(void)
{
  clock_t ct;
  long long ncalls = 0;
  taylor<NVAR,NDEG> res(0),inp(2.5,0);
  //double resd = 0;
  while ((ct = clock())/CLOCKS_PER_SEC < NSEC)
    {
      for (int i=0;i<INNER;i++)
	{
	  //resd += boys(2.5+int(ncalls/100 % 6) );
	  res += boys(inp);
	  ncalls++;
	}
    }
  cout << "taylor<" << NVAR << ", " << NDEG << "> coeffs/cpu s: " << binomial<NVAR+NDEG,NDEG>::value*(ncalls*double(CLOCKS_PER_SEC))/ct << endl;
  //return int(resd);
  return int(res[0]);
}
