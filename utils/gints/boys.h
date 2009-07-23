#ifndef BOYS_H
#define BOYS_H

#include <cmath>
#include "taylor.h"

namespace boys_internal
{
static double fac(int n)
{
  double f = 1;
  while (n > 1)
    f *= n--;
  return f;
}

static double boys_taylor(double x, int n)
{
  long double xn=1, sum = 0;
  int k = 0;
  while (1)
    {
      long double term = xn/(2*n+1+2*k);
      sum += term;
      if (fabs(term/sum) < 1e-16)
	return sum;
      k++;
      xn = xn*(-x)/k;
    }
}


// This method is accurate for the first 16 terms
// (relative accuracy better than 10^-14)
// for x0 >= 10
// It is _very_ in-accurate for x0 < 1 due to cancellation.
template<int Ndeg>
static void boys_taylor_largex(taylor<1,Ndeg>& res, tnum_t x0)
{
  taylor<1,Ndeg> sqrtt;
  pow_taylor(sqrtt,x0,0.5);
  res = (0.5*sqrt(TAYLOR_PI))/sqrtt;
  res *= erf(sqrtt);
}

// Boys function for small arguments, based on downward
// recursion.
template<int Ndeg>
static void boys_taylor_smallx(taylor<1,Ndeg>& res, double x)
{
  double expmx = exp(-x);
  int n = Ndeg;
  res[Ndeg] = boys_taylor(x,Ndeg);
  while (n>=0)
    {
      n--;
      if (n >= 0)
	res[n] = (2*x*res[n+1] + expmx)/(2*n+1);
    }
  for (n=0;n<=Ndeg;n++)
    res[n] /= pow(-1,n)*fac(n);
}

// At worst x (=5) Currently accurate to (relative) 
// 1e-12 for 15:th derivative, others better.
template<int Ndeg>
void boys_taylor(taylor<1,Ndeg>& res, tnum_t x0)
{
  if (x0 < 5) // 5 seems to give balanced accuracy, but smallx should be improved.
    boys_taylor_smallx(res,x0);
  else
    boys_taylor_largex(res,x0);
}

} // end of namespace boys_internal

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> boys(const taylor<Nvar,Ndeg>& arg)
{
  taylor<1,Ndeg> tmp;
  boys_internal::boys_taylor(tmp,arg[0]);
  taylor<Nvar,Ndeg> res;
  arg.compose(res,tmp);
  return res;
}

#endif
