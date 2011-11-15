#ifndef BOYS_H
#define BOYS_H

#include <cmath>
#include "taylor.h"
#include "boys_table.h"

namespace boys_internal
{

#define NTERMS 8 // Order of the interpolation in the table.

#define INVFAC_LEN 21
static const double invfac[INVFAC_LEN] = 
  {
    1,
    1,
    0.500000000000000,
    0.16666666666666667,
    0.041666666666666667,
    0.008333333333333333,
    0.001388888888888889,
    1.98412698412698e-04,
    2.48015873015873e-05,
    2.75573192239859e-06,
    2.75573192239859e-07,
    2.50521083854417e-08,
    2.08767569878681e-09,
    1.60590438368216e-10,
    1.14707455977297e-11,
    7.64716373181983e-13,
    4.77947733238739e-14,
    2.81145725434552e-15,
    1.56192069685862e-16,
    8.22063524662432e-18,
    4.11031762331216e-19,
  };

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

template<int Ndeg>
void boys_taylor_table(taylor<1,Ndeg>& val, double x)
{
  assert(NTERMS<=INVFAC_LEN);
  const int itab = int(x/BOYS_DELTA+0.5);
  assert(itab >= 0 && "Negative argument to boys()");
  assert(itab < BOYS_STEPS);
  double delta = x - itab*BOYS_DELTA;
  // Horner scheme
  for (int j=0;j<=Ndeg;j++)
    val[j] = boystab[itab][NTERMS-1+j]*invfac[NTERMS-1];
  for (int i=NTERMS-1;i>=0;i--)
    for (int j=0;j<=Ndeg;j++)
      val[j] = val[j]*delta + boystab[itab][i+j]*invfac[i];
  for (int j=0;j<=Ndeg;j++)
    val[j] *= invfac[j];
}

template<int Ndeg>
void boys_taylor(taylor<1,Ndeg>& res, tnum_t x0)
{
  if (x0 < BOYS_DELTA*(BOYS_STEPS-1)) 
    boys_taylor_table(res,x0);
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

template<int Nvar, int Ndeg, int Narg>
void boys_opt(taylor<Nvar,Ndeg>& res,
	      const taylor<Nvar,Narg>& arg)
{
  taylor<1,Ndeg> tmp;
  boys_internal::boys_taylor(tmp,arg[0]);
  arg.compose(res,tmp);
}

#endif
