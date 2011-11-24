#include <cmath>
#include <cassert>
#include <iostream>
#include <qd/qd_real.h>
#include "boys_table.h"

using namespace std;

#define NTERMS 8

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

double boys(double x, int n)
{
  assert(NTERMS<=INVFAC_LEN);
  assert(n>=0);
  assert(n<BOYS_MAXDEG-NTERMS+1);
  const int itab = int(x/BOYS_DELTA+0.5);
  assert(itab >= 0 && "Negative argument to boys()");
  assert(itab < BOYS_STEPS);
  double delta = x - itab*BOYS_DELTA;
  // Horner scheme
  double val = boystab[itab][n+NTERMS-1]*invfac[NTERMS-1];
  for (int i=NTERMS-1;i>=0;i--)
    val = val*delta + boystab[itab][n+i]*invfac[i];
  return val;
}

//Comput taylor coefficients for the boys function at x.
template<int Ndeg>
void boys_taylor(double x, double c[])
{
  assert(NTERMS<=INVFAC_LEN);
  assert(n>=0);
  assert(n<BOYS_MAXDEG-NTERMS+1);
  const int itab = int(x/BOYS_DELTA+0.5);
  assert(itab >= 0 && "Negative argument to boys()");
  assert(itab < BOYS_STEPS);
  double delta = x - itab*BOYS_DELTA;
  // Horner scheme
  for (int j=0;j<=Ndeg;j++)
    val[j] = boystab[itab][n+NTERMS-1+j]*invfac[NTERMS-1];
  for (int i=NTERMS-1;i>=0;i--)
    for (int j=0;j<=Ndeg;j++)
      val[j] = val[j]*delta + boystab[itab][n+i+j]*invfac[i];
  for (int j=0;j<=Ndeg;j++)
    val[j] *= invfac[j];
}

static qd_real boys_from_taylor(qd_real x, int n)
{
  qd_real xn = 1, sum = 0;
  int k = 0;
  while (1)
    {
      qd_real term = xn/(2*n+1+2*k);
      sum += term;
      if (fabs(term/sum) < 1e-53)
	return pow(-1,n)*sum;
      k++;
      xn = xn*(-x)/k;
    }
}

int main(void)
{
  while (!cin.fail())
    {
      double x;
      cin >> x;
      cout << (boys(x,8) - boys_from_taylor(x,8))/boys_from_taylor(x,8) << endl;
    }
  return 0;
}
