// By Ulf Ekstrom 2009
#ifndef TAYLOR_H
#define TAYLOR_H

/*TODO:

  * Implement shifting the Taylor series by some amount.

  * Right now the Boys benchmark is 50x slower for <3,6> than <3,5>,
  because <3,6> uses the generic code. Make the generic code a
  template, and make the recursion terminate at the tables instead of
  continuing all the way down.
*/


#include <cmath>

typedef double tnum_t; // Type of coefficients
typedef long double tlongnum_t; // Type used to accumulate sums

#include "taylor_multab.h"

// Template to evaluate binomial coefficients at compile time.
template <int n, int k>
struct binomial
{
    enum { value = binomial<n-1,k-1>::value + binomial<n-1,k>::value };
};

template <int n>
struct binomial<n,n>
{
    enum { value = 1 };
};

template <int n>
struct binomial<n,0>
{
    enum { value = 1 };
};

template <int k>
struct binomial<0,k>
{
    enum { value = 1 };
};

// Number of terms in a polynomial of nvar variables
// up to and including order ndeg.
static inline int polylen(int nvar, int ndeg)
{
  if (ndeg < 0)
    return 0;
  int len = 1;
  for (int k=1;k<=nvar;k++)
    {
      len *= ndeg + k;
      len /= k;
    }
  return len;
}

// Return the index of the term with exponents
// This is a slow operation
int poly_index(int nvar, const int exponents[]);

// Multiply polynomials with coefficients p1 and p2,
// and _add_ the result to dst.
void taylor_mul(double dst[], int nvar, int ndeg,
		const double p1[], const double p2[]);

// A polynomial in Nvar variables of order Ndeg (all terms
// x^iy^j..z^k: i+j+..+k <= Ndeg). Multiplication
// truncates all higher order terms.
template<int Nvar,int Ndeg>
class taylor
{
 public:
  taylor(void) {}
  taylor(tnum_t c0) // Set the constant term to c0 and all other to 0.
    {
      c[0] = c0;
      for (int i=1;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	c[i] = 0;
    }
  // Set the constant term to c0 and the first order term of
  // variable var to 1.
  taylor(tnum_t c0, int var)
    {
      c[0] = c0;
      for (int i=1;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	c[i] = 0;
      assert(Ndeg>=0);
      assert(Nvar>=1);
      assert(var>=0);
      assert(var<Nvar);
      c[1+var] = 1;
    }
  // Evaluate polynomial or derivative of polynomial at x.
  // This is not a fast or accurate function.
  tnum_t eval(tnum_t x, int derivative_order = 0)
    {
      assert(Nvar == 1);
      tlongnum_t xn = 1, res = 0, a = 1;
      for (int i=2;i<=derivative_order;i++)
	a *= i;
      for (int i=derivative_order;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	{
	  res += a*xn*c[i];
	  xn *= x;
	  a *= (i + 1)/(tlongnum_t)(i+1-derivative_order);
	}
      return res;
    }
  // Evaluate the taylor series using the sumalt algorithm.
  // Note: This is for the value only, not derivatives!
  tnum_t eval_sumalt(tnum_t x, tnum_t e = 1.0)
    {
      assert(Nvar == 1); // TODO: general case
      tlongnum_t d = pow(2*e+1+2*sqrt(e*(e+1)),binomial<Ndeg+Nvar,Ndeg>::value);
      d = (d + 1/d)/2;
      tlongnum_t b = 1;
      tlongnum_t cc = d;
      tlongnum_t s = 0;
      tlongnum_t xn = 1;
      for (int k=0;k<binomial<Ndeg+Nvar,Ndeg>::value;k++)
	{
	  cc -= b;
	  s += cc*xn*c[k];
	  xn *= x;
	  b *= (2*(binomial<Ndeg+Nvar,Ndeg>::value+k)*(binomial<Ndeg+Nvar,Ndeg>::value-k))/tnum_t((2*k+1)*(k+1))*e;
	}
      return s/d;
    }

  /* Shift the polynomial so that x:= x + d[0], y:= y + d[1] etc*/
  void shift(const tnum_t d[])
  {
  }
  void integrate(void)
    {
      assert(Nvar == 1 && "Implement this..");
      for (int i=Ndeg;i>0;i--)
	c[i] = c[i-1]/i;
      c[0] = 0;
    }
  // Set x := alpha*x, y = :alpha*y etc
  void stretch(tnum_t alpha)
    {
      tnum_t an = alpha;
      int k = 1;
      for (int i=1;i<=Ndeg;i++)
	{
	  for (;k<polylen(Nvar,i);k++)
	    c[k] *= an;
	  an *= alpha;
	}
    }
  void zero(void)
  {
    for (int i=0;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	c[i] = 0;
  }
  int size(void) const
    {
      return binomial<Ndeg+Nvar,Ndeg>::value;
    }
  // Lexicographically ordered access to terms
  tnum_t operator[](int i) const
    {
      assert(i>=0);
      assert(i<size());
      return c[i];
    }
  tnum_t& operator[](int i)
    {
      assert(i>=0);
      assert(i<size());
      return c[i];
    }
  void operator-=(const taylor<Nvar,Ndeg>& t)
    {
      for (int i=0;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	c[i] -= t.c[i];
    }
  void operator+=(const taylor<Nvar,Ndeg>& t)
    {
      for (int i=0;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	c[i] += t.c[i];
    }
  void operator*=(const tnum_t& scale)
    {
      for (int i=0;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	c[i] *= scale;
    }
  void operator*=(const taylor<Nvar,Ndeg>& t)
    {
      if(Ndeg <= TAYLOR_MULTAB_NDEG &&
	 Nvar <= TAYLOR_MULTAB_NVAR)
	{
	  taylor_multiplier<Nvar,Ndeg>::acc(c,t.c);
	}
      else
	{
	  tnum_t tmp[binomial<Ndeg+Nvar,Ndeg>::value];
	  taylor_mul(tmp,Nvar,Ndeg,t.c,c);
	  for (int i=0;i<binomial<Ndeg+Nvar,Ndeg>::value;i++)
	    c[i] = tmp[i];
	}
    }
  taylor<Nvar,Ndeg> operator-(void)
    {
      taylor<Nvar,Ndeg> res = *this;
      for (int i=0;i<res.size();i++)
	res[i] *= -1;
      return res;
    }
  /* Put sum_i coeff[i]*(this - this[0])^i in res,
     used when evaluating analytical functions of this
  */
  void compose(taylor<Nvar,Ndeg>& res, const taylor<1,Ndeg> &coeff) const
    {
      if(Ndeg <= TAYLOR_MULTAB_NDEG &&
	 Nvar <= TAYLOR_MULTAB_NVAR)
	{
	  taylor_multiplier<Nvar,Ndeg>::compose(res.c,coeff.c,c);
	}
      else
	{
	  // Slow generic code
	  taylor<Nvar,Ndeg> tmp = *this;
	  tmp[0] = 0;
	  res.zero();
	  for (int i=Ndeg;i>0;i--)
	    {
	      res[0] += coeff[i];
	      res *= tmp;
	    }
	  res[0] += coeff[0];
	}
    }
  // protected: //Need taylor<1,Ndeg> be friends with taylor<Nvar,Ndeg>, but how??
  tnum_t c[binomial<Ndeg+Nvar,Ndeg>::value];
};

// Operators and stuff

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator*(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp *= x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator*(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp *= x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator*(const taylor<Nvar,Ndeg>& t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<Nvar,Ndeg> tmp = t1;
  tmp *= t2;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator+(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator+(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp[0] += x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator+(const taylor<Nvar,Ndeg>& t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<Nvar,Ndeg> tmp = t1;
  tmp += t2;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator-(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp *= -1;
  tmp[0] += x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator-(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp[0] -= x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator-(const taylor<Nvar,Ndeg>& t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<Nvar,Ndeg> tmp = t1;
  tmp -= t2;
  return tmp;
}

template<int Ndeg1, int Nvar2, int Ndeg2>
void taylor_compose(taylor<Nvar2, Ndeg2>& res, const taylor<1,Ndeg1>& t1, const taylor<Nvar2,Ndeg2>& t2)
{
  // Horner scheme
  res.zero();
  for (int i=Ndeg1;i>0;i--)
    {
      res[0] += t1[i];
      res *= t2;
    }
  res[0] += t1[0];
}

// Elementary functions

//Taylor series of 1/(a+x)
template<int N>
void inv_taylor(taylor<1,N>& t, const tnum_t &a)
{
  assert(a != 0 && "1/(a+x) not analytic at a = 0");
  t[0] = 1/a;
  for (int i=1;i<=N;i++)
    t[i] = -t[i-1]/a;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator/(const tnum_t &x, const taylor<Nvar,Ndeg>& t)
{
  taylor<1,Ndeg> tmp;
  inv_taylor(tmp,t[0]);
  tmp*=x;
  taylor<Nvar, Ndeg> res;
  t.compose(res,tmp);
  return res;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator/(const taylor<Nvar,Ndeg>& t, const tnum_t &x)
{
  taylor<Nvar,Ndeg> tmp = t;
  tmp *= 1/x;
  return tmp;
}

template<int Nvar, int Ndeg>
taylor<Nvar, Ndeg> operator/(const taylor<Nvar,Ndeg>&t1, const taylor<Nvar,Ndeg>& t2)
{
  taylor<1,Ndeg> tmp;
  inv_taylor(tmp,t2[0]);
  taylor<Nvar, Ndeg> res;
  t2.compose(res,tmp);
  res*=t1;
  return res;
}

// Evaluate the taylor series of exp(x0+x)=exp(x0)*exp(x)
template<int Ndeg>
void exp_taylor(taylor<1,Ndeg> &t, const tnum_t &x0)
{
  tnum_t ifac = 1;
  t[0] = exp(x0);
  for (int i=1;i<=Ndeg;i++)
    {
      ifac *= i;
      t[i] = t[0]/ifac;
    }
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> exp(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  exp_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

// Log series log(a+x) = log(1+x/a) + log(a)
template<int N>
void log_taylor(taylor<1,N> &t, const tnum_t &x0)
{
  assert(x0 != 0 && "log(x) not analytic at x = 0");
  t[0] = log(x0);
  tnum_t xn = x0;
  for (int i=1;i<=N;i++)
    {
      t[i] = (2*(i & 1)-1)/(i*xn);
      xn *= x0;
    }
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> log(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  log_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

/* Use that (x0+x)^a=x0^a*(1+x/x0)^a */
template<int N>
void pow_taylor(taylor<1,N>& t, const tnum_t &x0, const tnum_t &a)
{
  assert(x0 != 0 && "pow(x,a) not analytic at x = 0");
  t[0] = pow(x0,a);
  for (int i=1;i<=N;i++)
    t[i] = t[i-1]*(a-i+1)/(x0*i);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> pow(const taylor<Nvar,Ndeg> &t, const tnum_t &a)
{
  taylor<1,Ndeg> tmp;
  pow_taylor(tmp,t[0],a);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> sqrt(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  pow_taylor(tmp,t[0],0.5);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

//Integer exponent version is analytical at t[0] = 0
// TODO: speed this up a bit by using compose
template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> pow(const taylor<Nvar,Ndeg> &t, int n)
{
  if (n < 1)
    return pow(t,tnum_t(n));
  taylor<Nvar,Ndeg> res(1), tpow = t;
  while (n>0)
    {
      if (n & 1)
	res *= tpow;
      tpow *= tpow;
      n = n/2;
    }
  return res;
}

// Use that d/dx atan(x) = 1/(1 + x^2),
// Taylor expand in x^2 and integrate.
template<int Ndeg>
void atan_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  // Calculate taylor expansion of 1/(1+a^2+x)
  taylor<1,Ndeg> invt,x;
  inv_taylor(invt,1+a*a);
  //insert x = 2*a*x + x^2
  x = taylor<1,Ndeg>(0);
  if (Ndeg > 0)
    x[1] = 2*a;
  if (Ndeg > 1)
    x[2] = 1;
  taylor_compose(t,invt,x);
  //Integrate each term and set the constant
  t.integrate();
  t[0] = atan(a);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> atan(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  atan_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


/* Taylor expansion of exp(-(a+x)^2) =
   exp(-a^2-2a*x)*exp(-x^2)

 */
template<int Ndeg>
void gauss_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  exp_taylor(t,-a*a);
  t.stretch(-2*a);
  taylor<1,Ndeg> g;
  g.zero();
  g[0] = 1;
  for (int i=1;i<=Ndeg/2;i++)
    g[2*i] = -g[2*(i-1)]/i;
  t*=g;
}

// Use that d/dx erf(x) = 2/sqrt(pi)*exp(-x^2),
// Taylor expand in x^2 and integrate.
template<int Ndeg>
void erf_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  gauss_taylor(t,a);
  t*=2/sqrt(M_PI);
  //Integrate each term and set the constant
  t.integrate();
  t[0] = erf(a);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> erf(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  erf_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}


// hyperbolic arcsin function. d/dx asinh(x) = 1/sqrt(1+x^2)
// 1 + (a+x)^2 = 1+a^2 + 2ax + x^2
template<int Ndeg>
void asinh_taylor(taylor<1,Ndeg>& t, const tnum_t &a)
{
  taylor<1,Ndeg> tmp(1+a*a);
  if (Ndeg>0)
    tmp[1] = 2*a;
  if (Ndeg>1)
    tmp[2] = 1;
  t = pow(tmp,-0.5);
  t.integrate();
  t[0] = asinh(a);
}

template<int Nvar, int Ndeg>
taylor<Nvar,Ndeg> asinh(const taylor<Nvar,Ndeg> &t)
{
  taylor<1,Ndeg> tmp;
  asinh_taylor(tmp,t[0]);

  taylor<Nvar,Ndeg> res;
  t.compose(res,tmp);
  return res;
}

#endif
