#ifndef POLYMUL_H
#define POLYMUL_H
#include <cassert>

/*
  By Ulf Ekstrom March-June 2009.

  For 5 var 5 deg template taylor multiplication is 
  about 25% slower than explicit table.
 */

/*
  Multiply polynomials p1 and p2, both polynomials of nvar variables,
  of degree n1 and n2 respectively.  The result is a polynomial of
  degree n1+n2, which is _added_ to dst.

  A degree N polynomial of K variables has
  /N + K\
  \  N  /   terms.

  A polynomial P[K,N] is here organized as a sequence of monomials in
  K variables:

  M[K,0], M[K,1] .. M[K,N]  (for example 1, x y, x^2 xy y^2, ... )

  However, for each monomial all terms are of the same degree and the
  exponents of one variable can be inferred from the other
  exponents. We therefore have

  M[K,N] ~ P[K-1,N] (example: x^2 xy y^2 ~ 1 y y^2).

  To multiply monomials we can then use the ordinary polymul()
  function with lower nvar value.
*/


// Template to evaluate binomial coefficients at compile time.
template <int n, int k>
struct poly_bin
{
    enum { value = poly_bin<n-1,k-1>::value + poly_bin<n-1,k>::value };
};

template <int n>
struct poly_bin<n,n>
{
    enum { value = 1 };
};

template <int n>
struct poly_bin<n,0>
{
    enum { value = 1 };
};

template <int k>
struct poly_bin<0,k>
{
    enum { value = 1 };
};

template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
class polynomial_multiplier
{
 public:
  // _add_ the product between p1 and p2 to dst.
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  { 
    polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2-1>::mul(dst,p1,p2);
    polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2>::mul_monomial(dst,p1,p2+poly_bin<Nvar+Ndeg2-1,Ndeg2-1>::value);
  }
  // m2 is a monomial in Nvar variables and of order Ndeg2
  static void mul_monomial(numtype dst[], const numtype p1[], const numtype m2[])
  { 
    polynomial_multiplier<numtype,Nvar-1,Ndeg1,Ndeg2>::mul(dst+poly_bin<Nvar+Ndeg1+Ndeg2-1,Ndeg1+Ndeg2-1>::value,
							   p1 +poly_bin<Nvar+Ndeg1-1,Ndeg1-1>::value,
							   m2);
    polynomial_multiplier<numtype,Nvar,Ndeg1-1,Ndeg2>::mul_monomial(dst,p1,m2);
  }
};

template<class numtype, int Ndeg1, int Ndeg2>
  class polynomial_multiplier<numtype, 1, Ndeg1, Ndeg2>
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<=Ndeg1;i++)
      for (int j=0;j<=Ndeg2;j++)
	dst[i+j] += p1[i]*p2[j];
  }
  static void mul_monomial(numtype dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=0;i<Ndeg1;i++)
      dst[i+Ndeg2] += m2[0]*p1[i];
  }
};

template<class numtype, int Nvar, int Ndeg2>
  class polynomial_multiplier<numtype, Nvar, 0, Ndeg2>
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<poly_bin<Nvar+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*p2[i];
  }
  static void mul_monomial(numtype dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=poly_bin<Nvar+Ndeg2-1,Ndeg2-1>::value;i<poly_bin<Nvar+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*m2[i-poly_bin<Nvar+Ndeg2-1,Ndeg2-1>::value];
  }
};

template<class numtype, int Ndeg2>
  class polynomial_multiplier<numtype, 1, 0, Ndeg2>
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<poly_bin<1+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*p2[i];
  }
  static void mul_monomial(numtype dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=poly_bin<1+Ndeg2-1,Ndeg2-1>::value;i<poly_bin<1+Ndeg2,Ndeg2>::value;i++)
      dst[i] += p1[0]*m2[i-poly_bin<1+Ndeg2-1,Ndeg2-1>::value];
  }
};

template<class numtype, int Nvar, int Ndeg1>
  class polynomial_multiplier<numtype, Nvar, Ndeg1, 0>
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<poly_bin<Nvar+Ndeg1,Ndeg1>::value;i++)
      dst[i] += p1[i]*p2[0];
  }
  static void mul_monomial(numtype dst[], const numtype p1[], const numtype m2[])
  {
    for (int i=0;i<poly_bin<Nvar+Ndeg1,Ndeg1>::value;i++)
      dst[i] += p1[i]*m2[0];
  }
};

template<class numtype, int Nvar>
  class polynomial_multiplier<numtype, Nvar, 0, 0>
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    dst[0] += p1[0]*p2[0];
  }
  static void mul_monomial(numtype dst[], const numtype p1[], const numtype m2[])
  {
    dst[0] += p1[0]*m2[0];
  }
};


// Like polymul but truncates the result to Ndeg1 order order.
// Ndeg2 _must_ be less or equal to Ndeg1, otherwise the recursion
// will not terminate.
template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
class taylor_multiplier
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    polynomial_multiplier<numtype,Nvar,Ndeg1-Ndeg2,Ndeg2>::mul_monomial(dst,p1,p2+poly_bin<Nvar+Ndeg2-1,Ndeg2-1>::value);   
    taylor_multiplier<numtype,Nvar,Ndeg1,Ndeg2-1>::mul(dst,p1,p2);
  }
};

template<class numtype, int Nvar, int Ndeg1>
  class taylor_multiplier<numtype, Nvar,Ndeg1,0>
{
 public:
  static void mul(numtype dst[], const numtype p1[], const numtype p2[])
  {
    for (int i=0;i<poly_bin<Ndeg1+Nvar,Ndeg1>::value;i++)
      dst[i] += p1[i]*p2[0];
  }
};


template<class numtype, int Nvar, int Ndeg>
  class polynomial
{
 public:
  polynomial(void) {}
  polynomial(numtype c0) 
    { 
      c[0] = c0;
      for (int i=1;i<this->size();i++)
	c[i] = 0;
    }
  numtype operator[](int i) const
  {
    assert(i>=0);
    assert(i<this->size());
    return c[i];
  }
  numtype &operator[](int i)
  {
    assert(i>=0);
    assert(i<this->size());
    return c[i];
  }
  int size(void) const { return poly_bin<Nvar+Ndeg,Ndeg>::value; }
  void zero(void)
  {
    for (int i=0;i<size();i++)
      c[i] = 0;
  }
  // This is a _very slow_ function to get the exponents
  // of a particular term. 
  static void exponents(int term, int exponents[])
  {
    if (Nvar == 1)
      {
	exponents[0] = term;
	return;
      }
    for (int i=0;i<Nvar;i++)
      exponents[i] = 0;
    assert(term >= 0);
    if (term >= poly_bin<Nvar+Ndeg,Ndeg>::value)
      {
	assert(0 && "term < poly_bin<Nvar+Ndeg,Ndeg>::value");
      }
    for (int i=0;i<term;i++)
      polynomial<numtype,Nvar,Ndeg>::next_exponents(Nvar,exponents);
  }
  // Return the index of the term with certain exponents
  static int term_index(const int exponents[])
  {
    int N = 0;
    for (int i=0;i<Nvar;i++)
      N += exponents[i];
    int i = 0, idx = 0;
    N--;
    while (N >= 0)
      {
	idx += polynomial<numtype,Nvar,Ndeg>::polylen(Nvar-i,N);
	N -= exponents[i];
	i++;
      }
    return idx;
  }

  numtype c[poly_bin<Nvar+Ndeg,Ndeg>::value];

 protected:
  static void next_exponents(int nvar, int m[])
  {
    int k = 0;
    for (int i=0;i<nvar-1;i++)
      k += m[i];
    if (k == 0)
      {
	m[0] = m[nvar-1] + 1;
	m[nvar-1] = 0;
	return;
      }
    if (m[nvar-2] > 0)
      {
	m[nvar-1]++;
	m[nvar-2]--;
      }
    else
      {
	next_exponents(nvar-1,m);
	for (int i=nvar-2;i>=0;i--)
	  {
	    if (m[i] > 0)
	      {
		m[i] += m[nvar-1];
		break;
	      }
	  }
	m[nvar-1] = 0;
      }
  }
  // = binomial(nvar+ndeg,ndeg), but this one
  // can be evaluated at run time.
  static int polylen(int nvar, int ndeg)
  {
    int len = 1;
    for (int k=1;k<=nvar;k++)
      {
	len *= ndeg + k;
	len /= k;
      }
    return len;
  }
};

template<class numtype, int Nvar, int Ndeg1, int Ndeg2>
  void polymul(polynomial<numtype, Nvar,Ndeg1+Ndeg2> &dst,
	       const polynomial<numtype, Nvar,Ndeg1> &p1,
	       const polynomial<numtype, Nvar,Ndeg2> &p2)
{
  polynomial_multiplier<numtype,Nvar,Ndeg1,Ndeg2>::mul(dst.c,p1.c,p2.c);
}

template<class numtype, int Nvar, int Ndeg>
  void taylormul(polynomial<numtype, Nvar,Ndeg> &dst,
		 const polynomial<numtype, Nvar,Ndeg> &p1,
		 const polynomial<numtype, Nvar,Ndeg> &p2)
{
  taylor_multiplier<numtype,Nvar,Ndeg,Ndeg>::mul(dst.c,p1.c,p2.c);
};


#endif
