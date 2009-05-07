#include <cassert>
#include "taylor.h"
#include <iostream>
using namespace std;

// Fortran interface to taylor_len
extern "C"
int taylor_len_(int *nvar, int *ndeg)
{
  return taylor_len(*nvar,*ndeg);
}

int taylor_index(int nvar, const int exponents[])
{
  int N = 0;
  for (int i=0;i<nvar;i++)
    N += exponents[i];
  int i = 0, idx = 0;
  N--;
  while (N >= 0)
    {
      idx += taylor_len(nvar-i,N);
      N -= exponents[i];
      i++;
    }
  return idx;
}

// Fortran version
extern "C"
int taylor_index_(const int *nvar, const int exponents[])
{
  return 1+taylor_index(*nvar,exponents);
}

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
void polymul(tnum_t dst[], int nvar,
	     const tnum_t p1[], int n1,
	     const tnum_t p2[], int n2)
{
  assert(nvar >= 1);
  assert(n1 >= 0);
  assert(n2 >= 0);
  if (nvar == 1)
    {
      for (int i=0;i<=n1;i++)
	for (int j=0;j<=n2;j++)
	  dst[i+j] += p1[i]*p2[j];
    }
  else if (n1 == 0)
    {
      for (int i=0;i<=taylor_len(nvar,n2);i++)
	dst[i] += p1[0]*p2[i];
    }
  else if (n2 == 0)
    {
      for (int i=0;i<=taylor_len(nvar,n1);i++)
	dst[i] += p2[0]*p1[i];
    }
  else
    {
      for (int i=0;i<=n1;i++)
	for (int j=0;j<=n2;j++) // dst monomial i+j.
	  polymul(dst+taylor_len(i+j-1+nvar,i+j-1),nvar-1,
		  p1+taylor_len(nvar,i-1),i,
		  p2+taylor_len(nvar,j-1),j);
    }
}


/* Multiply polynomials p1 and p2, truncating all terms of degree
   higher than ndeg. Dst is zeroed before accumulating the product. */
void taylor_mul(tnum_t dst[], int nvar, int ndeg,
		const tnum_t p1[], const tnum_t p2[])
{
  for (int i=0;i<taylor_len(nvar,ndeg);i++)
    dst[i] = 0;
  if (nvar == 1)
    {
      for (int i=0;i<=ndeg;i++)
	for (int j=0;j<=ndeg-i;j++)
	  dst[i+j] += p1[i]*p2[j];
    }
  else
    {
      for (int i=0;i<=ndeg;i++)
	for (int j=0;j<=ndeg-i;j++)
	  polymul(dst+taylor_len(nvar,i+j-1),nvar-1,
		  p1+taylor_len(nvar,i-1),i,
		  p2+taylor_len(nvar,j-1),j);
    }
}

