#ifndef SPECMATH_H
#define SPECMATH_H

// Some math-related functions useful in many places
#include "taylor.h"

template<class T>
T pow2(const T &t)
{
  return t*t;
}

template<class T>
T pow3(const T &t)
{
  return t*t*t;
}

template<class T, class T2>
T poly(const T &x, int ndeg, const T2 coeffs[])
{
  // Horner rule
  T res = coeffs[--ndeg];
  while (ndeg)
    {
      res *= x;
      res += coeffs[--ndeg];
    }
  return res;
}


// Rational function P(x)/Q(x), for P and Q polynomials
// of degree np and nq with coefficients pc and qc.
template<class T, class T2>
T ratfun(const T &x, int np, const T2 pc[], int nq, const T2 qc[])
{
  return poly(x,np,pc)/poly(x,nq,qc);
}

template<class T, int Nvar, int Ndeg, class T2>
taylor<T,Nvar,Ndeg> ratfun(const taylor<T,Nvar,Ndeg> &x, 
			   int np, const T2 pc[], int nq, const T2 qc[])
{
  taylor<T,1,Ndeg> t(x[0],0);
  return x.compose(poly(t,np,pc)/poly(t,nq,qc));
}


template<class T, int Nvar, int Ndeg>
taylor<T,Nvar,Ndeg> preexpand(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &),
			      const taylor<T,Nvar,Ndeg> &t)
{
  return t.compose(f(taylor<T,1,Ndeg>(t[0],0)));
}

template<class T, int Nvar, int Ndeg>
taylor<T,Nvar,Ndeg> preexpand(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &,
						    const T *),
			      const taylor<T,Nvar,Ndeg> &t,
			      const T *param)
{
  return t.compose(f(taylor<T,1,Ndeg>(t[0],0),param));
}



#endif
