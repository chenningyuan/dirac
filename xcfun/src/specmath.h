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

template<class T, int Nvar, int Ndeg>
taylor<T,Nvar,Ndeg> preexpand(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &),
			      const taylor<T,Nvar,Ndeg> &t)
{
  taylor<T,Nvar,Ndeg> res;
  t.compose(res,f(taylor<T,1,Ndeg>(t[0],0)));
  return res;
}

template<class T, int Nvar, int Ndeg>
taylor<T,Nvar,Ndeg> preexpand(taylor<T,1,Ndeg> (*f)(const taylor<T,1,Ndeg> &,
						    const T *),
			      const taylor<T,Nvar,Ndeg> &t,
			      const T *param)
{
  taylor<T,Nvar,Ndeg> res;
  t.compose(res,f(taylor<T,1,Ndeg>(t[0],0),param));
  return res;
}

#endif
