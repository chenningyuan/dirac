#ifndef FUN_HELPERS_H
#define FUN_HELPERS_H

#include "taylor.h"
#include "xc_fun.h"

template<int N, int Nvar>
taylor<Nvar,N> get_rhoa(const taylor<Nvar,N> &R, 
			const taylor<Nvar,N> &S)
{
  return 0.5*(R+S);
}

template<int N, int Nvar>
taylor<Nvar,N> get_rhob(const taylor<Nvar,N> &R, 
			const taylor<Nvar,N> &S)
{
  return 0.5*(R-S);
}

template<int N, int Nvar>
taylor<Nvar,N> get_grada2(const taylor<Nvar,N> &Z,
			  const taylor<Nvar,N> &X,
			  const taylor<Nvar,N> &Y)
{
  return 0.25*(Z+X+2*Y);
}

template<int N, int Nvar>
taylor<Nvar,N> get_gradab(const taylor<Nvar,N> &Z,
			  const taylor<Nvar,N> &X,
			  const taylor<Nvar,N> &Y)
{
  return 0.25*(Z-X);
}

template<int N, int Nvar>
taylor<Nvar,N> get_gradb2(const taylor<Nvar,N> &Z,
			  const taylor<Nvar,N> &X,
			  const taylor<Nvar,N> &Y)
{
  return 0.25*(Z+X-2*Y);
}

#endif
