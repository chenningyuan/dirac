#ifndef AUTO_FUNCTIONALS_H
#define AUTO_FUNCTIONALS_H

#include "taylor.h"

// The AUTO_FUNCTIONAL_SCALING number is a scale factor
// applied to each derivative order. For example a factor
// of 0.5 would scale first derivatives by 0.5, second 
// derivatives by 0.5^2  = 0.25, and so on. This can be used
// to avoid very large numbers for high derivatives.
#define AUTO_FUNCTIONAL_SCALING 1.0

//Functionals:

template<int N>
taylor<2,N> vwn_correlation(tnum_t R0);

template<int N>
taylor<2,N> slater_exchange(tnum_t R0);

template<int N>
taylor<5,N> lyp_correlation(tnum_t R0, tnum_t Z0);

template<int N>
taylor<5,N> becke_exchange(tnum_t R0, tnum_t Z0);

template<int N>
taylor <5,N> pbe_exchange(tnum_t R0, tnum_t Z0);

template<int N>
taylor <5,N> pbe_correlation(tnum_t R0, tnum_t Z0);


// Helpers

template<int N>
inline void get_rhoa(taylor<5,N> &rhoa, tnum_t R0)
{
  rhoa = R0/2;
  if (N>0)
    {
      rhoa[1] = 0.5;
      rhoa[2] = 0.5;
    }
}
template<int N>
inline void get_rhob(taylor<5,N> &rhob, tnum_t R0)
{
  rhob = R0/2;
  if (N>0)
    {
      rhob[1] = 0.5;
      rhob[2] = -0.5;
    }
}

template<int N>
inline void get_grada2(taylor<5,N> &ga2, tnum_t Z0)
{
  ga2 = Z0/4;
  if (N>0)
    {
      ga2[3] = 0.25; // Z
      ga2[4] = 0.25; // X
      ga2[5] = 0.5;  // Y
    }
}

template<int N>
inline void get_gradab(taylor<5,N> &gab, tnum_t Z0)
{
  gab = Z0/4;
  if (N>0)
    {
      gab[3] = 0.25; // Z
      gab[4] = -0.25; // X
      gab[5] = 0;  // Y
    }
}

template<int N>
inline void get_gradb2(taylor<5,N> &gb2, tnum_t Z0)
{
  gb2 = Z0/4;
  if (N>0)
    {
      gb2[3] = 0.25; // Z
      gb2[4] = 0.25; // X
      gb2[5] = -0.5;  // Y
    }
}


#endif
