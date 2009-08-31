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
taylor<2,N> vwn_correlation(const double d[2]);

template<int N>
taylor<2,N> slater_exchange(const double d[2]);

template<int N>
taylor<5,N> lyp_correlation(const double d[5]);

template<int N>
taylor<5,N> becke_exchange(const double d[5]);

template<int N>
taylor <5,N> pbe_exchange(const double d[5]);

template<int N>
taylor <5,N> pbe_correlation(const double d[5]);


// Helpers

template<int N>
inline void get_R(taylor<2,N> &res, const double d[2])
{
  res = d[0];
  if (N>0)
    res[1] = 1;
}

template<int N>
inline void get_R(taylor<5,N> &res, const double d[5])
{
  res = d[0];
  if (N>0)
    res[1] = AUTO_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_S(taylor<2,N> &res, const double d[2])
{
  res = d[1];
  if (N>0)
    res[2] = AUTO_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_S(taylor<5,N> &res, const double d[5])
{
  res = d[1];
  if (N>0)
    res[2] = AUTO_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_Z(taylor<5,N> &res, const double d[5])
{
  res = d[2];
  if (N>0)
    res[3] = AUTO_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_X(taylor<5,N> &res, const double d[5])
{
  res = d[3];
  if (N>0)
    res[4] = AUTO_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_Y(taylor<5,N> &res, const double d[5])
{
  res = d[4];
  if (N>0)
    res[5] = AUTO_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_res(taylor<2,N> &res, const double d[2])
{
  res = (d[0]+d[1])/2;
  if (N>0)
    {
      res[1] = 0.5*AUTO_FUNCTIONAL_SCALING;
      res[2] = 0.5*AUTO_FUNCTIONAL_SCALING;
    }
}
template<int N>
inline void get_rhob(taylor<2,N> &rhob, const double d[2])
{
  rhob = (d[0]-d[1])/2;
  if (N>0)
    {
      rhob[1] = 0.5*AUTO_FUNCTIONAL_SCALING;
      rhob[2] = -0.5*AUTO_FUNCTIONAL_SCALING;
    }
}

template<int N>
inline void get_rhoa(taylor<5,N> &rhoa, const double d[5])
{
  rhoa = (d[0]+d[1])/2;
  if (N>0)
    {
      rhoa[1] = 0.5*AUTO_FUNCTIONAL_SCALING;
      rhoa[2] = 0.5*AUTO_FUNCTIONAL_SCALING;
    }
}
template<int N>
inline void get_rhob(taylor<5,N> &rhob, const double d[5])
{
  rhob = (d[0]-d[1])/2;
  if (N>0)
    {
      rhob[1] = 0.5*AUTO_FUNCTIONAL_SCALING;
      rhob[2] = -0.5*AUTO_FUNCTIONAL_SCALING;
    }
}

template<int N>
inline void get_grada2(taylor<5,N> &ga2, const double d[5])
{
  ga2 = (d[2]+d[3]+2*d[4])/4;
  if (N>0)
    {
      ga2[3] = 0.25*AUTO_FUNCTIONAL_SCALING; // Z
      ga2[4] = 0.25*AUTO_FUNCTIONAL_SCALING; // X
      ga2[5] = 0.5*AUTO_FUNCTIONAL_SCALING;  // Y
    }
}

template<int N>
inline void get_gradab(taylor<5,N> &gab, const double d[5])
{
  gab = (d[2] - d[3])/4;
  if (N>0)
    {
      gab[3] = 0.25*AUTO_FUNCTIONAL_SCALING; // Z
      gab[4] = -0.25*AUTO_FUNCTIONAL_SCALING; // X
      gab[5] = 0;  // Y
    }
}

template<int N>
inline void get_gradb2(taylor<5,N> &gb2, const double d[5])
{
  gb2 = (d[2]+d[3]-2*d[4])/4;
  if (N>0)
    {
      gb2[3] = 0.25*AUTO_FUNCTIONAL_SCALING; // Z
      gb2[4] = 0.25*AUTO_FUNCTIONAL_SCALING; // X
      gb2[5] = -0.5*AUTO_FUNCTIONAL_SCALING;  // Y
    }
}


#endif
