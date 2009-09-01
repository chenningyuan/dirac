#ifndef FUN_HELPERS_H
#define FUN_HELPERS_H

#include "taylor.h"
#include "xc_fun.h"

template<int N>
inline void get_R(taylor<2,N> &res, const double d[2])
{
  res = d[XC_R];
  if (N>0)
    res[1+XC_R] = 1;
}

template<int N>
inline void get_R(taylor<5,N> &res, const double d[5])
{
  res = d[XC_R];
  if (N>0)
    res[1+XC_R] = XC_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_S(taylor<2,N> &res, const double d[2])
{
  res = d[XC_S];
  if (N>0)
    res[1+XC_S] = XC_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_S(taylor<5,N> &res, const double d[5])
{
  res = d[XC_S];
  if (N>0)
    res[1+XC_S] = XC_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_Z(taylor<5,N> &res, const double d[5])
{
  res = d[XC_Z];
  if (N>0)
    res[1+XC_Z] = XC_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_X(taylor<5,N> &res, const double d[5])
{
  res = d[XC_X];
  if (N>0)
    res[1+XC_X] = XC_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_Y(taylor<5,N> &res, const double d[5])
{
  res = d[XC_Y];
  if (N>0)
    res[1+XC_Y] = XC_FUNCTIONAL_SCALING;
}

template<int N>
inline void get_rhoa(taylor<2,N> &res, const double d[2])
{
  res = (d[XC_R]+d[XC_S])/2;
  if (N>0)
    {
      res[1+XC_R] = 0.5*XC_FUNCTIONAL_SCALING;
      res[1+XC_S] = 0.5*XC_FUNCTIONAL_SCALING;
    }
}
template<int N>
inline void get_rhob(taylor<2,N> &rhob, const double d[2])
{
  rhob = (d[XC_R]-d[XC_S])/2;
  if (N>0)
    {
      rhob[1+XC_R] =  0.5*XC_FUNCTIONAL_SCALING;
      rhob[1+XC_S] = -0.5*XC_FUNCTIONAL_SCALING;
    }
}

template<int N>
inline void get_rhoa(taylor<5,N> &rhoa, const double d[5])
{
  rhoa = (d[XC_R]+d[XC_S])/2;
  if (N>0)
    {
      rhoa[1+XC_R] = 0.5*XC_FUNCTIONAL_SCALING;
      rhoa[1+XC_S] = 0.5*XC_FUNCTIONAL_SCALING;
    }
}
template<int N>
inline void get_rhob(taylor<5,N> &rhob, const double d[5])
{
  rhob = (d[XC_R]-d[XC_S])/2;
  if (N>0)
    {
      rhob[1+XC_R] = 0.5*XC_FUNCTIONAL_SCALING;
      rhob[1+XC_S] = -0.5*XC_FUNCTIONAL_SCALING;
    }
}

template<int N>
inline void get_grada2(taylor<5,N> &ga2, const double d[5])
{
  ga2 = (d[XC_Z]+d[XC_X]+2*d[XC_Y])/4;
  if (N>0)
    {
      ga2[1+XC_Z] = 0.25*XC_FUNCTIONAL_SCALING; 
      ga2[1+XC_X] = 0.25*XC_FUNCTIONAL_SCALING; 
      ga2[1+XC_Y] = 0.5*XC_FUNCTIONAL_SCALING;  
    }
}

template<int N>
inline void get_gradab(taylor<5,N> &gab, const double d[5])
{
  gab = (d[XC_Z] - d[XC_X])/4;
  if (N>0)
    {
      gab[1+XC_Z] =  0.25*XC_FUNCTIONAL_SCALING;
      gab[1+XC_X] = -0.25*XC_FUNCTIONAL_SCALING;
      gab[1+XC_Y] = 0;  
    }
}

template<int N>
inline void get_gradb2(taylor<5,N> &gb2, const double d[5])
{
  gb2 = (d[XC_Z]+d[XC_X]-2*d[XC_Y])/4;
  if (N>0)
    {
      gb2[1+XC_Z] = 0.25*XC_FUNCTIONAL_SCALING; 
      gb2[1+XC_X] = 0.25*XC_FUNCTIONAL_SCALING; 
      gb2[1+XC_Y] = -0.5*XC_FUNCTIONAL_SCALING; 
    }
}


#endif
