#ifndef AUTO_FUNCTIONALS_H
#define AUTO_FUNCTIONALS_H

#include "taylor.h"

template<int N>
taylor<2,N> vwn_correlation(double R0);

template<int N>
taylor<2,N> slater_exchange(double R0);

template<int N>
taylor<5,N> lyp_correlation(double R0, double Z0);

template<int N>
taylor<5,N> becke_exchange(double R0, double Z0);

#endif
