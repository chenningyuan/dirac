#ifndef AUTO_FUNCTIONALS_H
#define AUTO_FUNCTIONALS_H

#include "taylor.h"

// The AUTO_FUNCTIONAL_SCALING number is a scale factor
// applied to each derivative order. For example a factor
// of 0.5 would scale first derivatives by 0.5, second 
// derivatives by 0.5^2  = 0.25, and so on. This can be used
// to avoid very large numbers for high derivatives.
#define AUTO_FUNCTIONAL_SCALING 1.0

template<int N>
taylor<2,N> vwn_correlation(tnum_t R0);

template<int N>
taylor<2,N> slater_exchange(tnum_t R0);

template<int N>
taylor<5,N> lyp_correlation(tnum_t R0, tnum_t Z0);

template<int N>
taylor<5,N> becke_exchange(tnum_t R0, tnum_t Z0);

#endif
