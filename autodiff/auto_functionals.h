#ifndef AUTO_FUNCTIONALS_H
#define AUTO_FUNCTIONALS_H

#include "taylor.h"

// The AUTO_FUNCTIONAL_SCALING number is a scale factor
// applied to each derivative order. For example a factor
// of 0.5 would scale first derivatives by 0.5, second 
// derivatives by 0.5^2  = 0.25, and so on. This can be used
// to avoid very large numbers for high derivatives.
#define AUTO_FUNCTIONAL_SCALING 1.0

// Helper functions:

template<int N>
void get_rhoa(taylor<5,N> &rhoa, tnum_t R0);
template<int N>
void get_rhob(taylor<5,N> &rhob, tnum_t R0);

template<int N>
void get_grada2(taylor<5,N> &ga2, tnum_t Z0);
template<int N>
void get_gradab(taylor<5,N> &gab, tnum_t Z0);
template<int N>
void get_gradb2(taylor<5,N> &gb2, tnum_t Z0);

//Functionals:

template<int N>
taylor<2,N> vwn_correlation(tnum_t R0);

template<int N>
taylor<2,N> slater_exchange(tnum_t R0);

template<int N>
taylor<5,N> lyp_correlation(tnum_t R0, tnum_t Z0);

template<int N>
taylor<5,N> becke_exchange(tnum_t R0, tnum_t Z0);


#endif
