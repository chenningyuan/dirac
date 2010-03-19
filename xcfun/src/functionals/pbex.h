#ifndef PBEX_H
#define PBEX_H
#include "pw9xx.h"

namespace pbex
{
  // enhancement factor F(S), common to PBEx and REVPBEx
  template<class num>
    static num enhancement(const double &R, 
			   const num &rho,
			   const num &grad)
    {
      using pw91_like_x_internal::S;
      
      // double mu = 0.066725*M_PI*M_PI/3.0;
      // ulfek: mu from Daresbury implementation
      static const double mu = 0.2195149727645171;
      num st = S(rho,grad);
      num t1 = 1 + mu*st*st/R;   
      return 1 + R - R/t1; 
    }
}

#endif
