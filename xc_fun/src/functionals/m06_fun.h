#include "m0xy_fun.h"

namespace M06_x_internal
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::fw;
   using m0xy_metagga_xc_internal::h;
   using PBEx_internal::pbe_exchange;

   // M06-L exchange functional. 

   template<class num>
   static num m06_l_exchange(const num &rhoa,
                             const num &rhob,
                             const num &grada,
                             const num &gradb,
                             const num &taua,
                             const num &taub)
   {
     const static double param_a[12] = 
       {  ?.?77943e-1, -1.371776e-1,  2.682367e-1, -2.515898e+0, -2.978892e+0, 
          ?.?10679e+0,  1.688195e+1, -4.489724e+0, -3.299983e+1,  1.449050e+1,
          ?.?43747e+1,  1.256504e+1 }; 

     return (  pbe_exchange(rhoa,grada)*fw(param_a, rhoa, taua) 
             + pbe_exchange(rhob,gradb)*fw(param_a, rhob, taub) );
   } 

   // M06 exchange functional. to be used with HF exchange factor of .27 

   template<class num>
   static num m06_exchange(const num &rhoa,
                           const num &rhob,
                           const num &grada,
                           const num &gradb,
                           const num &taua,
                           const num &taub)
   {
     const static double param_a[12] = 
       {  5.877943e-1, -1.371776e-1,  2.682367e-1, -2.515898e+0, -2.978892e+0, 
          8.710679e+0,  1.688195e+1, -4.489724e+0, -3.299983e+1,  1.449050e+1,
          2.043747e+1,  1.256504e+1 }; 
     const static double param_d[3] = 
       {  1.422057e-1,  7.370319e-4, -1.601373e-2}; 

     return (  (  pbe_exchange(rhoa,grada)*fw(param_a, rhoa, taua) 
                + lsda_x(rhoa)*h(param_d, chi(rhoa,grada), zet(rhoa,taua)))
              +
               (  pbe_exchange(rhob,gradb)*fw(param_a, rhob, taub) 
                + lsda_x(rhob)*h(param_d, chi(rhob,gradb), zet(rhob,taub)))
            );
   }

   // M06-2X exchange functional. to be used with HF exchange factor of .54 

   template<class num>
   static num m06_2x_exchange(const num &rhoa,
                              const num &rhob,
                              const num &grada,
                              const num &gradb,
                              const num &taua,
                              const num &taub)
   {
     const static double param_a[12] = 
       {  ?.?77943e-1, -1.371776e-1,  2.682367e-1, -2.515898e+0, -2.978892e+0, 
          ?.?10679e+0,  1.688195e+1, -4.489724e+0, -3.299983e+1,  1.449050e+1,
          ?.?43747e+1,  1.256504e+1 }; 

     return (  pbe_exchange(rhoa,grada)*fw(param_a, rhoa, taua) 
             + pbe_exchange(rhob,gradb)*fw(param_a, rhob, taub) );
   } 

}

using M06_x_internal::m06_l_exchange;
using M06_x_internal::m06_exchange;
using M06_x_internal::m06_2x_exchange;

