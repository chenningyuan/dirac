#include "m0xy_fun.h"

namespace M05_x_internal
{
   using m0xy_metagga_xc_internal::fw;
   using PBEx_internal::pbe_exchange;

   // M05 exchange functional. to be used with HF exchange factor of .28 

   template<class num>
   static num m05_exchange(const num &rhoa,
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

     return (  pbe_exchange(rhoa,grada)*fw(param_a, rhoa, taua) 
             + pbe_exchange(rhob,gradb)*fw(param_a, rhob, taub) );
   }

   // M05-2X exchange functional. to be used with HF exchange factor of .56 

   template<class num>
   static num m05_2x_exchange(const num &rhoa,
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

     return (  pbe_exchange(rhoa,grada)*fw(param_a, rhoa, taua) 
             + pbe_exchange(rhob,gradb)*fw(param_a, rhob, taub) );
   }

}

using M05_x_internal::m05_exchange;
using M05_x_internal::m05_2x_exchange;

