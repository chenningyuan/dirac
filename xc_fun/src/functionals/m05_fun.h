#include "m0xy_fun.h"


// Exchange functionals for the different M05 flavors 

namespace M05_x_internal
{
   using m0xy_metagga_xc_internal::fw;
   using PBEx_internal::pbe_exchange;

   // M05 exchange functional. to be used with HF exchange factor of .28 

   template<class num>
   static num m05_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  1.000000e+00,  8.151000e-02, -4.395600e-01, -3.224220e+00,  2.018190e+00, 
          8.794310e+00, -2.950000e-03,  9.820290e+00, -4.823510e+00, -4.817574e+01,
          3.648020e+00,  3.402248e+01 }; 

     return (  pbe_exchange(d.a,d.gaa)*fw(param_a, d.a, d.taua) 
             + pbe_exchange(d.b,d.gbb)*fw(param_a, d.b, d.taub) );
   }

   // M05-2X exchange functional. to be used with HF exchange factor of .56 

   template<class num>
   static num m05_2x_exchange(const densvars<num> &d)
   {
     const static double param_a[12] = 
       {  1.000000e+00, -5.683300e-01, -1.300570e+00,  5.500700e+00,  9.064020e+00, 
         -3.221075e+01, -2.373298e+01,  7.022996e+01,  2.988614e+01, -6.025778e+01,
         -1.322205e+01,  1.523694e+01 }; 

     return (  pbe_exchange(d.a,d.gaa)*fw(param_a, d.a, d.taua) 
             + pbe_exchange(d.b,d.gbb)*fw(param_a, d.b, d.taub) );
   }

}

using M05_x_internal::m05_exchange;
using M05_x_internal::m05_2x_exchange;


namespace M05_c_internal
{
   using pw91_like_x_internal::chi;
   using m0xy_metagga_xc_internal::zet;
   using m0xy_metagga_xc_internal::m05_c_anti;
   using m0xy_metagga_xc_internal::ueg_c_anti;
   using m0xy_metagga_xc_internal::m05_c_para;
   using m0xy_metagga_xc_internal::ueg_c_para;

   // M05 correlation functional. 

   template<class num>
   static num m05_correlation(const densvars<num> &d)
   {
     // parameters for anti-parallel spin contributions
     const static double param_c_anti[5] =
       {  1.000000e+00,  3.785690e+00, -1.415261e+01, -7.465890e+00,  1.794491e+00 };

     // parameters for parallel spin contributions
     const static double param_c_para[5] =
       {  1.000000e+00,  3.773440e+00, -2.604463e+01,  3.069913e+01, -9.226950e+00 };

     num chi_a = chi(d.a, d.gaa);
     num chi_b = chi(d.b, d.gbb);
     num zet_a = zet(d.a, d.taua);
     num zet_b = zet(d.b, d.taub);

     num Ec_ab = ueg_c_anti(d)  *m05_c_anti(param_c_anti[5],chi_a,chi_b); 
     num Ec_aa = ueg_c_para(d.a)*m05_c_para(param_c_para[5],chi_a,zet_a);
     num Ec_bb = ueg_c_para(d.b)*m05_c_para(param_c_para[5],chi_b,zet_b);

     return Ec_ab + Ec_aa + Ec_bb;
   }

   // M05-2X correlation functional. 

   template<class num>
   static num m05_2x_correlation(const densvars<num> &d)
   {
     // parameters for anti-parallel spin contributions
     const static double param_c_anti[5] =
       {  1.000000e+00,  1.092970e+00, -3.791710e+00,  2.828100e+00, -1.058909e+01 };

     // parameters for parallel spin contributions
     const static double param_c_para[5] =
       {  1.000000e+00, -3.054300e+00,  7.618540e+00,  1.476650e+00, -1.192365e+01 };

     num chi_a = chi(d.a, d.gaa);
     num chi_b = chi(d.b, d.gbb);
     num zet_a = zet(d.a, d.taua);
     num zet_b = zet(d.b, d.taub);

     num Ec_ab = ueg_c_anti(d)  *m05_c_anti(param_c_anti[5],chi_a,chi_b); 
     num Ec_aa = ueg_c_para(d.a)*m05_c_para(param_c_para[5],chi_a,zet_a);
     num Ec_bb = ueg_c_para(d.b)*m05_c_para(param_c_para[5],chi_b,zet_b);

     return Ec_ab + Ec_aa + Ec_bb;
   }
}

using M05_c_internal::m05_correlation;
using M05_c_internal::m05_2x_correlation;

