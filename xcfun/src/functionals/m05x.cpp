#include "functional.h"
#include "pbex.h"
#include "m0xy_fun.h"

// M05 exchange functional. to be used with HF exchange factor of .28 

template<class num>
static num energy (const densvars<num> &d)
{
   using m0xy_metagga_xc_internal::fw;

   const static double param_a[12] =
     {  1.000000e+00,  8.151000e-02, -4.395600e-01, -3.224220e+00,  2.018190e+00,
        8.794310e+00, -2.950000e-03,  9.820290e+00, -4.823510e+00, -4.817574e+01,
        3.648020e+00,  3.402248e+01 };

   return (  pbex::energy_pbe_ab(pbex::R_pbe,d.a,d.gaa)*fw(param_a, d.a, d.taua)
           + pbex::energy_pbe_ab(pbex::R_pbe,d.b,d.gbb)*fw(param_a, d.b, d.taub) );
}

SETUP_FUNCTIONAL(functional &f)
{
  f.describe("m05x",XC_MGGA,
	     "M05 exchange",
             "M05 Meta-GGA Exchange Functional\n"
             "Y Zhao, N. E. Schultz and D. G. Truhlar, J. Chem. Theory Comput. 2, 364 (2006)\n"
             "Implemented by Andre Gomes\n");

  SET_MGGA_ENERGY_FUNCTION(f,energy);

}
