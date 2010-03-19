#include "functional.h"
#include "pbex.h"

template<class num>
static num energy(const densvars<num> &d)
{
  static const double R = 1.245;
  using pw91_like_x_internal::prefactor;
  return prefactor(d.a,d.gaa)*pbex::enhancement(R,d.a,d.gaa) 
    + prefactor(d.b,d.gbb)*pbex::enhancement(R,d.b,d.gbb);
}


SETUP_FUNCTIONAL(functional &f)
{
  f.describe("revpbex",XC_GGA,
	     "Revised PBE Exchange Functional",
	     "Revised PBE Exchange Functional\n"
	     "Y. Zhang and W., Phys. Rev. Lett 80, 890 (1998)\n"
	     "Implemented by Ulf Ekstrom and Andre Gomes\n");
  SET_GGA_ENERGY_FUNCTION(f,energy);
}
