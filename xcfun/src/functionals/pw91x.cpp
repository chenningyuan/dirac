#include "functional.h"
#include "pw9xx.h"

template<class num>
static num energy(const densvars<num> &d)
{
  const static double param_AB[6] = 
    { 0.19645, 7.7956, 0.2743, 0.1508, 100.0, 0.004};
  using pw91_like_x_internal::prefactor;
  using pw91_like_x_internal::pw91xk_enhancement;
  return prefactor(d.a,d.gaa)*pw91xk_enhancement(param_AB,d.a,d.gaa)
    + prefactor(d.b,d.gbb)*pw91xk_enhancement(param_AB,d.b,d.gbb);
}

SETUP_FUNCTIONAL(functional &f)
{
  f.describe("pw91x",XC_LDA,
	     "Perdew-Wang 1991 GGA Exchange Functional",
	     "Perdew-Wang 1991 GGA Exchange Functional\n"
	     "J. P. Perdew, J. A. Chevary, S. H. Vosko, "
	     "K. A. Jackson, M. R. Pederson, and C. Fiolhais, "
	     "Phys. Rev. B 46, 6671 (1992)\n"
	     "Implemented by Andre Gomes.\n"); 
  SET_LDA_ENERGY_FUNCTION(f,energy);
}
