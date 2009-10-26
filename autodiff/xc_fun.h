#ifndef XC_FUN
#define XC_FUN

/* 
   Main entry point to DFT exchange-correlation functionals 
   by ulfek 2009 
*/

// Note that these parameters have be in sync with those in
// xc_fun_module.f90 (Fortran version)
#define XC_MAX_ORDER 5 // Update xc_eval* if this changes
#define XC_FUNCTIONAL_SCALING 1.0 // Scale derivatives by this factor

enum xc_vars
  {
    XC_R = 0, // R = \rho_\alpha + \rho_\beta
    XC_S,     // S = \rho_\alpha - \rho_\beta
    XC_Z,     // Z = \nabla R \cdot \nabla R
    XC_X,     // X = \nabla S \cdot \nabla S
    XC_Y      // Y = \nabla R \cdot \nabla S
  };

enum xc_funparams
  {
    XC_SLATER_EXCHANGE = 0,
    XC_VWN_CORRELATION,
    XC_BECKE_EXCHANGE,
    XC_LYP_CORRELATION,
    XC_PBE_EXCHANGE,
    XC_PBE_CORRELATION,
    XC_OLD_PBE_CORRELATION,
    XC_RPBE_EXCHANGE,
    XC_PW92_CORRELATION,
    XC_PARAMS_LEN // This has to be last!
  };

// Set the functional as defined by the parameters in the
// params array. This function also has a fortran version.
void xc_set_functional(const double params[XC_PARAMS_LEN]);

// Evaluate the exchange-correlation functional defined by
// a previous call to xc_set_functional() with the density
// parameters in densvar, with derivatives up to the given
// order. The result is put in the result array, in graded
// lexicographical ordering. This function also has a fortran 
// version.
void xc_eval_lda(double result[], int order, const double densvars[2]);
void xc_eval_gga(double result[], int order, const double densvars[5]);

// Return a string describing the functional and its implementation.
const char *xc_fun_reference(enum xc_funparams functions);

#endif
