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

enum xc_mode
  {
    XC_ORDER_AB = 0, //Using variables in xc_vars_ab
    XC_ORDER_RS      //Using variables in xc_vars_rs
  };

enum xc_fun_type
  {
    XC_LDA = 0, // 2 variables
    XC_GGA      // 5 variables
                // May add unpolarized versions here, with 1 and 2 variables.
  };

// Variables in the total density/spin polarization representation
enum xc_vars_rs
  {
    XC_R = 0, // R = \rho_\alpha + \rho_\beta
    XC_S,     // S = \rho_\alpha - \rho_\beta
    XC_Z,     // Z = \nabla R \cdot \nabla R
    XC_X,     // X = \nabla S \cdot \nabla S
    XC_Y      // Y = \nabla R \cdot \nabla S
  };

// Variables in the alpha/beta representation
enum xc_vars_ab
  {
    XC_A = 0, // \rho_\alpha
    XC_B,     // \rho_\beta
    XC_GAA,   // \nabla rho_a \cdot \nabla rho_a
    XC_GBB,   // \nabla rho_b \cdot \nabla rho_b
    XC_GAB    // \nabla rho_a \cdot \nabla rho_b
  };

enum xc_funparams
  {
    XC_SLATER_EXCHANGE = 0,
    XC_VWN5_CORRELATION,
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
// params array. 
void xc_set_functional(enum xc_fun_type funtype, 
		       const double params[XC_PARAMS_LEN]);

// Use alpha/beta or R/S variables
void xc_set_mode(enum xc_mode mode);

//Return the total number of terms calculated by xc_eval
int xc_len(enum xc_fun_type funtype, int order); 

// Evaluate the exchange-correlation functional defined by
// a previous call to xc_set_functional() with the density
// parameters in densvar, with derivatives up to the given
// order. The result is put in the result array, in graded
// lexicographical ordering. This function also has a fortran 
// version.
// The _rs functions calculate derivatives with respect to R,S,..
// while the _ab versions uses rho_a, rho_b, etc
void xc_eval(double result[], int order, const double densvars[]);

//Return the index (into result) for a given derivative, i.e.
//{0,2,1,0} gives the d^3/dS^2dZ term index.
int xc_index(enum xc_fun_type funtype, const int derivative[]);

// Return a string describing the functional and its implementation.
const char *xc_fun_reference(enum xc_funparams functions);

#ifdef WITH_QD
#include <qd/qd_real.h>
void xc_eval_qd(qd_real *result, int order, const qd_real densvars[]);
#endif

#endif
