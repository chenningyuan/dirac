#ifndef XC_FUN_INTERFACE_CPP
#define XC_FUN_INTERFACE_CPP
#endif
#include <xc_fun_interface.h>

#ifndef XC_FUN
#define XC_FUN

/* 
   Main entry point to DFT exchange-correlation functionals 
   by ulfek 2009 
*/

// If XC_PRECONDITION is defined then the seed polynomials
// will be scaled by the density variables, and then scaled
// back at the end of the calculation. This feature is experimental!
// #define XC_PRECONDITION

int xc_parse_functional(const char *functional_name);

// Set the functional as defined by the parameters in the
// params array. 
void xc_set_functional(enum xc_mode mode, 
		       const double params[XC_PARAMS_LEN]);

//Return the total number of terms calculated by xc_eval
int xc_len(enum xc_mode mode, int order); 

// Evaluate the exchange-correlation functional defined by
// a previous call to xc_set_functional() with the density
// parameters in densvar, with derivatives up to the given
// order. The result is put in the result array, in graded
// lexicographical ordering. This function also has a fortran 
// version.
void xc_eval(double result[], int order, const double densvars[]);
#ifdef WITH_SINGLE
void xc_eval_single(float *result, int order, const float densvars[]);
#endif
//Return the index (into result) for a given derivative, i.e.
//{0,2,1,0} gives the d^3/dS^2dZ term index.
int xc_index(enum xc_mode mode, const int derivative[]);

// Return a string describing the functional and its implementation.
const char *xc_fun_reference(enum xc_funparams functions);

#ifdef WITH_QD
#include <qd/qd_real.h>
void xc_eval_qd(qd_real *result, int order, const qd_real densvars[]);
#endif

#endif
