#ifndef XC_FUN
#define XC_FUN

/* 
   Main entry point to DFT exchange-correlation functionals 
   by ulfek 2009 
*/

#define XC_PARAMS_LEN 5 // Length of the xc_funparams enum
#define XC_MAX_ORDER 5 // Update xc_eval* if this changes

enum xc_funparams
  {
    XC_SLATER_EXCHANGE = 0,
    XC_VWN_CORRELATION,
    XC_BECKE_EXCHANGE,
    XC_LYP_CORRELATION,
    XC_PBE_EXCHANGE,
    XC_PBE_CORRELATION,
  };

// Set the functional as defined by the parameters in the
// params array. 
void xc_set_functional(const double params[XC_PARAMS_LEN]);

// Evaluate the exchange-correlation functional defined by
// a previous call to xc_set_functional() with the density
// parameters in densvar, with derivatives up to the given
// order. The result is put in the result array, in graded
// lexicographical ordering.
void xc_eval_lda(double result[], int order, const double densvars[2]);
void xc_eval_gga(double result[], int order, const double densvars[5]);


/*
  Use the variable ordering R,S,Z,X,Y; which are
  defined as

      R = \rho_\alpha + \rho_\beta
      S = \rho_\alpha - \rho_\beta

      Z = \nabla R \cdot \nabla R
      X = \nabla S \cdot \nabla S
      Y = \nabla R \cdot \nabla S

      and conversely

      na = (R + S)/2
      nb = (R - S)/2

      Dna^2 = (X+Z+2*Y)/4
      Dnb^2 = (X+Z-2*Y)/4
      DnaDnb = (Z - X)/4
*/

#endif
