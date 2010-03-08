#ifdef XC_FUN_INTERFACE_CPP /* C++ part */

#ifndef XC_FUN_INTERFACE /* trap */
#define XC_FUN_INTERFACE

// IMPORTANT:
// If you change parameters in this file please make sure
// they are in sync in all interface definitions (C++ and Fortran versions)

// The following variables are used in the code:
// A = \rho_\alpha
// B = \rho_\beta
// F = \nabla A \cdot \nabla A  sometimes called \sigma_{\alpha\alpha}
// G = \nabla B \cdot \nabla B                   \sigma_{\beta \beta}
// H = \nabla A \cdot \nabla B                   \sigma_{\alpha\beta}
//
// Alternatively you can use:
// R = \rho_\alpha + \rho_\beta
// S = \rho_\alpha - \rho_\beta
// Z = \nabla R \cdot \nabla R  sometimes called \sigma
// X = \nabla S \cdot \nabla S
// Y = \nabla R \cdot \nabla S

// Highest order of derivatives - 
// note that the generated code becomes huge if this is too large
#define XC_MAX_ORDER 5

enum xc_mode
  {
    XC_A = 0, // Fully spin polarized LDA
    XC_R,     // Unpolarized LDA
    XC_AB,    // LDA with alpha and beta densities
    XC_RS,    // LDA with total and spin densitise
    XC_AF,    // Fully spin polarized GGA
    XC_RZ,    // Unpolarized GGA
    XC_ABFGH, // Generic GGA (alpha and beta density variables)
    XC_RSZXY, // Generic GGA (total and spin density variables)
  };

enum xc_funparams
  {
    XC_VWN5_CORRELATION = 0,
    XC_SLATER_EXCHANGE,
    XC_LYP_CORRELATION,
    XC_PW91_EXCHANGE,
    XC_PBE_CORRELATION,
    XC_BECKE_EXCHANGE,           // Not including the Slater part
    XC_PBE_EXCHANGE,
    XC_REVPBE_EXCHANGE_UNTESTED,
    XC_PW92_CORRELATION,
    KIN_TF_UNTESTED,             // (orbital-free) kinetic energy functionals next
    KIN_PW91_UNTESTED,
    XC_KT_EXCHANGE,
    XC_PARAMS_LEN                // This has to be last!
  };

// This must match the order of enum xc_funparams! 
// This is defined as an array of char* to avoid relying
// on the C++ runtime (although the C++ string class
// would be better to use).

static const char *functional_names[] =
  {
    "vwn5c",
    "slaterx",
    "lypc",
    "pw91x",
    "pbec",
    "beckex",
    "pbex",
    "revpbex_untested",
    "pw92c_untested",
    "tfk_untested",
    "pw91k_untested",
    "ktx",
  };

#endif /* trap */
#endif /* ifndef XC_FUN_INTERFACE_CPP */

#ifdef XC_FUN_INTERFACE_F90 /* F90 part */

  integer, parameter ::  XC_MAX_ORDER = 5
      
! index of each functional parameter in the argument of set_functional
  integer, parameter :: XC_VWN5_CORRELATION         =  1
  integer, parameter :: XC_SLATER_EXCHANGE          =  2
  integer, parameter :: XC_LYP_CORRELATION          =  3
  integer, parameter :: XC_PW91_EXCHANGE            =  4
  integer, parameter :: XC_PBE_CORRELATION          =  5
  integer, parameter :: XC_BECKE_EXCHANGE           =  6
  integer, parameter :: XC_PBE_EXCHANGE             =  7
  integer, parameter :: XC_REVPBE_EXCHANGE_UNTESTED =  8
  integer, parameter :: XC_PW92_CORRELATION         =  9
  integer, parameter :: KIN_TF_UNTESTED             = 10
  integer, parameter :: KIN_PW91_UNTESTED           = 11
  integer, parameter :: XC_KT_EXCHANGE              = 12

! make sure xc_params_len matches the number of parameters above
  integer, parameter :: XC_PARAMS_LEN = 12
  
! make sure the order matches the order above
! definitions are the same
  integer, parameter :: XC_A     = 0
  integer, parameter :: XC_R     = 1
  integer, parameter :: XC_AB    = 2
  integer, parameter :: XC_RS    = 3
  integer, parameter :: XC_AF    = 4
  integer, parameter :: XC_RZ    = 5
  integer, parameter :: XC_ABFGH = 6
  integer, parameter :: XC_RSZXY = 7

! list of gga functionals
  integer            :: gga_functionals(6) = (/                    &
                                               XC_BECKE_EXCHANGE,  &
                                               XC_LYP_CORRELATION, &
                                               XC_PBE_EXCHANGE,    &
                                               XC_PBE_CORRELATION, &
                                               KIN_PW91_UNTESTED,  &
                                               XC_KT_EXCHANGE      &
                                              /)

#endif /* ifndef XC_FUN_INTERFACE_F90 */
