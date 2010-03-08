#ifdef XC_FUN_INTERFACE_CPP /* C++ part */

#ifndef XC_FUN_INTERFACE /* trap */
#define XC_FUN_INTERFACE

// Note that these parameters have be in sync with those in
// xc_fun_module.f90 (Fortran version)

// Highest order of derivatives. Note that the generated code
// becomes huge if this is too large
#define XC_MAX_ORDER 5

/*
  The following variables are used in the code:
  A = \rho_\alpha
  B = \rho_\beta
  F = \nabla A \cdot \nabla A (sometimes called \sigma_\alpha\alpha)
  G = \nabla B \cdot \nabla B (\sigma_\beta\beta)
  H = \nabla A \cdot \nabla B (\sigma_\alpha\beta)

  Alternatively you can use;
  R = \rho_\alpha + \rho_\beta
  S = \rho_\alpha - \rho_\beta
  Z = \nabla R \cdot \nabla R (sometimes called \sigma)
  X = \nabla S \cdot \nabla S
  Y = \nabla R \cdot \nabla S
*/

enum xc_mode
  {
    XC_A = 0, // Fully spin polarized LDA
    XC_R, // Unpolarized LDA
    XC_AB, // LDA with alpha/beta densities
    XC_RS, // LDA with total density/polarization density
    XC_AF, // Fully spin polarized GGA
    XC_RZ, // Unpolarized GGA
    XC_ABFGH, // Generic GGA alpha,beta variables 
    XC_RSZXY, // Generic GGA R,S variables
  };

enum xc_funparams
  {
    XC_VWN5_CORRELATION = 0,
    XC_SLATER_EXCHANGE,
    XC_LYP_CORRELATION,
    XC_PW91_EXCHANGE,
    XC_PBE_CORRELATION,
    XC_BECKE_EXCHANGE, //Not including the Slater part.
    XC_PBE_EXCHANGE,
    XC_REVPBE_EXCHANGE_UNTESTED,
    XC_PW92_CORRELATION,
    KIN_TF_UNTESTED,       // (orbital-free) kinetic energy functionals next
    KIN_PW91_UNTESTED,
    XC_PARAMS_LEN // This has to be last!
  };

#endif /* trap */
#endif /* ifndef XC_FUN_INTERFACE_CPP */

#ifdef XC_FUN_INTERFACE_F90 /* F90 part */

      integer, parameter ::  XC_MAX_ORDER = 5
      
!     Index of each functional parameter in the argument
!     of set_functional
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
!     Make sure XC_PARAMS_LEN matches the number of parameters above
      integer, parameter ::  XC_PARAMS_LEN = 11
      
!     The different modes (density variables) xc_fun can work with,
!     these are passed as the first argument to set_functional()
       integer, parameter ::   XC_A  = 0 !Fully spin polarized LDA
       integer, parameter ::   XC_R  = 1 !Unpolarized LDA
       integer, parameter ::   XC_AB = 2 !LDA with alpha/beta densities
       integer, parameter ::   XC_RS = 3 !LDA with total density/polarization
       integer, parameter ::   XC_AF = 4 !Fully spin polarized GGA
       integer, parameter ::   XC_RZ = 5 !Unpolarized GGA
       integer, parameter ::   XC_ABFGH = 6 !Generic GGA alpha,beta variables 
       integer, parameter ::   XC_RSZXY = 7 !Generic GGA R,S variables

#endif /* ifndef XC_FUN_INTERFACE_F90 */
