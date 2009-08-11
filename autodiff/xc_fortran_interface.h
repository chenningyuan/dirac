
#if 0
// This is a header file for both Fortran and C++ with symbolic
// constants for DFT functionals.
#endif

#define LDA_VWN    0
#define LDA_SLATER 1
#define GGA_BECKE  2
#define GGA_LYP    3
#define GGA_PBEX   4
#define GGA_PBEC   5

// If you change MAX_ORDER to be very high you might have to
// update dispatch.h to create a macro for that order.
#define MAX_LDA_ORDER 8
#define MAX_GGA_ORDER 5
