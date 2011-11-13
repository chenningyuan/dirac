!
!     File: aovec.h
!
!     MXCONT = maximum number of contracted functions in an AO input block
!              (in "MOLECULE.INP" input file)
!     MXAOVC = used for MXCONT in a few places
!
!     IF you change MXCONT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXCONT, MXAOVC
      PARAMETER (MXCONT = 25, MXAOVC = MXCONT)
