!
!     file: mxcent.h
!
!     MXATOM = max number of atomic types
!     MXCENT = max number of nuclei + point charges + ghost orbital centers
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXATOM, MXCENT, MXCOOR
#if ! defined (LARGE_MOLECULES)
      PARAMETER (MXATOM=10, MXCENT = 200, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXATOM=20, MXCENT = 600, MXCOOR = 3*MXCENT)
#endif
