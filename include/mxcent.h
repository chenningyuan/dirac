!
!     file: mxcent.h
!
!     MXCENT = max number of nuclei + point charges + ghost orbital centers
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXCENT, MXCOOR
#if ! defined (LARGE_MOLECULES)
      PARAMETER (MXCENT = 60, MXCOOR = 3*MXCENT)
#else
      PARAMETER (MXCENT = 600, MXCOOR = 3*MXCENT)
#endif


