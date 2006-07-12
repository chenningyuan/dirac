!
! File: maxorb.h
!
! hjaaj July 2000:
! MXSHEL = max number of symmetry independent contracted shells
! MXPRIM = max number of symmetry independent primitive shells
!  (e.g. a p shell contains 3 components)
! MXCORB = max total number of contracted functions = max NBASIS
! MXBSETS = max number of basis sets (typically 4 : large comp., small comp. + 
!           density fit sets for both)
!
!     IF you change any of these parameters you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT, MXBSETS
cmi   PARAMETER (MXSHEL = 750, MXPRIM = 2000, MXCORB = 5000,            &
#if ! defined (LARGE_MOLECULES)
      PARAMETER (MXSHEL = 950, MXPRIM = 2000, MXCORB = 5000,            &
#else
      PARAMETER (MXSHEL = 2375,MXPRIM = 5000, MXCORB = 12500,           &
#endif
     &           MXORBT = MXCORB*(MXCORB+1)/2, MXBSETS = 4)
