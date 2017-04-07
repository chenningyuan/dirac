!
! File: maxorb.h
!
! hjaaj July 2000:
! MXSHEL = max number of symmetry independent contracted shells
! MXPRIM = max number of symmetry independent primitive shells
!  (e.g. a p shell contains 3 components)
! MXCORB = max total number of contracted functions = max NBASIS
!
!     IF you change any of these parameters you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT
      PARAMETER (MXSHEL = 20000,MXPRIM = 40000, MXCORB = 40000,          &
     &           MXORBT = MXCORB*(MXCORB+1)/2)
