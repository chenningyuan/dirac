! hjaaj July 2000:
! MXSHEL = max number of symmetry independent contracted shells
! MXPRIM = max number of symmetry independent primitive shells
!  (e.g. a p shell contains 3 components)
! MXCORB = max total number of contracted functions = max NBASIS
! MXBSETS = max number of basis sets (typically 4 : large comp., small comp. + 
!           density fit sets for both)
      INTEGER MXSHEL, MXPRIM, MXCORB, MXORBT, MXBSETS
      PARAMETER (MXSHEL = 500, MXPRIM = 2000, MXCORB = 5000, MXORBT =   &
     &           MXCORB*(MXCORB+1)/2, MXBSETS = 4)
