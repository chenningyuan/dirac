!
!     file: mxcent.h
!
!     MXATOM  = max number of atomic types
!     MXCENT  = max number of nuclei + point charges + ghost orbital centers
!     MXBSETS = max number of basis sets (typically 4 : large comp., small comp. + 
!               density fit sets for both)
!
!     IF you change MXCENT you should do a "make depend"
!     and then rebuild the program using the command "make".
!
      INTEGER MXATOM, MXCENT, MXCOOR, MXBSETS
      PARAMETER (MXATOM=20, MXCENT = 300, MXCOOR = 3*MXCENT, MXBSETS=4)
