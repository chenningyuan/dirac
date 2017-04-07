!
!     file: mxcent.h
!
!     MXATOM  = max number of atomic types in the molecule
!     MXCENT  = max number of nuclei + point charges + ghost orbital centers in the molecule
!     MXBSETS = max number of basis sets (p.t. 5 : large comp., small comp. + 
!               density fit sets for both + extended Huckel basis)
!
      INTEGER MXATOM, MXCENT, MXCOOR, MXBSETS
      PARAMETER (MXATOM=20, MXCENT = 300, MXCOOR = 3*MXCENT, MXBSETS=5)
