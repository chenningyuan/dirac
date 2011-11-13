      LOGICAL ADDSTO

!     Worst case scenario for HUCEXC; a user running STO-3G (Hueckel basis
!     half the size of the total basis set size
!
      REAL*8 HUCCNT, HUCEXC(MXSHEL)
      INTEGER NHUCAO(8), NHUCBA, IHUCPT(MXSHEL)
      COMMON /HUCKEL/ HUCCNT, HUCEXC, ADDSTO, NHUCAO, NHUCBA, IHUCPT

