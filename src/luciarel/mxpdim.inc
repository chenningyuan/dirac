*
      INTEGER IMAXLBLKSZ
#ifdef LUCIAREL_STANDALONE
      PARAMETER (MXPWRD = 8 000 000)
#else
      COMMON /CMXPWRD/ MXPWRD, IMAXLBLKSZ
#endif
!
!     this file contains all required PARAMETERS for 
!     static dimensioning in LUCIAREL
!
!     max. number of roots per irrep
      parameter (MXROOT = 64  )

!     MXPSTT == max. number of symmetry type-types
!     decrease/increase when running smaller/larger systems
      PARAMETER (MXPSTT = 8000 )

!     MXINT_TP == max. number of SO-integral types 
      INTEGER MXINT_TP
      PARAMETER (MXINT_TP = 10000)
!     MXINT_DIA == max. number of diagonal integrals  
!     decrease/increase when running smaller/larger systems
      INTEGER MXINT_DIA
      PARAMETER (MXINT_DIA = 130 000)

      PARAMETER (MXPICI = 10  )
!
!     max. number of GASpaces
      PARAMETER (MXPNGAS = 8)
!     max. no. of double group irreps, limited to D32H: 128
      parameter (MXNDGIRR = 128 )
!     max. no. of point group irreps, limited to D32H: MXNDGIRR/2
      PARAMETER (MXPIRR = MXNDGIRR/2)
!     max. no. of orbital symmetries - to be consistent with the code (as it is): MXNDGIRR
      PARAMETER (MXPOBS = MXNDGIRR  )
!     max. no. of spinors
      PARAMETER (MXPORB = 800)
!     whatever the original meaning of MXPCSM was, it is now used as equivalent to MXNDGIRR, 
!     therefore:  (discovered while implementing linear symmetry)
      PARAMETER (MXPCSM = MXNDGIRR)
!
!     ALLOWD NUMBER OF MK SPACES
      PARAMETER(MXPNMK=100)
!     preliminary for SIGDEN route (this seems like a random guess to me - FIXME?) 
!                S. Knecht March 2011
      parameter (MXNOP    =  32)    ! Max. no. of exc. operators
      parameter (MXPLCCOP =  MXNOP)
