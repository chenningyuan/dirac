*
      INTEGER IMAXLBLKSZ
#ifdef LUCIAREL_STANDALONE
      PARAMETER (MXPWRD = 8 000 000)
#else
      COMMON /CMXPWRD/ MXPWRD, IMAXLBLKSZ
#endif
!
!     this file contains all PARAMETERS defining LUCIAREL
!
!
!     max. no. of NR orbitals per type and symmetry (now global)
      parameter (MXTSOB = 40)
!
      PARAMETER (MXPIRR = 20  )
      PARAMETER (MXPOBS = 30  )
      PARAMETER (MXPR4T =  0  )
      PARAMETER (MXPR0T = 1   )
      PARAMETER (MXPORB = 510 )
      PARAMETER (MXPXOT = 9   )
      PARAMETER (MXPXST = 100 )
      PARAMETER (MXPL   = 20  )
      PARAMETER (MXPXT  = 25  )
      PARAMETER (MXPICI = 10  )
      parameter (MXROOT = 64  )
!     MXPSTT == max. number of symmetry type-types
#if defined (LARGE_MOLECULES)
      PARAMETER (MXPSTT = 80000)
#else
      PARAMETER (MXPSTT = 8000 )
#endif
!     MXINT_TP == max. number of SO-integral types 
      INTEGER*8 MXINT_TP
      PARAMETER (MXINT_TP = 10000)
!     MXINT_DIA == max. number of diagonal integrals  
!     could be decreased when running smaller systems
      INTEGER*8 MXINT_DIA
#if defined (LARGE_MOLECULES)
      PARAMETER (MXINT_DIA = 500 000)
#else
      PARAMETER (MXINT_DIA = 125 000)
#endif
!
      PARAMETER (MXPCSM = 20)
      PARAMETER (MXPCTP = 30)
      PARAMETER (MXCNSM = 8)
      PARAMETER (MTYP = 30)
!     Required in order to handle GAS and RAS within /LUCINP/   
      PARAMETER (MXPNGAS = 8)
      PARAMETER (MXPNSMST = 8)
!     Largest allowed division of space for perturbation operator
      PARAMETER (MXPPTSPC=20)
!     contains all parameters defining PICASSO
      parameter (MXNSPIR = 8)
!     max. no. of double group irreps, limited to D2H: 12
      parameter (MXNDGIRR = 12)
      parameter (MXNATOMS = 100)
      parameter (MXNSPIPDMGRP = 250)
!
!     Has to be increased when larger systems are being treated !! (e.g. PbH)
!
      parameter (MXLISTO = 100)
      parameter (MXBLKO = 1000)
!     ALLOWD NUMBER OF MK SPACES
      PARAMETER(MXPNMK=100)
!     for MOLCAS compatibility
      parameter (MXTIT = 10)
!     preliminary for SIGDEN route:
      parameter (MXPLCCOP = 20)  ! preliminary for dimensioning
      parameter (MXNOP = 100)    ! Max. no. of exc. operators
