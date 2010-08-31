*
#ifdef LUCIAREL_STANDALONE
      PARAMETER (MXPWRD = 8 000 000)
#else
      COMMON /CMXPWRD/ MXPWRD
#endif
*
* contains all PARAMETERS defining LUCIA
      parameter (MXTSOB =  40)
*  Max. no. of NR orbitals per type and symmetry (now global)
*
      PARAMETER (MXPIRR = 20)
      PARAMETER (MXPOBS = 30)
C     PARAMETER (MXPR4T = -2 )
      PARAMETER (MXPR4T =  0 )
      PARAMETER (MXPR0T = 1 )
* Max number of orbitals
      PARAMETER (MXPORB = 510)
* Max number of open shells
      PARAMETER (MXOPEN = 22)
      PARAMETER (MXPXOT = 9)
      PARAMETER (MXPXST = 100 )
      PARAMETER (MXPL = 20 )
      PARAMETER (MXPXT = 25)
      PARAMETER (MXPICI = 2)
      parameter (MXROOT = 30)
      PARAMETER (MXPSTT = 8416)
*  !! Has to be increased when larger systems are being treated !! (e.g. Tl)
*
      PARAMETER (MXPCSM = 20)
      PARAMETER (MXPCTP = 30)
      PARAMETER (MXCNSM = 8)
      PARAMETER (MTYP = 30)
*. Note : MXPNGAS = MXPR4T+6 !!
*. Required in order to handle GAS and RAS within /LUCINP/   
      PARAMETER (MXPNGAS = 6)
      PARAMETER (MXPNSMST = 8)
*. Added by Lasse. Question to Timo does this parameter allready excist?
      PARAMETER (MXCMBSPCK = 3)
*. Largest allowed division of space for perturbation operator
      PARAMETER (MXPPTSPC=20)
Contains all parameters defining PICASSO
      parameter (MXNSPIR = 8)
      parameter (MXNDGIRR = 12)
* mx. no. of double group irreps, limited to D2H: 12
      parameter (MXNATOMS = 100)
      parameter (MXNSPIPDMGRP = 250)
*  !! Has to be increased when larger systems are being treated !! (e.g. PbH)
*
      parameter (MXLISTO = 100)
      parameter (MXBLKO = 1000)
*. ALLOWD NUMBER OF MK SPACES
      PARAMETER(MXPNMK=100)
* Preliminary for SIGDEN route:
      parameter (MXPLCCOP = 20)  ! preliminary for dimensioning
      parameter (MXNOP = 100)    ! Max. no. of exc. operators
