*
      COMMON /CMXPWRD/ MXPWRD
*
* contains all PARAMETERS defining LUCIA
      parameter (MXTSOB =  70)
*  Max. no. of NR orbitals per type and symmetry (now global)
*
      PARAMETER (MXPIRR = 20)
      PARAMETER (MXPOBS = 70)
* Max number of orbitals
      PARAMETER (MXPORB = 200)
* Max number of open shells
      PARAMETER (MXOPEN = 22)
      PARAMETER (MXPICI = 2)
      PARAMETER (MXPSTT = 10000)
*  !! Has to be increased when larger systems are being treated !! (e.g. Tl)
*
      PARAMETER (MXPCSM = 20)
*. Note : MXPNGAS = MXPR4T+6 !!
*. Required in order to handle GAS and RAS within /LUCINP/   
      PARAMETER (MXPNGAS = 6)
      PARAMETER (MXPNSMST = 8)
Contains all parameters defining PICASSO
      parameter (MXNSPIR = 8)
      parameter (MXNDGIRR = 12)
* mx. no. of double group irreps, limited to D2H: 12
* Preliminary for SIGDEN route:
      parameter (MXPLCCOP = 20)  ! preliminary for dimensioning
      parameter (MXNOP = 100)    ! Max. no. of exc. operators
