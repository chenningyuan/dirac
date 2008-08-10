!
! VERSION : $Revision: 6986 $
! DATE    : $Date: 2008-06-03 12:20:27 +0200 (Di, 03 Jun 2008) $
! FILE    : luciprop.h
!
!     This file conatins all information needed for 
!     property calculations using the LUCI* modules in 
!     DIRAC.
!
!     LUCIPROP: activate property calculation with LUCI* 
!               modules
!     NLUCIPRP: total number of one-electron operators
!     LLUCIPRP: pointer list to general common block list 
!               of one-electron operators
!     IPROPPRINT: print level
!
!     ================================================================
!
      PARAMETER(MXLUCIPRP=60)
      PARAMETER(IPROPPRINT=0)
	INTEGER NLUCIPRP, LLUCIPRP
	COMMON/LUCIIPRP/ NLUCIPRP, LLUCIPRP(MXLUCIPRP)
	
      LOGICAL LUCIPROP
	COMMON/LUCILPRP/ LUCIPROP
