!
! VERSION : $Revision: 6986 $
! DATE    : $Date: 2008-06-03 12:20:27 +0200 (Di, 03 Jun 2008) $
! FILE    : luciprop.h
!
!     This file conatins all information needed for 
!     property calculations inside the KR-CI module in 
!     DIRAC.
!
!     DOPROP_KRCI: activate property calculation with LUCI* 
!                  modules
!     NPROP_KRCI:  total number of one-electron operators
!     LPROP_KRCI:  pointer list to general common block list 
!                  of one-electron operators
!     NPROP_ROOTS_KRCI: 
!                  total number of states (roots) in property 
!                  calculation
!     IPROPPRINT:  print level
!
!     ================================================================
!
      PARAMETER(MXPROP_KRCI=60)
      PARAMETER(IPROPPRINT=0)
!
	INTEGER NPROP_KRCI, LPROP_KRCI, NPROP_ROOTS_KRCI
	COMMON/LUCIIPRP/ NPROP_KRCI, LPROP_KRCI(MXPROP_KRCI), NPROP_ROOTS_KRCI
!	
      LOGICAL DOPROP_KRCI, TWOEPROP, DOOSCILLST
	COMMON/LUCILPRP/ DOPROP_KRCI, TWOEPROP, DOOSCILLST
