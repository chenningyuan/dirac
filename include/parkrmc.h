!
! VERSION : $Revision: 6936 $
! DATE    : $Date: 2008-05-19 20:27:23 +0200 (Mo, 19 Mai 2008) $
! FILE    : parkrmc.h
!
!                     PARALLEL INFORMATION
!
!             common blocks file for parallel KR-MCSCF
!     ======================================================
!
!
      PARAMETER (IKRMCPARF = 8)
!
!     integer block         
!     -------------
      COMMON/KRMCPARICB/NUMCPU,LMCPRI
!
!     double precision block
!     ----------------------
!     DOUBLE PRECISION  
!     COMMON/KRMCPARRCB/ 
!
!     logical block
!     -------------
      LOGICAL PARMCSCF
      COMMON/KRMCPARLCB/PARMCSCF
!
!     character block
!     ---------------
!     CHARACTER 
!     COMMON/KRMCPARCCB/ 
