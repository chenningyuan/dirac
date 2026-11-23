!-----------------------------------------
! Module with variables previously defined
! at src/dirrci/iterat.inc include file (COMMON /DAVIDC/, /EIVALS/, /CMTRIX/)
!-----------------------------------------
!  SOURCE  : $Source$
   module iterat_mod
   use param_dirrci_mod, only: N2, N8
   implicit none
   real*8 :: EIGEN1(N8) = 0.0D0   ! initialized in BLOCK DATA DIRDAT before
   real*8 :: EIGEN(N2), EIGEN2(N8), CONV(N8), RESIDU(N8), REFVEC(N8), &
             REFRDU(N8), ECORE, CMTR(N2,N2), CMTI(N2,N2), EVECR(N2,N2), EVECI(N2,N2)
   logical :: DCORR(N8)
   end module iterat_mod
