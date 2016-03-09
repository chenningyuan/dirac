module mdavi
!
!***********************************************************************
!   ------ this module contains general informations for Davidson -----
!=======================================================================
!   IR,IW,IP(I4) = logical unit for read, write and punch
!   PGEO(C8) = flag of the ab initio program which created the MO's
!   TIMLIM,TI,TX,TIM,WTI,WTX,WTIM (R8) = timing info
!-----------------------------------------------------------------------
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
   implicit none
!
   integer,save :: itmax
   double precision :: ftol,gtol,stol
end module mdavi
