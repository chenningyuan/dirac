module mgenrl
!
!***********************************************************************
!   ------ this module contains general informations -----
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
   character (len=80), save :: check,run
   character (len=8), save :: pgeo
   integer, save :: ir=4,iw=6,ip=7
   integer, save :: labsiz
   integer, save :: oadim,omdim,actdim
   double precision :: tstcnf,tprint
!  double precision, save :: timlim,ti,tx,tim,wti,wtx,wtim
end module mgenrl
