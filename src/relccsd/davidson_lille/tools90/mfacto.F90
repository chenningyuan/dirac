module mfacto
!
!***********************************************************************
!   ----- this module contains the factorials -----
!-----------------------------------------------------------------------
!  NOTE: if lfact(i) is FALSE, fact1(i) has not yet been computed
!           it is then computed and stored in tfact(i)
!        if lfact(i) is TRUE, fact(i) is recovered from tfact(i)
!=======================================================================
!     written by: J. P. Flament and V. Vallet
!                 PhLAM
!                 Bat. P5
!                 Universite de Lille-1
!                 F-59655 Villeneuve d'Ascq Cedex
!***********************************************************************
!
  use mparam
  implicit none
!
!     ----- variables in common -----
!
  double precision, save, dimension(0:MAXACT) :: tfact
  logical, save, dimension(0:MAXACT) :: lfact
end module mfacto
