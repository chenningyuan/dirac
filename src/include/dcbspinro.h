! ********************************************************************
! ****************            NSR module              ****************
! ********************************************************************
! 
!  LINDET - LINDET is used to save that LINEAR symmetry is used in
!           pamprp.F, even if in further steps LINEAR could be
!           changed as occur for UKB, spin-free and ZORA calculations.
!           (see subroutine GLINSM, in src/dirac/dirone.F)
! 
! ********************************************************************
      LOGICAL LINDET
      COMMON LINDET
