! -- FILE : dcbnrt.h
!
! FILE : dcbnrt.h
!
! For freezing ("NO ROTATION") orbitals in MCSCF
!
!     Dependencies: primitives.F90
!
      INTEGER NOROT(MXCORB), NOROTC(MXCORB)
      COMMON /DCINRT/ NOROT, NOROTC

      LOGICAL LNOROT
      COMMON /DCLNRT/ LNOROT
! -- end of dcbnrt.h --
