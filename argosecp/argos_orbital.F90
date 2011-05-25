MODULE RECP_ORBITAL

CONTAINS

SUBROUTINE RECP_ORBITAL_CONT_COEFF_NORM(mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)
! normalize the contraction coefficients
! Y. C. Park

  IMPLICIT NONE

! variables : local 
  INTEGER :: mrcru,mconu ! from common (parmi)
  INTEGER :: ncons    ! I  Number of Contractions??
  INTEGER :: nrcr(*)  ! J
  INTEGER :: ncon(*)  ! K
  INTEGER :: lmnp1(*) !
  REAL(8) :: eta(mrcru,mconu,*)
  REAL(8) :: zet(mconu,*)
! variables : local 
  INTEGER :: I, J, K, L
  REAL(8) :: prtint, temp

  DO I = 1, ncons
    DO J = 1, nrcr(I)

!     make eta=1 if ncon=1
      if ( ncon(I) .eq. 1 ) then
        eta(1,1,I) = 1.0d0
        GOTO 312               
      endif

      prtint = eta(J,1,I)**2
      DO K = 2, ncon(I)
        DO L = 1, (K-1)
          temp = (zet(K,I) + zet(L,I)) * 0.5d0
          temp = zet(K,I)/temp * zet(L,I)/temp
          prtint = prtint + 2.0d0 * eta(J,K,I) * eta(J,L,I) * sqrt(temp**lmnp1(I) * sqrt(temp))
        ENDDO
        prtint = prtint + eta(J,K,I)**2
      ENDDO

      prtint = sqrt(prtint)
      DO K = 1, ncon(I)
        eta(J,K,I) = eta(J,K,I) / prtint
      ENDDO

      312 CONTINUE !yc-check
    ENDDO
  ENDDO

END SUBROUTINE RECP_ORBITAL_CONT_COEFF_NORM


!SUBROUTINE ARGOS_ORBITAL_SYM_ORB_NORM
! normalization of symmetry orbitals.
! deleted Y. C. Park
!END SUBROUTINE ARGOS_ORBITAL_SYM_ORB_NORM

END MODULE RECP_ORBITAL 
