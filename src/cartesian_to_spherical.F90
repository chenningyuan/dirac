      module cartesian_to_spherical

!     SPHCOM and two more routines taken from DIRAC, 
!     slightly modernized by replacing work array 
!     by dynamic allocation and removing dependencies
!     Lucas Visscher, october 2024

      private

      interface get_sphcom
           module procedure get_sphcom
      end interface

      public get_sphcom

      contains

      subroutine get_sphcom(lval1,tramat,morder,minteg)
         integer, intent(in) :: lval1
         integer, intent(in), optional :: morder, minteg
         real(8) :: tramat(:,:)
         integer :: nlm, nxyz, morder1, minteg1

         nxyz = size(tramat,1)
         nlm  = size(tramat,2)

         ! sanity checks
         if (nxyz /= nfunctions(lval1,1)) stop 'incorrect dimension in get_sphcom'
         if (nlm  /= nfunctions(lval1,2)) stop 'incorrect dimension in get_sphcom'

         if (present(morder)) then
            morder1 = morder
         else
            morder1 = 0
         endif

         if (present(minteg)) then
            minteg1 = minteg
         else
            minteg1 = 0
         endif

         call sphcom(lval1-1,tramat,nlm,nxyz,morder1,minteg1)

      end subroutine get_sphcom

! Below the original Fortran77 routines by Trygve Helgaker

!  /* Deck sphcom */
      SUBROUTINE SPHCOM(LVAL,TRAMAT,NLM,NXYZ,MORDER,MINTEG)

!     This routine generates coefficients for transforming from
!     Cartesian to spherical components. It is based on notes by
!     P. Wormer, September 90.

!     MORDER eq 0 : M_l order 0, +1, -1, +2, -2, ..., +LVAL, -LVAL
!     MORDER ne 0 : M_l order -LVAL, ..., -1, 0, 1, ..., +LVAL

!     MINTEG eq 0 : old Hermit normalization
!     MINTEG eq 1 : the smallest coefficient is one (abs.val.).
!     MINTEG eq 2 : the spherical components are normalized to unity

!     tuh 11.10.90

      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (D0 = 0.D0, D1 = 1.D0, D2 = 2.D0, D3 = 3.0D0)
      DIMENSION TRAMAT(NXYZ,NLM), COSSIN(0:LVAL,0:LVAL), PL(0:LVAL)

      LVAL1 = LVAL + 1

!     Legendre coefficients
!     ---------------------

      PL = D0
      DO 100 K = 0, LVAL/2
         PL(LVAL-2*K) = (dble((-1)**K)/dble(2**LVAL))*BINOM(LVAL,K)*BINOM(2*(LVAL-K),LVAL)
  100 CONTINUE

!     Cosine and sine coefficients
!     ----------------------------

      COSSIN = D0
      DO 200 M = 0, LVAL
         COSSIN(M,0) = D1
         DO 210 K = 1, M
            COSSIN(M,K) = COSSIN(M-1,K-1)*dble((-1)**(K-1))
            IF (M .GT. K) COSSIN(M,K) = COSSIN(M,K) + COSSIN(M-1,K)
  210    CONTINUE
  200 CONTINUE

!     Transformation coefficients
!     ---------------------------

      TRAMAT = D0
      DO 300 M = 0, LVAL
         CM = SQRT(D2*FACULT(LVAL-M)/FACULT(LVAL+M))
         IF (M .EQ. 0) CM = D1
         IF (MINTEG.EQ.2) CM = CM/SQRT(FACUL2(2*LVAL-1))
         DO 400 K = MOD(LVAL - M,2), LVAL - M, 2
            IF (M .GT. 0) PL(K) = dble((K+1))*PL(K+1)
            CMK = CM*PL(K)
            DO 500 I = 0, (LVAL - K - M)/2
               CMKI = CMK*BINOM((LVAL - K - M)/2,I)
               DO 600 J = 0, I
                  CMKIJ = CMKI*BINOM(I,J)
                  DO 700 N = 0, M
                    IX = LVAL - 2*J - M + N
                    IX = IX*(IX + 1)/2 + LVAL1 - M - 2*I
                    IF (MORDER .EQ. 0) THEN
                       ILM = MAX(1,2*M + MOD(N,2))
                    ELSE
                       IF (MOD(N,2) .EQ. 1) THEN
                          ILM = 1 + LVAL - M
                       ELSE
                          ILM = 1 + LVAL + M
                       END IF
                    END IF
                    TRAMAT(IX,ILM) = TRAMAT(IX,ILM) + CMKIJ*COSSIN(M,N)
  700             CONTINUE
  600          CONTINUE
  500       CONTINUE
  400    CONTINUE
  300 CONTINUE

!     Renormalize if requested with MINTEG

      IF (MINTEG.EQ.1) THEN
         IOFF = LVAL*LVAL
         DO 800 I = 1, 2*LVAL + 1
            TMIN = TRAMAT(IDAMAX(NXYZ,TRAMAT(1,I),1),I)
            DO 810 J = 1, NXYZ
               TJI = ABS(TRAMAT(J,I))
               IF ((TJI.GT.D0).AND.(TJI.LT.TMIN)) TMIN = TJI
  810       CONTINUE
            TMIN = D1 / TMIN
            CALL DSCAL(NXYZ,TMIN,TRAMAT(1,I),1)
  800    CONTINUE
      END IF
      RETURN
      END
!  /* Deck facult */
      FUNCTION FACULT(N)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer            :: LUPRI   = 6       
      PARAMETER (D1=1.D0)
      IF (N .LT. 0) THEN
         WRITE (LUPRI,'(/,A,I10,/A)') ' Argument less than zero in FACULT:',N,' Program cannot continue.'
         ERROR STOP 'Illegal argument in FACULT'
      ELSE
        FACULT = D1
        DO 100 I = 1, N
           FACULT = FACULT*dble(I)
  100   CONTINUE
      END IF
      RETURN
      END
!  /* Deck facul2 */
      FUNCTION FACUL2(N)

!     tuh

      IMPLICIT REAL*8 (A-H,O-Z)
      integer            :: LUPRI   = 6 
      PARAMETER (D1=1.D0)

!     N < 0

      IF (N .LT. 0) THEN
         FACUL2 = dble(N + 2)
         DO I = N + 4, 1, 2
            FACUL2 = FACUL2*dble(I)
         END DO
         IF (FACUL2 .EQ. dble(0)) THEN
            WRITE (LUPRI,'(/,A,I10,/A)') ' Double factorial undefined for ',N,' Program cannot continue.'
            ERROR STOP 'Illegal argument in FACUL2'
         ELSE
            FACUL2 = D1/FACUL2
         END IF

!     N = 0

      ELSE IF (N.EQ.0) THEN
         FACUL2 = D1

!     N > 0

      ELSE
        FACUL2 = dble(N)
        DO I = N - 2, 1, -2
           FACUL2 = FACUL2*dble(I)
        END DO
      END IF
      RETURN
      END
!  /* Deck binom */
      REAL*8 FUNCTION BINOM(I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer            :: LUPRI   = 6       
      PARAMETER (D1=1.D0)
      IF (I .LT. J) THEN
         WRITE (LUPRI,'(/,A,2I5,/A)') ' Second argument larger than first argument in BINOM:',I,J,' Program cannot continue.'
         ERROR STOP 'Illegal arguments in BINOM'
      ELSE
        BINOM = FACULT(I)/(FACULT(I-J)*FACULT(J))
      END IF
      RETURN
      END

      integer function nfunctions(l_value,basis_angular)
      ! Computes number of Cartesian/spherical functions for a given l-value
          integer, intent(in) :: l_value
          integer, intent(in) :: basis_angular
          if (basis_angular == 1) then
!            Cartesian
             nfunctions = l_value * ( l_value + 1 ) / 2
          else if (basis_angular == 2) then
!            spherical
             nfunctions = 2 * (l_value - 1) + 1
          end if
     end function nfunctions


      end module cartesian_to_spherical
