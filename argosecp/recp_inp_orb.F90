MODULE RECP_INP_ORBITAL
CONTAINS

SUBROUTINE RECP_INP_INTCHECK(nu,ns,ncrs,mcrs,lcr,lls,inam)
! determine the number of types of one-electron integrals
!   inam=4 : arep
!   inam=5 : sorep
  IMPLICIT NONE
  INTEGER :: I,nu,ns,ncrs,mcrs(*),lcr(*),lls(*),inam(5)

  nu = 3 !initialize nu for s(*),t(*),v(*)

  if ( ncrs .gt. 0 ) then
    do I = 1,ns
      if ( mcrs(I) .ne. 0 ) then
        if ( lcr(mcrs(I)) .gt. 0 ) then
          nu = 4
          inam(4) = 4
          go to 200
        endif
      endif
    enddo
    200 continue
    do I = 1, ns
      if ( mcrs(I) .ne. 0 ) then
        if ( lls(mcrs(I)) .gt. 0 ) then
          nu = nu + 1
          inam(nu) = 5
          go to 208
        endif
      endif
    enddo
    208 continue
  endif
END SUBROUTINE RECP_INP_INTCHECK

SUBROUTINE RECP_INP_ORBITAL_CCNORM(mrcru,mconu,ncons,nrcr,ncon,lmnp1,eta,zet)
! normalize the contraction coefficients
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
          prtint = prtint + 2.0d0 * eta(J,K,I)*eta(J,L,I)       &
                            * sqrt(temp**lmnp1(I) * sqrt(temp))
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

! DO I = 1, ncons
!    print *,'From RECP_INP(CCNORM)'
!    DO J = 1, ncon(I)
!       print *, zet(J,I),(eta(K,J,I),K=1,nrcr(I))
!    ENDDO
! ENDDO
END SUBROUTINE RECP_INP_ORBITAL_CCNORM



SUBROUTINE RECP_INP_ORBITAL_OFFSET(nsopr,nso,nblpr,nst)
! compute contracted orbital and lower-triangle-packed matrix symmetry offsets.
  USE RECP_FUNCTION
  IMPLICIT NONE
! variable : global
  INTEGER  :: nsopr(*),nso(*),nblpr(*),nst
! variable : local
  INTEGER  :: i 

  nsopr(1) = 0
  nblpr(1) = 0
  do i = 2,nst
     nsopr(i) = nsopr(i-1) + nso(i-1)
     nblpr(i) = nblpr(i-1) + ( (nso(i-1)+1)*((nso(i-1)+1)-1)/2 )
  enddo
END SUBROUTINE RECP_INP_ORBITAL_OFFSET


SUBROUTINE RECP_INP_ORBANGSYM(icu,itl,itu,ics,igcs,lai,ixyzir)
! # determine angular momentum symmetries.
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER :: icu,itl,itu,ics,igcs,lai,ixyzir(3)
! local variables
  REAL*8  :: c_sum 
  INTEGER :: I1,I2,I3,ict,ic,it
  LOGICAL :: L1,L2,L3
  DO I1 = 1, 3
     I2 = mod(I1,  3) + 1
     I3 = mod(I1+1,3) + 1
     c_sum = 0.0d0
     ict = 0
     do ic = 1, icu
        do it = itl,itu
           ict = ict + 1
           L1 = mod(IPT_LMNV(I1,it),2).eq.1
           L2 = mod(IPT_LMNV(I2,it),2).eq.0
           L3 = mod(IPT_LMNV(I3,it),2).eq.0
           IF(L1.AND.L2.AND.L3) c_sum=c_sum+IPT_AO2SO(ict,ics,igcs) 
        enddo
     enddo

     IF ( (c_sum.ne.0.0d0).and.(lai.ne.ixyzir(I1)) ) THEN 
        IF ( ixyzir(I1) .ne. 0 ) then
           WRITE (RECP_OUT,'(3X,A)') &
                 'axis choice is unsuitable for spin-orbit integrals'
           call bummer('RECP_INP_ORBANGSYM: ixyzir=',ixyzir(I1),2)
        ENDIF
        ixyzir(I1) = lai
     ENDIF
  ENDDO
ENd SUBROUTINE RECP_INP_ORBANGSYM

END MODULE RECP_INP_ORBITAL 
