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
#include "inc_print.h"
! variables : local 
  INTEGER :: mrcru,mconu ! from common (parmi)
  INTEGER :: ncons    ! I  Number of basis-set blocks
  INTEGER :: nrcr(*)  ! J  Number of contraction coefficient
  INTEGER :: ncon(*)  ! K  Number of primitive
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

! print
! IF (RECP_DBG.GE.5) THEN
!    DO I = 1, ncons
!       print *,'From RECP_INP(CCNORM)'
!       DO J = 1, ncon(I)
!          print *, zet(J,I),(eta(K,J,I),K=1,nrcr(I))
!       ENDDO
!    ENDDO
! ENDIF
END SUBROUTINE RECP_INP_ORBITAL_CCNORM



SUBROUTINE RECP_INP_ORBITAL_OFFSET(nsopr,nso,nblpr,nst,FILE_NUM)
! compute contracted orbital and lower-triangle-packed matrix symmetry offsets.
  USE RECP_FUNCTION
  IMPLICIT NONE
#include "inc_print.h"
! variable : global
  INTEGER  :: nsopr(*),nso(*),nblpr(*),nst,FILE_NUM(*)
! variable : local
  INTEGER  :: i 

  IF(RECP_DBG.GE.1) WRITE(FILE_NUM(1),'(/,A)') ' *** RECP_INP_ORBITAL_OFFSET ***'
  IF(RECP_DBG.GE.1) WRITE(FILE_NUM(1),'(6X,2A)') 'I',':  nsopr  nblpr'
  nsopr(1) = 0
  nblpr(1) = 0
  DO I = 2,nst
     nsopr(i) = nsopr(i-1) + nso(i-1)
     nblpr(i) = nblpr(i-1) + ( (nso(i-1)+1)*((nso(i-1)+1)-1)/2 )
     IF(RECP_DBG.GE.1) WRITE(FILE_NUM(1),'(I7,A,2I7)') I,':',nsopr(i),nblpr(i)
  ENDDO
END SUBROUTINE RECP_INP_ORBITAL_OFFSET


SUBROUTINE RECP_INP_ORBITAL_ANGSYM(icu,itl,itu,ICS,igcs,lai,ixyzir)
! # determine angular momentum symmetries.
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
  INTEGER :: icu,itl,itu,ICS,igcs,lai,ixyzir(3)
! local variables
  REAL*8  :: c_sum 
  INTEGER :: I1,I2,I3,ict,ic,it
  LOGICAL :: L1,L2,L3
  DO I1 = 1, 3
     I2 = mod(I1,  3) + 1
     I3 = mod(I1+1,3) + 1
     c_sum = 0.0d0
     ict = 0
     DO ic = 1, icu
        DO it = itl,itu
           ict = ict + 1
           L1 = mod(IPT_LMNV(I1,it),2).eq.1
           L2 = mod(IPT_LMNV(I2,it),2).eq.0
           L3 = mod(IPT_LMNV(I3,it),2).eq.0
           IF(L1.AND.L2.AND.L3) c_sum=c_sum+IPT_AO2SO(ict,ICS,igcs) 
        ENDDO
     ENDDO

     IF ( (c_sum.NE.0.0d0).AND.(lai.NE.ixyzir(I1)) ) THEN 
        IF ( ixyzir(I1).NE.0 ) THEN
           WRITE (RECP_OUT,'(3X,A)') 'axis choice is unsuitable for spin-orbit integrals'
           CALL bummer('RECP_INP_ORBANGSYM: ixyzir=',ixyzir(I1),2)
        ENDIF
        ixyzir(I1) = lai
     ENDIF
  ENDDO
END SUBROUTINE RECP_INP_ORBITAL_ANGSYM


! --------------- SOCFPD ---------------


SUBROUTINE RECP_SOCFPD_CX0(npr,icx,icxa,mics,mjcs,mau,ict,jct,igcs,jgcs, &
           FILE_NUM,easfb,msfbct)
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
! global variable
  INTEGER npr,icx,icxa,mics(*),mjcs(*),mau(*),ict,jct,igcs,jgcs,FILE_NUM(*)
  LOGICAL easfb,msfbct
! local variable
  INTEGER I,J,ICS,JCS
  REAL(8) CX1

  
  DO I = 1, npr

!    set icx/IPT_CX0 (icxa/IPT_CX0)
!    ------------------------------
     icx = icx + 1
     IPT_CX0(icx) = 0.0d0
     if (easfb) then
       icxa = icxa + 1
       IPT_CX0(icxa) = 0.0d0
     endif

!    set ICS/JCS
!    -----------
     ICS = mics(I)
     JCS = mjcs(I)

!    write file for debug
!    --------------------
     IF (RECP_DBG.GE.1) THEN
        WRITE (FILE_NUM(2),'(A,4I10)') 'I/icx/ICS/JCS',I,icx,ICS,JCS
        IF (easfb) THEN
           WRITE (FILE_NUM(2),'(A,4I10)') 'I/icxa/ICS/JCS',I,icxa,ICS,JCS
        ENDIF
     ENDIF

!    get IPT_CX1
!    -----------
     DO J = 1, mau(I)
        ICS = ICS + 1
        JCS = JCS + 1

        CX1 = (IPT_AO2SO(ict,ICS,igcs) * IPT_AO2SO(jct,JCS,jgcs))
        IPT_CX0(icx) = IPT_CX0(icx) + CX1 
        IF (easfb) IPT_CX0(icxa) = IPT_CX0(icxa) + CX1

        IF (.NOT.msfbct) THEN
           CX1 = (IPT_AO2SO(jct,ICS,igcs) * IPT_AO2SO(ict,JCS,jgcs))
           IPT_CX0(icx) = IPT_CX0(icx) + CX1 
           IF (easfb) IPT_CX0(icxa) = IPT_CX0(icxa) - CX1
        ENDIF

!       write file for debug
        IF (RECP_DBG.GE.1) THEN
           WRITE (FILE_NUM(2),'(A,I7,A,F15.7)')'IPT_CX0(',icx,') = ',IPT_CX0(icx)
           IF (easfb) THEN
              WRITE (FILE_NUM(3),'(A,I7,A,F15.7)')'IPT_CX0(',icxa,') = ',IPT_CX0(icxa)
           ENDIF
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE RECP_SOCFPD_CX0

!SUBROUTINE RECP_INP_ORBITAL_SHIFTCX(mblu,mcxu,mcxu1,mcxu2,icx,icxu,ijsf,icxast,ea)
SUBROUTINE RECP_SOCFPD_SHIFTCX(mblu,mcxu,mcxu1,mcxu2,icx,icxu,ijsf,icxast,ea)
! # move antisymmetric coefficient products down.

! #       offsets: 0               mcxu1  mblu
! # before: a(*)-->................c(*)
! #                0          mcxu mcxu1  mblu
! # now:    a(*)-->cx(*)......x(*),c(*)
! #                0    mcxu2 mcxu       (mblu+mcxu)
! # after:  a(*)-->cx(*),x(*),........
  USE RECP_IPT
  IMPLICIT NONE
#include "inc_print.h"
! global variable
  INTEGER mblu,mcxu,mcxu1,mcxu2,icx,icxu,ijsf,icxast(*)
  LOGICAL ea
! local variable
  INTEGER I,icxl
  if (ea) then
!    # on return, mcxu2 is the offset for the s.o. coefficients.
     mcxu2 = icx
     icxl  = icx + 1
     icxu  = mcxu - icx
     mcxu  = mcxu1 - icxu
!    print *,'socfpd :',icxl, mcxu
     do icx = icxl, mcxu
       IPT_CX0(icx) = IPT_CX0(icx+icxu)
!      print *,'socfpd:cx(',icx,')',IPT_CX0(icx)
     enddo

!    # make pointers relative to mcxu2 instead of 0.
     do I = 1, ijsf
       icxast(I) = icxast(I) - icxu
     enddo

!    # on return, mcxu and mblu divide the workspace as:
!    # a(*)-->[cx//x](1:mcxu),w(1:mblu)
     mblu = mblu - mcxu
     IF (RECP_DBG.GE.2) WRITE(RECP_OUT,'(A,I9)') ' * socfpd: mcxu2=',mcxu2

  else
!    # on return, mcxu and mblu divide the workspace as:
!    # a(*)-->cx(1:mcxu),w(1:mblu)
     mcxu = icx
     mblu = mblu - mcxu
  endif

  IF (RECP_DBG.GE.2) THEN
     write(RECP_OUT,'(A,I9)') ' * socfpd: mcxu=',mcxu
     write(RECP_OUT,'(A,I9)') ' * socfpd: left=',mblu
  ENDIF
END SUBROUTINE RECP_SOCFPD_SHIFTCX


END MODULE RECP_INP_ORBITAL 
