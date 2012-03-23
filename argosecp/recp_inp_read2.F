      SUBROUTINE RECP_INP_NUCREP
     &         ( ic,nc,chg,x,y,z )

!     # calculate the nuclear repulsion energy.
#include "inc_mxvalue.h"
#include "inc_print.h"

!     Global variables ------------------------------------
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst

      INTEGER ic, nc(*), chg(*)  !ic->I_CENT
      REAL*8  x(mcu,*), y(mcu,*), z(mcu,*)

!     Local variables -------------------------------------
      INTEGER I_SYM, J_SYM,
     &        J_CENT, N_J_CENT, !(Number of) J'th atomic center
     &        temp 
      REAL*8  NUCREP, ! nuclear repulsion
     &        CHGPRD  ! charge product

      NUCREP = 0.0d0 
      DO I_SYM = 1, ns  !number of symmetry?
        ic = nc(I_SYM)
        temp = (ic) * chg(I_SYM)
!       # rearrangement theorem eliminates the ic loop in favor
!       # of the overall (ic) factor.
        DO J_SYM = 1, I_SYM
          CHGPRD = temp * chg(J_SYM)
          if ( J_SYM .eq. I_SYM ) then
!           # repulsions within a set are double-counted, so
!           # compensate with a factor of 1/2.
            CHGPRD = 0.5d0 * CHGPRD
!           # don't include nc(I_SYM)=nc(J_SYM) self-repulsions.
            N_J_CENT = nc(J_SYM) - 1
          else
            N_J_CENT = nc(J_SYM)
          endif
!         ---------------------------------------------------------
!
!                      |    ic * chg(I_SYM) * chg(J_SYM)          |
!         NUCREP = SUM | ---------------------------------------- |
!                      | { (x-x0)^2 + (y-y0)^2 + (z-z0)^2 }^(1/2) |
!
!         ---------------------------------------------------------
          DO J_CENT = 1, N_J_CENT
            NUCREP = NUCREP 
     &             + CHGPRD / sqrt( (x(ic,I_SYM)-x(J_CENT,J_SYM))**2
     &                            + (y(ic,I_SYM)-y(J_CENT,J_SYM))**2 
     &                            + (z(ic,I_SYM)-z(J_CENT,J_SYM))**2 )
          ENDDO
        ENDDO
      ENDDO

      IF (PRINT_LEVEL.GE.2) WRITE (RECP_OUT,'(/,a,f14.8)')
     &         ' * nuclear repulsion energy', NUCREP
      END

      SUBROUTINE RECP_INP_PRINTCSO
     &  (is,mcrs,lcr,nkcrl,nkcru,ncr,nklsl,nklsu,lls,
     &   zcr,ccr,mtype,lblsh)
#include "inc_print.h"
!     Global variables 
      INTEGER is, mcrs(*), lcr(*),nkcrl(6,*), nkcru(6,*),
     &        ncr(*),nklsl(4,*), nklsu(4,*), lls(*)
      REAL*8  zcr(*), ccr(*)
      character*3 mtype(*)
      character*1 lblsh(21)
!     Local variables 
      INTEGER icrs,kcrl,kcru,kcr,AREP1,AREP2,SOREP1,SOREP2 

      IF (PRINT_LEVEL.LE.1) RETURN

      icrs = mcrs(is)
      if ( icrs .eq. 0 ) goto 1000
      AREP2 = lcr(icrs)
      if ( AREP2 .lt. 0 ) goto 1000
      kcrl = nkcrl(1,icrs)
      kcru = nkcru(1,icrs)

      WRITE (RECP_OUT,'(/,10X,60A)') ('-', I=1,60)

!     RECP : r^(-ncr) * exp( zcr * ccr ) 

!     Print AREP information (highest angular momentum)
!     -------------------------------------------------
      WRITE (RECP_OUT,'(/,31X,A3,A)') mtype(is),' core potential'
      WRITE (RECP_OUT,'(/,35X,A1,A)') lblsh(AREP2+1),' potential'
      WRITE (RECP_OUT,'(19X,A,6X,A,4X,A)') 
     &      'powers','exponentials','coefficients'
      WRITE (RECP_OUT,'(50(I22,4X,2(F16.7),/))') 
     &      (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)

!     Print AREP information
!     ----------------------
      do AREP1 = 1,AREP2
         kcrl = nkcrl(AREP1+1,icrs)
         kcru = nkcru(AREP1+1,icrs)
         WRITE (RECP_OUT,'(33X,A1,A,A1,A)')
     &         lblsh(AREP1),' - ',lblsh(AREP2+1),' potential'
         WRITE (RECP_OUT,'(19X,A,6X,A,4X,A)')
     &         'powers','exponentials','coefficients'
         WRITE (RECP_OUT,'(50(I22,4X,2(F16.7),/))')
     &         (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
      enddo

!     Print SOREP information
!     -----------------------
      SOREP2 = lls(icrs)
      if ( SOREP2 .le. 0 ) goto 1000
      write (RECP_OUT,820) mtype(is)
  820 format(//28x,a3,' spin-orbit potential')
      do SOREP1 = 1,SOREP2
         kcrl = nklsl(SOREP1,icrs)
         kcru = nklsu(SOREP1,icrs)
         WRITE (RECP_OUT,'(35X,A1,A)') lblsh(SOREP1+1),' potential'
         WRITE (RECP_OUT,'(19X,A,6X,A,4X,A)')
     &         'powers','exponentials','coefficients'
         WRITE (RECP_OUT,'(50(I22,4X,2(F16.7),/))')
     &         (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)

      enddo

 1000 CONTINUE
      RETURN
      END

!------------------------------------------------

      SUBROUTINE RECP_INP_CARTEXP
!     # generate cartesian gaussian exponent array.
!     # IPT_LMNV(*,*) = exponents of the cartesian gaussian basis functions.
!     #             s  p  d   f   g   h   i
!     #       lmn = 0  1  2   3   4   5   6
!     #    numxyz = 1, 3, 6, 10, 15  21  28      =((lmn+1)*(lmn+2))/2
      USE RECP_IPT
      USE RECP_FUNCTION
      IMPLICIT NONE
#include "inc_print.h"
      INTEGER X1,Y1,Z1,ANG0,NUM0,LMNV_MAX

!     Set up variable IPT_LMNV
      LMNV_MAX = 84  ! sum of numxyz
      ALLOCATE(IPT_LMNV(3,LMNV_MAX))
      CALL RECP_SETZERO_I2(IPT_LMNV,3,LMNV_MAX)

!     new ordering of exponents
!     =========================
      NUM0 = 2
      DO ANG0 = 1, 6  !From p(1) to i(6) cartesian fn 
!        print *,'new ordering L=',ANG0
         DO X1 = ANG0, 0, -1
            DO Y1 = (ANG0-X1), 0, -1
               Z1 = (ANG0-X1-Y1)
               IPT_LMNV(1,NUM0) = X1 
               IPT_LMNV(2,NUM0) = Y1
               IPT_LMNV(3,NUM0) = Z1
               IF (RECP_DBG.GE.5) THEN
                  WRITE (RECP_OUT,'(3X,A,I3,A,3I3)')  
     &                  'LMNV No.',NUM0,'  X1,Y1,Z1 : ',X1,Y1,Z1
               ENDIF
               NUM0 = NUM0 + 1
            ENDDO
         ENDDO
      ENDDO
      END

