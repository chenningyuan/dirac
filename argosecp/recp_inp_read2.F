      SUBROUTINE RECP_INP_NUCREP
     &         ( ic, nc, chg, NLIST, PRINT_LEVEL, x, y, z )

!     # calculate the nuclear repulsion energy.
#include "inc_mxvalue.h"

!     Global variables ------------------------------------
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst

      INTEGER ic, nc(*), chg(*)  !ic->I_CENT
      INTEGER NLIST, PRINT_LEVEL 
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
c       # rearrangement theorem eliminates the ic loop in favor
c       # of the overall (ic) factor.
        DO J_SYM = 1, I_SYM
          CHGPRD = temp * chg(J_SYM)
          if ( J_SYM .eq. I_SYM ) then
c           # repulsions within a set are double-counted, so
c           # compensate with a factor of 1/2.
            CHGPRD = 0.5d0 * CHGPRD
c           # don't include nc(I_SYM)=nc(J_SYM) self-repulsions.
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

      IF (PRINT_LEVEL .GE. 2) 
     &   WRITE (NLIST,'(/,a,f14.8)')
     &         ' * nuclear repulsion energy', NUCREP
      END


      SUBROUTINE ARGOS_INPUT_PRINT_RECP
     &         ( is, mcrs, lcr,nkcrl, nkcru,
     &           NLIST, ncr,nklsl, nklsu, lls,
     &           zcr, ccr, mtype, lblsh )

!     Global variables ------------------------------------
      INTEGER is, mcrs(*), lcr(*),nkcrl(6,*), nkcru(6,*),
     &        NLIST, ncr(*),nklsl(4,*), nklsu(4,*), lls(*)
      REAL*8  zcr(*), ccr(*)
      character*3 mtype(*)
      character*1 lblsh(21)
!     Local variables -------------------------------------
      INTEGER icrs,kcrl,kcru,kcr, 
     &        I_AREP, NUM_AREP, I_SOREP, NUM_SOREP

      icrs = mcrs(is)
      if ( icrs .eq. 0 ) goto 1000
      NUM_AREP = lcr(icrs)
      if ( NUM_AREP .lt. 0 ) goto 816
      kcrl = nkcrl(1,icrs)
      kcru = nkcru(1,icrs)

      WRITE (NLIST,'(/,10X,A)')
     &'-----------------------------------------------------------'

!     RECP : r^(-ncr) * exp( zcr * ccr ) ??
      write (nlist,804) mtype(is),lblsh(NUM_AREP+1),
     &      (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  804 format(//31x,a3,' core potential'//35x,a1,' potential'/
     &      19x,'powers',6x,'exponentials',4x,'coefficients'/
     &      (i22,7x,2(1pg16.7)))

!     Print AREP information
      if ( NUM_AREP .eq. 0 ) go to 1000
      do I_AREP = 1,NUM_AREP
         kcrl = nkcrl(I_AREP+1,icrs)
         kcru = nkcru(I_AREP+1,icrs)
         write (nlist,808) lblsh(I_AREP), lblsh(NUM_AREP+1),
     &         (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  808    format(/33x,a1,' - ',a1,' potential'/
     &         19x,'powers',6x,'exponentials',4x,'coefficients'/
     &         (i22,7x,2(1pg16.7)))
      enddo

!     Print SOREP information
  816 continue
      NUM_SOREP = lls(icrs)
      if ( NUM_SOREP .le. 0 ) goto 1000
      write (nlist,820) mtype(is)
  820 format(//28x,a3,' spin-orbit potential')
      do I_SOREP = 1,NUM_SOREP
         kcrl = nklsl(I_SOREP,icrs)
         kcru = nklsu(I_SOREP,icrs)
         write (nlist,824) lblsh(I_SOREP+1),
     &   (ncr(kcr), zcr(kcr), ccr(kcr), kcr=kcrl,kcru)
  824    format(/35x,a1,' potential'/
     &         19x,'powers',6x,'exponentials',4x,'coefficients'/
     &         (i22,7x,2(1pg16.7)))
      enddo

 1000 CONTINUE
      RETURN
      END

!------------------------------------------------

      subroutine mklmnv(lmn1u)
!     # generate cartesian gaussian exponent array.
!     # IPT_LMNV(*,*) = exponents of the cartesian gaussian basis functions.
!     #             s  p  d   f   g   h   i
!     #       lmn = 0  1  2   3   4   5   6
!     #    numxyz = 1, 3, 6, 10, 15  21  28      =((lmn+1)*(lmn+2))/2
      USE RECP_IPT
      implicit logical (a-z)
      integer i, j, k, lmn, lmn1u, ndx
      INTEGER X1,Y1,Z1,ANG0,NUM0

      ndx = 0
      do 40 lmn = 0,lmn1u-1
        do 30 i = lmn,(lmn+2)/3,-1
          do 20 j = min(i,lmn-i),(lmn-i+1)/2,-1
            k = lmn - i -j
            if( i.eq.j .and. j.eq.k ) then
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = i
              IPT_LMNV(2,ndx) = i
              IPT_LMNV(3,ndx) = i
            elseif( i.eq.j .and.j.gt.k ) then
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = i
              IPT_LMNV(2,ndx) = i
              IPT_LMNV(3,ndx) = k
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = i
              IPT_LMNV(2,ndx) = k
              IPT_LMNV(3,ndx) = i
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = k
              IPT_LMNV(2,ndx) = i
              IPT_LMNV(3,ndx) = i
            elseif( i.gt.j .and. j.eq.k ) then
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = i
              IPT_LMNV(2,ndx) = j
              IPT_LMNV(3,ndx) = j
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = j
              IPT_LMNV(2,ndx) = i
              IPT_LMNV(3,ndx) = j
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = j
              IPT_LMNV(2,ndx) = j
              IPT_LMNV(3,ndx) = i
            elseif( i.gt.j .and. j.gt.k ) then
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = i
              IPT_LMNV(2,ndx) = j
              IPT_LMNV(3,ndx) = k
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = i
              IPT_LMNV(2,ndx) = k
              IPT_LMNV(3,ndx) = j
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = j
              IPT_LMNV(2,ndx) = i
              IPT_LMNV(3,ndx) = k
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = k
              IPT_LMNV(2,ndx) = i
              IPT_LMNV(3,ndx) = j
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = j
              IPT_LMNV(2,ndx) = k
              IPT_LMNV(3,ndx) = i
              ndx = ndx + 1
              IPT_LMNV(1,ndx) = k
              IPT_LMNV(2,ndx) = j
              IPT_LMNV(3,ndx) = i
            endif
   20     continue
   30   continue
   40 continue

!     print *,'mklmnv : ndx',ndx
!     do i = 1, 100 
!        print *,i,'IPT_LMNV',IPT_LMNV(1,i),IPT_LMNV(2,i),IPT_LMNV(3,i)
!     enddo

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
!              print *,NUM0,'X1,Y1,Z1 : ',X1,Y1,Z1 
               NUM0 = NUM0 + 1
            ENDDO
         ENDDO
      ENDDO
      END

