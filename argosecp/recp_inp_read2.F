      SUBROUTINE ECPINT_INPUT_READ_RECP
     &      ( ncru, lproju, ncrs, ECP_INPUT,
     &        lcr,lls,nkcrl, nkcru,
     &        nklsl,nklsu, ncr, zcr, ccr )

#include "argos_mxvalue.h"
#include "argos_error.h"

! ### read in expansions for core and spin-orbit potentials.
      INTEGER ncru, lproju, ncrs,ECP_INPUT,
     &        lcr(*),lls(*),nkcrl(6,*), nkcru(6,*),
     &        nklsl(4,*),nklsu(4,*), ncr(*)
      REAL*8  zcr(*), ccr(*)
!     local
      INTEGER I_AREP,l,N_RECP_IN_BLOCK,kcr, kcru,kcrl, icrs,
     &        NUM_AREP, NUM_SOREP

!     # Skip if ncrs=0
      if ( ncrs .eq. 0 ) goto 1000

      ncru = -2
      lproju = 0
      kcru = 0
      do icrs = 1,ncrs

!      # read number of core and spin-orbit
        read (ECP_INPUT,*) NUM_AREP, NUM_SOREP
        lproju = max(lproju, NUM_AREP-1, NUM_SOREP)

        if ( NUM_AREP .gt. 5 ) then
          call bummer('NUM_AREP too large ',NUM_AREP,faterr)
        elseif ( NUM_SOREP .gt. 4 ) then
          call bummer('NUM_SOREP too large ',NUM_SOREP,faterr)
        endif
        lcr(icrs) = NUM_AREP
        lls(icrs) = NUM_SOREP

!      # Read AREP 
        do I_AREP = 1,(NUM_AREP+1)
          read (ECP_INPUT,*) N_RECP_IN_BLOCK
          kcrl = kcru + 1
          kcru = kcru + N_RECP_IN_BLOCK
          nkcrl(I_AREP,icrs) = kcrl
          nkcru(I_AREP,icrs) = kcru

          do kcr = kcrl,kcru
            read (ECP_INPUT,*) ncr(kcr), zcr(kcr), ccr(kcr)
            if ( ncr(kcr) .lt. 0 ) then
              call bummer('ncr(kcr) too small ',ncr(kcr),faterr)
            endif
            ncru = max(ncru,ncr(kcr))
          enddo
        enddo             

!      # Read SOREP 
!      # Be careful with spin-orbit factor
        do l = 1,NUM_SOREP
          read (ECP_INPUT,*) N_RECP_IN_BLOCK
          kcrl = kcru + 1
          kcru = kcru + N_RECP_IN_BLOCK
          nklsl(l,icrs) = kcrl
          nklsu(l,icrs) = kcru

          do kcr = kcrl,kcru
            read (ECP_INPUT,*) ncr(kcr), zcr(kcr), ccr(kcr)
            ccr(kcr)=ccr(kcr)*0.5
            if ( ncr(kcr) .lt. -2 ) then
              call bummer('ncr(kcr) too small ',ncr(kcr),faterr)
            endif
            ncru = max(ncru,ncr(kcr))
          enddo
        enddo

      enddo

      if ( kcru .gt. mcrup ) then
        call bummer('change mcrup (one place) to ',kcru,faterr)
      endif

 1000 CONTINUE
      RETURN
      END


      SUBROUTINE ARGOS_INPUT_CALC_NUCREP
     &         ( ic, nc, chg, NLIST, PRINT_LEVEL, x, y, z )

!     # calculate the nuclear repulsion energy.
#include "argos_mxvalue.h"
#include "argos_error.h"

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

      RETURN
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

!--------------------------------------------------------------------
      SUBROUTINE ECPINT_INPUT_READ_BASISSET
     &          (ECP_INPUT, lmnp1, nrcr, ncon, zet, eta, ncons)

#include "argos_mxvalue.h"
#include "argos_error.h"

!     global viariables
      INTEGER ncons     !Number of basis-set blocks
      INTEGER lmnp1(*)  !
      INTEGER ECP_INPUT !
      INTEGER nrcr(*)   !
      INTEGER ncon(*)   !
      REAL*8  zet(mconup,*), eta(mrcrup,mconup,*)
!     local viariables
      INTEGER I_BLOCK, J_EXP, N_EXP, J_COEFF, N_COEFF

      DO I_BLOCK = 1,ncons

!        # Read numbers to read in each block
         read (ECP_INPUT,*) N_EXP, lmnp1(I_BLOCK), N_COEFF
         if ( N_COEFF .eq. 0 ) N_COEFF = 1

!        # Check maximum value
         if ( N_COEFF .gt. mrcrup ) then
            call bummer('Number of coeff. is exceed ',N_COEFF,faterr)
         elseif ( N_EXP .gt. mconup ) then
            call bummer('Number of exp. is exceed ',N_EXP,faterr)
         endif

!        # Save to global variables
         nrcr(I_BLOCK) = N_COEFF
         ncon(I_BLOCK) = N_EXP

!        # Read exponent and coefficient in each block
         do J_EXP = 1,N_EXP
            read (ECP_INPUT,*) 
     &           zet(J_EXP,I_BLOCK), 
     &          (eta(J_COEFF,J_EXP,I_BLOCK), J_COEFF=1,N_COEFF)
         enddo
      ENDDO
      RETURN
      END

!--------------------------------------------------------------------
      SUBROUTINE RECP_INP_READ_IRREP(ECP_INPUT,idp,nd,ityp)
!     # read the number of product triples.

#include "argos_mxvalue.h"
#include "argos_error.h"

      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
!     global viariables
      INTEGER ECP_INPUT,idp(mstu,mstu,*),nd(*)
      character*3 ityp(*) 

!     local viariables
      INTEGER I_DPT, N_DPT, I,J,K, I2,J2,K2, iprd 
      INTEGER ndptsm(0:9), iprdsm(3,7)

!     the number of product triples for each nst.
      data (ndptsm(i),i=0,9)
     &    / -1, -1, 0, -1, 1, -1, -1, -1, 7, -1 /
!     nst=4 and nst=8 product triples.
      data iprdsm / 2,3,4, 2,5,6, 2,7,8, 3,5,7, 3,6,8, 4,5,8, 4,6,7 /
      do I2 = 1, mstu
         do J2 = 1, mstu
            do K2 = 1, mstu
!           do K2 = 1, nst
               idp(I2,J2,K2) = 0
            enddo
         enddo
      enddo
    
!     # read number of irreps, degeneracy, and character labels.
      read (ECP_INPUT,'(i3,12(i3,a3))') nst,(nd(I),ityp(I),I=1,nst)
      if ( nst .gt. mstup ) then
        call bummer('change mstup (two places) to ',nst,faterr)
      endif

!     print *,'mstu',mstu,'nst',nst
!     assign the totally symmetric irrep binary products.
      do I = 1,nst
        idp(I,I,1) = 1
        idp(I,1,I) = 1
        idp(1,I,I) = 1
      enddo  

!     # read the number of product triples.
      READ(ECP_INPUT,*) N_DPT
 
      IF ( N_DPT .GE. 0 ) THEN

!        read in the totally symmetric irrep product
!        triples from the input file.
         DO I_DPT = 1, N_DPT
            read (ECP_INPUT,*) I, J, K
            idp(I,J,K) = 1
            idp(I,K,J) = 1
            idp(J,I,K) = 1
            idp(J,K,I) = 1
            idp(K,I,J) = 1
            idp(K,J,I) = 1
         ENDDO

      else 
        N_DPT = ndptsm( min( max( 0, nst), 9 ) )
        if ( N_DPT .lt. 0 ) then
          call bummer('syminp: unknown N_DPT defaults for nst=',
     &     nst, faterr )
        else
          do iprd = 1, N_DPT
            I = iprdsm(1,iprd)
            J = iprdsm(2,iprd)
            K = iprdsm(3,iprd)
            idp(I,J,K) = 1
            idp(I,K,J) = 1
            idp(J,I,K) = 1
            idp(J,K,I) = 1
            idp(K,I,J) = 1
            idp(K,J,I) = 1
          enddo
        endif
      endif

      END

!--------------------------------------------------------------------
      SUBROUTINE RECP_INP_READ_AO2SO 
     &          (naords,ngcs,maords,L_IRREP,nir,ECP_INPUT,CX_TEMP)
!     # Read AO to SO transformation matrices 

      USE RECP_IPT
#include "argos_mxvalue.h"
#include "argos_error.h"

!     global viariables
      INTEGER naords, ngcs, maords(*)
      INTEGER L_IRREP(mrup,*)  ! la -> L_IRREP 
      INTEGER nir(*),ECP_INPUT
!     local viariables
      INTEGER AORDS1  ! Ith AO reduction set
      INTEGER N_IRREP, I_IRREP, IRREP_MAX
      INTEGER SO1, N_SO, AO1, N_AO, CX_TEMP
      REAL*8  AO2SO_TEMP1, AO2SO_TEMP2

!     # Read number(N_IRREP) and list(L_IRREP) of irrep 
!       in Ith AO reduction set
      IRREP_MAX = 0
      DO AORDS1 = 1, naords
        read (ECP_INPUT,*) N_IRREP, 
     &                 (L_IRREP(I_IRREP,AORDS1), I_IRREP=1,N_IRREP)
        if ( N_IRREP .gt. mrup ) then
          call bummer('change mrup (one place) to ',N_IRREP,faterr)
        endif

        IF ( N_IRREP.GE.IRREP_MAX) IRREP_MAX = N_IRREP
        nir(AORDS1) = N_IRREP
      ENDDO

!     Allocate 
!     print *,'allocate IPT_AO2SO(',IRREP_MAX,IRREP_MAX,naords,')'
      ALLOCATE( IPT_AO2SO(IRREP_MAX,IRREP_MAX,naords) )

      CX_TEMP = IRREP_MAX*IRREP_MAX*naords 
      CX_TEMP = CX_TEMP*CX_TEMP
!     print *,'allocate IPT_CX_TEMP(',CX_TEMP,')'
      ALLOCATE( IPT_CX_TEMP(CX_TEMP) )
      
! NGCS  Number of transformation matrices relating AOs to
!       SOs to be read in.  Usually this will be the same
!       as the number of AO reduction sets.

      do AORDS1 = 1,ngcs
        read (ECP_INPUT,*) N_SO, N_AO, maords(AORDS1)
!       if ( N_SO .gt. mcsup ) then
!         call bummer('change mcsup (one place) to ',N_SO,faterr)
!       elseif ( N_AO .gt. mctup ) then
!         call bummer('change mctup (one place) to ',N_AO,faterr)
!       endif
 
        do SO1 = 1,N_SO
          read (ECP_INPUT,*) (IPT_AO2SO(AO1,SO1,AORDS1), AO1=1,N_AO)

          DO AO1 = 1,N_AO
            AO2SO_TEMP1 = IPT_AO2SO(AO1,SO1,AORDS1)
            AO2SO_TEMP2 = ABS(IPT_AO2SO_TEMP1)
            IF (AO2SO_TEMP2.NE.0.0d0 .AND. AO2SO_TEMP2.NE.1.0d0 )
     &         IPT_AO2SO(AO1,SO1,AORDS1) 
     &              = SIGN(SQRT(AO2SO_TEMP2),AO2SO_TEMP1)
!           PRINT *,'IPT_AO2SO(',AO1,SO1,AORDS1,')',
!    &              IPT_AO2SO(AO1,SO1,AORDS1)
          ENDDO
        enddo
      enddo

      RETURN
      END

!--------------------------------------------------------------------

      SUBROUTINE ECPINT_INPUT_READ_GEO
     &           (mcu, msu, ng, nbft, ECP_INPUT, NUC_N, BASIS_BLOCK_N,
     &            NUC_OTHER_N, ica, mcons, mgcs, maords, lmnp1, lmn1u,
     &            nrcr, nir, nso, nt, ntl, ntu, la, lb, mcrs,
     &            NUC_CHG, X, Y, Z, NUC_NAME, mru, ncrs,nst)

      USE RECP_FUNCTION
      IMPLICIT NONE

#include "argos_mxvalue.h"
#include "argos_error.h"

!     global
      INTEGER     mru, mcu, msu, ng 
      INTEGER     nbft            ! nbft <- isfr
      INTEGER     ECP_INPUT
      INTEGER     NUC_N           ! NUC_N : ns
      INTEGER     BASIS_BLOCK_N(*) ! BASIS_BLOCK_N : nf
      INTEGER     NUC_OTHER_N(*) ! NUC_OTHER_N : nc
      INTEGER     ica(mcu,msu,*)
      INTEGER     mcons(*), mgcs(*)
      INTEGER     maords(*)  ! from ECPINT_INPUT_READ_AO2SO
      INTEGER     lmnp1(*)   ! (l+1)
      INTEGER     lmn1u      ! maximum value of lmnp1
      INTEGER     nrcr(*), nir(*), nso(*)
      INTEGER     nt(*), ntl(*), ntu(*)
      INTEGER     la(mru,*), lb(*)
      INTEGER     ncrs, mcrs(*)     ! ECP
      REAL*8      NUC_CHG(*)      ! NUC_CHG : chg
      REAL*8      X(mcu,*),Y(mcu,*),Z(mcu,*)
      CHARACTER*3 NUC_NAME(*)     ! NUC_NAME : mtype
      INTEGER     nst
!     local
      INTEGER     I, J, K, I2
      INTEGER     NUC_I          ! NUC_I(is)
      INTEGER     L_QUANT        ! Angular quantum number
      INTEGER     BASIS_BLOCK

      I2 = 0
      lmn1u = 0

!     nbft : total number of symmetry orbital basis functions.
      nbft = 0 

!     CHECKME
      DO I = 1, nst
         nso(i) = 0 
      ENDDO

      DO NUC_I = 1,NUC_N           ! NUC_I(is) -> NUC_N(ns)

!       Read nuclear information
        READ (ECP_INPUT,*) NUC_NAME(NUC_I),BASIS_BLOCK_N(NUC_I), 
     &                     NUC_OTHER_N(NUC_I),NUC_CHG(NUC_I)
!       Error if the number of nuclear is more than mcu
        IF ( NUC_OTHER_N(NUC_I) .GT. mcu ) THEN
          call bummer('change mcup (one place) to ',
     &                NUC_OTHER_N(NUC_I),faterr)
        ENDIF
!       Read nuclear geometry
        DO I = 1,NUC_OTHER_N(NUC_I)
          READ (ECP_INPUT,*) X(I,NUC_I), Y(I,NUC_I), Z(I,NUC_I)
        ENDDO
!       Read 
        IF ( NUC_OTHER_N(NUC_I) .GT. 1 ) THEN
          DO I = 2, ng
            READ (ECP_INPUT,*) (ica(J,NUC_I,I), 
     &            J=1,NUC_OTHER_N(NUC_I) )
          ENDDO  
        ENDIF

!       Match basis-set block to l-quantum number
        DO I = 1,BASIS_BLOCK_N(NUC_I)
          I2 = I2+1
!         Read each 'basis block number' and corresponding 'ang. quant. number'
          read (ECP_INPUT,*) BASIS_BLOCK, L_QUANT

          mcons(I2) = BASIS_BLOCK      ! save 
          mgcs(I2)  = L_QUANT

!         lmn1u : maximum value of lmnp1
          lmn1u     = max(lmn1u,lmnp1(BASIS_BLOCK))

!         nt : 
          nt(I2)=(lmnp1(BASIS_BLOCK)+1)*((lmnp1(BASIS_BLOCK)+1)-1)/2
          ntu(I2)   = (nt(I2)*(lmnp1(BASIS_BLOCK)+2))/3
          ntl(I2)   = ntu(I2)-nt(I2)+1
!         nct(I2)   = NUC_OTHER_N(NUC_I)*nt(I2)

          DO J = 1,nrcr(BASIS_BLOCK)
            DO K = 1,nir(maords(L_QUANT))
              nbft = nbft+1
              nso(la(K,maords(L_QUANT))) 
     &           = nso(la(K,maords(L_QUANT)))+1
              lb(nbft) = nso(la(K,maords(L_QUANT)))
!             print *,'la:',la(K,maords(L_QUANT)),
!    &                'nso:',nso(la(K,maords(L_QUANT)))
            ENDDO
          ENDDO
        ENDDO

        if ( ncrs .gt. 0 ) then
          read(ECP_INPUT,*) mcrs(NUC_I)
        endif
      ENDDO

!     check error
      if ( I2 .gt. msfup ) then
        call bummer('change msfup (two places) to ',I2,faterr)
      elseif ( nbft .gt. msfrup ) then
        call bummer('change msfrup (one place) to ',nbft,faterr)
      endif
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
civan-s
c  *********************
c     new ordering of 
c      d exponents
c  *********************
      IPT_LMNV(1,5) = 2
      IPT_LMNV(2,5) = 0
      IPT_LMNV(3,5) = 0
      IPT_LMNV(1,6) = 1
      IPT_LMNV(2,6) = 1
      IPT_LMNV(3,6) = 0
      IPT_LMNV(1,7) = 1
      IPT_LMNV(2,7) = 0
      IPT_LMNV(3,7) = 1
      IPT_LMNV(1,8) = 0
      IPT_LMNV(2,8) = 2
      IPT_LMNV(3,8) = 0
      IPT_LMNV(1,9) = 0
      IPT_LMNV(2,9) = 1
      IPT_LMNV(3,9) = 1
      IPT_LMNV(1,10) = 0
      IPT_LMNV(2,10) = 0
      IPT_LMNV(3,10) = 2
c  *********************
c     new ordering of
c      f exponents  
c  *********************
      IPT_LMNV(1,11) = 3
      IPT_LMNV(2,11) = 0
      IPT_LMNV(3,11) = 0
      IPT_LMNV(1,12) = 2
      IPT_LMNV(2,12) = 1
      IPT_LMNV(3,12) = 0
      IPT_LMNV(1,13) = 2
      IPT_LMNV(2,13) = 0
      IPT_LMNV(3,13) = 1
      IPT_LMNV(1,14) = 1
      IPT_LMNV(2,14) = 2
      IPT_LMNV(3,14) = 0
      IPT_LMNV(1,15) = 1
      IPT_LMNV(2,15) = 1
      IPT_LMNV(3,15) = 1
      IPT_LMNV(1,16) = 1
      IPT_LMNV(2,16) = 0
      IPT_LMNV(3,16) = 2
      IPT_LMNV(1,17) = 0
      IPT_LMNV(2,17) = 3
      IPT_LMNV(3,17) = 0
      IPT_LMNV(1,18) = 0
      IPT_LMNV(2,18) = 2
      IPT_LMNV(3,18) = 1
      IPT_LMNV(1,19) = 0
      IPT_LMNV(2,19) = 1
      IPT_LMNV(3,19) = 2
      IPT_LMNV(1,20) = 0
      IPT_LMNV(2,20) = 0
      IPT_LMNV(3,20) = 3
c  *********************
c     new ordering of
c      g exponents  
c  *********************
      IPT_LMNV(1,21) = 4
      IPT_LMNV(2,21) = 0
      IPT_LMNV(3,21) = 0
      IPT_LMNV(1,22) = 3
      IPT_LMNV(2,22) = 1
      IPT_LMNV(3,22) = 0
      IPT_LMNV(1,23) = 3
      IPT_LMNV(2,23) = 0
      IPT_LMNV(3,23) = 1
      IPT_LMNV(1,24) = 2
      IPT_LMNV(2,24) = 2
      IPT_LMNV(3,24) = 0
      IPT_LMNV(1,25) = 2
      IPT_LMNV(2,25) = 1
      IPT_LMNV(3,25) = 1
      IPT_LMNV(1,26) = 2
      IPT_LMNV(2,26) = 0
      IPT_LMNV(3,26) = 2
      IPT_LMNV(1,27) = 1
      IPT_LMNV(2,27) = 3
      IPT_LMNV(3,27) = 0
      IPT_LMNV(1,28) = 1
      IPT_LMNV(2,28) = 2
      IPT_LMNV(3,28) = 1
      IPT_LMNV(1,29) = 1
      IPT_LMNV(2,29) = 1
      IPT_LMNV(3,29) = 2
      IPT_LMNV(1,30) = 1
      IPT_LMNV(2,30) = 0
      IPT_LMNV(3,30) = 3
      IPT_LMNV(1,31) = 0
      IPT_LMNV(2,31) = 4
      IPT_LMNV(3,31) = 0
      IPT_LMNV(1,32) = 0
      IPT_LMNV(2,32) = 3
      IPT_LMNV(3,32) = 1
      IPT_LMNV(1,33) = 0
      IPT_LMNV(2,33) = 2
      IPT_LMNV(3,33) = 2
      IPT_LMNV(1,34) = 0
      IPT_LMNV(2,34) = 1
      IPT_LMNV(3,34) = 3
      IPT_LMNV(1,35) = 0
      IPT_LMNV(2,35) = 0
      IPT_LMNV(3,35) = 4
civan-e
      return
      end

