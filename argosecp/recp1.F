!*******************************************************************
!
!     ECPINT : RECP integral routine
!
!     This program is based on argos integral routine
!     with the permission of R. M. Pitzer
!     and maintained by Young Choon Park (KAIST)
!
!*******************************************************************
      SUBROUTINE RECP

      USE RECP_OUTPUT
#include "implicit.h"
#include "inc_print.h"

!     Enter ECPINT
      CALL QENTER('ECPINT')
      CALL TITLER('RECP integral start','*',126)
      CALL OUTPUT_TITLE   

!     Basis setup
      PRINT_LEVEL = 1 !10 
      RECP_DBG    = 0
      RECP_OUT    = 6 ! DIRAC output

!     RECP_INPUT generation
      CALL ARGOSIN

!     Main part
      CALL ARGOS_DRIVER
 
!     Exit ECPINT
      CALL TITLER('RECP integral end','*',126)
      CALL QEXIT('ECPINT')
      END

! ------------------------------------------------

      SUBROUTINE ARGOS_DRIVER
      USE RECP_IPT
      implicit logical(a-z)

#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h"
#include "inc_print.h"

      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst

      INTEGER     LCORED
      PARAMETER ( LCORED = 4500000 )
      INTEGER LCORE

      integer
     & nf(msup),      nc(msup),      ncon(mconsp),    nrcr(mconsp),
     & lmnp1(mconsp), ntl(msfup),    ntu(msfup),      nt(msfup),
     & mcons(msfup),  iprst(mnsfup), npair(2,mnsfup), icxst(2,mnsfup),
     & nfct(mccup),   ngw(mccup),    icb(4,24,mccup),
     & ica(mcup,msup,mgup),          nprir(2,mstup,mnsfup)

!     # ipt(*) is used for workspace stack pointers.
      integer ipt(55)

      real*8
     & zet(mconup,mconsp), eta(mrcrup,mconup,mconsp),
     & COOR_X(mcup,msup), COOR_Y(mcup,msup), COOR_Z(mcup,msup)

      LCORE = LCORED
      msu   = msup
      mstu  = mstup
      mrcru = mrcrup
      mconu = mconup
      mcu   = mcup
      mccu  = mccup

      ALLOCATE (IPT_A(10))
      ipt(1) = 1    ! allocate space for lmnv

      call recp_one( IPT_A, eta, ica,    icb,    icxst,
     &   iprst,  ipt,    lmnp1, mcons,  nc, ncon,
     &   nf,     nfct,   ngw,    npair,     nprir,  nrcr,   nt,
     &   ntl,    ntu,    COOR_X, COOR_Y, COOR_Z,     zet,    lcore )

      DEALLOCATE(IPT_A)
      end


      subroutine recp_one( a, eta, ica, icb, icxst,
     &   iprst,  ipt,    lmnp1,       mcons,  nc,     ncon,
     &   nf,     nfct,   ngw,    npair,     nprir,  nrcr,   nt, 
     &   ntl,    ntu,    x,      y,         z,      zet,    lcore )
!  read the user input and compute the 1-e integrals.
 
      USE RECP_IPT
      USE RECP_OUTPUT
      implicit logical(a-z)

#include "argos_mxvalue.h"
#include "argos_files.h"
#include "argos_error.h" 
#include "inc_print.h" 

      integer INT_AREP
 
      integer
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst 
      common /parmi/
     & mccu, mconu, mcu,   mpru, nbft, nnbft,
     & mrcru, mstu, msu,   ng,   ns,   nst 
 
      integer
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,  
     & lxyzir,         inam,        nnam, mdum
      common /ntgr/
     & mconsu, mru,   mcru, msfu, 
     & msfru,  mnru,   ngcs,  nu,   
     & lxyzir(3),      inam(5),     nnam, mdum(32)
!     # dummy:
      integer lcore, mgcsu
      integer
     & ica(mcu,msu,*), icb(4,24,*), icxst(2,*), iprst(*),
     & ipt(*),         lmnp1(*),    
     & mcons(*),       nc(*),       ncon(*),    nf(*),
     & nfct(*),        ngw(*),      npair(2,*), nprir(2,mstu,*),
     & nrcr(*),        nt(*),       ntl(*),     ntu(*)

      real*8
     & a(*),     eta(mrcru,mconu,*), x(mcu,*), y(mcu,*),
     & z(mcu,*), zet(mconu,*)

!     # local:
      integer i, n
      integer
     & lcr(msup),       lls(msup),      nkcrl(6,msup), nkcru(6,msup),
     & nklsl(4,msup),   nklsu(4,msup),  mcrs(msup),    nir(kaordp),
     & maords(mgcsup),  icxsv1(mgcsup), nd(mstup),     nso(mstup),
     & nsopr(mstup),    nblpr(mstup),   nopir(mstup),  ilxyz(3,mstup),
     & la(mrup,kaordp), ncr(mcrup),     nct(msfup),    mgcs(msfup),
     & lb(msfrup),      ms(msfrup),     mnl(msfrup),   mau(mnrup),
     & mics(mnrup),     mjcs(mnrup),    icxast(mnsfup),
     & icxsv2(mgcsup,mgcsup),           idp(mstup,mstup,mstup),
     & lmn1u,           ncru,           lproju,        ndfac,
     & lmax
      integer mblu,mcxu,mcxu2,CX_TEMP
      real*8 chg(msup), zcr(mcrup), ccr(mcrup), eps
      character*3 ityp(mstup), mtype(msup)
      integer  atebyt, forbyt
      external atebyt, forbyt

      mconsu = mconsp
      mgcsu  = mgcsup
      mru    = mrup
      mcru   = mcrup
      msfu   = msfup
      msfru  = msfrup
      mnru   = mnrup
      mblu   = lcore   !lcore = total available workspace

      do i = 1, 3
        inam(i) = i
      enddo

      ipt(55) = mblu + 1
      mcxu = mblu 

!     # input processing.
      call RECP_INP_MAIN(   a,      ccr,    chg,    eta,
     &   ica,    idp,    ipt,    ityp,   la,     lb,     lcr,
     &   lls,    lmnp1,    lmn1u,  lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,     mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,    nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,  nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,      y,      z,      zcr,
     &   zet,    lcore,  N1MAX,  L1REC,CX_TEMP)

!     # allocate space for il(*)
      ALLOCATE(IPT_IL(nnbft))
!     # mcxu is the remaining space.
      mcxu = ipt(55) - ipt(2)

!     Compute symmetry coefficient products.
      call recp_socfpd( mcxu,mcxu2,mblu,mgcsu,icxast,icxst,
     &   icxsv1, icxsv2, idp,  iprst, la,
     &   lb,     maords, mau,  mcons, mgcs,   mics,
     &   mjcs,   nblpr,  nc,   nd,    nf,     nir,
     &   npair,  nprir,  nrcr, nsopr, nt,     ntl,
     &   ntu, CX_TEMP )
!     Allocate variables
      CALL RECP_IPTA_CORTAB(ipt,ndfac,lmn1u,lproju,ncru,lmax)
      eps = 1.0d-12
      CALL CORTAB(eps,lmax,lmn1u,lproju,ndfac)

!     # compute the one-electron integrals.
      INT_AREP = 0
      IF (NU.eq.4) INT_AREP = 1

      do n = 4, nu
        nnam  = inam(n)
        CALL RECP_ONE_MAIN(
     &        a,         ccr,        chg,       eta,        ica,
     &        icb,       icxast,     icxst,     ilxyz,
     &        iprst,     ipt,        lcr,       lls,        lmnp1,
     &        lproju,    mcons,      mcrs,      nc,         ncon,
     &        ncr,       nf,         nfct,      ngw,        nkcrl,
     &        nkcru,     nklsl,      nklsu,     nopir,      npair,
     &        nprir,     nrcr,       nt,        ntl,        ntu,
     &        x,         y,          z,         zcr,        zet,
     &        nso, INT_AREP, N1MAX,  L1REC,mblu,mcxu,mcxu2 )
      enddo

!     Deallocate variables
      DEALLOCATE (IPT_CX)
      CALL RECP_IPTD_CORTAB
      DEALLOCATE(IPT_IL)
      DEALLOCATE(IPT_LMNV)
      end
