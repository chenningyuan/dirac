!*******************************************************************
!
!     RECP : RECP integral routine
!
!     This program is based on argos integral routine
!     with the permission of R. M. Pitzer
!     and maintained by Young Choon Park (KAIST)
!
!*******************************************************************
      SUBROUTINE RECP

      USE RECP_OUTPUT
      USE RECP_INP_READ
      IMPLICIT NONE
#include "priunit.h"
#include "argoscom.h"
#include "inc_print.h"

!     Enter RECP
      CALL QENTER('RECP')
      CALL TITLER('RECP integral start','*',126)

!     Parameters from menu (RECP_LNK_RDMENU)
      PRINT_LEVEL = ECPINT_PRINT ! .PRINT
      RECP_DBG    = ECPINT_DEBUG ! .DEBUG
      RECP_OUT    = LUPRI        ! DIRAC output

!     Print title
      CALL OUTPUT_TITLE   

!     Main part
      CALL RECP_ONE

!     Exit RECP
      CALL TITLER('RECP integral end','*',126)
      CALL QEXIT('RECP')
      END

! ------------------------------------------------


      SUBROUTINE RECP_ONE 
!  read the user input and compute the 1-e integrals.
      USE RECP_IPT
      USE RECP_NTR
      USE RECP_OUTPUT
      implicit logical(a-z)
#include "inc_mxvalue.h"
#include "inc_print.h" 
      integer
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      common /parmi/
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      integer
     & mconsu,mru,mcru,msfu,msfru,mnru,ngcs,nu,lxyzir,inam,nnam,mdum
      common /ntgr/
     & mconsu,mru,mcru,msfu,msfru,mnru,ngcs,nu,   
     & lxyzir(3),inam(5),nnam,mdum(32)

!     # dummy:
      integer
     & ica(mcup,msup,mgup),icb(4,24,mccup),icxst(2,mnsfup),
     & iprst(mnsfup),lmnp1(mconsp),mcons(msfup),nc(msup),
     & ncon(mconsp),nf(msup),nfct(mccup),ngw(mccup),npair(2,mnsfup),
     & nprir(2,mstup,mnsfup),nrcr(mconsp),nt(msfup),ntl(msfup),
     & ntu(msfup)
      real*8
     & zet(mconup,mconsp), eta(mrcrup,mconup,mconsp),
     & x(mcup,msup), y(mcup,msup), z(mcup,msup)
!     # local:
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
      integer mgcsu,CX_TEMP,i,n,INT_AREP
      real*8 chg(msup), zcr(mcrup), ccr(mcrup), eps
      character*3 ityp(mstup), mtype(msup)

      CALL RECP_SETPARM
      mgcsu  = mgcsup

      ALLOCATE (IPT_A(10))

!     # input processing.
      CALL RECP_INP_MAIN( IPT_A,      ccr,    chg,    eta,
     &   ica,    idp,    ityp,   la,     lb,     lcr,
     &   lls,    lmnp1,  lmn1u,  lproju, maords, mcons,
     &   mcrs,   mgcs,   mnl,    ms,     mtype,  nblpr,  nc,
     &   ncon,   ncr,    ncru,   nct,    nd,     nf,     nir,
     &   nkcrl,  nkcru,  nklsl,  nklsu,  nrcr,   nso,    nsopr,
     &   nt,     ntl,    ntu,    x,      y,      z,      zcr,
     &   zet,    N1MAX,  L1REC,CX_TEMP)

      ALLOCATE(IPT_IL(nnbft))  !allocate space for il(*)

!     Compute symmetry coefficient products.
      call recp_socfpd( mgcsu,icxast,icxst,
     &   icxsv1, icxsv2, idp,  iprst, la,
     &   lb,     maords, mau,  mcons, mgcs,   mics,
     &   mjcs,   nblpr,  nc,   nd,    nf,     nir,
     &   npair,  nprir,  nrcr, nsopr, nt,     ntl,
     &   ntu, CX_TEMP )
!     Allocate variables
      CALL RECP_IPTA_CORTAB(ndfac,lmn1u,lproju,ncru,lmax)
      eps = 1.0d-12
      CALL CORTAB(eps,lmax,lmn1u,lproju,ndfac)

!     # compute the one-electron integrals.
      INT_AREP = 0
      IF (nu.eq.4) INT_AREP = 1

      do n = 4, nu
        nnam  = inam(n)
        CALL RECP_ONE_MAIN(
     &        IPT_A,     ccr,        chg,       eta,        ica,
     &        icb,       icxast,     icxst,     ilxyz,
     &        iprst,     lcr,       lls,        lmnp1,
     &        lproju,    mcons,      mcrs,      nc,         ncon,
     &        ncr,       nf,         nfct,      ngw,        nkcrl,
     &        nkcru,     nklsl,      nklsu,     nopir,      npair,
     &        nprir,     nrcr,       nt,        ntl,        ntu,
     &        x,         y,          z,         zcr,        zet,
     &        nso, INT_AREP, N1MAX,  L1REC )
      enddo

!     Deallocate variables
      DEALLOCATE (IPT_CX)
      CALL RECP_IPTD_CORTAB
      DEALLOCATE(IPT_IL)
      DEALLOCATE(IPT_LMNV)
      DEALLOCATE(IPT_A)
      IF (ALLOCATED(RECP_SET)) DEALLOCATE(RECP_SET)
      end


      SUBROUTINE RECP_SETPARM
      IMPLICIT NONE
#include "inc_mxvalue.h"
      integer
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      common /parmi/
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      integer
     & mconsu,mru,mcru,msfu,msfru,mnru,ngcs,nu,lxyzir,inam,nnam,mdum
      common /ntgr/
     & mconsu,mru,mcru,msfu,msfru,mnru,ngcs,nu,
     & lxyzir(3),inam(5),nnam,mdum(32)
      integer i

!     parmi values
!     ------------
      mccu  = mccup
      mconu = mconup
      mcu   = mcup
!     mpru,nbft,nnbft
      mrcru = mrcrup
      mstu  = mstup
      msu   = msup
!     ng,ns,nst


!     ntgr values
!     -----------
      mconsu = mconsp
      mru    = mrup
      mcru   = mcrup
      msfu   = msfup
      msfru  = msfrup
      mnru   = mnrup
!     ngcs,nu,lxyzir
      do i = 1, 3
        inam(i) = i
      enddo
!     nnam,mdum
      END SUBROUTINE RECP_SETPARM

