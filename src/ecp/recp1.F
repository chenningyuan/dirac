!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

!*******************************************************************
!
!     RECP : RECP integral routine
!
!     This program is based on argos integral routine
!     with the permission of R. M. Pitzer
!     and maintained by Yoon Sup Lee and Young Choon Park (KAIST)
!
!*******************************************************************
      SUBROUTINE RECP

      USE RECP_OUTPUT
      USE RECP_CHK
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

      CALL OUTPUT_TITLE   ! Print title
      CALL RECP_CHK_MAIN  ! Check consistency
      CALL RECP_1E        ! Main part

!     Exit RECP
      CALL TITLER('RECP integral end','*',126)
      CALL QEXIT('RECP')
      END

! ------------------------------------------------

      SUBROUTINE RECP_1E 
!  read the user input and compute the 1-e integrals.
      USE RECP_CSO_TABLE
      USE RECP_INP
      USE RECP_IPT
      USE RECP_NTR
      USE RECP_ONE
      USE RECP_OUTPUT
      USE RECP_SOCFPD
      implicit logical(a-z)
#include "inc_mxvalue.h"
#include "inc_print.h" 
      integer
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      common /parmi/
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      integer
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,inam,nnam,mdum
      common /ntgr/
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,   
     &     lxyzir(3),inam(5),nnam,mdum(32)
      integer, allocatable :: ica(:,:,:)
      integer, allocatable :: icb(:,:,:)
      integer, allocatable :: icxst(:,:)
      integer, allocatable :: iprst(:)
      integer, allocatable :: lmnp1(:)
      integer, allocatable :: mcons(:)
      integer, allocatable :: nc(:)
      integer, allocatable :: ncon(:)
      integer, allocatable :: nf(:)
      integer, allocatable :: nfct(:)
      integer, allocatable :: ngw(:)
      integer, allocatable :: npair(:,:)
      integer, allocatable :: nprir(:,:,:)
      integer, allocatable :: nrcr(:)
      integer, allocatable :: nt(:)
      integer, allocatable :: ntl(:)
      integer, allocatable :: ntu(:)    
      real(kind=8), allocatable :: zet(:,:)
      real(kind=8), allocatable :: eta(:,:,:)
      real(kind=8), allocatable :: x(:,:), y(:,:), z(:,:)
!     # local:
      integer, allocatable :: lcr(:), lls(:)
      integer, allocatable :: nkcrl(:,:), nkcru(:,:)
      integer, allocatable :: nklsl(:,:), nklsu(:,:)
      integer, allocatable :: mcrs(:), nir(:)
      integer, allocatable :: maords(:)
      integer, allocatable :: nd(:), nso(:)
      integer, allocatable :: nsopr(:), nblpr(:), nopir(:)
      integer, allocatable :: ilxyz(:,:)
      integer, allocatable :: la(:,:)
      integer, allocatable :: ncr(:), nct(:), mgcs(:)
      integer, allocatable :: lb(:), ms(:), mnl(:)
      integer, allocatable :: mau(:)
      integer, allocatable :: icxast(:)
      integer, allocatable :: idp(:,:,:)

      integer lmn1u,ncru,lproju,ndfac,lmax
      integer mgcsu,MAXCX,i,n,INT_AREP
      real*8 chg(msup), zcr(mcrup), ccr(mcrup), eps
      character*3 ityp(mstup), mtype(msup)
!
!     Dynamic allocations
!
      allocate( ica(mcup, msup, mgup) )
      allocate( icb(4, 24, mccup) )
      allocate( icxst(2, mnsfup) )
      allocate( iprst(mnsfup) )
      allocate( lmnp1(mconsp) )
      allocate( mcons(msfup) )
      allocate( nc(msup) )
      allocate( ncon(mconsp) )
      allocate( nf(msup) )
      allocate( nfct(mccup) )
      allocate( ngw(mccup) )
      allocate( npair(2, mnsfup) )
      allocate( nprir(2, mstup, mnsfup) )
      allocate( nrcr(mconsp) )
      allocate( nt(msfup) )
      allocate( ntl(msfup) )
      allocate( ntu(msfup) )
      allocate( zet(mconup, mconsp) )
      allocate( eta(mrcrup, mconup, mconsp) )
      allocate( x(mcup, msup) )
      allocate( y(mcup, msup) )
      allocate( z(mcup, msup) )
      allocate( lcr(msup) )
      allocate( lls(msup) )
      allocate( nkcrl(6, msup) )
      allocate( nkcru(6, msup) )
      allocate( nklsl(4, msup) )
      allocate( nklsu(4, msup) )
      allocate( mcrs(msup) )
      allocate( nir(kaordp) )
      allocate( maords(mgcsup) )
      allocate( nd(mstup) )
      allocate( nso(mstup) )
      allocate( nsopr(mstup) )
      allocate( nblpr(mstup) )
      allocate( nopir(mstup) )
      allocate( ilxyz(3, mstup) )
      allocate( la(mrup, kaordp) )
      allocate( ncr(mcrup) )
      allocate( nct(msfup) )
      allocate( mgcs(msfup) )
      allocate( lb(msfrup) )
      allocate( ms(msfrup) )
      allocate( mnl(msfrup) )
      allocate( mau(mnrup) )
      allocate( icxast(mnsfup) )
      allocate( idp(mstup, mstup, mstup) )
      
      CALL RECP_SETPARM
      mgcsu  = mgcsup

!     ================
!     input processing
!     ================
      CALL RECP_INP_MAIN(ccr,chg,eta,ica,idp,ityp,la,lb,lcr,
     &     lls,lmnp1,lmn1u,lproju,maords,mcons,
     &     mcrs,mgcs,mnl,ms,mtype,nblpr,nc,
     &     ncon,ncr,ncru,nct,nd,nf,nir,
     &     nkcrl,nkcru,nklsl,nklsu,nrcr,nso,nsopr,
     &     nt,ntl,ntu,x,y,z,zcr,zet,MAXCX)

      ALLOCATE(IPT_IL(nnbft))  !allocate space for il(*)

!     Compute symmetry coefficient products.
      CALL RECP_SOCFPD_MAIN(icxast,icxst,idp,iprst,la,
     &     lb,maords,mau,mcons,mgcs,nblpr,nc,nd,nf,nir,
     &     npair,nprir,nrcr,nsopr,nt,ntl,ntu,MAXCX)

!     Allocate variables
      eps = 1.0d-12
      CALL RECP_IPTA_CORTAB(ndfac,lmn1u,lproju,ncru,lmax)
      CALL RECP_CSO_TABLE0(eps,lmax,lmn1u,lproju,ndfac)

!     =========================
!     compute the ecp integrals
!     =========================
      INT_AREP = 0
      IF (nu.eq.4) INT_AREP = 1

      ALLOCATE (IPT_A(10))
      DO n = 4, nu
         nnam  = inam(n)
         CALL RECP_ONE_MAIN(IPT_A,ccr,chg,eta,ica,
     &        icb,icxast,icxst,ilxyz,iprst,lcr,lls,lmnp1,
     &        lproju,mcons,mcrs,nc,ncon,ncr,nf,nfct,ngw,nkcrl,
     &        nkcru,nklsl,nklsu,nopir,npair,nprir,nrcr,nt,ntl,
     &        ntu,x,y,z,zcr,zet,nso,INT_AREP)
      ENDDO

!     Deallocate variables
      deallocate( ica )
      deallocate( icb )
      deallocate( icxst )
      deallocate( iprst )
      deallocate( lmnp1 )
      deallocate( mcons )
      deallocate( nc )
      deallocate( ncon )
      deallocate( nf )
      deallocate( nfct )
      deallocate( ngw )
      deallocate( npair )
      deallocate( nprir )
      deallocate( nrcr )
      deallocate( nt )
      deallocate( ntl )
      deallocate( ntu )
      deallocate( zet )
      deallocate( eta )
      deallocate( x )
      deallocate( y )
      deallocate( z )
      deallocate( lcr )
      deallocate( lls )
      deallocate( nkcrl )
      deallocate( nkcru )
      deallocate( nklsl )
      deallocate( nklsu )
      deallocate( mcrs )
      deallocate( nir )
      deallocate( maords )
      deallocate( nd )
      deallocate( nso )
      deallocate( nsopr )
      deallocate( nblpr )
      deallocate( nopir )
      deallocate( ilxyz )
      deallocate( la )
      deallocate( ncr )
      deallocate( nct )
      deallocate( mgcs )
      deallocate( lb )
      deallocate( ms )
      deallocate( mnl )
      deallocate( mau )
      deallocate( icxast )
      deallocate( idp )
      
      DEALLOCATE (IPT_CX)
      CALL RECP_IPTD_CORTAB
      DEALLOCATE(IPT_IL)
      DEALLOCATE(IPT_LMNV)
      DEALLOCATE(IPT_A)
      IF (ALLOCATED(RECP_SET)) DEALLOCATE(RECP_SET)
      END


      SUBROUTINE RECP_SETPARM
      IMPLICIT NONE
#include "inc_mxvalue.h"
#include "inc_print.h" 
      integer
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      common /parmi/
     & mccu,mconu,mcu,mpru,nbft,nnbft,mrcru,mstu,msu,ng,ns,nst
      integer
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,lxyzir,inam,nnam,mdum
      common /ntgr/
     & mconsu,mru,mcru,msfu,msfru,ngcs,nu,
     & lxyzir(3),inam(5),nnam,mdum(32)
      integer i

      DATA RECP_FNUM / 71, 72, 73, 74, 75, 76, 77, 78, 79, 80 /
      DATA RECP_FNAM / 'RECP_XXXXX', 'RECP_XXXXX', 'RECP_XXXXX',
     &                 'RECP_H2___', 'RECP_XXXXX', 'RECP_XXXXX',
     &                 'RECP_XXXXX', 'RECP_XXXXX', 'RECP_XXXXX',
     &                 'RECP_XXXXX' /

!     parmi values
!     ------------
      mccu  = mccup
      mconu = mconup
      mcu   = mcup
      mrcru = mrcrup
      mstu  = mstup
      msu   = msup

!     ntgr values
!     -----------
      mconsu = mconsp
      mru    = mrup
      mcru   = mcrup
      msfu   = msfup
      msfru  = msfrup

      do i = 1, 3
        inam(i) = i
      enddo
      END SUBROUTINE RECP_SETPARM

