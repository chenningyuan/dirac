module func_PipekMezey
!================================================================================!
!                  This module is used in optimization module                    !
!           it will provide information about Pipek-Mezey functional             !
!================================================================================!
use memory_allocator
implicit none

public  :: initialize_func, get_func, exp_func, get_gradient, get_hessian

private :: ORBLP1, ORBLO1

type :: PipekMezey_flags
  integer :: iprloc, itrloc
  real*8  :: thfull, thdiag, thgrad
  logical :: prjloc, lgfull, lgdiag, lggrad, lgchck
  character*4 :: hesloc
end type PipekMezey_flags

type :: PipekMezey_variales
  integer :: lupri, lucoef
  integer :: nfbas(2,0:2), n2bbasx, n2bbasxq
  integer :: nish(2), norb(2), npsh(2), nesh(2), norbt
  real*8  :: ssmtrc
  logical :: twocompbss, levyle
  integer :: nz, ipqtoq(4,0:7)
  integer :: nucdep
  real*8  :: thrzer
  integer :: nrefs, iprprj
  integer :: maxref
  real*8  :: prothr
  logical :: ownbas
  character*6  :: reffil(15)    ! reffil(maxref)    -->  see dcbprj.h
  character*72 :: vecref(2,15)  ! vecref(2,maxref)  -->  see dcbprj.h
end type PipekMezey_variales

type( PipekMezey_flags     ), private :: loc
type( PipekMezey_variales  ), private :: var

! Some global parameters
real*8, parameter, private :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

contains

  subroutine initialize_func
  !==============================================================================!
  !            Initialize func_PipekMezey module for Dirac program               !
  !==============================================================================!
  implicit none

#include "dcbloc.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "mxcent.h"
#include "nuclei.h"
#include "thrzer.h"
#include "dcbprj.h"

  ! Internal variables
  integer :: i

  !------------------------------------------------------------------------------!
  ! dcbloc.h
  loc%iprloc     = iprloc
  loc%itrloc     = itrloc
  loc%thfull     = thfull
  loc%thdiag     = thdiag
  loc%thgrad     = thgrad
  loc%prjloc     = prjloc
  loc%lgfull     = lgfull
  loc%lgdiag     = lgdiag
  loc%lggrad     = lggrad
  loc%lgchck     = lgchck
  loc%hesloc     = hesloc
  ! priunit.h
  var%lupri      = lupri
  ! dcbgen.h
  var%lucoef     = lucoef
  ! dcbbas.h
  var%nfbas(:,:) = nfbas(:,:)
  var%n2bbasx    = n2bbasx
  var%n2bbasxq   = n2bbasxq
  ! dcborb.h
  var%nish(:)    = nish(:)
  var%norb(:)    = norb(:)
  var%npsh(:)    = npsh(:)
  var%nesh(:)    = nesh(:)
  var%norbt      = norbt
  ! dcbham.h
  var%ssmtrc     = ssmtrc
  var%twocompbss = twocompbss
  var%levyle     = levyle
  ! dgroup.h
  var%ipqtoq(:,:)= ipqtoq(:,:)
  var%nz         = nz
  ! nuclei.h
  var%nucdep     = nucdep
  ! thrzer.h
  var%thrzer     = thrzer
  ! dcbprj.h
  var%nrefs      = nrefs
  var%iprprj     = iprprj
  var%maxref     = min(maxref,15)
  var%prothr     = prothr
  var%ownbas     = ownbas

  if(maxref /= 15)then
    write(lupri,'(a)')' WARNING: The maxref dimension in localization module and Dirac differ'
    write(lupri,'(a)')'          See type external_variales in localization.f90 '
    write(lupri,'(a,I4)')' maxref(Dirac) =',maxref
    write(lupri,'(a,I4)')' maxref(Localization module) = 15'
  endif

  do i=1,var%maxref
    var%reffil(i)   = reffil(i)
    var%vecref(1,i) = vecref(1,i)
    var%vecref(2,i) = vecref(2,i)
  enddo

  end subroutine initialize_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_func(smat,scoef,pmat,nfnuc,dinvmg)
  !==============================================================================!
  !                       Get value of your functional                           !
  !==============================================================================!
  implicit none

  ! External variables
  integer, intent(in )   :: nfnuc(*)

  real*8,  intent(  out) :: dinvmg
  real*8,  intent(in   ) :: scoef(*)
  real*8,  intent(inout) :: smat(*)
  real*8,  intent(  out) :: pmat(var%nish(1),var%nish(1),4,var%nucdep)

  !------------------------------------------------------------------------------!
  if(loc%prjloc)then ! localization with projection analysis
    call gtovlx(smat,var%ssmtrc)
    call orblp1(smat,scoef,pmat,var%nfbas(1,1),var%nfbas(1,2),var%nish(1), &
                var%nucdep,var%maxref,dinvmg)
  else               ! localization with Mulliken analysis
    call orblo1(smat,scoef,pmat,var%nfbas(1,1),var%nfbas(1,2),var%nish(1), &
                var%nucdep,nfnuc,dinvmg)
  endif

  !------------------------------------------------------------------------------!
  end subroutine get_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_hessian(pmat,hmat,full,totdiag,totnond)
  !******************************************************************************!
  !                                                                              !
  !  Construction of full hessian                                                !
  !                                                                              ! 
  !     HMAT         A = d^2 G / dK* dK     B = d^2 G / dK dK          K         !
  !  / A   B  \       / A_uu   A_ub \        / B_uu   B_ub \      / K(unbared) \ !
  !  \ B*  A* /       \ A_bu   A_bb /        \ B_bu   B_bb /      \ K(bared)   / !
  !                                                                              ! 
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: full
  real*8,  intent(out) :: totdiag, totnond

  real*8,  intent(in ), dimension(:,:,:,:) :: pmat
  real*8,  intent(out), pointer :: hmat(:,:)

  ! Internal variables
  integer :: i, j, k, l, m, n, ij, kl
  integer :: nucdep, ni, ndimh, ni12, lupri
  real*8  :: pij1, pij2, pkml, p2kml, p2lmk, buf, buf1, buf2

  integer, allocatable :: ijn(:,:), iju(:), ijl(:)

  real*8,  allocatable, dimension(:,:,:) :: p1, p2, p3, p4
  real*8,  allocatable, dimension(:,:)   :: pd

  !------------------------------------------------------------------------------!
  lupri = var%lupri

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ni12   = ( ni*ni - ni ) / 2
  ndimh  = 4*ni12

  !------------------------------------------------------------------------------!
  ! Calculate real representation of hessian and get the diagonal dominance info !
  !------------------------------------------------------------------------------!
  totdiag = d0
  totnond = d0

  if(.not.full)then  ! label 111

    if(associated(hmat)) call deallocp(hmat,id='hmat in get_hessian')
    call allocp(hmat,ndimh,1,id='hmat in get_hessian')
    hmat(:,:) = d0

    ! unbared-unbared A^R + B^R contribution
    do m = 1,nucdep
      ij = 0
      do j=1,ni
        do i=j+1,ni
          ij = ij + 1
          pij1 = pmat(i,j,1,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d4*pij1*pij1 - buf*buf )
        enddo
      enddo
    enddo

    ! bared-bared A^R + B^R contribution
    do m = 1,nucdep
      ij = ni12
      do j=1,ni
        do i=j+1,ni
          ij = ij + 1
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d2*(pij1*pij1 + pij2*pij2) - buf*buf )
        enddo
      enddo
    enddo

    ! unbared-unbared A^R - B^R contribution
    do m = 1,nucdep
      ij = 2*ni12
      do j=1,ni
        do i=j+1,ni
          ij = ij + 1
          pij2 = pmat(i,j,2,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d4*pij2*pij2 - buf*buf )
        enddo
      enddo
    enddo

    ! bared-bared A^R + B^R contribution
    do m = 1,nucdep
      ij = 3*ni12
      do j=1,ni
        do i=j+1,ni
          ij = ij + 1
          pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
          buf  = pmat(i,i,1,m) - pmat(j,j,1,m)
          hmat(ij,1) = hmat(ij,1) + d2 * ( d2*(pij1*pij1 + pij2*pij2) - buf*buf )
        enddo
      enddo
    enddo

    do i=1,ndimh
      totdiag = totdiag + dabs( hmat(i,1) )
      hmat(i,1) = - hmat(i,1)  ! we are looking for the minimum
    enddo

    if(loc%iprloc.ge.3)then
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R + B^R)_uu'
      write(lupri,'(8d12.4)') (hmat(i,1),i=1,ni12)
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R + B^R)_bb'
      write(lupri,'(8d12.4)') (hmat(i,1),i=ni12+1,2*ni12)
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R - B^R)_uu'
      write(lupri,'(8d12.4)') (hmat(i,1),i=2*ni12+1,3*ni12)
      write(lupri,'(/a,2f10.5)') 'Print diagonal (A^R - B^R)_bb'
      write(lupri,'(8d12.4)') (hmat(i,1),i=3*ni12+1,4*ni12)
    endif


  else ! label 111

    if(associated(hmat)) call deallocp(hmat,id='hmat in get_hessian')
    call allocp(hmat,ndimh,ndimh,id='hmat in get_hessian')
    hmat(:,:) = d0

    call alloc(ijn,ni,ni,      id='ijn in get_hessian')
    call alloc(iju,ni,         id='iju in get_hessian')
    call alloc(ijl,ni,         id='ijl in get_hessian')
    call alloc(pd,ni,nucdep,   id='pd  in get_hessian')
    call alloc(p1,ni,ni,nucdep,id='p1  in get_hessian')
    call alloc(p2,ni,ni,nucdep,id='p2  in get_hessian')
    call alloc(p3,ni,ni,nucdep,id='p3  in get_hessian')
    call alloc(p4,ni,ni,nucdep,id='p4  in get_hessian')

    ! prepare some arrays for effective construction of the hessian
    do m=1,nucdep
      do i=1,ni
        pd(i,m) = pmat(i,i,1,m)
        do j=1,ni
          p1(j,i,m) = pmat(j,i,1,m)
          p2(j,i,m) = pmat(j,i,2,m)
          p3(j,i,m) = pmat(j,i,3,m)
          p4(j,i,m) = pmat(j,i,4,m)
        enddo
      enddo
    enddo

    ijn(:,:) = 0
    kl = 0
    do l=1,ni
      do k=l+1,ni
        kl = kl + 1
        ijn(k,l) = kl
        ijn(l,k) = ijn(k,l)
      enddo
      iju(l) = l - 1
      ijl(l) = l + 1
    enddo

    !----------------------------------------------------------------------------!
    !                       Construction of the HESSIAN                          !
    !----------------------------------------------------------------------------!
    !                            kl = < 1 , ni12 >                               !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = 0
      do l=1,ni
        do k=l+1,ni
          kl = kl + 1
          buf  = d4*p1(k,l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)
 
          ! (A^R + B^R)_uu = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (32) }                 + Real{ (36) }
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) + buf*p1(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) + buf*p1(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) - buf*p1(j,l,m) - (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) - buf*p1(i,k,m) - (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !                       kl = < ni12 + 1 , 2*ni12 >                           !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = ni12
      do l=1,ni
        do k=l+1,ni
          kl = kl + 1
          buf1 = d2*p3(k,l,m)
          buf2 = d2*p4(k,l,m)
          buf  = d2*buf1
          pkml =    pd(k,m) - pd(l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)

          ! (A^R + B^R)_ub = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (33) }                 + Real{ (37) }
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) + buf*p1(j,k,m) + (p2kml-pd(j,m))*p3(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) + buf*p1(i,l,m) - (p2lmk-pd(i,m))*p3(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) - buf*p1(j,l,m) + (p2lmk-pd(j,m))*p3(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) - buf*p1(i,k,m) - (p2kml-pd(i,m))*p3(i,l,m)
          enddo

          ! (A^R + B^R)_bb = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (35) }                 + Real{ (39) }
          n = ni12
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) &
                               - buf1*p3(j,k,m) - buf2*p4(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) &
                               + buf1*p3(i,l,m) + buf2*p4(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) &
                               + buf1*p3(j,l,m) + buf2*p4(j,l,m) + (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) &
                               - buf1*p3(i,k,m) - buf2*p4(i,k,m) + (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !                       kl = < 2*ni12 + 1 , 3*ni12 >                         !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = 2*ni12
      do l=1,ni
        do k=l+1,ni
          kl = kl + 1
          buf  = d4*p2(k,l,m)
          pkml =    pd(k,m) - pd(l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)

          ! (B^I - A^I)_uu = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (36) }                - Imag{ (32) }
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) - buf*p1(j,k,m) - (p2kml-pd(j,m))*p2(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) - buf*p1(i,l,m) + (p2lmk-pd(i,m))*p2(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) + buf*p1(j,l,m) - (p2lmk-pd(j,m))*p2(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) + buf*p1(i,k,m) + (p2kml-pd(i,m))*p2(i,l,m)
          enddo

          ! (B^I - A^I)_bu = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (38) }                - Imag{ (34) }
          n = ni12
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) + buf*p3(j,k,m) - (p2kml-pd(j,m))*p4(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) - buf*p3(i,l,m) - (p2lmk-pd(i,m))*p4(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) - buf*p3(j,l,m) - (p2lmk-pd(j,m))*p4(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) + buf*p3(i,k,m) - (p2kml-pd(i,m))*p4(i,l,m)
          enddo

          ! (A^R - B^R)_uu = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (32) }                 - Real{ (36) }
          n = 2*ni12
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) - buf*p2(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) + buf*p2(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) + buf*p2(j,l,m) + (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) - buf*p2(i,k,m) + (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !                    kl = < 3*ni12 + 1 , 4*ni12 >                            !
    !----------------------------------------------------------------------------!
    do m = 1,nucdep
      kl = 3*ni12
      do l=1,ni
        do k=l+1,ni
          kl = kl + 1
          buf1 = d2*p3(k,l,m)
          buf2 = d2*p4(k,l,m)
          buf  = d2*buf2
          pkml =    pd(k,m) - pd(l,m)
          p2kml= d2*pd(k,m) - pd(l,m)
          p2lmk= d2*pd(l,m) - pd(k,m)

          ! (B^I - A^I)_ub = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (37) }                - Imag{ (33) }
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k),kl) = hmat(ijn(j,k),kl) - buf*p1(j,k,m) - (p2kml-pd(j,m))*p4(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l),kl) = hmat(ijn(i,l),kl) - buf*p1(i,l,m) + (p2lmk-pd(i,m))*p4(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l),kl) = hmat(ijn(j,l),kl) + buf*p1(j,l,m) - (p2lmk-pd(j,m))*p4(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k),kl) = hmat(ijn(i,k),kl) + buf*p1(i,k,m) + (p2kml-pd(i,m))*p4(i,l,m)
          enddo

          ! (B^I - A^I)_bb = Imag{ d^2 G / dK*_ij dK*_kl } - Imag{ d^2 G / dK*_ij dK_kl }
          !                = - Imag{ (39) }                - Imag{ (35) }
          n = ni12
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) &
                               - buf1*p4(j,k,m) + buf2*p3(j,k,m) + (p2kml-pd(j,m))*p2(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) &
                               + buf1*p4(i,l,m) - buf2*p3(i,l,m) + (p2lmk-pd(i,m))*p2(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) &
                               + buf1*p4(j,l,m) - buf2*p3(j,l,m) + (p2lmk-pd(j,m))*p2(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) &
                               - buf1*p4(i,k,m) + buf2*p3(i,k,m) + (p2kml-pd(i,m))*p2(i,l,m)
          enddo

          ! (A^R - B^R)_ub = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (33) }                 - Real{ (37) }
          n = 2*ni12
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) - buf*p2(j,k,m) + (p2kml-pd(j,m))*p3(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) + buf*p2(i,l,m) + (p2lmk-pd(i,m))*p3(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) + buf*p2(j,l,m) + (p2lmk-pd(j,m))*p3(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) - buf*p2(i,k,m) + (p2kml-pd(i,m))*p3(i,l,m)
          enddo

          ! (A^R - B^R)_bb = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
          !                = Real{ (35) }                 - Real{ (39) }
          n = 3*ni12
          do j=1,iju(k)   ! i==k
           hmat(ijn(j,k)+n,kl) = hmat(ijn(j,k)+n,kl) &
                               - buf1*p3(j,k,m) - buf2*p4(j,k,m) + (p2kml-pd(j,m))*p1(j,l,m)
          enddo

          do i=ijl(l),ni  ! j==l
           hmat(ijn(i,l)+n,kl) = hmat(ijn(i,l)+n,kl) &
                               + buf1*p3(i,l,m) + buf2*p4(i,l,m) + (p2lmk-pd(i,m))*p1(i,k,m)
          enddo

          do j=1,iju(l)   ! i==l
           hmat(ijn(j,l)+n,kl) = hmat(ijn(j,l)+n,kl) &
                               + buf1*p3(j,l,m) + buf2*p4(j,l,m) + (p2lmk-pd(j,m))*p1(j,k,m)
          enddo

          do i=ijl(k),ni  ! j==k
           hmat(ijn(i,k)+n,kl) = hmat(ijn(i,k)+n,kl) &
                               - buf1*p3(i,k,m) - buf2*p4(i,k,m) + (p2kml-pd(i,m))*p1(i,l,m)
          enddo

        enddo
      enddo
    enddo

    !----------------------------------------------------------------------------!
    !            Calculate diagonal dominance and deallocate arrays              !
    !----------------------------------------------------------------------------!
    do j=1,ndimh
      do i=1,j
        if(i==j)then
          totdiag = totdiag + dabs( hmat(i,j) )
        else
          totnond = totnond + dabs( hmat(i,j) )
        endif
        hmat(i,j) = - hmat(i,j)   ! we are looking for the minimum
      enddo
    enddo
    totnond = d2*totnond

    !----------------------------------------------------------------------------!
    !                              PRINT section                                 !
    !----------------------------------------------------------------------------!
    if(loc%iprloc.ge.3)then
      write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_uu'
      do i=1,ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=1,ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_uu'
      do i=1,ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=2*ni12+1,3*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_uu'
      do i=2*ni12+1,3*ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=2*ni12+1,3*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_ub'
      do i=1,ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=ni12+1,2*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_ub'
      do i=1,ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*ni12+1,4*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_ub'
      do i=2*ni12+1,3*ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*ni12+1,4*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_bb'
      do i=ni12+1,2*ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=ni12+1,2*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_bb'
      do i=ni12+1,2*ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*ni12+1,4*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_bb'
      do i=3*ni12+1,4*ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=3*ni12+1,4*ni12)
      enddo

      write(lupri,'(/a,2f10.5)') 'Print (B^I - A^I)_bu'
      do i=ni12+1,2*ni12
        write(lupri,'(i4,15f10.5)') i,(hmat(i,j),j=2*ni12+1,3*ni12)
      enddo
    endif

    ! deallocate arrays
    call dealloc(ijn,id='ijn in get_hessian')
    call dealloc(iju,id='iju in get_hessian')
    call dealloc(ijl,id='ijl in get_hessian')
    call dealloc(pd, id='pd  in get_hessian')
    call dealloc(p1, id='p1  in get_hessian')
    call dealloc(p2, id='p2  in get_hessian')
    call dealloc(p3, id='p3  in get_hessian')
    call dealloc(p4, id='p4  in get_hessian')

  endif ! label 111

  ! Calculate mean value
  totdiag = totdiag / dfloat(ndimh)
  totnond = totnond / dfloat(ndimh*ndimh-ndimh)

  end subroutine get_hessian

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_gradient(pmat,gradient,dnrmg)
  !******************************************************************************!
  !                                                                              !
  !                     Construction of the gradient                             !
  !                                                                              ! 
  !******************************************************************************!
  implicit none

  ! External variables
  real*8, intent(out) :: dnrmg

  real*8, intent(in ), dimension(:,:,:,:) ::  pmat
  real*8, intent(out), dimension(:)       ::  gradient

  ! Internal variables
  integer :: i, j, k, ij1, ij2, ij3, ij4
  integer :: nucdep, ni, ndimg, ni12, lupri

  real*8  :: dtmp, dnrm2

  !------------------------------------------------------------------------------!
  lupri = var%lupri
  gradient(:) = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimg  = size(gradient,1)
  ni12   = ndimg / 4

  !------------------------------------------------------------------------------!
  !               Calculate real representation of the gradient                  !
  !------------------------------------------------------------------------------!
  do k = 1,nucdep
    ij1 = 0
    do j=1,ni
      do i=j+1,ni
        ij1 = ij1 + 1
        ij2 = ij1 + ni12
        ij3 = ij2 + ni12
        ij4 = ij3 + ni12
        dtmp = d2*(pmat(j,j,1,k) - pmat(i,i,1,k))
        gradient(ij1) = gradient(ij1) + dtmp * pmat(i,j,1,k)
        gradient(ij2) = gradient(ij2) + dtmp * pmat(i,j,3,k)
        gradient(ij3) = gradient(ij3) - dtmp * pmat(i,j,2,k)
        gradient(ij4) = gradient(ij4) - dtmp * pmat(i,j,4,k)
      enddo
    enddo
  enddo

  do i=1,ndimg
    gradient(i) = - gradient(i)  ! we are looking for the minimum
  enddo

  if(loc%iprloc.ge.2)then
    write(lupri,'(/a,2f10.5)') 'Print gradient^R_u'
    write(lupri,'(8d12.4)') (gradient(i),i=1,ni12)
    write(lupri,'(/a,2f10.5)') 'Print gradient^R_b}'
    write(lupri,'(8d12.4)') (gradient(i),i=1+ni12,2*ni12)
    write(lupri,'(/a,2f10.5)') 'Print gradient^I_u'
    write(lupri,'(8d12.4)') (gradient(i),i=1+2*ni12,3*ni12)
    write(lupri,'(/a,2f10.5)') 'Print gradient^I_b'
    write(lupri,'(8d12.4)') (gradient(i),i=1+3*ni12,4*ni12)
  endif

  ! Norm of the gradient
  dnrmg = dnrm2(ndimg,gradient,1)

  end subroutine get_gradient

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine exp_func(pkappa,tgrad,eval,ni,dinvmq)
  !******************************************************************************!
  !                                                                              !
  !       Calculate value of the quadratic approximation to the functional       !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  integer, intent(in ) :: ni
  real*8,  intent(out) :: dinvmq

  real*8,  intent(in), pointer :: pkappa(:), tgrad(:), eval(:)

  ! Internal variables
  integer :: i

  !------------------------------------------------------------------------------!
  dinvmq = d0
  do i=1,size(eval,1)
    dinvmq = dinvmq + d2*pkappa(i)*tgrad(i) + pkappa(i)*eval(i)*pkappa(i)
  enddo
  dinvmq = dinvmq / dfloat(ni)

  end subroutine exp_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLP1(SMAT,SCOEF,PMAT,NBL,NBS,NI,NUCDEP,MAXREF,DINVM)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
      USE SELECT_COEFF_OWN_BAS
      IMPLICIT NONE

      ! External variables
      INTEGER, INTENT(IN ) :: NBL, NBS, NI, NUCDEP, MAXREF
      REAL*8,  INTENT(OUT) :: DINVM

      REAL*8,  INTENT(IN   ) :: SCOEF(*)
      REAL*8,  INTENT(INOUT) :: SMAT(*)
      REAL*8,  INTENT(  OUT) :: PMAT(NI,NI,4,NUCDEP)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      REAL*8, PARAMETER :: DM1=-1.0D0

      LOGICAL   :: TOBE

      INTEGER :: I, J, K, IZ, IER, IFRP, ICFS, NB, NQ, NCDIM, NBDIM, NMOL
      INTEGER :: IVEC, JOFF, NREF, N2REF, LUPRI, NPSH, NESH, NORB, NREFS
      INTEGER :: NSTR(2,0:2,0:MAXREF), NVECS(2,0:2), IOFF(2), IDUMARR(2), IPQTOQ(4,0:7)

      REAL*8    :: DDUMARR(2), DINV, QLOC2, DLOC

      INTEGER, ALLOCATABLE, DIMENSION(:) :: KVECTMP, JREF, IPIV
      REAL*8,  ALLOCATABLE, DIMENSION(:) :: RCOEF, CMO, BV, GVEC, AMAT, WMAT

      TYPE :: PTR2ARRAY
        INTEGER, DIMENSION(:), POINTER :: P
      END TYPE PTR2ARRAY

      TYPE(PTR2ARRAY), ALLOCATABLE :: KVEC(:,:)

  !------------------------------------------------------------------------------!
      IFRP   = 1
      LUPRI= VAR%LUPRI
      NPSH = VAR%NPSH(IFRP)
      NESH = VAR%NESH(IFRP)
      NORB = VAR%NORB(IFRP)
      NREFS= VAR%NREFS
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)
      VECREF(:,:) = VAR%VECREF(:,:)
      REFFIL(:)   = VAR%REFFIL(:)

      ALLOCATE( KVEC(2,NREFS), STAT=IER )
      IF(IER /= 0)THEN
        WRITE(LUPRI,*) "Error in allocation of KVEC array of pointers in OBLP1"
        WRITE(LUPRI,*) "STAT=",IER
        WRITE(0,*) "Error in allocation ==> see output file "
        STOP
      ENDIF
      DO I=1,2
      DO J=1,NREFS
        NULLIFY(KVEC(I,J)%P)
      ENDDO
      ENDDO

      NB     = NBL + NBS
!     DO IFRP = 1,NFSYM
        NVECS(IFRP,1) = 0
        CALL ALLOC(KVECTMP,NORB,ID='KVECTMP IN ORBLP1')
        DO I = 1,NREFS
          NSTR(IFRP,0,I) = 1
          CALL NUMLST(VECREF(IFRP,I),KVECTMP, &
               NORB,-NPSH,NESH,               &
               IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(KVECTMP,NSTR(IFRP,0,I), &
               NPSH,NESH,                     &
               NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NVECS(IFRP,1) = NVECS(IFRP,1) + NSTR(IFRP,0,I)
          CALL ALLOCP(KVEC(IFRP,I)%P,NSTR(IFRP,0,I),ID='KVEC%P IN ORBLP1')
          KVEC(IFRP,I)%P(1:NSTR(IFRP,0,I)) = KVECTMP(1:NSTR(IFRP,0,I))
        ENDDO
        CALL DEALLOC(KVECTMP,ID='KVECTMP IN ORBLP1')
!     ENDDO
      NCDIM = NB*NVECS(IFRP,1)*4
      CALL ALLOC(RCOEF,NCDIM,ID="RCOEF in ORBLP1")

!     COEFFICIENT SELECTION
!     =====================

      ICFS = 6
      IF(VAR%IPRPRJ.GE.5) ICFS = 7
      CALL ALLOC(CMO,VAR%N2BBASXQ,ID="CMO in ORBLP1")

      IOFF(1) = 0
      VAR%IPRPRJ = 4
      DO I = 1,NREFS
        INQUIRE(FILE=REFFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
          WRITE(LUPRI,'(A,A6,A)') 'Reference coefficient file ',REFFIL(I),' not present'
          CALL QUIT('PRJANA: No reference coefficients !')
        ENDIF
        CALL OPNFIL(VAR%LUCOEF,REFFIL(I),'OLD','PRJANA')
        IF(VAR%OWNBAS)THEN
          CALL SELOWN_DRV(2,CMO,DDUMARR,IDUMARR,               &
               I,1,VAR%LUCOEF,REFFIL(I),RCOEF,DDUMARR,IDUMARR, &
               NVECS(1,1),IOFF,KVEC(IFRP,I)%P,                 &
               NSTR(1,0,I))
          IF(LOC%IPRLOC.GE.5)THEN
            WRITE(LUPRI,'(A,A6)') '* Selected coefficients for reference ',REFFIL(I)
            NQ = NB*IOFF(IFRP)+1
            CALL PRQMAT(RCOEF(NQ),NB,NSTR(IFRP,0,I),NB,NSTR(IFRP,0,I),4,IPQTOQ(1,0),LUPRI)
          ENDIF
!         DO IFRP = 1, NFSYM
            IOFF(IFRP) = IOFF(IFRP) + NSTR(IFRP,0,I)
!         ENDDO
        ELSE
          WRITE(LUPRI,*) 'Current version work only with .OWNBAS option'
          WRITE(0,*) 'Current version work only with .OWNBAS option'
          STOP
        ENDIF
        CLOSE(VAR%LUCOEF,STATUS='KEEP')
      ENDDO

!     Throw away the full coefficients and eigenvalues
!     ================================================

      CALL DEALLOC(CMO,ID='CMO IN ORBLP1')

!     Do projection analysis
!     ======================

!     NBDIM = 0
!     DO IFRP = 1,NFSYM
        NBDIM = NVECS(IFRP,1)*NI*4
!     ENDDO
      CALL ALLOC(BV,NBDIM,ID="BV in ORBLP1")
!     DO IFRP = 1,NFSYM
!     CALL HEADER('Fermion ircop '//FREP(IFRP),-1)
      NMOL   = NI
      NREF   = NVECS(IFRP,1)
      N2REF  = NREF*NREF
      CALL ALLOC(GVEC,(NREF+1)*2,ID="GVEC in ORBLP1")
      CALL ALLOC(AMAT,N2REF*4,ID="AMAT in ORBLP1")
      CALL ALLOC(JREF,NREF*3,ID="JREF in ORBLP1")
      CALL ALLOC(IPIV,NREF,ID="IPIV in ORBLP1")
      IVEC = 0
      DO I = 1,NREFS
        CALL PRJINF(JREF,IFRP,I,IVEC,NSTR,KVEC(IFRP,I)%P)
        IVEC = IVEC + NSTR(IFRP,0,I)
      ENDDO
      CALL ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NI,NVECS(IFRP,1),BV,NSTR,IPIV,MAXREF)

! Memory allocation for the matrix W = Ac

      CALL ALLOC(WMAT,NREF*NI*4,ID="WMAT in ORBLP1")

! Construction of W = Ac

      IF(LOC%IPRLOC.GE.5)THEN
        WRITE(LUPRI,'(/A)') '* Matrice BV'
        CALL PRQMAT(BV,NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A)') '* Matrice AMAT'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NREF,NI,NREF,D1,                &
           'N','N',IPQTOQ(1,0),AMAT,NREF,NREF,4, &
           'N','N',IPQTOQ(1,0),BV,NREF,NI,4,     &
           D0,IPQTOQ(1,0),WMAT,NREF,NI,4)

      IF(LOC%IPRLOC.GE.5)THEN
        WRITE(LUPRI,'(/A)') '* Matrice W'
        CALL PRQMAT(WMAT,NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF

! Construct P_A for each center A

      JOFF = 0
      DO I = 1,NUCDEP
        CALL QGEMM(NI,NI,NSTR(IFRP,0,I),D1,                    &
                   'H','N',IPQTOQ(1,0),BV(JOFF+1),NREF,NI,4,   &
                   'N','N',IPQTOQ(1,0),WMAT(JOFF+1),NREF,NI,4, &
                   D0,IPQTOQ(1,0),PMAT(1,1,1,I),NI,NI,4)
        JOFF = JOFF + NSTR(IFRP,0,I)
      ENDDO

! P_A matrices symmetrization

      DO I = 1,NUCDEP
        CALL FULMAT('S',NI,NI,PMAT(1,1,1,I))
        DO IZ = 2,4
          CALL FULMAT('A',NI,NI,PMAT(1,1,IZ,I))
        ENDDO
      ENDDO

! Multiply pmat with 1/2 just to avoid confusion in the code below
      DO IZ=1,4
        DO I=1,NI
          DO J=1,NI
            DO K=1,NUCDEP
              PMAT(I,J,IZ,K) = PMAT(I,J,IZ,K) / D2
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      IF(LOC%IPRLOC.GE.2)THEN
        DO IZ=1,4
          WRITE(LUPRI,*)
          WRITE(LUPRI,'(A,I3)') 'Print Q^A matrices, IZ = ',IZ
          DO I=1, NUCDEP
            WRITE(LUPRI,'(A,I4)') ' Nucleus #',I
            DO J=1,NI
              WRITE(LUPRI,'(I4,10ES9.1)') J,(PMAT(J,K,IZ,I),K=1,NI)
            ENDDO
            WRITE(LUPRI,*)
          ENDDO
        ENDDO
      ENDIF

      IF(LOC%IPRLOC.GE.5)THEN
        DO I = 1,NUCDEP
          WRITE(LUPRI,'(/A,I3)') '* ORBLP1: Matrice P  Atom ', I
          WRITE(LUPRI,'(A)') '  ============================'
          CALL PRQMAT(PMAT(1,1,1,I),NI,NI,NI,NI,4,IPQTOQ(1,0),LUPRI)
        ENDDO
      ENDIF

! Calculation and print of the localisation criterion

    WRITE(LUPRI,'(/20X,A)') '** Pipek-Mezey localization criteria **' 
    WRITE(LUPRI,'(20X,A)')  '  ===================================  '
    DINV = D0
    DO I = 1,NI
      QLOC2 = D0
      DO J = 1,NUCDEP
         QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)
      ENDDO
      DLOC = D1/QLOC2
      DINV = DINV + QLOC2
      WRITE(LUPRI,'(11X,A,I4,A,F7.4)')' Orbital',I,' Pipek-Mezey localization measure d = ',DLOC
    ENDDO
    DINVM = DINV/DFLOAT(NI)
    WRITE(LUPRI,'(/11X,A,F21.15)')' Inverse mean delocalization = ',DINVM
    DINVM = - DINVM  ! we are looking for the minimum

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
    DO I=1,NREFS
      CALL DEALLOCP(KVEC(1,I)%P,ID='KVEC%P IN ORBLP1')
    ENDDO
    DEALLOCATE( KVEC )
    CALL DEALLOC(RCOEF,ID='RCOEF IN ORBLP1')
    CALL DEALLOC(BV,ID='BV IN ORBLP1')
    CALL DEALLOC(GVEC,ID='GVEC IN ORBLP1')
    CALL DEALLOC(AMAT,ID='AMAT IN ORBLP1')
    CALL DEALLOC(JREF,ID='JREF IN ORBLP1')
    CALL DEALLOC(IPIV,ID='IPIV IN ORBLP1')
    CALL DEALLOC(WMAT,ID='WMAT IN ORBLP1')

  END SUBROUTINE ORBLP1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

 SUBROUTINE ORBLP2(SMAT,SCOEF,RCOEF,AMAT,NB,NMOL,NREF,BVEC,NSTR,IPIVOT,MAXREF)
!***********************************************************************
!
!     Projection analysis for fermion ircop IFRP
!
!     Written by T.Saue March 26 1999
!
!***********************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER, INTENT(IN ) :: NB, NMOL, NREF, MAXREF
      INTEGER, INTENT(IN ) :: NSTR(2,0:2,0:MAXREF)
      INTEGER, INTENT(OUT) :: IPIVOT(*)

      REAL*8,  INTENT(IN   ) :: SCOEF(*), RCOEF(*)
      REAL*8,  INTENT(INOUT) :: SMAT(*)
      REAL*8,  INTENT(  OUT) :: AMAT(NREF,NREF,4), BVEC(NREF,NMOL,4)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      INTEGER      :: I, J, IFRP, JOFF, KOFF, JOB, NEFF, NREFS, LUPRI
      INTEGER      :: IPQTOQ(4,0:7)

      REAL*8 :: PTOL
      REAL*8, PARAMETER :: DC=1.0D-2

      LOGICAL DOPOP

  !------------------------------------------------------------------------------!
      IFRP = 1
      NREFS = VAR%NREFS

      DOPOP = NREFS.GT.1
      PTOL  = VAR%PROTHR*DC
      LUPRI = VAR%LUPRI
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)
      VECREF(:,:) = VAR%VECREF(:,:)
      REFFIL(:)   = VAR%REFFIL(:)

!     Make projection vector

      CALL QTRANS90('AOMO','S',D0,NB,NB,NREF,NMOL, &
                    SMAT,NB,NB,1,IPQTOQ(1,0),      &
                    BVEC,NREF,NMOL,4,IPQTOQ(1,0),  &
                    RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                    SCOEF,NB,NMOL,4,IPQTOQ(1,0),VAR%IPRPRJ)
      IF(NREFS.GT.1) THEN

!       Make overlap matrix

        CALL QTRANS90('AOMO','S',D0,NB,NB,NREF,NREF, &
                      SMAT,NB,NB,1,IPQTOQ(1,0),      &
                      AMAT,NREF,NREF,4,IPQTOQ(1,0),  &
                      RCOEF,NB,NREF,4,IPQTOQ(1,0),   &
                      RCOEF,NB,NREF,4,IPQTOQ(1,0),VAR%IPRPRJ)

!       Print overlaps

        IF(LOC%IPRLOC.GE.5)THEN
           JOFF = 1
           DO J = 1,NREFS
              KOFF = 1
              DO I = 1,(J-1)
                 WRITE(LUPRI,'(A,2(3X,A6))') '*Overlaps :',REFFIL(I),REFFIL(J)
                 CALL PRQMAT(AMAT(KOFF,JOFF,1),                  &
                      NSTR(IFRP,0,I),NSTR(IFRP,0,J),NREF,NREF,4, &
                      IPQTOQ(1,0),LUPRI)
                 KOFF = KOFF + NSTR(IFRP,0,I)
              ENDDO
              JOFF = JOFF + NSTR(IFRP,0,J)
           ENDDO
        ENDIF

!       Solve linear system by Cholesky decomposition
!       with full pivoting
!        --> notice that overlap matrix SMAT is destroyed !

        JOB = 1
        CALL QCHOLD(AMAT,NREF,4,NREF,NREF,SMAT,PTOL,NEFF,JOB,IPIVOT)
        CALL QCHOLS(AMAT,NREF,NEFF,NMOL,4,NREF,NREF,SMAT,BVEC, &
             NREF,NMOL,BVEC,NREF,NMOL,JOB,IPIVOT,SMAT(1+NREF))
        IF(NEFF.LT.NREF) THEN
          WRITE(LUPRI,'(A,A)') '* WARNING: ', &
           ' linear dependencies detected in Cholesky decomposition'
          WRITE(LUPRI,'(11X,A,I5,A,I5)') ' Reduced dimensionality: ',NREF,' --> ',NEFF
        ENDIF

!       Restore overlap matrix

        CALL QHMRST(AMAT,NREF,VAR%NZ,NREF,NREF)
      ENDIF

      END SUBROUTINE ORBLP2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO1(SMAT,SCOEF,PMAT,NBL,NBS,NI,NUCDEP,NFNUC,DINVM)
!***********************************************************************
!
!     Written by S.Dubillard Jan 21 2004
!     Last revision:
!
!     Construction of the matrices P_A=c_A^+Sc from Mulliken analysis
!
!     Calculation of the localisation criteria DINV
!
!***********************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER, INTENT(IN ) :: NBL, NBS, NI, NUCDEP
      REAL*8,  INTENT(OUT) :: DINVM

      INTEGER, INTENT(IN ) :: NFNUC(*)
      REAL*8,  INTENT(IN ) :: SMAT(*), SCOEF(*)
      REAL*8,  INTENT(OUT) :: PMAT(NI,NI,4,NUCDEP)

      ! Internal variables
      INTEGER :: I, J, K, IZ, NB, IOFFL, IOFFS, LUPRI
      INTEGER :: IPQTOQ(4,0:7)
      REAL*8  :: DINV, QLOC2, DLOC

      REAL*8, ALLOCATABLE, DIMENSION(:) :: WLMAT, WSMAT

  !------------------------------------------------------------------------------!
      NB = NBL + NBS

      LUPRI       = VAR%LUPRI
      IPQTOQ(:,:) = VAR%IPQTOQ(:,:)

! Memory allocation for the matrix W = Sc

      CALL ALLOC(WLMAT,NBL*NI*4,ID='WLMAT IN ORBLO1')
      CALL ALLOC(WSMAT,NBS*NI*4,ID='WSMAT IN ORBLO1')

! Construction of W = Sc (large and small components)

      CALL QGEMM(NBL,NI,NBL,D1,               &
           'N','N',IPQTOQ(1,0),SMAT,NB,NB,1,  &
           'N','N',IPQTOQ(1,0),SCOEF,NB,NI,4, &
           D0,IPQTOQ(1,0),WLMAT,NBL,NI,4)
      IF(.NOT.VAR%TWOCOMPBSS) THEN
         CALL QGEMM(NBS,NI,NBS,D1,                            &
              'N','N',IPQTOQ(1,0),SMAT(1+NB*NBL+NBL),NB,NB,1, &
              'N','N',IPQTOQ(1,0),SCOEF(1+NBL),NB,NI,4,       &
              D0,IPQTOQ(1,0),WSMAT,NBS,NI,4)
      ENDIF

      IF(LOC%IPRLOC.GE.5)THEN
         WRITE(LUPRI,'(/A)') '* WL matrix :'
         CALL PRQMAT(WLMAT,NBL,NI,NBL,NI,4,IPQTOQ(1,0),LUPRI)
         IF(.NOT.VAR%TWOCOMPBSS) THEN
            WRITE(LUPRI,'(/A)') '* WS matrix :'
            CALL PRQMAT(WSMAT,NBS,NI,NBS,NI,4,IPQTOQ(1,0),LUPRI)
         ENDIF
      ENDIF

! Construct P_A for each center A

      IOFFL = 0
      IOFFS = 0
      DO I = 1,NUCDEP
         CALL QGEMM(NI,NI,NFNUC(I),D1,                     &
              'H','N',IPQTOQ(1,0),SCOEF(1+IOFFL),NB,NI,4,  &
              'N','N',IPQTOQ(1,0),WLMAT(1+IOFFL),NBL,NI,4, &
              D0,IPQTOQ(1,0),PMAT(1,1,1,I),NI,NI,4)
         IF(.NOT.VAR%TWOCOMPBSS) THEN
            CALL QGEMM(NI,NI,NFNUC(I+NUCDEP),D1,                 &
                 'H','N',IPQTOQ(1,0),SCOEF(1+NBL+IOFFS),NB,NI,4, &
                 'N','N',IPQTOQ(1,0),WSMAT(1+IOFFS),NBS,NI,4,    &
                 D1,IPQTOQ(1,0),PMAT(1,1,1,I),NI,NI,4)
            IOFFS = IOFFS + NFNUC(NUCDEP+I)
         ENDIF
         IOFFL = IOFFL + NFNUC(I)
      ENDDO

! P_A matrices symetrisation

      DO I = 1,NUCDEP
         CALL FULMAT('S',NI,NI,PMAT(1,1,1,I))
         DO IZ = 2,4
            CALL FULMAT('A',NI,NI,PMAT(1,1,IZ,I))
         ENDDO
      ENDDO

! Multiply pmat with 1/2 just to avoid confusion in the code below                                
      DO IZ=1,4                                                                                   
        DO I=1,NI                                                                                 
          DO J=1,NI                                                                               
            DO K=1,NUCDEP                                                                         
              PMAT(I,J,IZ,K) = PMAT(I,J,IZ,K) / D2                                                
            ENDDO                                                                                 
          ENDDO                                                                                   
        ENDDO                                                                                     
      ENDDO

      IF(LOC%IPRLOC.GE.2)THEN                                                                     
        DO IZ=1,4                                                                                 
          WRITE(LUPRI,*)                                                                          
          WRITE(LUPRI,'(A,I3)') 'Print Q^A matrices, IZ = ',IZ                                    
          DO I=1, NUCDEP                                                                          
            WRITE(LUPRI,'(A,I4)') ' Nucleus #',I                                                  
            DO J=1,NI                                                                             
              WRITE(LUPRI,'(I4,10ES9.1)') J,(PMAT(J,K,IZ,I),K=1,NI)                               
            ENDDO                                                                                 
            WRITE(LUPRI,*)                                                                        
          ENDDO                                                                                   
        ENDDO                                                                                     
      ENDIF

      IF(LOC%IPRLOC.GE.5)THEN
         DO I = 1,NUCDEP
            WRITE(LUPRI,'(/A,I3)') '* ORBLO1: Matrice P  Atom ', I
            WRITE(LUPRI,'(A)') '  ============================'
            CALL PRQMAT(PMAT(1,1,1,I),NI,NI,NI,NI,4,IPQTOQ(1,0),LUPRI)
         ENDDO
      ENDIF

! Calculation and print of the localisation criterion

    WRITE(LUPRI,'(/20X,A)') '** Pipek-Mezey localization criteria **' 
    WRITE(LUPRI,'(20X,A)')  '  ===================================  '
    DINV = D0
    DO I = 1,NI
      QLOC2 = D0
      DO J = 1,NUCDEP
         QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)
      ENDDO
      DLOC = D1/QLOC2
      DINV = DINV + QLOC2
      WRITE(LUPRI,'(11X,A,I4,A,F7.4)')' Orbital',I,' Pipek-Mezey localization measure d = ',DLOC
    ENDDO
    DINVM = DINV/DFLOAT(NI)
    WRITE(LUPRI,'(/11X,A,F21.15)')' Inverse mean delocalization = ',DINVM
    DINVM = - DINVM  ! we are looking for the minimum

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
    CALL DEALLOC(WLMAT,ID='WLMAT IN ORBLO1')
    CALL DEALLOC(WSMAT,ID='WSMAT IN ORBLO1')

  END SUBROUTINE ORBLO1

end module func_PipekMezey
