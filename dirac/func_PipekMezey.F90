module func_PipekMezey
!================================================================================!
!!! This is not a true module. It depends on external subroutines and .h files !!!
!                                                                                !
!                  This module is used in optimization module                    !
!           it will provide information about Pipek-Mezey functional             !
!                                                                                !
!   Two forms of the Pipek-Mezey functional are available:                       !
!   1) based on Mulliken population analysis  -->  orblo1                        !
!   2) based on projection analysis           -->  ORBLP1                        !
!                                                                                !
!   [1] S. Dubillard, J.-B. Rota, and T. Saue, J.Chem.Phys. 124, 154307 (2006)   !
!================================================================================!
use memory_allocator
implicit none

public  :: init_func, quit_func, mv2point, get_func, get_gradient, get_hessian, get_date

private :: ORBLP1, ORBLP2, orblo1, ORBLO4, ORBLO5, ROTCMO_MOD

type :: PipekMezey_flags
  integer :: iprloc
  logical :: prjloc
end type PipekMezey_flags

type :: PipekMezey_variales
  integer :: lupri, lucoef
  integer :: nfbas(2,0:2), n2bbasx, n2bbasxq
  integer :: nish(2), norb(2), npsh(2), nesh(2), norbt
  real*8  :: ssmtrc
  logical :: h2comp, levyle
  integer :: nz, ipqtoq(4,0:7)
  integer :: nucdep
  real*8  :: thrzer
  integer :: nrefs, iprprj
  integer :: maxref
  real*8  :: prothr
  logical :: ownbas
  character*6  :: reffil(15)    ! reffil(maxref)    -->  see dcbprj.h
  character*72 :: vecref(2,15)  ! vecref(2,maxref)  -->  see dcbprj.h
  real*8  :: timstr
end type PipekMezey_variales

type( PipekMezey_flags     ), private :: PM
type( PipekMezey_variales  ), private :: var

! Some global parameters
real*8, parameter, private :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

! functional value of our functional (stored for output)
real*8, save :: dinvmg4out

! Global arrays
real*8,  allocatable, save, private :: pmat(:,:,:,:), smat(:), scoef(:), storec(:)
integer, allocatable, save, private :: nfnuc(:)

contains

  subroutine init_func(ndim)
  !==============================================================================!
  !            Initialize module func_PipekMezey for Dirac program               !
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

  ! External variables
  integer :: ndim

  ! Internal variables
  logical :: tobe
  integer :: i, ifrp, idum
  real*8  :: toterg, second

  integer, allocatable, dimension(:) :: isel, indao, indaos, indnao
  real*8,  allocatable, dimension(:) :: cmo

  !------------------------------------------------------------------------------!
  ! Start localization of Molecular orbitals
  call titler('Molecular orbital localization','*',122)

  !------------------------------------------------------------------------------!
  !                       Initialize global variables                            !
  !------------------------------------------------------------------------------!
  ! dcbloc.h
  PM%iprloc      = iprloc
  PM%prjloc      = prjloc
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

  ! two component calculation
  var%h2comp = twocomp .or. x2c

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

  ! starting time
  var%timstr = second()

  ! dimension of your space
  ifrp = 1
  ndim = 4 * ( nish(ifrp)*nish(ifrp) - nish(ifrp) ) / 2

  !------------------------------------------------------------------------------!
  !                             Print debug info                                 !
  !------------------------------------------------------------------------------!
  if(PM%iprloc.ge.4)then
    write(lupri,*) 'lupri    =',lupri
    write(lupri,*) 'nish     =',nish(ifrp)
    write(lupri,*) 'norb     =',norb(ifrp)
    write(lupri,*) 'nucdep   =',nucdep
    write(lupri,*) 'maxref   =',var%maxref
    write(lupri,*) 'thrzer   =',thrzer
    write(lupri,*) 'nfbas(0) =',nfbas(ifrp,0)
    write(lupri,*) 'nfbas(1) =',nfbas(ifrp,1)
    write(lupri,*) 'nfbas(2) =',nfbas(ifrp,2)
  endif

  !------------------------------------------------------------------------------!
  !                              Allocate arrays                                 !
  !------------------------------------------------------------------------------!
  ! allocate global arrays
  call alloc(scoef ,var%nish(ifrp)*var%nfbas(ifrp,0)*4,        id='scoef  in init_func')
  call alloc(storec,var%nish(ifrp)*var%nfbas(ifrp,0)*4,        id='storec in init_func')
  call alloc(pmat  ,var%nish(ifrp),var%nish(ifrp),4,var%nucdep,id='pmat   in init_func')
  call alloc(smat  ,var%n2bbasx,                               id='smat   in init_func')
  call alloc(nfnuc ,2*var%nucdep,                              id='nfnuc  in init_func')
  scoef = d0; storec = d0; pmat = d0; smat = d0; nfnuc = 0

  ! allocate local arrays
  call alloc(cmo   ,var%n2bbasxq,       id='cmo    in init_func')

  call alloc(isel  ,var%nish(ifrp),     id='isel   in init_func')
  call alloc(indao ,var%nfbas(ifrp,0),  id='indao  in init_func')
  call alloc(indaos,var%nfbas(ifrp,0),  id='indaos in init_func')
  call alloc(indnao,2*var%nucdep+2,     id='indnao in init_func')
  cmo = d0; indao = 0; indaos = 0; indnao = 0

  !------------------------------------------------------------------------------!
  !                      Precalculate some global arrays                         !
  !------------------------------------------------------------------------------!
  !----------------------------------------------------!
  ! -- nfnuc -- Get the number of functions per center !
  !----------------------------------------------------!
  call indaoc(indao,indaos,indnao)
  call orblo4(indnao,var%nucdep)

  !-------------------------------!
  ! -- smat -- Get overlap matrix !
  !-------------------------------!
  call gtovlx(smat,var%ssmtrc)
  if(PM%iprloc.ge.6)then
    write(lupri,'(/a)') '* orbloc: Overlap matrix' 
    call prqmat(smat,var%nfbas(ifrp,0),var%nfbas(ifrp,0),var%nfbas(ifrp,0),var%nfbas(ifrp,0), &
                1,ipqtoq(1,0),lupri)
  endif

  !-----------------------------------------------------!
  ! -- scoef -- Electronic orbital coeffients selection !
  !-----------------------------------------------------!
  ! Check if coefficients are on file
  inquire(file='DFCOEF',exist=tobe)
  if(.not.tobe)then
    write(lupri,'(a)') 'init_func: Coefficient file not found !'
    call quit('init_func: Coefficients not found !')
  endif

  ! Read coefficients
  call opnfil(var%lucoef,'DFCOEF','OLD','PAMANA')
  call reacmo_no_work(var%lucoef,'DFCOEF',cmo,idum,idum,toterg,2)
  close(var%lucoef,status='keep')

  ! select orbitals
  call occupy(isel,var%nish(ifrp))
  call selcfs(cmo,ifrp,scoef,var%nish(ifrp),isel,0, &
              var%nish(ifrp),var%nfbas(ifrp,0),var%norb(ifrp))
  if(PM%iprloc.ge.5)then
    write(lupri,'(/a)') '* orblo1: Selected molecular coefficients'
    call prqmat(scoef,var%nfbas(ifrp,0),var%nish(ifrp),var%nfbas(ifrp,0), &
                var%nish(ifrp),4,ipqtoq(1,0),lupri)
  endif

  !------------------------------------------------------------------------------!
  !                         Deallocate local arrays                              !
  !------------------------------------------------------------------------------!
  call dealloc(cmo   ,id='cmo    in init_func')

  call dealloc(isel  ,id='isel   in init_func')
  call dealloc(indao ,id='indao  in init_func')
  call dealloc(indaos,id='indaos in init_func')
  call dealloc(indnao,id='indnao in init_func')

  end subroutine init_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine mv2point(restore,rejected,xtmp)
  !==============================================================================!
  !                    Move to the new expansion point                           !
  !                      --> get new MO coefficients                             !
  !                                                                              !
  !   Note:  Deffinition of kappa matrix in rotcmo_mod is different as in        !
  !          calculation of gradient and hessian in get_gradient and get_hessian !
  !          --> Indexes are exchanged and j quaterion is with minus sign        !
  !              Qkappa_lk = kappa_kl - bar{kappa}_kl^* j                        !
  !              It follows from definition of bared unitary matrix:             !
  !              tilde{Psi_i} = U_ij Psi_j + bar{U}_ij Psi_j                     !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8,  intent(inout) :: scoef(*), storec(*)
  ! real*8,  intent(  out) :: pmat(nish,nish,4,nucdep)

  ! External variables
  logical, intent(in) :: restore, rejected
  real*8,  intent(in) :: xtmp(:)

  ! Internal variables
  integer :: i, j, ij1, ij2, ij3, ij4, iz, nish, ni12, ifrp

  real*8, allocatable :: xkappa(:,:,:)

  !------------------------------------------------------------------------------!
  ifrp = 1
  nish = size(pmat,1)
  ni12 = size(xtmp,1) / 4

  ! allocate local variable
  call alloc(xkappa,nish,nish,4,id='xkappa in mv2point')

  ! Unpack X to kappa matrix
  xkappa(:,:,:) = d0
  ij1 = 0
  do j = 1,nish
    do i = j+1,nish
      ij1 = ij1 + 1                                                                             
      ij2 = ij1 + ni12
      ij3 = ij2 + ni12                                                                          
      ij4 = ij3 + ni12
      xkappa(j,i,1) =   xtmp(ij1)
      xkappa(i,j,1) = - xkappa(j,i,1)

      xkappa(j,i,2) =   xtmp(ij3)
      xkappa(i,j,2) =   xkappa(j,i,2)

      xkappa(j,i,3) = - xtmp(ij2)
      xkappa(i,j,3) =   xkappa(j,i,3)

      xkappa(j,i,4) =   xtmp(ij4)
      xkappa(i,j,4) =   xkappa(j,i,4)
    enddo
  enddo

  if(PM%iprloc.ge.2)then
    do iz=1,4
      write(var%lupri,'(/a,i4)') 'Print xkappa matrix, IZ =',iz
      do i=1,nish
        write(var%lupri,'(i4,10es9.1)') i,(xkappa(i,j,iz),j=1,nish)
      enddo
    enddo
  endif

  ! get new MO coefficients(scoef)
  if(restore)then
    scoef(:) = storec(:)  ! restore MO coefficients
  elseif(rejected)then
    scoef(:) = storec(:)  ! restore MO coefficients; the step was rejected
  elseif(.not.rejected)then
    storec(:) = scoef(:)  ! store MO coefficients
  endif
  call rotcmo_mod(var%nfbas(ifrp,0),nish,4,xkappa,nish,var%ipqtoq(ifrp,0))

  ! deallocate local variable
  call dealloc(xkappa,id='xkappa in mv2point')

  !------------------------------------------------------------------------------!
  end subroutine mv2point

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  real*8 function get_func()
  !==============================================================================!
  !                       Get value of your functional                           !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8,  intent(inout) :: smat(*)
  ! real*8,  intent(  out) :: pmat(nish,nish,4,nucdep)

  ! Internal variables
  integer :: nish, nucdep
  real*8  :: dinvmg

  !------------------------------------------------------------------------------!
  nish   = size(pmat,1)
  nucdep = size(pmat,4)

  if(PM%prjloc)then ! localization with projection analysis
    call gtovlx(smat,var%ssmtrc)
    call orblp1(var%nfbas(1,1),var%nfbas(1,2),nish,nucdep,var%maxref,dinvmg)
  else               ! localization with Mulliken analysis
    call orblo1(var%nfbas(1,1),var%nfbas(1,2),nish,nucdep,dinvmg)
  endif
  dinvmg4out = dinvmg   ! global variable
  get_func   = dinvmg

  !------------------------------------------------------------------------------!
  end function get_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine quit_func(hmat)
  !==============================================================================!
  !                  Uninitialize module func_PipekMezey                         !
  !==============================================================================!
  implicit none

  ! External variables
  real*8,  intent(inout), pointer :: hmat(:,:)

  ! Internal variables
  character :: resst*50
  integer   :: ifrp
  real*8    :: timend, second

  ifrp = 1
  ! Put the localised coefficients in DFCOEF (scoef)
  call orblo5(var%nfbas(ifrp,0),var%nish(ifrp),var%norb(ifrp))

  ! deallocate global arrays
  call dealloc(storec,id='storec in quit_func')
  call dealloc(scoef ,id='scoef  in quit_func')
  call dealloc(pmat  ,id='pmat   in quit_func')
  call dealloc(smat  ,id='smat   in quit_func')
  call dealloc(nfnuc ,id='nfnuc  in quit_func')

  ! deallocate pointers
  call deallocp(hmat ,id='hmat   in quit_func')

  ! print final functional value
  write(resst,'(a,f21.15)') 'Inverse mean delocalization =',dinvmg4out
  call titler(resst,'*',104)

  ! calculate time spend on localization
  timend = second()
  call timtxt('          --- Time used in localization is :',timend-var%timstr,var%lupri)

  end subroutine quit_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_hessian(hmat,full,totdiag,totnond)
  !==============================================================================!
  !                                                                              !
  !  Construction of full hessian                                                !
  !                                                                              ! 
  !     HMAT         A = d^2 G / dK* dK     B = d^2 G / dK dK          K         !
  !  / A   B  \       / A_uu   A_ub \        / B_uu   B_ub \      / K(unbared) \ !
  !  \ B*  A* /       \ A_bu   A_bb /        \ B_bu   B_bb /      \ K(bared)   / !
  !                                                                              ! 
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8, intent(in) :: pmat(nish,nish,4,nucdep)

  ! External variables
  logical, intent(in ) :: full
  real*8,  intent(out) :: totdiag, totnond

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
      hmat(i,1) = - hmat(i,1) / dfloat(ni) ! minus sign --> we are looking for the minimum
      totdiag = totdiag + dabs( hmat(i,1) )
    enddo

    if(PM%iprloc.ge.3)then
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
        hmat(i,j) = - hmat(i,j) / dfloat(ni)  ! minus sign --> we are looking for the minimum
        if(i==j)then
          totdiag = totdiag + dabs( hmat(i,j) )
        else
          totnond = totnond + dabs( hmat(i,j) )
        endif
      enddo
    enddo
    totnond = d2*totnond

    !----------------------------------------------------------------------------!
    !                              PRINT section                                 !
    !----------------------------------------------------------------------------!
    if(PM%iprloc.ge.3)then
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

  subroutine get_gradient(gradient,dnrmg)
  !==============================================================================!
  !                     Construction of the gradient                             !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8, intent(in) :: pmat(nish,nish,4,nucdep)

  ! External variables
  real*8, intent(out) :: dnrmg

  real*8, intent(out), dimension(:) :: gradient

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
    gradient(i) = - gradient(i) / dfloat(ni)  ! minus sign --> we are looking for the minimum
  enddo

  if(PM%iprloc.ge.2)then
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

  subroutine get_date(daytid)
  implicit none
  character, intent(out) :: daytid*24

  call gtinfo(daytid)

  end subroutine get_date
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLP1(NBL,NBS,NI,NUCDEP,MAXREF,DINVM)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
      USE SELECT_COEFF_OWN_BAS
      IMPLICIT NONE

      ! Global variables
      ! REAL*8,  INTENT(  OUT) :: PMAT(NI,NI,4,NUCDEP)

      ! External variables
      INTEGER, INTENT(IN ) :: NBL, NBS, NI, NUCDEP, MAXREF
      REAL*8,  INTENT(OUT) :: DINVM

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
          IF(PM%IPRLOC.GE.5)THEN
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
      CALL ORBLP2(RCOEF,AMAT,NB,NI,NVECS(IFRP,1),BV,NSTR,IPIV,MAXREF)

! Memory allocation for the matrix W = Ac

      CALL ALLOC(WMAT,NREF*NI*4,ID="WMAT in ORBLP1")

! Construction of W = Ac

      IF(PM%IPRLOC.GE.5)THEN
        WRITE(LUPRI,'(/A)') '* Matrice BV'
        CALL PRQMAT(BV,NREF,NI,NREF,NI,4,IPQTOQ(1,0),LUPRI)
        WRITE(LUPRI,'(/A)') '* Matrice AMAT'
        CALL PRQMAT(AMAT,NREF,NREF,NREF,NI,4,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL QGEMM(NREF,NI,NREF,D1,                &
           'N','N',IPQTOQ(1,0),AMAT,NREF,NREF,4, &
           'N','N',IPQTOQ(1,0),BV,NREF,NI,4,     &
           D0,IPQTOQ(1,0),WMAT,NREF,NI,4)

      IF(PM%IPRLOC.GE.5)THEN
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

      IF(PM%IPRLOC.GE.2)THEN
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

      IF(PM%IPRLOC.GE.5)THEN
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

  SUBROUTINE ORBLP2(RCOEF,AMAT,NB,NMOL,NREF,BVEC,NSTR,IPIVOT,MAXREF)
!***********************************************************************
!
!     Projection analysis for fermion ircop IFRP
!
!     Written by T.Saue March 26 1999
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! REAL*8,  INTENT(IN   ) :: SCOEF(*)
      ! REAL*8,  INTENT(INOUT) :: SMAT(*)

      ! External variables
      INTEGER, INTENT(IN ) :: NB, NMOL, NREF, MAXREF
      INTEGER, INTENT(IN ) :: NSTR(2,0:2,0:MAXREF)
      INTEGER, INTENT(OUT) :: IPIVOT(*)

      REAL*8,  INTENT(IN   ) :: RCOEF(*)
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

        IF(PM%IPRLOC.GE.5)THEN
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

  subroutine orblo1(nbl,nbs,ni,nucdep,dinvm)
  !==============================================================================!
  !                                                                              !
  !    Construction of Q^A_i matrices from Mulliken analysis (eq. 3 in [1])      !
  !                         Q^A_i == PMAT(i,i,4,A) =                             !
  !     = 1/2 ( SUM_{nu=A} SUM_{mu} C_{nu,i}^* S_{nu,mu} C_{mu,i} + c.c. )       !
  !                                                                              !
  !              Calculation of the inverse mean delocalization                  !
  !                  DINVM = 1/M  SUM_i^M SUM_A ( Q^A_i )^2                      !
  !==============================================================================!
  implicit none

  ! Global variables
  ! real*8,  intent(in ) :: smat(*), scoef(*)
  ! real*8,  intent(out) :: pmat(ni,ni,4,nucdep)
  ! integer, intent(in ) :: nfnuc(*)

  ! External variables
  integer, intent(in ) :: nbl, nbs, ni, nucdep
  real*8,  intent(out) :: dinvm

  ! Internal variables
  integer :: i, j, k, iz, nb, ioffl, ioffs, lupri
  integer :: ipqtoq(4,0:7)
  real*8  :: dinv, qloc2, dloc

  real*8, allocatable, dimension(:) :: wlmat, wsmat

  !------------------------------------------------------------------------------!
  nb = nbl + nbs

  lupri       = var%lupri
  ipqtoq(:,:) = var%ipqtoq(:,:)

  ! Memory allocation for the matrix W = Sc
  call alloc(wlmat,nbl*ni*4,id='wlmat in orblo1')
  call alloc(wsmat,nbs*ni*4,id='wsmat in orblo1')

  ! Construction of W = Sc (large and small components)
  call qgemm(nbl,ni,nbl,d1,                     &
             'N','N',ipqtoq(1,0),smat,nb,nb,1,  &
             'N','N',ipqtoq(1,0),scoef,nb,ni,4, &
              d0,ipqtoq(1,0),wlmat,nbl,ni,4)
  if(.not.var%h2comp)then
    call qgemm(nbs,ni,nbs,d1,                                  &
               'N','N',ipqtoq(1,0),smat(1+nb*nbl+nbl),nb,nb,1, &
               'N','N',ipqtoq(1,0),scoef(1+nbl),nb,ni,4,       &
                d0,ipqtoq(1,0),wsmat,nbs,ni,4)
  endif

  if(pm%iprloc.ge.5)then
    write(lupri,'(/a)') '* WL matrix :'
    call prqmat(wlmat,nbl,ni,nbl,ni,4,ipqtoq(1,0),lupri)
    if(.not.var%h2comp)then
      write(lupri,'(/a)') '* WS matrix :'
      call prqmat(wsmat,nbs,ni,nbs,ni,4,ipqtoq(1,0),lupri)
    endif
  endif

  ! Construct P_A for each center A
  ioffl = 0
  ioffs = 0
  do i = 1,nucdep
    call qgemm(ni,ni,nfnuc(i),d1,                           &
               'H','N',ipqtoq(1,0),scoef(1+ioffl),nb,ni,4,  &
               'N','N',ipqtoq(1,0),wlmat(1+ioffl),nbl,ni,4, &
                d0,ipqtoq(1,0),pmat(1,1,1,i),ni,ni,4)
    if(.not.var%h2comp)then
      call qgemm(ni,ni,nfnuc(i+nucdep),d1,                       &
                 'H','N',ipqtoq(1,0),scoef(1+nbl+ioffs),nb,ni,4, &
                 'N','N',ipqtoq(1,0),wsmat(1+ioffs),nbs,ni,4,    &
                 d1,ipqtoq(1,0),pmat(1,1,1,i),ni,ni,4)
      ioffs = ioffs + nfnuc(nucdep+i)
    endif
    ioffl = ioffl + nfnuc(i)
  enddo

  ! P_A matrices symetrisation
  do i = 1,nucdep
    call fulmat('S',ni,ni,pmat(1,1,1,i))
    do iz = 2,4
      call fulmat('A',ni,ni,pmat(1,1,iz,i))
    enddo
  enddo

  ! Multiply pmat with 1/2 just to avoid confusion in the code below
  do iz=1,4
    do i=1,ni
      do j=1,ni
        do k=1,nucdep
          pmat(i,j,iz,k) = pmat(i,j,iz,k) / d2
        enddo
      enddo
    enddo
  enddo

  if(pm%iprloc.ge.2)then
    do iz=1,4
      write(lupri,'(/a,i3)') 'Print Q^A matrices, IZ = ',iz
      do i=1,nucdep
        write(lupri,'(a,i4)') ' Nucleus #',i
        do j=1,ni
          write(lupri,'(i4,10es9.1)') j,(pmat(j,k,iz,i),k=1,ni)
        enddo
        write(lupri,*)
      enddo
    enddo
  endif

  if(pm%iprloc.ge.5)then
     do i = 1,nucdep
        write(lupri,'(/a,i3)') '* ORBLO1: Matrice P  Atom ',i
        write(lupri,'(a)') '  ============================'
        call prqmat(pmat(1,1,1,i),ni,ni,ni,ni,4,ipqtoq(1,0),lupri)
     enddo
  endif

  ! Calculation and print of the localisation criterion
  write(lupri,'(/20x,a)') '** Pipek-Mezey localization criteria **' 
  write(lupri,'(20x,a)')  '  ===================================  '
  dinv = d0
  do i = 1,ni
    qloc2 = d0
    do j = 1,nucdep
       qloc2 = qloc2 + pmat(i,i,1,j)*pmat(i,i,1,j)
    enddo
    dloc = d1/qloc2
    dinv = dinv + qloc2
    write(lupri,'(11x,a,i4,a,f7.4)')' Orbital',i,' Pipek-Mezey localization measure d = ',dloc
  enddo
  dinvm = dinv/dfloat(ni)
  write(lupri,'(/11x,a,f21.15)')' Inverse mean delocalization = ',dinvm
  dinvm = - dinvm  ! we are looking for the minimum

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
  call dealloc(wlmat,id='wlmat in orblo1')
  call dealloc(wsmat,id='wsmat in orblo1')

  end subroutine orblo1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO4(INDNAO,NUCDEP)
!***********************************************************************
!
!     Written by S.Dubillard Jan 21 2004
!     Last revision:
!
!     Construction of the array NFNUC which contains the number of
!     large components basis functions for each center and then the 
!     number of small components basis functions for each center
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! INTEGER, INTENT(OUT) :: NFNUC(2*NUCDEP)

      ! External variables
      INTEGER, INTENT(IN ) :: NUCDEP
      INTEGER, INTENT(IN ) :: INDNAO(2,*)

      ! Internal variables
      INTEGER :: IFRP, I

  !------------------------------------------------------------------------------!
      IFRP = 1
      DO I = 1,NUCDEP-1
         NFNUC(I)=INDNAO(1,I+1)-INDNAO(1,I)
      ENDDO
      NFNUC(NUCDEP)=INDNAO(2,1)-INDNAO(1,NUCDEP)
      DO I = 1,NUCDEP-1
         NFNUC(I+NUCDEP)=INDNAO(2,I+1)-INDNAO(2,I)
      ENDDO
      NFNUC(2*NUCDEP)=VAR%NFBAS(IFRP,0)-INDNAO(2,NUCDEP)

  END SUBROUTINE ORBLO4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO5(NB,NI,NO)
!***********************************************************************
!
!     Written by S.Dubillard Apr 13 2004
!     Last revision:
!
!     Put the localised coefficients into the DFCOEF file (1/2)
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! REAL*8  :: SCOEF(*) ! SCOEF(NB,NI,4)

      ! External variables
      INTEGER :: NB, NI, NO

      ! Internal variables
      INTEGER :: I, J, IZ, L, IJK
      REAL*8  :: TOTERG

      INTEGER, ALLOCATABLE :: IBEIG(:)
      REAL*8,  ALLOCATABLE :: CMO(:,:,:), EIG(:)

  !------------------------------------------------------------------------------!
      ! allocate arrays
      CALL ALLOC(CMO  ,NB,NO,4  ,ID='CMO   IN ORBLO5')
      CALL ALLOC(EIG  ,VAR%NORBT,ID='EIG   IN ORBLO5')
      CALL ALLOC(IBEIG,VAR%NORBT,ID='IBEIG IN ORBLO5')

      ! open DFCOEF file and read EIG, IBEIG and TOTERG
      CALL OPNFIL(VAR%LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO_no_work(VAR%LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)

      ! reorganize MO
      DO IZ = 1,4
        DO J = 1,NI
          L = J + NO/2
          IF (VAR%H2COMP .OR. VAR%LEVYLE) L = J
          DO I = 1,NB
            IJK = NB*NI*(IZ-1) + NB*(J-1) + I
            CMO(I,L,IZ) = SCOEF(IJK)
          ENDDO
        ENDDO
      ENDDO
      IF (PM%IPRLOC.GE.5)THEN
         WRITE(VAR%LUPRI,'(/A)') '*ORBLO6 : New CMO matrix'
         CALL PRQMAT(CMO,NB,NO,NB,NO,4,VAR%IPQTOQ(1,0),VAR%LUPRI)
      ENDIF

      ! write MO to file and close him
      CALL WRICMO(VAR%LUCOEF,CMO,EIG,IBEIG,TOTERG)
      CLOSE(VAR%LUCOEF,STATUS='KEEP')

      ! deallocate arrays
      CALL DEALLOC(CMO  ,ID='CMO   IN ORBLO5')
      CALL DEALLOC(EIG  ,ID='EIG   IN ORBLO5')
      CALL DEALLOC(IBEIG,ID='IBEIG IN ORBLO5')

  END SUBROUTINE ORBLO5

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ROTCMO_MOD(NBASI,NORBI,NZ,XKAPPA,NORBT,IP)
!***********************************************************************
!
!     Rotate coefficients
!
!     SCOEF_{new} = SCOEF_{old} exp( -\kappa )
!     We use a trick a la Sostrup notes, chapter 3.1.5
!
!     Input:
!        SCOEF   - old coefficients
!        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
!        IPRINT  - print level
!
!     Output:
!        SCOEF   - new coefficients
!
!     Written by J. Thyssen - Nov 24 1998
!
!***********************************************************************
      IMPLICIT NONE

      ! Global variables
      ! REAL*8,  INTENT(INOUT) :: SCOEF(*)  ! SCOEF(NBASI,NORBI,NZ)

      ! External variables
      INTEGER, INTENT(IN   ) :: NBASI, NORBI, NZ, NORBT, IP(4)
      REAL*8,  INTENT(IN   ) :: XKAPPA(NORBT,NORBT,NZ)

      ! Internal variables
      REAL*8, PARAMETER :: DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00
      REAL*8, PARAMETER :: D1P5 = 1.5D00, THREQL = 1.0D-12, DM1 = -D1

      CHARACTER SECTID*12,CPUTID*12,WALLTID*12

      INTEGER :: I, J, NORBXQ, IERR, IADDR, IOFF, I_MAX, J_MAX, LUPRI
      REAL*8  :: CPU1, WALL1, CPU2, WALL2, CPU3, WALL3, CPU4, WALL4, EIG, SQEIG, D_MAX
      REAL*8  :: ABS, SQRT, COS, SIN

      REAL*8, ALLOCATABLE, DIMENSION(:) :: REIG, RSIN, RCOS, RV, RTMP

  !------------------------------------------------------------------------------!
      LUPRI  = VAR%LUPRI

      CALL GETTIM(CPU1,WALL1)
      IF (PM%IPRLOC.GE.5) CALL HEADER('Output from ROTCMO_MOD',-1)

      NORBXQ = NORBI*NORBI*NZ
      CALL ALLOC(REIG,NORBI ,ID='REIG IN ROTCMO_MOD')
      CALL ALLOC(RSIN,NORBXQ,ID='RSIN IN ROTCMO_MOD')
      CALL ALLOC(RV  ,NORBXQ,ID='KV   IN ROTCMO_MOD')
      CALL ALLOC(RTMP,NORBXQ,ID='RTMP IN ROTCMO_MOD')

      CALL ALLOC(RCOS,NORBI*NBASI*NZ,ID='RCOS IN ROTCMO_MOD')

!     Calculate -\kappa^2 = \kappa^{\dagger} \kappa
!     (use RCOS as temporary array)
!     ---------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,D1,             &
                 'H','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D0,IP,RCOS,NORBI,NORBI,NZ)

      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: -kappa^2',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Diagonalize -\kappa^2
!     ---------------------------------------------

      CALL QDIAG90(NZ,NORBI,RCOS,NORBI,NORBI, &
                   REIG,1,RV,NORBI,NORBI,IERR)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V matrix (eigenvector)',-1)
         CALL PRQMAT(RV,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO_MOD: eigenvalues',-1)
         CALL PRQMAT(REIG,1,NORBI,1,NORBI,1,IP,LUPRI)
      END IF

!     Make cos(d) and sin(d)/d ; d = sqrt(eigenvalues)
!     ------------------------------------------------

      CALL GETTIM(CPU2,WALL2)
      RCOS(1:NORBXQ) = D0
      RSIN(1:NORBXQ) = D0
      DO I = 1,NORBI
         IADDR = (I-1) + (I-1)*NORBI
         EIG = REIG(I)
         IF (EIG .LT. D0) THEN
            IF (PM%IPRLOC.GE.3 .AND. ABS(EIG).GT.VAR%THRZER) WRITE(LUPRI,9010) EIG
            EIG = D0
            REIG(I)   = D0
         END IF
         IF (EIG .LT. DTHRS) THEN

!           Use Taylor expansion to avoid division by zero and/or
!           to speed up evaluation (avoid call to sin/cos/sqrt).
!           For EIG < 1.0D-8 the 3. order Taylor expansion is correct
!             to machine precision.

            RCOS(1+IADDR) = D1-DP5*EIG
            RSIN(1+IADDR) = D1-DSIXTH*EIG
         ELSE
            SQEIG = SQRT(EIG)
            RCOS(1+IADDR) = COS(SQEIG)
            RSIN(1+IADDR) = SIN(SQEIG)/SQEIG
         END IF
      END DO
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: sin(d)/d',-1)
         CALL PRQMAT(RSIN,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
         CALL HEADER('ROTCMO_MOD: cos(d)',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU3,WALL3)
 9010 FORMAT('WARNING from ROTCMO_MOD: Negative eigenvalue ',1P,D15.6)

!     Make V cos(d) V^{\dagger}
!     (use KTMP as temporatory storage)
!     ------------------------------------

      CALL QTRANS90('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI, &
                    RTMP,NORBI,NORBI,NZ,IP,                &
                    RCOS,NORBI,NORBI,NZ,IP,                &
                    RV,NORBI,NORBI,NZ,IP,                  &
                    RV,NORBI,NORBI,NZ,IP,PM%IPRLOC)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V cos(d) V^{dagger}',-1)
         CALL PRQMAT(RTMP,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V sin(d)/d V^{\dagger}
!     (use RCOS as temporatory storage)
!     ------------------------------------

      CALL QTRANS90('MOAO','S',D0,NORBI,NORBI,NORBI,NORBI, &
                    RCOS,NORBI,NORBI,NZ,IP,                &
                    RSIN,NORBI,NORBI,NZ,IP,                &
                    RV,NORBI,NORBI,NZ,IP,                  &
                    RV,NORBI,NORBI,NZ,IP,PM%IPRLOC)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
!     ---------------------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,DM1,            &
                 'N','N',IP,RCOS,NORBI,NORBI,NZ,   &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D1,IP,RTMP,NORBI,NORBI,NZ)
      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: exp(-kappa)',-1)
         CALL PRQMAT(RTMP,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Check if unitary

      CALL QGEMM(NORBI,NORBI,NORBI,D1,           &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ, &
                 'H','N',IP,RTMP,NORBI,NORBI,NZ, &
                 D0,IP,RCOS,NORBI,NORBI,NZ)
      IOFF = 0
      I_MAX = 0
      J_MAX = 0
      D_MAX = d0
      DO I = 1,NORBI
      DO J = 1,NORBI
         IF (I.EQ.J) THEN
            IF (ABS(RCOS(1+IOFF)-D1).GT.D_MAX) THEN
               I_MAX = I
               J_MAX = J
               D_MAX = ABS(RCOS(1+IOFF)-D1)
            END IF
         ELSE
            IF (RCOS(1+IOFF).GT.D_MAX) THEN
               I_MAX = I
               J_MAX = J
               D_MAX = ABS(RCOS(1+IOFF))
            END IF
         END IF
         IOFF = IOFF + 1
      END DO
      END DO
      IF(PM%IPRLOC.GE.1)THEN
        IF(D_MAX.GT.VAR%THRZER)THEN
           WRITE(LUPRI,'(/A)')' WARNING: ROTCMO_MOD: maximum deviation from unitary matrix'
           WRITE(LUPRI,*)'WARNING: Element ',J_MAX,I_MAX,' value :',D_MAX
        ENDIF
      ENDIF

!     Calculate new orbitals: SCOEF_{old} exp(-kappa)
!     ---------------------------------------------

      IF (PM%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: old orbitals',-1)
         CALL PRQMAT(SCOEF,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
!SK   CALL DZERO(RCOS,NORBXQ)
      CALL QGEMM(NBASI,NORBI,NORBI,D1,             &
                 'N','N',IP,SCOEF,NBASI,NORBI,NZ,  &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ,   &
                 D0,IP,RCOS,NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,RCOS,1,SCOEF,1)
      IF (PM%IPRLOC.GE.5)THEN
         CALL HEADER('ROTCMO_MOD: new orbitals',-1)
         CALL PRQMAT(SCOEF,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (PM%IPRLOC.GE.2)THEN
         CPUTID = SECTID(CPU4-CPU1)
         WALLTID = SECTID(WALL4-WALL1)
         WRITE(LUPRI,9020) CPUTID,WALLTID
         CPUTID = SECTID(CPU3-CPU2)
         WALLTID = SECTID(WALL3-WALL2)
         WRITE(LUPRI,9021) CPUTID,WALLTID
      END IF
 9020 FORMAT(/'ROTCMO_MOD: total CPU (WALL) time       : ',A12,' (',A12,')')
 9021 FORMAT( 'ROTCMO_MOD: sin/sqrt/cos CPU (WALL) time: ',A12,' (',A12,')')

      CALL DEALLOC(REIG,ID='REIG IN ROTCMO_MOD')
      CALL DEALLOC(RSIN,ID='RSIN IN ROTCMO_MOD')
      CALL DEALLOC(RCOS,ID='RCOS IN ROTCMO_MOD')
      CALL DEALLOC(RV  ,ID='KV   IN ROTCMO_MOD')
      CALL DEALLOC(RTMP,ID='RTMP IN ROTCMO_MOD')

  END SUBROUTINE ROTCMO_MOD

end module func_PipekMezey
