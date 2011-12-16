module localization

use memory_allocator
use func_PipekMezey
implicit none

public  :: orbloc, initialize_localization

private :: ORBLP2, ORBLO4, ORBLO5, ROTCMO_MOD
private :: get_kappa, hdiag
private :: get_mu, updateT

type :: localization_flags
  integer :: iprloc, itrloc
  real*8  :: thfull, thdiag, thgrad
  logical :: prjloc, lgfull, lgdiag, lggrad, lgchck
  character*4 :: hesloc
end type localization_flags

type :: external_variales
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
end type external_variales

type( localization_flags ), private :: loc
type( external_variales  ), private :: var

! Some global parameters
real*8, parameter, private :: thresh = 1.0d-13, trsh_mu = 1.0d-4
real*8, parameter, private :: d0 = 0.0d0, d1 = 1.0d0, d2 = 2.0d0, d4 = 4.0d0

contains

  subroutine initialize_localization
  !==============================================================================!
  !               Initialize localization for Dirac program                      !
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

  end subroutine initialize_localization

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine orbloc
  !==============================================================================!
  !                                                                              !
  !     Molecular orbital localization                                           !
  !                                                                              !
  !==============================================================================!
  implicit none

  ! Internal variables
  character*12 :: sectid, cputid, walltid
  character    :: daytid*24, hessian*4, chstart*23
  logical      :: tobe, full, rejected, Iminmin, trdef, flat, nstep, restore

  integer :: i, ifrp, ndimh, idum, itr, ist, lupri, nish, norb, nucdep, maxref
  integer :: nfbas(0:2), ipqtoq(4,0:7), nep, nem, nez, nref, ndim

  real*8  :: dinvmo, dinvmg, dinvmq, dnrmk, dnrmg, diffg, diffq, diffa
  real*8  :: timend, time, timstr, second, toterg, totdiag, totnond, shift, tradius, tpred

  integer, allocatable, dimension(:)       :: nfnuc, isel, indao, indaos, indnao
  real*8,  allocatable, dimension(:)       :: scoef, storec, smat, cmo, grad, etmp
  real*8,  allocatable, dimension(:,:,:)   :: xkappa
  real*8,  allocatable, dimension(:,:,:,:) :: pmat

  integer, pointer :: iorig(:)
  real*8,  pointer :: pkappa(:), tgrad(:), eval(:), hmat(:,:)

  !------------------------------------------------------------------------------!
  call initialize_func

  ! Set some parameters
  ifrp   = 1
  lupri  = var%lupri
  nish   = var%nish(ifrp)
  norb   = var%norb(ifrp)
  nucdep = var%nucdep
  maxref = var%maxref

  nfbas(:)    = var%nfbas(ifrp,:)
  ipqtoq(:,:) = var%ipqtoq(:,:)

  ! Check thresholds
  if(loc%thdiag < thresh)then
    write(6,'(/a)')      ' WARNING: The threshold (THDIAG) from the input is too small '
    write(6,'(a,d8.1)')  ' WARNING: It will be changed to the minimal reasonable value: ',thresh
    loc%thdiag = thresh
  endif
  if(loc%thfull < thresh)then
    write(6,'(/a)')      ' WARNING: The threshold (THFULL) from the input is too small '
    write(6,'(a,d8.1)')  ' WARNING: It will be changed to the minimal reasonable value: ',thresh
    loc%thfull = thresh
  endif

  ! Set defaults
  full = .false.
  if(loc%hesloc == 'FULL') full = .true.
  hessian = 'DIAG'
  if(full) hessian = 'FULL'

  rejected = .false.
  Iminmin  = .false.
  chstart = 'START ITERATION NO.    '

  tradius = 10.0d0
  trdef   = .true.
  flat    = .false.
  nstep   = .true.
  restore = .false.

  if(loc%iprloc.ge.4)then
    write(lupri,*) 'lupri =',lupri
    write(lupri,*) 'nish  =',nish
    write(lupri,*) 'norb  =',norb
    write(lupri,*) 'nucdep=',nucdep
    write(lupri,*) 'maxref=',maxref
    write(lupri,*) 'full  =',full
    write(lupri,*) 'thdiag=',loc%thdiag
    write(lupri,*) 'thfull=',loc%thfull
    write(lupri,*) 'thgrad=',loc%thgrad
    write(lupri,*) 'lgdiag=',loc%lgdiag
    write(lupri,*) 'lgfull=',loc%lgfull
    write(lupri,*) 'lggrad=',loc%lggrad
    write(lupri,*) 'lgchck=',loc%lgchck
    write(lupri,*) 'hesloc=',loc%hesloc
    write(lupri,*) 'thrzer=',var%thrzer
  endif

  call titler('Molecular orbital localization','*',122)

  ! Check if coefficients are on file
  inquire(file='DFCOEF',exist=tobe)
  if(.not.tobe)then
    write(lupri,'(a)') 'orbloc: Coefficient file not found !'
    call quit('orbloc: Coefficients not found !')
  endif

  ! Memory allocation
  ndimh = 4 * ( nish*nish - nish ) / 2
  call alloc(scoef ,nish*nfbas(0)*4,   id='scoef  in orbloc')
  call alloc(storec,nish*nfbas(0)*4,   id='storec in orbloc')
  call alloc(smat  ,var%n2bbasx,       id='smat   in orbloc')
  call alloc(cmo   ,var%n2bbasxq,      id='cmo    in orbloc')
  call alloc(grad  ,ndimh,             id='grad   in orbloc')
  call alloc(etmp  ,ndimh,             id='etmp   in orbloc')
  call alloc(xkappa,nish,nish,4,       id='xkappa in orbloc')
  call alloc(pmat  ,nish,nish,4,nucdep,id='pmat   in orbloc')

  call alloc(nfnuc ,2*nucdep,  id='nfnuc  in orbloc')
  call alloc(isel  ,nish,      id='isel   in orbloc')
  call alloc(indao ,nfbas(0),  id='indao  in orbloc')
  call alloc(indaos,nfbas(0),  id='indaos in orbloc')
  call alloc(indnao,2*nucdep+2,id='indnao in orbloc')

  ! nullify pointers for later allocation using memory_allocator module
  nullify(iorig )
  nullify(tgrad )
  nullify(eval  )
  nullify(pkappa)
  nullify(hmat  )

  ! Read coefficients
  call opnfil(var%lucoef,'DFCOEF','OLD','PAMANA')
  call reacmo_no_work(var%lucoef,'DFCOEF',cmo,idum,idum,toterg,2)
  close(var%lucoef,status='keep')

  ! Get overlap matrix 
  call gtovlx(smat,var%ssmtrc)
  if(loc%iprloc.ge.6)then
    write(lupri,'(/a)') '* orbloc: Molecular coefficients:'
    call prqmat(cmo,nfbas(0),norb,nfbas(0),norb,4,ipqtoq(1,0),lupri)
    write(lupri,'(/a)') '* orbloc: Overlap matrix' 
    call prqmat(smat,nfbas(0),nfbas(0),nfbas(0),nfbas(0),1,ipqtoq(1,0),lupri)
  endif

  ! Occupied electronic orbital coeffients selection
  call occupy(isel,nish)
  call selcfs(cmo,ifrp,scoef,nish,isel,d0,nish,nfbas(0),norb)
  if(loc%iprloc.ge.5)then
    write(lupri,'(/a)') '* orblo1: Selected molecular coefficients'
    call prqmat(scoef,nfbas(0),nish,nfbas(0),nish,4,ipqtoq(1,0),lupri)
  endif

  ! Get the number of functions per center
  call indaoc(indao,indaos,indnao)
  call orblo4(indnao,nfnuc,nucdep)

  !------------------------------------------------------------------------------!
  !                           Main iteration procedure                           !
  !------------------------------------------------------------------------------!
  call titler(' REFERENCE POINT ','=',122)

  timstr = second()
  call get_func(smat,scoef,pmat,nfnuc,dinvmg)
  dinvmo = dinvmg

  write(lupri,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 4x,a, 3x,a, 4x,a, 5x,a, 6x,a)') &
    '*T','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'
  write(lupri,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 2x,a, 3x,a, 5x,a, 6x,a, 6x,a)') &
    '*C','ITR','DIFFA','|GRAD|','|KAPPA|','H(nonD/D)','e: +,-,0','H','DIFFG','DIFFQ'

  ! zero step
  call get_hessian(pmat,hmat,full,totdiag,totnond)
  call get_gradient(pmat,grad,dnrmg)
  call hdiag(hmat,grad,etmp,tgrad,eval,iorig,nep,nem,nez,Iminmin,full)

  shift  = d0
  dnrmk  = d0
  dinvmq = d0; diffq = d0
  dinvmg = d0; diffg = d0
  itr    = 0
  ist    = 0
  write(lupri,'(a,2x,i4,3d11.2,f10.5,1x,3i4,3x,a,2d11.2)') &
    '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq

  do
    ist = ist + 1

    ! start i-th iteration
    if( .not.rejected .and. .not.restore )then
      write(chstart(20:23),'(i4)') itr+1
      call titler(chstart,'=',122)
      call gtinfo(daytid)
      write(lupri,'(/6x,a,a24)')'--- ',daytid
    endif

    ! get new kappa
    shift = d0
    if(rejected) call get_mu(eval,tgrad,tradius,nref,shift,nstep,tpred)
    call get_kappa(xkappa,pkappa,tgrad,eval,hmat,iorig,nstep,full,tradius,shift,dnrmk)
    call exp_func(pkappa,tgrad,eval,nish,dinvmq)

    if( (rejected) .and. (dabs(tpred-dnrmk) > trsh_mu*tradius) )then
      write(lupri,'(a)') ' WARNING: Predicted and actual norm of kappa are too different'
      write(lupri,'(a,i4,a,d12.4,a,d12.4,a,2d12.4)')&
        ' WARNING: itr:',itr,' Predicted:',tpred,' Actual:',dnrmk,' Trust-radius:',tradius
    endif

    ! get new MO coefficients
    if(restore)then
      scoef(:) = storec(:)  ! restore MO coefficients
    elseif(rejected)then
      scoef(:) = storec(:)  ! restore MO coefficients; the step was rejected
    elseif(.not.rejected)then
      storec(:) = scoef(:)  ! store MO coefficients
    endif
    call rotcmo_mod(scoef,nfbas(0),nish,4,xkappa,nish,ipqtoq(ifrp,0))

    ! get functional value in the new point (new MO)
    call get_func(smat,scoef,pmat,nfnuc,dinvmg)

    ! Update Trust-radius
    diffg = dinvmg - dinvmo
    diffq = dinvmq
    call updateT(dnrmk,diffg,diffq,tradius,rejected,trdef,flat,nstep,restore,full,itr,nem,&
                 pkappa,tgrad,eval,nref,dnrmg)

    ! Move to the new point if the step was not rejected or old step is going to be restored
    if( .not.rejected .and. .not.restore )then
      !---------------------------------------------!
      ! Calculate hessian and gradient in new point !
      !---------------------------------------------!
      call get_hessian(pmat,hmat,full,totdiag,totnond)
      call get_gradient(pmat,grad,dnrmg)
      call hdiag(hmat,grad,etmp,tgrad,eval,iorig,nep,nem,nez,Iminmin,full)

      itr = itr + 1

      if((diffg > d0).and.(diffg > thresh)) &
        write(lupri,'(a,i4)')' WARNING: dinvm_new - dinvm_old > 0; itr =',itr
      diffa = dabs(dinvmg - dinvmo)
      hessian = 'DIAG'
      if(full) hessian = 'FULL'
      dinvmo = dinvmg

      write(lupri,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 2x,a, 3x,a, 5x,a, 6x,a, 6x,a)') &
        '  ','ITR','DIFFA','|GRAD|','|KAPPA|','H(nonD/D)','e: +,-,0','H','DIFFG','DIFFQ'
      write(lupri,'(a,2x,i4,3d11.2,f10.5,1x,3i4,3x,a,2d11.2)') &
        '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq

      !----------------------------!
      ! check convergence criteria !
      !----------------------------!
      if(loc%hesloc == 'DIAG')then
        if( .not.loc%lggrad .and.   .not.loc%lgdiag    .and. Iminmin ) exit
        if(      loc%lggrad .and. (dnrmg < loc%thgrad) .and. Iminmin ) exit
        if(      loc%lgdiag .and. (diffa < loc%thdiag) .and. Iminmin ) exit
      endif

      if(loc%hesloc == 'FULL')then
        if( .not.loc%lggrad .and. .not.loc%lgfull      .and. Iminmin ) exit
        if(      loc%lggrad .and. (dnrmg < loc%thgrad) .and. Iminmin ) exit
        if(      loc%lgfull .and. (diffa < loc%thfull) .and. Iminmin ) exit
      endif

      if(loc%hesloc == 'COMB')then
        if(full)then
          if( .not.loc%lggrad .and. .not.loc%lgfull      .and. Iminmin ) exit
          if(      loc%lggrad .and. (dnrmg < loc%thgrad) .and. Iminmin ) exit
          if(      loc%lgfull .and. (diffa < loc%thfull) .and. Iminmin ) exit
        else
          if( .not.loc%lgdiag                            .and. Iminmin ) full = .true.
          if(      loc%lgdiag .and. (dnrmg < loc%thdiag) .and. Iminmin ) full = .true.
          if(full)then
            call get_hessian(pmat,hmat,full,totdiag,totnond)
            call get_gradient(pmat,grad,dnrmg)
            call hdiag(hmat,grad,etmp,tgrad,eval,iorig,nep,nem,nez,Iminmin,full)
            hessian = 'FULL'
            write(lupri,'(a,2x,i4,3d11.2,f10.5,1x,3i4,3x,a,2d11.2)') &
              '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq
            if(loc%lgchck) exit
          endif
        endif
      endif

      if(itr > loc%itrloc)exit
    endif

    if(ist > 2*loc%itrloc)exit

  enddo

  write(lupri,'(//11x,a)')      '************************************************************'
  if(Iminmin)then
    write(lupri,'(11x,a,i4,a,i4,a)') &
         '**** MAXIMUM reached after',itr,' iterations and',ist,' steps ****'
  elseif( full .and. ((diffa < loc%thfull).and.(dnrmg < loc%thgrad)) )then
    write(lupri,'(11x,a,i4,a,i4,a)') &
         '****    Convergence in',itr,' iterations and',ist,' steps     ****'
  else
    write(lupri,'(11x,a,i4,a,i4,a)') &
         '**** No convergence after',itr,' iterations and',ist,' steps  ****'
  endif
  write(lupri,'(11x,a,f21.15,a)')'**** Inverse mean delocalization =',dinvmo,' ****'
  write(lupri,'(11x,a//)')      '************************************************************'

  timend = second()
  time   = timend - timstr
  call timtxt('          --- Time used in localization is :',time,lupri)

  ! Put the localised coefficients in DFCOEF
  call orblo5(scoef,nfbas(0),nish,norb)

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
  call dealloc(scoef ,id='scoef  in orbloc')
  call dealloc(storec,id='storec in orbloc')
  call dealloc(smat  ,id='smat   in orbloc')
  call dealloc(pmat  ,id='pmat   in orbloc')
  call dealloc(xkappa,id='xkappa in orbloc')
  call dealloc(grad  ,id='grad   in orbloc')
  call dealloc(etmp  ,id='etmp   in orbloc')
  call dealloc(cmo   ,id='cmo    in orbloc')

  call dealloc(nfnuc ,id='nfnuc  in orbloc')
  call dealloc(isel  ,id='isel   in orbloc')
  call dealloc(indao ,id='indao  in orbloc')
  call dealloc(indaos,id='indaos in orbloc')
  call dealloc(indnao,id='indnao in orbloc')

  call deallocp(iorig ,id='iorig  in orbloc')
  call deallocp(pkappa,id='pkappa in orbloc')
  call deallocp(tgrad ,id='tgrad  in orbloc')
  call deallocp(eval  ,id='eval   in orbloc')
  call deallocp(hmat  ,id='hmat   in orbloc')

  end subroutine orbloc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_kappa(xkappa,pkappa,tgrad,eval,hmat,iorig,nstep,full,tradius,shift,dnrmk)
  !******************************************************************************!
  !                                                                              !
  !           Construction of the XKAPPA matrix using full hessian               !
  !          The linear set of equations is solved using information             !
  !                   from diagonalization of the hessian                        !
  !                                                                              !
  !   (H-shift) * X = -g     -->     (tilde{H}-shift) * tilde{X} = - tilde{g}    !
  !                                                                              !
  !   tilde{H} = U^T * H * U = eval                 U = hmat (iorig)             !
  !   tilde{X} = U^T * X     = pkappa               X = xkappa                   !
  !   tilde{g} = U^T * g     = tgrad                dnrmk = |X|                  !
  !                                                                              !
  ! Notes: 1. Deffinition of kappa matrix in rotcmo_mod is different as in       !
  !           calculation of gradient and hessian in get_gradient and get_hessian!
  !           --> Indexes are exchanged and j quaterion is with minus sign       !
  !               Qkappa_lk = kappa_kl - bar{kappa}_kl^* j                       !
  !               It follows from definition of bared unitary matrix:            !
  !               tilde{Psi_i} = U_ij Psi_j + bar{U}_ij Psi_j                    !
  !        2. Be also careful about minus in imaginary part of gradient in       !
  !                                H * K = - g                                   !
  !                                 H                     K            g         !
  !  / A   B  \   -->   / A^R + B^I   B^I - A^I \      / K^R \      / g^R \      !
  !  \ B*  A* /   -->   \ B^I + A^I   A^R - B^I /      \ K^I /      \-g^I /      !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in   ) :: full
  logical, intent(inout) :: nstep
  real*8,  intent(in ) :: shift, tradius
  real*8,  intent(out) :: dnrmk

  real*8, intent(in ), dimension(:,:)   :: hmat
  real*8, intent(out), dimension(:,:,:) :: xkappa

  integer, intent(in ), pointer :: iorig(:)
  real*8,  intent(in ), pointer :: tgrad(:), eval(:)
  real*8,  intent(out), pointer :: pkappa(:)

  ! Internal variables
  integer :: i, j, ij1, ij2, ij3, ij4, iz, ni, ndimh, ni12, nred, lupri
  real*8  :: dnrm2, dtmp

  real*8, allocatable, dimension(:) :: xtmp

  !------------------------------------------------------------------------------!
  lupri  = var%lupri

  ! Get some dimensions
  ni    = size(xkappa,1)
  ndimh = size(hmat,1)
  ni12  = ndimh / 4
  nred  = size(eval,1)

  ! allocate arrays
  if(associated(pkappa)) call deallocp(pkappa,id='pkappa in get_kappa')
  call allocp(pkappa,nred ,id='pkappa in get_kappa')
  call alloc(xtmp   ,ndimh,id='xtmp   in get_kappa')

  ! calculate tilde{X}
  if(nstep)then
    do i=1,nred
      pkappa(i) = - tgrad(i) / (eval(i) - shift)
    enddo
  else
    dtmp = d0
    do i=1,nred
      if(eval(i) > d0) exit
      dtmp = dtmp + eval(i)*eval(i)
    enddo
    dtmp = dsqrt(dtmp)  ! norm of the negative eigenvalues
    ! make step in the direction of negative eigenvalues
    ! eigenvalues are weighted with their relative size 
    pkappa(:) = d0
    do i=1,nred
      if(eval(i) > d0) exit
      pkappa(i) = tradius * dabs(eval(i)) / dtmp
    enddo
  endif

  ! calculate X
  xtmp(:) = d0
  do j=1,nred
    if(full)then
      do i=1,ndimh
        xtmp(i) = xtmp(i) + hmat(i,j)*pkappa(j)
      enddo
    else
      xtmp(iorig(j)) = pkappa(j)
    endif
  enddo

  ! Norm of the X
  dnrmk = dnrm2(ndimh,xtmp,1)
  if(loc%iprloc.ge.1)then
    write (lupri,'(/a,e12.2)') ' Norm of the kappa: ',dnrmk
  endif

  ! Unpack X to kappa matrix
  xkappa(:,:,:) = d0
  ij1 = 0
  do j = 1,ni
    do i = j+1,ni
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

  if(loc%iprloc.ge.2)then
    do iz=1,4
      write(lupri,'(/a,i4)') 'Print xkappa matrix, IZ =',iz
      do i=1,ni
        write(lupri,'(i4,10es9.1)') i,(xkappa(i,j,iz),j=1,ni)
      enddo
    enddo
  endif

  call dealloc(xtmp,id='xtmp in get_kappa')

  end subroutine get_kappa

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLO4(INDNAO,NFNUC,NUCDEP)
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

      ! External variables
      INTEGER, INTENT(IN ) :: NUCDEP

      INTEGER, INTENT(IN ) :: INDNAO(2,*)
      INTEGER, INTENT(OUT) :: NFNUC(2*NUCDEP)

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

  SUBROUTINE ORBLO5(SCOEF,NB,NI,NO)
!***********************************************************************
!
!     Written by S.Dubillard Apr 13 2004
!     Last revision:
!
!     Put the localised coefficients into the DFCOEF file (1/2)
!
!***********************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER :: NB, NI, NO
      REAL*8  :: SCOEF(NB,NI,4)

      ! Internal variables
      INTEGER :: I, J, K, L
      REAL*8  :: TOTERG

      INTEGER, ALLOCATABLE :: IBEIG(:)
      REAL*8,  ALLOCATABLE :: CMO(:,:,:), EIG(:)

  !------------------------------------------------------------------------------!
      ! allocate arrays
      CALL ALLOC(CMO  ,NB,NO,4,ID='CMO   IN ORBLO5')
      CALL ALLOC(EIG  ,VAR%NORBT  ,ID='EIG   IN ORBLO5')
      CALL ALLOC(IBEIG,VAR%NORBT  ,ID='IBEIG IN ORBLO5')

      ! open DFCOEF file and read EIG, IBEIG and TOTERG
      CALL OPNFIL(VAR%LUCOEF,'DFCOEF','OLD','PAMANA')
      CALL REACMO_no_work(VAR%LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)

      ! reorganize MO
      DO I = 1,NB
        DO J = 1,NI
          DO K = 1,4
            L = J + NO/2
            IF (VAR%LEVYLE) L = J
            CMO(I,L,K) = SCOEF(I,J,K)
          ENDDO
        ENDDO
      ENDDO
      IF (LOC%IPRLOC.GE.5)THEN
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

  SUBROUTINE ROTCMO_MOD(CMO,NBASI,NORBI,NZ,XKAPPA,NORBT,IP)
!***********************************************************************
!
!     Rotate coefficients
!
!     CMO_{new} = CMO_{old} exp( -\kappa )
!     We use a trick a la Sostrup notes, chapter 3.1.5
!
!     Input:
!        CMO     - old coefficients
!        XKAPPA  - Kappa matrix (Kappa is overwritten with junk)
!        IPRINT  - print level
!
!     Output:
!        CMO     - new coefficients
!
!     Written by J. Thyssen - Nov 24 1998
!
!***********************************************************************
      IMPLICIT NONE

      ! External variables
      INTEGER, INTENT(IN   ) :: NBASI, NORBI, NZ, NORBT, IP(4)
      REAL*8,  INTENT(IN   ) :: XKAPPA(NORBT,NORBT,NZ)
      REAL*8,  INTENT(INOUT) :: CMO(NBASI,NORBI,NZ)

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
      IF (LOC%IPRLOC.GE.5) CALL HEADER('Output from ROTCMO_MOD',-1)

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

      IF (LOC%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: -kappa^2',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Diagonalize -\kappa^2
!     ---------------------------------------------

      CALL QDIAG90(NZ,NORBI,RCOS,NORBI,NORBI, &
                   REIG,1,RV,NORBI,NORBI,IERR)
      IF (LOC%IPRLOC.GE.10)THEN
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
            IF (LOC%IPRLOC.GE.3 .AND. ABS(EIG).GT.VAR%THRZER) WRITE(LUPRI,9010) EIG
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
      IF (LOC%IPRLOC.GE.10)THEN
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
                    RV,NORBI,NORBI,NZ,IP,LOC%IPRLOC)
      IF (LOC%IPRLOC.GE.10)THEN
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
                    RV,NORBI,NORBI,NZ,IP,LOC%IPRLOC)
      IF (LOC%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: V sin(d)/d V^{dagger}',-1)
         CALL PRQMAT(RCOS,NORBI,NORBI,NORBI,NORBI,NZ,IP,LUPRI)
      END IF

!     Make V cos(d) V^{\dagger} - V sin(d)/d V^{\dagger} \kappa
!     ---------------------------------------------------------

      CALL QGEMM(NORBI,NORBI,NORBI,DM1,            &
                 'N','N',IP,RCOS,NORBI,NORBI,NZ,   &
                 'N','N',IP,XKAPPA,NORBT,NORBT,NZ, &
                 D1,IP,RTMP,NORBI,NORBI,NZ)
      IF (LOC%IPRLOC.GE.10)THEN
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
      IF(LOC%IPRLOC.GE.1)THEN
        IF(D_MAX.GT.VAR%THRZER)THEN
           WRITE(LUPRI,'(/A)')' WARNING: ROTCMO_MOD: maximum deviation from unitary matrix'
           WRITE(LUPRI,*)'WARNING: Element ',J_MAX,I_MAX,' value :',D_MAX
        ENDIF
      ENDIF

!     Calculate new orbitals: CMO_{old} exp(-kappa)
!     ---------------------------------------------

      IF (LOC%IPRLOC.GE.10)THEN
         CALL HEADER('ROTCMO_MOD: old orbitals',-1)
         CALL PRQMAT(CMO,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
!SK   CALL DZERO(RCOS,NORBXQ)
      CALL QGEMM(NBASI,NORBI,NORBI,D1,           &
                 'N','N',IP,CMO,NBASI,NORBI,NZ,  &
                 'N','N',IP,RTMP,NORBI,NORBI,NZ, &
                 D0,IP,RCOS,NBASI,NORBI,NZ)
      CALL DCOPY(NBASI*NORBI*NZ,RCOS,1,CMO,1)
      IF (LOC%IPRLOC.GE.5)THEN
         CALL HEADER('ROTCMO_MOD: new orbitals',-1)
         CALL PRQMAT(cmo,NBASI,NORBI,NBASI,NORBI,NZ,IP,LUPRI)
      END IF
      CALL GETTIM(CPU4,WALL4)
      IF (LOC%IPRLOC.GE.2)THEN
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine hdiag(hmat,gradient,etmp,tgrad,eval,iorig,nep,nem,nez,Iminmin,full)
  !******************************************************************************!
  !                                                                              !
  !                Diagonalize real representation of hessian                    !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: full
  logical, intent(out) :: Iminmin
  integer, intent(out) :: nep, nem, nez
  real*8,  intent(in ), dimension(:)   :: gradient
  real*8,  intent(out), dimension(:)   :: etmp

  integer, intent(  out), pointer :: iorig(:)
  real*8,  intent(inout), pointer :: hmat(:,:)
  real*8,  intent(  out), pointer :: tgrad(:), eval(:)

  ! Internal variables
  integer :: i, j, k, ndimh, nred, ier, lupri, lwork, ilaenv, nb, icur
  real*8  :: acur

  real*8,  allocatable, dimension(:)   :: work

  !------------------------------------------------------------------------------!
  lupri = var%lupri

  ! Get some dimensions
  ndimh = size(hmat,1)

  if(full)then

    nb = ilaenv(1,'dsytrd','u',ndimh,-1,-1,-1)
    lwork = (nb+2)*ndimh  ! see the comment in dsyev lapack routine
  
    call alloc(work,lwork,id='work in hdiag')
  
    ! get eigenvalues (etmp) and eigenvectors (hmat)
    call dsyev('V','U',ndimh,hmat,ndimh,etmp,work,lwork,ier)
    if(ier /= 0)then
      write(lupri,*) ' ERROR: The dsyev LAPACK subroutine in hdiag finished unsuccessfully'
      write(lupri,*) ' info =',ier
      stop
    endif

    call dealloc(work,id='work in hdiag')

  else
  
    if(associated(iorig)) call deallocp(iorig,id='iorig in hdiag')
    call allocp(iorig,ndimh,id='tgrad in hdiag')

    ! get eigenvalues (etmp)
    ! We use diagonal approximation so it will just sort the diagonal of the hessian
    ! iorig will contain information how the eigenvalues were sorted
    do i=1,ndimh
      etmp(i) = hmat(i,1)
      iorig(i)= i
    enddo
    do i=1,ndimh
      acur = etmp(i)
      icur = i
      do j=i+1,ndimh
        if(acur > etmp(j))then
          acur = etmp(j)
          icur = j
        endif
      enddo
      etmp(icur)  = etmp(i)
      etmp(i)     = acur
      j = iorig(icur)
      iorig(icur) = iorig(i)
      iorig(i)    = j
    enddo

  endif

  ! allocate unnamed arrays (dimension is number of eigenvalues over some threshold)
  nred = 0
  do i=1,ndimh
    if(dabs(etmp(i)) > thresh) nred = nred + 1
  enddo
  if(associated(eval))then
    call deallocp(tgrad,id='tgrad in hdiag')
    call deallocp(eval ,id='eval in hdiag')
  endif
  call allocp(tgrad,nred,id='tgrad in hdiag')
  call allocp(eval ,nred,id='eval in hdiag')

  ! project out redundant eigevalues and corresponding eigenvectors
  k = 0
  do i=1,ndimh
    if(dabs(etmp(i)) > thresh)then
      k = k + 1
      eval(k) = etmp(i)
      if(full)then
        do j=1,ndimh
          hmat(j,k) = hmat(j,i)
        enddo
      else
        iorig(k) = iorig(i)
      endif
    endif
  enddo

  ! calculate transformed gradient: tilde{g} = U^T * g
  tgrad(:) = d0
  do i=1,nred
    if(full)then
      do j=1,ndimh
        tgrad(i) = tgrad(i) + hmat(j,i) * gradient(j)
      enddo
    else
      tgrad(i) = gradient(iorig(i))
    endif
  enddo

  ! count number of positive, negative and zero eigenvalues
  nep=0; nem=0; nez=0
  do i=1,ndimh
    if(dabs(etmp(i)) > thresh)then
      if(etmp(i) < d0)nem = nem + 1
      if(etmp(i) > d0)nep = nep + 1
    else
      nez = nez + 1
    endif
  enddo
  if(nem == 0)Iminmin = .true.  ! we have reached the minimum of the function
  if(nem /= 0)Iminmin = .false.

  if(loc%iprloc.ge.1)then
    write(lupri,'(/a,3i4)')' Eigenvalues of the hessian;  +,-,0  =  ',nep,nem,nez
    write(lupri,'(8d12.4)') (etmp(i),i=1,ndimh)
  endif

  if(loc%iprloc.ge.4)then
    write(lupri,'(/a,3i4)')' Reduced eigenvalues of the hessian:'
    write(lupri,'(8d12.4)') (eval(i),i=1,nred)
    write(lupri,'(/a,3i4)')' Transformed gradient:'
    write(lupri,'(8d12.4)') (tgrad(i),i=1,nred)
  endif

  end subroutine hdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_mu(eval,tgrad,tradius,nref,shift,nstep,tpred)
  !******************************************************************************!
  !                                                                              !
  !                   Calculate level-shift parameter mu                         !
  !                                                                              !
  !          tradius = sqrt{ sum{ (tgrad)^2_i / (eval_i-shift)^2  }  }           !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: nstep
  integer, intent(in ) :: nref
  real*8,  intent(in ) :: tradius
  real*8,  intent(out) :: shift, tpred

  real*8, intent(in), pointer :: eval(:), tgrad(:)

  ! Internal variables
  integer :: i, j, k, ndim, lupri, imax
  real*8  :: tcurr, xcurr, cutoff, xmax, xmin, step

  !------------------------------------------------------------------------------!
  lupri = var%lupri

  ! If we are not performing Newton step quit.
  if(.not.nstep)then
    tpred = tradius
    return
  endif

  ! Get some dimensions
  ndim = size(eval,1)

  ! initialization
  imax  = 10000
  cutoff= trsh_mu * tradius
  step  = 0.1d0
  xmax  = eval(nref)
  xmin  = eval(nref) - step
  tcurr = f(xmin)

  ! find lower boundary if necessary
  i = 0
  if(tcurr > tradius)then
    do
      i = i + 1
      if(i > imax)then
        write(lupri,'(a)')' WARNING: Maximum number of iteration reached in get_mu(1)!!!'
        write(lupri,'(10x,a,d12.4,a,d12.4)')'f(xmin) =',tcurr,' tradius =',tradius
        exit
      endif
      step = xmax - xmin
      xmax = xmin
      xmin = xmin - d2*step
      tcurr = f(xmin)
      if(tcurr < tradius)exit
    enddo
  endif 

  i = 0
  do
  i = i + 1
  if(i > imax)then
    write(lupri,'(a)')' WARNING: Maximum number of iteration reached in get_mu(2)!!!'
    write(lupri,'(10x,a,d12.4,a,d12.4)')&
      'h(desired) =',tradius,' |h(desired)-h(current)| =',dabs(tradius - tcurr)
    exit
  endif

  if(dabs(tcurr-tradius) < cutoff) exit

  xcurr = (xmax+xmin)/d2
  tcurr = f(xcurr)

  if(tcurr > tradius) xmax = xcurr
  if(tcurr < tradius) xmin = xcurr

  enddo

  shift = xcurr
  tpred = f(xcurr)

  !------------------------------------------------------------------------------!

  contains

    function f(x)
    implicit none
    real*8, intent(in ) :: x
    real*8  :: f
    integer :: ii

    f = d0
    do ii=1,ndim
      if(dabs(tgrad(ii)) > thresh)then
        f = f + tgrad(ii) * tgrad(ii) / ((eval(ii) - x) * (eval(ii) - x))
      endif
    enddo
    f = dsqrt(f)

    end function f

  end subroutine get_mu

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine updateT(dnrmk,diffg,diffq,tradius,rejected,trdef,flat,nstep,restore,full,itr,nem,&
                     pkappa,tgrad,eval,nref,dnrmg)
  !******************************************************************************!
  !                                                                              !
  !                           Update trust radius                                !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  integer, intent(in ) :: itr, nem
  integer, intent(out) :: nref
  real*8,  intent(in   ) :: dnrmk, diffg, diffq, dnrmg
  real*8,  intent(inout) :: tradius
  logical, intent(in   ) :: full
  logical, intent(inout) :: rejected, trdef, flat, nstep, restore

  real*8,  intent(in), pointer :: pkappa(:), tgrad(:), eval(:)

  ! Internal variables
  character*6 :: cond
  integer :: i, lupri, nred
  integer, save :: nrun  = 0

  real*8            :: R, tradold, dmaxa, dmaxn
  real*8, save      :: diffg1, diffg2
  real*8, parameter :: d1q = 0.25d0, d3q = 0.75d0, d5q = 1.25d0

  !------------------------------------------------------------------------------!
  !                              Initialization                                  !
  !------------------------------------------------------------------------------!
  lupri   = var%lupri
  nred    = size(eval,1)
  tradold = tradius
  R       = diffg / diffq

  !------------------------------------------------------------------------------!
  !                            Trus-radius update                                !
  !------------------------------------------------------------------------------!
  if( dabs(diffg)<thresh .or. dabs(diffq)<thresh )then
    cond = '0     '
    if(rejected)then
      if(diffq > d0)then
        cond = '0r q>0'
        tradius  = tradius / 10.0d0
      else
        cond = '0r    '
        rejected = .false.
      endif
    endif
  else
    if(rejected)then
      if(diffq > d0)then
        cond = 'r  q>0'
        tradius  = tradius / 10.0d0
      else
        cond = 'r     '
        if(R < d0 ) cond = 'r  R<0'
        if(R < d1q) cond = 'r  R<1'
        if(R > d3q) cond = 'r  R>3'
        rejected = .false.
        if(R < d0 ) rejected = .true.
        if(R < d1q) tradius  = d1q* tradius
        if(R > d3q) tradius  = d2 * tradius
      endif
    else
      if(dnrmk > tradius)then
        if(full)then
          if(diffq < d0)then
            cond = 'Tq    '
            if(R < d0 ) cond = 'Tq R<0'
            if(R < d1q) cond = 'Tq R<1'
            if(R > d5q) cond = 'Tq R>5'
            tradius = dnrmk
            if(R < d0 ) rejected = .true.
            if(R < d1q) tradius  = d1q*tradius
            if(R > d5q) tradius  = d2*tradius
          else 
            cond = 'TQ    '
            if(R > d0 ) cond = 'TQ R>0'
            if((R < d1q).or.(R > d1+d3q)) cond = 'TQ R>7'
            if((R > d3q).and.(R < d5q))   cond = 'TQ R<5'
            tradius = dnrmk
            if(R > d0 ) rejected = .true.
            if((R < d1q).or.(R > d1+d3q)) tradius = d1q*tradius
            if((R > d3q).and.(R < d5q))   tradius = d2*tradius
          endif
        else
          cond = 'T  dia'
          rejected = .true.
        endif
      else
        if(diffq < d0)then
          cond = 'tq     '
          if(R < d0 ) cond = 'tq R<0'
          if(R < d1q) cond = 'tq R<1'
          if(R > d5q) cond = 'tq R>5'
          tradius = tradold
          if(trdef) tradius = dnrmk
          if(R < d0 ) rejected = .true.
          if(R < d1q) tradius  = d1q * dnrmk
          if(R > d5q) tradius  = max(d2*dnrmk,tradold)
        else 
          cond = 'tQ    '
          if(R > d0 ) cond = 'tQ R>0'
          if((R < d1q).or.(R > d1+d3q)) cond = 'tQ R>7'
          if((R > d3q).and.(R < d5q))   cond = 'tQ R<5'
          tradius = tradold
          if(trdef) tradius = dnrmk
          if(R > d0 ) rejected = .true.
          if((R < d1q).or.(R > d1+d3q)) tradius = d1q * dnrmk
          if((R > d3q).and.(R < d5q))   tradius = max(d2*dnrmk,tradold)
        endif
      endif
    endif
  endif

  !------------------------------------------------------------------------------!
  !                       Are we close to the minimum?                           !
  !            If both dmaxa, dmaxn are small we are near minimum.               !
  !------------------------------------------------------------------------------!
  dmaxn = d0
  do i=1,nred
    if(eval(i) > d0) exit
    dmaxn = dmaxn + dabs(d2*eval(i)*dnrmk)
  enddo
  if(nem == 0)then
    dmaxn = d1
  else
    dmaxn = dmaxn / dfloat(nem)
  endif
  dmaxa = dnrmg*dnrmk
  if(loc%iprloc.ge.2)then
    write(lupri,'(a,i4,2d14.6)') ' itr,dmaxa,dmaxn: ',itr,dmaxa,dmaxn
  endif

  !------------------------------------------------------------------------------!
  !         If all gradients in the direction of negative eigenvalues            !
  !              are below threshold, choose the best step from                  !
  !         (the purpose of this condition is to replace ineffective             !
  !                   rejection of Newton step - get_mu )                        !
  ! Newton step,                                        -->  rejected = false    !
  ! step on the Trust-region                            -->  rejected = true     !
  ! and step in the direction of positive eigenvalues.  -->  nstep    = false    !
  !------------------------------------------------------------------------------!
  nref = 0
  do i=1,nred
    if(dabs(tgrad(i)) > thresh)then
      nref = i
      exit
    endif
  enddo
  if(nref == 0)then
    if(nem == 0)then
      write(lupri,'(a)')" WARNING: We are in minimum. Why I'm here??? "
      write(lupri,'(a)')" WARNING:    ==>  STOP "
      stop
    else
      nstep = .false.
    endif
  else
    if((eval(nref) > d0).and.(nem/=0)) nstep = .false.
    if( dmaxa<10.0d0*thresh .or. dmaxn<10.0d0*thresh ) nstep = .true.
  endif

  if( .not.nstep )then
    cond(3:3) = '!'
    if(restore)then          ! Newton step was restored
      rejected = .false.
      nstep    = .true.
      restore  = .false.
    elseif(nrun == 0)then    ! try side step
      diffg1   = diffg
      tradius  = tradold
      rejected = .true.
      nstep    = .false.
      nrun     = 1
    elseif( nrun == 1 .and. nref == 0 )then
      if( diffg1<diffg )then  ! Newton step                           
        tradius  = tradold                                                                    
        rejected = .false.                                                                    
        nstep    = .true.                                                                     
        restore  = .true.                                                                     
        nrun     = 0                                                                          
      else                    ! side step                             
        tradius  = tradold                                                                    
        rejected = .true.                                                                     
        nstep    = .false.                                                                    
        nrun     = 3                                                                          
      endif
    elseif(nrun == 1)then    ! try step on the Trust-region
      diffg2   = diffg
      tradius  = tradold
      rejected = .true.
      nstep    = .true.
      nrun     = 2
    elseif(nrun == 2)then    ! choose the best step
      if( diffg<d0 .or. diffg1<d0 .or. diffg2<d0 )then
        if( diffg<diffg1 .and. diffg<diffg2 )then       ! step on the Trust-region
          rejected = .false.
          nstep    = .true.
          nrun     = 0
        elseif( diffg1<diffg .and. diffg1<diffg2 )then  ! Newton step
          tradius  = tradold
          rejected = .false.
          nstep    = .true.
          restore  = .true.
          nrun     = 0
        elseif( diffg2<diffg .and. diffg2<diffg1 )then  ! side step
          tradius  = tradold
          rejected = .true.
          nstep    = .false.
          nrun     = 3
        endif
      else
        write(lupri,'(a)')" WARNING: Both regular Newton step and side step were rejected! "
        write(lupri,'(a)')" WARNING:    ==>  STOP "
        stop
      endif
    elseif(nrun == 3)then
      rejected = .false.
      nstep    = .true.
      nrun     = 0
    endif
  endif

  !------------------------------------------------------------------------------!
  !                       Check if we are in the flat area                       !
  !------------------------------------------------------------------------------!
  if( dabs(diffg)<thresh .and. dabs(diffq)<thresh .and. nem/=0 .and. &
      ( dmaxa>10.0d0*thresh .and. dmaxn>10.0d0*thresh ) .and. cond(3:3) /= '!' )then
    rejected = .true.
    if(flat) tradius = 10.0d0*tradius
    flat     = .true.
    cond(3:3) = '.'
  else
    flat     = .false.
  endif

  !------------------------------------------------------------------------------!
  !                Check if we are in the saddle point or maximum                !
  !------------------------------------------------------------------------------!
  if(dmaxa<thresh .and. nem/=0 .and. cond(3:3) /= '!')then
    rejected  = .true.
    if(cond(3:3) == ' ') cond(3:3) = 's'
    if(cond(3:3) == '.') cond(3:3) = 'd'
  endif

  if(tradius /= tradold) trdef = .false.

  write(lupri,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 4x,a, 3x,a, 4x,a, 5x,a, 6x,a)') &
    '  ','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'
  write(lupri,'(a,2x,i4,3d11.2,f11.5,3x,l1,3x,a,4d11.2)') &
    '*T',itr+1,tradold,tradius,dnrmk,R,rejected,cond,diffg,diffq


  end subroutine updateT

end module localization
