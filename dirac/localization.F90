module localization

use memory_allocator
implicit none

public :: ORBLOC, initialize_localization

private :: ORBLP1, ORBLP2, ORBLO1, ORBLO4, ORBLO5, ROTCMO_MOD
private :: get_kappa, get_hessian, get_gradient, exp_func, hdiag
private :: get_mu, updateT

type :: localization_flags
  integer :: iprloc, itrloc
  real*8  :: thfull, thdiag
  logical :: prjloc, lgfull, lgdiag
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
  integer :: thrzer
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
real*8, parameter, private :: tresh = 1.0d-13, trsh_mu = 1.0d-4
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
  loc%prjloc     = prjloc
  loc%lgfull     = lgfull
  loc%lgdiag     = lgdiag
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
  character    :: daytid*24, hessian*4 = 'DIAG', chstart*23
  logical      :: tobe, full, rejected, Iminmax, diagconv

  integer :: ifrp, ndimh, ndimhr, idum, itr, ist, lupri, nish, norb, nucdep, maxref
  integer :: nfbas(0:2), ipqtoq(4,0:7), nep, nem, nez

  real*8  :: dinvmo, dinvmg, dinvmq, dnrmk, dnrmg, diffg, diffq, diffa
  real*8  :: timend, time, timstr, second, toterg, totdiag, totnond, shift, tradius, tpred

  integer, allocatable, dimension(:)       :: nfnuc, isel, indao, indaos, indnao
  real*8,  allocatable, dimension(:)       :: scoef, storec, smat, cmo, grad, etmp
  real*8,  allocatable, dimension(:,:)     :: hmat, evec
  real*8,  allocatable, dimension(:,:,:)   :: xkappa
  real*8,  allocatable, dimension(:,:,:,:) :: pmat

  real*8, pointer :: pkappa(:), tgrad(:), eval(:)

  !------------------------------------------------------------------------------!
  call qenter('orbloc')

  ! Set some parameters
  ifrp   = 1
  lupri  = var%lupri
  nish   = var%nish(ifrp)
  norb   = var%norb(ifrp)
  nucdep = var%nucdep
  maxref = var%maxref

  nfbas(:)    = var%nfbas(ifrp,:)
  ipqtoq(:,:) = var%ipqtoq(:,:)

  full = .false.
  if(loc%hesloc == 'FULL') full = .true.

  rejected = .false.
  Iminmax  = .false.
  chstart = 'START ITERATION NO.    '

  tradius = 10.0d0

  if(loc%iprloc.ge.4)then
    write(lupri,*) 'lupri =',lupri
    write(lupri,*) 'nish  =',nish
    write(lupri,*) 'norb  =',norb
    write(lupri,*) 'nucdep=',nucdep
    write(lupri,*) 'maxref=',maxref
    write(lupri,*) 'full  =',full
    write(lupri,*) 'thdiag=',loc%thdiag
    write(lupri,*) 'thfull=',loc%thfull
    write(lupri,*) 'lgdiag=',loc%lgdiag
    write(lupri,*) 'lgfull=',loc%lgfull
    write(lupri,*) 'hesloc=',loc%hesloc
  endif

  call titler('Molecular orbital localization','*',122)

  ! Check if coefficients are on file
  inquire(file='DFCOEF',exist=tobe)
  if(.not.tobe)then
    write(lupri,'(a)') 'orbloc: Coefficient file not found !'
    call quit('orbloc: Coefficients not found !')
  endif

  ! Memory allocation
  ndimh = 4 * ( nish*nish + nish ) / 2
  ndimhr= ndimh - nish
  call alloc(scoef ,nish*nfbas(0)*4,   id='scoef  in orbloc')
  call alloc(storec,nish*nfbas(0)*4,   id='storec in orbloc')
  call alloc(smat  ,var%n2bbasx,       id='smat   in orbloc')
  call alloc(cmo   ,var%n2bbasxq,      id='cmo    in orbloc')
  call alloc(hmat  ,ndimhr,ndimhr,     id='hmat   in orbloc')
  call alloc(grad  ,ndimhr,            id='grad   in orbloc')
  call alloc(etmp  ,ndimhr,            id='etmp   in orbloc')
  call alloc(evec  ,ndimhr,ndimhr,     id='evec   in orbloc')
  call alloc(xkappa,nish,nish,4,       id='xkappa in orbloc')
  call alloc(pmat  ,nish,nish,4,nucdep,id='pmat   in orbloc')

  call alloc(nfnuc ,2*nucdep,  id='nfnuc  in orbloc')
  call alloc(isel  ,nish,      id='isel   in orbloc')
  call alloc(indao ,nfbas(0),  id='indao  in orbloc')
  call alloc(indaos,nfbas(0),  id='indaos in orbloc')
  call alloc(indnao,2*nucdep+2,id='indnao in orbloc')

  ! nullify pointers for later allocation using memory_allocator module
  nullify(tgrad )
  nullify(eval  )
  nullify(pkappa)

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

  timstr = second()
  if(loc%prjloc) then
    ! PRJLOC = TRUE localization with projection analysis
    call orblp1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,maxref,dinvmg)
  else
    ! PRJLOC = FALSE localization with Mulliken analysis
    call orblo1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,nfnuc,dinvmg)
  endif
  write(lupri,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 2x,a, 3x,a, 5x,a, 6x,a, 6x,a)') &
    '*C','ITR','DIFFA','|GRAD|','|KAPPA|','H(nonD/D)','e: +,-,0','H','DIFFG','DIFFQ'
  write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 4x,a, 5x,a, 6x,a)') &
    '*T','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'

  ! SCF cycle
  dinvmq = d0
  itr    = 0
  ist    = 0
  do

    ist = ist + 1
    ! get new kappa
    if(.not.rejected)then
      call get_hessian(pmat,hmat,totdiag,totnond)
      call get_gradient(pmat,grad,dnrmg)
      call hdiag(hmat,grad,etmp,tgrad,eval,evec,nep,nem,nez,Iminmax,full)
      shift = d0
    else
      call get_mu(eval,tgrad,tradius,shift,tpred)
    endif
    call get_kappa(xkappa,pkappa,tgrad,eval,evec,shift,dnrmk)
    call exp_func(pkappa,tgrad,eval,nish,dinvmq)
    if( (rejected) .and. (dabs(tpred-dnrmk) > trsh_mu*tradius) )then
      write(lupri,'(a)') ' WARNING: Predicted and actual norm of kappa are too different'
      write(lupri,'(a,i4,a,d12.4,a,d12.4,a,2d12.4)')&
        ' WARNING: itr:',itr,' Predicted:',tpred,' Actual:',dnrmk,' Trust-radius:',tradius
    endif

    if(.not.rejected) storec(:) = scoef(:)  ! store MO coeff. in the case you reject the step
    if(     rejected) scoef(:)  = storec(:) ! restore MO coeff.; previous step was rejected
    call rotcmo_mod(scoef,nfbas(0),nish,4,xkappa,nish,ipqtoq(ifrp,0))
    if(loc%prjloc) then
      ! PRJLOC = TRUE localization with projection analysis
      call gtovlx(smat,var%ssmtrc)
      call orblp1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,maxref,dinvmg)
    else
      ! PRJLOC = FALSE localization with Mulliken analysis
      call orblo1(smat,scoef,pmat,nfbas(1),nfbas(2),nish,nucdep,nfnuc,dinvmg)
    endif

    ! Update Trust-radius
    diffg = dinvmg - dinvmo
    diffq = dinvmq
    call updateT(dnrmk,diffg,diffq,tradius,rejected,full,itr)

    if(.not.rejected)then
      itr = itr + 1
      call gtinfo(daytid)

      write(chstart(20:23),'(i4)') itr
      call titler(chstart,'=',122)
      write(lupri,'(/6x,a,a24)')'>>> ',daytid

      if (diffg < d0) write(lupri,'(a,i4)')' WARNING: dinvm_new - dinvm_old < 0; itr =',itr
      diffa = dabs(dinvmg - dinvmo)
      if(full) hessian = 'FULL'
      dinvmo = dinvmg

      write(lupri,'(/a,2x,a4, 4x,a, 6x,a, 5x,a, 2x,a, 3x,a, 5x,a, 6x,a, 6x,a)') &
        '  ','ITR','DIFFA','|GRAD|','|KAPPA|','H(nonD/D)','e: +,-,0','H','DIFFG','DIFFQ'
      write(lupri,'(a,2x,i4,3d11.2,f10.4,1x,3i4,3x,a,2d11.2)') &
        '*C',itr,diffa,dnrmg,dnrmk,totnond/totdiag,nep,nem,nez,hessian,diffg,diffq

      ! Exit conditions
      if(loc%hesloc == 'DIAG')then
        if(Iminmax) exit
        if(loc%lgdiag.and.(diffa < loc%thdiag)) exit
      endif

      if(loc%hesloc == 'FULL')then
        if(Iminmax) exit
        if(loc%lgfull.and.(diffa < loc%thfull)) exit
      endif

      if(loc%hesloc == 'COMB')then
        if(full)then
          if(Iminmax) exit
          if(loc%lgfull.and.(diffa < loc%thfull)) exit
        else
          if(.not.loc%lgdiag.and.Iminmax) full = .true.
          if(     loc%lgdiag.and.(diffa < loc%thdiag)) full = .true.
        endif
      endif

      if(itr > loc%itrloc)exit
      ! Exit conditions

    endif

  enddo

  write(lupri,'(//11x,a)')      '************************************************************'
  if(Iminmax)then
    write(lupri,'(11x,a,i4,a,i4,a)') &
         '>>>> MAXIMUM reached after',itr,' iterations and',ist,' steps <<<<'
  elseif((full.and.(diffa < loc%thfull)).or.(.not.full.and.(diffa < loc%thdiag)))then
    write(lupri,'(11x,a,i4,a,i4,a)') &
         '>>>>    Convergence in',itr,' iterations and',ist,' steps     <<<<'
  else
    write(lupri,'(11x,a,i4,a,i4,a)') &
         '>>>> No convergence after',itr,' iterations and',ist,' steps  <<<<'
  endif
  write(lupri,'(11x,a,f21.15,a)')'>>>> Inverse mean delocalization =',dinvmo,' <<<<'
  write(lupri,'(11x,a//)')      '************************************************************'

  timend = second()
  time   = timend - timstr
  call timtxt('          >>> Time used in localization is :',time,lupri)

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
  call dealloc(hmat  ,id='hmat   in orbloc')
  call dealloc(grad  ,id='grad   in orbloc')
  call dealloc(etmp  ,id='etmp   in orbloc')
  call dealloc(evec  ,id='evec   in orbloc')
  call dealloc(cmo   ,id='cmo    in orbloc')

  call dealloc(nfnuc ,id='nfnuc  in orbloc')
  call dealloc(isel  ,id='isel   in orbloc')
  call dealloc(indao ,id='indao  in orbloc')
  call dealloc(indaos,id='indaos in orbloc')
  call dealloc(indnao,id='indnao in orbloc')

  call deallocp(pkappa,id='pkappa in orbloc')
  call deallocp(tgrad ,id='tgrad  in orbloc')
  call deallocp(eval  ,id='eval   in orbloc')

  call qexit('orbloc')

  end subroutine orbloc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE ORBLP1(SMAT,SCOEF,PMAT,NBL,NBS,NI,NUCDEP,MAXREF,DINVM)
!***********************************************************************
!
!     Based on T.Saue March 26 1999 : SUBROUTINE ANAPRJ
!     Modify by S. Dubillard - Apr 21 2004
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      USE SELECT_COEFF_OWN_BAS
      IMPLICIT NONE

      ! External variables
      INTEGER*4, INTENT(IN ) :: NBL, NBS, NI, NUCDEP, MAXREF
      REAL*8,    INTENT(OUT) :: DINVM

      REAL*8,    INTENT(IN   ) :: SCOEF(*)
      REAL*8,    INTENT(INOUT) :: SMAT(*)
      REAL*8,    INTENT(  OUT) :: PMAT(NI,NI,4,NUCDEP)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      REAL*8, PARAMETER :: DM1=-1.0D0

      LOGICAL   :: TOBE

      INTEGER*4 :: I, J, K, IZ, IER, IFRP, ICFS, NB, NQ, NCDIM, NBDIM, NMOL
      INTEGER*4 :: IVEC, JOFF, NREF, N2REF, LUPRI, NPSH, NESH, NORB, NREFS
      INTEGER*4 :: NSTR(2,0:2,0:MAXREF), NVECS(2,0:2), IOFF(2), IDUMARR(2), IPQTOQ(4,0:7)

      REAL*8    :: DDUMARR(2), DINV, QLOC2, DLOC

      INTEGER*4, ALLOCATABLE, DIMENSION(:) :: KVECTMP, JREF, IPIV
      REAL*8,    ALLOCATABLE, DIMENSION(:) :: RCOEF, CMO, BV, GVEC, AMAT, WMAT

      TYPE :: PTR2ARRAY
        INTEGER*4, DIMENSION(:), POINTER :: P
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

      CALL QENTER('ORBLP1')
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
        WRITE(LUPRI,'(11X,A,I3,A,2F7.4)') &
        '* Orbital ', I ,' Pipek-Mezey localization measure d = ',DLOC
      ENDDO
      DINVM = DINV/DFLOAT(NI)
      WRITE(LUPRI,'(/11X,A,F21.15)') '* Inverse mean delocalization = ',DINVM

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

      CALL QEXIT('ORBLP1')

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
      INTEGER*4, INTENT(IN ) :: NB, NMOL, NREF, MAXREF
      INTEGER*4, INTENT(IN ) :: NSTR(2,0:2,0:MAXREF)
      INTEGER*4, INTENT(OUT) :: IPIVOT(*)

      REAL*8,    INTENT(IN   ) :: SCOEF(*), RCOEF(*)
      REAL*8,    INTENT(INOUT) :: SMAT(*)
      REAL*8,    INTENT(  OUT) :: AMAT(NREF,NREF,4), BVEC(NREF,NMOL,4)

      ! Internal variables
      CHARACTER*6  :: REFFIL(MAXREF)
      CHARACTER*72 :: VECREF(2,MAXREF)
      INTEGER*4    :: I, J, IFRP, JOFF, KOFF, JOB, NEFF, NREFS, LUPRI
      INTEGER*4    :: IPQTOQ(4,0:7)

      REAL*8    :: PTOL
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
      USE MEMORY_ALLOCATOR
      IMPLICIT NONE

      ! External variables
      INTEGER*4, INTENT(IN ) :: NBL, NBS, NI, NUCDEP
      REAL*8,    INTENT(OUT) :: DINVM

      INTEGER*4, INTENT(IN ) :: NFNUC(*)
      REAL*8,    INTENT(IN ) :: SMAT(*), SCOEF(*)
      REAL*8,    INTENT(OUT) :: PMAT(NI,NI,4,NUCDEP)

      ! Internal variables
      INTEGER*4 :: I, J, IZ, NB, IOFFL, IOFFS, LUPRI
      INTEGER*4 :: IPQTOQ(4,0:7)
      REAL*8    :: DINV, QLOC2, DLOC

      REAL*8, ALLOCATABLE, DIMENSION(:) :: WLMAT, WSMAT

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLO1')

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
!
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
            QLOC2 = QLOC2 + PMAT(I,I,1,J)*PMAT(I,I,1,J)/D4
         ENDDO
         DLOC = D1/QLOC2
         DINV = DINV + QLOC2
         WRITE(LUPRI,'(11X,A,I3,A,F7.4)') &
         '* Orbital ', I ,' Pipek-Mezey localization measure d = ', DLOC
      ENDDO
      DINVM = DINV/DFLOAT(NI)
      WRITE(LUPRI,'(/11X,A,F21.15)') '* Inverse mean delocalization = ',DINVM

  !------------------------------------------------------------------------------!
  !                             Deallocate arrays                                !
  !------------------------------------------------------------------------------!
      CALL DEALLOC(WLMAT,ID='WLMAT IN ORBLO1')
      CALL DEALLOC(WSMAT,ID='WSMAT IN ORBLO1')

      CALL QEXIT('ORBLO1')

  END SUBROUTINE ORBLO1

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_kappa(xkappa,pkappa,tgrad,eval,evec,shift,dnrmk)
  !******************************************************************************!
  !                                                                              !
  !           Construction of the XKAPPA matrix using full hessian               !
  !          The linear set of equations is solved using information             !
  !                   from diagonalization of the hessian                        !
  !                                                                              !
  !   (H-shift) * X = -g     -->     (tilde{H}-shift) * tilde{X} = - tilde{g}    !
  !                                                                              !
  !   tilde{H} = U^T * H * U = eval                 U = evec                     !
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
  real*8, intent(in ) :: shift
  real*8, intent(out) :: dnrmk

  real*8, intent(in ), dimension(:,:)     :: evec
  real*8, intent(out), dimension(:,:,:)   :: xkappa

  real*8, intent(in ), pointer :: tgrad(:), eval(:)
  real*8, intent(out), pointer :: pkappa(:)

  ! Internal variables
  integer :: i, j, ij, iz, ni, nhtmp, ni12, nred, lupri
  real*8  :: dnrm2

  real*8, allocatable, dimension(:) :: xtmp

  !------------------------------------------------------------------------------!
  call qenter('get_kappa')

  lupri  = var%lupri

  ! Get some dimensions
  ni    = size(xkappa,1)
  nhtmp = size(evec,1)
  ni12  = (nhtmp+ni) / 4
  nred  = size(eval,1)

  ! allocate arrays
  if(associated(pkappa)) call deallocp(pkappa,id='pkappa in get_kappa')
  call allocp(pkappa,nred ,id='pkappa in get_kappa')
  call alloc(xtmp  ,nhtmp,id='xtmp   in get_kappa')

  ! calculate tilde{X}
  do i=1,nred
    pkappa(i) = - tgrad(i) / (eval(i) - shift)
  enddo

  ! calculate X
  xtmp(:) = d0
  do j=1,nred
    do i=1,nhtmp
      xtmp(i) = xtmp(i) + evec(i,j)*pkappa(j)
    enddo
  enddo

  ! Norm of the X
  dnrmk = dnrm2(nhtmp,xtmp,1)
  if(loc%iprloc.ge.1)then
    write (lupri,'(/a,e12.2)') ' Norm of the kappa: ',dnrmk
  endif

  ! Unpack X to kappa matrix
  xkappa(:,:,:) = d0
  ij = 0
  do j = 1,ni
    do i = j,ni
      if(i==j)cycle
      ij = ij + 1
      xkappa(j,i,1) =   xtmp(ij)
      xkappa(i,j,1) = - xkappa(j,i,1)
    enddo
  enddo

  ij = ni12-ni
  do j = 1,ni
    do i = j,ni
      ij = ij + 1
      xkappa(j,i,2) =   xtmp(ni12+ij)
      xkappa(i,j,2) =   xkappa(j,i,2)

      xkappa(j,i,3) = - xtmp(ij)
      xkappa(i,j,3) =   xkappa(j,i,3)

      xkappa(j,i,4) =   xtmp(2*ni12+ij)
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

  !------------------------------------------------------------------------------!
  call dealloc(xtmp,id='xtmp in get_kappa')

  call qexit('get_kappa')

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
      INTEGER*4, INTENT(IN ) :: NUCDEP

      INTEGER*4, INTENT(IN ) :: INDNAO(2,*)
      INTEGER*4, INTENT(OUT) :: NFNUC(2*NUCDEP)

      ! Internal variables
      INTEGER*4 :: IFRP, I

  !------------------------------------------------------------------------------!
      IFRP = 1
      CALL QENTER('ORBLO4')
      DO I = 1,NUCDEP-1
         NFNUC(I)=INDNAO(1,I+1)-INDNAO(1,I)
      ENDDO
      NFNUC(NUCDEP)=INDNAO(2,1)-INDNAO(1,NUCDEP)
      DO I = 1,NUCDEP-1
         NFNUC(I+NUCDEP)=INDNAO(2,I+1)-INDNAO(2,I)
      ENDDO
      NFNUC(2*NUCDEP)=VAR%NFBAS(IFRP,0)-INDNAO(2,NUCDEP)
      CALL QEXIT('ORBLO4')

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
      USE MEMORY_ALLOCATOR
      IMPLICIT NONE

      ! External variables
      INTEGER*4 :: NB, NI, NO
      REAL*8    :: SCOEF(NB,NI,4)

      ! Internal variables
      INTEGER*4 :: I, J, K, L
      REAL*8    :: TOTERG

      INTEGER*4, ALLOCATABLE :: IBEIG(:)
      REAL*8,    ALLOCATABLE :: CMO(:,:,:), EIG(:)

  !------------------------------------------------------------------------------!
      CALL QENTER('ORBLO5')

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

      CALL QEXIT('ORBLO5')

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
      USE MEMORY_ALLOCATOR
      IMPLICIT NONE

      ! External variables
      INTEGER*4, INTENT(IN   ) :: NBASI, NORBI, NZ, NORBT, IP(4)
      REAL*8,    INTENT(IN   ) :: XKAPPA(NORBT,NORBT,NZ)
      REAL*8,    INTENT(INOUT) :: CMO(NBASI,NORBI,NZ)

      ! Internal variables
      REAL*8, PARAMETER :: DTHRS = 1.0D-8, DSIXTH = D1/6.00D00, DP5 = 0.50D00
      REAL*8, PARAMETER :: D1P5 = 1.5D00, THREQL = 1.0D-12, DM1 = -D1

      CHARACTER SECTID*12,CPUTID*12,WALLTID*12

      INTEGER*4 :: I, J, NORBXQ, IERR, IADDR, IOFF, I_MAX, J_MAX, LUPRI
      REAL*8    :: CPU1, WALL1, CPU2, WALL2, CPU3, WALL3, CPU4, WALL4, EIG, SQEIG, D_MAX
      REAL*8    :: ABS, SQRT, COS, SIN

      REAL*8, ALLOCATABLE, DIMENSION(:) :: REIG, RSIN, RCOS, RV, RTMP

  !------------------------------------------------------------------------------!
      CALL QENTER('ROTCMO_MOD')

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

      CALL QEXIT('ROTCMO_MOD')

  END SUBROUTINE ROTCMO_MOD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_hessian(pmat,hmat,totdiag,totnond)
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
  real*8, intent(out) :: totdiag, totnond

  real*8, intent(in ), dimension(:,:,:,:) :: pmat
  real*8, intent(out), dimension(:,:)     :: hmat

  ! Internal variables
  real*8, parameter :: dhalf = 0.5d0

  real*8  :: pij1, pji1, pkl1, pjl1, plj1, pki1, pik1, pjk1, pil1, pli1, pkj1
  real*8  :: pij2, pji2, pkl2, pjl2, plj2, pki2, pik2, pjk2, pil2, pli2, pkj2
  real*8  :: pii, pjj, pkk, pll, sumdiag, sumnond, sumtmp
  real*8  :: dtmp1, dtmp2, dtmp3, dtmp4, dtmp5, dtmp6

  integer :: i, j, i1, j1, k, l, m, ij, kl, nucdep, ni, ndimh, ni12, lupri

  !------------------------------------------------------------------------------!
  call qenter('get_hessian')

  lupri = var%lupri

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimh  = size(hmat,1)
  ni12   = (ndimh+ni) / 4

  !------------------------------------------------------------------------------!
  ! Calculate real representation of hessian and get the diagonal dominance info !
  !------------------------------------------------------------------------------!
  hmat(:,:) = d0
  totdiag = d0
  totnond = d0

  !------------------------------------------------------------------------------!
  !                 Unbared-unbared contribution to the hessian                  !
  !------------------------------------------------------------------------------!
  ! A^R + B^R = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
  !           = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK_ij dK_kl }
  !           = Real{ (32) }                 + Real{ (36) }
  sumdiag = d0
  sumnond = d0
  ij = 0
  do j=1,ni
    do i=j,ni
      if(i==j) cycle
      ij = ij + 1
      kl = 0
      do l=1,ni
        do k=l,ni
          if(k==l) cycle
          kl = kl + 1
          do m = 1,nucdep
            pij1 = pmat(i,j,1,m); pkl1 = pmat(k,l,1,m)
            pjl1 = pmat(j,l,1,m); pil1 = pmat(i,l,1,m)
            pki1 = pmat(k,i,1,m); pkj1 = pmat(k,j,1,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d4*pij1*pkl1

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * pjl1
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp1 - d2 * ( pjj - dhalf*(pii + pll) ) * pil1
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp1 - d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
          enddo
          if(ij==kl)then
            sumdiag = sumdiag + dabs( hmat(ij,kl) )
          else
            sumnond = sumnond + dabs( hmat(ij,kl) )
          endif
        enddo
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_uu'
    do j=1,ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=1,ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^R + B^R)_uu elements (diag. vs. non-diag): ',sumdiag,sumnond
  endif
  totdiag = totdiag + sumdiag
  totnond = totnond + sumnond

  ! A^I + B^I = Imag{ d^2 G / dK*_ij dK_kl } + Imag{ d^2 G / dK*_ij dK*_kl }
  !           = Imag{ d^2 G / dK*_ij dK_kl } - Imag{ d^2 G / dK_ij dK_kl }
  !           = Imag{ (32) }                 - Imag{ (36) }
  sumnond = d0
  ij = 2*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          if(k==l) cycle
          kl = kl + 1
          do m = 1,nucdep
            pij2 = pmat(i,j,2,m); pkl1 = pmat(k,l,1,m)
            pjl2 = pmat(j,l,2,m); pil2 = pmat(i,l,2,m)
            pki2 = pmat(k,i,2,m); pkj2 = pmat(k,j,2,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp2= -d4*pij2*pkl1

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pii - dhalf*(pjj + pll) ) * pjl2
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki2
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pjj - dhalf*(pii + pll) ) * pil2
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj2
          enddo
          sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        enddo
      enddo
    enddo
  enddo

  ! A^I + B^I =   Imag{ d^2 G / dK*_kk dK_ij } + Imag{ d^2 G / dK*_kk dK*_ij }
  !           =   Imag{ d^2 G / dK*_kk dK_ij } - Imag{ d^2 G / dK_kk dK_ij }
  !           = - Imag{ (40) }                 - Imag{ (44) }
  ij = 0
  do j=1,ni
    do i=j,ni
      if(i==j) cycle
      ij = ij + 1
      kl = 2*ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij2 = pmat(i,j,2,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(kl,ij) = hmat(kl,ij) + pij2 * ( pii - pjj )
          if(j==k) hmat(kl,ij) = hmat(kl,ij) - pij2 * ( pii - pjj )
        enddo
        sumnond = sumnond + d2*dabs( hmat(kl,ij) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^I + B^I)_uu'
    do j=2*ni12-ni+1,3*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=1,ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^I + B^I)_uu elements (non-diag): ',sumnond
  endif
  totnond = totnond + sumnond

  ! A^R - B^R = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
  !           = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK_ij dK_kl }
  !           = Real{ (32) }                 - Real{ (36) }
  sumdiag = d0
  sumnond = d0
  ij = 2*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 2*ni12 - ni
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij2 = pmat(i,j,2,m); pkl2 = pmat(k,l,2,m)
            pjl1 = pmat(j,l,1,m); pil1 = pmat(i,l,1,m)
            pki1 = pmat(k,i,1,m); pkj1 = pmat(k,j,1,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d4*pij2*pkl2

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * pjl1
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pil1
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
          enddo
          if(ij==kl)then
            sumdiag = sumdiag + dabs( hmat(ij,kl) )
          else
            sumnond = sumnond + dabs( hmat(ij,kl) )
          endif
        enddo
      enddo
    enddo
  enddo

  ! A^R - B^R = Real{ d^2 G / dK*_ij dK_kk } - Real{ d^2 G / dK*_ij dK*_kk }
  !           = Real{ d^2 G / dK*_ij dK_kk } - Real{ d^2 G / dK_ij dK_kk }
  !           = Real{ (40) }                 - Real{ (44) }
  ! A^R - B^R = Real{ d^2 G / dK*_kk dK_ij } - Real{ d^2 G / dK*_kk dK*_ij }
  !           = Real{ d^2 G / dK*_kk dK_ij } - Real{ d^2 G / dK_kk dK_ij }
  !           = Real{ (40) }                 - Real{ (44) }
  ij = 2*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 2*ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,1,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(ij,kl) = hmat(ij,kl) + pij1 * ( pii - pjj )
          if(j==k) hmat(ij,kl) = hmat(ij,kl) - pij1 * ( pii - pjj )
        enddo
        hmat(kl,ij) = hmat(ij,kl)
        sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_uu'
    do j=2*ni12-ni+1,3*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=2*ni12-ni+1,3*ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^R - B^R)_uu elements (diag. vs. non-diag): ',sumdiag,sumnond
  endif
  totdiag = totdiag + sumdiag
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                   Bared-unbared contribution to the hessian                  !
  !------------------------------------------------------------------------------!
  ! A^R + B^R = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
  !           = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK_ij dK_kl }
  !           = Real{ (34) }                 + Real{ (38) }
  sumnond = d0
  ij = ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          if(k==l) cycle
          kl = kl + 1
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pkl1 = pmat(k,l,1,m)
            plj1 = pmat(l,j,3,m); pkj1 = pmat(k,j,3,m)
            pli1 = pmat(l,i,3,m); pki1 = pmat(k,i,3,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d4*pij1*pkl1

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * plj1
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp1 - d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pli1
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp1 - d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
          enddo
          sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        enddo
      enddo
    enddo
  enddo

  ! A^R + B^R = Real{ d^2 G / dK*_kk dK_ij } + Real{ d^2 G / dK*_kk dK*_ij }
  !           = Real{ d^2 G / dK*_kk dK_ij } + Real{ d^2 G / dK_kk dK_ij }
  !           = Real{ (41) }                 + Real{ (45) }
  ij = 0
  do j=1,ni
    do i=j,ni
      if(i==j) cycle
      ij = ij + 1
      kl = ni12 -ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(kl,ij) = hmat(kl,ij) - pij1 * ( pii - pjj )
          if(j==k) hmat(kl,ij) = hmat(kl,ij) + pij1 * ( pii - pjj )
        enddo
        sumnond = sumnond + d2*dabs( hmat(kl,ij) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_bu'
    do j=ni12-ni+1,2*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=1,ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^R + B^R)_bu elements (non-diag): ',sumnond
  endif
  totnond = totnond + sumnond

  ! A^I + B^I = Imag{ d^2 G / dK*_ij dK_kl } + Imag{ d^2 G / dK*_ij dK*_kl }
  !           = Imag{ d^2 G / dK*_ij dK_kl } - Imag{ d^2 G / dK_ij dK_kl }
  !           = Imag{ (34) }                 - Imag{ (38) }
  sumnond = d0
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 0
      do l=1,ni
        do k=l,ni
          if(k==l) cycle
          kl = kl + 1
          do m = 1,nucdep
            pij2 = pmat(i,j,4,m); pkl1 = pmat(k,l,1,m)
            plj2 = pmat(l,j,4,m); pkj2 = pmat(k,j,4,m)
            pli2 = pmat(l,i,4,m); pki2 = pmat(k,i,4,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp2= -d4*pij2*pkl1

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp2 - d2 * ( pii - dhalf*(pjj + pll) ) * plj2
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki2
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp2 - d2 * ( pjj - dhalf*(pii + pll) ) * pli2
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj2
          enddo
          sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        enddo
      enddo
    enddo
  enddo

  ! A^I + B^I =   Imag{ d^2 G / dK*_kk dK_ij } + Imag{ d^2 G / dK*_kk dK*_ij }
  !           =   Imag{ d^2 G / dK*_kk dK_ij } - Imag{ d^2 G / dK_kk dK_ij }
  !           = - Imag{ (41) }                 - Imag{ (45) }
  ij = 0
  do j=1,ni
    do i=j,ni
      if(i==j) cycle
      ij = ij + 1
      kl = 3*ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(kl,ij) = hmat(kl,ij) + pij2 * ( pii - pjj )
          if(j==k) hmat(kl,ij) = hmat(kl,ij) - pij2 * ( pii - pjj )
        enddo
        sumnond = sumnond + d2*dabs( hmat(kl,ij) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^I + B^I)_bu'
    do j=3*ni12-ni+1,4*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=1,ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^I + B^I)_bu elements (non-diag): ',sumnond
  endif
  totnond = totnond + sumnond

  ! A^R - B^R = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
  !           = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK_ij dK_kl }
  !           = Real{ (34) }                 - Real{ (38) }
  sumnond = d0
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 2*ni12 - ni
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij2 = pmat(i,j,4,m); pkl2 = pmat(k,l,2,m)
            plj1 = pmat(l,j,3,m); pkj1 = pmat(k,j,3,m)
            pli1 = pmat(l,i,3,m); pki1 = pmat(k,i,3,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d4*pij2*pkl2

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * plj1
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pli1
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
          enddo
          sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        enddo
      enddo
    enddo
  enddo

  ! A^R - B^R = Real{ d^2 G / dK*_ij dK_kk } - Real{ d^2 G / dK*_ij dK*_kk }
  !           = Real{ d^2 G / dK*_ij dK_kk } - Real{ d^2 G / dK_ij dK_kk }
  !           = Real{ (42) }                 - Real{ (46) }
  ! A^R - B^R = Real{ d^2 G / dK*_kk dK_ij } - Real{ d^2 G / dK*_kk dK*_ij }
  !           = Real{ d^2 G / dK*_kk dK_ij } - Real{ d^2 G / dK_kk dK_ij }
  !           = Real{ (41) }                 - Real{ (45) }
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 2*ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,3,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(ij,kl) = hmat(ij,kl) + pij1 * ( pii - pjj )
          if(j==k) hmat(ij,kl) = hmat(ij,kl) + pij1 * ( pii - pjj )
        enddo
        hmat(kl+ni12,ij-ni12) = - hmat(ij,kl)
        sumnond = sumnond + d4*dabs( hmat(ij,kl) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_bu'
    do j=3*ni12-ni+1,4*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=2*ni12-ni+1,3*ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^R - B^R)_bu elements (non-diag): ',sumnond
  endif
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                   Bared-bared contribution to the hessian                    !
  !------------------------------------------------------------------------------!
  ! A^R + B^R = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK*_ij dK*_kl }
  !           = Real{ d^2 G / dK*_ij dK_kl } + Real{ d^2 G / dK_ij dK_kl }
  !           = Real{ (35) }                 + Real{ (39) }
  sumdiag = d0
  sumnond = d0
  ij = ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 - ni
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
            pkl1 = pmat(k,l,3,m); pkl2 = pmat(k,l,4,m)
            pki1 = pmat(k,i,1,m); pkj1 = pmat(k,j,1,m)
            pli1 = pmat(l,i,1,m); plj1 = pmat(l,j,1,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 + pij2*pkl2)

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * plj1
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pli1
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
          enddo
        if(ij==kl)then
          sumdiag = sumdiag + dabs( hmat(ij,kl) )
        else
          sumnond = sumnond + dabs( hmat(ij,kl) )
        endif
        enddo
      enddo
    enddo
  enddo

  ! A^R + B^R = Real{ d^2 G / dK*_ij dK_kk } + Real{ d^2 G / dK*_ij dK*_kk }
  !           = Real{ d^2 G / dK*_ij dK_kk } + Real{ d^2 G / dK_ij dK_kk }
  !           = Real{ (43) }                 + Real{ (47) }
  ! A^R + B^R = Real{ d^2 G / dK*_kk dK_ij } + Real{ d^2 G / dK*_kk dK*_ij }
  !           = Real{ d^2 G / dK*_kk dK_ij } + Real{ d^2 G / dK_kk dK_ij }
  !           = Real{ (43) }                 + Real{ (47) }
  ij = ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,1,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(ij,kl) = hmat(ij,kl) + pij1 * ( pii - pjj )
          if(j==k) hmat(ij,kl) = hmat(ij,kl) - pij1 * ( pii - pjj )
        enddo
        hmat(kl,ij) = hmat(ij,kl)
        sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^R + B^R)_bb'
    do j=ni12-ni+1,2*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=ni12-ni+1,2*ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^R + B^R)_bb elements (diag. vs. non-diag): ',sumdiag,sumnond
  endif
  totdiag = totdiag + sumdiag
  totnond = totnond + sumnond

  ! A^I + B^I = Imag{ d^2 G / dK*_ij dK_kl } + Imag{ d^2 G / dK*_ij dK*_kl }
  !           = Imag{ d^2 G / dK*_ij dK_kl } - Imag{ d^2 G / dK_ij dK_kl }
  !           = Imag{ (35) }                 - Imag{ (39) }
  sumnond = d0
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 - ni
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
            pkl1 = pmat(k,l,3,m); pkl2 = pmat(k,l,4,m)
            pki2 = pmat(k,i,2,m); pkj2 = pmat(k,j,2,m)
            pli2 = pmat(l,i,2,m); plj2 = pmat(l,j,2,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp2= d2*(pij1*pkl2 - pij2*pkl1)

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pii - dhalf*(pjj + pll) ) * plj2
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki2
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pjj - dhalf*(pii + pll) ) * pli2
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj2
          enddo
        sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        enddo
      enddo
    enddo
  enddo

  ! A^I + B^I =   Imag{ d^2 G / dK*_ij dK_kk } + Imag{ d^2 G / dK*_ij dK*_kk }
  !           =   Imag{ d^2 G / dK*_ij dK_kk } - Imag{ d^2 G / dK_ij dK_kk }
  !           =   Imag{ (43) }                 - Imag{ (47) }
  ! A^I + B^I =   Imag{ d^2 G / dK*_kk dK_ij } + Imag{ d^2 G / dK*_kk dK*_ij }
  !           =   Imag{ d^2 G / dK*_kk dK_ij } - Imag{ d^2 G / dK_kk dK_ij }
  !           = - Imag{ (43) }                 - Imag{ (47) }
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij2 = pmat(i,j,2,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(ij,kl) = hmat(ij,kl) + pij2 * ( pii - pjj )
          if(j==k) hmat(ij,kl) = hmat(ij,kl) + pij2 * ( pii - pjj )
        enddo
        hmat(kl+2*ni12,ij-2*ni12) = - hmat(ij,kl)
        sumnond = sumnond + d4*dabs( hmat(ij,kl) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^I + B^I)_bb'
    do j=3*ni12-ni+1,4*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=ni12-ni+1,2*ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^I + B^I)_bb elements (non-diag): ',sumnond
  endif
  totnond = totnond + sumnond

  ! A^R - B^R = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK*_ij dK*_kl }
  !           = Real{ d^2 G / dK*_ij dK_kl } - Real{ d^2 G / dK_ij dK_kl }
  !           = Real{ (35) }                 - Real{ (39) }
  sumdiag = d0
  sumnond = d0
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 3*ni12 - ni
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij1 = pmat(i,j,3,m); pij2 = pmat(i,j,4,m)
            pkl1 = pmat(k,l,3,m); pkl2 = pmat(k,l,4,m)
            pki1 = pmat(k,i,1,m); pkj1 = pmat(k,j,1,m)
            pli1 = pmat(l,i,1,m); plj1 = pmat(l,j,1,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp1= d2*(pij1*pkl1 + pij2*pkl2)
 
            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pii - dhalf*(pjj + pll) ) * plj1
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp1 + d2 * ( pjj - dhalf*(pii + pkk) ) * pki1
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pjj - dhalf*(pii + pll) ) * pli1
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp1 + d2 * ( pii - dhalf*(pjj + pkk) ) * pkj1
          enddo
        if(ij==kl)then
          sumdiag = sumdiag + dabs( hmat(ij,kl) )
        else
          sumnond = sumnond + dabs( hmat(ij,kl) )
        endif
        enddo
      enddo
    enddo
  enddo
 
  ! A^R - B^R = Real{ d^2 G / dK*_ij dK_kk } - Real{ d^2 G / dK*_ij dK*_kk }
  !           = Real{ d^2 G / dK*_ij dK_kk } - Real{ d^2 G / dK_ij dK_kk }
  !           = Real{ (43) }                 - Real{ (47) }
  ! A^R - B^R = Real{ d^2 G / dK*_kk dK_ij } - Real{ d^2 G / dK*_kk dK*_ij }
  !           = Real{ d^2 G / dK*_kk dK_ij } - Real{ d^2 G / dK_kk dK_ij }
  !           = Real{ (43) }                 - Real{ (47) }
  ij = 3*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 3*ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij1 = pmat(i,j,1,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
 
          if(i==k) hmat(ij,kl) = hmat(ij,kl) + pij1 * ( pii - pjj )
          if(j==k) hmat(ij,kl) = hmat(ij,kl) - pij1 * ( pii - pjj )
        enddo
        hmat(kl,ij) = hmat(ij,kl)
        sumnond = sumnond+d2*dabs( hmat(ij,kl) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^R - B^R)_bb'
    do j=3*ni12-ni+1,4*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=3*ni12-ni+1,4*ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^R - B^R)_bb elements (diag. vs. non-diag): ',sumdiag,sumnond
  endif
  totdiag = totdiag + sumdiag
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                 Unbared-bared contribution to the hessian                    !
  !------------------------------------------------------------------------------!
  ! A^I + B^I = Imag{ d^2 G / dK*_ij dK_kl } + Imag{ d^2 G / dK*_ij dK*_kl }
  !           = Imag{ d^2 G / dK*_ij dK_kl } - Imag{ d^2 G / dK_ij dK_kl }
  !           = Imag{ (33) }                 - Imag{ (37) }
  sumnond = d0
  ij = 2*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 - ni
      do l=1,ni
        do k=l,ni
          kl = kl + 1
          if(k==l) cycle
          do m = 1,nucdep
            pij2 = pmat(i,j,2,m); pkl1 = pmat(k,l,3,m)
            pik2 = pmat(i,k,4,m); pjk2 = pmat(j,k,4,m)
            pil2 = pmat(i,l,4,m); pjl2 = pmat(j,l,4,m)
            pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)
            pkk  = pmat(k,k,1,m); pll  = pmat(l,l,1,m)
            dtmp2= - d4*pij2*pkl1

            if(i==k) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pii - dhalf*(pjj + pll) ) * pjl2
            if(j==l) hmat(ij,kl) = hmat(ij,kl) + dtmp2 + d2 * ( pjj - dhalf*(pii + pkk) ) * pik2
            if(j==k) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pjj - dhalf*(pii + pll) ) * pil2
            if(i==l) hmat(ij,kl) = hmat(ij,kl) - dtmp2 + d2 * ( pii - dhalf*(pjj + pkk) ) * pjk2
          enddo
          sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        enddo
      enddo
    enddo
  enddo

  ! A^I + B^I =   Imag{ d^2 G / dK*_ij dK_kk } + Imag{ d^2 G / dK*_ij dK*_kk }
  !           =   Imag{ d^2 G / dK*_ij dK_kk } - Imag{ d^2 G / dK_ij dK_kk }
  !           =   Imag{ (41) }                 - Imag{ (45) }
  ij = 2*ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pji2 = pmat(j,i,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(ij,kl) = hmat(ij,kl) + pji2 * ( pii - pjj )
          if(j==k) hmat(ij,kl) = hmat(ij,kl) + pji2 * ( pii - pjj )
        enddo
        sumnond = sumnond + d2*dabs( hmat(ij,kl) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  ! A^I + B^I =   Imag{ d^2 G / dK*_kk dK_ij } + Imag{ d^2 G / dK*_kk dK*_ij }
  !           =   Imag{ d^2 G / dK*_kk dK_ij } - Imag{ d^2 G / dK_kk dK_ij }
  !           = - Imag{ (42) }                 - Imag{ (46) }
  ij = ni12 - ni
  do j=1,ni
    do i=j,ni
      ij = ij + 1
      if(i==j) cycle
      kl = 2*ni12 - ni + 1
      do k=1,ni
        do m = 1,nucdep
          pij2 = pmat(i,j,4,m)
          pii  = pmat(i,i,1,m); pjj  = pmat(j,j,1,m)

          if(i==k) hmat(kl,ij) = hmat(kl,ij) + pij2 * ( pii - pjj )
          if(j==k) hmat(kl,ij) = hmat(kl,ij) + pij2 * ( pii - pjj )
        enddo
        sumnond = sumnond + d2*dabs( hmat(kl,ij) )
        kl = kl + ni - (k-1)  ! diagonal
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.3)then
    write(lupri,'(/a,2f10.5)') 'Print (A^I + B^I)_ub'
    do j=2*ni12-ni+1,3*ni12-ni
      write(lupri,'(i4,15f10.5)') j,(hmat(j,k),k=ni12-ni+1,2*ni12-ni)
    enddo
  endif
  if(loc%iprloc.ge.1)then
    write(lupri,'(a,2f10.5)') '(A^I + B^I)_ub elements (non-diag): ',sumnond
  endif
  totnond = totnond + sumnond

  !------------------------------------------------------------------------------!
  !                     Get the rest of the hessian matrix                       !
  !------------------------------------------------------------------------------!
  ! 2,1  -->  1,2
  do i=1*ni12-ni+1,2*ni12-ni
    do j=1,ni12-ni
      hmat(j,i) = hmat(i,j)
    enddo
  enddo

  ! 43  -->  34
  do i=3*ni12-ni+1,4*ni12-ni
    do j=2*ni12-ni+1,3*ni12-ni
      hmat(j,i) = hmat(i,j)
    enddo
  enddo

  ! 34,12  -->  12,34
  do i=2*ni12-ni+1,4*ni12-ni
    do j=1,2*ni12-ni
      hmat(j,i) = hmat(i,j)
    enddo
  enddo

  call qexit('get_hessian')

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
  integer :: i, j, k, ij, ij1, ij2, ij3, ij4
  integer :: nucdep, ni, ndimg, ni12, lupri

  real*8  :: dtmp, dnrm2

  !------------------------------------------------------------------------------!
  call qenter('get_gradient')

  lupri = var%lupri
  gradient(:) = d0

  ! Get some dimensions
  nucdep = size(pmat,4)
  ni     = size(pmat,1)
  ndimg  = size(gradient,1)
  ni12   = (ndimg+ni) / 4

  !------------------------------------------------------------------------------!
  !               Calculate real representation of the gradient                  !
  !------------------------------------------------------------------------------!
  ! Do not calculate redundant parameters (diagonal real parameters of unbared kappa)
  ij = 0; ij1 = 0
  do j=1,ni
    do i=j,ni
      ij  = ij  + 1
      ij2 = ij  + ni12-ni
      ij3 = ij2 + ni12
      ij4 = ij3 + ni12
      if(i==j) cycle
      ij1 = ij1 + 1
      do k = 1,nucdep
        dtmp = d2*(pmat(j,j,1,k) - pmat(i,i,1,k))
        gradient(ij1) = gradient(ij1) + dtmp * pmat(i,j,1,k)
        gradient(ij2) = gradient(ij2) + dtmp * pmat(i,j,3,k)
        gradient(ij3) = gradient(ij3) - dtmp * pmat(i,j,2,k)
        gradient(ij4) = gradient(ij4) - dtmp * pmat(i,j,4,k)
      enddo
    enddo
  enddo

  if(loc%iprloc.ge.2)then
    write(lupri,'(/a,2f10.5)') 'Print gradient(:,1)  -->  gradient'
    write(lupri,'(9d12.4)') (gradient(j),j=1,ni12-ni)
    write(lupri,'(/a,2f10.5)') 'Print gradient(:,2)  -->  gradient'
    write(lupri,'(9d12.4)') (gradient(j),j=1+2*ni12-ni,3*ni12-ni)
    write(lupri,'(/a,2f10.5)') 'Print gradient(:,1)  -->  /bar{gradient}'
    write(lupri,'(9d12.4)') (gradient(j),j=1+ni12-ni,2*ni12-ni)
    write(lupri,'(/a,2f10.5)') 'Print gradient(:,2)  -->  /bar{gradient}'
    write(lupri,'(9d12.4)') (gradient(j),j=1+3*ni12-ni,4*ni12-ni)
  endif

  ! Norm of the gradient
  dnrmg = dnrm2(ndimg,gradient,1)

  call qexit('get_gradient')

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
  call qenter('exp_func')

  dinvmq = d0
  do i=1,size(eval,1)
    dinvmq = dinvmq + d2*pkappa(i)*tgrad(i) + pkappa(i)*eval(i)*pkappa(i)
  enddo
  dinvmq = dinvmq / dfloat(ni)

  call qexit('exp_func')

  end subroutine exp_func

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine hdiag(hmat,gradient,etmp,tgrad,eval,evec,nep,nem,nez,Iminmax,full)
  !******************************************************************************!
  !                                                                              !
  !                Diagonalize real representation of hessian                    !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  logical, intent(in ) :: full
  logical, intent(out) :: Iminmax
  integer, intent(out) :: nep, nem, nez
  real*8,  intent(in ), dimension(:)   :: gradient
  real*8,  intent(in ), dimension(:,:) :: hmat
  real*8,  intent(out), dimension(:)   :: etmp
  real*8,  intent(out), dimension(:,:) :: evec
  real*8,  intent(out), pointer :: tgrad(:), eval(:)

  ! Internal variables
  integer :: i, j, k, nhtmp, nred, ier, lupri, lwork, ilaenv, nb, icur
  real*8  :: acur

  integer, allocatable, dimension(:)   :: iorig
  real*8,  allocatable, dimension(:)   :: work

  !------------------------------------------------------------------------------!
  call qenter('hdiag')

  lupri = var%lupri

  ! Get some dimensions
  nhtmp = size(hmat,1)

  if(full)then

    nb = ilaenv(1,'dsytrd','u',nhtmp,-1,-1,-1)
    lwork = (nb+2)*nhtmp  ! see the comment in dsyev lapack routine
  
    call alloc(work,lwork,id='work in hdiag')
  
    ! prepare evec
    evec(1:nhtmp,1:nhtmp) = hmat(1:nhtmp,1:nhtmp)
  
    ! get eigenvalues (etmp) and eigenvectors (evec)
    call dsyev('V','U',nhtmp,evec,nhtmp,etmp,work,lwork,ier)
    if(ier /= 0)then
      write(lupri,*) ' ERROR: The dsyev LAPACK subroutine in hdiag finished unsuccessfully'
      write(lupri,*) ' info =',ier
      stop
    endif

    call dealloc(work,id='work in hdiag')

  else
  
    call alloc(iorig,nhtmp,id='iorig in hdiag')

    ! get eigenvalues (etmp) and eigenvectors (evec)
    ! We use diagonal approximation so it will just sort the diagonal of the hessian
    do i=1,nhtmp
      etmp(i) = hmat(i,i)
      iorig(i)= i
    enddo
    evec(:,:) = d0
    do i=1,nhtmp
      acur = etmp(i)
      icur = i
      do j=i+1,nhtmp
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

    ! eigenvectors will contain information how the eigenvalues were sorted
    do i=1,nhtmp
      evec(iorig(i),i) = d1
    enddo

    call dealloc(iorig,id='iorig in hdiag')

  endif

  ! allocate unnamed arrays (dimension is number of eigenvalues over some threshold)
  nred = 0
  do i=1,nhtmp
    if(dabs(etmp(i)) > tresh) nred = nred + 1
  enddo
  if(associated(eval))then
    call deallocp(tgrad,id='tgrad in hdiag')
    call deallocp(eval ,id='eval in hdiag')
  endif
  call allocp(tgrad,nred,id='tgrad in hdiag')
  call allocp(eval ,nred,id='eval in hdiag')

  ! project out redundant eigevalues and corresponding eigenvectors
  k = 0
  do i=1,nhtmp
    if(dabs(etmp(i)) > tresh)then
      k = k + 1
      eval(k) = etmp(i)
      do j=1,nhtmp
        evec(j,k) = evec(j,i)
      enddo
    endif
  enddo

  ! calculate transformed gradient: tilde{g} = U^T * g
  tgrad(:) = d0
  do i=1,nred
    do j=1,nhtmp
      tgrad(i) = tgrad(i) + evec(j,i) * gradient(j)
    enddo
  enddo

  ! count number of positive, negative and zero eigenvalues
  nep=0; nem=0; nez=0
  do i=1,nhtmp
    if(dabs(etmp(i)) > tresh)then
      if(etmp(i) < d0)nem = nem + 1
      if(etmp(i) > d0)nep = nep + 1
    else
      nez = nez + 1
    endif
  enddo
  if(nep == 0)Iminmax = .true.  ! we reached the maximum of the function
  if(nep /= 0)Iminmax = .false.

  if(loc%iprloc.ge.1)then
    write(lupri,'(/a,3i4)')' Eigenvalues of the hessian;  +,-,0  =  ',nep,nem,nez
    write(lupri,'(8d12.4)') (etmp(i),i=1,nhtmp)
  endif

  call qexit('hdiag')

  !------------------------------------------------------------------------------!

  end subroutine hdiag

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine get_mu(eval,tgrad,tradius,shift,tpred)
  !******************************************************************************!
  !                                                                              !
  !                   Calculate level-shift parameter mu                         !
  !                                                                              !
  !          tradius = sqrt{ sum{ (tgrad)^2_i / (eval_i-shift)^2  }  }           !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  real*8, intent(in ) :: tradius
  real*8, intent(out) :: shift, tpred

  real*8, intent(in), pointer :: eval(:), tgrad(:)

  ! Internal variables
  integer :: i, j, k, ndim, lupri, imax
  real*8  :: tcurr, xcurr, cutoff, xmax, xmin, step

  !------------------------------------------------------------------------------!
  call qenter('get_mu')

  lupri = var%lupri

  ! Get some dimensions
  ndim = size(eval,1)

  ! initialization
  imax  = 10000
  cutoff= trsh_mu * tradius
  step  = 0.1d0
  xmax  = eval(ndim) + step
  xmin  = eval(ndim)
  tcurr = f(xmax)

  ! find upper boundary if necessary
  i = 0
  if(tcurr > tradius)then
    do
      i = i + 1
      if(i > imax)then
        write(lupri,'(a)')' WARNING: Maximum number of iteration reached in get_mu(1)!!!'
        write(lupri,'(10x,a,d12.4,a,d12.4)')'f(xmax) =',tcurr,' tradius =',tradius
        exit
      endif
      step = xmax - xmin
      xmin = xmax
      xmax = xmax + d2*step
      tcurr = f(xmax)
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

  if(tcurr > tradius) xmin = xcurr
  if(tcurr < tradius) xmax = xcurr

  enddo

  shift = xcurr
  tpred = f(xcurr)

  !------------------------------------------------------------------------------!

  call qexit('get_mu')

  !------------------------------------------------------------------------------!

  contains

    function f(x)
    implicit none
    real*8, intent(in ) :: x
    real*8  :: f
    integer :: ii

    f = d0
    do ii=1,ndim
      if(dabs(eval(ii) - x) < 1.0d-14)then
        write(lupri,'(a)')' WARNING: mu is too close to one of the eigenvalues!!!'
        write(lupri,'(a,i4,a,d12.4)')' WARNING: i=',ii,' mu-e_i =',x-eval(ii)
      endif
      f = f + tgrad(ii) * tgrad(ii) / ((eval(ii) - x) * (eval(ii) - x))
    enddo
    f = dsqrt(f)

    end function f

  end subroutine get_mu

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine updateT(dnrmk,diffg,diffq,tradius,rejected,full,itr)
  !******************************************************************************!
  !                                                                              !
  !                           Update trust radius                                !
  !                                                                              !
  !******************************************************************************!
  implicit none

  ! External variables
  integer, intent(in   ) :: itr
  real*8,  intent(in   ) :: dnrmk, diffg, diffq
  real*8,  intent(inout) :: tradius
  logical, intent(in   ) :: full
  logical, intent(inout) :: rejected

  ! Internal variables
  character*6 :: cond
  integer :: lupri
  real*8, parameter :: d1q = 0.25d0, d3q = 0.75d0, d5q = 1.25d0
  real*8 :: R, tradius_old

  !------------------------------------------------------------------------------!
  call qenter('updateT')

  lupri  = var%lupri
  R = diffg / diffq
  tradius_old = tradius
  cond = '??????'

  if(rejected)then
    cond = 'r     '
    if(R < d0 ) cond = 'r  R<0'
    if(R < d1q) cond = 'r  R<1'
    if(R > d3q) cond = 'r  R>3'
  else
    if(dnrmk > tradius)then
      cond = '   k>T'
    else
      if(diffq > d0)then
        cond = 'qP q>0 '
        if(R < d0 ) cond = 'qP R<0'
        if(R < d1q) cond = 'qP R<1'
        if(R > d5q) cond = 'qP R>5'
      else 
        cond = 'qN    '
        if(R > d0 ) cond = 'qN R>0'
        if((R < d1q).or.(R > d1+d3q)) cond = 'qN R>7'
        if((R > d3q).and.(R < d5q))   cond = 'qN R<5'
      endif
    endif
  endif

  if(rejected)then
    if(diffq < d0)then
      write(lupri,'(a)')' STOP: The value of the quadratic approximation of the functional '
      write(lupri,'(a)')'       decreased although the level-shifted Newton step was used '
      stop
    endif
    rejected = .false.
    if(R < d0 ) rejected = .true.
    if(R < d1q) tradius  = d1q* tradius
    if(R > d3q) tradius  = d2 * tradius
  else
    if(dnrmk > tradius)then
!     if(diffg < d0) rejected = .true.
      rejected = .true.
    else
      if(diffq > d0)then
        if(full.and.(R < d1q)) tradius = dnrmk
        if(R < d0 ) rejected = .true.
        if(R < d1q) tradius  = d1q * tradius
        if(R > d5q) tradius  = max(d2*dnrmk,tradius_old)
      else 
        if(R > d0 ) rejected = .true.
        if(full) tradius = dnrmk
        if((R < d1q).or.(R > d1+d3q)) tradius = d1q * tradius
        if((R > d3q).and.(R < d5q))   tradius = max(d2*dnrmk,tradius_old)
      endif
    endif
  endif

  write(lupri,'(/a,2x,a4,2x, 2x,a,4x, 2x,a,3x, 2x,a,2x, 2x,a,2x, 1x,a, 4x,a, 5x,a, 6x,a)') &
    '  ','ITR','T_old','T_new ','|KAPPA|','  R   ','REJ','CON','DIFFG','DIFFQ'
  write(lupri,'(a,2x,i4,3d11.2,f11.5,a,l,a,2x,a,4d11.2)') &
    '*T',itr,tradius_old,tradius,dnrmk,R,'  ',rejected,' ',cond,diffg,diffq

  call qexit('updateT')

  !------------------------------------------------------------------------------!

  end subroutine updateT

end module localization
