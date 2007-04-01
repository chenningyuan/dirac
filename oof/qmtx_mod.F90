!
!
! This module contains tools for working with quaternion matrixes
!
!
MODULE QMatrixes

! for allocation/deallocations
use Memory  

#include "iratdef.h"

!     print parameter
integer(kind=4), private :: iprqmtx=0

type, public :: irp
private
 integer (kind=8/IRAT) :: irep 
end type irp

type, public :: nzp
private
! integer (kind=4) :: nz
 integer (kind=8/IRAT) :: nz
end type nzp
 
! data type for symmetry packed quaternion matrixes
type, public :: qmtx_pck
private
 integer (kind=8/IRAT) :: nr1,nc1,nr2,nc2
end type

! data type for symmetry unpacked quaternion SQUARE matrixes
type, public :: qmtx_unpck
private
 integer (kind=8/IRAT) :: n1,n2
end type

!
! universal structure for quaternion matrix, either symmetry packed, symmetry unpacked
! or simple unsymmetric
!
type QMtx
private
! symm.distingiushed dimensions (1:NFSYM)
!  integer (kind=4), dimension(2)  :: nr, nc
  integer (kind=8/IRAT), dimension(2)  :: nr, nc
!       ...  IREP and NZ
!  integer (kind=4) :: irep,nz
  integer (kind=8/IRAT) :: irep,nz
  integer (kind=4) :: kstart ! contains position on the PW (WORK) array !
! ... offsets for symmetry distinguished matrixes
  integer (kind=8) :: ioff_pk = 0
  integer (kind=8), dimension(2,2) :: ioff_upk = 0
!  
  logical ::  simple_mtx=.false.
  logical ::  symm_packed=.false.
  logical ::  is_alloc = .false.
end type QMtx

! common name for contructors
interface MakeQMtx
   module procedure MakeQMtxPck, MakeQMtxUnPck
end interface 

!interface assignment ( = )
!  module procedure assign_qmatrix ! for QMatrix
!end interface 

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_iprqmtx(iprint)
!=================================================
!            set up print variable
!=================================================
#include "priunit.h"
integer (kind=8/IRAT), intent(in) :: iprint
CALL QENTER('QMatrixes: set_iprqmt')
  iprqmtx = iprint
!  write(LUPRI,*) 'set_iprqmtx: set iprint=',iprqmtx
CALL QEXIT('QMatrixes: set_iprqmt')
end subroutine set_iprqmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function set_irp (irep_in) result (irep_out)
!================================================================
!        define setting of irp
!================================================================
#include "priunit.h"
type(irp) :: irep_out
!integer (kind=4), intent(in) :: irep_in
integer (kind=8/IRAT), intent(in) :: irep_in
CALL QENTER('QMatrixes: set_irp')
if (.NOT.(irep_in .ge. 0 .and. irep_in .le. 8)) then
 write(lupri,*) 'QMatrixes: set_irp, entering irep_in=',irep_in
 call quit('QMatrixes: set_irp - wrong value of irep !!!')
endif
irep_out%irep = irep_in
CALL QEXIT('QMatrixes: set_irp')
end function set_irp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function set_nz (nz_in) result (nz_out)
!
! define setting of NZ
!
#include "priunit.h"
type(nzp) :: nz_out
!integer (kind=4), intent(in) :: nz_in
integer (kind=8/IRAT), intent(in) :: nz_in
CALL QENTER('QMatrixes: set_nz')
  if (.NOT.(nz_in .eq. 1 .or. nz_in .eq. 2  .or. nz_in .eq. 4)) then
    write(lupri,*) 'QMatrixes: set_nz, entering nz_in=',nz_in
    call quit('QMatrixes: set_nz - wrong value of nz_in (no 1/2/4) !!!')
  endif
  ! nz_out%nz = nz_in
  nz_out = nzp(nz_in)
CALL QEXIT('QMatrixes: set_nz')
end function set_nz

!=======================================================================
function set_qmtx_pck (nr1in,nc1in,nr2in,nc2in) result(qmtx_pck_out)
!
! set dimensions (NFSYM=2) of a symmetry packed quaternion matrix
!
#include "priunit.h"
type(qmtx_pck) :: qmtx_pck_out
integer (kind=8/IRAT), intent(in) :: nr1in,nc1in
integer (kind=8/IRAT), intent(in), optional :: nr2in,nc2in
CALL QENTER('QMatrixes: set_qmtx_pck')
qmtx_pck_out%nr1 = nr1in; qmtx_pck_out%nc1 = nc1in 
! check entering dimensions !
if (present(nr2in).and.present(nc2in)) then
   qmtx_pck_out%nr2 = nr2in
   qmtx_pck_out%nc2 = nc2in
else
   qmtx_pck_out%nr2 = 0
   qmtx_pck_out%nc2 = 0
 ! ... must be taken care off  
endif
CALL QEXIT('QMatrixes: set_qmtx_pck')
end function set_qmtx_pck

function set_qmtx_unpck (n1in,n2in) result(qmtx_unpck_out)
!==============================================================================
!
!   set dimensions (NFSYM=2) of a symmetry unpacked SQUARE quaternion matrix
!
!==============================================================================
#include "priunit.h"
type(qmtx_unpck) :: qmtx_unpck_out
integer (kind=8/IRAT), intent(in) :: n1in
integer (kind=8/IRAT), intent(in), optional :: n2in
CALL QENTER('QMatrixes: set_qmtx_unpck')
qmtx_unpck_out%n1 = n1in
! check entering dimensions !
if (present(n2in)) then
   qmtx_unpck_out%n2 = n2in
else
   qmtx_unpck_out%n2 = 0
endif
CALL QEXIT('QMatrixes: set_qmtx_unpck')
end function set_qmtx_unpck

!======================================================
!         print out the quaternion matrix
!======================================================
SUBROUTINE QMtxInfo(qm)
#include "priunit.h"
#include "dgroup.h"
type(QMtx),intent(in) :: qm
integer(kind=4) :: ifsym,nt
integer (kind=8) :: ioff
CALL QENTER('QMtxInfo')
write(LUPRI,'(/,4x,a)') 'QMtx: info '

if (.not.qm%simple_mtx.and.qm%symm_packed) then
 write(LUPRI,*) 'fermion symmetry packed matrix !!!'
 !ioff = 1
 ioff = qm%kstart
 do ifsym = 1, NFSYM
   write(LUPRI,*) 'fermion symmetry ',ifsym,'/',NFSYM,  &
  'nrows,ncols,nz:',qm%nr(ifsym),qm%nc(ifsym),qm%nz
   write(LUPRI,*) 'offset=',qm%ioff_pk
  if (qm%is_alloc) then
! ... print QMatrix
   CALL PRQMAT(PW(ioff),qm%nr(ifsym),qm%nc(ifsym),qm%nr(ifsym),qm%nc(ifsym),qm%nz,  &
              IPQTOQ(1,qm%irep),LUPRI)
  else
    write(LUPRI,*) 'matrix is not allocated,no print out'
  endif
  ioff = ioff + qm%ioff_pk
 enddo

else if (.not.qm%simple_mtx.and..not.qm%symm_packed) then
! ... symmetry unpacked (full) quaternion matrix
  if (qm%is_alloc) then
! ... print TOTAL QMatrix
    nt = ( qm%nr(1) + qm%nr(2) )
    CALL PRQMAT(PW(qm%kstart),nt,nt,nt,nt,qm%nz,IPQTOQ(1,qm%irep),LUPRI) 
  else
    write(LUPRI,*) 'this matrix is not allocated, so no print out'
  endif

else if (qm%simple_mtx) then
! ... print out of a trivial quaternion matrix

else
 call quit('QMtxInfo: blind "if" branch ')
endif

 write(LUPRI,*) 'qm%kstart=',qm%kstart

CALL QEXIT('QMtxInfo')
END SUBROUTINE QMtxInfo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION MakeQMtxPck(qm_pk_in, irp_in, nz_in) result(qm)
!-------------------------------------------------------------------
!
!    basic constructor for symmetry packed quaternion matrix
!
!-------------------------------------------------------------------
#include "priunit.h"
#include "dgroup.h"
#include "iratdef.h"
type(QMtx) :: qm
type(qmtx_pck), intent(in) :: qm_pk_in
type(irp), optional, intent(in) ::  irp_in
type(nzp), optional, intent(in) ::  nz_in
integer (kind=8) :: isize
CALL QENTER('QMatrixes: MakeQMtxPck')

if (present(irp_in)) then
  qm%irep = irp_in%irep
else
  qm%irep = 0 
endif

if (present(nz_in)) then
 qm%nz = nz_in%nz
else
 qm%nz = NZ ! from dgroup.h
endif

qm%nr(1) = qm_pk_in%nr1
qm%nc(1) = qm_pk_in%nc1
qm%nr(2) = qm_pk_in%nr2
qm%nc(2) = qm_pk_in%nc2
if (NFSYM == 2) then
  qm%ioff_pk = qm%nr(1) * qm%nc(1) * qm%nz
else 
  qm%ioff_pk = 0
endif
qm%simple_mtx  = .false.
qm%symm_packed = .true.
!  array size of the symmetry packed matrix
isize = qm%nr(1)*qm%nc(1)*qm%nz + qm%nr(2)*qm%nc(2)*qm%nz

! check the array
!CALL MEMCHK('QMatrixes: MakeQMtx',PW,1)
if (iprqmtx >=5 ) then
 call print_mem_param()
! write(lupri,*) 'KFREE,LFREE=',KFREE,LFREE
 write(lupri,*) 'isize=',isize
endif

! allocate array for the matrix
!CALL MEMGET('REAL',qm%kstart,isize,PW,KFREE,LFREE)
CALL mymemget('REAL',qm%kstart,isize)

qm%is_alloc = .true. ! array is allocated

CALL QEXIT('QMatrixes: MakeQMtxPck')
END FUNCTION MakeQMtxPck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION MakeQMtxUnPck(qm_upk_in, irp_in, nz_in) result(qm)
!---------------------------------------------------------------------
!
!  Make quaternion full (unpacked) matrix.
!
!---------------------------------------------------------------------
#include "priunit.h"
#include "dgroup.h"
#include "iratdef.h"
type(QMtx) :: qm
type(irp), optional, intent(in) ::  irp_in
type(nzp), optional, intent(in) ::  nz_in
type(qmtx_unpck),  intent(in) :: qm_upk_in
integer (kind=8) :: isize
logical, save :: was_called = .false.
CALL QENTER('QMatrixes: MakeQMtxUnPck')

if (present(irp_in)) then
  qm%irep = irp_in%irep
else
  qm%irep = 0 
endif

if (present(nz_in)) then
 qm%nz = nz_in%nz
else
  qm%nz = NZ ! from dgroup.h
endif

qm%nr(1) = qm_upk_in%n1
qm%nc(1) = qm_upk_in%n1
qm%nr(2) = qm_upk_in%n2
qm%nc(2) = qm_upk_in%n2
qm%ioff_upk(1,1) = 0
if (NFSYM == 2) then
  qm%ioff_upk(1,2) = qm_upk_in%n1
  qm%ioff_upk(2,1) = (qm_upk_in%n1+qm_upk_in%n2)*qm_upk_in%n2
  qm%ioff_upk(2,2) = (qm_upk_in%n1+qm_upk_in%n2)*qm_upk_in%n2 + qm_upk_in%n1
else
  qm%ioff_upk(1,2) = 0; qm%ioff_upk(2,1) = 0; qm%ioff_upk(2,2) = 0 
endif
qm%simple_mtx  = .false.
qm%symm_packed = .false.

! size of array to be allocated
isize = (qm%nr(1)+qm%nr(2))*(qm%nr(1)+qm%nr(2))*qm%nz

if (iprqmtx >=5 ) then
 call print_mem_param()
! write(lupri,*) 'KFREE,LFREE=',KFREE,LFREE
 write(lupri,*) 'isize=',isize
endif

! allocate array for the matrix
!CALL MEMGET('REAL',qm%kstart,isize,PW,KFREE,LFREE)
 CALL mymemget('REAL',qm%kstart,isize)

qm%is_alloc = .true.

CALL QEXIT('QMatrixes: MakeQMtxUnPck')
END FUNCTION MakeQMtxUnPck

SUBROUTINE DeAllocQMtx(qm)
#include "priunit.h"
type(QMtx) :: qm
call QENTER('QMatrixes: DeAllocQMtx')

!call MEMREL('QMatrixes: DestroyQMtx',PW,qm%kstart,qm%kstart,KFREE,LFREE)
call mymemrel('QMatrixes: DestroyQMtx',qm%kstart,qm%kstart)

qm%is_alloc = .false.

call QEXIT('QMatrixes: DeAllocQMtx')
END SUBROUTINE DeAllocQMtx


END MODULE QMatrixes
