!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!
!     This module contains tools for working with quaternion matrixes
!
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
MODULE QMatrixes

! for allocation/deallocations
use Memory  
! deta type for strings
use mystring

! neede for 64bit machines and -i8 flags
#include "iratdef.h"

! common printing parameter in the module 
integer(kind=4), private :: iprint_qmtx=0

! own data type for IREP value
type, public :: irp
private
 integer (kind=8/IRAT) :: irep 
end type irp

! own data type for NZ value
type, public :: nzp
private
! integer (kind=4) :: nz
 integer (kind=8/IRAT) :: nz
end type nzp
 
! data type - size of matrixes in 2 ferm.irreps - for symmetry packed quaternion matrixes
type, public :: qmtx_pck
private
 integer (kind=8/IRAT) :: nr1,nc1,nr2,nc2
end type

! data type for symmetry unpacked quaternion SQUARE matrixes
type, public :: qmtx_unpck
private
 integer (kind=8/IRAT) :: n1,n2
end type

! universal structure for the quaternion matrix, either symmetry packed, symmetry unpacked
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
!  the total size on the array
  integer (kind=8) :: array_size = 0
  logical ::  simple_mtx=.false.
  logical ::  symm_packed=.false.
! ...important flag to check if matrix is allocated !
  logical ::  is_alloc = .false.
! matrix may be real or integer
  logical ::  is_real = .true.
! short text attached to the matrix
  !character (len=30) :: qmtx_text 
  type(string), dimension(1) :: qmtx_text 
end type QMtx

! common name for contructors
interface MakeQMtx
   module procedure MakeQMtxPck, MakeQMtxUnPck
end interface 

!interface assignment ( = )
!  module procedure assign_qmatrix ! for QMatrix
!end interface 

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function compatible_2qmtx(amat,bmat) result(compatible)
!---------------------------------------------------------------------------
!
!  Logical function returning compatibility confirmation of two-matrixes
!
!------------------------------------------------------------------
type(QMtx), intent(in) :: amat,bmat
logical (kind=1) :: compatible
call QENTER('compatible_2qmtx')
! check matrix compatibilities
    compatible =                          &
                 (amat%is_alloc)   .AND.  &
                 (bmat%is_alloc)   .AND.  & 
                 (amat%nr(1) == bmat%nr(1)) .AND.  &
                 (amat%nr(2) == bmat%nr(2)) .AND.  &
                 (amat%nc(1) == bmat%nc(1)) .AND.  &  
                 (amat%nc(2) == bmat%nc(2)) .AND.  & 
                 (amat%array_size == bmat%array_size)  .AND.  & 
                 (.not.xor(amat%simple_mtx,bmat%simple_mtx))   .AND.  &   
                  (.not.xor(amat%symm_packed,bmat%symm_packed) ) .AND.  &
                 (.not.xor(amat%is_real,bmat%is_real))   .AND.  &   
                 (amat%nz == bmat%nz)
!
call QEXIT('compatible_2qmtx')
end function compatible_2qmtx


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function compatible_2qmtx_part(amat,bmat) result(compatible)
!---------------------------------------------------------------------------
!
!  Logical function returning compatibility of two-matrixes after their comarison
!
! NOTE: the nz may not be identical !
!
!------------------------------------------------------------------
type(QMtx), intent(in) :: amat,bmat
logical (kind=1) :: compatible
call QENTER('compatible_2qmtx_part')
! check matrix compatibilities
    compatible =                          &
                 (amat%is_alloc)   .AND.  &
                 (bmat%is_alloc)   .AND.  & 
                 (amat%nr(1) == bmat%nr(1)) .AND.  &
                 (amat%nr(2) == bmat%nr(2)) .AND.  &
                 (amat%nc(1) == bmat%nc(1)) .AND.  &  
                 (amat%nc(2) == bmat%nc(2)) .AND.  & 
                 (.not.xor(amat%simple_mtx,bmat%simple_mtx))   .AND.  &   
                  (.not.xor(amat%symm_packed,bmat%symm_packed) ) .AND.  &
                 (.not.xor(amat%is_real,bmat%is_real))  
!
call QEXIT('compatible_2qmtx_part')
end function compatible_2qmtx_part

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_isalloc(qmat)
!-----------------------------------------------------------------------
!
! set the is_alloc attribute to .TRUE.   
!
!-----------------------------------------------------------------------
type(QMtx), intent(inout) :: qmat
call QENTER('set_isalloc')
 qmat%is_alloc = .true.
call QEXIT('set_isalloc')
end subroutine set_isalloc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine unset_isalloc(qmat)
!-----------------------------------------------------------------------
!
! set the is_alloc attribute to .FALSE.   
!
!-----------------------------------------------------------------------
type(QMtx), intent(inout) :: qmat
call QENTER('unset_isalloc')
 qmat%is_alloc = .false.
call QEXIT('unset_isalloc')
end subroutine unset_isalloc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_qmtx_kstart(qmat) result (kstart_out)
!----------------------------------------------------------------------
!
!      returns kstart value of entering quaternion matrix "qmtx"
!
!----------------------------------------------------------------------
#include "priunit.h"
type(QMtx), intent(in) :: qmat
integer(kind=4) :: kstart_out
!integer :: kstart_out
call QENTER('get_qmtx_kstart')
  if (.not.qmat%is_alloc) then
    call QUIT('get_qmtx_kstart: can not return "kstart" position &
                when memory for matrix is NOT allocated !')
  else
    kstart_out =  ( qmat%kstart )
  endif
call QEXIT('get_qmtx_kstart')
end function get_qmtx_kstart

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_iprqmtx(iprint)
!=================================================
!            set up print variable
!=================================================
#include "priunit.h"
integer (kind=8/IRAT), intent(in) :: iprint
CALL QENTER('QMatrixes: set_iprqmt')
 iprint_qmtx = iprint
CALL QEXIT('QMatrixes: set_iprqmt')
end subroutine set_iprqmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_iprint_qmtx() result(iprint)
!----------------------------------------------------------------
!
!        get the (protected) printing variable
!
!----------------------------------------------------------------
!integer (kind=8/IRAT) :: iprint
integer (kind=4) :: iprint
CALL QENTER('QMatrixes: get_iprqmtx')
 iprint = iprint_qmtx
CALL QEXIT('QMatrixes: get_iprqmtx')
end function get_iprint_qmtx

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
!-----------------------------------------------
!
!        define setting of NZ
!
!-----------------------------------------------
#include "priunit.h"
type(nzp) :: nz_out
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function set_qmtx_pck (nr1in,nc1in,nr2in,nc2in) result(qmtx_pck_out)
!-------------------------------------------------------------------------
!
!   set dimensions (NFSYM=2) of a symmetry packed quaternion matrix
!
!-------------------------------------------------------------------------
#include "priunit.h"
type(qmtx_pck) :: qmtx_pck_out
integer (kind=8/IRAT), intent(in) :: nr1in,nc1in
integer (kind=8/IRAT), intent(in), optional :: nr2in,nc2in
CALL QENTER('QMatrixes: set_qmtx_pck')
! set dimensions of fermion irep 1
qmtx_pck_out%nr1 = nr1in; qmtx_pck_out%nc1 = nc1in 
! check entering dimensions !
if (present(nr2in).and.present(nc2in)) then
   qmtx_pck_out%nr2 = nr2in
   qmtx_pck_out%nc2 = nc2in
else
! ...  if both dimensions of second irrep are not present, set them to zero
   qmtx_pck_out%nr2 = 0
   qmtx_pck_out%nc2 = 0
endif
CALL QEXIT('QMatrixes: set_qmtx_pck')
end function set_qmtx_pck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE AddTextQMtx(qm,text)
!-------------------------------------------------------------------------
!
!    Attach text(string) attribute to the entering quaternion matrix
! If there is a previous text, replace it; 
!
! if there is no "text" variable in the input, clear it
!
!-------------------------------------------------------------------------
#include "priunit.h"
type(QMtx),intent(inout) :: qm
character(*), intent(in), optional :: text
call QENTER('QMatrixes: AddTextQMtx')

if (len(text).ge.1) then
! if there is some string, attach it, otherwise do nothing
  if (is_string(qm%qmtx_text)) then
    call cancel_string(qm%qmtx_text)
    call make_string(qm%qmtx_text,text)
  else
    call make_string(qm%qmtx_text,text)
  endif
else
! if no string is attached, cancel o
    call cancel_string(qm%qmtx_text)
    call set_no_string(qm%qmtx_text)
endif

call QEXIT('QMatrixes: AddTextQMtx')
END SUBROUTINE 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE set_qmtx_nz(qm,nz_in)
!-------------------------------------------------------------------------
!
!               set nz_in into "qm"
!
!-------------------------------------------------------------------------
type(QMtx),intent(inout) :: qm
integer(kind=8/IRAT),intent(in) :: nz_in
CALL QENTER('set_qmtx_nz')
  qm%nz = nz_in
CALL QEXIT('set_qmtx_nz')
END SUBROUTINE set_qmtx_nz

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE QMtxInfo(qm)
!=========================================================================
!
!    print info (basic variables etc) of the quaternion matrix,
!  plus the matrix alone, if allocated
!
!=========================================================================
#include "priunit.h"
#include "dgroup.h"
type(QMtx),intent(in) :: qm
integer(kind=4) :: ifsym,nt
integer (kind=8) :: ioff
!real (kind=8), pointer, dimension(:) :: my_work_p
CALL QENTER('QMtxInfo')

CALL HEADER('Output from QMtxInfo',-1)

!call get_mem_pointer(my_work_p)

if (is_string(qm%qmtx_text)) then
   !write(LUPRI,'(4x,a,a)') 'Attached text 1:',get_string(qm%qmtx_text)
   !write(LUPRI,*) 'Attached text 1:',get_string(qm%qmtx_text)
  !  write(LUPRI,*) get_string(qm%qmtx_text)
   ! print *,get_string(qm%qmtx_text)
   write(LUPRI,'(4x,a,$)') 'Attached text:'
   call print_string(qm%qmtx_text)
endif

write(LUPRI,'(4x,a,i1,1x,a,i1)') 'IREP=',qm%irep,'NZ=',qm%nz
write(LUPRI,'(2x,a,i5,a,i5)') 'start position on the big array:',qm%kstart,   & 
    ' size of the allocated array:',qm%array_size

if (.not.qm%simple_mtx.and.qm%symm_packed) then
 write(LUPRI,*) 'fermion symmetry packed matrix '
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
  write(LUPRI,*) 'full (symmetry unpacked) quaternion matrix '
  if (qm%is_alloc) then
! ... print TOTAL QMatrix
    nt = ( qm%nr(1) + qm%nr(2) )
    CALL PRQMAT(PW(qm%kstart),nt,nt,nt,nt,qm%nz,IPQTOQ(1,qm%irep),LUPRI) 
  else
    write(LUPRI,*) 'no printing out - this matrix is not allocated'
  endif

else if (qm%simple_mtx) then
! ... print out of a trivial quaternion matrix
   write(LUPRI,*) 'triavial quaternion matrix'
    call quit('QMtxInfo: needs to be implemented ')
else
 call quit('QMtxInfo: blind "if" branch ')
endif

 call FLSHFO(LUPRI)

CALL QEXIT('QMtxInfo')
END SUBROUTINE QMtxInfo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION MakeQMtxPck(qm_pk_in, irp_in, nz_in, text) result(qm)
!-------------------------------------------------------------------
!
!   basic constructor for symmetry packed quaternion matrix
!
!-------------------------------------------------------------------
#include "priunit.h"
#include "dgroup.h"
type(QMtx) :: qm
type(qmtx_pck), intent(in) :: qm_pk_in
type(irp), optional, intent(in) ::  irp_in
type(nzp), optional, intent(in) ::  nz_in
character(*), optional, intent(in) :: text
!integer (kind=8) :: isize
integer (kind=8/IRAT) :: isize
integer(kind=4) :: good_alloc
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

if (present(text)) then
   call make_string(qm%qmtx_text, text)
else
   call set_no_string(qm%qmtx_text)
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
! store the size
qm%array_size = isize

! allocate array for the matrix
!CALL MEMGET('REAL',qm%kstart,isize,PW,KFREE,LFREE)
CALL mymemget('REAL',qm%kstart,isize)
qm%is_real  = .true.
qm%is_alloc = .true. ! array is allocated

!if (iprint_qmtx >=5 ) then
! call print_mem_param()
!endif

CALL QEXIT('QMatrixes: MakeQMtxPck')
END FUNCTION MakeQMtxPck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION MakeQMtxUnPck(qm_upk_in, irp_in, nz_in, text) result(qm)
!---------------------------------------------------------------------
!
!             Make quaternion full (unpacked) matrix.
!
!   On input: 
!            qm_upk_in  - structure containing dimensions of matrixe
!            irp_in, nz_in - structures with IREP and NZ 
!            text - optional text attached to the matrix
!
!---------------------------------------------------------------------
#include "priunit.h"
#include "dgroup.h"
type(QMtx) :: qm
type(irp), optional, intent(in) ::  irp_in
type(nzp), optional, intent(in) ::  nz_in
type(qmtx_unpck),  intent(in) :: qm_upk_in
character(*), optional, intent(in) :: text
!integer (kind=8) :: isize
integer (kind=8/IRAT) :: isize
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

if (present(text)) then
  call make_string(qm%qmtx_text,text)
else
  call set_no_string(qm%qmtx_text)
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
! store the size
qm%array_size = isize
!if (iprint_qmtx >=5 ) then
! call print_mem_param()
!endif

! allocate array for the matrix
CALL mymemget('REAL',qm%kstart,isize)
!if (iprint_qmtx >=5 ) then
! call print_mem_param()
!endif
qm%is_alloc = .true.
qm%is_real  = .true.
CALL QEXIT('QMatrixes: MakeQMtxUnPck')
END FUNCTION MakeQMtxUnPck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE DeAllocQMtx(qm)
!-----------------------------------------------------------------------------
!
!   Deallocate quaternion matrix; 
!
!-----------------------------------------------------------------------------
#include "priunit.h"
type(QMtx) :: qm
call QENTER('QMatrixes: DeAllocQMtx')

!call MEMREL('QMatrixes: DestroyQMtx',PW,qm%kstart,qm%kstart,KFREE,LFREE)
call mymemrel('QMatrixes: DestroyQMtx',qm%kstart,qm%kstart)
qm%is_alloc = .false.
qm%array_size = 0
call cancel_string(qm%qmtx_text) ! cancel the string text
! shall I zero all attribites ?
call QEXIT('QMatrixes: DeAllocQMtx')
END SUBROUTINE DeAllocQMtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function getdiff_qmtx(amat,bmat) result(diff)
!-------------------------------------------------------------------------
!
! returns the average difference of two quatern.matrixes
!
!-------------------------------------------------------------------------
type(QMtx), intent(in) :: amat,bmat
real(kind=8) :: diff
logical(kind=1) :: compatible
integer(kind=4) :: i
call QENTER('getdiff_qmtx')
! check matrix compatibilities
compatible = compatible_2qmtx(amat,bmat)

if (.not.compatible) then
  call print_qmtx_attrib(amat)
  call print_qmtx_attrib(bmat)
  call QUIT('getdiff_qmtx: matrixes are not compatible !')
endif

diff = 0.D0
do i=0,amat%array_size -1
 diff = abs(PW(amat%kstart+i) - PW(bmat%kstart+i)) + diff
enddo
diff = diff / dfloat(amat%array_size)

call QEXIT('getdiff_qmtx')
end function getdiff_qmtx


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function getdiff_part_qmtx(amat,bmat) result(diff)
!-------------------------------------------------------------------------
!
! returns the average difference of two quatern.matrixes
!
!-------------------------------------------------------------------------
type(QMtx), intent(in) :: amat,bmat
real(kind=8) :: diff
logical(kind=1) :: compatible
integer(kind=8) :: size_min
integer(kind=4) :: i
call QENTER('getdiff_part_qmtx')
! check matrix compatibilities
compatible = compatible_2qmtx_part(amat,bmat)

if (.not.compatible) then
  call print_qmtx_attrib(amat)
  call print_qmtx_attrib(bmat)
  call QUIT('getdiff_qmtx: matrixes are not compatible !')
endif

size_min = min(amat%array_size,bmat%array_size)

diff = 0.D0
do i=0, size_min -1
 diff = abs(PW(amat%kstart+i) - PW(bmat%kstart+i)) + diff
enddo
diff = diff / dfloat(amat%array_size)

call QEXIT('getdiff_part_qmtx')
end function getdiff_part_qmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine substract_qmtx(amat,bmat,fac_a,fac_b)
!--------------------------------------------------------------------
!
!  Substract bmat from amat
!
!   DO: amat := fac_a*amat - fac_b*bmat;
!
!  returns result in "amat"
!
!--------------------------------------------------------------------
type(QMtx), intent(inout) :: amat
type(QMtx), intent(in) :: bmat
real(kind=8),intent(in),optional :: fac_a,fac_b
real(kind=8) :: temp
logical(kind=1) :: compatible
integer(kind=4) :: i
call QENTER('substract_qmtx')
! check matrix compatibilities
compatible = compatible_2qmtx(amat,bmat)

if (.not.compatible) then
 call print_qmtx_attrib(amat)
 call print_qmtx_attrib(bmat)
 call QUIT('getdiff_qmtx: matrixes are not compatible !')
endif

!
! todo: replace later by DAXPY  !!!
!
if (present(fac_a).and.present(fac_b)) then

do i=0,amat%array_size -1
 temp = PW(amat%kstart+i)
 PW(amat%kstart+i)  =  fac_a*temp - fac_b*PW(bmat%kstart+i)
enddo

else if (present(fac_a).and..not.present(fac_b)) then

do i = 0, amat%array_size -1
 temp = PW(amat%kstart+i)
 PW(amat%kstart+i)  =  fac_a*temp - PW(bmat%kstart+i)
enddo

else if (.not.present(fac_a).and..not.present(fac_b)) then

do i=0,amat%array_size -1
 temp = PW(amat%kstart+i)
 PW(amat%kstart+i)  =  temp - PW(bmat%kstart+i)
enddo

else
 call QUIT('substract_qmtx: dead branch !!!')
endif

call QEXIT('substract_qmtx')
end subroutine substract_qmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine substract_part_qmtx(amat,bmat,fac_a,fac_b)
!--------------------------------------------------------------------
!
!    DO: amat := fac_a*amat - fac_b*bmat;
!  returns results in "amat"
!
!  Substract partially bmat from amat :  find the lowest NZ first 
!  to set the size
!
!--------------------------------------------------------------------
type(QMtx), intent(inout) :: amat
type(QMtx), intent(in) :: bmat
real(kind=8),intent(in),optional :: fac_a,fac_b
real(kind=8) :: temp
logical(kind=1) :: compatible
integer(kind=4) :: i
integer(kind=8) :: size_min
call QENTER('substract_part_qmtx')
! check matrix compatibilities
compatible = compatible_2qmtx_part(amat,bmat)

if (.not.compatible) then
 call print_qmtx_attrib(amat)
 call print_qmtx_attrib(bmat)
 call QUIT('getdiff_part_qmtx: matrixes are not compatible !')
endif

size_min = min(amat%array_size,bmat%array_size)

!
! todo: replace later by DAXPY  !!!
!
if (present(fac_a).and.present(fac_b)) then

do i=0,size_min -1
 temp = PW(amat%kstart+i)
 PW(amat%kstart+i)  =  fac_a*temp - fac_b*PW(bmat%kstart+i)
enddo

else if (present(fac_a).and..not.present(fac_b)) then

do i = 0, size_min -1
 temp = PW(amat%kstart+i)
 PW(amat%kstart+i)  =  fac_a*temp - PW(bmat%kstart+i)
enddo

else if (.not.present(fac_a).and..not.present(fac_b)) then

do i=0,size_min -1
 temp = PW(amat%kstart+i)
 PW(amat%kstart+i)  =  temp - PW(bmat%kstart+i)
enddo

else
 call QUIT('substract_qmtx: dead branch !!!')
endif

call QEXIT('substract_part_qmtx')
end subroutine substract_part_qmtx


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_qmtx_attrib(amat)
!-------------------------------------------------------------------
!
!   print out only attributes of the matrix data type (not whole matrixes)
!
!-------------------------------------------------------------------
#include "priunit.h"
type(QMtx), intent(in) :: amat
call QENTER('print_qmtx_attrib')
 write(LUPRI,'(/2x,a)')  & 
 '---------- Basic attributes of the QMtx matrix data type -----------'
 write(LUPRI,'(a,l1,$)') 'is_alloc:',amat%is_alloc
 write(LUPRI,'(a,l1,$)') '  simple_mtx:',amat%simple_mtx
 write(LUPRI,'(a,l1,$)') '  symmetry packed:',amat%symm_packed
 write(LUPRI,'(a,l1)') '  is real:',amat%is_real
 if (amat%is_alloc) then
  write(LUPRI,'(a,i7)') 'size of array:',amat%array_size
 endif
 write(LUPRI,'(a,2i3,a,2i3,a,i1,a,i2)') 'nr(1),nc(1);nr(2),nc(2:', &
       amat%nr(1),amat%nc(1),';',amat%nr(2),amat%nc(2),  & 
       ' --- NZ:',amat%nz,' IREP:',amat%irep
 if (is_string(amat%qmtx_text)) then
  write(LUPRI,'(4x,a,$)') '...attached text:'
  call print_string(amat%qmtx_text)
 endif
 call FLSHFO(LUPRI)
call QEXIT('print_qmtx_attrib')
end subroutine print_qmtx_attrib
END MODULE QMatrixes

