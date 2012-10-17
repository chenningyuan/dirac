!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!
!       This module contains tools for working with quaternion matrixes
!
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
MODULE QMatrixes

! for allocation/deallocations
!use Memory  
! deta type for strings
use mystring

! neede for 64bit machines and -i8 flags
#include "iratdef.h"

! common printing parameter in the module 
integer, private :: iprint_qmtx=0

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

!----------------------------------------------------------------------------------------------
!
!   general data type for the quaternion matrix, either symmetry packed, symmetry unpacked
! or simple unsymmetric
!
!----------------------------------------------------------------------------------------------
type QMtx
private
!public

! symm.distingiushed dimensions (1:NFSYM)
!  integer (kind=4), dimension(2)  :: nr, nc
 ! integer (kind=8/IRAT), dimension(2)  :: nr, nc
  integer, dimension(2)  :: nr, nc
!       ...  IREP and NZ
!  integer (kind=4) :: irep,nz


  integer (kind=8/IRAT) :: irep,nz
 ! integer (kind=4) :: kstart ! contains position on the PW (WORK) array !

  integer  :: kstart ! contains position on the PW (WORK) array !

! ... offsets for symmetry distinguished matrixes
  integer (kind=8) :: ioff_pk = 0

  integer (kind=8), dimension(2,2) :: ioff_upk = 0

!  the total size on the array

!  integer (kind=8) :: array_size = 0
  integer :: array_size = 0

  logical ::  simple_mtx=.false.

  logical ::  symm_packed=.false.

! ...important flag to check if matrix is allocated !
  logical ::  is_alloc = .false.

! matrix may be real or integer
  logical ::  is_real = .true.

! short text attached to the matrix
  type(string), dimension(1) :: qmtx_text 

! ... pointers
  real(kind=8), dimension(:,:,:), pointer :: p_rA3D1 => NULL(), p_rA3D2 => NULL()
  !integer,      dimension(:,:,:), pointer :: p_iA3D => NULL()
  !logical,      dimension(:,:,:), pointer :: p_lA3D => NULL()
end type QMtx

! common name for contructors
interface MakeQMtx
   module procedure MakeQMtxPck
   module procedure MakeQMtxUnPck
end interface 

!interface assignment ( = )
!  module procedure assign_qmatrix ! for QMatrix
!end interface 

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function compatible_2qmtx(amat,bmat) result(compatible)
!---------------------------------------------------------------------------
!
!    Logical function returning compatibility confirmation 
!  of two quaternion matrixes
!
!  If matrixes are fully identical, returns .true.
!
!---------------------------------------------------------------------------
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
                      compatible_2qmtx_logical(amat,bmat) .AND. & 
              !  (.not.xor(amat%simple_mtx,bmat%simple_mtx))   .AND.  &   
              !  (.not.xor(amat%symm_packed,bmat%symm_packed) ) .AND.  &
              !  (.not.xor(amat%is_real,bmat%is_real))   .AND.  &   
                 (amat%nz == bmat%nz)
!
call QEXIT('compatible_2qmtx')
end function compatible_2qmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function compatible_2qmtx_logical(amat,bmat) result (log_comp)
!---------------------------------------------------------------------------
!
! .... instead of XOR I have to use new function - as g95 has no XOR (only IEOR)
!
!---------------------------------------------------------------------------
type(QMtx), intent(in) :: amat,bmat
logical (kind=1) :: log_comp
logical (kind=1) :: simple_mtx2, symm_packed2, is_real2
call QENTER('compatible_2qmtx_logical')

           ! ......... doing this: xor in NOT in g95 !!!
              !  (.not.xor(amat%simple_mtx,bmat%simple_mtx))   .AND.  &   
              !  (.not.xor(amat%symm_packed,bmat%symm_packed) ) .AND.  &
              !  (.not.xor(amat%is_real,bmat%is_real))   .AND.  &   
  
     simple_mtx2 = (amat%simple_mtx.AND.bmat%simple_mtx) .OR. &
                   (.NOT.amat%simple_mtx.AND..NOT.bmat%simple_mtx)

     symm_packed2 = (amat%symm_packed.AND.bmat%symm_packed) .OR. &
                    (.NOT.amat%symm_packed.AND..NOT.bmat%symm_packed)

     is_real2    = (amat%is_real.AND.bmat%is_real) .OR.   &
                   (.NOT.amat%is_real.AND..NOT.bmat%is_real) 

  log_comp = simple_mtx2 .AND. symm_packed2 .AND. is_real2 ! resulting

call QEXIT('compatible_2qmtx_logical')
end function compatible_2qmtx_logical


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
                      compatible_2qmtx_logical(amat,bmat) 
          !       (.not.xor(amat%simple_mtx,bmat%simple_mtx))   .AND.  &   
          !       (.not.xor(amat%symm_packed,bmat%symm_packed) ) .AND.  &
          !       (.not.xor(amat%is_real,bmat%is_real))  
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
 qmat%is_alloc = .true.
end subroutine set_isalloc

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine unset_isalloc(qmat)
!-----------------------------------------------------------------------
!
! set the is_alloc attribute to .FALSE.   
!
!-----------------------------------------------------------------------
type(QMtx), intent(inout) :: qmat
 qmat%is_alloc = .false.
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
!integer(kind=4) :: kstart_out ! returning integer*4 !!!
integer :: kstart_out ! returning integer*4 !!!
  if (.not.qmat%is_alloc) then
    call QUIT('get_qmtx_kstart: can not return "kstart" position' & 
                //'when memory for matrix is NOT allocated !')
  else
    kstart_out =  ( qmat%kstart )
  endif
end function get_qmtx_kstart

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_iprqmtx(iprint)
!=================================================
!            set up print variable
!=================================================
#include "priunit.h"
integer, intent(in) :: iprint
iprint_qmtx = iprint
if (iprint_qmtx >= 6) then
  write(LUPRI,'(2x,a,i4)')   &  
  'QMatrixes/set_iprqmt - print level "iprint_qmtx" set to',iprint_qmtx
endif
end subroutine set_iprqmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_iprint_qmtx() result(iprint)
!----------------------------------------------------------------
!
!           get the (protected) printing variable
!
!----------------------------------------------------------------
integer  :: iprint
 iprint = iprint_qmtx
end function get_iprint_qmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function set_irp (irep_in) result (irep_out)
!================================================================
!        define setting of irp
!================================================================
#include "priunit.h"
type(irp) :: irep_out
integer, intent(in) :: irep_in
if (.NOT.(irep_in .ge. 0 .and. irep_in .le. 8)) then
 write(lupri,*) 'QMatrixes: set_irp, entering irep_in=',irep_in
 call quit('QMatrixes: set_irp - wrong value of irep !!!')
endif
irep_out%irep = irep_in
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
!integer (kind=8/IRAT), intent(in) :: nz_in
integer, intent(in) :: nz_in
  if (.NOT.(nz_in .eq. 1 .or. nz_in .eq. 2  .or. nz_in .eq. 4)) then
    write(lupri,*) 'QMatrixes: set_nz, entering nz_in=',nz_in
    call quit('QMatrixes: set_nz - wrong value of nz_in (no 1/2/4) !!!')
  endif
  ! nz_out%nz = nz_in
  nz_out = nzp(nz_in)
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
integer,  intent(in) :: nr1in, nc1in
integer, intent(in), optional :: nr2in,nc2in
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
!integer (kind=8/IRAT), intent(in) :: n1in
integer, intent(in) :: n1in
!integer (kind=8/IRAT), intent(in), optional :: n2in
integer, intent(in), optional :: n2in
qmtx_unpck_out%n1 = n1in
! check entering dimensions !
if (present(n2in)) then
   qmtx_unpck_out%n2 = n2in
else
   qmtx_unpck_out%n2 = 0
endif
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

if (len(text).ge.1) then
! if there is some string, attach it, otherwise do nothing
  if (is_string(qm%qmtx_text)) then
    call cancel_string(qm%qmtx_text)
    call make_string(qm%qmtx_text,text)
  else
    call make_string(qm%qmtx_text,text)
  endif
else
! if no string is attached, cancel string ...
    call cancel_string(qm%qmtx_text)
    call set_no_string(qm%qmtx_text)
endif

END SUBROUTINE 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE set_qmtx_nz(qm,nz_in)
!-------------------------------------------------------------------------
!
!               set nz_in into "qm"
!
!-------------------------------------------------------------------------
type(QMtx),intent(inout) :: qm
!integer(kind=8/IRAT),intent(in) :: nz_in
integer,intent(in) :: nz_in
   qm%nz = nz_in
END SUBROUTINE set_qmtx_nz

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE QMtxInfo(qm, print_qmtx)
!=========================================================================
!
!    print info (basic variables etc) of the quaternion matrix,
!  plus the matrix alone, if allocated
!
!=========================================================================
use Memory
#include "priunit.h"
#include "dgroup.h"
type(QMtx),intent(in) :: qm
logical, optional, intent(in) :: print_qmtx
logical :: doprint = .true.
!integer(kind=4) :: ifsym, ntr, ntc
integer :: ifsym, ntr, ntc
!integer (kind=8) :: ioff
integer :: ioff
!real (kind=8), pointer, dimension(:) :: my_work_p

CALL HEADER('Output from QMtxInfo',-1)

write(LUPRI,'(2x,a,i3)') 'print out parameter, iprint_qmtx:',iprint_qmtx

if (is_string(qm%qmtx_text)) then
   !write(LUPRI,'(4x,a,a)') 'Attached text 1:',get_string(qm%qmtx_text)
   !write(LUPRI,*) 'Attached text 1:',get_string(qm%qmtx_text)
  !  write(LUPRI,*) get_string(qm%qmtx_text)
   ! print *,get_string(qm%qmtx_text)
   write(LUPRI,'(4x,a,$)') 'Attached text to QMtx object:'
   call print_string(qm%qmtx_text)
endif

write(LUPRI,'(4x,a,i1,1x,a,i1)') 'IREP=',qm%irep,'NZ=',qm%nz
call FLSHFO(LUPRI)

if (present(print_qmtx)) doprint = print_qmtx

!
! ... print out only when the matrix is allocated
!
if (qm%is_alloc) then

write(LUPRI,'(2x,a,i5,a,i5)') 'start position on the big array:',qm%kstart,   & 
    ' size of the allocated array:',qm%array_size
call FLSHFO(LUPRI)

if (.not.qm%simple_mtx.and.qm%symm_packed) then

 write(LUPRI,'(2x,a)') 'fermion symmetry packed matrix '
 !ioff = 1
 ioff = qm%kstart
 do ifsym = 1, NFSYM
   write(LUPRI,*) 'fermion symmetry ',ifsym,'/',NFSYM,  &
  'nrows, ncols, nz:',qm%nr(ifsym),qm%nc(ifsym),qm%nz
   write(LUPRI,*) 'offset=',qm%ioff_pk
!  ... print QMatrix
  ! CALL PRQMAT(PW(ioff),qm%nr(ifsym),qm%nc(ifsym),qm%nr(ifsym),qm%nc(ifsym),qm%nz,  &
   if (ifsym == 1) then
     if (.not.associated(qm%p_rA3D1)) then
        call QUIT('QMtxInfo: qm%p_rA3D1 pointer is NOT associated !')
     endif
     
     if (doprint) then
      CALL PRQMAT(qm%p_rA3D1,qm%nr(ifsym),qm%nc(ifsym),qm%nr(ifsym),qm%nc(ifsym),qm%nz,  &
                   IPQTOQ(1,qm%irep),LUPRI )
     endif

   else if (ifsym == 2) then
     if (.not.associated(qm%p_rA3D2)) then
        call QUIT('QMtxInfo: qm%p_rA3D2 pointer is NOT associated !')
     endif
     if (doprint) then
       CALL PRQMAT(qm%p_rA3D2,qm%nr(ifsym),qm%nc(ifsym),qm%nr(ifsym),qm%nc(ifsym),qm%nz,  &
                   IPQTOQ(1,qm%irep),LUPRI )
     endif
   else
     call QUIT('QMtxInfo: ifsym <> 1 or 2 !')
   endif
  ! ... offset for the symmetry packed matrix - cycle over ifsym/NFSYM is needed !
  ioff = ioff + qm%ioff_pk
 enddo

else if (.not.qm%simple_mtx.and..not.qm%symm_packed) then
! ... symmetry unpacked (full) quaternion matrix
   write(LUPRI,'(2X,A)') 'this is full (symmetry unpacked) quaternion matrix '
   call FLSHFO(LUPRI)
   if (.not.associated(qm%p_rA3D1)) then
       call QUIT('QMtxInfo: qm%p_rA3D1 pointer for (symmetry unpacked)' &
                //' quaternion matrix is NOT associated !')
   endif
   ! ... print TOTAL QMatrix
   ntr = ( qm%nr(1) + qm%nr(2) )
   ntc = ( qm%nc(1) + qm%nc(2) )
   if (ntr .ne. ntc) call QUIT('QMtxInfo: no square matrix !')
   !CALL PRQMAT(PW(qm%kstart),ntr,ntc,ntr,ntc,qm%nz,IPQTOQ(1,qm%irep),LUPRI) 

   if (doprint) then
     CALL PRQMAT(qm%p_rA3D1(1,1,1),ntr,ntc,ntr,ntc,qm%nz,IPQTOQ(1,qm%irep),LUPRI) 
   endif
else if (qm%simple_mtx) then
! ... print out of a trivial quaternion matrix
   write(LUPRI,*) 'trivial quaternion matrix'
   call quit('QMtxInfo: trivial quaternion matrix needs to be implemented ')
else
   call quit('QMtxInfo: blind "if" branch ')
endif

  write(LUPRI,'(2x,a)')  &
  'no matrix print out - this matrix is not allocated'
  write(LUPRI,*) 'qm%simple_mtx=',qm%simple_mtx,' qm%symm_packed=',qm%symm_packed

endif

 call FLSHFO(LUPRI)

END SUBROUTINE QMtxInfo

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION MakeQMtxPck(qm_pk_in, irp_in, nz_in, text) result(qm)
!-------------------------------------------------------------------
!
!    basic constructor for the symmetry packed quaternion matrix
!
!-------------------------------------------------------------------
use Memory
#include "priunit.h"
#include "dgroup.h"
type(QMtx) :: qm
type(qmtx_pck), intent(in) :: qm_pk_in
type(irp), optional, intent(in) ::  irp_in
type(nzp), optional, intent(in) ::  nz_in
character(*), optional, intent(in) :: text
integer  :: isize
real*8, dimension(:), pointer :: temp

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

! ... array size of the symmetry packed matrix
isize = qm%nr(1)*qm%nc(1)*qm%nz + qm%nr(2)*qm%nc(2)*qm%nz
! save the size
qm%array_size = isize

! allocate array for the matrix
CALL mymemget('REAL',qm%kstart,isize)

temp => AccessWork(qm%kstart,isize) 
CALL DZERO(temp,isize)

qm%is_real  = .true.
qm%is_alloc = .true. ! array is allocated

! ... assign the pointer
qm%p_rA3D1 => NULL(); qm%p_rA3D2 => NULL()
call AssignQmtxPointer(qm)

!if (iprint_qmtx >=5 ) then
! call print_mem_param()
!endif

END FUNCTION MakeQMtxPck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignQmtxPointer(qm)
!----------------------------------------------------------------------
!
!           Assign pointers for a quaternion matrix
!
!----------------------------------------------------------------------
use Memory
use GeneralPointers
#include "priunit.h"
#include "dgroup.h"
type(QMtx),intent(inout) :: qm
!integer(kind=4) :: ifsym, ntr, ntc
integer :: ifsym, ntr, ntc
!integer (kind=8) :: ioff ! assumed bigg number 
integer :: ioff ! assumed bigg number 
real*8, dimension(:), pointer :: temp

if (.not.qm%is_alloc) then
  call QUIT('AssignQmtxPointer: matrix(array) NOT allocated !')
endif

if (qm%is_real) then

  if (.not.qm%symm_packed.and..not.qm%simple_mtx) then

    ntr = ( qm%nr(1) + qm%nr(2) ) ! square matrix assumed !
    ntc = ( qm%nc(1) + qm%nc(2) ) ! square matrix assumed !
   
    !if (ntr <> ntc) then
    if (ntr .ne. ntc) then
      call print_qmtx_attrib(qm)
      call QUIT('AssignQMtxArrayPointer: unpacked matrix is NOT square !')
    endif

   ! call AssignQMtx3D_real(PW(qm%kstart),iprint_qmtx,ntr,ntc,qm%nz)
    temp => AccessWork(qm%kstart)
   ! call AssignQMtx3D_real(AccessWork(qm%kstart),iprint_qmtx,ntr,ntc,qm%nz)
    call AssignQMtx3D_real(temp,iprint_qmtx,ntr,ntc,qm%nz)

    qm%p_rA3D1 => p3_real8;  p3_real8=> NULL()
 !   qm%p_rA3D1 => Access3Dpointer(qm%kstart)

    if (.not.associated(qm%p_rA3D1)) then
      call QUIT('AssignQMtxArrayPointer: p_rA3D1 not associated !')
    endif

  else if (qm%symm_packed.and..not.qm%simple_mtx) then

   ioff = qm%kstart ! offset
   do ifsym = 1, NFSYM
    ! call AssignQMtx3D_real(PW(ioff), iprint_qmtx, qm%nr(ifsym), qm%nc(ifsym), qm%nz )
  ! ...does not work !!! needed another method to reach pointer !
    temp => AccessWork(qm%kstart)
   ! call AssignQMtx3D_real(AccessWork(ioff), iprint_qmtx, qm%nr(ifsym), qm%nc(ifsym), qm%nz )
    call AssignQMtx3D_real( temp, iprint_qmtx, qm%nr(ifsym), qm%nc(ifsym), qm%nz )
     if (ifsym == 1) then
       qm%p_rA3D1 => p3_real8;  p3_real8=> NULL()
     else if (ifsym == 2) then
       qm%p_rA3D2 => p3_real8;  p3_real8=> NULL()
     else
       call QUIT('AssignQMtxArrayPointer: ifsym out of bounds (neither 1 nor 2) !')
     endif
     ioff = ioff + qm%ioff_pk
   enddo
  else
    call QUIT('AssignQMtxArrayPointer: blind branch !')
  endif 
else
  call QUIT('QMatrixes: AssignQmtxPointer - only real quaternion matrixes !')
endif

end subroutine AssignQmtxPointer

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
use Memory
#include "priunit.h"
#include "dgroup.h"
type(QMtx) :: qm
type(irp), optional, intent(in) ::  irp_in
type(nzp), optional, intent(in) ::  nz_in
type(qmtx_unpck),  intent(in) :: qm_upk_in
character(*), optional, intent(in) :: text
!integer (kind=8) :: isize
!integer (kind=8/IRAT) :: isize
integer :: isize
logical, save :: was_called = .false.

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

! ... total size of array to be allocated
isize = (qm%nr(1)+qm%nr(2))*(qm%nc(1)+qm%nc(2))*qm%nz
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

! ... assign the pointer
qm%p_rA3D1 => NULL(); qm%p_rA3D2 => NULL()
call AssignQmtxPointer(qm)

END FUNCTION MakeQMtxUnPck

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE DeAllocQMtx(qm)
!-----------------------------------------------------------------------------
!
!   Deallocate quaternion matrix; 
!
!-----------------------------------------------------------------------------
use Memory
#include "priunit.h"
type(QMtx) :: qm

!call MEMREL('QMatrixes: DestroyQMtx',PW,qm%kstart,qm%kstart,KFREE,LFREE)
call mymemrel('QMatrixes: DestroyQMtx',qm%kstart,qm%kstart)
qm%is_alloc = .false.
qm%array_size = 0
call cancel_string(qm%qmtx_text) ! cancel the string text
! shall I zero all attribites ?
END SUBROUTINE DeAllocQMtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function getdiff_qmtx(amat,bmat) result(diff)
!-------------------------------------------------------------------------
!
! returns the average difference of two FULLY COMPATIBLE quatern.matrixes
!
!-------------------------------------------------------------------------
use Memory
type(QMtx), intent(in) :: amat,bmat
real(kind=8) :: diff
logical(kind=1) :: compatible
integer(kind=4) :: i
! check matrix compatibilities
compatible = compatible_2qmtx(amat,bmat)

if (.not.compatible) then
  call print_qmtx_attrib(amat)
  call print_qmtx_attrib(bmat)
  call QUIT('getdiff_qmtx: matrixes are not compatible !')
endif

diff = 0.D0
do i=0,amat%array_size -1
! diff = abs(PW(amat%kstart+i) - PW(bmat%kstart+i)) + diff
 diff = abs(GetWorkElem(amat%kstart+i) - GetWorkElem(bmat%kstart+i)) + diff
enddo
! using size of the allocated array - depends on NZ !
diff = diff / dfloat(amat%array_size)

end function getdiff_qmtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function getdiff_part_qmtx(amat,bmat) result(diff)
!-------------------------------------------------------------------------
!
! returns the average difference of two quatern.matrixes
!
!-------------------------------------------------------------------------
use Memory
type(QMtx), intent(in) :: amat,bmat
real(kind=8) :: diff
logical(kind=1) :: compatible
integer(kind=8) :: size_min
integer(kind=4) :: i
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
 !diff = abs(PW(amat%kstart+i) - PW(bmat%kstart+i)) + diff
 diff = abs( GetWorkElem(amat%kstart+i) - GetWorkElem(bmat%kstart+i)) + diff
enddo
diff = diff / dfloat(amat%array_size)

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
use Memory  
type(QMtx), intent(inout) :: amat
type(QMtx), intent(in) :: bmat
real(kind=8),intent(in),optional :: fac_a,fac_b
real(kind=8) :: temp, temp1
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
 !temp = PW(amat%kstart+i)
 temp = GetWorkElem(amat%kstart+i)
! PW(amat%kstart+i)  =  fac_a*temp - fac_b*PW(bmat%kstart+i)
 temp1  =  fac_a*temp - fac_b*GetWorkElem(bmat%kstart+i)
 call PutWorkElem(amat%kstart+i,temp1) 
enddo

else if (present(fac_a).and..not.present(fac_b)) then

do i = 0, amat%array_size -1
! temp = PW(amat%kstart+i)
 temp = AccessWorkElem(amat%kstart+i)
! PW(amat%kstart+i)  =  fac_a*temp - PW(bmat%kstart+i)
 temp1  =  fac_a*temp - AccessWorkElem(bmat%kstart+i)
 call PutWorkElem(amat%kstart+i,temp1) 
enddo

else if (.not.present(fac_a).and..not.present(fac_b)) then

do i=0,amat%array_size -1
 !temp = PW(amat%kstart+i)
 temp = GetWorkElem(amat%kstart+i)
! PW(amat%kstart+i)  =  temp - PW(bmat%kstart+i)
 !GetWorkElem(amat%kstart+i)  =  temp - GetWorkElem(bmat%kstart+i)
 call PutWorkElem( amat%kstart+i, temp - GetWorkElem(bmat%kstart+i) )
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
use Memory  
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
 !temp = PW(amat%kstart+i)
 temp = GetWorkElem(amat%kstart+i)
 !PW(amat%kstart+i)  =  fac_a*temp - fac_b*PW(bmat%kstart+i)
 call PutWorkElem( amat%kstart+i,   fac_a*temp - fac_b*GetWorkElem(bmat%kstart+i) )
enddo

else if (present(fac_a).and..not.present(fac_b)) then

do i = 0, size_min -1
 !temp = PW(amat%kstart+i)
 temp = GetWorkElem(amat%kstart+i)
 !PW(amat%kstart+i)  =  fac_a*temp - PW(bmat%kstart+i)
 call PutWorkElem (amat%kstart+i,  fac_a*temp - GetWorkElem(bmat%kstart+i) )
enddo

else if (.not.present(fac_a).and..not.present(fac_b)) then

do i=0,size_min -1
 !temp = PW(amat%kstart+i)
 temp = GetWorkElem(amat%kstart+i)
 !PW(amat%kstart+i)  =  temp - PW(bmat%kstart+i)
 call PutWorkElem (amat%kstart+i, temp -  GetWorkElem(bmat%kstart+i) )
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

 write(LUPRI,'(/2x,a)')  & 
 '---------- Basic attributes of the QMtx matrix data type -----------'
 write(LUPRI,'(a,l1,$)') 'is_alloc:',amat%is_alloc
 write(LUPRI,'(a,l1,$)') '  simple_mtx:',amat%simple_mtx
 write(LUPRI,'(a,l1,$)') '  symmetry packed:',amat%symm_packed
 write(LUPRI,'(a,l1)') '  is real:',amat%is_real
 if (amat%is_alloc) then
  write(LUPRI,'(a,i7)') 'size of array:',amat%array_size
 endif
 write(LUPRI,'(a,2i3,a,2i3,a,i1,a,i2)') &
       'dimensions nr(1),nc(1); nr(2),nc(2):', &
       amat%nr(1),amat%nc(1),';',amat%nr(2),amat%nc(2),  & 
       ' --- NZ:',amat%nz,' IREP:',amat%irep
 call FLSHFO(LUPRI)
 if (is_string(amat%qmtx_text)) then
  write(LUPRI,'(2x,a,$)') '...attached text:'
  call print_string(amat%qmtx_text)
 endif

end subroutine print_qmtx_attrib
END MODULE QMatrixes

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine AssignQMtx3D_real(QMtx,iprint,n1,n2,nz)
!---------------------------------------------------------------------
!
! Assign general pointer to the 3-dim. array - quaternion matrix
!
!---------------------------------------------------------------------
use GeneralPointers
#include "priunit.h"
!integer(kind=4), intent(in) :: iprint,n1,n2,nz
integer, intent(in) :: iprint,n1,n2,nz
real(kind=8), target :: QMtx(n1,n2,nz)

!write(LUPRI,*) 'in AssignQMtx3Dreal1...'
!call FLSHFO(LUPRI)

p3_real8 => QMtx
if (iprint >= 10) then
  write(LUPRI,*) 'AssignQMtx3Dreal:1 real 3dim pointer to the matrix was assigned.'
  call FLSHFO(LUPRI)
endif

!write(LUPRI,*) 'in AssignQMtx3Dreal1 2...'
!call FLSHFO(LUPRI)

end subroutine AssignQMtx3D_real



