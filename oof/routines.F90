!
! Module contains basic variables of the DIRAC memmory array
!
MODULE Memmory

integer :: KWORK, KFREE, LFREE, LWORK ! basic variables of memmory handling
double precision,  pointer :: PW(:) ! points to the big WORK array

  contains

   ! print out basic memmory variables
   subroutine print_mem_param 
#include "priunit.h"
   write(LUPRI,*) 'Memmory parameters: KWORK=',KWORK,' KFREE=',KFREE,' LFREE=',LFREE, &
   ' LWORK=',LWORK 

   if (associated(PW)) then
     write(LUPRI,*) 'Pointer PW is associated with the WORK array.'
   else
     write(LUPRI,*) 'Pointer PW is NOT associated with the WORK array.'
  endif
   end subroutine print_mem_param

SUBROUTINE Stick_into_array(WORK,LLWORK)
#include "priunit.h"
 integer, intent(in) :: LLWORK
 double precision, intent(in),  target :: WORK(LLWORK)
 CALL QENTER('Stick_into_array')
 if (associated(PW)) then
   write(LUPRI,*) 'Stick_into_array: PW is already associated - Exit !'
   CALL QUIT('Stick_into_array: PW is already associated !')
 endif
 PW => WORK ! PW points now to the big WORK array ... 
 if (size(PW) /= LLWORK) then
    write(LUPRI,*) 'somme error in array sticking !'
    CALL QUIT('Stick_into_array: error in length !')
 endif
 LWORK = LLWORK ! story into module variable
 CALL QEXIT('Stick_into_array')
END SUBROUTINE Stick_into_array

END MODULE Memmory

!
!
! This module contains tools for working with quaternion matrixes
!
MODULE QMatrixes

use Memmory  

type QMatrix
private
  integer::  nrows=-1,ncols=-1,nz=-1 ! 
  integer::  IREP=-1
  integer::  kstart = -1 ! contains position on PW (WORK) array !
  logical::  isalloc=.false.
end type QMatrix 

interface assignment ( = )
  module procedure assign_qmatrix
end interface 

contains

!
!    Constructor for a simple quaternion matrix.
!
SUBROUTINE Allocate_QMATRIX(qmat,nrows,ncols,nz,IREP)
#include "priunit.h"
type(QMatrix),intent(inout) :: qmat
integer, intent(in) :: nrows,ncols,nz,IREP
integer :: isize
logical, save :: was_called = .false.

  CALL QENTER('Allocate_QMATRIX')

! Check is PW pointer is associated with the DIRAC WORK array
  if (.NOT.associated(PW)) call QUIT('PW pointer MUST be associated with the &
 DIRAC WORK array!')

  if (qmat%kstart.ne.-1) then
   call QUIT('The matrix is already allocated !')
  endif

! Check dimensions of the quaternion matrix
  if (nz.ne.1.and.nz.ne.2.and.ne.ne.4) then
   write(LUPRI,*) 'Allocate_QMATRIX: entering NZ=',nz
   CALL QUIT('Allocate_QMATRIX: NZ  must be 1 or 2 or 4 !')
  endif
  if (.not.(nrows>0.and.ncols>0)) then
   write(LUPRI,*) 'Allocate_QMATRIX: entering nrows,nclos=',nrows,ncols
   CALL QUIT('Allocate_QMATRIX: entering nrows,ncols must be greater than zero !')
  endif

!   ...  assign entering parameters  ...
  qmat%nrows = nrows
  qmat%ncols = ncols
  qmat%nz = nz
  qmat%IREP = IREP

 isize = qmat%nrows * qmat%ncols * qmat%nz  ! set up isize...

 write(LUPRI,*) 'size to allocate, isize=',isize
 write(LUPRI,*) ' PW (WORK)array, size=',size(PW)
 write(LUPRI,*) 'LWORK from Memory:',LWORK

if (.not.was_called) then
#include "memint.h"
was_called=.true.
write(LUPRI,*) '1.calling!'
else
write(LUPRI,*) 'not 1.calling!'
endif

    CALL print_mem_param()

  CALL MEMGET('REAL',qmat%kstart,isize,PW,KFREE,LFREE) ! WORK(kstart)
  qmat%isalloc = .true.
  write(LUPRI,*) 'memmory is allocated ...KSTART=',qmat%kstart
     CALL print_mem_param()
 
! by default, fill the matrix with zeroes 
  CALL DZERO(PW(qmat%kstart),isize) 

  CALL QEXIT('Allocate_QMATRIX')

END SUBROUTINE Allocate_QMATRIX

!
! Destructor deallocates the memmory.
!
SUBROUTINE DeAllocate_QMATRIX(qmat)
type(QMatrix),intent(inout) :: qmat

#include "priunit.h"

  CALL QENTER('DeAllocate_QMATRIX')

     write(LUPRI,*) 'before deallocation: qmat%kstart=',qmat%kstart
    CALL print_mem_param()
  CALL MEMREL('DeAllocate_QMATRIX',PW,qmat%kstart,qmat%kstart,KFREE,LFREE)
     write(LUPRI,*) 'after deqllocation:'
    CALL print_mem_param()

  qmat%nrows = -1
  qmat%ncols = -1
  qmat%nz = -1
  qmat%IREP = -1
  qmat%isalloc = .false.

  CALL QEXIT('DeAllocate_QMATRIX')

END SUBROUTINE DeAllocate_QMATRIX

SUBROUTINE MakeUnitQMATRIX(qmat)
type(QMatrix),intent(inout) :: qmat
#include "priunit.h"

call QENTER('MakeUnitQMATRIX')

 if (.not.qmat%isalloc) then
   call QUIT('MakeUnitQMATRIX: matrix is not allocated !')
 endif

 if (qmat%nrows /= qmat%ncols) then
   call QUIT('MakeUnitQMATRIX: matrix is not SQUARE for making it unit !')
 endif

 ! call it to make (:,:,:) dimension structure of the PW matrix array !
 CALL MakeUnitQMATRIX2(qmat,PW(qmat%kstart))

call QEXIT('MakeUnitQMATRIX')

END SUBROUTINE MakeUnitQMATRIX

SUBROUTINE  MakeUnitQMATRIX2(qmat,QM)
type(QMatrix),intent(in) :: qmat
double precision, intent(inout) :: QM(qmat%nrows,qmat%ncols,qmat%nz) 
 CALL QENTER('MakeUnitQMATRIX2')

!  call basic routine from gp/gpsaue.F
   CALL DUNIT2(QM,qmat%nrows,qmat%nrows,qmat%ncols,qmat%nz) 

 CALL QEXIT('MakeUnitQMATRIX2')
END SUBROUTINE  MakeUnitQMATRIX2

SUBROUTINE PrintQMATRIX(qmat)
type(QMatrix),intent(in) :: qmat
#include "priunit.h"
 CALL QENTER('PrintQMATRIX')

! if (qmat%isalloc == .false.) then
 if (.NOT.qmat%isalloc) then
   call QUIT('PrintQMATRIX: matrix is not allocated !')
 endif

 ! call it to make (:,:,:) dimension structure of the PW matrix array !
  Call PrintQMATRIX2(qmat,PW(qmat%kstart))

CALL QEXIT('PrintQMATRIX')

END SUBROUTINE PrintQMATRIX

SUBROUTINE PrintQMATRIX2(qmat,QM) 
#include "priunit.h"
#include "dgroup.h"
type(QMatrix),intent(in) :: qmat
double precision, intent(inout) :: QM(qmat%nrows,qmat%ncols,qmat%nz)
CALL QENTER('PrintQMATRIX2')
! this is the lower layer routine from qp/
CALL PRQMAT(QM,qmat%nrows,qmat%ncols,qmat%nrows,qmat%ncols,qmat%nz,  &
            IPQTOQ(1,qmat%IREP),LUPRI)
CALL QEXIT('PrintQMATRIX2')
END SUBROUTINE PrintQMATRIX2

Subroutine assign_qmatrix(pmat,qmat) ! pmat: = qmat 
#include "priunit.h"
type(QMatrix),intent(inout) :: pmat
type(QMatrix),intent(in) :: qmat
! First check 
!   write(LUPRI,*) 'rows of pmat,qmat:',pmat%nrows,qmat%nrows    
 if (.not.check_qmatrixes(pmat,qmat)) then
   write(LUPRI,*) 'rows of pmat,qmat:',pmat%nrows,qmat%nrows    
   write(LUPRI,*) 'columns of pmat,qmat:',pmat%ncols,qmat%ncols
   write(LUPRI,*) 'NZ of pmat,qmat:',pmat%nz,qmat%nz
   call QUIT('entering quaternion matrixes do not have same dimensions !')
 endif
 isize = pmat%nrows*pmat%ncols*pmat%nz
 write(lupri,*) 'isize for copy:',isize
 write(lupri,*) 'qmat%kstart,pmat%kstart:',qmat%kstart,pmat%kstart
! copy matrixes-parts of PW(WORK) arrays - using lower layer library routine
 CALL DCOPY(isize,PW(qmat%kstart),1,PW(pmat%kstart),1)
END subroutine assign_qmatrix

Function check_qmatrixes(pmat,qmat) result (good)
type(QMatrix), intent(in) :: pmat,qmat
logical :: good
CALL QENTER('check_qmatrixes')
good = ( (pmat%nrows == qmat%nrows) .AND. &
         (pmat%ncols == qmat%ncols) .AND. &
         (pmat%kstart /= qmat%kstart) .AND. &
         (pmat%nz == qmat%nz) .AND. pmat%isalloc .AND. qmat%isalloc ) 
CALL QEXIT('check_qmatrixes')
END Function check_qmatrixes

END MODULE QMatrixes

