!
! Module contains basi variables of the memmory initialization
!
MODULE Memmory
   integer :: KWORK, KFREE, LFREE, LWORK
   double precision,  pointer :: PW(:) ! points to the big WORK array
  ! double precision,  pointer :: PW ! points to the big WORK array
  contains
   subroutine print_mem_param() ! print out basic memmory variables
     #include "priunit.h"
     write(LUPRI,*) 'Memmory parameters KWORK, KFREE, LFREE, LWORK:', &
     KWORK, KFREE, LFREE, LWORK
   end subroutine print_mem_param
END MODULE Memmory

MODULE Memmory_handling

use Memmory

CONTAINS

SUBROUTINE Stick_array(WORK,IWORK)

 #include "priunit.h"

 integer, intent(in) :: IWORK

 double precision, intent(in),  target :: WORK(IWORK)

 CALL QENTER('Stick_array')

  PW => WORK ! PW points to the big WORK array ... 

!  write(LUPRI,*) 'Array WORK(',IWORK,') is sticked into PW!'
!  write(LUPRI,*) 'length=',size(WORK)

   if (size(PW) /= IWORK) then
      write(LUPRI,*) 'somme error in array sticking !'
      CALL QUIT('Stick_array: error in length !')
   endif

  LWORK = IWORK ! part of memmory
  CALL QEXIT('Stick_array')

END SUBROUTINE Stick_array

END MODULE Memmory_handling

MODULE QMatrixes

type QMatrix
private
  integer::  nrows=-1,ncols=-1,nz=-1 ! 
  integer::  IREP=-1
  integer::  kstart = -1
  logical::  isalloc=.false.
  double precision, allocatable :: QM(:,:,:) !...problematic
  ! double precision, pointer :: QM(:,:,:) !...problematic
  ! double precision, pointer :: QM(:) !...problematic
end type QMatrix 

contains

SUBROUTINE Allocate_QMATRIX(qmat,nrows,ncols,nz,IREP)

use Memmory  

#include "priunit.h"
type(QMatrix),intent(inout) :: qmat
integer, intent(in) :: nrows,ncols,nz,IREP
integer :: isize

!      ...  assign entering parameters  ...
  CALL QENTER('Allocate_QMATRIX')

!      ...  assign entering parameters  ...
  qmat%nrows = nrows
  qmat%ncols = ncols
  qmat%nz = nz
  qmat%IREP = IREP

 isize = qmat%nrows * qmat%ncols * qmat%nz  ! set up isize...

 write(LUPRI,*) 'size to allocate, isize=',isize
 write(LUPRI,*) ' PW (WORK)array, size=',size(PW)
 write(LUPRI,*) 'LWORK from Memory:',LWORK

#include "memint.h"

!     CALL print_mem_param()

  CALL MEMGET('REAL',qmat%kstart,isize,PW,KFREE,LFREE) ! WORK(kstart)
  qmat%isalloc = .true.
 ! write(LUPRI,*) 'memmory is allocated !.'

!    CALL print_mem_param()

!.... does not work !!!! ...
!  qmat%QM = PW(qmat%kstart)
  ! qmat%QM(1) = PW(qmat%kstart)
!   qmat%QM = PW(qmat%kstart)
 !  qmat%QM(1,1,1) = PW(qmat%kstart)

!  write(LUPRI,*)  'pointer active ....'

  CALL QEXIT('Allocate_QMATRIX')

END SUBROUTINE Allocate_QMATRIX

!
!
!
!
SUBROUTINE DeAllocate_QMATRIX(qmat)

 use Memmory  

type(QMatrix),intent(inout) :: qmat

#include "priunit.h"
!      ...  assign entering parameters  ...
  CALL QENTER('DeAllocate_QMATRIX')

!    CALL print_mem_param()
!#include "memint.h"

! isize = qmat%nrows * qmat%ncols * qmat%nz  ! set up isize...

!    CALL print_mem_param()
 CALL MEMREL('DeAllocate_QMATRIX',PW,qmat%kstart,qmat%kstart,KFREE,LFREE)

  qmat%nrows = -1
  qmat%ncols = -1
  qmat%nz = -1
  qmat%IREP = -1
  qmat%isalloc = .false.

!nullify (qmat%QM)

  CALL QEXIT('DeAllocate_QMATRIX')

END SUBROUTINE DeAllocate_QMATRIX

END MODULE QMatrixes

