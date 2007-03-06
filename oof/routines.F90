!
!
MODULE Memmory
  integer :: KWORK, KFREE, LFREE, LWORK
!  integer, save :: KWORK, KFREE, LFREE, LWORK
!  integer, save :: LWORK
!  integer :: KWORK, KFREE, LFREE
  double precision, pointer :: PW(:)
END MODULE Memmory

MODULE Memmory_handling
! IMPLICIT NONE

!PRIVATE
!PUBLIC:
CONTAINS

!
SUBROUTINE Stick_array(WORK,IWORK)

! use Memmory, ONLY: PW ! import pointer from PW
use Memmory

#include "priunit.h"

! integer, intent(in) :: LWORK
 integer, intent(inout) :: IWORK
 integer :: KSTART

! double precision, target, intent(in) :: WORK(LWORK)
 !double precision, target, intent(inout) :: WORK(LWORK)
double precision, intent(inout),  target :: WORK(IWORK)

 CALL QENTER('Stick_array')

   PW => WORK ! PW points to the big WORK array ... 

  write(LUPRI,*) 'Array WORK(',IWORK,') is sticked into PW!'
  write(LUPRI,*) 'length=',size(WORK)

   if (size(PW) /= IWORK) then
      write(LUPRI,*) 'somme error in array sticking !'
      CALL QUIT('Stick_array: error !')
   endif

   LWORK = IWORK

#include "memint.h"
      write(LUPRI,*) 'LFREE,LWORK', LFREE, LWORK

!   CALL MEMGET('REAL',KSTART,10,WORK,KFREE,LFREE) ! WORK(kstart)
    CALL MEMGET('REAL',KSTART,10,PW,KFREE,LFREE) ! WORK(kstart)

    write(LUPRI,*) 'this memget works, LWORK in Memmory=!',LWORK

! written PW prints the whole array  WORK !!

 CALL QEXIT('Stick_array')
  END SUBROUTINE Stick_array

  SUBROUTINE Allocate_memmory

! IMPLICIT NONE
 
!#include "implicit.h"
#include "priunit.h"

!  double precision :: a

  CALL QENTER('Allocate_memmory')

   write(lupri,*) 'memmory allocation on...' 
!#include "memint.h"

  CALL QEXIT('Allocate_memmory')

  END SUBROUTINE Allocate_memmory

END MODULE Memmory_handling

MODULE QMatrixes

type QMatrix
private
  integer::  nrows,ncols,nz ! 
  integer::  IREP
  integer::  kstart = -1
  logical::  isalloc=.false.
  double precision, pointer :: QM(:,:,:)
end type QMatrix 

contains

SUBROUTINE Allocate_QMATRIX(qmat,nrows,ncols,nz,IREP)

use Memmory  

type(QMatrix),intent(inout) :: qmat
integer, intent(in) :: nrows,ncols,nz,IREP
integer :: isize

#include "priunit.h"
!      ...  assign entering parameters  ...
  CALL QENTER('Allocate_QMATRIX')

  qmat%nrows = nrows
  qmat%ncols = ncols
  qmat%nz = nz
  qmat%IREP = IREP

!#include "memint.h"

 isize = qmat%nrows * qmat%ncols * qmat%nz  ! set up isize...

 write(LUPRI,*) 'size to allocate, isize=',isize
 write(LUPRI,*) ' PW (WORK)array, size=',size(PW)
 write(LUPRI,*) 'LWORK from Memory:',LWORK

#include "memint.h"
      write(LUPRI,*) 'after memint: KFREE,LFREE,KWORK', KFREE, LFREE, KWORK
CALL MEMGET('REAL',qmat%kstart,isize,PW,KFREE,LFREE) ! WORK(kstart)
  qmat%isalloc = .true.

  write(LUPRI,*) 'memmory allocated !.'

 ! qmat%QM = PW(qmat%kstart)
!  qmat%QM(1,1,1) = PW(qmat%kstart)

  write(LUPRI,*)  'pointer active ....'

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

  qmat%nrows = -1
  qmat%ncols = -1
  qmat%nz = -1
  qmat%IREP = -1
  qmat%isalloc = .false.

#include "memint.h"

isize = qmat%nrows * qmat%ncols * qmat%nz  ! set up isize...

 CALL MEMREL('PLOT',PW,qmat%kstart,qmat%kstart,KFREE,LFREE)

nullify (qmat%QM)

write(LUPRI,*)  'pointer active ....'

  CALL QEXIT('DeAllocate_QMATRIX')

END SUBROUTINE DeAllocate_QMATRIX


END MODULE QMatrixes

