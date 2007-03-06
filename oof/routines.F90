!
! Module contains basic variables of the DIRAC memmory array
!
MODULE Memmory

   integer :: KWORK, KFREE, LFREE, LWORK
   double precision,  pointer :: PW(:) ! points to the big WORK array

  contains

   subroutine print_mem_param() ! print out basic memmory variables
     #include "priunit.h"
     write(LUPRI,*) 'Memmory parameters KWORK, KFREE, LFREE, LWORK:', &
     KWORK, KFREE, LFREE, LWORK
   end subroutine print_mem_param

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
  ! double precision, allocatable :: QM(:,:,:) !...problematic
  ! double precision, pointer :: QM(:,:,:) !...problematic
  ! double precision, pointer :: QM(:) !...problematic
end type QMatrix 

contains

SUBROUTINE Allocate_QMATRIX(qmat,nrows,ncols,nz,IREP)

!use Memmory  

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
  write(LUPRI,*) 'memmory is allocated ...KSTART=',qmat%kstart

  CALL QEXIT('Allocate_QMATRIX')

END SUBROUTINE Allocate_QMATRIX

!
!
SUBROUTINE DeAllocate_QMATRIX(qmat)
type(QMatrix),intent(inout) :: qmat

#include "priunit.h"
!      ...  assign entering parameters  ...
  CALL QENTER('DeAllocate_QMATRIX')

!    CALL print_mem_param()
  CALL MEMREL('DeAllocate_QMATRIX',PW,qmat%kstart,qmat%kstart,KFREE,LFREE)

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

 if (qmat%isalloc == .false.) then
   call QUIT('MakeUnitQMATRIX: matrix is not allocated !')
 endif

 if (qmat%nrows /= qmat%ncols) then
   call QUIT('MakeUnitQMATRIX: matrix is not SQUARE !')
 endif

 CALL MakeUnitQMATRIX2(qmat,PW(qmat%kstart))

call QEXIT('MakeUnitQMATRIX')

END SUBROUTINE MakeUnitQMATRIX

SUBROUTINE  MakeUnitQMATRIX2(qmat,QM)
type(QMatrix),intent(in) :: qmat
double precision, intent(inout) :: QM(qmat%nrows,qmat%ncols,qmat%nz)

!  call basic routine from gp/gpsaue.F
   CALL DUNIT2(QM,qmat%nrows,qmat%nrows,qmat%ncols,qmat%nz) 

END SUBROUTINE  MakeUnitQMATRIX2

SUBROUTINE PrintQMATRIX(qmat)
type(QMatrix),intent(in) :: qmat
#include "priunit.h"

 CALL QENTER('PrintQMATRIX')

 if (qmat%isalloc == .false.) then
   call QUIT('PrintQMATRIX: matrix is not allocated !')
 endif
!
  Call PrintQMATRIX2(qmat,PW(qmat%kstart))

 CALL QEXIT('PrintQMATRIX')

END SUBROUTINE PrintQMATRIX

SUBROUTINE PrintQMATRIX2(qmat,QM) 
#include "priunit.h"
#include "dgroup.h"
type(QMatrix),intent(in) :: qmat
double precision, intent(inout) :: QM(qmat%nrows,qmat%ncols,qmat%nz)

  CALL QENTER('PrintQMATRIX2')

  CALL PRQMAT(QM,qmat%nrows,qmat%ncols,qmat%nrows,qmat%ncols,qmat%nz,  &
           IPQTOQ(1,qmat%IREP),LUPRI)

  CALL QEXIT('PrintQMATRIX2')

END SUBROUTINE PrintQMATRIX2

END MODULE QMatrixes

