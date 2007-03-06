!
!
! Tests of Object oriented Fortran90/95 features in DIRAC
!
!
  SUBROUTINE TEST_OOF(work,lwork)

   use Memmory_handling
!  use Memmory
   use QMatrixes

#include "priunit.h"

 double precision :: work(lwork)
 type(QMatrix) :: qmat 
call qenter('TEST_OOF')

  write(lupri,*)
  write(lupri,*) 'in test_oof routine'

!      WORK is already in, continue
   call Stick_array(WORK,LWORK)
   call Allocate_QMATRIX(qmat,3,3,1,IREP)
   call DeAllocate_QMATRIX(qmat)

 call qexit('TEST_OOF')

! call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF


