!
!
! Tests of Object oriented Fortran90/95 features in DIRAC
!
!
  SUBROUTINE TEST_OOF(work,lwork)

   use Memmory
   use QMatrixes

#include "priunit.h"

 double precision, intent(inout) :: work(lwork) ! for testing...
 type(QMatrix) :: qmat 
call qenter('TEST_OOF')

  write(lupri,*)
  write(lupri,*) 'in test_oof routine'

   call Stick_array(WORK,LWORK) ! take in WORK(LWORK) array ...
!      WORK is already in, continue
   call Allocate_QMATRIX(qmat,3,3,1,0) ! allocate array...

   call MakeUnitQMATRIX (qmat)
   call PrintQMATRIX (qmat)

   call DeAllocate_QMATRIX(qmat)

 call qexit('TEST_OOF')
! call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF

