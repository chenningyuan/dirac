!
!
!    Tests of Object oriented Fortran90/95 features in DIRAC
!
!
  SUBROUTINE TEST_OOF(work,llwork)

   use Memmory
   use QMatrixes

#include "priunit.h"

 double precision, intent(inout) :: work(llwork) ! for testing...
 type(QMatrix) :: qmat  ! derived data type from the "QMATRIXES" module
call qenter('TEST_OOF')

  write(lupri,*)
  write(lupri,*) 'in test_oof routine'

   call Stick_into_array(work,llwork) ! take in WORK(LWORK) array ...

   call Allocate_QMATRIX(qmat,3,3,1,0) ! allocate quaternion matrix

   call MakeUnitQMATRIX (qmat)
   call PrintQMATRIX (qmat)

   call DeAllocate_QMATRIX(qmat)

 call qexit('TEST_OOF')
! call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF

