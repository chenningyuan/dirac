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
 type(QMatrix) :: pmat,qmat,rmat  ! derived data type from the "QMATRIXES" module
call qenter('TEST_OOF')

  write(lupri,*)
  write(lupri,*) 'in test_oof routine'

   call Stick_into_array(work,llwork) ! take in WORK(LWORK) array ...

   call Allocate_QMATRIX(qmat,3,3,1,0) ! allocate quaternion matrix
   call Allocate_QMATRIX(pmat,3,3,1,0) ! allocate quaternion matrix

   call MakeUnitQMATRIX (qmat)
   call PrintQMATRIX (qmat)

   call MakeUnitQMATRIX (pmat)
   call PrintQMATRIX (pmat)

   call Allocate_QMATRIX(rmat,3,3,1,0) ! allocate quaternion matrix

! assign: rmat = pmat
   rmat = pmat
   write(LUPRI,*) 'Assigned rmat = pmat:'
   call PrintQMATRIX (rmat)

   !.. correct order of deallocation = mirror order of allocatiosn -  is needed !!!!!!
   call DeAllocate_QMATRIX(rmat)
   call DeAllocate_QMATRIX(pmat)
   call DeAllocate_QMATRIX(qmat)

 call qexit('TEST_OOF')

 call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF

