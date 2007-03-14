!
!
!    Tests of Object oriented Fortran90/95 features in DIRAC
!
!
  SUBROUTINE TEST_OOF(work,llwork)
!  ... this is very simple test of
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

   !.. correct order of deallocation = mirror order of allocatins -  is needed !!!!!!
   call DeAllocate_QMATRIX(rmat)
   call DeAllocate_QMATRIX(pmat)
   call DeAllocate_QMATRIX(qmat)

 call qexit('TEST_OOF')

! call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF

!
!  Handle the one-electron Fock matrix, VMAT - on entry
!
SUBROUTINE HANDLE_MATRIXES(VMAT1,TMAT1,WORK,LLWORK)
! VMAT - Fock SAAO matrix, TMAT - symm.packed transf. matrix
!
  use Memmory
  use QMatrixes
#include "priunit.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
 double precision, intent(inout) :: work(llwork) 
 double precision, intent(in) ::  vmat1(*),tmat1(*)
 type(QMtx) :: tmat

call qenter('HANDLE_MATRIXES')
! make new allocation for VMAT
    
   call UnStickArray ! release the previous work array
   call Stick_into_array(work,llwork) ! plug into entering WORK(LWORK) array

    call MakeQMtxPacked(tmat,NFBAS(1,0),NORB(1),NFBAS(2,0),NORB(2),NZ,0)
    call QMtxInfo(tmat)

call qexit('HANDLE_MATRIXES')
END SUBROUTINE HANDLE_MATRIXES

