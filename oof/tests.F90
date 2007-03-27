!
!    Tests of Object oriented Fortran90/95 features in DIRAC
!
  SUBROUTINE TEST_OOF(work,llwork)
!  ... this is very simple test of
!   use Memory
   use QMatrixes

#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"

#ifdef VARX
#include "dcbham.h"
#endif

#ifdef VARY
#include "dcbdhf.h"
#endif

 double precision, intent(inout) :: work(llwork) ! for testing...
 integer*4, intent(in) :: llwork
 type(QMtx) :: smat,pmat
call qenter('TEST_OOF')

CALL HEADER( 'in test_oof routine',-1)

call set_iprqmt(5)

call Stick_into_array(work,llwork)  ! take in WORK(LWORK) array ...

smat = MakeQMtx( set_qmtx_pck(NFBAS(1,0),NFBAS(1,0),NFBAS(2,0),NFBAS(2,0)) )
call QMtxInfo(smat)

!pmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(0),set_nz(1))
pmat = MakeQMtx(  set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0))  )
call QMtxInfo(pmat)

 call DeAllocQMtx(pmat)
 call DeAllocQMtx(smat)

 call test_oof1

  call UnStickArray

!   call Allocate_QMATRIX(qmat,3,3,1,0) ! allocate quaternion matrix
!   call Allocate_QMATRIX(pmat,3,3,1,0) ! allocate quaternion matrix
!   call MakeUnitQMATRIX (qmat)
!   call PrintQMATRIX (qmat)
!   call MakeUnitQMATRIX (pmat)
!   call PrintQMATRIX (pmat)
!   call Allocate_QMATRIX(rmat,3,3,1,0) ! allocate quaternion matrix
!   rmat = pmat
!   write(LUPRI,*) 'Assigned rmat = pmat:'
!   call PrintQMATRIX (rmat)
   !.. correct order of deallocation = mirror order of allocatins -  is needed !!!!!!
!   call DeAllocate_QMATRIX(rmat)
!   call DeAllocate_QMATRIX(pmat)
!   call DeAllocate_QMATRIX(qmat)

!  call UnStickArray

 call qexit('TEST_OOF')

! call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF

 subroutine test_oof1
!   use Memmory
   use QMatrixes
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
 type(QMtx) :: pmat

CALL HEADER( 'in test_oof1 routine',-1)

pmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(0),set_nz(1))
call QMtxInfo(pmat)

 end subroutine test_oof1

!
!  Handle the one-electron Fock matrix, VMAT - on entry
!
SUBROUTINE HANDLE_MATRIXES(VMAT1,TMAT1,WORK,LLWORK)
! VMAT - Fock SAAO matrix, TMAT - symm.packed transf. matrix
!
! use Memory
  use QMatrixes
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
 integer*4, intent(in) :: llwork
 double precision, intent(inout) :: work(llwork) 
 double precision, intent(in) ::  vmat1(*),tmat1(*)
! type(QMtx) :: tmat

call qenter('HANDLE_MATRIXES')
! make new allocation for VMAT
    
!   call UnStickArray ! release the previous work array
!   call Stick_into_array(work,llwork) ! plug into entering WORK(LWORK) array
!   call MakeQMtxPacked(tmat,NFBAS(1,0),NORB(1),NFBAS(2,0),NORB(2),NZ,0) 
!   call QMtxInfo(tmat)
!   call UnStickArray ! release the previous work array

call qexit('HANDLE_MATRIXES')
END SUBROUTINE HANDLE_MATRIXES

