!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

!------------------------------------------------------------------
!
!    Tests of Object oriented Fortran90/95 features in DIRAC
!
!------------------------------------------------------------------


SUBROUTINE TEST_OO_MEMORY(WORK,LWORK)

 use Memory

#include "implicit.h"
#include "priunit.h"

 double precision, intent(inout) :: work(lwork) !

call set_iprmem(20) 
call StickMemory(WORK,LWORK)  ! take in WORK(LWORK) array ...

call print_mem_param

call UnStickMemory


!call quit('in TEST_OO_MEMORY')

END SUBROUTINE TEST_OO_MEMORY



  SUBROUTINE TEST_OOF(work,llwork)
 ! SUBROUTINE TEST_OOF

   use Memory
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
! integer*4, intent(in) :: llwork
 type(QMtx) :: smat,pmat
 logical :: domtx

 !!!  CALL TEST_MEMGET_CHAR(work,llwork)

!domtx = .false.
domtx = .true.

if (domtx) then

CALL HEADER( 'in test_oof routine',-1)
call StickMemory(work,llwork)  ! take in WORK(LWORK) array ...
call set_iprqmtx(100)
call set_iprmem(100)

!call StickMemory(work,llwork)  ! take in WORK(LWORK) array ...

! call print_mem_param

write(LUPRI,*) 'test_oof: smat called ... '
!smat = MakeQMtx( set_qmtx_pck(NFBAS(1,0),NFBAS(1,0),NFBAS(2,0),NFBAS(2,0)) )
smat = MakeQMtx( set_qmtx_pck(NFBAS(1,0),NFBAS(1,0)), text='testing smat' )
!call QMtxInfo(smat,  print_qmtx = .false.)
call QMtxInfo(smat)
write(LUPRI,*) 'test_oof: after smat called ... '

!pmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(0),set_nz(1))
!pmat = MakeQMtx(  set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0))  )
!pmat = MakeQMtx(  set_qmtx_unpck(NFBAS(1,0))  )
!call QMtxInfo(pmat, print_qmtx = .false.)

write(LUPRI,*) 'test_oof: after pmat called ... '

! call DeAllocQMtx(pmat)
 call DeAllocQMtx(smat)

 call test_oof1

 call UnStickMemory

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

 endif

! call quit('TEST_OOF')   

  END SUBROUTINE TEST_OOF

  subroutine test_oof1
use QMatrixes
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
 type(QMtx) :: pmat

CALL HEADER( 'in test_oof1 routine !',-1)

!   make a new matrix ...
pmat = MakeQMtx( set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)), set_irp(0), set_nz(1) )
call QMtxInfo(pmat)

 end subroutine test_oof1

!
!  Handle the one-electron Fock matrix, VMAT - on entry
!
SUBROUTINE HANDLE_MATRIXES(VMAT1,TMAT1,WORK,LLWORK)
!
! VMAT - Fock SAAO matrix, TMAT - symm.packed transf. matrix
!
! use Memory
  use QMatrixes
#include "implicit.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
! integer*4, intent(in) :: llwork
 integer, intent(in) :: llwork
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


subroutine test_oof_basis_fu(WORK,LWORK)
!------------------------------------------------------------
!
!
!
! control print out of atomic basis function; exponents and contraction coefficients from non-abacus
! 
! part of the DIRAC
!
!
!
!------------------------------------------------------------
use memory
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "dcbgen.h"
#include "cbirea.h"
#include "molinp.h"
#include "ccom.h"
#include "nuclei.h"
#include "frame.h"
#include "primit.h"
#include "shells.h"
#include "symmet.h"
#include "chrsgn.h"
#include "aosotr.h"
#include "infpar.h"


real*8 WORK(LWORK)
logical SPHRA


call StickMemory(WORK,LWORK)
 
  call HEADER('test_oof_basis_fu',-1)

      write(lupri,*) 'NLARGE=',NLARGE,' NSMALL=',NSMALL

      write(lupri,*) 'NUCIND=',NUCIND
      write(lupri,*) 'NLRGSH=',NLRGSH
      write(lupri,*) 'NSMLSH=',NSMLSH
      write(lupri,*) 'KMAX=',KMAX
      write(lupri,*) 'NORBS=',NORBS

      write(lupri,*) 'charge pre nucleus:',(charge(i),i=1,nucind)

      write(lupri,*) 'nuclear degeneracy:',(nucdeg(i),i=1,nucind)

      write(lupri,*) 'Total number of centers:',NUCDEP
      write(lupri,*)


   ! run over shells...
     ! DO ISHELL = 1, NLRGSH
      DO ISHELA = 1, KMAX

     !  L     = NHKT(ISHELL)
     !  ICENT = NCENT(ISHELL)
     !  NDEG  = NUCDEG(ICENT)
     !  INDX=NUMCF(ISHELL)

         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)

         SPHRA = SPHR(ISHELA) 
         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)

         NCENTA = NCENT(ISHELA)
         ICENTA = NUCNUM(NCENTA,1)
         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)

         NUCDEGA  = NUCDEG(NCENTA)       
       !     ... offset for normalized coefficients ...
         JSTA   = JSTRT(ISHELA)
       !    ... coordinates
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)

   !    write(LUPRI,*) 'ISHELL=',ISHELL,' L=',L
   !    write(LUPRI,*) 'ICENT=',ICENT,' NDEG=',NDEG
   !    write(LUPRI,*) 'KHKT(ISHELL)=',KHKT(ISHELL)
   !    write(LUPRI,*) 'KCKT(ISHELL)=',KCKT(ISHELL)
   !    write(LUPRI,*) 'ICA(ISHELL)=',LCLASS(ISHELL)
   !    write(LUPRI,*) 'SPHR(ISHELL)=',SPHR(ISHELL)
   !    write(LUPRI,*) 'NCENT(ISHELL)=',NCENT(ISHELL)
   !    write(LUPRI,*) 'NUCNUM(NCENT(ISHELL),1)=',NUCNUM(NCENT(ISHELL),1)
   !    write(LUPRI,*) 'MULT(ISTBAO(ISHELL))=',MULT(ISTBAO(ISHELL))
   !    write(LUPRI,*) 'NUCO(ISHELL)=',NUCO(ISHELL)
   !    write(LUPRI,*) 'NUMCF(ISHELL)=',NUMCF(ISHELL)
   !    write(LUPRI,*) 'JSTRT(ISHELL)=',JSTRT(ISHELL)

        write(LUPRI,*)
        write(LUPRI,'(2X,A,I3,A,I3,A,I2)') 'ISHELA/KMAX=',ISHELA,'/',KMAX,' large(1)/small(2) comp., ICA:',ICA
        write(LUPRI,'(2X,A,I3,A,I3)') 'NHKTA(alias L)=',NHKTA,' cartesian degeneracy (KHKTA)=',KHKTA

    !    write(LUPRI,'(2X,A,I3,A,I3,A,I3)') 'ISHELA=',ISHELA,' NUCA=',NUCA,' NUCMFA=',NUMCFA
        write(LUPRI,'(2X,A,I3,A,I3)') 'contracted shell (NUMCFA) #/ no f of uncotracted shells (NUCA)=',NUMCFA,'/',NUCA
        write(LUPRI,'(2X,A,I3,A,I3)') 'symmetric center (NCENTA) / total nu of symm.nucleis (NUCIND)=',NCENTA,'/',NUCIND

        write(LUPRI,'(2X,A,L1)')      'logical SPHRA=',SPHRA 

        write(LUPRI,'(2X,A,I3,A,I3,A,F10.3,A,A4,A,I3)')  &
   'ICENTA=',ICENTA,' NCENTA=',NCENTA,' charge:',CHARGE(NCENTA),' name:',NAMN(NCENTA), 'deg:',NUCDEGA

        write(LUPRI,'(2X,A,3I3)')    'L,M,N=',LVALUA,MVALUA,NVALUA
        write(LUPRI,'(2X,A,3F10.4)')  'coordinates:',CORAX,CORAY,CORAZ
 
        write(LUPRI,'(2X,A)') 'Exponents and normalized contract.coefficients...'
        DO  IPRIMA = 1,NUCA
           JPRIMA = JSTA + IPRIMA
           CONTA = PRICCF(JPRIMA,NUMCFA)
           EXPA = PRIEXP(JPRIMA)
           write(LUPRI,'(2X,A,F16.6,A,F20.12)') ' exp:',EXPA,'  normalized contr.coeff:',CONTA
         ENDDO

    !   write(LUPRI,*)

      ENDDO


 ! call quit('test_oof_basis_fu')

  call UnstickMemory

end subroutine test_oof_basis_fu


  SUBROUTINE TEST_MEMGET_CHAR(work,lwork)
! -------------------------------------------------------------
!
!
!
! -------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
  DIMENSION WORK(LWORK)

  CALL HEADER( 'in TEST_MEMGET_CHAR routine',-1)

#include "memint.h"
  IS1_1 = 3; IS1_2=10
  CALL MEMGET('CHAR',KCHAR1,IS1_1*IS1_2,WORK,KFREE,LFREE)
  IS2_1 = 4; IS2_2 = 5; IS2_3 = IS1_2
  CALL MEMGET('CHAR',KCHAR2,IS2_1*IS2_2*IS2_3,WORK,KFREE,LFREE)

  write(LUPRI,*) 'allocated !!!!'

  CALL TEST_CHAR_ARR(IS1_1,IS1_2,WORK(KCHAR1),IS2_1,IS2_2,IS2_3,WORK(KCHAR2))

  CALL MEMREL('TEST_MEMGET_CHAR',WORK,KWORK,KCHAR1,KFREE,LFREE)

! CALL QUIT('TEST_MEMGET_CHAR')
  END SUBROUTINE TEST_MEMGET_CHAR

SUBROUTINE TEST_CHAR_ARR(N1_1,N1_2,ARR1,N2_1,N2_2,N2_3,ARR2)
! ---------------------------------------------------------------------
!
!
! ---------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
CHARACTER*3 ARR1(N1_2)
!CHARACTER ARR1(N1_1,N1_2)
!CHARACTER*4 ARR2(5,10)
CHARACTER*4 ARR2(N2_2,N2_3)


DO I=1,N1_2
  ARR1(I) = 'abcdefghi'
  write(LUPRI,*) I,ARR1(I)
  DO J = 1, N2_2
    ARR2(I,J)='xyz123'
    write(LUPRI,*) I,J,ARR2(I,J)
  ENDDO
ENDDO

CALL QEXIT('TEST_CHAR_ARR')
RETURN
END

!===============================================================
!
!===============================================================
SUBROUTINE  TEST_CHAR2(ARR2,N2_2,N2_3)
#include "implicit.h"
#include "priunit.h"
CHARACTER*4 ARR2(N2_2,N2_3),Z

 DO I=1, N2_3
 DO J=1, N2_2
   Z = ARR2(I,J); write(LUPRI,*) 'Z=',Z
   ARR2(I,J)='nove'
   write(LUPRI,*) I,J,Z,' new:',ARR2(I,J)
 ENDDO
 ENDDO

RETURN
END

