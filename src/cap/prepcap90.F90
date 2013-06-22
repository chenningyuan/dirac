!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine CAP_DIPVEL_COMP
!----------------------------------------------------------------------------
!
! compare CAP and DIRAC dipole velocity integrals, which are now 
! adapted and written  on the AOPROPER files
!
! called from CAP_STUFF
!
!----------------------------------------------------------------------------
 ! use Operators
 ! use Operators
   use QMatrixes
#include "implicit.h"
#include "priunit.h"
#include "dcbham.h"
#include "dcbcap.h"
 ! type(DirOper) :: cap_dipvel, dir_dipvel
  type(QMtx) :: cap_dipvel, dir_dipvel 
! adjust the threshold
 double precision, parameter :: diff_thrs = 1.0D-5
 call QENTER('CAP_DIPVEL_COMP')

! call add_operator(cap_dipvel,IPCAP_DV); call add_operator(dir_dipvel,IPKINERG)

 call get_prop_saao(IPCAP_DV,cap_dipvel)
 call AddTextQMtx(cap_dipvel,'CAP dipole velocity integrals')

 call get_prop_saao(IPKINERG,dir_dipvel)
 call AddTextQMtx(dir_dipvel,'DIRAC dipole velocity integrals')

!  compare smat & capovlp quaternion matrixes
 ! diff = getdiff_part_qmtx(cap_dipvel,dir_dipvel)
  diff = getdiff_qmtx(cap_dipvel,dir_dipvel)
  write(LUPRI,'(/A)') '*************************************************************'
  write(LUPRI,'(A,D22.10)') '*** comparison of CAP-DIRAC dipole velocity operators: ',diff
  write(LUPRI,'(A)')  '*************************************************************'
  CALL FLSHFO(LUPRI)

 if (diff > diff_thrs) then
! ... print individual matrixes, plus difference between them
  call QMTxInfo(cap_dipvel)
  call QMTxInfo(dir_dipvel)

  call  AddTextQMtx(cap_dipvel,'******** difference  CAP - DIRAC ***********')
  call substract_part_qmtx(cap_dipvel,dir_dipvel)
  call QMTxInfo(cap_dipvel)
  write(LUPRI,'(A,D22.10)') 'absolute difference between CAP-DIRAC dipole velocity operators: ',diff
  CALL FLSHFO(LUPRI)
  call quit('CAP_DIPVEL_COMP: big differences between'//  &
            ' CAP/DIRAC dipole velocity integrals !')
 endif

 call QEXIT('CAP_DIPVEL_COMP')

 end subroutine CAP_DIPVEL_COMP

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine cap_ovlp_comp
!---------------------------------------------------------------------------
!
!    compare overlap integrals from DIRAC and from Sajeev's CAP integral
! generator
!
!  called from CAP_STUFF
!
!----------------------------------------------------------------------------
    use QMatrixes
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
type(QMtx) :: smat,capovl
double precision, parameter :: diff_thrs = 1.0D-7
call QENTER('cap_ovlp_comp')

 call set_iprqmtx(IPRCAP);
! write(LUPRI,*) 'iprqmtx=',get_iprint_qmtx()

 call GetUnpckOvlpMtx(smat)
! call QMTxInfo(smat)

 call GetCapOvlp(capovl)
 !call QMTxInfo(capovl)

!  compare smat & capovlp quaternion matrixes
  diff = getdiff_part_qmtx(smat,capovl)
  write(LUPRI,'(/A)') '*************************************************************'
  write(LUPRI,'(A,D22.10)') '*** comparison of CAP-DIRAC overlap matrixes: ',diff
  write(LUPRI,'(A)')  '*************************************************************'
  CALL FLSHFO(LUPRI)

 if (diff > diff_thrs) then
! ... print individual matrixes, plus difference between them
  call QMTxInfo(smat)
  call QMTxInfo(capovl)

  call  AddTextQMtx(smat,'******** difference smat - capovlp ***********')
  call substract_part_qmtx(smat,capovl)
  call QMTxInfo(smat)

  call quit('CAP_DIPVEL_COMP: big differences between'//  &
            ' CAP/DIRAC overlap integrals !')
 endif

call QEXIT('cap_ovlp_comp')
end subroutine cap_ovlp_comp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GetCapOvlp(capovlp_mtx)
!--------------------------------------------------------------------
!
! get CAP overlap matrix into "capovlp_mtx" matrix data type
!
!--------------------------------------------------------------------
   use QMatrixes
#include "priunit.h"
!#include "implicit.h"
!#include "dcbbas.h"
!#include "dcborb.h"
#include "dcbcap.h"
!#include "dcbham.h"
 type(QMtx), intent(inout) :: capovlp_mtx
call QENTER('GetCapOvlp')

 call get_prop_saao(IPCAP_S,capovlp_mtx) 
 call AddTextQMtx(capovlp_mtx,'CAP overlap integrals')

 if (get_iprint_qmtx() >= 7) then 
! if (get_iprint_qmtx() >= 0) then 
  write(LUPRI,*) 'GetCapOvlp: CAP overlap matrix...'
  call print_qmtx_attrib(capovlp_mtx)
  call QMtxInfo(capovlp_mtx)
 endif

call QEXIT('GetCapOvlp')
end subroutine GetCapOvlp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GetUnpckOvlpMtx(smat)
!===========================================================================
!
!
!     get unpacked overlap matrix, "smat" - allocate inside, then return
!
!
!===========================================================================
   use Memory
   use QMatrixes
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcbham.h"
type(QMtx), intent(inout) :: smat
call QENTER('GetUnpckOvlpMtx')

 smat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(0),set_nz(1), & 
                 'DIRAC unpacked overlap SA-AO matrix' )
! call QMtxInfo(smat)
! call print_mem_param()

 if (get_iprint_qmtx() >= 7) then
   write(LUPRI,*) 'N2BAST=',N2BAST
   write(LUPRI,*) 'get_kfree, get_lfree=',get_kfree(),get_lfree()
   write(LUPRI,*) 'get_qmtx_kstart(smat)=',get_qmtx_kstart(smat)
   !write(LUPRI,*) 'SSMTRC=',SSMTRC, ' but I set 1.0'
   write(LUPRI,*) 'SSMTRC=',SSMTRC
 endif

 ILFREE =  get_lfree() ! must be due to integre conversion !!!

! call lower level routine
 !call GTOVLX( PW(get_qmtx_kstart(smat)), SSMTRC, PW(get_kfree()), ILFREE )
 call GTOVLX( AccessWork(get_qmtx_kstart(smat)), SSMTRC)

 !call GTOVLX( PW(get_qmtx_kstart(smat)), SSMTRC, PW(get_kfree()), get_lfree() )
 !call GTOVLX( PW(get_qmtx_kstart(smat)), 1.0D0, PW(get_kfree()), get_lfree() )

 if (get_iprint_qmtx() >= 5) then
   write(LUPRI,*) 'GetUnpckOvlpMtx: going to print info on smat!'
   call print_qmtx_attrib(smat)
   call QMtxInfo(smat)
 endif

  call QEXIT('GetUnpckOvlpMtx')

end subroutine GetUnpckOvlpMtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine CAP_STUFF(WORK,LWORK)
!subroutine CAP_STUFF
!---------------------------------------------------------------------------
!
!
!   Called from 'main' dirac flow - dirset.F/PAMINP
!
!
!---------------------------------------------------------------------------
!use Vectors
!use Memory
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
#include "mxcent.h"
#include "nuclei.h"
#include "iratdef.h"

INTERFACE 
    SUBROUTINE CAP_DIRAC_AO(icabfo,iobf1_2) 
       integer, allocatable :: icabfo(:), iobf1_2(:)
   END SUBROUTINE CAP_DIRAC_AO
END INTERFACE

!type(Vect) :: icabfo, iobf1_2
integer, allocatable :: icabfo(:), iobf1_2(:)

LOGICAL IS_FILE, IS_OPEN, CAPDIPVEL_SAVE

DIMENSION WORK(LWORK)
!INTEGER(KIND=8/IRAT) :: ILFREE

!  ... control test... 
 !   call test_oof_basis_fu(WORK,LWORK)
 !call set_vect_iprint(IPRCAP)
 !write(LUPRI,*) 'before stick into memory...LWORK=',LWORK
! ... stick the memory  ... 
 !call StickMemory(WORK,LWORK)

   ! ... new testing procedure ...
   if (CAPNEW) call CAP_new_stuff

  ! icabfo = MakeVector('INTE',NBASIS,'CAP-DIRAC correspondence')   
   allocate(icabfo(NBASIS))
   write(LUPRI,*) 'CAP_STUFF: icabfo allocated, size=',NBASIS
   !iobf1_2 = MakeVector('INTE',NBASIS,'AO function types for obtaining proper factors') 
   allocate(iobf1_2(NBASIS))
   write(LUPRI,*) 'CAP_STUFF: iobf1_2 allocated !'

!    CALL CAP_DIRAC_AO(icabfo,iobf1_2)  ! fill vectors

!  .. 1st contact with the PEMAT file
     LUPEMAT=14
     INQUIRE(FILE='PEMAT',EXIST=IS_FILE,OPENED=IS_OPEN)
     IF (IS_FILE) THEN
      CALL OPNFIL(LUPEMAT,'PEMAT','OLD','CAP_STUFF')
     ELSE
      CALL OPNFIL(LUPEMAT,'PEMAT','UNKNOWN','CAP_STUFF')
      CALL NEWLAB('EOFLABEL',LUPEMAT,LUPRI)
     ENDIF

     CAPDIPVEL_SAVE = CAPDIPVEL
     CAPDIPVEL = .FALSE.

     !call print_mem_param()
     !ILFREE=get_lfree();  write(LUPRI,*) '.....ILFREE=',ILFREE

     IF (CAPOVLP) THEN
         write(LUPRI,*) 'calling DOSAJEEV..1'; CALL  FLSHFO(LUPRI)
         !CALL DOSAJEEV(AccessWork(get_kfree()), ILFREE)
         CALL DOSAJEEV(WORK,LWORK)
         CALL CAP_DIRAC_AO(icabfo,iobf1_2) 
           !  ILFREE=get_lfree()
         !CALL CAPADAPT(AccessWork(get_vect_kstart(icabfo)),AccessWork(get_vect_kstart(iobf1_2)),AccessWork(get_kfree()),ILFREE)
         CALL CAPADAPT(icabfo,iobf1_2,WORK,LWORK)
         CALL CAP_OVLP_COMP
         call quit('CAP_STUFF: only integral comparison !')
     ELSE
       IF (.NOT.NOOVLP) THEN
!      .... first compare DIRAC-CAP overlap integrals
         CAPOVLP=.TRUE.
        ! CALL DOSAJEEV(PW(get_kfree()),get_lfree())
         ! ILFREE=get_lfree(); write(LUPRI,*) '.....ILFREE 2=',ILFREE
        ! CALL DOSAJEEV(PW(get_kfree()),ILFREE)
         !CALL DOSAJEEV(AccessWork(get_kfree()),ILFREE)
         CALL DOSAJEEV(WORK,LWORK)
        ! write(LUPRI,*) 'CAP_STUFF: CAP_DIRAC_AO done 0.1'
         CALL CAP_DIRAC_AO(icabfo,iobf1_2) 
        ! write(LUPRI,*) 'CAP_STUFF: CAP_DIRAC_AO done 1'
          
         !CALL CAPADAPT(AccessWork(get_vect_kstart(icabfo)),AccessWork(get_vect_kstart(iobf1_2)),AccessWork(get_kfree()),ILFREE)
         CALL CAPADAPT(icabfo,iobf1_2,WORK,LWORK)

         CALL CAP_OVLP_COMP
       ENDIF
        CAPOVLP=.FALSE.
        IF (CAPDIPVEL_SAVE.AND..NOT.DONRCAP) THEN
          CAPDIPVEL = .TRUE.
         !  write(LUPRI,*) 'calling DOSAJEEV..2'
         ! ILFREE=get_lfree(); write(LUPRI,*) '.....ILFREE 2=',ILFREE
          !CALL DOSAJEEV(AccessWork(get_kfree()),ILFREE)
          CALL DOSAJEEV(WORK,LWORK)
          CALL CAP_DIRAC_AO(icabfo,iobf1_2)
          !CALL CAPADAPT(AccessWork(get_vect_kstart(icabfo)),AccessWork(get_vect_kstart(iobf1_2)),AccessWork(get_kfree()),ILFREE) 
          CALL CAPADAPT(icabfo,iobf1_2,WORK,LWORK) 
          CALL CAP_DIPVEL_COMP
          call QUIT('CAP_STUFF: comparison')
        ENDIF
        CAPDIPVEL = .FALSE.
         ! write(LUPRI,*) 'calling DOSAJEEV..3'
         ! ILFREE=get_lfree(); write(LUPRI,*) '.....ILFREE 4=',ILFREE
        !CALL DOSAJEEV(AccessWork(get_kfree()),ILFREE)
        CALL DOSAJEEV(WORK,LWORK)
        CALL CAP_DIRAC_AO(icabfo,iobf1_2) 
       ! write(LUPRI,*) 'CAP_STUFF: CAP_DIRAC_AO done 2'
        !CALL CAPADAPT(AccessWork(get_vect_kstart(icabfo)),AccessWork(get_vect_kstart(iobf1_2)),AccessWork(get_kfree()),ILFREE)
        CALL CAPADAPT(icabfo,iobf1_2,WORK,LWORK)
        CAPDIPVEL = CAPDIPVEL_SAVE
     ENDIF

! at the end run 'new' CAP implementation...
! if (CAPNEW) call CAP_new_stuff

! release the memory
 !call UnStickMemory

CALL QEXIT('CAP_STUFF')
end subroutine CAP_STUFF

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine CAP_DIRAC_AO(icapbfo,iobf1_2)
!---------------------------------------------------------------------------
!
! returns filled  icabfo,iobf1_2
!
! called from prepcap90.F90/CAP_STUFF
!
!---------------------------------------------------------------------------
!use Memory
!use Vectors
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
!type(Vect), intent(inout) :: icapbfo, iobf1_2
integer, intent(inout) :: icapbfo(:), iobf1_2(:)
!integer :: KIOBF1,KIOBF2,KILC1,KILC2
integer, allocatable :: IOBF1(:,:),IOBF2(:,:),ILC1(:),ILC2(:)

 !call mymemget('INTE',KIOBF1,3*NBASIS)
 !call mymemget('INTE',KIOBF2,3*NBASIS)
 !call mymemget('INTE',KILC1,NUCDEP)
 !call mymemget('INTE',KILC2,NUCDEP)
 allocate(IOBF1(3,NBASIS));  !call izero(IOBF1,3*NBASIS)
 allocate(IOBF2(3,NBASIS));  !call izero(IOBF2,3*NBASIS)
 allocate(ILC1(NUCDEP));     call izero(ILC1,NUCDEP)
 allocate(ILC2(NUCDEP));     call izero(ILC2,NUCDEP)

write(LUPRI,*) "IOBF1, IOBF2....NBASIS=",NBASIS
write(LUPRI,*) "ILC1, ILC2....NUCDEP=",NUCDEP

 !CALL IZERO(AccessWork(KIOBF1),3*NBASIS)
 !CALL IZERO(AccessWork(KIOBF2),3*NBASIS)
 !CALL IZERO(AccessWork(KILC1),NUCDEP)
 !CALL IZERO(AccessWork(KILC2),NUCDEP)

 !CALL GET_DIR_CAP_AO_MATCH( AccessWork(get_vect_kstart(icapbfo)),AccessWork(get_vect_kstart(iobf1_2)),AccessWork(KIOBF1), AccessWork(KIOBF2),AccessWork(KILC1),AccessWork(KILC2) )

  IOBF1(1,1)=0; IOBF2(1,1)=0
! this routine hurts allocated arrays !
 CALL GET_DIR_CAP_AO_MATCH( icapbfo,iobf1_2,IOBF1, IOBF2,ILC1,ILC2 )

 deallocate(ILC2)
 deallocate(ILC1) ! some problems with the deallocation
 deallocate(IOBF2)
 deallocate(IOBF1)

 call quit('stop here...')

 !call mymemrel('CAP_DIRAC_AO',KILC2,KILC2)
 !call mymemrel('CAP_DIRAC_AO',KILC1,KILC1)
 !call mymemrel('CAP_DIRAC_AO',KIOBF2,KIOBF2)
 !call mymemrel('CAP_DIRAC_AO',KIOBF1,KIOBF1)
 !write(LUPRI,*) 'CAP_DIRAC_AO: after mymemrel !'

! if (get_vect_iprint() >= 4) then
!   print out icapbfo and iobf1_2
  ! call print_vector(icapbfo); 
  ! call print_vector(iobf1_2)
! endif

end subroutine CAP_DIRAC_AO
