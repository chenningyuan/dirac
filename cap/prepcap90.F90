!=======================================================================
!
! file contains pure Fortran90 classes related to CAP project
!
!=======================================================================

  module cap

! private ! must be public that subroutines can be used in F77 parts !
!#include "dcbcap.h"

!  logical, private :: cap_ovlp
!   integer(kind=4), private :: iprint
!   integer(kind=4), private :: iprcap
!  type cap_int_param
!  end type cap_int_param
  contains

  subroutine cap_init
use QMatrixes

#include "implicit.h"
#include "priunit.h"

!#include "dcbcap.h"

! iprcap = IPRCAP
! cap_ovlp = CAP_OVLP

  if (IPRCAP >= 2) then
   write(LUPRI,*) 'in init routine !'
   write(LUPRI,*) 'module iprcap =',IPRCAP
   write(LUPRI,*) 'module cap_ovlp =',CAPOVLP
!   write(LUPRI,*) 'CPRIMU ',CPRIMU
!   write(LUPRI,*) 'ALPHA: ',ALPHA
!   write(LUPRI,*) 'NUC: ',NUC

  endif

! call quit('in cap_init...')

   end subroutine cap_init

 end module cap

subroutine test_files
use QMatrixes

#include "implicit.h"
#include "priunit.h"

#include "mxcent.h"
#include "maxorb.h"

#include "aovec.h"
#include "maxaqn.h"
#include "primit.h"
!C#include "cbirea.h"
#include "chrxyz.h"
#include "nuclei.h"
#include "symmet.h"
#include "dummy.h"
#include "infpar.h"

!  write(LUPRI,*) 'CPRIMU ',CPRIMU
  write(LUPRI,*) 'ALPHA: ',ALPHA
!   write(LUPRI,*) 'NUC: ',NUC

  call quit('end in test files...')
end subroutine test_files

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
 endif

call QEXIT('cap_ovlp_comp')
end subroutine cap_ovlp_comp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
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

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! subroutine get_prop_saao(INDXPR,propmat,nz_in)
 subroutine get_prop_saao(INDXPR,propmat)
!----------------------------------------------------------------
!
!   general purpose routine; 
!
! get the SA-AO unpacked matrix "propmat" of the property INDXPR
!
!----------------------------------------------------------------
      use QMatrixes
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbxpr.h"
#include "dcbprl.h"
#include "dummy.h"

type(QMtx),intent(inout) :: propmat

integer, intent(in) :: INDXPR
!integer(kind=8/IRAT),intent(in),optional :: nz_in
!integer(kind=4),intent(in),optional :: nz_in
!integer,intent(in),optional :: nz_in

logical :: is_square
integer(kind=4) :: KPTRI, KONEIN, KFIRST

 call QENTER('get_prop_saao')

!     Extract symmetry information - ISYM, IREP, ITIM
!  =====================================================
!
  ISYM = IPRPSYM(INDXPR)
  IREP = ISYM - 1
  ITIM = IPRPTIM(INDXPR)

!  write(LUPRI,*) 'property index=',INDXPR,' IREP,ITIM=',IREP,ITIM

  propmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(IREP),set_nz(NZ) )

  IPRINT = get_iprint_qmtx()
  !IPRINT = 12

  if (ITIM == 0) then 
      is_square = .true.
      IDIMENS = N2BBASX
      CALL MYMEMGET('REAL',KPTRI  ,NZ*IDIMENS)
      CALL MYMEMGET('REAL',KONEIN ,   IDIMENS)
  else
      is_square = .false.
      IDIMENS = NNBBASX
      CALL MYMEMGET('REAL',KPTRI  ,NZ*IDIMENS)
      CALL MYMEMGET('REAL',KONEIN ,   IDIMENS) 
  endif

!   check it if it's not needed for property !!!
  CALL MYMEMGET('LOGI',KFIRST ,NZ        )

  IF (IPRINT >= 7) THEN
      write(LUPRI,*) 'get_prop_saao  memmory parameters'
     call print_mem_param
  ENDIF

#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED', &
              FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED', &
              FILE='AOPROPER')
#endif

     IF (IS_SQUARE) THEN
! call 'old' F77 routines
       CALL PRPMSAO(INDXPR,IDIMENS,DUMMY,PW(get_qmtx_kstart(propmat)),PW(KONEIN),PW(KIRST), &
                    PW(get_KFREE()),get_LFREE,IPRINT)
     ELSE
       CALL PRPMSAO(INDXPR,IDIMENS,PW(KPTRI),PW(get_qmtx_kstart(propmat)),PW(KONEIN),PW(KFIRST), &
                    PW(get_KFREE()),get_LFREE,IPRINT )
     ENDIF

     CLOSE(LU1INT,STATUS = 'KEEP')

     IF (IPRINT >= 7) THEN
      write(LUPRI,*) 'get_prop_saao 1 memmory parameters'
      call print_mem_param
     ENDIF
!     release all the memmory
     call mymemrel('get_prop_saao',KFIRST,KFIRST)
!     call print_mem_param
     call mymemrel('get_prop_saao',KONEIN,KONEIN)
!     call print_mem_param
     call mymemrel('get_prop_saao',KPTRI,KPTRI)

     IF (IPRINT >= 7) THEN
      write(LUPRI,*) 'get_prop_saao 2 memmory parameters'
      call print_mem_param
     ENDIF

 call QEXIT('get_prop_saao')
 end subroutine get_prop_saao

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GetUnpckOvlpMtx(smat)
!===========================================================================
!
!
!    get unpacked overlap matrix, "smat" - allocate inside, then return
!
!
!===========================================================================
   use QMatrixes
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
   write(LUPRI,*) 'SSMTRC=',SSMTRC, ' but I set 1.0'
 endif

! call lower level routine
 !call GTOVLX( PW(get_qmtx_kstart(smat)), SSMTRC, PW(get_kfree()), get_lfree() )
 call GTOVLX( PW(get_qmtx_kstart(smat)), 1.0D0, PW(get_kfree()), get_lfree() )

 if (get_iprint_qmtx() >= 5) then
   call print_qmtx_attrib(smat)
   call QMtxInfo(smat)
 endif

  call QEXIT('GetUnpckOvlpMtx')
end subroutine GetUnpckOvlpMtx

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine CAP_STUFF(WORK,LWORK)
!---------------------------------------------------------------------------
!
!
!    Called from 'main' dirac flow - dirset.F/PAMINP
!
!
!---------------------------------------------------------------------------
use Vectors
use Memory
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
#include "mxcent.h"
#include "nuclei.h"

type(Vect) :: icabfo, iobf1_2
!LOGICAL*1 IS_FILE,IS_OPEN
LOGICAL IS_FILE,IS_OPEN
DIMENSION WORK(LWORK)

 CALL QENTER('CAP_STUFF')

 call set_vect_iprint(IPRCAP)

 call StickMemory(WORK,LWORK)

     icabfo = MakeVector('INTE',NBASIS,'CAP-DIRAC correspondence')   
     iobf1_2 = MakeVector('INTE',NBASIS,'AO function types for obatining factors') 

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

     IF (CAPOVLP) THEN
         CALL DOSAJEEV(PW(get_kfree()),get_lfree())
         CALL CAP_DIRAC_AO(icabfo,iobf1_2) 
         CALL CAPADAPT(PW(get_vect_kstart(icabfo)),PW(get_vect_kstart(iobf1_2)),PW(get_kfree()),get_lfree())
         CALL CAP_OVLP_COMP
         call quit('CAP_STUFF: only integral comparison !')
     ELSE
       IF (.NOT.NOOVLP) THEN
!      .... first compare DIRAC-CAP overlap integrals
         CAPOVLP=.TRUE.
         CALL DOSAJEEV(PW(get_kfree()),get_lfree())
         CALL CAP_DIRAC_AO(icabfo,iobf1_2) 
         CALL CAPADAPT(PW(get_vect_kstart(icabfo)),PW(get_vect_kstart(iobf1_2)),PW(get_kfree()),get_lfree())
         CALL CAP_OVLP_COMP
       ENDIF
        CAPOVLP=.FALSE.
        CALL DOSAJEEV(PW(get_kfree()),get_lfree())
        CALL CAP_DIRAC_AO(icabfo,iobf1_2) 
        CALL CAPADAPT(PW(get_vect_kstart(icabfo)),PW(get_vect_kstart(iobf1_2)),PW(get_kfree()),get_lfree())
     ENDIF

! release the memory
 call UnStickMemory

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
use Memory
use Vectors
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
type(Vect), intent(inout) :: icapbfo, iobf1_2
integer(kind=4) :: KIOBF1,KIOBF2,KILC1,KILC2

  CALL QENTER('CAP_DIRAC_AO')

 call mymemget('INTE',KIOBF1,3*NBASIS)
 call mymemget('INTE',KIOBF2,3*NBASIS)
 call mymemget('INTE',KILC1,NUCDEP)
 call mymemget('INTE',KILC2,NUCDEP)

 CALL IZERO(PW(KIOBF1),3*NBASIS)
 CALL IZERO(PW(KIOBF2),3*NBASIS)
 CALL IZERO(PW(KILC1),NUCDEP)
 CALL IZERO(PW(KILC2),NUCDEP)

 CALL GET_DIR_CAP_AO_MATCH(PW(get_vect_kstart(icapbfo)),PW(get_vect_kstart(iobf1_2)),PW(KIOBF1), &
           PW(KIOBF2),PW(KILC1),PW(KILC2))

 call mymemrel('CAP_DIRAC_AO',KILC2,KILC2)
 call mymemrel('CAP_DIRAC_AO',KILC1,KILC1)
 call mymemrel('CAP_DIRAC_AO',KIOBF2,KIOBF2)
 call mymemrel('CAP_DIRAC_AO',KIOBF1,KIOBF1)

 if (get_vect_iprint() >= 4) then
!   print out icapbfo and iobf1_2
   call print_vector(icapbfo); 
   call print_vector(iobf1_2)
 endif

  CALL QEXIT('CAP_DIRAC_AO')

end subroutine CAP_DIRAC_AO

