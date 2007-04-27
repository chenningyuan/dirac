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
!   called from CAP_STUFF
!
!----------------------------------------------------------------------------
   use QMatrixes
#include "implicit.h"
#include "priunit.h"
type(QMtx) :: smat,capovlp
double precision, parameter :: diff_thrs = 1.0D-7
call QENTER('cap_ovlp_comp')

 call GetUnpckOvlpMtx(smat)

 call GetCapOvlp(capovlp)

!  compare smat & capovlp quaternion matrixes
 diff = getdiff_part_qmtx(smat,capovlp)
  write(LUPRI,*) 
  write(LUPRI,*) 'comparison of the 2 matrixes: ',diff
  write(LUPRI,*) 

 if (diff > diff_thrs) then
  call  AddTextQMtx(smat,'******** difference smat - capovlp ***********')
  call substract_part_qmtx(smat,capovlp)
  call QMTxInfo(smat)
 endif

call QEXIT('cap_ovlp_comp')
end subroutine cap_ovlp_comp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GetCapOvlp(capovlp_mtx)
!--------------------------------------------------------------------
!
! get CAP overlap matrix
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

 call get_prop_saao(IPCAP_S,capovlp_mtx,1) 
! call set_qmtx_nz(capovlp,1)
 call AddTextQMtx(capovlp_mtx,'CAP overlap integrals')

 call QMtxInfo(capovlp_mtx)

call QEXIT('GetCapOvlp')
end subroutine GetCapOvlp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine get_prop_saao(INDXPR,propmat,nz_in)
!----------------------------------------------------------------
!
!   general purpose routine; 
!
! get the SA-AO unpacked matrix "propmat" of the property INDXPR
!
! nz_in : optional parameter
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
integer(kind=8/IRAT),intent(in),optional :: nz_in
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

  write(LUPRI,*) 'property index=',INDXPR,' IREP,ITIM=',IREP,ITIM

  propmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(IREP),set_nz(NZ) )

     IPRINT = 20

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

     call print_mem_param

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

     call print_mem_param
    
!     release all the memmory
     call mymemrel('get_prop_saao',KFIRST,KFIRST)
!     call print_mem_param
     call mymemrel('get_prop_saao',KONEIN,KONEIN)
!     call print_mem_param
     call mymemrel('get_prop_saao',KPTRI,KPTRI)

     call print_mem_param

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
!   use Memory
#include "priunit.h"
#include "dcbbas.h"
#include "dcbham.h"
 type(QMtx), intent(inout) :: smat
call QENTER('GetUnpckOvlpMtx')

 call set_iprqmtx(10)
 call set_iprmem(10)
 smat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(0),set_nz(1), & 
                 'DIRAC unpacked overlap SA-AO matrix' )
! call QMtxInfo(smat)
! call print_mem_param()

 write(LUPRI,*) 'N2BAST=',N2BAST
 write(LUPRI,*) 'get_kfree, get_lfree=',get_kfree(),get_lfree()
 write(LUPRI,*) 'get_qmtx_kstart(smat)=',get_qmtx_kstart(smat)
 write(LUPRI,*) 'SSMTRC=',SSMTRC

! call lower level routine
 call GTOVLX( PW(get_qmtx_kstart(smat)), SSMTRC, PW(get_kfree()), get_lfree() )

 call QMtxInfo(smat)

  call QEXIT('GetUnpckOvlpMtx')
end subroutine GetUnpckOvlpMtx

subroutine CAP_STUFF(WORK,LWORK)
!
!
!  Called from dirset.F/PAMINP
!
!
!
use Vectors
use Memory
#include "implicit.h"
#include "priunit.h"
#include "dcbcap.h"
#include "mxcent.h"
#include "nuclei.h"

type(Vect) :: icabfo, iobf1_2

DIMENSION WORK(LWORK)
CALL QENTER('CAP_STUFF')


 call StickMemory(WORK,LWORK)

     icabfo = MakeVector('INTE',NBASIS,'CAP-DIRAC correspondence')   

     iobf1_2 = MakeVector('INTE',NBASIS,'AO function types for obatining factors') 

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

subroutine CAP_DIRAC_AO(icapbfo,iobf1_2)
!
! returns:  icabfo,iobf1_2
!
! called from CAP_STUFF
!
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

  call mymemrel('..',KILC2,KILC2)
  call mymemrel('..',KILC1,KILC1)
  call mymemrel('..',KIOBF2,KIOBF2)
  call mymemrel('..',KIOBF1,KIOBF1)

! print out icapbfo and iobf1_2
   call print_vector(icapbfo); 
   call print_vector(iobf1_2)

CALL QEXIT('CAP_DIRAC_AO')

end subroutine CAP_DIRAC_AO

!
!
!
!

