!=======================================================================
!
! file contains pure Fortran90 classes related to CAP project
!
!=======================================================================

  module cap

! private ! must be public that subroutines can be used in F77 parts !
#include "dcbcap.h"

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

 subroutine cap_ovlp
use Memory
use cap
#include "priunit.h"
! check if not make check of DIRAC/CAP overlap integrals
 CALL QENTER('cap_ovlp')

 write(LUPRI,*) 'cap_ovlp, CAPOVLP=',CAPOVLP,' IPRCAP=',IPRCAP
 call print_mem_param
 
 if (CAPOVLP) then
  call HEADER('Output form cap_ovlp',-1)

  call cap_ovlp_comp

  call quit('cap_ovlp')
 endif
 CALL QEXIT('cap_ovlp')
 end subroutine cap_ovlp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine cap_ovlp_comp
!---------------------------------------------------------------------------
!
!  compare overlap integrals from DIRAC and from Sajeev's CAP integral
! generator
!
!----------------------------------------------------------------------------
   use QMatrixes
#include "implicit.h"
#include "priunit.h"
type(QMtx) :: smat,capovlp

call QENTER('cap_ovlp_comp')

 call GetUnpckOvlpMtx(smat)

 call GetCapOvlp(capovlp)

!  compare smat & capovlp quaternion matrixes
   diff = getdiff_qmtx(smat,capovlp)
   write(LUPRI,*) 
   write(LUPRI,*) 'comparison of 2 matrixes: ',diff
   write(LUPRI,*) 

  call  AddTextQMtx(smat,'difference smat - capovlp ')
 ! call  AddTextQMtx(smat,'')
  call substract_qmtx(smat,capovlp)
  call QMTxInfo(smat)

call QEXIT('cap_ovlp_comp')
end subroutine cap_ovlp_comp

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine GetCapOvlp(capovlp)
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
!#include "dcbcap.h"
#include "dcbham.h"
 type(QMtx), intent(inout) :: capovlp
call QENTER('GetCapOvlp')

 call get_prop_saao(IPCAPR,capovlp,1) 
! call set_qmtx_nz(capovlp,1)
 call AddTextQMtx(capovlp,'CAP overlap integrals')

 call QMtxInfo(capovlp)

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

integer, intent(in) :: INDXPR
integer(kind=8/IRAT),intent(in),optional :: nz_in
!integer(kind=4),intent(in),optional :: nz_in
!integer,intent(in),optional :: nz_in

type(QMtx),intent(inout) :: propmat

logical :: is_square
integer(kind=4) :: KPTRI, KONEIN, KFIRST
integer(kind=8/IRAT) :: NZX

 call QENTER('get_prop_saao')

!
!     Extract symmetry information - ISYM, IREP, ITIM
!     =============================
!
      ISYM = IPRPSYM(INDXPR)
      IREP = ISYM - 1
      ITIM = IPRPTIM(INDXPR)

     write(LUPRI,*) 'property index=',INDXPR,' IREP,ITIM=',IREP,ITIM

    if (present(nz_in)) then
      NZX = nz_in
    else
      NZX = NZ
    endif

!         allocate array for matrix
   ! propmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(IREP),set_nz(NZX) )
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
       CALL PRPMSAO(INDXPR,IDIMENS,DUMMY,PW(get_qmtx_kstart(propmat)),PW(KONEIN),PW(KIRST), &
                    PW(get_KFREE),get_LFREE,IPRINT)
     ELSE
       CALL PRPMSAO(INDXPR,IDIMENS,PW(KPTRI),PW(get_qmtx_kstart(propmat)),PW(KONEIN),PW(KFIRST), &
                     PW(get_KFREE),get_LFREE,IPRINT )
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
   use Memory
#include "priunit.h"
#include "dcbbas.h"
#include "dcbham.h"
 type(QMtx), intent(inout) :: smat
call QENTER('GetUnpckOvlpMtx')

 call set_iprqmtx(10)
 call set_iprmem(10)
 smat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(0),set_nz(1), & 
                 'DIRAC unpacked overlap SA-AO matrix')
! call QMtxInfo(smat)
! call print_mem_param()

 write(LUPRI,*) 'N2BAST=',N2BAST
 write(LUPRI,*) 'get_kfree, get_lfree=',get_kfree(),get_lfree()
 write(LUPRI,*) 'get_qmtx_kstart(smat)=',get_qmtx_kstart(smat)
 write(LUPRI,*) 'SSMTRC=',SSMTRC

! call lower level routine
 call GTOVLT(PW(get_qmtx_kstart(smat)), SSMTRC, PW(get_kfree), get_lfree)

 call QMtxInfo(smat)

  call QEXIT('GetUnpckOvlpMtx')
 end subroutine GetUnpckOvlpMtx
