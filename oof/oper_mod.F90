!
!
!
 MODULE Operators

  use mystring

 type DirOper
  private
 ! ... index of the operator
  integer (kind=4) :: indxpr=-1
 ! ... is on the operator list ?
  logical (kind=1) :: is_on_list=.false. 
 ! ... attached text
  type(string), dimension(1) :: operator_text 
 end type DirOper
 
 contains
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! function is_operator_on_list(self) result(is)
 function is_operator_on_list(oper_indx) result(is)
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
 ! type(DiOper), intent(in) :: self
  integer (kind=4), intent(in) :: oper_indx
  logical (kind=1) :: is
 call QENTER('Operators:is_operator_on_list')

 call QEXIT('Operators:is_operator_on_list')
 end function is_operator_on_list



 END MODULE Operators

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! subroutine get_prop_saao(INDXPR,propmat,nz_in)
 subroutine get_prop_saao(INDXPR,propmat)
!----------------------------------------------------------------
!
!
! general purpose F90 routine; 
!
! return the SA-AO unpacked matrix "propmat" of the property INDXPR
!
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


