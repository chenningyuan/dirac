!dirac_copyright_start
!      Copyright (c) 2013 by the authors of DIRAC.
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
 function operator_on_dirac_list(oper_indx) result(is)
!---------------------------------------------------------------------------
!
!
! Check if this operator is present in the DIRAC list of operators
!
!
!---------------------------------------------------------------------------
 ! integer (kind=4), intent(in) :: oper_indx
  integer , intent(in) :: oper_indx
  logical (kind=1) :: is

    is = .false.

 end function operator_on_dirac_list

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine add_operator(self,oper_indx_temp)
!---------------------------------------------------------------------------
!
!   add operator to the wrapping structure
!
!---------------------------------------------------------------------------
 type(DirOper), intent(inout):: self
 integer, intent(in) :: oper_indx_temp
 !integer (kind=4)  :: oper_indx
 integer :: oper_indx
 logical (kind=1) :: is

 call QENTER('Operators:add_operator')

  oper_indx = oper_indx_temp

 if (.not.operator_on_dirac_list(oper_indx)) then

   self%is_on_list = .TRUE.
 else
   write(LUPRI,*) 'add_operator: operator',oper_indx,' is not in the DIRAC list!'
   self%is_on_list = .FALSE.
 endif

 call QEXIT('Operators:add_operator')
 end subroutine add_operator


 END MODULE Operators

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! subroutine get_prop_saao(INDXPR,propmat,nz_in)
 subroutine get_prop_saao(INDXPR,propmat)
!----------------------------------------------------------------
!
! general purpose F90 routine; 
!
! return the SA-AO unpacked matrix "propmat"
!  of the property operator INDXPR
!
! TODO: verify presence of INDXPR on the DIRAC operator list !
!
! Written by Miro ILIAS, Tel Aviv, 2007
!
!----------------------------------------------------------------
      use QMatrixes
      use Memory
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
!integer(kind=4) :: KPTRI, KONEIN, KFIRST
integer :: KPTRI, KONEIN, KFIRST
!integer(kind=4) :: IKFREE, ILFREE,IKSTART
integer :: IKFREE, ILFREE,IKSTART

 call QENTER('get_prop_saao')

!     Extract symmetry information - ISYM, IREP, ITIM
!  =====================================================
!
  ISYM = IPRPSYM(INDXPR)
  IREP = ISYM - 1
  ITIM = IPRPTIM(INDXPR)

!  write(LUPRI,*) 'property index=',INDXPR,' IREP,ITIM=',IREP,ITIM

! ... allocate the matrix
  propmat = MakeQMtx(set_qmtx_unpck(NFBAS(1,0),NFBAS(2,0)),set_irp(IREP),set_nz(NZ) )

  IPRINT = get_iprint_qmtx()
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

   !mi: do not call, causing mem.error
   !CALL MEMCHK('get_prop_saao',PW,1)

   write(LUPRI,*) 'get_prop_saao: KFREE,LFREE,LWORK=',get_KFREE(),get_LFREE(),get_LWORK()
   write(LUPRI,*) 'get_qmtx_kstart(propmat)=',get_qmtx_kstart(propmat)
   call  print_mem_param
     
  ILFREE = get_LFREE(); IKFREE = get_KFREE(); IKSTART=get_qmtx_kstart(propmat)

  CALL PRPMSAO(INDXPR,AccessWork(get_qmtx_kstart(propmat)),AccessWork(KFIRST), &
               AccessWork(get_KFREE()), ILFREE, IPRINT )
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

