!
! Module contains basic variables of the DIRAC memmory array
!
MODULE Memory

!integer :: KWORK, KFREE, LFREE, LWORK ! basic variables of memmory handling
integer, private :: KWORK, KFREE, LFREE, LWORK ! basic variables of memmory handling

! print level
integer, private :: iprmem  = 0

! pointer to the working array; public status is needed !
real (kind=8),  pointer :: PW(:) 
!real (kind=8), private,  pointer :: PW(:) ! points to the big WORK array

! logical variable on if memory is sticked or not
logical (kind=1), private :: is_mem_sticked = .false.

  contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine mymemget( datatype, kbase, length )
!************************************************************************
!
!  public procedure to replace the common MEMGET procedure
! with hidding of private variables KFREE,LFREE
!
! Parameters:  datatype : 'REAL','INTE','LOGI'
!              kbase    : starting position on superarray
!              length   : length of allocated array 
!
!************************************************************************
#include "priunit.h"
#include "iratdef.h"
character*(*), intent(in) :: datatype
integer(kind=4), intent(inout) :: kbase
!integer(kind=8/IRAT), intent(inout) :: kbase
!integer, intent(inout) :: kbase
!integer (kind=8), intent(in) :: length
integer (kind=8/IRAT), intent(in) :: length
integer :: kbase_temp ! needs to fit into internal MEM_stuff 
call QENTER('mymemget')

 CALL MEMCHK('mymemget',PW,1) 
  
 kbase_temp = kbase

 CALL MEMGET(datatype,kbase_temp,length,PW,KFREE,LFREE)

 kbase = kbase_temp

 CALL MEMCHK('mymemget',PW,1) 

call QEXIT('mymemget')
end subroutine mymemget

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine mymemrel (datatype, kfirst,krel)
!************************************************************************
!
!  Public procedure to replace common MEMREL subroutine to hide 
! KFREE,LFREE variables
!
!************************************************************************
#include "priunit.h"
character*(*), intent(in) :: datatype
integer(kind=4), intent (inout) :: kfirst,krel
integer :: kfirst_temp, krel_temp

call QENTER('mymemrel')
 ! this substitution is needed due to old DIRAC memory stuff
 kfirst_temp = kfirst; krel_temp  = krel
 call MEMREL(datatype,PW,kfirst_temp,krel_temp,KFREE,LFREE)
 kfirst = kfirst_temp;  krel  =  krel_temp
call QEXIT('mymemrel')
end subroutine mymemrel
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_mem_param 
!************************************************************************
!
!        print out the DIRAC memory parameters
!
!************************************************************************
#include "priunit.h"
CALL QENTER('print_mem_param')
write(LUPRI,*) '*********************************************************************'
write(LUPRI,'(3(A,I12,/),A,I12)')                &
   '  Actual memory parameters:   KWORK=',KWORK, &
   '----------------------------  KFREE=',KFREE, &
   '                              LFREE=',LFREE, &
   '                              LWORK=',LWORK 
if (associated(PW)) then
   !write(LUPRI,'(A)') '----> Pointer PW is associated with the WORK array.'
   if (size(PW) /= LWORK) then
     write(LUPRI,*) ' size(PW), LWORK:',size(PW), LWORK
     call quit('print_mem_param: LWORK .ne. size(PW) !')
   endif
   CALL MEMCHK('print_mem_param',PW,1)
  ! CALL MEMCHK('print_mem_param',PW(KFREE),LFREE) 
   write(LUPRI,'(A)') '----> Pointer PW is associated with the WORK array, memory checked.'
else
   write(LUPRI,'(A)') '----> Pointer PW is NOT associated with the WORK array.'
endif
write(LUPRI,*) '*********************************************************************'
CALL QEXIT('print_mem_param')
end subroutine print_mem_param

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_iprmem(i)
!------------------------------------------------------------------------
!
!                  set up the module printing parameter
!
!------------------------------------------------------------------------
integer, intent(in) :: i
call QENTER('set_iprmem')
if (i<0) call QUIT('set_iprmem: print level must be greater than 0 !')
iprmem = i
call QEXIT('set_iprmem')
end subroutine set_iprmem

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE StickMemory(WORK,LWORK_IN)
!------------------------------------------------------------------------
!
! Stick into the WORK array as the primary source of the memory
!
!                     On input: WORK(LWORK_IN)
!
!------------------------------------------------------------------------
#include "priunit.h"
#include "iratdef.h"
!integer*4, intent(in) :: LWORK_IN
integer(kind=8/IRAT), intent(in) :: LWORK_IN
real*8, intent(in), target :: WORK(LWORK_IN)
CALL QENTER('StickMemory')
if (iprmem > 10) then
  call HEADER('Output from StickMemory',-1)
endif

if (associated(PW)) then
  if (iprmem > 10) then
    write(LUPRI,*) 'PW is already associated - Exit !'
  endif
  CALL QUIT('PW is already associated !')
else
  if (iprmem > 10) then
    write(LUPRI,*) 'PW is not associated'
  endif
endif
! ..................................................
PW => WORK ! PW points now to the big WORK array ... 
! ..................................................
is_mem_sticked = .true. ! update logical variable
LWORK = LWORK_IN ! actualize module variable
if (size(PW) /= LWORK_IN) then
 if (iprmem > 10) then
   write(LUPRI,*) ' size(PW)=',size(PW), &
                  ' while entering LWORK=',LWORK_IN
   write(LUPRI,*) 'somme error in array sticking !'
 endif
 CALL QUIT('StickMemory: error in length !')
endif

! initialize the memmory
#include "memint.h"

! check it
CALL MEMCHK('StickMemory',PW,1)

if (iprmem > 10 ) then
  CALL print_mem_param()
  write(LUPRI,*) 'StickMemory: PW sticked into WORK, memmory checked !'
endif
CALL QEXIT('StickMemory')
END SUBROUTINE StickMemory

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE UnStickMemory
!--------------------------------------------------------------
!
! release the PW pointer sticked into WORK array
!
!--------------------------------------------------------------
#include "priunit.h"
integer(kind=4), parameter :: print_default = 10
CALL QENTER('UnStickMemory')
if (associated(PW)) then
  if ( iprmem > print_default ) then
    write(LUPRI,*) 'Entering PW is already associated'
  endif
else
  if ( iprmem > print_default ) then
    write(LUPRI,*) 'Entering PW is NOT associated !'
  endif
endif
nullify(PW);
is_mem_sticked = .false. ! update logical variable
!  ... 
if ( iprmem > print_default ) then
  write(LUPRI,'(A)') 'PW pointer was deassociated!'
endif
CALL QEXIT('UnStickMemory')
END SUBROUTINE UnStickMemory

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION IsMemorySticked() result(is_sticked)
!------------------------------------------------------------------------
!
! returns the (private) 'is_mem_sticked' value which is indicator
! of the presence of the large memory block PW => WORK
!
!------------------------------------------------------------------------
#include "priunit.h"
logical(kind=1) :: is_sticked
CALL QENTER('IsMemorySticked')
 is_sticked = is_mem_sticked
CALL QEXIT('IsMemorySticked')
END FUNCTION IsMemorySticked

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION get_kfree() result(kfree_out)
!------------------------------------------------------------------------
! returns the (private) KFREE value
!------------------------------------------------------------------------
#include "priunit.h"
!integer :: kfree_out
integer(kind=4) :: kfree_out
CALL QENTER('get_kfree')
 kfree_out = KFREE
 if (iprmem >= 15) then
  write(LUPRI,*) 'get_kfree: KFREE=',KFREE,'kfree_out=',kfree_out
 endif
CALL QEXIT('get_kfree')
END FUNCTION get_kfree

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION get_lfree() result(lfree_out)
!------------------------------------------------------------------------
!
!    returns the (protected) LFREE value
!
!------------------------------------------------------------------------
#include "priunit.h"
!integer :: lfree_out
integer(kind=4) :: lfree_out
CALL QENTER('get_lfree')
 lfree_out = LFREE
 if (iprmem >= 15) then
  write(LUPRI,'(2x,a,i8,a,i8)')  &
  'Memory/get_lfree: LFREE:',LFREE,' = lfree_out:',lfree_out
 endif
CALL QEXIT('get_lfree')
END FUNCTION get_lfree

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION get_lwork() result(lwork_out)
!------------------------------------------------------------------------
!
! return the (module protected) LWORK value
!
!------------------------------------------------------------------------
#include "priunit.h"
!integer :: lwork_out
integer (kind=4) :: lwork_out
CALL QENTER('get_lwork')
 lwork_out = LWORK
 if (iprmem >= 15) then
  write(LUPRI,*) 'get_lwork: LWORK=',LWORK,'lwork_out=',lwork_out
 endif
CALL QEXIT('get_lwork')
END FUNCTION get_lwork

END MODULE Memory
