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

!----------------------------------------------------------------------------------
!
!
!  Module for handling the DIRAC memory.
!
!  Contains basic variables and methods to operate upon WORK array.
!  Encapsulates MEMGER/MEMREL/MEMCHK F77 routines.
!
!----------------------------------------------------------------------------------
MODULE Memory

integer, private :: KWORK, KFREE, LFREE, LWORK ! basic variables for memmory handling

! print level
integer, private :: iprmem  = 0

!------------------------------------------------------------------------------
! ... private status 
real (kind=8), private,  pointer :: PW(:) ! points to the big WORK array
!------------------------------------------------------------------------------

! logical variable on if memory is sticked or not
logical (kind=1), private :: is_mem_sticked = .false.

  contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine mymemget( datatype, kbase, length )
!************************************************************************
!
!  public procedure to replace the common MEMGET procedure
!  by hidding private variables KFREE,LFREE
!
! Parameters:  datatype : 'REAL','INTE','LOGI'
!              kbase    : starting position on superarray
!              length   : length of allocated array 
!
!************************************************************************
#include "priunit.h"
character*(*), intent(in) :: datatype
integer, intent(inout) :: kbase
integer, intent(in) :: length
integer :: kbase_temp ! needs to fit into internal MEM_stuff 

 if (iprmem > 25 ) then
   write(LUPRI,'(/,2x,a)') 'mymemget: Before MEMGET:'
   CALL print_mem_param()
 endif

 CALL MEMCHK('mymemget',PW,1) 
  
 kbase_temp = kbase

 ! ... call own memory allocation routine 
 CALL MEMGET(datatype,kbase_temp,length,PW,KFREE,LFREE)

 kbase = kbase_temp

 if (iprmem > 25 ) then
   write(LUPRI,'(/,2x,a)') 'mymemget: After MEMGET:'
   CALL print_mem_param()
 endif

 CALL MEMCHK('mymemget',PW,1) 

end subroutine mymemget

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine mymemrel (datatype, kfirst,krel)
!************************************************************************
!
!   Public procedure to replace common MEMREL subroutine; hides 
! KFREE,LFREE variables.
!
!************************************************************************
#include "priunit.h"
character*(*), intent(in) :: datatype
integer, intent (inout) :: kfirst,krel
integer :: kfirst_temp, krel_temp

 if (iprmem > 25 ) then
   write(LUPRI,'(/,2x,a)') 'mymemrel: Before MEMREL:'
   CALL print_mem_param()
 endif

 ! this substitution is needed due to old DIRAC memory stuff
 kfirst_temp = kfirst; krel_temp  = krel

 call MEMREL(datatype,PW,kfirst_temp,krel_temp,KFREE,LFREE)

 kfirst = kfirst_temp;  krel  =  krel_temp

 if (iprmem > 25 ) then
   write(LUPRI,'(/,2x,a)') 'mymemrel: After MEMREL:'
   CALL print_mem_param()
 endif

end subroutine mymemrel
 
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_mem_param 
!************************************************************************
!
!         print out the DIRAC memory parametersa
!
!************************************************************************
#include "priunit.h"
write(LUPRI,*) '*********************************************************************'
write(LUPRI,'(3(A,I12,/),A,I12)')                &
   '  Actual memory parameters:   KWORK=',KWORK, &
   '----------------------------  KFREE=',KFREE, &
   '                              LFREE=',LFREE, &
   '                              LWORK=',LWORK 
if (associated(PW)) then
   !write(LUPRI,'(A)') '----> Pointer PW is associated with the WORK array.'

   if (size(PW) /= LWORK) then
     write(LUPRI,'(2x,a,2i9)')  &
    ' size(PW) vs LWORK:',size(PW),LWORK
     call quit('print_mem_param: LWORK .ne. size(PW) !')
   endif

   CALL MEMCHK('print_mem_param',PW,1)
  ! CALL MEMCHK('print_mem_param',PW(KFREE),LFREE) 
   write(LUPRI,'(A)') '----> Pointer PW is associated with the WORK array, memory checked.'
else
   write(LUPRI,'(A)') '----> Pointer PW is NOT associated with the WORK array.'
endif
write(LUPRI,*) '*********************************************************************'
call FLSHFO(LUPRI)
end subroutine print_mem_param

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_iprmem(i)
!------------------------------------------------------------------------
!
!     set up the private module printing parameter
!  
! Written by:  Miro ILIAS, Tel Aviv 2007,2008  
!  
!------------------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: i
if (i<0) call QUIT('set_iprmem: print level must be greater than 0 !')
iprmem = i
if (iprmem >= 8) then
  write(LUPRI,'(/,2x,a,i4,/)')  & 
  '**** memory set_iprmem: iprmem set to ',iprmem
endif
end subroutine set_iprmem

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE StickMemory(WORK,LWORK_IN)
!------------------------------------------------------------------------
!
! Stick into the WORK array as the primary source of the memory
!
!                On input: WORK(LWORK_IN)
!
! Note: we use LWORK_IN because LWORK is the private parameter of the module.
!
!------------------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: LWORK_IN
real*8, intent(in), target :: WORK(LWORK_IN)
if (iprmem > 10) then
  call HEADER('Output from StickMemory',-1)
endif

if (associated(PW)) then
  if (iprmem > 10) then
    write(LUPRI,*) 'StickMemory: PW is already associated - Error Exit !'
    call FLSHFO(LUPRI)
  endif
  CALL QUIT('StickMemory: PW is already associated !')
else
  if (iprmem > 10) then
    write(LUPRI,*) 'StickMemory: PW is not associated'
    call FLSHFO(LUPRI)
  endif
endif
! ..................................................
PW => WORK ! PW points now to the big WORK array ... 
! ..................................................
is_mem_sticked = .true. ! update logical variable

LWORK = LWORK_IN ! actualize module variable

if (size(PW) /= LWORK_IN) then
! if (iprmem > 10) then
   write(LUPRI,*) ' size(PW)=',size(PW), &
                  ' while entering LWORK=',LWORK_IN
   write(LUPRI,*) 'somme error in array sticking !'
   call FLSHFO(LUPRI)
! endif
 CALL QUIT('StickMemory: error in length !')
endif

! initialize the memmory
#include "memint.h"

! check it
CALL MEMCHK('StickMemory',PW,1)

if (iprmem > 10 ) then
  CALL print_mem_param()
  write(LUPRI,'(/,2x,a,/)') & 
  '*** StickMemory: PW sticked into WORK, memmory checked ***'
  call FLSHFO(LUPRI)
endif
END SUBROUTINE StickMemory

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION  AccessWork(ipos,end_pos) result(pw_out)
!--------------------------------------------------------------
!
! returns pointer to the chunk of array starting from WORK(ipos)
!
!   Function used as parameter for subroutine/function 
! instead of direct PW(i)
!
!--------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: ipos
integer, intent(in), optional :: end_pos
integer :: this_end_pos
real*8, pointer, dimension(:) :: pw_out 

! write(LUPRI,*) 'in AccessWork...'

if (.not.present(end_pos)) then 
   this_end_pos = -1
else
   this_end_pos = end_pos
endif

! ... check ipos ...  
if ( CheckArrayPos(ipos=ipos,end_pos=this_end_pos) ) then
  pw_out => PW(ipos:)
else
  write(LUPRI,'(2x,a)')  'AccessWork error !'
  call print_mem_param
  call quit('AccessWork: wrong array ! Error exit')
endif

!write(LUPRI,*) 'in AccessWork...end'
!call FLSHFO(LUPRI)

END FUNCTION  AccessWork

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION  AccessWorkElem(ipos) result(pwi_out)
!--------------------------------------------------------------
!
! returns pointer of the element WORK(ipos)
!
! returning pointer enables change of WORK(ipos) elem
!
! Used as parameter for subroutine/function instead of direct PW(i)
!
!--------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: ipos
real*8, pointer :: pwi_out 

! ... check ipos ...  
if ( CheckArrayPos(ipos) ) then
  pwi_out => PW(ipos)
else
  write(LUPRI,'(/,2x,a)') 'AccessWorkElem: error ! '
  call print_mem_param
  call quit('AccessWorkElem: bad array !')
endif

END FUNCTION  AccessWorkElem


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION GetWorkElem(ipos) result(WorkElem)
!--------------------------------------------------------------
!
!  returns element from the big Work array
!  
!--------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: ipos
real*8 :: WorkElem

if ( CheckArrayPos(ipos) ) then

 WorkElem = PW(ipos)

else
 call print_mem_param
 call quit('GetWorkElem: wrong array or pos !')
endif

END FUNCTION GetWorkElem

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE  PutWorkElem(ipos,WorkElem)
!--------------------------------------------------------------
!
! put real*8 element into WORK(ipos) 
!
!--------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: ipos
real*8, intent(in) :: WorkElem

if ( CheckArrayPos(ipos) ) then

 PW(ipos) = WorkElem

else
   call print_mem_param
   call quit('PutWorkElem: error !')
endif


END SUBROUTINE  PutWorkElem

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION  CheckArrayPos(ipos, end_pos) result(is_good)
!--------------------------------------------------------------------
!
! Checks array WORK (through poiner PW) and position on it "ipos"
!
!--------------------------------------------------------------------
#include "priunit.h"
integer, intent(in) :: ipos
integer, intent(in), optional :: end_pos
logical :: is_good

 if (ipos <= LFREE .and. ipos >= 1 .and. associated(PW) .and. is_mem_sticked ) then

   if (size(PW) /= LWORK) then
     write(LUPRI,'(2x,a,2i9)')  &
     'CheckArrayPos:  size(PW) vs LWORK:',size(PW),LWORK
     call FLSHFO(LUPRI)
     call quit('CheckArrayPos: LWORK .ne. size(PW) !')
   endif

   if (present(end_pos)) then
     if (end_pos > LFREE) then
       write(LUPRI,*) 'CheckArrayPos: LFREE, end_pos=',LFREE, end_pos
       call FLSHFO(LUPRI)
       call quit('CheckArrayPos: end_pos > LFREE !')
     endif
   endif

   ! check the memory ...
   CALL MEMCHK('CheckArrayPos',PW,1)

   is_good = .true.
 else
   is_good = .false.
 endif

END FUNCTION  CheckArrayPos

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
SUBROUTINE UnStickMemory
!--------------------------------------------------------------
!
!      release the PW pointer sticked into WORK array
!
!--------------------------------------------------------------
#include "priunit.h"
integer, parameter :: print_default = 10
if (associated(PW)) then
  if ( iprmem > print_default ) then
    write(LUPRI,*) 'The PW pointer is already associated, all right'
  endif
else
   write(LUPRI,'(/,2x,a,/)') 'UnStickMemory Error: The PW pointer is NOT associated !'
   call quit('UnStickMemory Error: The PW pointer is NOT associated !')
endif
!  ...
nullify(PW);
is_mem_sticked = .false. ! update logical variable
!  ... 
if ( iprmem > print_default ) then
  write(LUPRI,'(/,A)') 'UnStickMemory: PW pointer was deassociated!'
  call FLSHFO(LUPRI)
endif
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
 is_sticked = is_mem_sticked
END FUNCTION IsMemorySticked

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION get_kfree() result(kfree_out)
!------------------------------------------------------------------------
! returns the (private) KFREE value
!------------------------------------------------------------------------
#include "priunit.h"
integer :: kfree_out
 kfree_out = KFREE
 if (iprmem >= 15) then
  write(LUPRI,*) 'get_kfree: KFREE=',KFREE,'= returning kfree_out:',kfree_out
  call FLSHFO(LUPRI)
 endif
END FUNCTION get_kfree

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION get_lfree() result(lfree_out)
!------------------------------------------------------------------------
!
!    returns the (protected) LFREE value
!
!------------------------------------------------------------------------
#include "priunit.h"
integer :: lfree_out
 lfree_out = LFREE
 if (iprmem >= 15) then
   write(LUPRI,'(2x,a,i8,a,i8)')  &
   'Memory/get_lfree: LFREE:',LFREE,' = lfree_out:',lfree_out
   call FLSHFO(LUPRI)
 endif
END FUNCTION get_lfree

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
FUNCTION get_lwork() result(lwork_out)
!------------------------------------------------------------------------
!
!         returns the (module protected) LWORK value
!
!------------------------------------------------------------------------
#include "priunit.h"
integer :: lwork_out
 lwork_out = LWORK
 if (iprmem >= 15) then
  write(LUPRI,*) 'get_lwork: LWORK=',LWORK,'lwork_out=',lwork_out
  call FLSHFO(LUPRI)
 endif
END FUNCTION get_lwork

END MODULE Memory
