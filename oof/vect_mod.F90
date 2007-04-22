 MODULE Vectors

!use Memory
! deta type for strings
use mystring

!   print parameter
integer(kind=4), private :: iprint_vect=0

type Vect
private
! contains position on the PW (WORK) array 
  integer (kind=4) :: kstart 
!   the total size on the array
  integer (kind=8) :: array_size = 0
! ...important flag to check if vector is allocated !
  logical ::  is_alloc = .false.
! may be real or integer 
  logical ::  is_real = .true.
  logical ::  is_logi = .false.
! short text attached to the vector
  type(string), dimension(1) :: vect_text
end type Vect

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_isalloc_vect(self)
!-----------------------------------------------------------------------
!
! set the is_alloc attribute to .TRUE.
!
!-----------------------------------------------------------------------
type(Vect), intent(inout) :: self
call QENTER('set_isalloc_vect')
 self%is_alloc = .true.
call QEXIT('set_isalloc_vect')
end subroutine set_isalloc_vect

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine unset_isalloc_vect(self)
!-----------------------------------------------------------------------
!
!      set the is_alloc attribute to .FALSE.
!
!-----------------------------------------------------------------------
type(Vect), intent(inout) :: self
call QENTER('unset_isalloc_vect')
 self%is_alloc = .false.
call QEXIT('unset_isalloc_vect')
end subroutine unset_isalloc_vect

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_vect_kstart(self) result (kstart_out)
!----------------------------------------------------------------------
!
!      returns kstart value of entering vector data type "self"
!
!----------------------------------------------------------------------
#include "priunit.h"
type(Vect), intent(in) :: self
integer(kind=4) :: kstart_out
call QENTER('get_qmtx_kstart')
  if (.not.self%is_alloc) then
    call QUIT('get_vect_kstart: can not return "kstart" position &
               when memory for vector is NOT allocated !')
  else
    kstart_out = self%kstart
  endif
call QEXIT('get_vect_kstart')
end function get_vect_kstart

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_vect_iprint(iprint)
!=================================================
!            set up print variable
!=================================================
#include "priunit.h"
#include "iratdef.h"
integer (kind=8/IRAT), intent(in) :: iprint
CALL QENTER('set_vect_iprint')
 iprint_qmtx = iprint
CALL QEXIT('set_vect_iprint')
end subroutine set_vect_iprint

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_vect_attrib(self)
!-------------------------------------------------------------------
!
!   print out only attributes of the vector data type
!
!-------------------------------------------------------------------
#include "priunit.h"
type(Vect), intent(in) :: self
call QENTER('print_vect_attrib')
 write(LUPRI,'(/2x,a)')  &
 '---------- Basic attributes of the  "vector" data type -----------'
 write(LUPRI,'(a,l1,$)') 'is_alloc:',self%is_alloc
 write(LUPRI,'(a,l1)') '  is real:',self%is_real
 if (self%is_alloc) then
  write(LUPRI,'(a,i7)') 'size of array:',self%array_size
 endif
 if (is_string(self%vect_text)) then
  write(LUPRI,'(4x,a,$)') '...attached text:'
  call print_string(self%vect_text)
 endif
call QEXIT('print_vect_attrib')
end subroutine print_vect_attrib


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function MakeVector(char_type,isize,text) result(self)
!--------------------------------------------------------------------------
!
!   basic constructor for s simple vector: real,integer or logical
!
!  On entry:
!===============
!  char_type  'REAL,'INTE','LOGI' for real, integer and logical variables
!
!
!
!---------------------------------------------------------------------------
use Memory
#include "iratdef.h"
#include "priunit.h"
character(*), intent(in) :: char_type
integer(kind=8/IRAT),intent(in) :: isize
character(*), optional, intent(in) :: text
type(Vect) :: self
call QENTER('MakeVector')

! add the text if is present
if (present(text)) then
   call make_string(self%vect_text, text)
else
   call set_no_string(self%vect_text)
endif

! allocate memmory
if (char_type(1:4) .EQ. 'REAL') then
  CALL mymemget('REAL',self%kstart,isize)
  self%is_real  = .true.
  self%is_logi  = .false.
else if (char_type(1:4) .EQ. 'INTE') then
  CALL mymemget('INTE',self%kstart,isize)
  self%is_real  = .false.
  self%is_logi  = .false.
else if (char_type(1:4) .EQ. 'LOGI') then
  self%is_real  = .false.
  self%is_logi  = .true.
else
  write(LUPRI,*) '=== Output from MakeVector ==='
  write(LUPRI,*) 'Entering vector TYPE: ',char_type(1:4)
  call QUIT('MakeVector - neither REAL/INTE/LOGI ! Dead branch, check it ! ')
endif
self%array_size = isize
self%is_alloc = .true. ! array is allocated

 Call ZeroVector(self)

call QEXIT('MakeVector')
end function MakeVector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine ZeroVector(self)
!-----------------------------------------------------------------------
!
!         Fill vector with zeroes !
!
! For logical vector set them to .FALSE.
!
!-----------------------------------------------------------------------
use Memory
#include "priunit.h"
type(Vect), intent(inout) :: self
call QENTER('ZeroVector')
if (self%is_alloc) then
  if (self%is_real.and..not.self%is_logi) then
     call DZERO( PW(self%kstart), self%array_size )
  else if (.not.self%is_real.and..not.self%is_logi) then
     call IZERO( PW(self%kstart), self%array_size )
  else if (.not.self%is_real.and.self%is_logi) then
     call LZERO( PW(self%kstart), self%array_size )
  else
     call print_vect_attrib(self)
     call QUIT('ZeroVector: dead branch !')
  endif
endif
call QEXIT('ZeroVector')
end subroutine ZeroVector

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine DeAllocVect(self)
!-----------------------------------------------------------------------
!
!     deallocate the vector
!
!-----------------------------------------------------------------------
use Memory
#include "priunit.h"
type(Vect),intent(inout) :: self
call QENTER('DeAllocVect')
if (self%is_alloc) then
 call mymemrel('Vectors: DeAllocVect',self%kstart,self%kstart)
 self%is_alloc = .false.
 self%array_size = 0
 call cancel_string(self%vect_text) ! cancel the string text
endif
call QEXIT('DeAllocVect')
end subroutine DeAllocVect

 END MODULE Vectors
