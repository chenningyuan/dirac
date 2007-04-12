 module mystring
!--------------------------------------------------------------------
!
!       Module for handling simple dynamical strings
!
!--------------------------------------------------------------------
 type string
  private
   character*1, dimension (:), pointer :: strp  => NULL()
   logical :: is_alloc = .false. ! the default does not work for older F90 compilars !!
   integer (kind=4) :: str_length = 0
 end type string

 contains 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine make_string(self, text)
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
#include "priunit.h"
 type(string), dimension(:), intent(inout) :: self
 character(*), intent(in) :: text
 integer(kind=4) :: i, alloc_stat
 call QENTER('make_string')

 if (.not.self(1)%is_alloc) then
   allocate( self(1)%strp(len(text)), stat = alloc_stat )
   ! ... store the actual size
   self(1)%str_length = size(self(1)%strp)
   if (self(1)%str_length.ne.len(text)) then
     write(LUPRI,*) 'new string lenghts:',len(text),self(1)%str_length
     call QUIT('make_string: error in string length !')
   endif
   if ( alloc_stat > 0 ) then
     write(LUPRI,*) 'alloc_stat erro !'
     call QUIT('mystring/make_string: erro allocation !')
     stop 
   endif

  ! copy the text ....
  do i = 1, len(text)
    self(1)%strp(i) = text(i:i)
  enddo
  self(1)%is_alloc=.true.
 else
  write(LUPRI,*) 'mystring/make_string: this string is already allocated, see:',self(1)%strp
 endif

 call QEXIT('make_string')
 end subroutine make_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine cancel_string(self)
!----------------------------------------------------------------
!
!   cancel the string
!
!----------------------------------------------------------------
#include "priunit.h"
 type(string), dimension(:), intent(inout) :: self
 !type(string), intent(inout) :: self
 integer :: dealloc_stat
 call QENTER('mystring:cancel_string')

 if ( self(1)%is_alloc) then
     deallocate(self(1)%strp, stat = dealloc_stat)
     if (dealloc_stat > 0) then
       write(LUPRI,*) 'dellocation error !'
       call QUIT('mystring/cancel_string: dellocation error !')
     endif
     nullify( self(1)%strp )
     self(1)%is_alloc=.false.
     self(1)%str_length = 0
 endif
 call QEXIT('mystring:cancel_string')
 end subroutine cancel_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_string(self)
!----------------------------------------------------------------------------
!
!  print out the string
!
!----------------------------------------------------------------------------
#include "priunit.h"
type(string), dimension(:), intent(in) :: self
 call QENTER('print_string')
  if (self(1)%is_alloc) then
    write(LUPRI,*) (self(1)%strp)
  endif
 call QEXIT('print_string')
end subroutine print_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 function is_string(self) result(is)
!---------------------------------------------------------------------------
!
!  Returns true/false depending on whether was string allocted/deallocated
!
!----------------------------------------------------------------------------
   type(string), dimension(:), intent(in) :: self
   logical (kind=1) :: is
   call QENTER('is_string')
     is = self(1)%is_alloc
   call QEXIT('is_string')
 end function is_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_string(self) result(str)
!-------------------------------------------------------------------
!
!         returns the string
!
!-------------------------------------------------------------------
  type(string), dimension(:), intent(in) :: self
  character(*) :: str
  integer(kind=4) :: i
  call QENTER('get_string')
  if (self(1)%is_alloc) then
   do i = 1, self(1)%str_length
    str(i:i) = self(1)%strp(i)
   enddo
  else
    str=''
  endif
  call QEXIT('get_string')
end function get_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine set_no_string(self)
!------------------------------------------------------------------------
!
!     sets logical flag indicated allocated string to .FALSE.
!
!------------------------------------------------------------------------
   type(string), dimension(:), intent(inout) :: self
   call QENTER('mystring/set_no_string')
     self(1)%is_alloc = .FALSE.
     self(1)%str_length = 0
   call QEXIT('mystring/set_no_string')
 end subroutine set_no_string

 end module mystring

