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

 module mystring
!--------------------------------------------------------------------
!
!     Module for handling simple dynamical string
!
!--------------------------------------------------------------------
 type string
  private
   character*1, dimension (:), pointer :: strp  => NULL()
   logical :: is_alloc = .false. ! this default does not work for older F90 compilars !!
   integer (kind=4) :: str_length = 0
 end type string

 integer, private :: iprint_string = 0  ! control print-out variable

 contains 

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine make_string(self, text)
!-----------------------------------------------------------------------
!
!       create (allocate) the string
!
!-----------------------------------------------------------------------
#include "priunit.h"
 type(string), dimension(:), intent(inout) :: self
 character(*), intent(in) :: text
 integer(kind=4) :: i, alloc_stat

 if (.not.self(1)%is_alloc) then
   allocate( self(1)%strp(len(text)), stat = alloc_stat )
   ! ... store the actual size
   self(1)%str_length = size(self(1)%strp)
   if (self(1)%str_length.ne.len(text)) then
     write(LUPRI,*) 'new string lenghts:',len(text),self(1)%str_length
     call FLSHFO(LUPRI)
     call QUIT('make_string: error in string length !')
   endif
   if ( alloc_stat > 0 ) then
     write(LUPRI,*) 'alloc_stat erro !'
     call FLSHFO(LUPRI)
     call QUIT('mystring/make_string: erro allocation !')
     stop 
   endif

  ! copy the text ....
  do i = 1, len(text)
    self(1)%strp(i) = text(i:i)
  enddo
  self(1)%is_alloc=.true.
 else
   ! fix for older F90 compilers which do not assign default values 
 !  call FLSHFO(LUPRI)
   !write(lupri,*) 'make_string hallooo ! size(self(1)%strp)=',size(self(1)%strp)
  ! write(lupri,*) 'make_string hallooo ! size(self(1)%strp)=',SIZE(self)
 !  call FLSHFO(LUPRI)
 !  if ( size(self(1)%strp) > 0 ) then
     write(LUPRI,*) 'mystring/make_string: this string '// &
     'is already allocated, see: ',self(1)%strp
  ! endif
 endif

 end subroutine make_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine cancel_string(self)
!----------------------------------------------------------------
!
!   cancel (deallocate) the string
!
!----------------------------------------------------------------
#include "priunit.h"
 type(string), dimension(:), intent(inout) :: self
 !type(string), intent(inout) :: self
! must be kind=4 otherwise generates error !
 integer(kind=4) :: dealloc_stat

 if ( self(1)%is_alloc) then
 !    write(LUPRI,*) 'mystring:cancel_string entering string: ',self(1)%strp
 !    write(LUPRI,*) ' length of string is:',self(1)%str_length
     deallocate(self(1)%strp, stat = dealloc_stat)
     if (dealloc_stat /= 0) then
       write(LUPRI,*) 'dellocation error !, dealloc_stat=',dealloc_stat
       write(LUPRI,*) ' length of the delloc. string was:',self(1)%str_length
       call FLSHFO(LUPRI)
       call QUIT('mystring/cancel_string: deallocation error !')
     endif
     nullify( self(1)%strp )
     self(1)%is_alloc=.false.
     self(1)%str_length = 0
 endif
 end subroutine cancel_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_string(self)
!----------------------------------------------------------------------------
!
!       print out the string (LUPRI channel)
!
!----------------------------------------------------------------------------
#include "priunit.h"
type(string), dimension(:), intent(in) :: self
  if (self(1)%is_alloc) then
    write(LUPRI,*) (self(1)%strp)
  endif
  call FLSHFO(LUPRI)
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
     is = self(1)%is_alloc
 end function is_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function Xget_string(self) result(str_out)
!-------------------------------------------------------------------
!
!        should return the string
!
!-------------------------------------------------------------------
  type(string), dimension(:), intent(in) :: self
  character :: str_out  !  fix for gfortran
  !character(*) :: str_out
  integer(kind=4) :: i
  call QENTER('get_string')
  if (self(1)%is_alloc) then
   do i = 1, self(1)%str_length
    str_out(i:i) = self(1)%strp(i)
   enddo
  else
    str_out=' '
  endif
  call QEXIT('get_string')
end function Xget_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function get_string(self) result(str_out)
!-------------------------------------------------------------------
!
!             should return the string
!
!-------------------------------------------------------------------
  type(string), dimension(:), intent(in) :: self
  !character :: str_out  !  fix for gfortran
  !character(len=20) :: str_out ! the rest of length is filled with rubbish
  !character(*), allocatable(:),save :: str_out
  !character*1, dimension (:), pointer, save :: str_out  => NULL()
  !character*1, dimension (:), pointer :: str_out  => NULL()
  character*1, dimension (:), pointer :: str_out  

  integer(kind=4) :: i


  !if ( allocated(str_out(1)) ) then
  !if ( associated(str_out(1)) ) then
  if ( associated(str_out) ) then
   deallocate(str_out); nullify(str_out)
  endif
  !allocate( str_out(self(1)%str_length) )

  if (self(1)%is_alloc) then

   allocate( str_out(self(1)%str_length) )

   do i = 1, self(1)%str_length
    str_out(i:i) = self(1)%strp(i)
   enddo
  else
    str_out=' '
  endif

end function get_string

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine copy_strings(source,dest)
!-----------------------------------------------------------------------
!
!  copy allocated "source" -> "dest"; allocate "dest" first
!
!-----------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
!type(string), dimension(:), intent(inout) :: source, dest
!type(string), dimension(1), intent(inout) :: source, dest
type(string), dimension(1), intent(inout) :: source, dest

if (.not.source(1)%is_alloc) then
   call Quit('copy_strings: "source" is NOT allocated !')
endif

if (associated(dest(1)%strp)) then
  call cancel_string(dest)
endif

allocate( dest(1)%strp(source(1)%str_length), stat =  ierr )
if (ierr /= 0) then
  call Quit('copy_strings: allocation error of "dest(1)%strp" !')
endif

dest(1)%is_alloc   = .true.
dest(1)%str_length = source(1)%str_length

!...copy the text
do i = 1, dest(1)%str_length
 dest(1)%strp(i) = source(1)%strp(i)
enddo

end subroutine copy_strings

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine set_no_string(self)
!------------------------------------------------------------------------
!
!    sets logical flag indicated allocated string to .FALSE.
!
!------------------------------------------------------------------------
type(string), dimension(:), intent(inout) :: self
     self(1)%is_alloc = .FALSE.
     self(1)%str_length = 0
end subroutine set_no_string

end module mystring
