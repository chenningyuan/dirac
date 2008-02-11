Module Indexed_Quantities

use mystring
use control_print

type indexed_term
! assigned unique ID, for instance 1 - Fij, 2 - Tij, 3 - Vpqrs, 4 - Tijab ... 
integer :: id = -1
! ... store the selected name of undexed quantity (like T_ij^ab, F_pq,...)
!character(len=20) :: it_name
type(string), dimension(1) :: it_name
! ... number of indexes (2,4,6...) - must be even number
integer :: nindx 
! ... array for indexes of endexes !
integer, allocatable, dimension(:) :: indxs

end type indexed_term

! ... common name
interface fill_indexes
  module  procedure fill_indexes_2
  module  procedure fill_indexes_4
!  module  procedure fill_indexes_6
!  module  procedure fill_indexes_8
  module  procedure fill_indexes_array
end interface

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
function make_it(id_in, nindx_in, it_name_in ) result(this)
!----------------------------------------------------------------------------
!
!                 Returns ALLOCATED indexed quantity
!
! Entering: id_in - unique ID
!           nindx_in - number of indexes (even number)
!           it_name_in - string to be attached
!
! Like  T_ij, Zijab, Xpq, Y_abi ....
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
integer, intent(in) :: id_in, nindx_in
character(len=*),optional, intent(in) :: it_name_in
type(indexed_term), pointer :: this

call QEnter('make_it')

if (get_iprint() > 70) then
 write(LUPRI,'(2x,a,i3,a,i3)')   &
  'make_it: entering id=',id_in,' entering nindx=',nindx_in
   if ( present(it_name_in) ) then
      write(LUPRI,'(2x,a,a)')   &
         'make_it: also entering name=',it_name_in
    endif
endif

nullify(this)
allocate(this, stat=ierr)

if (ierr /= 0) then
  write(LUPRI,'(2x,a)') 'make_it: error in allocation !'
  call Quit('make_it: error in allocation !')
endif

! ... check even number of indexes
if ( MOD(nindx_in,2) /= 0) then
  write(LUPRI,'(2x,a,i3)')  &
  'number of indexes NOT even number !, nindx_in=',nindx_in
  call Quit('make_it: number of indexes NOT even number !')
endif
this%nindx = nindx_in

this%id = id_in

allocate(this%indxs(nindx),stat=ierr)
if (ierr /= 0) then
  write(LUPRI,'(2x,a)') 'error in allocation !'
  call Quit('make_it: error in allocation !')
endif

this%indxs = -1 

if ( present(it_name_in) ) then
   call cancel_string(this%it_name)
   call make_string(this%it_name,it_name_in)
endif

call QExit('make_it')
end function  make_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine cancel_it(this)
!----------------------------------------------------------------------------
!
!                   completely deallocates indexed term
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
call QEnter('cancel_it')

if (.not.associated(this)) then
  write(LUPRI,'(2x,a)') 'cancel_it: entering pointer NOT associated !  '
  call Quit('cancel_it: entering pointer NOT associated !')
endif

if (is_string(this%it_name)) then
  ! ... deallocate the attached string first
 call cancel_string(this%it_name)
endif

if (allocated(this%indxs)) then
  ! ... deallocate indexes
  deallocate(this%indxs,stat=ierr)
  if (ierr /= 0) then
    call Quit('cancel_it: error in deallocation !')
  endif
endif

! ... now deallocate whole cell indicated with pointer
deallocate(this,stat=ierr)
if (ierr /= 0) then
 call Quit('cancel_it: error in deallocation of "this" !')
endif

call QExit('cancel_it')
end subroutine cancel_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes_array(this,indxs_array)
!----------------------------------------------------------------------------
!
!      fills indexed term indexes numbers from array
!
!  entering "indxs_array" contains indexes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indxs_array(:)
call QEnter('fill_indexes_array')

if (.not.associated(this)) then
 call Quit('fill_indexes: pointer NOT associated !')
endif

if (size(indxs_array) /= this%nindx) then
 call Quit('fill_indexes: no indexes !')
endif

do i = 1, this%nindx
   this%indxs(i) = indxs_array(i)
   if (indxs_array(i).le.0) then
     write(LUPRI,'(2x,a)') 'something wrong !'
   endif
enddo

if (get_iprint() > 70) then
  write(LUPRI,'(2x,a,i3)') 'indexes numbers: ',(this%indxs(i),i=1,this%nindx )
endif

call QExit('fill_indexes_array')
end subroutine fill_indexes_array

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes_2(this,indx1,indx2)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indx1,indx2
call QEnter('fill_indexes_2')

if (.not.associated(this)) then
 call Quit('fill_indexes_2: entering pointer NOT associated !')
endif

if (this%nindx /= 2) then
 call Quit('fill_indexes_2: wrong number of indexes !')
endif

this%indxs(1) = indx1
this%indxs(2) = indx2

call QExit('fill_indexes_2')
end subroutine fill_indexes_2

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes_4(this,indx1,indx2,indx3,indx4)
!----------------------------------------------------------------------------
!
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indx1,indx2,indx3,indx4
call QEnter('fill_indexes_4')

if (.not.associated(this)) then
 call Quit('fill_indexes_4: entering pointer NOT associated !')
endif

if (this%nindx /= 4) then
 call Quit('fill_indexes_4: wrong number of indexes !')
endif

this%indxs(1) = indx1
this%indxs(2) = indx2
this%indxs(3) = indx3
this%indxs(4) = indx4

call QExit('fill_indexes_4')
end subroutine fill_indexes_4

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_it(this)
!----------------------------------------------------------------------------
!
!        prints allocated scalar indexed quantity - all attributes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
call QEnter('print_it')

if (.not.associated(this)) then
 call Quit('fill_indexes: pointer of "type(indexed_term)" NOT associated !')
else
 write(LUPRI,'(3x,a,i1,a,i3,$)')  'indexed_quantity (ID:',this%id,  &
 '); # of indexes',this%nindx
 write(LUPRI,'(1x,a,6i2)') '; indexes:',(this%indxs(i),i=1,this%nindx )

 if (is_string(this%it_name)) then
   write(LUPRI,'(3x,a,$)') 'indexed_quantity text description:'
   call print_string(this%it_name)
 endif
 
endif

call QExit('print_it')
end subroutine print_it

End Module Indexed_Quantities
