Module Indexed_Quantities

use mystring

use control_print
!use mystring

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
!               deallocate indexed term
!
!----------------------------------------------------------------------------
use mystring
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
!integer, intent(in) :: indxs_array(:)
call QEnter('cancel_it')

if (.not.associated(this)) then
  write(LUPRI,'(2x,a)') 'cancel_it: entering pointer NOT associated !  '
  call Quit('cancel_it: entering pointer NOT associated !')
endif

!if (is_string(this%it_name)) then
  ! ... 
  call cancel_string(this%it_name)
!endif


call QExit('cancel_it')
end subroutine cancel_it

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine fill_indexes(this,indxs_array)
!----------------------------------------------------------------------------
!
! fills indexes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this
integer, intent(in) :: indxs_array(:)
call QEnter('fill_indexes')

if (.not.associated(this)) then
 call Quit('fill_indexes: pointer NOT associated !')
endif

if (size(indxs_array) /= this%nindx) then
 call Quit('fill_indexes: no indexes !')
endif

do i = 1, this%nindx
   this%indxs(i) = indxs_array(i)
enddo

if (get_iprint() > 70) then
  write(LUPRI,'(2x,a,i3)') 'indexes numbers: ',(this%indxs(i),i=1,this%nindx )
endif

call QExit('fill_indexes')
end subroutine fill_indexes

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
subroutine print_it(this)
!----------------------------------------------------------------------------
!
! a   prints allocated scalar indexed quantity - all attributes
!
!----------------------------------------------------------------------------
#include "implicit.h"
#include "priunit.h"
type(indexed_term), pointer :: this

if (.not.associated(this)) then
 call Quit('fill_indexes: pointer NOT associated !')
else
 write(LUPRI,'(2x,a,i3,a,i3)')  'indexed_quantity ID=',this%id,  &
 ' number of indexes',this%nindx
 
endif


end subroutine print_it

End Module Indexed_Quantities

