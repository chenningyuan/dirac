!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
module print_level

 integer, private :: iprint=0
 logical, private :: print_level_read = .false.

contains

subroutine extract_print_level
! extract print level from common block variabla
#include "dat.h"
 iprint = iopt(1)
 print_level_read = .true.! print variable is extracted
end subroutine extract_print_level

function getiprint() result(ipr)
! get private module value
  integer :: ipr

 if (.not.print_level_read) then
   call quit('getiprint: print variable NOT extracted !')
 else
  ipr = iprint
  return
 endif

end function getiprint

end module print_level
