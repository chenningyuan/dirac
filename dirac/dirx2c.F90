module dirx2c

! stefan: - this module will contain all necessary functionality 
!           to perform the transformation from a four-component 
!           to a two-component framework 

! use so_far_nothing

  implicit none

  public say_hello_dirx2c

  private

#include "dgroup.h"

contains

  subroutine say_hello_dirx2c()

    write(6,*) ' *** Welcome to the new X2C interface in DIRAC ***'

  end subroutine

end module
