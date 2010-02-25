module pammag_info

  use interface_for_pammag 

  implicit none

  save

  integer, parameter, public :: nr_fock    = 3    ! x/y/z
!from old twolon:
  integer, parameter, public :: ncomp      = 3    ! x/y/z
  integer, parameter, public :: ndmat      = 4
  
!  public f77_functions
  public lbit


!  interface to f77_functions

  interface lbit
    logical function lbit(i, n)
      integer :: i, n
    end function
  end interface


  private

!contains


end module
