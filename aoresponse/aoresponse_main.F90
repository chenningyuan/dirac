module aoresponse_main_m

  use kinds_m
  use memory_allocator
  use dirac_matrix_operations_m
  use dirac_interface_m
  use propint
  use pertdens
  use proptest
  use aoresponse_input_m

  implicit none

! public something

! private

contains

  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(matrix)               :: P
    type(matrix)               :: C
    type(matrix)               :: D, D_mo
    type(matrix)               :: D1
    type(matrix)               :: S
    type(matrix)               :: DSD
    type(matrix)               :: F

    real(kind_dp), allocatable :: freq(:)
    integer(kind_i4)           :: nr_freq = 1

    call alloc(freq, nr_freq)
    freq = 0.0d0

!!!!!!!!!!
    freq(1) = 0.00d0
    freq(1) = 0.02d0
!!!!!!!!!!

    call get_P('ZDIPLEN', P)

    call get_C(C)
    call get_D_mo(D_mo)

    D = C*D_mo*dag(C)

    D1 = 1.0d-100*D

!   call response_solver(freq, P, D1)

    write(*, *) 'P**D   = ', P**D
    write(*, *) 'P**D1  = ', P**D1

    write(*, *) 'compare: -4.828661928957 at 0.00'
    write(*, *) 'compare: -4.833040214626 at 0.02'

    call get_F(F)

    write(*, *) 'F**D   = ', F**D

    call get_S(S)

    write(*, *) 'S**D   = ', S**D

    call dipmom(S, D, F)
    call polari(S, D, F, freq(1))

    P    = 0
    C    = 0
    D    = 0
    D_mo = 0
    D1   = 0
    S    = 0
!   DSD  = 0
    F    = 0

    call dealloc(freq)

  end subroutine aoresponse_driver

end module aoresponse_main_m
