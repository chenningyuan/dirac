module aoresponse_main_m

  use kinds_m
  use memory_allocator
  use dirac_matrix_operations_m
  use dirac_interface_m
  use aoresponse_test_m
  use aoresponse_xc_m

  implicit none

contains

  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(matrix)               :: S, C, D_mo, D, F
    real(kind_dp), allocatable :: freq(:)
    integer                    :: nr_freq = 4

    call alloc(freq, nr_freq)
    freq = 0.0d0

    call get_S(S)

    call get_C(C)

    call mat_init(D, ntbas(0), ntbas(0))
    D%elements = 0.0d0

    call genden(D%elements, C%elements, 1, 0)

    C    = 0

!   call get_D_mo(D_mo)
!   D = C*D_mo*dag(C)
!   D_mo = 0

    D = D/2.0d0
    D%is_density_matrix = .true.

    call get_F(F)

    call nr_electrons(S, D, F)
    call energy(S, D, F)
    call mu(S, D, F)

!   freq(1) = 0.00d0
!   call alpha(S, D, F, freq(1))

!   freq(1) = 0.02d0
!   call alpha(S, D, F, freq(1))

!   freq(1) = 0.07d0
!   call alpha(S, D, F, freq(1))








    freq(1) = -0.07d0
    freq(2) =  0.07d0
    call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))

    freq(1) = -0.14d0
    freq(2) =  0.07d0
    call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))

    freq(1) = -0.21d0
    freq(2) =  0.14d0
    call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))


!   freq(3) =  0.00d0
!   freq    =  0.00d0

!   call gamma(S, D, F, (/freq(1:3), -sum(freq(1:3))/))

    S    = 0
    D    = 0
    F    = 0

    call dealloc(freq)

  end subroutine aoresponse_driver

  subroutine test_matrix_stat()

    type(matrix) :: A, B, C

                            call matrixstats('entering subroutine   ')

    call mat_init(A, 5, 5); call matrixstats('call mat_init(A, 5, 5)')
    A = A                 ; call matrixstats('A = A                 ')
    B = 2.0d0*A           ; call matrixstats('B = 2.0d0*A           ')
    C = A + B             ; call matrixstats('C = A + B             ')
    C = C/2.0d0           ; call matrixstats('C = C/2.0d0           ')

    call test_matrix_stat_sub()
    call test_matrix_stat_sub()
    call test_matrix_stat_sub()

    A = 0                 ; call matrixstats('A = 0                 ')
    B = 0                 ; call matrixstats('B = 0                 ')
    C = 0                 ; call matrixstats('C = 0                 ')

                            call matrixstats('leaving subroutine    ')

  end subroutine test_matrix_stat

  subroutine test_matrix_stat_sub()

    type(matrix) :: A, B, C

                            call matrixstats('entering subroutine   ')

    call mat_init(A, 5, 5); call matrixstats('call mat_init(A, 5, 5)')
    A = A                 ; call matrixstats('A = A                 ')
    B = 2.0d0*A           ; call matrixstats('B = 2.0d0*A           ')
    C = A + B             ; call matrixstats('C = A + B             ')
    C = C/2.0d0           ; call matrixstats('C = C/2.0d0           ')

    A = 0                 ; call matrixstats('A = 0                 ')
    B = 0                 ; call matrixstats('B = 0                 ')
    C = 0                 ; call matrixstats('C = 0                 ')

                            call matrixstats('leaving subroutine    ')

  end subroutine test_matrix_stat_sub

end module aoresponse_main_m
