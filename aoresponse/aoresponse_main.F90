module aoresponse_main_m

  use kinds_m
  use memory_allocator
  use matrix_mod
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

    logical :: test_nr_electrons = .true.
    logical :: test_energy       = .true.
    logical :: test_mu           = .true.
    logical :: test_alpha        = .true.
    logical :: test_beta         = .false.
    logical :: test_gamma        = .false.

    call get_S(S)

    call get_C(C)
    call get_D_mo(D_mo)

    D = C*D_mo*dag(C)
    D = D/2.0d0

    D_mo = 0
    C    = 0

    call get_F(F)

    if (test_nr_electrons) then
      call nr_electrons(S, D)
    end if

    if (test_energy) then
      call energy(D, F)
    end if

    if (test_mu) then
      call mu(S, D, F)
    end if

    call alloc(freq, nr_freq)
    freq = 0.0d0

    if (test_alpha) then
      freq(1) =  0.00d0
      call alpha(S, D, F, (/freq(1), -freq(1)/))

      freq(1) = -0.07d0
      call alpha(S, D, F, (/freq(1), -freq(1)/))
    end if

    if (test_beta) then
      freq(1) = -0.07d0
      freq(2) =  0.07d0
      call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))

      freq(1) = -0.14d0
      freq(2) =  0.07d0
      call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))

      freq(1) = -0.21d0
      freq(2) =  0.14d0
      call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))
    end if

    if (test_gamma) then
      freq(1) = -0.14d0
      freq(2) =  0.07d0
      freq(3) =  0.07d0
      call gamma(S, D, F, (/freq(1:3), -sum(freq(1:3))/))
    end if

    S = 0
    D = 0
    F = 0

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
