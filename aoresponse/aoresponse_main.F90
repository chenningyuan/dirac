module aoresponse_main_m

  use kinds_m
  use memory_allocator
  use matrix_mod
  use dirac_interface_m
  use aoresponse_input_m
  use aoresponse_test_m
  use aoresponse_xc_m

  implicit none

contains

  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(matrix)  :: S, C, D_mo, D, F
    real(kind_dp) :: freq(4)

    freq(1) = aoresponse_freq(1)
    freq(2) = aoresponse_freq(2)
    freq(3) = aoresponse_freq(3)
    freq(4) = aoresponse_freq(4)

    call get_S(S)

    call get_C(C)
    call get_D_mo(D_mo)

    D = C*D_mo*dag(C)
    D = D/2.0d0

    D_mo = 0
    C    = 0

    call get_F(F)

!   if (test_nr_electrons) then
      call nr_electrons(S, D)
!   end if

!   if (test_energy) then
      call energy(D, F)
!   end if

!   if (test_mu) then
      call mu(S, D, F)
!   end if

!   0.08855851799036d0 = 514.5 nm

    if (aoresponse_alpha) then
      call alpha(S, D, F, (/freq(1), -freq(1)/))
    end if

    if (aoresponse_beta) then
      call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))
    end if

    if (aoresponse_gamma) then
      call gamma(S, D, F, (/freq(1:3), -sum(freq(1:3))/))
!     call gamma_element(1, 2, 1, 2, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
!     call gamma_element(1, 1, 2, 2, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
!     call gamma_element(3, 1, 3, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
!     call gamma_element(1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
!     call gamma_element(1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
!     call gamma_element(3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
    end if

    S = 0
    D = 0
    F = 0

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
