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

    type(matrix)  :: S, C_i, D, F
    real(kind_dp) :: freq(4)
    real(kind_dp) :: add_to

    freq(1) = aoresponse_freq(1)
    freq(2) = aoresponse_freq(2)
    freq(3) = aoresponse_freq(3)
    freq(4) = aoresponse_freq(4)

    call get_S(S)

    call get_i_or_s_C(C_i, 'i')

    D   = C_i*dag(C_i)
    D   = D/2.0d0
    C_i = 0

    call get_F(F)

    call nr_electrons(S, D)
    call energy(D, F)
    call mu(S, D, F)

    if (aoresponse_alpha) then
      call alpha(S, D, F, (/freq(1), -freq(1)/))
    end if

    if (aoresponse_beta) then
      call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))
    end if

    if (aoresponse_beta_2np1) then
      call beta_2np1(S, D, F, (/freq(1:2), -sum(freq(1:2))/))
    end if

    if (aoresponse_gamma) then
      call gamma(S, D, F, (/freq(1:3), -sum(freq(1:3))/))
    end if

    if (aoresponse_dc_kerr) then

      add_to = 0.0d0

      freq(1) = -aoresponse_freq_process
      freq(2) =  aoresponse_freq_process
      freq(3) =  0.0d0
      freq(4) =  0.0d0

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(8.0d0, add_to, 1, 3, 3, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 3, 3, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '("@ averaged gamma dc Kerr (w = "4f7.3"):" f26.10)', freq(1), &
                                                                  freq(2), &
                                                                  freq(3), &
                                                                  freq(4), &
                                                                  add_to

    end if

    if (aoresponse_dc_shg) then

      add_to = 0.0d0

      freq(1) = -aoresponse_freq_process*2.0d0
      freq(2) =  aoresponse_freq_process
      freq(3) =  aoresponse_freq_process
      freq(4) =  0.0d0

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 1, 3, 3, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 3, 1, 1, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(4.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(4.0d0, add_to, 3, 3, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '("@ averaged gamma dc SHG (w = "4f7.3"):" f26.10)', freq(1), &
                                                                 freq(2), &
                                                                 freq(3), &
                                                                 freq(4), &
                                                                 add_to

    end if

    if (aoresponse_thg) then

      add_to = 0.0d0

      freq(1) = -aoresponse_freq_process*3.0d0
      freq(2) =  aoresponse_freq_process
      freq(3) =  aoresponse_freq_process
      freq(4) =  aoresponse_freq_process

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(6.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(6.0d0, add_to, 3, 3, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '("@ averaged gamma THG (w = "4f7.3"):" f26.10)', freq(1), &
                                                              freq(2), &
                                                              freq(3), &
                                                              freq(4), &
                                                              add_to

    end if

    if (aoresponse_idri) then

      add_to = 0.0d0

      freq(1) = -aoresponse_freq_process
      freq(2) =  aoresponse_freq_process
      freq(3) = -aoresponse_freq_process
      freq(4) =  aoresponse_freq_process

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(8.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(4.0d0, add_to, 1, 3, 1, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '("@ averaged gamma IDRI (w = "4f7.3"):" f26.10)', freq(1), &
                                                               freq(2), &
                                                               freq(3), &
                                                               freq(4), &
                                                               add_to

    end if

    if (aoresponse_beta_static) then

      add_to = 0.0d0
      freq   = 0.0d0

      call beta_element(2.0d0, add_to, 1, 1, 3, S, D, F, (/freq(1:2), -sum(freq(1:2))/))
      call beta_element(1.0d0, add_to, 3, 3, 3, S, D, F, (/freq(1:2), -sum(freq(1:2))/))

      add_to = 3.0d0*add_to/5.0d0

      print '("@ averaged static beta (w = "3f7.3"):" f26.10)', freq(1), &
                                                                freq(2), &
                                                                freq(3), &
                                                                add_to

    end if

    if (aoresponse_gamma_static) then

      add_to = 0.0d0
      freq   = 0.0d0

      call gamma_element( 8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(12.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element( 3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '("@ averaged static gamma (w = "4f7.3"):" f26.10)', freq(1), &
                                                                 freq(2), &
                                                                 freq(3), &
                                                                 freq(4), &
                                                                 add_to

    end if

    S = 0
    D = 0
    F = 0

  end subroutine

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

  end subroutine

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

  end subroutine

end module
