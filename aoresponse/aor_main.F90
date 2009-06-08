module aor_main

  use kinds_m
  use memory_allocator
  use matrix_defop
  use matrix_stat
  use dirac_interface
  use aor_cfg
  use aor_binomial
  use aor_property_definitions
  use aor_property_definitions_zzz

  implicit none

contains

  subroutine aoresponse_driver()

#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(matrix)  :: S, C_i, D, F
    real(kind_dp) :: freq(4)
    real(kind_dp) :: add_to
    integer       :: natoms_ifc

    call init_binomial()

    freq(1) = aor_cfg_freq(1)
    freq(2) = aor_cfg_freq(2)
    freq(3) = aor_cfg_freq(3)
    freq(4) = aor_cfg_freq(4)

    call get_S(S)

    call get_i_or_s_C(C_i, 'i')

    D   = C_i*dag(C_i)
    D   = 0.5d0*D

    call get_F(F)

!   call test_matrix_defop(S, D, F, C_i)
    C_i = 0

    call nr_electrons(S, D)
    call energy(D, F)
    call mu(S, D, F)

    if (aor_cfg_test_geo) then
!     call test_grcont(D%nrow, natoms_ifc(), .false., .true.)
      call debug_1grad(D, F, natoms_ifc())
      call debug_2grad(D,    natoms_ifc())
      call MolGradi(S,D,F,3*natoms_ifc())
    end if

    if (aor_cfg_alpha) then
      call alpha(S, D, F, (/freq(1), -freq(1)/))
    end if

    if (aor_cfg_beta) then
      call beta(S, D, F, (/freq(1:2), -sum(freq(1:2))/))
    end if

    if (aor_cfg_beta_2np1) then
      call beta_2np1(S, D, F, (/freq(1:2), -sum(freq(1:2))/))
    end if

    if (aor_cfg_gamma) then
      call gamma(S, D, F, (/freq(1:3), -sum(freq(1:3))/))
    end if

    if (aor_cfg_gamma_lin_dc_kerr) then

      add_to = 0.0d0

      freq(1) = -aor_cfg_freq_process
      freq(2) =  aor_cfg_freq_process
      freq(3) =  0.0d0
      freq(4) =  0.0d0

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(8.0d0, add_to, 1, 3, 3, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 3, 3, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma dc Kerr',       &
            freq(1),                          &
            freq(2),                          &
            freq(3),                          &
            freq(4),                          &
            add_to

    end if

    if (aor_cfg_gamma_lin_dc_shg) then

      add_to = 0.0d0

      freq(1) = -aor_cfg_freq_process*2.0d0
      freq(2) =  aor_cfg_freq_process
      freq(3) =  aor_cfg_freq_process
      freq(4) =  0.0d0

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 1, 3, 3, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(2.0d0, add_to, 3, 1, 1, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(4.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(4.0d0, add_to, 3, 3, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma dc SHG',        &
            freq(1),                          &
            freq(2),                          &
            freq(3),                          &
            freq(4),                          &
            add_to

    end if

    if (aor_cfg_gamma_lin_thg) then

      add_to = 0.0d0

      freq(1) = -aor_cfg_freq_process*3.0d0
      freq(2) =  aor_cfg_freq_process
      freq(3) =  aor_cfg_freq_process
      freq(4) =  aor_cfg_freq_process

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(6.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(6.0d0, add_to, 3, 3, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma dc THG',        &
            freq(1),                          &
            freq(2),                          &
            freq(3),                          &
            freq(4),                          &
            add_to

    end if

    if (aor_cfg_zzz) call aor_zzz(aor_cfg_zzz_order, S, D, F)

    if (aor_cfg_gamma_lin_idri) then

      add_to = 0.0d0

      freq(1) = -aor_cfg_freq_process
      freq(2) =  aor_cfg_freq_process
      freq(3) = -aor_cfg_freq_process
      freq(4) =  aor_cfg_freq_process

      call gamma_element(8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(8.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(4.0d0, add_to, 1, 3, 1, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged gamma IDRI',          &
            freq(1),                          &
            freq(2),                          &
            freq(3),                          &
            freq(4),                          &
            add_to

    end if



    if (aor_cfg_beta_lin_static) then

      add_to = 0.0d0
      freq   = 0.0d0

      call beta_element(2.0d0, add_to, 1, 1, 3, S, D, F, (/freq(1:2), -sum(freq(1:2))/))
      call beta_element(1.0d0, add_to, 3, 3, 3, S, D, F, (/freq(1:2), -sum(freq(1:2))/))

      add_to = 3.0d0*add_to/5.0d0

      print '(a30 " (w = "3f9.4"):" f26.10)', &
            '@ averaged static beta',         &
            freq(1),                          &
            freq(2),                          &
            freq(3),                          &
            add_to

    end if

    if (aor_cfg_gamma_lin_static) then

      add_to = 0.0d0
      freq   = 0.0d0

      call gamma_element( 8.0d0, add_to, 1, 1, 1, 1, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element(12.0d0, add_to, 1, 1, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))
      call gamma_element( 3.0d0, add_to, 3, 3, 3, 3, S, D, F, (/freq(1:3), -sum(freq(1:3))/))

      add_to = add_to/15.0d0

      print '(a30 " (w = "4f9.4"):" f26.10)', &
            '@ averaged static gamma',        &
            freq(1),                          &
            freq(2),                          &
            freq(3),                          &
            freq(4),                          &
            add_to

    end if

    S = 0
    D = 0
    F = 0

    call print_matrix_stat()

  end subroutine

end module
