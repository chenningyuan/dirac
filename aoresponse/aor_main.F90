module aor_main

  use kinds_m
  use memory_allocator
  use matrix_defop
  use matrix_stat
  use dirac_interface
  use vib_prop
  use aor_cfg
! use aor_binomial
  use aor_property_definitions
  use aor_property_definitions_zzz

  implicit none

contains

  subroutine aoresponse_driver(work, lwork)

#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"

    type(matrix)  :: S, C_i, D, F
    complex(8) :: freq(4)
    real(kind_dp) :: add_to
    integer       :: i, j
    integer, intent(in)    :: lwork
    real(8), intent(inout) :: work(lwork)
    complex(8),allocatable:: tsr(:) !ajt: scratch for property tensors

!   initialize fortran 77 work memory interface
    call di_set_f77_work(work, lwork)
    
!   for debugging, remove later
!   call print_prp_common_blocks()

    freq(1) = aor_cfg_freq(1)
    freq(2) = aor_cfg_freq(2)
    freq(3) = aor_cfg_freq(3)
    freq(4) = aor_cfg_freq(4)

    call get_S(S)

    call get_C(C_i, i=1.0d0, s=0.0d0, g=1.0d0, u=1.0d0)
    D = C_i*dag(C_i)

    call get_F(F)

!   call test_matrix_defop(S, D, F, C_i)
    C_i = 0

    call nr_electrons(S, D)
    call energy(D, F)
    call mu(S, D, F)

    if (aor_cfg_quadrupole) then
      call aor_quadrupole(S, D, F)
    end if

    if (aor_cfg_roa) then
           print '(//a,f11.8/)','linscf calls roa_pol_gpri_aten_gra w=',aor_cfg_freq(1)
           call get_natoms(j); j=3*j
           allocate(tsr(3*3+3*3+3*6+3*3*j+3*3*j+3*6*j)) !ajt: hack
           call roa_pol_gpri_aten_gra(S,D,F,j,aor_cfg_freq(1)*(1d0,0d0),tsr(1:9),tsr(10:18),tsr(19:36), &
                                      tsr(37:36+9*j),tsr(37+9*j:36+18*j),tsr(37+18*j:36+36*j))
           !ajt VIBCTL seems to want g-prime = <<mu,mag>>(w;-w) (opposite sign) (maybe I'm wrong)
           tsr(10:18) = -tsr(10:18)
           tsr(37+9*j:36+18*j) = -tsr(37+9*j:36+18*j)
           call VIBCTL_ifc(j,aor_cfg_freq(1),real(tsr(1:9)),real(tsr(10:18)),real(tsr(19:36)),     &
                           (/((0d0,0d0),i=1,3*j)/),real(tsr(37:36+9*j)),real(tsr(37+9*j:36+18*j)), &
                           real(tsr(37+18*j:36+36*j)),work,lwork)
           deallocate(tsr)
    end if

    if (aor_cfg_cars) then
           print '(//a,f11.8/)','linscf calls cars_shyp_polgra w=',aor_cfg_freq(1)
           call get_natoms(j); j=3*j
           allocate(tsr(3+3*3+3*3*3*3+3*3*j)) !ajt: hack
           call cars_shyp_polgra(S,D,F,j,aor_cfg_freq(1)*(1d0,0d0),tsr(1:3),tsr(4:12), &
                                 tsr(13:93),tsr(94:93+9*j))
           call print_shypol(aor_cfg_freq(1)*(/-1,1,-1,1/),real(tsr(1:3)), &
                             (/real(tsr(4:12)),(0d0,i=1,9*3)/),(/(0d0,i=1,27*6)/),real(tsr(13:93)),lupri)
           call vib_ana_polari(aor_cfg_freq(1)*(/-1,1,-1,1/),real(tsr(1:3)),j,(/(0d0,i=1,3*j*4)/), &
                               (/real(tsr(94:93+9*j)),(0d0,i=1,9*j),real(tsr(94:93+9*j)), &
                               real(tsr(94:93+9*j)),(0d0,i=1,9*j),real(tsr(94:93+9*j))/), &
                               (/(0d0,i=1,27*j*4)/),lupri)
           deallocate(tsr)
    end if

    if (aor_cfg_dg_energy) then
!     call test_grcont(D%nrow, natoms_ifc(), .false., .true.)
!     call debug_1grad(D, F, natoms_ifc())
!     call debug_2grad(D,    natoms_ifc())
      call dg_energy(S, D, F, 3*nr_atoms)
    end if

    if (aor_cfg_alpha) then
      call alpha(S, D, F, (/freq(1), -freq(1)/))
    end if

    if (aor_cfg_magnetizability) then
      call aor_magnetizability(S, D, F, freq(1))
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
            real(freq(1)),                    &
            real(freq(2)),                    &
            real(freq(3)),                    &
            real(freq(4)),                    &
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
            real(freq(1)),                    &
            real(freq(2)),                    &
            real(freq(3)),                    &
            real(freq(4)),                    &
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
            real(freq(1)),                    &
            real(freq(2)),                    &
            real(freq(3)),                    &
            real(freq(4)),                    &
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
            '@ averaged gamma dc IDRI',       &
            real(freq(1)),                    &
            real(freq(2)),                    &
            real(freq(3)),                    &
            real(freq(4)),                    &
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
            real(freq(1)),                          &
            real(freq(2)),                          &
            real(freq(3)),                          &
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
            real(freq(1)),                          &
            real(freq(2)),                          &
            real(freq(3)),                          &
            real(freq(4)),                          &
            add_to

    end if

    S = 0
    D = 0
    F = 0

!radovan: debug this later
!   call print_matrix_stat()

  end subroutine

end module
