!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

module input_reader

! radovan: - this module should not initialize anything - it is a reader
!          - when we have only modules and no common blocks the initialization
!            problem goes away since variables can be initialized at definition
!          - please feel free to copy/use/improve this inside DIRAC/Dalton
!            without asking

   use keyword
   use character_processing
   use dirac_cfg
   use interface_input

   implicit none

   public read_input
   public move_to_next_star

   private

#include "dgroup.h"

contains

   subroutine read_input()

!    radovan: - this routine should not include any common blocks
!               this should be done in called subroutines
!             - please keep ifndef traps inside subroutines for better readability
!               of this subroutine

!     --------------------------------------------------------------------------
      character(kw_length) :: word
      character(kw_length) :: kw_section
!     --------------------------------------------------------------------------

!     this is to catch keywords that appear before some section starts
      kw_section = '       '

      open(unit_in, file='DIRAC.INP')
      rewind unit_in

      do while (.true.)

         read(unit_in, '(a7)', end=1) word

         if (word == '       ') then
!           blank line
         else
            select case (word(1:1))

               case ('!', '#')
!                 comment

               case ('*')
!                 section
                  kw_section = uppercase(word)

               case default
!                 keyword
                  select case (kw_section)

!                    case ('**METHO', '**WAVE ')
!                       call read_input_method(word, kw_section)

                     case ('*DFT   ')
                        call read_input_dft(word, kw_section)

                     case ('**GRID ')
                        call read_input_grid(word, kw_section)

!                    *VISUAL is for backward compatibility
!                    *VISUAL used to be under **ANALYZE
                     case ('**VISUA', '*VISUAL')
                        call read_input_visual(word, kw_section)

                     case ('*OPENRS')
                        call read_input_openrsp(word, kw_section)

                     case ('*X2C')
                        call read_input_x2c(word, kw_section)

                     case ('**RELCC')
                        call read_input_relcc(word, kw_section)

                     case ('*CCENER')
                        call read_input_ccener(word, kw_section)

                     case ('*CCFOPR')
                        call read_input_ccfopr(word, kw_section)

                     case ('*CCFSPC')
                        call read_input_ccfspc(word, kw_section)

                     case ('*CCIH  ')
                        call read_input_ccih(word, kw_section)

                     case ('*END OF')
                        go to 1

!                    case ('       ')
!                       call quit('keyword '//word//' without section')

                     case default
!                       activate as soon as everything is merged to here
!                       call quit('section '//kw_section//' not recognized')

                  end select

            end select
         end if

      end do

1     close(unit_in)

   end subroutine

   subroutine move_to_next_star(word_io)

!     --------------------------------------------------------------------------
      character(kw_length), intent(inout) :: word_io
!     --------------------------------------------------------------------------
      character(kw_length)                :: word
!     --------------------------------------------------------------------------

      do while (.true.)
         read(unit_in, '(a7)') word
         if (word(1:1) == '*') then
            word_io = word
            return
         end if
      end do

   end subroutine

   subroutine read_input_grid(word, kw_section)

      use num_grid_cfg

!     ---------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     ---------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.RADINT')) then
         call kw_read(word, num_grid_cfg_radint)
      end if

      if (kw_matches(word, '.ANGINT')) then
         call kw_read(word, num_grid_cfg_angint)
      end if

      if (kw_matches(word, '.NOPRUN')) then
!        turn off pruning of angular grid
         num_grid_cfg_no_pruning = .true.
      end if

      if (kw_matches(word, '.ANGMIN')) then
         call kw_read(word, num_grid_cfg_angmin)
      end if

      if (kw_matches(word, '.ATSIZE')) then
!        estimate relative atomic sizes for use in the Becke
!        partitioning scheme from atomic contributions
         num_grid_cfg_estimate_radii = .true.
      end if

      if (kw_matches(word, '.IMPORT')) then
         num_grid_cfg_import_grid = .true.
         call kw_read(word, num_grid_cfg_gridfile)
      end if

!radovan: not documented on wiki
      if (kw_matches(word, '.ZIPGRI')) then
!        is true by default, keep for backward compatibimity
         num_grid_cfg_zipgrid = .true.
      end if

      if (kw_matches(word, '.NOZIP ')) then
         num_grid_cfg_zipgrid = .false.
      end if

      if (kw_matches(word, '.4CGRID')) then
!        include small component in the grid generation
!        also for 1- and 2-c dft calculations
         num_grid_cfg_force_4c_grid = .true.
      end if

      if (kw_matches(word, '.DEBUG ')) then
         num_grid_cfg_radint = 1.0d-3
         num_grid_cfg_angint = 10
      end if

      if (kw_matches(word, '.COARSE')) then
         num_grid_cfg_radint = 1.0d-11
         num_grid_cfg_angint = 35
      end if

      if (kw_matches(word, '.ULTRAF')) then
         num_grid_cfg_radint = 2.0d-15
         num_grid_cfg_angint = 64
      end if

      if (kw_matches(word, '.INTCHK')) then
         call kw_read(word, num_grid_cfg_integration_check_level)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_visual(word, kw_section)

      use visual_cfg
      use character_processing

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i
      real(8)                          :: s, v_12(3), v_13(3)
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      dirac_cfg_visual = .true.

      if (kw_matches(word, '.LIST  ')) then
         visual_cfg_list = .true.
         call kw_read(word, visual_cfg_nr_points_in_list)
         allocate(visual_cfg_xyz_list(visual_cfg_nr_points_in_list, 3))
         do i = 1, visual_cfg_nr_points_in_list
            call kw_read(word, visual_cfg_xyz_list(i, 1), &
                               visual_cfg_xyz_list(i, 2), &
                               visual_cfg_xyz_list(i, 3))
         end do
      end if

      if (kw_matches(word, '.LINE  ')) then
         visual_cfg_line = .true.
         call kw_read(word, visual_cfg_line_from(1), &
                            visual_cfg_line_from(2), &
                            visual_cfg_line_from(3))
         call kw_read(word, visual_cfg_line_to(1), &
                            visual_cfg_line_to(2), &
                            visual_cfg_line_to(3))
         call kw_read(word, visual_cfg_line_nr_steps)
      end if

      if (kw_matches(word, '.2D    ')) then
         visual_cfg_2d = .true.
         read(unit_in, *) visual_cfg_2d_p_origin(1:3)
         read(unit_in, *) visual_cfg_2d_p_right(1:3)
         read(unit_in, *) visual_cfg_2d_nr_right
         read(unit_in, *) visual_cfg_2d_p_top(1:3)
         read(unit_in, *) visual_cfg_2d_nr_top
!        check that angle between 1-2 and 1-3 is 90 degs
         v_12 = visual_cfg_2d_p_right - visual_cfg_2d_p_origin
         v_13 = visual_cfg_2d_p_top   - visual_cfg_2d_p_origin
         s = v_12(1)*v_13(1) + v_12(2)*v_13(2) + v_12(3)*v_13(3)
         if (s > tiny(0.0d0)) then
            call quit('.2D: lines connecting points 1-2 and 1-3 do not enclose 90 degrees')
         end if
      end if

      if (kw_matches(word, '.2D_INT')) then
         visual_cfg_2d_integration = .true.
         read(unit_in, *) visual_cfg_2d_integration_p_origin(1:3)
         read(unit_in, *) visual_cfg_2d_integration_p_right(1:3)
         read(unit_in, *) visual_cfg_2d_integration_nr_right
         read(unit_in, *) visual_cfg_2d_integration_p_top(1:3)
         read(unit_in, *) visual_cfg_2d_integration_nr_top
         read(unit_in, *) visual_cfg_2d_integration_order
!        check that angle between 1-2 and 1-3 is 90 degs
         v_12 = visual_cfg_2d_integration_p_right - visual_cfg_2d_integration_p_origin
         v_13 = visual_cfg_2d_integration_p_top   - visual_cfg_2d_integration_p_origin
         s = v_12(1)*v_13(1) + v_12(2)*v_13(2) + v_12(3)*v_13(3)
         if (s > tiny(0.0d0)) then
            call quit('.2D_INT: lines connecting points 1-2 and 1-3 do not enclose 90 degrees')
         end if
      end if

      if (kw_matches(word, '.3D    ')) then
         visual_cfg_3d = .true.
         call kw_read(word, visual_cfg_ncube(1), &
                            visual_cfg_ncube(2), &
                            visual_cfg_ncube(3))
      end if

      if (kw_matches(word, '.3D_ADD')) then
         call kw_read(word, visual_cfg_add_3d)
      end if

      if (kw_matches(word, '.3D_INT')) then
         visual_cfg_3d_integration = .true.
      end if

      if (kw_matches(word, '.SCALE ')) then
         call kw_read(word, visual_cfg_scale)
      end if

      if (kw_matches(word, '.DENSIT')) then
         call read_input_visual_helper(iq_density)
      end if

      if (kw_matches(word, '.ELF   ')) then
         call read_input_visual_helper(iq_elf)
      end if

      if (kw_matches(word, '.GAMMA5')) then
         call read_input_visual_helper(iq_gamma5)
      end if

      if (kw_matches(word, '.JDIA  ')) then
         call read_input_visual_helper(iq_jdia)
      end if

      if (kw_matches(word, '.J     ')) then
         call read_input_visual_helper(iq_j)
      end if

      if (kw_matches(word, '.DIVJ  ')) then
         call read_input_visual_helper(iq_divj)
      end if

      if (kw_matches(word, '.ROTJ  ')) then
         call read_input_visual_helper(iq_rotj)
      end if

      if (kw_matches(word, '.S     ')) then
         call read_input_visual_helper(iq_s)
      end if

      if (kw_matches(word, '.DIVS  ')) then
         call read_input_visual_helper(iq_divs)
      end if

      if (kw_matches(word, '.ROTS  ')) then
         call read_input_visual_helper(iq_rots)
      end if

      if (kw_matches(word, '.EDIPX ')) then
         call read_input_visual_helper(iq_edipx)
      end if
      if (kw_matches(word, '.EDIPY ')) then
         call read_input_visual_helper(iq_edipy)
      end if
      if (kw_matches(word, '.EDIPZ ')) then
         call read_input_visual_helper(iq_edipz)
      end if

      if (kw_matches(word, '.BDIPX ')) then
         call read_input_visual_helper(iq_bdipx)
      end if
      if (kw_matches(word, '.BDIPY ')) then
         call read_input_visual_helper(iq_bdipy)
      end if
      if (kw_matches(word, '.BDIPZ ')) then
         call read_input_visual_helper(iq_bdipz)
      end if

      if (kw_matches(word, '.NDIPX ')) then
         call read_input_visual_helper2(iq_ndipx)
      end if
      if (kw_matches(word, '.NDIPY ')) then
         call read_input_visual_helper2(iq_ndipy)
      end if
      if (kw_matches(word, '.NDIPZ ')) then
         call read_input_visual_helper2(iq_ndipz)
      end if

      if (kw_matches(word, '.BDIPXD')) then
         call read_input_visual_helper(iq_bdipxdia)
      end if
      if (kw_matches(word, '.BDIPYD')) then
         call read_input_visual_helper(iq_bdipydia)
      end if
      if (kw_matches(word, '.BDIPZD')) then
         call read_input_visual_helper(iq_bdipzdia)
      end if

      if (kw_matches(word, '.NDIPXD')) then
         call read_input_visual_helper2(iq_ndipxdia)
      end if
      if (kw_matches(word, '.NDIPYD')) then
         call read_input_visual_helper2(iq_ndipydia)
      end if
      if (kw_matches(word, '.NDIPZD')) then
         call read_input_visual_helper2(iq_ndipzdia)
      end if
      if (kw_matches(word, '.ESP   ')) then
         call read_input_visual_helper(iq_esp)
      end if
      if (kw_matches(word, '.ESPRHO')) then
         call read_input_visual_helper(iq_esprho)
      end if

      if (kw_matches(word, '.GAUGE ')) then
         call kw_read(word, visual_cfg_gauge_origin(1), &
                            visual_cfg_gauge_origin(2), &
                            visual_cfg_gauge_origin(3))
      end if

      if (kw_matches(word, '.RADIAL')) then
         visual_cfg_radial = .true.
      end if

      if (kw_matches(word, '.SMALLA')) then
         visual_cfg_force_small_ao = .true.
      end if

      if (kw_matches(word, '.ORBSTR')) then
         visual_cfg_use_orbital_string = .true.
         call kw_read(word, visual_cfg_orbital_string(1))
         if (nfsym > 1) then
            call kw_read(word, visual_cfg_orbital_string(2))
         end if
      end if

#ifdef MOD_UNRELEASED
      if (kw_matches(word, '.LONDON')) then
         visual_cfg_london = .true.
         read(unit_in, *) visual_cfg_london_component
      end if

      if (kw_matches(word, '.NONE  ')) then
         visual_cfg_london_none      = .true.
      end if

      if (kw_matches(word, '.NODIRE')) then
         visual_cfg_london_skip_direct = .true.
      end if

      if (kw_matches(word, '.NOREOR')) then
         visual_cfg_london_skip_ro = .true.
      end if

      if (kw_matches(word, '.NOKAPP')) then
         visual_cfg_london_skip_kappa = .true.
      end if
#endif

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_visual_helper(j)

      use visual_cfg

!     --------------------------------------------------------------------------
      integer, intent(in) :: j
!     --------------------------------------------------------------------------
      integer             :: i
      character(80)       :: line
!     --------------------------------------------------------------------------

      visual_cfg_nr_dmat = visual_cfg_nr_dmat + 1
      i = visual_cfg_nr_dmat
      visual_cfg_property(i) = j
      read(unit_in, *) line
      backspace(unit_in)
      if (word_contains(line, 'DFCOEF')) then
         read(unit_in, *) visual_cfg_dmat_file(i)
      else
         read(unit_in, *) visual_cfg_dmat_file(i), &
                          visual_cfg_dmat_file_record(i)
      end if

   end subroutine

   subroutine read_input_visual_helper2(j)

      use visual_cfg

!     --------------------------------------------------------------------------
      integer, intent(in) :: j
!     --------------------------------------------------------------------------
      integer             :: i
!     --------------------------------------------------------------------------

      visual_cfg_nr_dmat = visual_cfg_nr_dmat + 1
      i = visual_cfg_nr_dmat
      visual_cfg_property(i) = j
      read(unit_in, *) visual_cfg_ref_nucleus(i), &
                       visual_cfg_dmat_file(i),   &
                       visual_cfg_dmat_file_record(i)

   end subroutine

   subroutine read_input_openrsp(word, kw_section)

#ifdef MOD_OPENRSP
      use openrsp_cfg
#endif

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i
!     --------------------------------------------------------------------------

#ifndef MOD_OPENRSP
      call quit(kw_section//' not available in this version')
#else

      call reset_available_kw_list()

      if (kw_matches(word, '.THRESH')) then
         call kw_read(word, openrsp_cfg_threshold_response, openrsp_cfg_threshold_rhs)
      end if

      if (kw_matches(word, '.FREQ  ')) then
         call kw_read(word, openrsp_cfg_nr_freq)
         allocate(openrsp_cfg_freq(openrsp_cfg_nr_freq))
         do i = 1, openrsp_cfg_nr_freq
            call kw_read(word, openrsp_cfg_freq(i))
         end do
      end if

      if (kw_matches(word, '.FREQ_P')) then
!        give frequency for a process
         call kw_read(word, openrsp_cfg_freq_process)
      end if

      if (kw_matches(word, '.SKIPEE')) then
         openrsp_cfg_skip_pp = .true.
      end if

      if (kw_matches(word, '.SKIPEP')) then
         openrsp_cfg_skip_pn = .true.
      end if

      if (kw_matches(word, '.SKIP1E')) then
         openrsp_cfg_skip_1el = .true.
      end if

      if (kw_matches(word, '.ALPHA ')) then
!        frequency-dependent \alpha
         openrsp_cfg_alpha = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.0ALPHA')) then
!        static \alpha
         openrsp_cfg_alpha_lin_static = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.BETA  ')) then
!        frequency-dependent \beta by the n+1 rule
         openrsp_cfg_beta = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.BETA2N')) then
!        frequency-dependent \beta by the 2n+1 rule
         openrsp_cfg_beta_2np1 = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.0BETA ')) then
!        static \beta
         openrsp_cfg_beta_lin_static = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.GAMMA ')) then
!        frequency-dependent \gamma
         openrsp_cfg_gamma = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.0GAMMA')) then
!        static \gamma
         openrsp_cfg_gamma_lin_static = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.DC-SHG')) then
         openrsp_cfg_gamma_lin_dc_shg = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.THG   ')) then
         openrsp_cfg_gamma_lin_thg = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.IDRI  ')) then
         openrsp_cfg_gamma_lin_idri = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.DC-KER')) then
         openrsp_cfg_gamma_lin_dc_kerr = .true.
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.MAGNET')) then
         openrsp_cfg_magnetizability = .true.
         need_1el_o = .true.
         need_1el_b = .true.
      end if

      if (kw_matches(word, '.QUADRU')) then
         openrsp_cfg_quadrupole = .true.
         need_1el_q = .true.
      end if

      if (kw_matches(word, '.ROA   ')) then
         openrsp_cfg_roa = .true.
         need_1el_f = .true.
         need_1el_q = .true.
         need_1el_o = .true.
         need_1el_g = .true.
      end if

      if (kw_matches(word, '.CARS  ')) then
         openrsp_cfg_cars = .true.
         need_1el_f = .true.
         need_1el_q = .true.
         need_1el_o = .true.
         need_1el_g = .true.
      end if

      if (kw_matches(word, '.PNC   ')) then
         openrsp_cfg_pnc = .true.
         need_1el_g = .true.
      end if

      if (kw_matches(word, '.PVBETA')) then
         openrsp_cfg_pv_beta = .true.
         need_1el_f = .true.
         need_1el_g = .true.
      end if

      if (kw_matches(word, '.PVGAMM')) then
         openrsp_cfg_pv_gamma = .true.
         need_1el_f = .true.
         need_1el_g = .true.
      end if

      if (kw_matches(word, '.NOLAO ')) then
         openrsp_cfg_nolao = .true.
      end if

      if (kw_matches(word, '.EFGB  ')) then
         openrsp_cfg_efgb = .true.
         need_1el_f = .true.
         need_1el_q = .true.
         need_1el_v = .true.
         need_1el_o = .true.
         need_1el_b = .true.
      end if

      if (kw_matches(word, '.JONES ')) then
         openrsp_cfg_jones = .true.
         need_1el_f = .true.
         need_1el_q = .true.
         need_1el_o = .true.
         need_1el_b = .true.
      end if

      if (kw_matches(word, '.CME   ')) then
         openrsp_cfg_cme = .true.
         need_1el_f = .true.
         need_1el_q = .true.
         need_1el_o = .true.
         need_1el_b = .true.
      end if

      if (kw_matches(word, '.MOLGRD')) then
         openrsp_cfg_dg_energy = .true.
         need_1el_f = .true.
         need_1el_g = .true.
      end if

      if (kw_matches(word, '.INTEGR')) then
         openrsp_cfg_integrate = .true.
      end if

      if (kw_matches(word, '.VISUAL')) then
         openrsp_cfg_visual = .true.
      end if

      if (kw_matches(word, '.ZZZ   ')) then
!        calculate static derivatives E_{zzz...}
!        of given order
         openrsp_cfg_zzz = .true.
         call kw_read(word, openrsp_cfg_zzz_order)
         need_1el_f = .true.
      end if

      if (kw_matches(word, '.DELTA ')) then
         openrsp_cfg_test_delta = .true.
         need_1el_f = .true.
      end if

      call check_whether_kw_found(word, kw_section)
#endif /* ifndef MOD_OPENRSP */

   end subroutine

   subroutine read_input_relcc(word, kw_section)

      use relcc_cfg

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.NELEC ')) then
         nelec_input = .true.
         call kw_read(word, relcc_nelec, 2*nfsym)
      end if

      if (kw_matches(word, '.NEL_F1')) then
         call kw_read(word, relcc_nelec_f1, 2*nfsym)
      end if

      if (kw_matches(word, '.NEL_F2')) then
         call kw_read(word, relcc_nelec_f2, 2*nfsym)
      end if

      if (kw_matches(word, '.NFROZ ')) then
         call kw_read(word, relcc_nfroz, 2*nfsym)
      end if

      if (kw_matches(word, '.DEBUG ')) then
         relcc_debug = .true.
      end if

      if (kw_matches(word, '.CARITH')) then
         relcc_carith = .true.
      end if

      if (kw_matches(word, '.TIMING')) then
         relcc_timing = .true.
      end if

      if (kw_matches(word, '.NOSORT')) then
         relcc_do_sort = .false.
      end if

      if (kw_matches(word, '.ENERGY')) then
         relcc_do_energy = .true.
      end if

      if (kw_matches(word, '.GRADIE')) then
         relcc_do_gradient = .true.
      end if

      if (kw_matches(word, '.FOCKSP')) then
         relcc_do_fspc     = .true.
      end if

      if (kw_matches(word, '.PRINT ')) then
         call kw_read(word, relcc_print)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccener(word, kw_section)

      use relcc_cfg

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.DOMP2 ')) then
         relcc_do_mp2 = .true.
      end if

      if (kw_matches(word, '.NOSING')) then
         relcc_do_singles = .false.
      end if

      if (kw_matches(word, '.NODOUB')) then
         relcc_do_doubles = .false.
      end if

      if (kw_matches(word, '.NOSD  ')) then
         relcc_do_ccsd = .false.
      end if

      if (kw_matches(word, '.NOSDT ')) then
         relcc_do_ccsd_t = .false.
      end if

      if (kw_matches(word, '.MAXIT ')) then
         read(unit_in,'(i5)') relcc_ccener_max_iterations
      end if

      if (kw_matches(word, '.MAXDIM')) then
         read(unit_in,'(i5)') relcc_ccener_max_dimension_diis
      end if

      if (kw_matches(word, '.NTOL  ')) then
         read(unit_in,'(i5)') relcc_ccener_desired_convergence
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccfopr(word, kw_section)

      use relcc_cfg

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.MP2G   ')) then
         relcc_do_mp2gradient = .true.
      end if

      if (kw_matches(word, '.MP2G_O')) then
         relcc_do_oldmp2gradient = .true.
      end if

      if (kw_matches(word, '.CCSDG ')) then
         relcc_do_ccsdgradient = .true.
      end if

      if (kw_matches(word, '.CCSDTG')) then
         relcc_do_ccsdtgradient = .true.
      end if

      if (kw_matches(word, '.NATORB')) then
         relcc_do_naturalorbitals = .true.
      end if

      if (kw_matches(word, '.NOZG  ')) then
         relcc_do_zipgrid = .true.
      end if

      if (kw_matches(word, '.NEOPER')) then
         read(unit_in,'(i5)') relcc_ne_oper
      end if

      if (kw_matches(word, '.MAXIT ')) then
         read(unit_in,'(i5)') relcc_fopr_max_iterations
      end if

      if (kw_matches(word, '.MAXDIM')) then
         read(unit_in,'(i5)') relcc_fopr_max_dimension_diis
      end if

      if (kw_matches(word, '.NTOL  ')) then
         read(unit_in,'(i5)') relcc_fopr_desired_convergence
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccfspc(word, kw_section)

      use relcc_cfg

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------
      integer                          :: i
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.DOEA  ')) then
         relcc_fs_do_ea = .true.
      end if

      if (kw_matches(word, '.DOIE  ')) then
         relcc_fs_do_ie = .true.
      end if

      if (kw_matches(word, '.DOEA2 ')) then
         relcc_fs_do_ea2 = .true.
      end if

      if (kw_matches(word, '.DOIE2 ')) then
         relcc_fs_do_ie2 = .true.
      end if

      if (kw_matches(word, '.DOEXC ')) then
         relcc_fs_do_exc = .true.
      end if

      if (kw_matches(word, '.FSSECT ')) then
         read(unit_in,*) (relcc_fs_fssect(i), i=1,6)
      end if

      if (kw_matches(word, '.DOIH  ')) then
         relcc_fs_do_ih = .true.
      end if

      if (kw_matches(word, '.NACTH ')) then
         call kw_read(word, relcc_fs_nacth, 2*NFSYM)
      end if

      if (kw_matches(word, '.NACTP ')) then
         call kw_read(word, relcc_fs_nactp, 2*NFSYM)
      end if

      if (kw_matches(word, '.MAXIT ')) then
         read(unit_in,'(i5)') relcc_fs_max_iterations
      end if

      if (kw_matches(word, '.MXIT02')) then
         read(unit_in,'(i5)') relcc_fs_max02_iterations
      end if

      if (kw_matches(word, '.MXIT01')) then
         read(unit_in,'(i5)') relcc_fs_max01_iterations
      end if

      if (kw_matches(word, '.MXIT20')) then
         read(unit_in,'(i5)') relcc_fs_max20_iterations
      end if

      if (kw_matches(word, '.MXIT10')) then
         read(unit_in,'(i5)') relcc_fs_max10_iterations
      end if

      if (kw_matches(word, '.MXIT11')) then
         read(unit_in,'(i5)') relcc_fs_max11_iterations
      end if

      if (kw_matches(word, '.MXIT00')) then
         read(unit_in,'(i5)') relcc_fs_max00_iterations
      end if

      if (kw_matches(word, '.MAXDIM')) then
         read(unit_in,'(i5)') relcc_fs_max_dimension_diis
      end if

      if (kw_matches(word, '.NTOL  ')) then
         read(unit_in,'(i5)') relcc_fs_desired_convergence
      end if

      if (kw_matches(word, '.TSHOLD')) then
         read(unit_in,'(i5)') relcc_fs_tshold
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_ccih(word, kw_section)

      use relcc_cfg

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.IHSCHE')) then
         read(unit_in,'(i5)') relcc_ih_scheme
      end if

      if (kw_matches(word, '.SH_H11')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_h11
      end if

      if (kw_matches(word, '.SH_H12')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_h12
      end if

      if (kw_matches(word, '.SH_P11')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_p11
      end if

      if (kw_matches(word, '.SH_P12')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_p12
      end if

      if (kw_matches(word, '.SH_P2 ')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_p2
      end if

      if (kw_matches(word, '.SH_H2 ')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_h2
      end if

      if (kw_matches(word, '.SH_HP ')) then
         read(unit_in,'(f15.8)') relcc_ih_shift_hp
      end if

      if (kw_matches(word, '.AIH   ')) then
         read(unit_in,'(f15.8)') relcc_ih_aih
      end if

      if (kw_matches(word, '.EHMIN ')) then
         read(unit_in,'(f15.8)') relcc_ih_eh_min
      end if

      if (kw_matches(word, '.EHMAX ')) then
         read(unit_in,'(f15.8)') relcc_ih_eh_max
      end if

      if (kw_matches(word, '.EPMIN ')) then
         read(unit_in,'(f15.8)') relcc_ih_ep_min
      end if

      if (kw_matches(word, '.EPMAX ')) then
         read(unit_in,'(f15.8)') relcc_ih_ep_max
      end if

      if (kw_matches(word, '.NIH   ')) then
         read(unit_in,'(i5)') relcc_ih_nih
      end if

      if (kw_matches(word, '.NACTHI')) then
         call kw_read(word, relcc_ih_nacthi, 2*NFSYM)
      end if

      if (kw_matches(word, '.NACTPI')) then
         call kw_read(word, relcc_ih_nactpi, 2*NFSYM)
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

   subroutine read_input_x2c(word, kw_section)

      use x2cmod_cfg

!     --------------------------------------------------------------------------
      character(kw_length), intent(in) :: word
      character(kw_length), intent(in) :: kw_section
!     --------------------------------------------------------------------------

      call reset_available_kw_list()

      if (kw_matches(word, '.FGMx2c')) then
!        x2cmod_fragment_x2c = .true.
         x2cmod_fragment_x2c = .false.
      end if

      if (kw_matches(word, '.DEBUG ')) then
         x2cmod_debug = .true.
      end if

      if (kw_matches(word, '.hdirty')) then
         x2cmod_h1dirty = .true.
      end if

!     free-particle matrix as defining h1 hamiltonian
      if (kw_matches(word, '.freeh1')) then
         x2c_free_part_mtx_defh1 = .true.
      end if

!     use DIRACs in-house AMFI
!     (full keyword: .diramfi)
      if (kw_matches(word, '.diramf')) then
         x2c_add_amfi = 2
      end if
      if (kw_matches(word, '.fockao')) then
         x2c_fock_saao_basis = .true.
      end if

      call check_whether_kw_found(word, kw_section)

   end subroutine

end module
