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

module visual

  use memory_allocator
  use visual_cfg
  use codata
  use ao_eval
  use interface_ao
  use matrix_defop
  use matrix_genop, only: copy_matrix_sym
  use dirac_interface
#ifdef MOD_UNRELEASED
  use visual_london
  use xc_london_c1
#endif
  use visual_integrate
  use visual_in_point


  implicit none

  public run_visual

  private

#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "dcbxlr.h"



    integer, parameter :: file_unit = 137

contains

    subroutine run_visual()

!   ----------------------------------------------------------------------------
    type(matrix)         :: D_mo, D(max_para)
    type(matrix)         :: D_con_ao, D_con_mo, D_con, D_0, D_0_full
    type(matrix)         :: C, Cig, Ciu, Csg, Csu

    real(8), allocatable :: eigenvalues(:)
    real(8), allocatable :: occupation(:)
    real(8), allocatable :: ao(:)
    real(8), allocatable :: buffer(:)

    real(8)              :: energy, time_start, second
    integer              :: max_ao_g_order, max_ao_m_order
    integer              :: ioff, imat, iq, iz, m, ivector

    logical              :: debug_me = .false.

    integer :: i, j, off_ij, off_iz, off_ij_s, off_ij_e
    type(matrix) :: D_mo_dagger, TP, TM, T
!   ----------------------------------------------------------------------------

    max_ao_g_order = 0
    do imat = 1, visual_cfg_nr_dmat
       if (                                             &
           (visual_cfg_property(imat) == iq_elf)   .or. &
           (visual_cfg_property(imat) == iq_divs)  .or. &
           (visual_cfg_property(imat) == iq_rots)  .or. &
           (visual_cfg_property(imat) == iq_divj)  .or. &
           (visual_cfg_property(imat) == iq_rotj)       &
          ) then
          max_ao_g_order = 1
       end if
       if (levyle) then
!         here we obtain the paramagnetic part of the current density
!         via the nr paramagnetic current density expression
!         for this we need the cartesian gradient
          if (                                             &
              (visual_cfg_property(imat) == iq_j)     .or. &
              (visual_cfg_property(imat) == iq_ndipx) .or. &
              (visual_cfg_property(imat) == iq_ndipy) .or. &
              (visual_cfg_property(imat) == iq_ndipz) .or. &
              (visual_cfg_property(imat) == iq_bdipx) .or. &
              (visual_cfg_property(imat) == iq_bdipy) .or. &
              (visual_cfg_property(imat) == iq_bdipz)      &
             ) then
             max_ao_g_order = 1
          end if
       end if
    end do

    max_ao_m_order = 0
    if (visual_cfg_london) then
       max_ao_m_order = 1
    end if

    write(LUPRI,'(/,2x,a)') '*** Visualisation module ***'


!   find occupation
!   ===============

    call alloc(eigenvalues, norbt)
    eigenvalues = 0.0d0

    call opnfil(file_unit, 'DFCOEF', 'OLD', 'visual.F90')
    call reacmo_no_work(file_unit, 'DFCOEF', (/0.0d0/), eigenvalues, (/0/), energy, 4)

    close(file_unit, status = 'keep')
    write(LUPRI,'(2x,a)') 'Visualisation : DFCOEF file read'

    call fndocc(eigenvalues, .true.)
    call dealloc(eigenvalues)


!   get orbital string
!   ==================

    call alloc(occupation, norbt)
    call setup_orbital_string(occupation)


!   get coefficients

    call get_c(C,   i=1.0d0, s=1.0d0, g=1.0d0, u=1.0d0)
    call get_c(Cig, i=1.0d0, s=0.0d0, g=1.0d0, u=0.0d0)
    call get_c(Csg, i=0.0d0, s=1.0d0, g=1.0d0, u=0.0d0)
    if (debug_me) then
       call print_mat(C,   label='debug C')
       call print_mat(Cig, label='debug Cig')
       call print_mat(Csg, label='debug Csg')
    end if

    if (nfsym == 2) then
       call get_c(Ciu, i=1.0d0, s=0.0d0, g=0.0d0, u=1.0d0)
       call get_c(Csu, i=0.0d0, s=1.0d0, g=0.0d0, u=1.0d0)
       if (debug_me) then
          call print_mat(Ciu, label='debug Ciu')
          call print_mat(Csu, label='debug Csu')
       end if
    end if

    call init_mat(D_mo, norbt, norbt)

#ifdef MOD_UNRELEASED
    if (visual_cfg_london) then

       D_mo%elms = 0.0d0
       do m = 1, norbt
          D_mo%elms((m - 1)*norbt + m) = occupation(m)
       end do
       D(imat) = C*(D_mo*dag(C))

       D_0 = C*(D_mo*dag(C))
       D_0 = 2.0d0*D_0
       D_0_full = 2.0d0*C*dag(C)
       
       call init_mat(D_con_ao, ntbas(0), ntbas(0))
       D_con_ao%elms = 0.0d0
       
!      read connection matrices in AO basis from AOPROPER file
       call setup_fo_conmat_ao(D_con_ao, 1)
       
       
!      transform connection matrices to MO basis
       D_con_mo = dag(C)*(D_con_ao*C)

       if (.not. levyle) then
!         radovan: magic factor for levyle
!                  need to find the reason for it
          D_con_mo = 0.5d0*D_con_mo
       end if

       if (visual_cfg_london_symmetric) then
!         symmetric connection
       
          D_con = D_0_full*D_con_ao*D_0 &
                + D_0*D_con_ao*D_0_full
          D_con = 0.125d0*D_con
       
       else if (visual_cfg_london_natural) then
!         natural connection
       
          D_con = D_0_full*D_con_ao*D_0 &
                - D_0*dag(D_con_ao)*D_0_full
          D_con = 0.125d0*D_con
       
       else if (visual_cfg_london_none) then
!         none connection
       
          D_con = D_0*D_con_ao*D_0
          D_con = 0.25d0*D_con
       end if

       if (levyle) then
          D_con = 2.0d0*D_con
       end if
    end if
#endif

    do imat = 1, visual_cfg_nr_dmat

       if (visual_cfg_dmat_file(imat) == 'DFCOEF') then
!         get unperturbed density matrix
      
          D_mo%elms = 0.0d0
          do m = 1, norbt
             D_mo%elms((m - 1)*norbt + m) = occupation(m)
          end do

          D_mo%ih_sym = 1
          D_mo%tr_sym = 1
          D_mo%irep   = 0

          D(imat) = C*(D_mo*dag(C))
          if (debug_me) then
             call print_mat(D_mo,    label = 'debug MO D matrix')
             call print_mat(D(imat), label = 'debug AO D matrix')
          end if
      
       else
!         get perturbed density matrix
      
          call setup_fo_dmat_ao(imat, occupation, D_mo)

          if (visual_cfg_london .and. visual_cfg_london_none) then
             TP = D_mo + D_con_mo
             TM = D_mo - D_con_mo
          else
             TP = 1.0d0*D_mo
             TM = 1.0d0*D_mo
          end if

          if (nfsym == 2) then
             if (jbtof(D_mo%irep, 1) == 2) then
!               ungerade perturbation
                D(imat) = (Cig*(TP*dag(Csu))) &
                        + (Ciu*(TP*dag(Csg))) &
                        - (Csg*(TM*dag(Ciu))) &
                        - (Csu*(TM*dag(Cig)))
             else
!               gerade perturbation
                D(imat) = (Cig*(TP*dag(Csg))) &
                        + (Ciu*(TP*dag(Csu))) &
                        - (Csg*(TM*dag(Cig))) &
                        - (Csu*(TM*dag(Ciu)))
             end if
          else
             D(imat) = (Cig*(TP*dag(Csg))) &
                     - (Csg*(TM*dag(Cig)))
          end if
          call copy_matrix_sym(D(imat), D_mo)
      
       end if
    end do

!   insert half phases
    if (nz < 4) then
       do iz = 1, nz
         ioff = ntbas(0)*ntbas(0)*(iz - 1)
         if (visual_cfg_london) then
            iq   = ipqtoq(iz, D_0%irep)
            call q2bphase('D', iq, 1, D_0%elms(1 + ioff))
            iq   = ipqtoq(iz, D_con%irep)
            call q2bphase('D', iq, 1, D_con%elms(1 + ioff))
         end if
         do imat = 1, visual_cfg_nr_dmat
            iq   = ipqtoq(iz, D(imat)%irep)
            call q2bphase('D', iq, 1, D(imat)%elms(1 + ioff))
         end do
       end do
    end if

    do imat = 1, visual_cfg_nr_dmat
#ifdef MOD_UNRELEASED
       if (visual_cfg_london) then
          if (.not. visual_cfg_london_skip_ro) then
!            dirty trick add D_con on top of D
!            on trunk we will have to check that this only
!            happes to the current, mag field, and so on
             D(imat) = D(imat) + D_con
          end if
       end if
#endif

!      scale density matrix
       D(imat) = 2.0d0*D(imat)

!      collect matrices for identical densities to save time
       if (visual_cfg_property(imat) /= iq_elf) then
!         elf has nonlinear dependence
          j = visual_cfg_property_to_matrix(visual_cfg_property(imat))
          if (j == 0) then
             visual_cfg_property_to_matrix(visual_cfg_property(imat)) = imat
          else
             D(j) = D(j) + D(imat)
             if (visual_cfg_london) then
                T   = 1.0d0*D_0
                D_0 = D_0 + T
                T   = 0
             end if
             visual_cfg_skip(imat) = .true.
          end if
       end if
    end do


!   deallocate various things
!   =========================

    call dealloc(occupation)

    C    = 0
    Cig  = 0
    Csg  = 0

    if (nfsym == 2) then
      Ciu = 0
      Csu = 0
    end if


!   initialize ao module
!   ====================

    call ao_eval_init(max_ao_g_order, &
                      max_ao_m_order, &
                      .false.)

    call alloc(ao,     nr_ao_slices*nr_ao_cartesian)
    call alloc(buffer, nr_ao_slices*nr_ao_cartesian)

#ifdef MOD_UNRELEASED
    if (visual_cfg_london) then
       D_con_ao = 0
       D_con_mo = 0
    end if
#endif


!   plot and integrate
!   ==================

    time_start = second()

    call visual_driver(visual_cfg_nr_dmat, &
                       D,                  &
                       D_0,                &
                       ao,                 &
                       buffer)

    call timtxt('>>>> time used in plot:', (second() - time_start), lupri)

    time_start = second()
    if (visual_cfg_2d_integration .or. visual_cfg_3d_integration) then
       call visual_integration(visual_cfg_nr_dmat, &
                               D,                  &
                               D_0,                &
                               ao,                 &
                               buffer)
    end if
    call timtxt('>>>> time used in integration:', (second() - time_start), lupri)

    call dealloc(ao)
    call dealloc(buffer)


!   deallocate density matrices

    do imat = 1, visual_cfg_nr_dmat
       D(imat) = 0
    end do
#ifdef MOD_UNRELEASED
    if (visual_cfg_london) then
       D_con = 0
    end if
#endif

    end subroutine




   subroutine visual_driver(nr_dmat, &
                            D,       &
                            D_0,     &
                            ao,      &
                            buffer)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: nr_dmat
      type(matrix)         :: D(nr_dmat)
      type(matrix)         :: D_0
      real(8), intent(in)  :: ao(*)
      real(8), intent(in)  :: buffer(*)
!     --------------------------------------------------------------------------
      real(8), allocatable :: q(:, :)
      real(8), allocatable :: atom_xyz(:, :)
      real(8)              :: step_3d(3, 0:3)
      integer              :: i, j, k, l, m, istep
      integer              :: icube, jcube, kcube
      integer              :: ix, iy
      real(8)              :: vx(3), vx_norm
      real(8)              :: vy(3), vy_norm
      real(8)              :: r(3), v(3), s
      integer              :: lunit_2d_scalar
      integer              :: lunit_2d_vector
      integer              :: lunit_3d_cube
      integer              :: lunit_3d_scalar
      integer              :: lunit_3d_vector
      integer              :: lunit_line_scalar
      integer              :: lunit_line_vector
!     --------------------------------------------------------------------------

      i = 1
      if (visual_cfg_ncube(3) > 0) then
         i = visual_cfg_ncube(3)
      end if
      allocate(q(0:3, i))

      allocate(atom_xyz(4, nucdep))

      lunit_2d_scalar   = 14
      lunit_2d_vector   = 15
      lunit_3d_cube     = 16
      lunit_3d_scalar   = 17
      lunit_3d_vector   = 18
      lunit_line_scalar = 19
      lunit_line_vector = 20

      if (visual_cfg_2d) then
         open(lunit_2d_scalar,           &
              file   = 'plot.2d.scalar', &
              status = 'new',            &
              access = 'sequential',     &
              form   = 'formatted')
         open(lunit_2d_vector,           &
              file   = 'plot.2d.vector', &
              status = 'new',            &
              access = 'sequential',     &
              form   = 'formatted')
      end if

      if (visual_cfg_3d) then
         open(lunit_3d_cube,           &
              file   = 'plot.3d.cube', &
              status = 'new',          &
              access = 'sequential',   &
              form   = 'formatted')
         open(lunit_3d_vector,           &
              file   = 'plot.3d.vector', &
              status = 'new',            &
              access = 'sequential',     &
              form   = 'formatted')
         open(lunit_3d_scalar,           &
              file   = 'plot.3d.scalar', &
              status = 'new',            &
              access = 'sequential',     &
              form   = 'formatted')
      end if

      if (visual_cfg_line) then
         open(lunit_line_scalar,           &
              file   = 'plot.line.scalar', &
              status = 'new',              &
              access = 'sequential',       &
              form   = 'formatted')
         open(lunit_line_vector,           &
              file   = 'plot.line.vector', &
              status = 'new',              &
              access = 'sequential',       &
              form   = 'formatted')
      end if

      if (visual_cfg_list) then

         call around('Requested density at point(s)')

         write(lupri, *)
         write(lupri, '(a8, 2a14, a18, a21, 2a18)') &
            'r_x',                                  &
            'r_y',                                  &
            'r_z',                                  &
            'scalar',                               &
            'x-component',                          &
            'y-component',                          &
            'z-component'
         write(lupri, *)

         do i = 1, visual_cfg_nr_points_in_list
            call get_quantity_in_point(visual_cfg_nr_dmat,        &
                                       D,                         &
                                       D_0,                       &
                                       visual_cfg_xyz_list(i, 1), &
                                       visual_cfg_xyz_list(i, 2), &
                                       visual_cfg_xyz_list(i, 3), &
                                       ao,                        &
                                       buffer,                    &
                                       q)

            write(lupri, '(3e14.6, 4e18.10)') &
               visual_cfg_xyz_list(i, 1),     &
               visual_cfg_xyz_list(i, 2),     &
               visual_cfg_xyz_list(i, 3),     &
               q(0:3, 1)
         end do

         write(lupri, *)
         call prsymb(lupri, '-', 100, 0)
         write(lupri, *)
         write(lupri, *)

      end if

      if (visual_cfg_line) then

         v = visual_cfg_line_to - visual_cfg_line_from
         v = v/visual_cfg_line_nr_steps

         do istep = 0, visual_cfg_line_nr_steps

            r = visual_cfg_line_from + istep*v
            call get_quantity_in_point(visual_cfg_nr_dmat, &
                                       D,                  &
                                       D_0,                &
                                       r(1),               &
                                       r(2),               &
                                       r(3),               &
                                       ao,                 &
                                       buffer,             &
                                       q)

            write(lunit_line_scalar, '(4e18.10)') &
               r(1),                              &
               r(2),                              &
               r(3),                              &
               q(0, 1)
            write(lunit_line_vector, '(6e18.10)') &
               r(1),                              &
               r(2),                              &
               r(3),                              &
               q(1:3, 1)
        
         end do
      end if

      if (visual_cfg_2d) then

         vx = (visual_cfg_2d_p_right - visual_cfg_2d_p_origin)/visual_cfg_2d_nr_right
         vy = (visual_cfg_2d_p_top   - visual_cfg_2d_p_origin)/visual_cfg_2d_nr_top

         vx_norm = dsqrt(vx(1)*vx(1) + vx(2)*vx(2) + vx(3)*vx(3))
         vy_norm = dsqrt(vy(1)*vy(1) + vy(2)*vy(2) + vy(3)*vy(3))

         i = 0
         do ix = 0, visual_cfg_2d_nr_right
            do iy = 0, visual_cfg_2d_nr_top
               i = i + 1
               r = visual_cfg_2d_p_origin + ix*vx + iy*vy
               call get_quantity_in_point(visual_cfg_nr_dmat, &
                                          D,                  &
                                          D_0,                &
                                          r(1),               &
                                          r(2),               &
                                          r(3),               &
                                          ao,                 &
                                          buffer,             &
                                          q)

               s      = q(0,   1)
               v(1:3) = q(1:3, 1)

               write(lunit_2d_scalar, '(3e18.10)')               &
                  ix*vx_norm - vx_norm*visual_cfg_2d_nr_right/2, &
                  iy*vy_norm - vy_norm*visual_cfg_2d_nr_top/2,   &
                  s

               write(lunit_2d_vector, '(4e18.10)')                &
                  ix*vx_norm - vx_norm*visual_cfg_2d_nr_right/2,  &
                  iy*vy_norm - vy_norm*visual_cfg_2d_nr_top/2,    &
                  (v(1)*vx(1) + v(2)*vx(2) + v(3)*vx(3))/vx_norm, &
                  (v(1)*vy(1) + v(2)*vy(2) + v(3)*vy(3))/vy_norm
            end do
         end do
      end if

      if (visual_cfg_3d) then

!        step_3d(*, 0)  coordinates of initial point
!        step_3d(*, 1)  step vector for the "x"-direction
!        step_3d(*, 2)  step vector for the "y"-direction
!        step_3d(*, 3)  step vector for the "z"-direction

!        atom_xyz(1, i) x-coordinate of atom i
!        atom_xyz(2, i) y-coordinate of atom i
!        atom_xyz(3, i) z-coordinate of atom i
!        atom_xyz(4, i) charge of atom i

         step_3d = 0.0d0

         do i = 1, 3
            step_3d(i, 0) = cord(i, 1)
            step_3d(i, i) = cord(i, 1)
         enddo

         i = 0
         do j = 1, nucind
            k = istbnu(j)
            do l = 0, maxopr
               if (iand(l, k) == 0) then
                  i = i + 1
                  do m = 1, 3
                     atom_xyz(m, i) = pt(iand(isymax(m, 1), l))*cord(m, j)

                     step_3d(m, 0) = min(step_3d(m, 0), atom_xyz(m, i))
                     step_3d(m, m) = max(step_3d(m, m), atom_xyz(m, i))
                  end do
                  atom_xyz(4, i) = charge(j)
               end if
            end do
         end do

!        add some space, calculate step sizes and convert step vectors to angstrom

         do i = 1, 3
            step_3d(i, 0) =  step_3d(i, 0) -  visual_cfg_add_3d
            step_3d(i, i) = (step_3d(i, i) + (visual_cfg_add_3d*2.0d0))*xtang
            step_3d(i, i) = (step_3d(i, i) -  step_3d(i, 0))/(visual_cfg_ncube(i) - 1)
         end do

!        write header

         write(lunit_3d_cube, '(a)') 'cube file'
         write(lunit_3d_cube, '(a)') 'generated by DIRAC'

         write(lunit_3d_cube,'(i5, 4f12.6)') &
            nucdep,                          &
           (step_3d(j, 0), j = 1, 3)
         do i = 1, 3
           write(lunit_3d_cube,'(i5, 4f12.6)') &
              visual_cfg_ncube(i),             &
             (step_3d(j, i), j = 1, 3)
         end do

         do i = 1, nucdep
           write(lunit_3d_cube, '(i5, 4f12.6)') &
                nint(atom_xyz(4, i)),           &
                atom_xyz(4, i),                 &
               (atom_xyz(j, i), j = 1, 3)
         end do

         do icube = 1, visual_cfg_ncube(1)
            do jcube = 1, visual_cfg_ncube(2)
               do kcube = 1, visual_cfg_ncube(3)

                  r(1) = step_3d(1, 0) + (icube-1)*step_3d(1, 1) &
                                       + (jcube-1)*step_3d(1, 2) &
                                       + (kcube-1)*step_3d(1, 3)
                  r(2) = step_3d(2, 0) + (icube-1)*step_3d(2, 1) &
                                       + (jcube-1)*step_3d(2, 2) &
                                       + (kcube-1)*step_3d(2, 3)
                  r(3) = step_3d(3, 0) + (icube-1)*step_3d(3, 1) &
                                       + (jcube-1)*step_3d(3, 2) &
                                       + (kcube-1)*step_3d(3, 3)

                  call get_quantity_in_point(visual_cfg_nr_dmat, &
                                             D,                  &
                                             D_0,                &
                                             r(1),               &
                                             r(2),               &
                                             r(3),               &
                                             ao,                 &
                                             buffer,             &
                                             q(0, kcube))

                  write(lunit_3d_vector, '(6e13.5)') &
	             r(1:3),                         &
	             q(1:3, kcube)
                  write(lunit_3d_scalar, '(4e13.5)') &
	             r(1:3),                         &
	             q(0, kcube)
               end do
               write(lunit_3d_cube, '(6e13.5)')  &
                  (q(0, m), m = 1, visual_cfg_ncube(3))
            end do
         end do
      end if

      if (visual_cfg_2d)   close(lunit_2d_scalar, status = 'keep')
      if (visual_cfg_2d)   close(lunit_2d_vector, status = 'keep')
      if (visual_cfg_3d)   close(lunit_3d_cube,   status = 'keep')
      if (visual_cfg_3d)   close(lunit_3d_scalar, status = 'keep')
      if (visual_cfg_3d)   close(lunit_3d_vector, status = 'keep')
      if (visual_cfg_line) close(lunit_line_scalar,     status = 'keep')
      if (visual_cfg_line) close(lunit_line_vector,     status = 'keep')

      deallocate(q)
      deallocate(atom_xyz)

   end subroutine


      subroutine setup_orbital_string(occupation)

      real(8)              :: occupation(norbt)

      integer, allocatable :: kvec(:)
      integer              :: i, ifsym, nvec, ioff, ishift

      occupation = 0.0d0

      ioff = 0

      do ifsym = 1, nfsym

        ioff = ioff + npsh(ifsym)

        if (visual_cfg_use_orbital_string) then

          call alloc(kvec, norb(ifsym))
          kvec = 0
          nvec = 1

          call numlst(visual_cfg_orbital_string(ifsym), &
     &                kvec, &
     &                norb(ifsym), &
     &               -npsh(ifsym), &
     &                nesh(ifsym), &
     &                ifsym, &
     &                nvec)

          do i = 1, nvec
            ishift = kvec(i)
            if (ishift < 0) ishift = ishift + 1
            occupation(ioff + ishift) = 1.0d0
          end do

          call dealloc(kvec)

        else

          do i = 1, nish(ifsym)
            occupation(ioff + i) = 1.0d0
          end do

        end if

        ioff = ioff + nesh(ifsym)

      end do

      end subroutine

      SUBROUTINE SETUP_FO_DMAT_AO(IXVECTOR, &
     &                            STRING, D_mo)

      CHARACTER LAB1*16,&
     &          LAB2*16,&
     &          TYP1*2,&
     &          TYP2*2,&
     &          LBUF*16,&
     &          PBUF*1
      LOGICAL   FILE_EXISTS
      type(matrix) :: D_mo
      real(8) :: STRING(*)

      real(8), allocatable :: bvec(:)

      integer, allocatable :: orbital_rotation_indices(:, :)
      real(8), allocatable :: response_vector(:, :)
      integer              :: length
      character(2)         :: vector_type
      real(8)              :: ih_factor

      integer              :: nr_vectors, LUNIT_XVECTOR
      integer              :: ixvector, NBAS, NORBS, len
      integer              :: JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, IVECTOR
      real(8)              :: FREQ1, FREQ2, ERGRSP, unused

!gosia, i'll remove it later
      integer :: i, j
      character(10) :: pattern
      character(16) :: text

      INQUIRE(FILE = visual_cfg_dmat_file(IXVECTOR),EXIST = FILE_EXISTS)
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A,A6,A)') 'solution vector file ', &
     &    visual_cfg_dmat_file(IXVECTOR),' not found'
        CALL QUIT('ERROR in VISUAL')
      ENDIF
      LUNIT_XVECTOR = 17
      OPEN(LUNIT_XVECTOR, &
     &     FILE = visual_cfg_dmat_file(IXVECTOR), &
     &     STATUS = 'OLD', &
     &     ACCESS = 'SEQUENTIAL', &
     &     FORM = 'UNFORMATTED')

!     get number of vectors on file

      NR_VECTORS = 0
      REWIND LUNIT_XVECTOR
 1    CONTINUE
      READ(LUNIT_XVECTOR) LAB1
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 2
      ELSE
        NR_VECTORS = NR_VECTORS + 1
!       skip solution vector
        READ(LUNIT_XVECTOR)
!       skip orbital rotations
        READ(LUNIT_XVECTOR)
        GOTO 1
      ENDIF
 2    CONTINUE

!     errors

      if ((visual_cfg_dmat_file_record(ixvector) > nr_vectors) &
          .or. (visual_cfg_dmat_file_record(ixvector) < 1)) then
         write(lupri,'(a)') 'error in visual:'
         write(lupri,'(a)') 'you try to read an invalid record number from file '//visual_cfg_dmat_file(ixvector)
         call quit('error in visual')
      end if

!     move to the right record

      IVECTOR = 0
      REWIND LUNIT_XVECTOR
 3    CONTINUE
      READ(LUNIT_XVECTOR)
      IVECTOR = IVECTOR + 1
      IF(IVECTOR .EQ. visual_cfg_dmat_file_record(IXVECTOR)) THEN
        BACKSPACE(LUNIT_XVECTOR)
        GOTO 4
      ELSE
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        GOTO 3
      ENDIF
 4    CONTINUE


!     get type and dimensions of record IWHICH_XVECTOR(IXVECTOR)

      READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2, &
     &                    JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2, &
     &                    UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS



      visual_cfg_irep_xvector(IXVECTOR) = JSYMOP1 - 1
      length  = LEN/NZ
      D_mo%irep   = JSYMOP1 - 1


    select case (typ1(1:1))
      case ('E')
        vector_type = 'pp'
      case ('P')
        vector_type = 'pn'
      case default
        call quit('unexpected typ1(1:1)')
    end select

    select case (typ1(2:2))
      case ('+')
        D_mo%ih_sym =  1
        ih_factor   =  1.0d0
      case ('-')
        D_mo%ih_sym = -1
        ih_factor   = -1.0d0
      case default
        call quit('unexpected typ1(2:2)')
    end select


!   allocate and read response vector
    call alloc(response_vector, length, nz)
    call readt(lunit_xvector, length*nz, response_vector)
!   allocate and read orbital rotation array
    call alloc(orbital_rotation_indices, 2, length)
    call readi(lunit_xvector, length*2, orbital_rotation_indices)
!   close file of solution vectors
    close(lunit_xvector, status = 'keep')

!fixme: this is a good place to zero elements according to occupation

    D_mo%elms = 0.0d0
    call scatter_vector(length,                   &
                        orbital_rotation_indices, &
                        ih_factor,                &
                        response_vector,          &
                        D_mo%elms,                &
                        D_mo%irep)

    call dealloc(orbital_rotation_indices)
    call dealloc(response_vector)

  end subroutine

  subroutine scatter_vector(length,                   &
                            orbital_rotation_indices, &
                            h,                        &
                            response_vector,          &
                            matrix,                   &
                            irep)

!   ----------------------------------------------------------------------------
    integer, intent(in)  :: length
    integer, intent(in)  :: orbital_rotation_indices(2, length)
    real(8), intent(in)  :: h
    real(8), intent(in)  :: response_vector(length, nz)
    real(8), intent(out) :: matrix(norbt, norbt, nz)
    integer, intent(in)  :: irep
!   ----------------------------------------------------------------------------
    integer              :: i, s, is, iz
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do is = 1, length

      i = orbital_rotation_indices(1, is)
      s = orbital_rotation_indices(2, is)

      do iz = 1, nz

        if (ipqtoq(iz, irep) > 1) then
          f = -1.0
        else
          f =  1.0
        end if

!              row
!              ¦
!              ¦  column
!              ¦  ¦
        matrix(s, i, iz) = matrix(s, i, iz) &
                         +     response_vector(is, iz)
        matrix(i, s, iz) = matrix(i, s, iz) &
                         - f*h*response_vector(is, iz)
      end do
    end do

  end subroutine

end module
