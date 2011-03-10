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

module ao_eval

!   this is based on DIRAC and Dalton
!   routines written by T. Helgaker, P. Salek, T. Saue, and others

! - in order to calculate aos using this module
!   you need to initialize some data first
!   by calling ao_eval_init from each node
! - ao_eval_init will set symmetry and basis set information
!   that is needed in each point and then mpi_bcast this info to the nodes
! - we try to not recalculate data which does not change from point to point
!   (like symmetry)
! - *after* the module data is initialized you can allocate
!   the ao vector and a buffer, both with the length total_ao_vector_length
! - (de)allocate them only once (each iteration), not in every point
!   this is also the reason why the ao vector and buffer are not allocated here - this would
!   be very inefficient
! - aos are not dzeroed for more performance
!   if you change the code make sure that elements of these vectors
!   are not read before they are written for the first time
!   note that in levy-leblond or x2c4 mode some elements will
!   contain garbage but they are not (= should not be) used
! - this module is called dft_ao but should be usable also outside dft
!   this is why some information is passed through ao_eval_init and not queried
!   from dft_cfg, to keep this module generic, please dont "use dft_cfg" here

!   some "bad ideas" = ideas which did not improve things
!   so that future generations do not waste time on this:
! - ao to so transformation using dgemv: possible, but much slower than present scheme
! - derivatives of aos by numerical differentiation: again, turned out to be much slower

!   FAQ
!   ===
!
!   1. how is the ao array sorted on output of get_ao?
!      - it is sorted according to ibbas and nbbas
!        a convenient way to see the sorting of (symmetry adapted) aos
!        in the output is to search for 'Generating Lowdin matrix'
!        as an example, with a d2h calculation i got:
!
!                                 Generating Lowdin matrix:
!                                 -------------------------
!
!        L   Ag    * Deleted:          2(Proj:          2, Lindep:          0)
!        L   B1g   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B2g   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B3g   * Deleted:          0(Proj:          0, Lindep:          0)
!        S   B3u   * Deleted:          2(Proj:          2, Lindep:          0)
!        S   B2u   * Deleted:          2(Proj:          2, Lindep:          0)
!        S   B1u   * Deleted:          6(Proj:          6, Lindep:          0)
!        S   Au    * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B3u   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B2u   * Deleted:          0(Proj:          0, Lindep:          0)
!        L   B1u   * Deleted:          1(Proj:          1, Lindep:          0)
!        L   Au    * Deleted:          0(Proj:          0, Lindep:          0)
!        S   Ag    * Deleted:          9(Proj:          9, Lindep:          0)
!        S   B1g   * Deleted:          0(Proj:          0, Lindep:          0)
!        S   B2g   * Deleted:          1(Proj:          1, Lindep:          0)
!        S   B3g   * Deleted:          1(Proj:          1, Lindep:          0)
!
!        so the order is:
!        | L Ag | L B1g | L B2g | L B3g | S B3u | ... | S B3g | and then the derivatives


   use basis_info
   use block_config
   use ao_xyz_power
 
   implicit none
 
   public get_ao
   public ao_eval_init

   private

   interface get_ao
      module procedure get_ao_point
      module procedure get_ao_block
   end interface

   integer, parameter :: max_order = 2
   integer, public    :: total_ao_vector_length = 0
   integer, public    :: ao_start     (0:max_order, 0:max_order) = 0
   integer, public    :: ao_end       (0:max_order, 0:max_order) = 0
   logical, public    :: need_ao_order(0:max_order, 0:max_order) = .false.
!  standard LDA calculation: nr_ao_slices = 1
!                            ao vector composition = | ao_g0 |
!                            meaning there are no geometric derivatives
!  standard GGA calculation: nr_ao_slices = 4
!                            ao vector composition = | ao_g0 | ao_gx | ao_gy | ao_gz |
!  with LAO there are more slices
!  also with higher geometric derivatives (soon) there are more and more slices
!  (higher derivatives and mixed derivatives)
!  slices start at ao_start and end at ao_end
!  and the total length of the ao vector is total_ao_vector_length
!  buffer has the same length for convenience (one could get away with less for the buffer)
   integer            :: nr_ao_slices = 0
   logical            :: use_ao_screening = .false.

contains

   subroutine get_ao_point(cx, &
                           cy, &
                           cz, &
                           ao, &
                           buffer)

!     --------------------------------------------------------------------------
      real(8), intent(in)    :: cx
      real(8), intent(in)    :: cy
      real(8), intent(in)    :: cz
      real(8), intent(inout) :: ao(*)
      real(8), intent(inout) :: buffer(*)
!     --------------------------------------------------------------------------

      call get_ao_point_sub(cx, &
                            cy, &
                            cz, &
                            ao)

!     without symmetry skip the rest
      if (nr_boson_sym == 1) return

!     transform from ao(bu) to so(bs)
      call ao_to_so_point(ao, &
                          buffer)

   end subroutine

   subroutine get_ao_block(block_length, &
                           cx,           &
                           cy,           &
                           cz,           &
                           ao,           &
                           buffer)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: cx(*)
      real(8), intent(in)    :: cy(*)
      real(8), intent(in)    :: cz(*)
      real(8), intent(inout) :: ao(block_length, *)
      real(8), intent(inout) :: buffer(block_length, *)
!     --------------------------------------------------------------------------

      call get_ao_block_sub(block_length, &
                            cx,           &
                            cy,           &
                            cz,           &
                            ao)

!     without symmetry skip the rest
      if (nr_boson_sym == 1) return

!     transform from ao(bu) to so(bs)
      call ao_to_so_block(block_length, &
                          ao,           &
                          buffer)

   end subroutine



  subroutine get_ao_point_sub(cx,   &
                    cy,   &
                    cz,   &
                    ao)
   use ao_eval_point
!   ----------------------------------------------------------------------------
    real(8), intent(in) :: cx
    real(8), intent(in) :: cy
    real(8), intent(in) :: cz
    real(8)             :: ao(*)
!   ----------------------------------------------------------------------------
    integer             :: i, iend, ishell, isym, ioff, j
    real(8)             :: px, py, pz, p2, f1, f2
    integer             :: ofg, ofm, off_power

!   ----------------------------------------------------------------------------

    iend = nr_shells


    i = 1
    do ishell = 1, iend
      do isym = 0, nr_boson_sym - 1
        if (iand(isym, shell_stabilizer(ishell)) == 0) then

          px = cx - bitstring_parity(iand(sym_axes(1, 1), isym))*shell_center_xyz(1, ishell)
          py = cy - bitstring_parity(iand(sym_axes(2, 1), isym))*shell_center_xyz(2, ishell)
          pz = cz - bitstring_parity(iand(sym_axes(3, 1), isym))*shell_center_xyz(3, ishell)

          p2 = px*px + py*py + pz*pz

!         note that in contrast to dalton:
!       - shell_extent contains already the half_cell_diagonal
!       - is squared (for more performance)
!         if (                             &
!             (.not. use_ao_screening)     &
!             .or.                         &
!             (p2 <= shell_extent(ishell)) &
!            ) then
!         else zero out!


            g_00_00_00_000 = 0.0d0
            if (need_ao_order(1, 0)) then
               g_00_00_00_100 = 0.0d0
               g_00_00_00_010 = 0.0d0
               g_00_00_00_001 = 0.0d0
            end if
            do j = 1, nr_nonzero_contraction_coefs(ishell)
               f1 = contraction_coef(ishell, j)*dexp(-primitive_exp(ishell, j)*p2)
               g_00_00_00_000 = g_00_00_00_000 + f1
               if (need_ao_order(1, 0)) then
                  f2 = -2.0d0*f1*primitive_exp(ishell, j)
                  g_00_00_00_100 = g_00_00_00_100 + px*f2
                  g_00_00_00_010 = g_00_00_00_010 + py*f2
                  g_00_00_00_001 = g_00_00_00_001 + pz*f2
               end if
            end do

!           this is because get_gao1 does both ao and ao1
            if (need_ao_order(0, 0) .and. .not. need_ao_order(1, 0)) then
              off_power = power_offset(l_quantum_nr(ishell))
              call get_gao(ishell,             &
                           ao(i),              &
                           px,                 &
                           py,                 &
                           pz,                 &
                           x_power(off_power), &
                           y_power(off_power), &
                           z_power(off_power))
            else if (need_ao_order(1, 0)) then
              ! ao_start(1, 0) indicate the place in ao where
              ! geom derivs should start
              ofg = ao_start(1, 0) + i - 1
              off_power = power_offset(l_quantum_nr(ishell))
              call get_gao1(ishell,             &
                            ao(i),              &
                            ao(ofg),            &
                            ao(ofg + nr_ao),    &
                            ao(ofg + 2*nr_ao),  &
                            px,                 &
                            py,                 &
                            pz,                 &
                            p2,                 &
                            x_power(off_power), &
                            y_power(off_power), &
                            z_power(off_power))
            end if

            if (need_ao_order(2, 0)) then
               ofg = ao_start(2, 0) + i - 1
               off_power = power_offset(l_quantum_nr(ishell))
               call get_gao2(ishell,             &
                             ao(ofg          ),  &
                             ao(ofg + 1*nr_ao),  &
                             ao(ofg + 2*nr_ao),  &
                             ao(ofg + 3*nr_ao),  &
                             ao(ofg + 4*nr_ao),  &
                             ao(ofg + 5*nr_ao),  &
                             px,                 &
                             py,                 &
                             pz,                 &
                             p2,                 &
                             x_power(off_power), &
                             y_power(off_power), &
                             z_power(off_power))
            end if

            if (need_ao_order(0, 1)) then
              ! ao_start(0, 1) indicate the place in ao where
              ! magn dervs should start
              ofm = ao_start(0, 1) + i - 1
              call get_ao_b1(px, py, pz,         &
                             ishell,             &
                             ao(i),           &
                             ao(ofm),         &
                             ao(ofm + nr_ao),  &
                             ao(ofm + 2*nr_ao))
              if (need_ao_order(1, 1)) then
                ofm = ao_start(1, 1) + i - 1
                ofg = ao_start(1, 0) + i - 1
                call get_ao_b1(px, py, pz,         &
                               ishell,                 &
                               ao(ofg),           &
                               ao(ofm),         &
                               ao(ofm + nr_ao),  &
                               ao(ofm + 2*nr_ao))
                ofm = ofm + nr_ao*3
                ofg = ofg + nr_ao
                call get_ao_b1(px, py, pz,         &
                               ishell,                 &
                               ao(ofg),           &
                               ao(ofm),         &
                               ao(ofm + nr_ao),  &
                               ao(ofm + 2*nr_ao))
                ofm = ofm + nr_ao*3
                ofg = ofg + nr_ao
                call get_ao_b1(px, py, pz,         &
                               ishell,                 &
                               ao(ofg),           &
                               ao(ofm),         &
                               ao(ofm + nr_ao),  &
                               ao(ofm + 2*nr_ao))
              end if
            end if

!>          if spherical and d or higher then transform
!> \fixme   spherical lao: nothing is done for them now
!> \todo    this should happen outside this routine
            if (shell_is_spherical(ishell) .and. (l_quantum_nr(ishell) > 2)) then
               call transform_c_to_s(ishell, ao(i))
               if (need_ao_order(1, 0)) then
                  ioff = ao_start(1, 0) + i - 1
                  call transform_c_to_s(ishell, ao(ioff          ))
                  call transform_c_to_s(ishell, ao(ioff +   nr_ao))
                  call transform_c_to_s(ishell, ao(ioff + 2*nr_ao))
               end if
            end if


          i = i + spherical_deg(ishell)

        end if
      end do
    end do


  end subroutine


   subroutine get_ao_block_sub(block_length, &
                               cx,           &
                               cy,           &
                               cz,           &
                               ao)

      use ao_eval_block
!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(in)    :: cx(*)
      real(8), intent(in)    :: cy(*)
      real(8), intent(in)    :: cz(*)
      real(8), intent(inout) :: ao(block_length, *)
!     --------------------------------------------------------------------------
      integer                :: i, ishell, isym, j, k
      real(8)                :: px(mxbllen)
      real(8)                :: py(mxbllen)
      real(8)                :: pz(mxbllen)
      real(8)                :: p2(mxbllen)
      real(8)                :: sx
      real(8)                :: sy
      real(8)                :: sz
      real(8)                :: f1, f2
    integer             :: ofg, ofm, off_power
!     --------------------------------------------------------------------------

      i = 1
      do ishell = 1, nr_shells
         do isym = 0, nr_boson_sym - 1
            if (iand(isym, shell_stabilizer(ishell)) == 0) then
      
               sx = bitstring_parity(iand(sym_axes(1, 1), isym))*shell_center_xyz(1, ishell)
               sy = bitstring_parity(iand(sym_axes(2, 1), isym))*shell_center_xyz(2, ishell)
               sz = bitstring_parity(iand(sym_axes(3, 1), isym))*shell_center_xyz(3, ishell)
    
               do k = 1, block_length
                  px(k) = cx(k) - sx
                  py(k) = cy(k) - sy
                  pz(k) = cz(k) - sz
                  p2(k) = px(k)*px(k) &
                        + py(k)*py(k) &
                        + pz(k)*pz(k)
               end do
    
               do k = 1, block_length
                  g_00_00_00_000(k) = 0.0d0
               end do
               if (need_ao_order(1, 0)) then
                  do k = 1, block_length
                     g_00_00_00_100(k) = 0.0d0
                     g_00_00_00_010(k) = 0.0d0
                     g_00_00_00_001(k) = 0.0d0
                  end do
               end if
               do j = 1, nr_nonzero_contraction_coefs(ishell)
                  do k = 1, block_length
                     f1 = contraction_coef(ishell, j)*dexp(-primitive_exp(ishell, j)*p2(k))
                     g_00_00_00_000(k) = g_00_00_00_000(k) + f1
                  end do
                  if (need_ao_order(1, 0)) then
                     do k = 1, block_length
                        !fixme, repetitive
                        f1 = contraction_coef(ishell, j)*dexp(-primitive_exp(ishell, j)*p2(k))
                        f2 = -2.0d0*f1*primitive_exp(ishell, j)
                        g_00_00_00_100(k) = g_00_00_00_100(k) + px(k)*f2
                        g_00_00_00_010(k) = g_00_00_00_010(k) + py(k)*f2
                        g_00_00_00_001(k) = g_00_00_00_001(k) + pz(k)*f2
                     end do
                  end if
               end do
    
               off_power = power_offset(l_quantum_nr(ishell))
               if (need_ao_order(1, 0)) then
                  call get_gao1_block(block_length, &
                                      ishell,       &
                                      ao(1, i),     &
                                      ao(1, i + nr_ao*1),     &
                                      ao(1, i + nr_ao*2),     &
                                      ao(1, i + nr_ao*3),     &
                                      px,           &
                                      py,           &
                                      pz,           &
                                      p2,           &
                                      x_power(off_power), &
                                      y_power(off_power), &
                                      z_power(off_power))
               else
                  call get_gao_block(block_length,       &
                                     ishell,             &
                                     ao(1, i),           &
                                     px,                 &
                                     py,                 &
                                     pz,                 &
                                     x_power(off_power), &
                                     y_power(off_power), &
                                     z_power(off_power))
               end if
    
               i = i + spherical_deg(ishell)
    
            end if
         end do
      end do

   end subroutine

   subroutine get_gao_block(block_length, &
                            ishell,       &
                            ao,           &
                            px,           &
                            py,           &
                            pz, &
                     x_pwr, &
                     y_pwr, &
                     z_pwr)

      use ao_eval_block
!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      integer, intent(in)    :: ishell
      real(8), intent(inout) :: ao(block_length, *)
      real(8), intent(in)    :: px(*)
      real(8), intent(in)    :: py(*)
      real(8), intent(in)    :: pz(*)
      integer, intent(in)    :: x_pwr(*)
      integer, intent(in)    :: y_pwr(*)
      integer, intent(in)    :: z_pwr(*)
!     --------------------------------------------------------------------------
      integer                :: ic, jc, k
!     --------------------------------------------------------------------------

      if (l_quantum_nr(ishell) < 8) then
!        faster explicit code
         call get_cartesian_ao_g0(block_length,         &
                                  l_quantum_nr(ishell), &
                                  px,                   &
                                  py,                   &
                                  pz,                   &
                                  ao)
      else
!        slower general-order code
         do ic = 1, cartesian_deg(ishell)
            do k = 1, block_length
               ao(k, ic) = (px(k)**x_pwr(ic)) &
                          *(py(k)**y_pwr(ic)) &
                          *(pz(k)**z_pwr(ic))*g_00_00_00_000(k)
            end do
         end do
      end if

   end subroutine

   subroutine get_gao1_block(block_length, &
                            ishell,       &
                            ao,           &
                            aox,           &
                            aoy,           &
                            aoz,           &
                            px,           &
                            py,           &
                            pz, &
                            p2, &
                     x_pwr, &
                     y_pwr, &
                     z_pwr)

      use ao_eval_block
!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      integer, intent(in)    :: ishell
      real(8), intent(inout) :: ao(block_length, *)
      real(8), intent(inout) :: aox(block_length, *)
      real(8), intent(inout) :: aoy(block_length, *)
      real(8), intent(inout) :: aoz(block_length, *)
      real(8), intent(in)    :: px(*)
      real(8), intent(in)    :: py(*)
      real(8), intent(in)    :: pz(*)
      real(8), intent(in)    :: p2(*)
      integer, intent(in)    :: x_pwr(*)
      integer, intent(in)    :: y_pwr(*)
      integer, intent(in)    :: z_pwr(*)
!     --------------------------------------------------------------------------
      real(8)                :: p0
      real(8)                :: tx
      real(8)                :: ty
      real(8)                :: tz
      real(8)                :: ux
      real(8)                :: uy
      real(8)                :: uz
      real(8)                :: f2
      integer                :: ic, jc, k, l, m, n, i
      real(8)                :: p0t
!     --------------------------------------------------------------------------

      if (l_quantum_nr(ishell) < 8) then
!        faster explicit code
         call get_cartesian_ao_g1(block_length,         &
                                  l_quantum_nr(ishell), &
                                  px,                   &
                                  py,                   &
                                  pz,                   &
                                  ao,                   &
                                  aox,                  &
                                  aoy,                  &
                                  aoz)
      else
!        slower general-order code
         do ic = 1, cartesian_deg(ishell)
        
            l = x_pwr(ic)
            m = y_pwr(ic)
            n = z_pwr(ic)

            do k = 1, block_length
               tx = 0.0d0
               ty = 0.0d0
               tz = 0.0d0
               
               if (l > 0) tx = dfloat(l)*(px(k)**(l-1))*(py(k)**m)*(pz(k)**n)
               if (m > 0) ty = dfloat(m)*(px(k)**l)*(py(k)**(m-1))*(pz(k)**n)
               if (n > 0) tz = dfloat(n)*(px(k)**l)*(py(k)**m)*(pz(k)**(n-1))
               
               p0 = (px(k)**l)*(py(k)**m)*(pz(k)**n)
               
               ao (k, ic) = 0.0d0
               aox(k, ic) = 0.0d0
               aoy(k, ic) = 0.0d0
               aoz(k, ic) = 0.0d0
               do i = 1, nr_nonzero_contraction_coefs(ishell)
                  f2 = contraction_coef(ishell, i)*dexp(-primitive_exp(ishell, i)*p2(k))
               
                  p0t = -2.0d0*primitive_exp(ishell, i)*p0
                  ux  = px(k)*p0t + tx
                  uy  = py(k)*p0t + ty
                  uz  = pz(k)*p0t + tz
               
                  ao (k, ic) = ao (k, ic) + p0*f2
                  aox(k, ic) = aox(k, ic) + ux*f2
                  aoy(k, ic) = aoy(k, ic) + uy*f2
                  aoz(k, ic) = aoz(k, ic) + uz*f2
               end do
            end do
         end do
      end if

   end subroutine

  subroutine get_gao(ishell,  &
                     gao,     &
                     px,      &
                     py,      &
                     pz,      &
                     x_pwr, &
                     y_pwr, &
                     z_pwr)

   use ao_eval_point

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8), intent(out)   :: gao(*)
    real(8), intent(in)    :: px
    real(8), intent(in)    :: py
    real(8), intent(in)    :: pz
    integer, intent(in)    :: x_pwr(*)
    integer, intent(in)    :: y_pwr(*)
    integer, intent(in)    :: z_pwr(*)
!   ----------------------------------------------------------------------------
    integer                :: i, ic
!   ----------------------------------------------------------------------------

    if (l_quantum_nr(ishell) < 12) then
!      faster explicit code
       call get_cartesian_ao_g0(l_quantum_nr(ishell), &
                                px, py, pz,           &
                                gao)
    else
!      slower general-order code
       do ic = 1, cartesian_deg(ishell)
          gao(ic) = (px**x_pwr(ic)) &
                   *(py**y_pwr(ic)) &
                   *(pz**z_pwr(ic))*g_00_00_00_000
       end do
    end if

  end subroutine

  subroutine get_gao1(ishell, &
                      gao,    &
                      gaox,   &
                      gaoy,   &
                      gaoz,   &
                      px,     &
                      py,     &
                      pz,     &
                      p2,     &
                      x_pwr,  &
                      y_pwr,  &
                      z_pwr)

   use ao_eval_point

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8)                :: gao(*)
    real(8)                :: gaox(*)
    real(8)                :: gaoy(*)
    real(8)                :: gaoz(*)
    real(8), intent(in)    :: px, py, pz, p2
    integer, intent(in)    :: x_pwr(*)
    integer, intent(in)    :: y_pwr(*)
    integer, intent(in)    :: z_pwr(*)
!   ----------------------------------------------------------------------------
    integer                :: i, ic
    integer                :: l, m, n
    real(8)                :: f2
    real(8)                :: tx, ty, tz, ux, uy, uz
    real(8)                :: p0, p0t
!   ----------------------------------------------------------------------------

    if (l_quantum_nr(ishell) < 12) then
!      faster explicit code
       call get_cartesian_ao_g1(l_quantum_nr(ishell), &
                                px, py, pz,           &
                                gao,                  &
                                gaox,                 &
                                gaoy,                 &
                                gaoz)


    else
!      slower general-order code
       do ic = 1, cartesian_deg(ishell)
      
          l = x_pwr(ic)
          m = y_pwr(ic)
          n = z_pwr(ic)
      
          tx = 0.0d0
          ty = 0.0d0
          tz = 0.0d0
      
          if (l > 0) tx = dfloat(l)*(px**(l-1))*(py**m)*(pz**n)
          if (m > 0) ty = dfloat(m)*(px**l)*(py**(m-1))*(pz**n)
          if (n > 0) tz = dfloat(n)*(px**l)*(py**m)*(pz**(n-1))
      
          p0 = (px**l)*(py**m)*(pz**n)
      
          gao (ic) = 0.0d0
          gaox(ic) = 0.0d0
          gaoy(ic) = 0.0d0
          gaoz(ic) = 0.0d0
          do i = 1, nr_nonzero_contraction_coefs(ishell)
             f2 = contraction_coef(ishell, i)*dexp(-primitive_exp(ishell, i)*p2)
       
             p0t = -2.0d0*primitive_exp(ishell, i)*p0
             ux  = px*p0t + tx
             uy  = py*p0t + ty
             uz  = pz*p0t + tz
       
             gao (ic) = gao (ic) + p0*f2
             gaox(ic) = gaox(ic) + ux*f2
             gaoy(ic) = gaoy(ic) + uy*f2
             gaoz(ic) = gaoz(ic) + uz*f2
          end do
       end do
    end if

  end subroutine

  subroutine get_gao2(ishell, &
                      gaoxx,  &
                      gaoxy,  &
                      gaoxz,  &
                      gaoyy,  &
                      gaoyz,  &
                      gaozz,  &
                      px,     &
                      py,     &
                      pz,     &
                      p2,     &
                      x_pwr,  &
                      y_pwr,  &
                      z_pwr)

   use ao_eval_point

!   ----------------------------------------------------------------------------
    integer, intent(in)    :: ishell
    real(8)                :: gaoxx(*)
    real(8)                :: gaoxy(*)
    real(8)                :: gaoxz(*)
    real(8)                :: gaoyy(*)
    real(8)                :: gaoyz(*)
    real(8)                :: gaozz(*)
    real(8), intent(in)    :: px, py, pz, p2
    integer, intent(in)    :: x_pwr(*)
    integer, intent(in)    :: y_pwr(*)
    integer, intent(in)    :: z_pwr(*)
!   ----------------------------------------------------------------------------
    integer                :: i, ic
    integer                :: l, m, n
    real(8)                :: f2
    real(8)                :: tx, ty, tz, ux, uy, uz
    real(8)                :: p0, p0t
    real(8)                :: talph, tapax, tapay, tapaz
    real(8)                :: px0, pxd, pxm, pxp
    real(8)                :: py0, pyd, pym, pyp
    real(8)                :: pz0, pzd, pzm, pzp
    real(8)                :: gaxx
    real(8)                :: gaxy
    real(8)                :: gaxz
    real(8)                :: gayy
    real(8)                :: gayz
    real(8)                :: gazz
    real(8)                :: alpha, p000
!   ----------------------------------------------------------------------------

       do ic = 1, cartesian_deg(ishell)
      
          l = x_pwr(ic)
          m = y_pwr(ic)
          n = z_pwr(ic)
      
          pxd = 0.0d0
          pyd = 0.0d0
          pzd = 0.0d0
      
          if (l > 1) pxd = dfloat(l*(l-1))*(px**(l-2))
          if (m > 1) pyd = dfloat(m*(m-1))*(py**(m-2))
          if (n > 1) pzd = dfloat(n*(n-1))*(pz**(n-2))
      
          pxm = 0.0d0
          pym = 0.0d0
          pzm = 0.0d0
      
          if (l > 0) pxm = dfloat(l)*(px**(l-1))
          if (m > 0) pym = dfloat(m)*(py**(m-1))
          if (n > 0) pzm = dfloat(n)*(pz**(n-1))

          px0 = px**l
          py0 = py**m
          pz0 = pz**n
          p000 = px0*py0*pz0
      
          gaoxx(ic) = 0.0d0
          gaoxy(ic) = 0.0d0
          gaoxz(ic) = 0.0d0
          gaoyy(ic) = 0.0d0
          gaoyz(ic) = 0.0d0
          gaozz(ic) = 0.0d0
          do i = 1, nr_nonzero_contraction_coefs(ishell)

             f2 = contraction_coef(ishell, i)*dexp(-primitive_exp(ishell, i)*p2)
             alpha = primitive_exp(ishell, i)

             talph = -2.0d0*alpha
             tapax = talph*px
             tapay = talph*py
             tapaz = talph*pz

             pxp = tapax*px0
             pyp = tapay*py0
             pzp = tapaz*pz0

             gaxx = (tapax**2 + talph*dfloat(2*l + 1))*p000 + pxd*py0*pz0
             gayy = (tapay**2 + talph*dfloat(2*m + 1))*p000 + px0*pyd*pz0
             gazz = (tapaz**2 + talph*dfloat(2*n + 1))*p000 + px0*py0*pzd
             gaxy = tapax*tapay*p000 + (pxp*pym + pxm*pyp + pxm*pym)*pz0
             gaxz = tapax*tapaz*p000 + (pxp*pzm + pxm*pzp + pxm*pzm)*py0
             gayz = tapay*tapaz*p000 + (pyp*pzm + pym*pzp + pym*pzm)*px0
       
             gaoxx(ic) = gaoxx(ic) + gaxx*f2
             gaoxy(ic) = gaoxy(ic) + gaxy*f2
             gaoxz(ic) = gaoxz(ic) + gaxz*f2
             gaoyy(ic) = gaoyy(ic) + gayy*f2
             gaoyz(ic) = gaoyz(ic) + gayz*f2
             gaozz(ic) = gaozz(ic) + gazz*f2
          end do
       end do

  end subroutine

  subroutine get_ao_b1(px,     &
                       py,     &
                       pz,     &
                       ishell, &
                       ao,     &
                       ao_bx,  &
                       ao_by,  &
                       ao_bz)

!   ----------------------------------------------------------------------------
    real(8), intent(in)  :: px
    real(8), intent(in)  :: py
    real(8), intent(in)  :: pz
    integer, intent(in)  :: ishell
    real(8), intent(in)  :: ao(*)
    real(8), intent(out) :: ao_bx(*)
    real(8), intent(out) :: ao_by(*)
    real(8), intent(out) :: ao_bz(*)
!   ----------------------------------------------------------------------------
    integer              :: ic
    real(8)              :: f
!   ----------------------------------------------------------------------------

    do ic = 1, cartesian_deg(ishell)

      f = 0.5d0*ao(ic)

      ao_bx(ic) = f*px
      ao_by(ic) = f*py
      ao_bz(ic) = f*pz

    end do

  end subroutine


!> \brief point-wise transform ao to symmetrized ao
   subroutine ao_to_so_point(ao, &
                             buffer)

!     --------------------------------------------------------------------------
      real(8), intent(inout) :: ao(*)
      real(8), intent(inout) :: buffer(*)
!     --------------------------------------------------------------------------
      integer                :: i, n, off
!     --------------------------------------------------------------------------

      do i = 1, total_ao_vector_length
         buffer(i) = 0.0d0
      end do
      do n = 1, nr_ao_slices
         off = nr_ao*(n - 1)
         do i = 1, nr_so_ao
            buffer(off + i_so_ao(1, i)) = buffer(off + i_so_ao(1, i)) &
                                            + ao(off + i_so_ao(2, i))*f_so_ao(i)
         end do
      end do
      do i = 1, total_ao_vector_length
         ao(i) = buffer(i)
      end do

   end subroutine

!> \brief block transform ao to symmetrized ao
   subroutine ao_to_so_block(block_length, &
                             ao,           &
                             buffer)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: block_length
      real(8), intent(inout) :: ao(block_length, total_ao_vector_length)
      real(8), intent(inout) :: buffer(block_length, total_ao_vector_length)
!     --------------------------------------------------------------------------
      integer                :: i, k, n, off
!     --------------------------------------------------------------------------

      do i = 1, total_ao_vector_length
         do k = 1, block_length
            buffer(k, i) = 0.0d0
         end do
      end do
      do n = 1, nr_ao_slices
         off = nr_ao*(n - 1)
         do i = 1, nr_so_ao
            do k = 1, block_length
               buffer(k, off + i_so_ao(1, i)) = buffer(k, off + i_so_ao(1, i)) &
                                                  + ao(k, off + i_so_ao(2, i))*f_so_ao(i)
            end do
         end do
      end do
      do i = 1, total_ao_vector_length
         do k = 1, block_length
            ao(k, i) = buffer(k, i)
         end do
      end do

   end subroutine

!> \brief initialize ao_eval module
   subroutine ao_eval_init(max_ao_g_order, &
                           max_ao_m_order, &
                           use_screening)

!     --------------------------------------------------------------------------
      integer, intent(in) :: max_ao_g_order
      integer, intent(in) :: max_ao_m_order
      logical, intent(in) :: use_screening
!     --------------------------------------------------------------------------
      integer             :: g_order
      integer             :: m_order
!     --------------------------------------------------------------------------

      use_ao_screening = use_screening
      need_ao_order    = .false.
      nr_ao_slices     = 0

     do g_order = 0, max_ao_g_order
        do m_order = 0, max_ao_m_order
           need_ao_order(g_order, m_order) = .true.
           ao_start     (g_order, m_order) = nr_ao_slices*nr_ao + 1
           nr_ao_slices = nr_ao_slices + (3**g_order)*(3**m_order)
           if (g_order > 1) then
             nr_ao_slices = nr_ao_slices - 3*(g_order - 1)
           end if
           if (m_order > 1) then
             nr_ao_slices = nr_ao_slices - 3*(m_order - 1)
           end if
           ao_end       (g_order, m_order) = nr_ao_slices*nr_ao
        end do
     end do

!     total_ao_vector_length = nr_ao_slices*nr_ao
!     radovan: the above can lead to "plus one" or "plus few" out-of-bounds
!              error if spherical gaussians are used
!              it is safe to do instead
      total_ao_vector_length = nr_ao_slices*nr_ao_cartesian
!              although it consumes more memory than needed
!              the allocation could be refined

   end subroutine

!> \brief point-wise transform cartesian to spherical
   subroutine transform_c_to_s(ishell, &
                               ao)

!     --------------------------------------------------------------------------
      integer, intent(in)    :: ishell
      real(8), intent(inout) :: ao(*)
!     --------------------------------------------------------------------------
!> \fixme allocate cartesian to spherical buffer dynamically
      real(8)                :: buffer(1000)
      integer                :: i, k
!     --------------------------------------------------------------------------

      do i = 1, cartesian_deg(ishell)
         buffer(i) = ao(i)
      end do

      call mxma(cartesian_to_spherical(spherical_address(l_quantum_nr(ishell))), &
                1,                                                               &
                spherical_deg(ishell),                                           &
                buffer,                                                          &
                1,                                                               &
                cartesian_deg(ishell),                                           &
                ao,                                                              &
                1,                                                               &
                spherical_deg(ishell),                                           &
                spherical_deg(ishell),                                           &
                cartesian_deg(ishell),                                           &
                1)

   end subroutine

end module
