!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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

module visual_in_point

  use memory_allocator
  use electrostatic_potential
  use interface_ao
  use visual_cfg
  use visual_london
  use codata
  use ao_eval
  use density_eval
  use matrix_defop_old
  use visual_london
  use xc_london_c1

  implicit none

  public get_quantity_in_point

  private

#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"

contains

   subroutine get_quantity_in_point(nr_dmat, &
                                    D,       &
                                    D_0,     &
                                    px,      &
                                    py,      &
                                    pz,      &
                                    ao,      &
                                    buffer,  &
                                    result)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: nr_dmat
      type(matrix)         :: D(nr_dmat)
      type(matrix)         :: D_0
      real(8), intent(in)  :: px
      real(8), intent(in)  :: py
      real(8), intent(in)  :: pz
      real(8)              :: ao(*)
      real(8)              :: buffer(*)
      real(8), intent(out) :: result(0:3)
!     --------------------------------------------------------------------------
      integer              :: imat, i
      real(8)              :: v(3), v2(3)
      real(8)              :: matrix_33(3, 3)
      real(8)              :: temp_33(3, 3)
      real(8)              :: elf_chi, tau
      real(8)              :: density, esp(1)
      real(8)              :: density_gradient(3)
      real(8)              :: density_gradient_squared
      real(8)              :: b, f, distance_to_origin
      real(8)              :: rx, ry, rz, r(3)
!     --------------------------------------------------------------------------

      call get_ao(1, (/px/), (/py/), (/pz/), ao, buffer)

      result = 0.0d0
      do imat = 1, nr_dmat
         if (.not. visual_cfg_skip(imat)) then

!           gamma5 density: \psi^\dag \gamma^5 \psi
            if (visual_cfg_property(imat) == iq_gamma5) then
               call get_gamma5(density,       &
                               D(imat)%irep,  &
                               D(imat)%nrow,  &
                               D(imat)%elms,  &
                               buffer,        &
                               ao)
               result(0) = result(0) + density
            end if

!           electron density
            if (visual_cfg_property(imat) == iq_density) then
               call get_n(density,       &
                          D(imat)%irep,  &
                          D(imat)%nrow,  &
                          D(imat)%elms,  &
                          buffer,        &
                          ao)
               result(0) = result(0) + density
            end if

!           elf: electron localization function
            if (visual_cfg_property(imat) == iq_elf) then
               call get_kin_tau(tau,          &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            buffer,       &
                            ao)
               call get_gn(density,          &
                           density_gradient, &
                           D(imat)%irep,     &
                           D(imat)%nrow,     &
                           D(imat)%elms,     &
                           buffer,           &
                           ao)

               density_gradient_squared = density_gradient(1)*density_gradient(1) &
                                        + density_gradient(2)*density_gradient(2) &
                                        + density_gradient(3)*density_gradient(3)

               b = density_gradient_squared/(8.0d0*density)

!                  2.87123400018819d0 = 0.3d0*((3.0d0*pi*pi)**(2.0d0/3.0d0))
               f = 2.87123400018819d0*(density**(5.0d0/3.0d0))
               elf_chi = (tau - b)/f

               result(0) = result(0) + 1.0d0/(1.0d0 + elf_chi*elf_chi)
            end if

!           electrostatic potential - total
            if ((visual_cfg_property(imat) == iq_esp) .or. &
                (visual_cfg_property(imat) == iq_esprho)) then
               call get_esp(1, esp,          &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            (/px, py, pz/))
            end if

!           electrostatic potential - electronic part
            if ((visual_cfg_property(imat) == iq_espe) .or. &
                (visual_cfg_property(imat) == iq_esperho)) then
               call get_esp(1, esp,          &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            (/px, py, pz/), &
                            include_nuc_part=.false., &
                            include_el_part=.true.)
            end if

!           electrostatic potential - nuclear part
            if ((visual_cfg_property(imat) == iq_espn) .or. &
                (visual_cfg_property(imat) == iq_espnrho)) then
               call get_esp(1, esp,          &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            (/px, py, pz/), &
                            include_nuc_part=.true., &
                            include_el_part=.false.)
            end if

            if ((visual_cfg_property(imat) == iq_esprho)  .or. &
                (visual_cfg_property(imat) == iq_esperho) .or. &
                (visual_cfg_property(imat) == iq_espnrho)) then
!              multiply electrostatic potential with density
               call get_n(density,       &
                          D(imat)%irep,  &
                          D(imat)%nrow,  &
                          D(imat)%elms,  &
                          buffer,        &
                          ao)
               result(0) = result(0) + density*esp(1)
            else if ((visual_cfg_property(imat) == iq_esp)  .or. &
                     (visual_cfg_property(imat) == iq_espe) .or. &
                     (visual_cfg_property(imat) == iq_espn)) then
!              just the electrostatic potential
               result(0) = result(0) +         esp(1)
            end if

            if (visual_cfg_property(imat) == iq_kin) then
               call get_kin(density,      &
                            cval,         &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            buffer,       &
                            ao)
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_kin_ls) then
               call get_kin_ls(density,      &
                               cval,         &
                               D(imat)%irep, &
                               D(imat)%nrow, &
                               D(imat)%elms, &
                               buffer,       &
                               ao)
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_kin_sl) then
               call get_kin_sl(density,      &
                               cval,         &
                               D(imat)%irep, &
                               D(imat)%nrow, &
                               D(imat)%elms, &
                               buffer,       &
                               ao)
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_kin_tau) then
               call get_kin_tau(density,      &
                                D(imat)%irep, &
                                D(imat)%nrow, &
                                D(imat)%elms, &
                                buffer,       &
                                ao)
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_kin_lap) then
               call get_kin_lap(density,      &
                                D(imat)%irep, &
                                D(imat)%nrow, &
                                D(imat)%elms, &
                                buffer,       &
                                ao)
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_kin_nr) then
               call get_kin_tau(density,      &
                                D(imat)%irep, &
                                D(imat)%nrow, &
                                D(imat)%elms, &
                                buffer,       &
                                ao)
               result(0) = result(0) + 0.5d0*density
               call get_kin_lap(density,      &
                                D(imat)%irep, &
                                D(imat)%nrow, &
                                D(imat)%elms, &
                                buffer,       &
                                ao)
               result(0) = result(0) + 0.5d0*density
            end if

!           current density: -ec \psi^\dag \alpha \psi
            if (visual_cfg_property(imat) == iq_j) then
               if (levyle) then
                  call get_gn(density,      &
                              v,            &
                              D(imat)%irep, &
                              D(imat)%nrow, &
                              D(imat)%elms, &
                              buffer,       &
                              ao)
                  v = 0.5d0*v
               else
                  call get_j(v,            &
                             cval,         &
                             D(imat)%irep, &
                             D(imat)%nrow, &
                             D(imat)%elms, &
                             buffer,       &
                             ao)
               end if
               if (visual_cfg_london .and. .not. visual_cfg_london_skip_direct) then
                  if (levyle) then
                     call get_jblao_ll(visual_cfg_london_component,D_0%elms,                      &
                                       ntbas(0),                      &
                                       ao, &
                                       (/px, py, pz/),                &
                                       v2)
                  else
                     call get_jblao(visual_cfg_london_component,cval,D_0%elms,                      &
                                    ntbas(0),                      &
                                    ao, &
                                    (/px, py, pz/),                &
                                    v2)
                  end if
                  v = v + v2
               end if
               result(1:3) = result(1:3) + v(1:3)
            end if

!           spin density: \psi^\dag \Sigma \psi
            if (visual_cfg_property(imat) == iq_s) then
               call get_s(v,            &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)
               result(1:3) = result(1:3) + v(1:3)
            end if

!           divergence or curl of spin density
            if ((visual_cfg_property(imat) == iq_divs) .or. &
                (visual_cfg_property(imat) == iq_rots)) then
               call get_gs(v,               &
                           matrix_33,       &
                           D(imat)%irep,    &
                           D(imat)%nrow,    &
                           D(imat)%elms,    &
                           buffer,          &
                           ao)
               if (visual_cfg_property(imat) == iq_divs) then
                  result(0) = result(0) + matrix_33(1, 1) + matrix_33(2, 2) + matrix_33(3, 3)
               end if
               if (visual_cfg_property(imat) == iq_rots) then
                  result(1) = result(1) - matrix_33(2, 3) + matrix_33(3, 2)
                  result(2) = result(2) - matrix_33(3, 1) + matrix_33(1, 3)
                  result(3) = result(3) - matrix_33(1, 2) + matrix_33(2, 1)
               end if
            end if

!           divergence or curl of current density
            if ((visual_cfg_property(imat) == iq_divj) .or. &
                (visual_cfg_property(imat) == iq_rotj)) then

               if (levyle) then
                  print *, 'error: Levy-Leblond divj and rotj not implemented'
                  stop
               end if

               call get_gj(v,            &
                           matrix_33,    &
                           cval,         &
                           D(imat)%irep, &
                           D(imat)%nrow, &
                           D(imat)%elms, &
                           buffer,       &
                           ao)

               if (visual_cfg_london .and. .not. visual_cfg_london_skip_direct) then
                  if (levyle) then
                     print *, 'error: Levy-Leblond divj and rotj not implemented'
                     stop
                  else
                     call get_gjblao(visual_cfg_london_component, &
                                     cval,                        &
                                     D_0%elms,                    &
                                     ntbas(0),                    &
                                     ao,                          &
                                     (/px, py, pz/),              &
                                     temp_33)
                  end if
                  matrix_33 = matrix_33 + temp_33
               end if

               if (visual_cfg_property(imat) == iq_divj) then
                  result(0) = result(0) + matrix_33(1, 1) + matrix_33(2, 2) + matrix_33(3, 3)
               end if
               if (visual_cfg_property(imat) == iq_rotj) then
                  result(1) = result(1) - matrix_33(2, 3) + matrix_33(3, 2)
                  result(2) = result(2) - matrix_33(3, 1) + matrix_33(1, 3)
                  result(3) = result(3) - matrix_33(1, 2) + matrix_33(2, 1)
               end if
            end if

!           nonrel diamagnetic current density (hardcoded for B along z)
            if (visual_cfg_property(imat) == iq_jdia) then
               call get_n(density,      &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)

               rx = px - visual_cfg_gauge_origin(1)
               ry = py - visual_cfg_gauge_origin(2)
               rz = pz - visual_cfg_gauge_origin(3)

               result(1) = result(1) + 0.5d0*ry*density
               result(2) = result(2) - 0.5d0*rx*density
               result(3) = result(3) - 0.0d0
            end if

!           responding operator is electric dipole
            if ((visual_cfg_property(imat) == iq_edipx) .or. &
                (visual_cfg_property(imat) == iq_edipy) .or. &
                (visual_cfg_property(imat) == iq_edipz)) then

               call get_n(density,      &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)
               if (visual_cfg_property(imat) == iq_edipx) then
                  result(0) = result(0) + px*density
               end if
               if (visual_cfg_property(imat) == iq_edipy) then
                  result(0) = result(0) + py*density
               end if
               if (visual_cfg_property(imat) == iq_edipz) then
                  result(0) = result(0) + pz*density
               end if
            end if

!           responding operator is magnetic dipole
            if ((visual_cfg_property(imat) == iq_bdipx) .or. &
                (visual_cfg_property(imat) == iq_bdipy) .or. &
                (visual_cfg_property(imat) == iq_bdipz)) then

               call get_j(v,            &
                          cval,         &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)

               rx = px - visual_cfg_gauge_origin(1)
               ry = py - visual_cfg_gauge_origin(2)
               rz = pz - visual_cfg_gauge_origin(3)

               f = -0.5d0

               if (visual_cfg_property(imat) == iq_bdipx) then
                  result(0) = result(0) + f*(ry*v(3) - rz*v(2))
               end if
               if (visual_cfg_property(imat) == iq_bdipy) then
                  result(0) = result(0) + f*(rz*v(1) - rx*v(3))
               end if
               if (visual_cfg_property(imat) == iq_bdipz) then
                  result(0) = result(0) + f*(rx*v(2) - ry*v(1))
               end if
            end if

!           responding operator is nuclear magnetic dipole
            if ((visual_cfg_property(imat) == iq_ndipx) .or. &
                (visual_cfg_property(imat) == iq_ndipy) .or. &
                (visual_cfg_property(imat) == iq_ndipz)) then

               if (levyle) then
                  call get_gn(density,      &
                              v,            &
                              D(imat)%irep, &
                              D(imat)%nrow, &
                              D(imat)%elms, &
                              buffer,       &
                              ao)
                  v = 0.5d0*v
               else
                  call get_j(v,            &
                             cval,         &
                             D(imat)%irep, &
                             D(imat)%nrow, &
                             D(imat)%elms, &
                             buffer,       &
                             ao)
               end if
               if (visual_cfg_london .and. .not. visual_cfg_london_skip_direct) then
                  if (levyle) then
                     call get_jblao_ll(visual_cfg_london_component,D_0%elms,                      &
                                       ntbas(0),                      &
                                       ao, &
                                       (/px, py, pz/),                &
                                       v2)
                  else
                     call get_jblao(visual_cfg_london_component,cval,D_0%elms,                      &
                                    ntbas(0),                      &
                                    ao, &
                                    (/px, py, pz/),                &
                                    v2)
                  end if
                  v = v + v2
               end if

               rx = px - cord(1, visual_cfg_ref_nucleus(imat))
               ry = py - cord(2, visual_cfg_ref_nucleus(imat))
               rz = pz - cord(3, visual_cfg_ref_nucleus(imat))

               f = -1.0d0/((rx*rx + ry*ry + rz*rz)**1.5d0)

               if (visual_cfg_property(imat) == iq_ndipx) then
                  result(0) = result(0) + f*(ry*v(3) - rz*v(2))
               end if
               if (visual_cfg_property(imat) == iq_ndipy) then
                  result(0) = result(0) + f*(rz*v(1) - rx*v(3))
               end if
               if (visual_cfg_property(imat) == iq_ndipz) then
                  result(0) = result(0) + f*(rx*v(2) - ry*v(1))
               end if
            end if

!           responding operator is nuclear magnetic dipole
            if ((visual_cfg_property(imat) == iq_ndipxdia) .or. &
                (visual_cfg_property(imat) == iq_ndipydia) .or. &
                (visual_cfg_property(imat) == iq_ndipzdia)) then

               call j_dia(matrix_33,               &
                          visual_cfg_london,       &
                          visual_cfg_gauge_origin, &
                          D(imat)%nrow,            &
                          D(imat)%elms,            &
                          ao,                      &
                          (/px, py, pz/))

               if (visual_cfg_property(imat) == iq_ndipxdia) i = 1
               if (visual_cfg_property(imat) == iq_ndipydia) i = 2
               if (visual_cfg_property(imat) == iq_ndipzdia) i = 3

               v(1) = matrix_33(1, i)
               v(2) = matrix_33(2, i)
               v(3) = matrix_33(3, i)

               r(1) = px - cord(1, visual_cfg_ref_nucleus(imat))
               r(2) = py - cord(2, visual_cfg_ref_nucleus(imat))
               r(3) = pz - cord(3, visual_cfg_ref_nucleus(imat))

               v2(1) = r(2)*v(3) - r(3)*v(2)
               v2(2) = r(3)*v(1) - r(1)*v(3)
               v2(3) = r(1)*v(2) - r(2)*v(1)

               f = -1.0d0/((r(1)*r(1) + r(2)*r(2) + r(3)*r(3))**1.5d0)
               result(0) = result(0) + f*v2(i)

            end if

         end if
      end do

      f = visual_cfg_scale
      if (visual_cfg_radial) then
         distance_to_origin = dsqrt(px*px + py*py + pz*pz)
         f = f*(distance_to_origin**2.0d0)*4.0d0*3.14159265358979323846d0
      end if
      result = f*result

   end subroutine

end module
