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

module visual_in_point

  use memory_allocator
  use vc_grid
  use visual_cfg
  use visual_london
  use codata
  use ao_eval
  use density_eval
  use matrix_defop
#ifdef MOD_UNRELEASED
  use visual_london
  use xc_london_c1
#endif

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
                                    q)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: nr_dmat
      type(matrix)         :: D(nr_dmat)
      type(matrix)         :: D_0
      real(8), intent(in)  :: px
      real(8), intent(in)  :: py
      real(8), intent(in)  :: pz
      real(8)              :: ao(*)
      real(8)              :: buffer(*)
      real(8), intent(out) :: q(0:3, nr_dmat)
!     --------------------------------------------------------------------------
      integer              :: imat, i
      real(8)              :: v(3), v2(3)
      real(8)              :: matrix_33(3, 3)
      real(8)              :: elf_chi, tau
      real(8)              :: density, esp
      real(8)              :: density_gradient(3)
      real(8)              :: density_gradient_squared
      real(8)              :: b, f, distance_to_origin
      real(8)              :: rx, ry, rz, r(3)
!     --------------------------------------------------------------------------

      call get_ao(px, py, pz, ao, buffer)

      q = 0.0d0
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
               q(0, imat) = density
            end if
            
!           electron density
            if (visual_cfg_property(imat) == iq_density) then
               call get_n(density,       &
                          D(imat)%irep,  &
                          D(imat)%nrow,  &
                          D(imat)%elms,  &
                          buffer,        &
                          ao)
               q(0, imat) = density
            end if
            
!           elf: electron localization function
            if (visual_cfg_property(imat) == iq_elf) then
               call get_tau(tau,          &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            buffer,       &
                            ao(ao_start(1, 0):ao_end(1, 0)))
               call get_gn(density,          &
                           density_gradient, &
                           D(imat)%irep,     &
                           D(imat)%nrow,     &
                           D(imat)%elms,     &
                           buffer,           &
                           ao,               &
                           ao(ao_start(1, 0):ao_end(1, 0)))
              
               density_gradient_squared = density_gradient(1)*density_gradient(1) &
                                        + density_gradient(2)*density_gradient(2) &
                                        + density_gradient(3)*density_gradient(3)
              
               b = density_gradient_squared/(8.0d0*density)
              
!                  2.87123400018819d0 = 0.3d0*((3.0d0*pi*pi)**(2.0d0/3.0d0))
               f = 2.87123400018819d0*(density**(5.0d0/3.0d0))
               elf_chi = (tau - b)/f
              
               q(0, imat) = 1.0d0/(1.0d0 + elf_chi*elf_chi)
            end if
            
!           total electrostatic potential
            if (visual_cfg_property(imat) == iq_esp) then
               call get_esp(esp,3,       &
                          D(imat)%irep,  &
                          D(imat)%nrow,  &
                          D(imat)%elms,  &
                          px,py,pz)
                         
               q(0, imat) = esp
            end if

!           density of Coulomb 2-electron contribution
            if (visual_cfg_property(imat) == iq_esprho) then
               call get_esp(esp,1,       &
                          D(imat)%irep,  &
                          D(imat)%nrow,  &
                          D(imat)%elms,  &
                          px,py,pz)
               call get_n(density,       &
                          D(imat)%irep,  &
                          D(imat)%nrow,  &
                          D(imat)%elms,  &
                          buffer,        &
                          ao)
               q(0, imat) = density*esp*0.5d0

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
                              ao,           &
                              ao(ao_start(1, 0):ao_end(1, 0)))
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
#ifdef MOD_UNRELEASED
               if (visual_cfg_london .and. .not. visual_cfg_london_skip_direct) then
                  if (levyle) then
                     call get_jblao_ll(D_0%elms,                      &
                                       ntbas(0),                      &
                                       ao(ao_start(0,0):ao_end(0,0)), &
                                       ao(ao_start(0,1):ao_end(0,1)), &
                                       ao(ao_start(1,0):ao_end(1,0)), &
                                       ao(ao_start(1,1):ao_end(1,1)), &
                                       (/px, py, pz/),                &
                                       v2)
                  else
                     call get_jblao(D_0%elms,                      &
                                    ntbas(0),                      &
                                    ao(ao_start(0,0):ao_end(0,0)), &
                                    ao(ao_start(0,1):ao_end(0,1)), &
                                    (/px, py, pz/),                &
                                    v2)
                  end if
                  v = v + v2
               end if
#endif
               q(1:3, imat) = v(1:3)
            end if
            
!           spin density: \psi^\dag \Sigma \psi
            if (visual_cfg_property(imat) == iq_s) then
               call get_s(v,            &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)
               q(1:3, imat) = v(1:3)
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
                           ao,              &
                           ao(ao_start(1, 0):ao_end(1, 0)))
               if (visual_cfg_property(imat) == iq_divs) then
                  q(0, imat) = matrix_33(1, 1) + matrix_33(2, 2) + matrix_33(3, 3)
               end if
               if (visual_cfg_property(imat) == iq_rots) then
                  q(1, imat) = matrix_33(3, 2) - matrix_33(2, 3)
                  q(2, imat) = matrix_33(1, 3) - matrix_33(3, 1)
                  q(3, imat) = matrix_33(2, 1) - matrix_33(1, 2)
               end if
            end if
            
!           divergence or curl of current density
            if ((visual_cfg_property(imat) == iq_divj) .or. &
                (visual_cfg_property(imat) == iq_rotj)) then
               call get_gj(v,            &
                           matrix_33,    &
                           cval,         &
                           D(imat)%irep, &
                           D(imat)%nrow, &
                           D(imat)%elms, &
                           buffer,       &
                           ao,           &
                           ao(ao_start(1, 0):ao_end(1, 0)))
               if (visual_cfg_property(imat) == iq_divj) then
                  q(0, imat) = matrix_33(1, 1) + matrix_33(2, 2) + matrix_33(3, 3)
               end if
               if (visual_cfg_property(imat) == iq_rotj) then
                  q(1, imat) = matrix_33(3, 2) - matrix_33(2, 3)
                  q(2, imat) = matrix_33(1, 3) - matrix_33(3, 1)
                  q(3, imat) = matrix_33(2, 1) - matrix_33(1, 2)
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
            
               q(1, imat) = +0.5d0*ry*density
               q(2, imat) = -0.5d0*rx*density
               q(3, imat) = -0.0d0           
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
                  q(0, imat) = px*density
               end if
               if (visual_cfg_property(imat) == iq_edipy) then
                  q(0, imat) = py*density
               end if
               if (visual_cfg_property(imat) == iq_edipz) then
                  q(0, imat) = pz*density
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
                  q(0, imat) = f*(ry*v(3) - rz*v(2))
               end if
               if (visual_cfg_property(imat) == iq_bdipy) then
                  q(0, imat) = f*(rz*v(1) - rx*v(3))
               end if
               if (visual_cfg_property(imat) == iq_bdipz) then
                  q(0, imat) = f*(rx*v(2) - ry*v(1))
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
                              ao,           &
                              ao(ao_start(1, 0):ao_end(1, 0)))
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
#ifdef MOD_UNRELEASED
               if (visual_cfg_london .and. .not. visual_cfg_london_skip_direct) then
                  if (levyle) then
                     call get_jblao_ll(D_0%elms,                      &
                                       ntbas(0),                      &
                                       ao(ao_start(0,0):ao_end(0,0)), &
                                       ao(ao_start(0,1):ao_end(0,1)), &
                                       ao(ao_start(1,0):ao_end(1,0)), &
                                       ao(ao_start(1,1):ao_end(1,1)), &
                                       (/px, py, pz/),                &
                                       v2)
                  else
                     call get_jblao(D_0%elms,                      &
                                    ntbas(0),                      &
                                    ao(ao_start(0,0):ao_end(0,0)), &
                                    ao(ao_start(0,1):ao_end(0,1)), &
                                    (/px, py, pz/),                &
                                    v2)
                  end if
                  v = v + v2
               end if
#endif
            
               rx = px - cord(1, visual_cfg_ref_nucleus(imat))
               ry = py - cord(2, visual_cfg_ref_nucleus(imat))
               rz = pz - cord(3, visual_cfg_ref_nucleus(imat))
               
               f = -1.0d0/((rx*rx + ry*ry + rz*rz)**1.5d0)
               
               if (visual_cfg_property(imat) == iq_ndipx) then
                  q(0, imat) = f*(ry*v(3) - rz*v(2))
               end if
               if (visual_cfg_property(imat) == iq_ndipy) then
                  q(0, imat) = f*(rz*v(1) - rx*v(3))
               end if
               if (visual_cfg_property(imat) == iq_ndipz) then
                  q(0, imat) = f*(rx*v(2) - ry*v(1))
               end if
            end if
            
!           responding operator is nuclear magnetic dipole
            if ((visual_cfg_property(imat) == iq_ndipxdia) .or. &
                (visual_cfg_property(imat) == iq_ndipydia) .or. &
                (visual_cfg_property(imat) == iq_ndipzdia)) then

#ifdef MOD_UNRELEASED
               call j_dia(matrix_33,               &
                          visual_cfg_london,       &
                          visual_cfg_gauge_origin, &
                          D(imat)%nrow,            &
                          D(imat)%elms,            &
                          ao,                      &
                          (/px, py, pz/))
#endif

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
               q(0, imat) = f*v2(i)

            end if

         end if
      end do

      f = visual_cfg_scale
      if (visual_cfg_radial) then
         distance_to_origin = dsqrt(px*px + py*py + pz*pz)
         f = f*(distance_to_origin**2.0d0)
      end if
      q = f*q

   end subroutine

end module
