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
  use visual_cfg
  use visual_london
  use codata
  use ao_eval
  use density_eval
  use basis_info
  use matrix_defop
#ifdef MOD_UNRELEASED
  use visual_london
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
      integer              :: imat
      real(8)              :: v(3), v2(3)
      real(8)              :: matrix_33(3, 3)
      real(8)              :: elf_chi, tau
      real(8)              :: density
      real(8)              :: density_gradient(3)
      real(8)              :: density_gradient_squared
      real(8)              :: b, f, distance_to_origin
      real(8)              :: rx, ry, rz
!     --------------------------------------------------------------------------

      call get_ao(px, py, pz, ao, buffer)

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
              
               result(0) = result(0) + 1.0d0/(1.0d0 + elf_chi*elf_chi)
            end if
            
!           current density: -ec \psi^\dag \alpha \psi
            if (visual_cfg_property(imat) == iq_j) then
               call get_j(v,            &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)
               if (visual_cfg_london) then
                  call get_jblao(imat,       &
                                 D_0,        &
                                 ao,         &
                                 buffer,     &
                                 px, py, pz, &
                                 v2)
                  v = v - v2
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
                           ao,              &
                           ao(ao_start(1, 0):ao_end(1, 0)))
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
               call get_gj(v,            &
                           matrix_33,    &
                           D(imat)%irep, &
                           D(imat)%nrow, &
                           D(imat)%elms, &
                           buffer,       &
                           ao,           &
                           ao(ao_start(1, 0):ao_end(1, 0)))
               if (visual_cfg_property(imat) == iq_divj) then
                  result(0) = result(0) + matrix_33(1, 1) + matrix_33(2, 2) + matrix_33(3, 3)
               end if
               if (visual_cfg_property(imat) == iq_rotj) then
                  result(1) = result(1) - matrix_33(2, 3) + matrix_33(3, 2)
                  result(2) = result(2) - matrix_33(3, 1) + matrix_33(1, 3)
                  result(3) = result(3) - matrix_33(1, 2) + matrix_33(2, 1)
               end if
            end if
            
!           nonrel paramagnetic current density (orbital part)
            if (visual_cfg_property(imat) == iq_jpara) then
               call get_gn(density,      &
                           v,            &
                           D(imat)%irep, &
                           D(imat)%nrow, &
                           D(imat)%elms, &
                           buffer,       &
                           ao,           &
                           ao(ao_start(1, 0):ao_end(1, 0)))
               result(1:3) = result(1:3) + 0.5d0*v(1:3)
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
            
               result(1) = result(1) - 0.5d0*ry*density
               result(2) = result(2) + 0.5d0*rx*density
               result(3) = result(3) + 0.0d0
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
            
               call get_j(v,            &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)
               if (visual_cfg_london) then
                  call get_jblao(imat,       &
                                 D_0,        &
                                 ao,         &
                                 buffer,     &
                                 px, py, pz, &
                                 v2)
                  v = v - v2
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

         end if
      end do

      f = visual_cfg_scale
      if (visual_cfg_radial) then
         distance_to_origin = dsqrt(px*px + py*py + pz*pz)
         f = f*(distance_to_origin**2.0d0)
      end if
      result = f*result

   end subroutine

end module



!radovan: rotjpara note for myself please do not delete
!     k = 1
!     l = 1
!     m = 2
!     n = 2
!     a = d1
!     b = d1
!     itransform = 0

!     do ixyz = 1,3

!       jxyz = ixyz + 1
!       kxyz = jxyz + 1

!       if(jxyz .gt. 3) jxyz = jxyz - 3
!       if(kxyz .gt. 3) kxyz = kxyz - 3

!       call get_density_sub(lget,.true.,
!    &                       vector(ixyz),
!    &                       a,b,k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao1(1+joff*(jxyz-1)),
!    &                       gao1(1+joff*(kxyz-1)),
!    &                       dmat,buf)
!       call get_density_sub(lget,.false.,
!    &                       vector(ixyz),
!    &                       -a,-b,k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao1(1+joff*(kxyz-1)),
!    &                       gao1(1+joff*(jxyz-1)),
!    &                       dmat,buf)

!     enddo

!radovan: rotvpara note for myself please do not delete
!     k = 1
!     l = 1
!     m = 2
!     n = 2
!     a = d1/density
!     b = d1/density
!     itransform = 0

!     do ixyz = 1,3

!       jxyz = ixyz + 1
!       kxyz = jxyz + 1

!       if(jxyz .gt. 3) jxyz = jxyz - 3
!       if(kxyz .gt. 3) kxyz = kxyz - 3

!       call get_density_sub(lget,.true.,
!    &                       vector(ixyz),
!    &                       a,
!    &                       b,
!    &                       k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao1(1+joff*(jxyz-1)),
!    &                       gao1(1+joff*(kxyz-1)),
!    &                       dmat,buf)
!       call get_density_sub(lget,.false.,
!    &                       vector(ixyz),
!    &                       -a,
!    &                       -b,
!    &                       k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao1(1+joff*(kxyz-1)),
!    &                       gao1(1+joff*(jxyz-1)),
!    &                       dmat,buf)

!       c = -density_gradient(jxyz)/(density*d2)
!       call get_density_sub(lget,.false.,
!    &                       vector(ixyz),
!    &                       a*c,
!    &                       b*c,
!    &                       k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao,
!    &                       gao1(1+joff*(kxyz-1)),
!    &                       dmat,buf)
!       call get_density_sub(lget,.false.,
!    &                       vector(ixyz),
!    &                       -a*c,
!    &                       -b*c,
!    &                       k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao1(1+joff*(kxyz-1)),
!    &                       gao,
!    &                       dmat,buf)

!       c = density_gradient(kxyz)/(density*d2)
!       call get_density_sub(lget,.false.,
!    &                       vector(ixyz),
!    &                       a*c,
!    &                       b*c,
!    &                       k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao,
!    &                       gao1(1+joff*(jxyz-1)),
!    &                       dmat,buf)
!       call get_density_sub(lget,.false.,
!    &                       vector(ixyz),
!    &                       -a*c,
!    &                       -b*c,
!    &                       k,l,m,n,
!    &                       1,itransform,irepop,
!    &                       gao1(1+joff*(jxyz-1)),
!    &                       gao,
!    &                       dmat,buf)

!     enddo
