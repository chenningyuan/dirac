!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful, 
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

module visual_in_point

  use memory_allocator
  use electrostatic_potential
  use interface_ao
  use interface_mo
  use visual_cfg
  use visual_london
  use dirac_ao_eval
  use density_eval
  use matrix_defop_old
  use visual_london
  use xc_london_c1

  use visual_datatypes

  implicit none

  public get_quantity_in_point
  public get_quantity_in_point_refactored
  public get_esp_functions_in_points

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
#include "ccom.h"

contains

! ------------------------------------------------------------------------------
  subroutine get_quantity_in_point_refactored(dmat_0,&
                                              ao,buffer,&
                                              grid_function,&
                                              grid_function_values,&
                                              px,py,pz,&
                                              ipoint)

    type(dmat_t), intent(in):: dmat_0
    real(8), intent(in):: ao(*)
    real(8), intent(in):: buffer(*)
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)
    real(8), intent(in):: px, py, pz
    integer, intent(in), optional  :: ipoint !thiswill not be optional

    integer:: mat_dim, mat_nz, mat_irep
    integer:: imat,i,j,itemp(3)
    real(8):: a, b, c, f, v1, v2(3), v3(3,3)
    real(8):: elf_chi, tau
    real(8):: eigenvalues(3),eigenvectors(3, 3)
    real(8):: temp(3), temp2(9)
  integer:: ioff

     !px = grid%points(1, ipoint)
     !py = grid%points(2, ipoint)
     !pz = grid%points(3, ipoint)
     !todo: call grid%points%get_cartesian_coords()

    mat_dim = dmat_0%mat_dim
    mat_nz  = dmat_0%mat_nz
    mat_irep= 0

    select case (grid_function%name)

      case ('gamma5')
!       gamma5 density: \psi^\dag \gamma^5 \psi
        v1 = 0.0d0
        call get_gamma5(v1,       &
                        mat_irep, &
                        mat_dim,  &
                        dmat_0%elements,   &
                        buffer,   &
                        ao)
        grid_function_values(1:1,ipoint:ipoint) = v1

      case ('density')
!       electron density
        v1 = 0.0d0
        call get_n(v1,       &
                   mat_irep, &
                   mat_dim,  &
                   dmat_0%elements,   &
                   buffer,   &
                   ao)
        grid_function_values(1:1,ipoint:ipoint) = v1

      case ('density_gradient')
!       electron density gradient
        v1 = 0.0d0
        v2 = 0.0d0
        call get_gn(v1,      &
                    v2,      &
                    mat_irep,&
                    mat_dim, &
                    dmat_0%elements,   &
                    buffer,  &
                    ao)               
        grid_function_values(1:1,ipoint:ipoint) = v2(1)
        grid_function_values(2:2,ipoint:ipoint) = v2(2)
        grid_function_values(3:3,ipoint:ipoint) = v2(3)

      case ('density_hessian')
!       gosia fixme:
!       density Hessian is implemented for hermitian density matrices (such as the unperturbed density matrix)
!       todo: add a warning
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        call get_hn(v1,       &
                    v2,       &
                    v3,       &
                    mat_irep, &
                    mat_dim,  &
                    dmat_0%elements,   &
                    buffer,   &
                    ao)               
        grid_function_values(1:1,ipoint:ipoint) = v3(1, 1)   ! xx
        grid_function_values(2:2,ipoint:ipoint) = v3(1, 2)   ! xy
        grid_function_values(3:3,ipoint:ipoint) = v3(1, 3)   ! xz
        grid_function_values(4:4,ipoint:ipoint) = v3(2, 1)   ! yx
        grid_function_values(5:5,ipoint:ipoint) = v3(2, 2)   ! yy
        grid_function_values(6:6,ipoint:ipoint) = v3(2, 3)   ! yz
        grid_function_values(7:7,ipoint:ipoint) = v3(3, 1)   ! zx
        grid_function_values(8:8,ipoint:ipoint) = v3(3, 2)   ! zy
        grid_function_values(9:9,ipoint:ipoint) = v3(3, 3)   ! zz

      case ('rdg')
!       reduced density gradient
        v1 = 0.0d0
        v2 = 0.0d0
        call get_gn(v1,       &
                    v2,       &
                    mat_irep, &
                    mat_dim,  &
                    dmat_0%elements,   &
                    buffer,   &
                    ao)

        b = sqrt(v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3))
        c = v1**(4.0d0/3.0d0)
!       0.161620459674d0 = 1/(2.0d0*((3.0d0*pi*pi)**(1.0d0/3.0d0))
        grid_function_values(1:1,ipoint:ipoint) = 0.161620459674d0*(b/c)

      case ('ed_sign_lambda2')
!       sign(lambda2), where lambda2=the second eigenvalue of the Hessian
!       of the electron density
!
!       gosia fixme:
!       density Hessian is implemented for hermitian density matrices (such as the unperturbed density matrix)
!       todo: add a warning
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        call get_hn(v1,       &
                    v2,       &
                    v3,       &
                    mat_irep, &
                    mat_dim,  &
                    dmat_0%elements,   &
                    buffer,   &
                    ao)               
        eigenvalues  = 0.0d0
        eigenvectors = 1.0d0
        temp         = 0.0d0
        itemp        = 0
        temp2 = reshape(v3, (/9/))

        call jaco(temp2, eigenvectors, 3, 3, 3, temp, itemp)

        temp(1)        = temp2(1)
        temp(2)        = temp2(3)
        temp(3)        = temp2(6)
        eigenvalues(1) = temp2(1)
        eigenvalues(2) = temp2(3)
        eigenvalues(3) = temp2(6)
        do i = 1, 3
           do j = 1, 3
              if (j .ne. i) then
                 if (eigenvalues(i) .lt. eigenvalues(j)) then
                    eigenvalues(i) = temp(j)
                 end if  
              end if  
           end do
        end do

        if (eigenvalues(2) .gt. 0.0d0) then
           grid_function_values(1:1,ipoint:ipoint) = v1
        else if (eigenvalues(2) .lt. 0.0d0) then
           grid_function_values(1:1,ipoint:ipoint) = -v1
        end if


      case ('ed_laplacian')
!
!       gosia fixme:
!       density Hessian is implemented for hermitian density matrices (such as the unperturbed density matrix)
!       todo: add a warning
        v1 = 0.0d0
        v2 = 0.0d0
        v3 = 0.0d0
        call get_hn(v1,       &
                    v2,       &
                    v3,       &
                    mat_irep, &
                    mat_dim,  &
                    dmat_0%elements,   &
                    buffer,   &
                    ao)               
        eigenvectors = 1.0d0
        temp         = 0.0d0
        itemp        = 0
        temp2 = reshape(v3, (/9/))

        call jaco(temp2, eigenvectors, 3, 3, 3, temp, itemp)

        grid_function_values(1:1,ipoint:ipoint) = temp2(1) + temp2(3) + temp2(6)

      case ('elf')
!       elf: electron localization function
        tau = 0.0d0
        v1 = 0.0d0
        v2 = 0.0d0
        call get_kin_tau(tau,      &
                         mat_irep, &
                         mat_dim,  &
                         dmat_0%elements,&
                         buffer,   &
                         ao)
        call get_gn(v1,&
                    v2,&
                    mat_irep, &
                    mat_dim,  &
                    dmat_0%elements,&
                    buffer,  &
                    ao)

        a = v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3)
        b = a/(8.0d0*v1)

!           2.87123400018819d0 = 0.3d0*((3.0d0*pi*pi)**(2.0d0/3.0d0))
        f = 2.87123400018819d0*(v1**(5.0d0/3.0d0))
        elf_chi = (tau-b)/f

        grid_function_values(1:1,ipoint:ipoint) = 1.0d0/(1.0d0+elf_chi*elf_chi)

      case ('kin')
!       kin: 
        v1 = 0.0d0
        call get_kin(v1,       &
                     visual_cfg_cval,&
                     mat_irep, &
                     mat_dim,  &
                     dmat_0%elements,&
                     buffer,   &
                     ao)

        grid_function_values(1:1,ipoint:ipoint) = v1


      case ('kinls')
!       kin: 
        v1 = 0.0d0
        call get_kin_ls(v1,       &
                        visual_cfg_cval,&
                        mat_irep, &
                        mat_dim,  &
                        dmat_0%elements,&
                        buffer,   &
                        ao)

        grid_function_values(1:1,ipoint:ipoint) = v1


      case ('kinsl')
!       kin: 
        v1 = 0.0d0
        call get_kin_sl(v1,       &
                        visual_cfg_cval,&
                        mat_irep, &
                        mat_dim,  &
                        dmat_0%elements,&
                        buffer,   &
                        ao)

        grid_function_values(1:1,ipoint:ipoint) = v1


      case ('kintau')
!       kin: 
        v1 = 0.0d0
        call get_kin_tau(v1,       &
                         mat_irep, &
                         mat_dim,  &
                         dmat_0%elements,&
                         buffer,   &
                         ao)

        grid_function_values(1:1,ipoint:ipoint) = v1


      case ('kinlap')
!       kin: 
        v1 = 0.0d0
        call get_kin_lap(v1,       &
                         mat_irep, &
                         mat_dim,  &
                         dmat_0%elements,&
                         buffer,   &
                         ao)

        grid_function_values(1:1,ipoint:ipoint) = v1


      case ('kinnr')
!       kin: 
        v1 = 0.0d0
        call get_kin_tau(v1,       &
                         mat_irep, &
                         mat_dim,  &
                         dmat_0%elements,&
                         buffer,   &
                         ao)
 
        grid_function_values(1:1,ipoint:ipoint) = 0.5d0*v1

        call get_kin_lap(v1,       &
                         mat_irep, &
                         mat_dim,  &
                         dmat_0%elements,&
                         buffer,   &
                         ao)

        grid_function_values(1:1,ipoint:ipoint) = &
                 grid_function_values(1:1,ipoint:ipoint)+0.5d0*v1

      case ('s')
!       spin density: \psi^\dag \Sigma \psi
!       FIXME!
        v2 = 0.0d0
        !write(*,'(a,2i4,3es16.8e2)') 'before get_s NEW dims: ', mat_irep, mat_dim,v2(1),v2(2),v2(3)
        !write(*,'(a,3es16.8e2)') 'before get_s NEW: ',ao(1), ao(20), ao(147)
        call get_s(v2,       &
                   mat_irep, &
                   mat_dim,  &
                   dmat_0%elements,&
                   buffer,   &
                   ao,r=(/px, py, pz/))
        grid_function_values(1:1,ipoint:ipoint) = v2(1)
        grid_function_values(2:2,ipoint:ipoint) = v2(2)
        grid_function_values(3:3,ipoint:ipoint) = v2(3)
        !write(*,'(a,3es16.8e2)') 'after get_s NEW: ', v2(1), v2(2), v2(3)

    end select

  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_esp_functions_in_points(my_rank,&
                                         dmat_0,&
                                         ao,buffer,&
                                         ndim,nr_points,grid_points,&
                                         grid_function,grid_function_values)

    integer, intent(in) :: my_rank
    type(dmat_t), intent(in):: dmat_0
    real(8), intent(in):: ao(*)
    real(8), intent(in):: buffer(*)
    integer, intent(in):: ndim,nr_points
    real(8), intent(in):: grid_points(ndim,nr_points)
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)

    integer  :: mat_dim, mat_nz, mat_irep

    if (grid_function%name .eq. 'espe') then
      mat_dim  = dmat_0%mat_dim
      mat_nz   = dmat_0%mat_nz
      mat_irep = 0
      call get_esp(nr_points,&
                   grid_function_values(1:1,1:nr_points),&
                   mat_irep, &
                   mat_dim,  &
                   dmat_0%elements,   &
                   grid_points(1:3,1:nr_points),&
                   include_nuc_part=.false., &
                   include_el_part=.true., my_rank=my_rank)
!     gosia: figure out why this should be scaled
      grid_function_values(1:1,1:nr_points) = 0.5d0*grid_function_values
    end if

  end subroutine



!  ----------- legacy ----------------------------------------------------------
   subroutine get_quantity_in_point(nr_dmat, &
                                    D,       &
                                    D_0,     &
                                    px,      &
                                    py,      &
                                    pz,      &
                                    ao,      &
                                    buffer,  &
                                    result_max_dim, &
                                    result,  &
                                    grid_function, &
                                    ipoint, &
                                    jb_input)

!     --------------------------------------------------------------------------
      integer, intent(in)  :: nr_dmat
      type(matrix)         :: D(nr_dmat)
      type(matrix)         :: D_0
      real(8), intent(in)  :: px
      real(8), intent(in)  :: py
      real(8), intent(in)  :: pz
      real(8)              :: ao(*)
      real(8)              :: buffer(*)
      integer, intent(in)  :: result_max_dim
      real(8), intent(out):: result(0:result_max_dim)
      real(8), intent(in), optional    :: jb_input(3)
  type(grid_function_t), intent(inout), optional:: grid_function !thiswill not be optional
  integer, intent(in), optional  :: ipoint !thiswill not be optional
  integer:: ioff
!     --------------------------------------------------------------------------
      integer              :: imat, i, j, ij, icomp
      real(8)              :: v(3), v2(3)
      real(8)              :: matrix_33(3, 3)
      real(8)              :: temp_33(3, 3)
      real(8)              :: elf_chi, tau
      real(8)              :: density, esp(1), esp_e(1), esp_n(1)
      real(8)              :: density_gradient(3)
      real(8)              :: density_hessian(3, 3)
      real(8)              :: eigenvalues(3)
      real(8)              :: eigenvectors(3, 3)
      real(8)              :: temp(3), temp2(9)
      integer              :: itemp(3)
      real(8)              :: density_gradient_squared
      real(8)              :: b, c, f, distance_to_origin
      real(8)              :: rx, ry, rz, r(3)

      logical:: gosia_debug_print
!     --------------------------------------------------------------------------
      gosia_debug_print = .false.

      !ao(1:nr_ao_slices*nr_ao_cartesian*max_batch_length) = 0.0d0
      !buffer(1:nr_ao_slices*nr_ao_cartesian*max_batch_length) = 0.0d0
      call get_ao(1, (/px/), (/py/), (/pz/), ao, buffer)

      result = 0.0d0
      do imat = 1, nr_dmat
         if (gosia_debug_print) then
           write(*, *) 'imat, visual_cfg_property(imat), result_max_dim: ', &
                        imat, visual_cfg_property(imat), result_max_dim
           write(*, *) 'result1 before: ', result(1)
           if (result_max_dim > 3) then
             write(*, *) 'result4 before: ', result(4)
           end if
         end if
         if (.not. visual_cfg_skip(imat)) then

!           components of current density: -ec \psi^\dag \alpha_q \psi 
            if (visual_cfg_property(imat) == iq_jx) then
               icomp = 1
               if (levyle) then
                  print *, 'error: Levy-Leblond jx not implemented'
                  stop
               else
                  call get_jcomp(density,  &
                             icomp,        &
                             cval,         &
                             D(imat)%irep, &
                             D(imat)%nrow, &
                             D(imat)%elms, &
                             buffer,       &
                             ao)
               end if
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_jy) then
               icomp = 2
               if (levyle) then
                  print *, 'error: Levy-Leblond jy not implemented'
                  stop
               else
                  call get_jcomp(density,  &
                             icomp,        &
                             cval,         &
                             D(imat)%irep, &
                             D(imat)%nrow, &
                             D(imat)%elms, &
                             buffer,       &
                             ao)
               end if
               result(0) = result(0) + density
            end if

            if (visual_cfg_property(imat) == iq_jz) then
               icomp = 3
               if (levyle) then
                  print *, 'error: Levy-Leblond jz not implemented'
                  stop
               else
                  call get_jcomp(density,  &
                             icomp,        &
                             cval,         &
                             D(imat)%irep, &
                             D(imat)%nrow, &
                             D(imat)%elms, &
                             buffer,       &
                             ao)
               end if
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
               elf_chi = (tau-b)/f

               result(0) = result(0) + 1.0d0/(1.0d0+elf_chi*elf_chi)
            end if

!           electrostatic potential-total
            if ((visual_cfg_property(imat) == iq_esp) .or. &
                (visual_cfg_property(imat) == iq_esprho)) then
               call get_esp(1, esp_e,     &
                            d(imat)%irep, &
                            d(imat)%nrow, &
                            d(imat)%elms, &
                            (/px, py, pz/), &
                            include_nuc_part=.false., &
                            include_el_part=.true.)
               call get_esp(1, esp_n,     &
                            d(imat)%irep, &
                            d(imat)%nrow, &
                            d(imat)%elms, &
                            (/px, py, pz/), &
                            include_nuc_part=.true., &
                            include_el_part=.false.)
               esp = esp_e*0.5d0+esp_n
            end if

!           electrostatic potential-electronic part
            if ((visual_cfg_property(imat) == iq_espe) .or. &
                (visual_cfg_property(imat) == iq_esperho)) then
               call get_esp(1, esp_e,     &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            (/px, py, pz/), &
                            include_nuc_part=.false., &
                            include_el_part=.true.)
               esp = esp_e*0.5d0
            end if

!           electrostatic potential-nuclear part
            if ((visual_cfg_property(imat) == iq_espn) .or. &
                (visual_cfg_property(imat) == iq_espnrho)) then
               call get_esp(1, esp_n,     &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            (/px, py, pz/), &
                            include_nuc_part=.true., &
                            include_el_part=.false.)
               esp = esp_n
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
                     call get_jblao_ll(visual_cfg_london_component, D_0%elms,                      &
                                       ntbas(0),                      &
                                       ao, &
                                       (/px, py, pz/),                &
                                       v2)
                  else
                     call get_jblao(visual_cfg_london_component, cval, D_0%elms,                      &
                                    ntbas(0),                      &
                                    ao, &
                                    (/px, py, pz/),                &
                                    v2)
                  end if
                  v(1:3) = v(1:3) + v2(1:3)
               end if
               result(1:3) = result(1:3) + v(1:3)
               if (gosia_debug_print) then
                 ! gosia: "v2" contribution depends on D_0-be careful when calculaing more than 1 density
                 ! at once (see 'collect matrices for identical densities to save time' in visual.F90)
                 write(*, *) 'after iq_j: visual_cfg_london_component, v2: ', visual_cfg_london_component, v2 
                 write(*, *) 'result1 after iq_j: ', result(1)
                 if (result_max_dim > 3) then
                   write(*, *) 'result4 after iq_j: ', result(4)
                 end if
               end if
            end if

!           spin density: \psi^\dag \Sigma \psi
            if (visual_cfg_property(imat) == iq_s) then
              v  = 0.0d0
              buffer(1:nr_ao_slices*nr_ao_cartesian) = 0.0d0
             !write(*,'(a,2i4,3es16.8e2)') 'before get_s OLD dims: ',D(imat)%irep,D(imat)%nrow,v(1),v(2),v(3)
             !write(*,'(a,3es16.8e2)') 'before get_s OLD: ',ao(1), ao(20), ao(147)
               call get_s(v,            &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          !D_0%irep, &
                          !D_0%nrow, &
                          !D_0%elms, &
                          buffer,       &
                          ao,r=(/px, py, pz/))
               result(1:3) = result(1:3) + v(1:3)
               !write(*,'(a,3es16.8e2)') 'after get_s OLD: ', v(1), v(2), v(3)
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

!           divergence or curl or gradient of current density
            if ((visual_cfg_property(imat) == iq_divj)  .or.  &
                (visual_cfg_property(imat) == iq_rotj)  .or.  &
                (visual_cfg_property(imat) == iq_gradj)) then

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
                     print *, 'error: Levy-Leblond divj and rotj and gradj not implemented'
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
                  matrix_33(1:3, 1:3) = matrix_33(1:3, 1:3) + temp_33(1:3, 1:3)
               end if

               ! the order of elements on matrix_33(1:3, 1:3) from get_gj*subroutines:
               !
               !  (1, 1) (1, 2) (1, 3)        d(jx)/dx d(jy)/dx d(jz)/dx
               !  (2, 1) (2, 2) (2, 3)  -->   d(jx)/dy d(jy)/dy d(jz)/dy
               !  (3, 1) (3, 2) (3, 3)        d(jx)/dz d(jy)/dz d(jz)/dz
               !

               if (visual_cfg_property(imat) == iq_divj) then
                  result(0) = result(0) + matrix_33(1, 1) + matrix_33(2, 2) + matrix_33(3, 3)
               end if
               if (visual_cfg_property(imat) == iq_rotj) then
                  ! gosia: i think the signs should be reversed here (to test):
                  result(1) = result(1) - matrix_33(2, 3) + matrix_33(3, 2)
                  result(2) = result(2) - matrix_33(3, 1) + matrix_33(1, 3)
                  result(3) = result(3) - matrix_33(1, 2) + matrix_33(2, 1)
               end if
               if (visual_cfg_property(imat) == iq_gradj) then
                  result(4)  = result(4)  + matrix_33(1, 1)  ! xx: d(jx)/dx
                  result(5)  = result(5)  + matrix_33(2, 1)  ! xy: d(jx)/dy
                  result(6)  = result(6)  + matrix_33(3, 1)  ! xz: d(jx)/dz
                  result(7)  = result(7)  + matrix_33(1, 2)  ! yx: d(jy)/dx
                  result(8)  = result(8)  + matrix_33(2, 2)  ! yy: d(jy)/dy
                  result(9)  = result(9)  + matrix_33(3, 2)  ! yz: d(jy)/dz
                  result(10) = result(10) + matrix_33(1, 3)  ! zx: d(jz)/dx
                  result(11) = result(11) + matrix_33(2, 3)  ! zy: d(jz)/dy
                  result(12) = result(12) + matrix_33(3, 3)  ! zz: d(jz)/dz
               end if
               if (gosia_debug_print) then
                ! its better not to run it with .J
                ! see debug tests in visual_lao_current
                 write(*, *) 'result1 after iq_gradj: ', result(1)
                 if (result_max_dim > 3) then
                   write(*, *) 'result4 after iq_gradj: ', result(4)
                 end if
               end if
            end if

!           nonrel diamagnetic current density (hardcoded for B along z)
            if (visual_cfg_property(imat) == iq_jdia) then
               call j_dia(matrix_33,               &
                          visual_cfg_london,       &
                          visual_cfg_gauge_origin, &
                          D(imat)%irep,            &
                          D(imat)%nrow,            &
                          D(imat)%elms,            &
                          ao,                      &
                          (/px, py, pz/))

               result(1) = matrix_33(1, 3)
               result(2) = matrix_33(2, 3)
               result(3) = matrix_33(3, 3)
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

               if (present(jb_input)) then
                 v = jb_input
               else
                 call get_j(v,            &
                            cval,         &
                            D(imat)%irep, &
                            D(imat)%nrow, &
                            D(imat)%elms, &
                            buffer,       &
                            ao)
               end if

               rx = px-visual_cfg_gauge_origin(1)
               ry = py-visual_cfg_gauge_origin(2)
               rz = pz-visual_cfg_gauge_origin(3)

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

!           full light-matter interaction: cosine part
            if (visual_cfg_property(imat) == iq_tcos) then

               call get_j(v,            &
                          cval,         &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)

               f = (visual_cfg_wave_vector(1)*px+   &
                    visual_cfg_wave_vector(2)*py+   &
                    visual_cfg_wave_vector(3)*pz)*visual_cfg_freq/cval
               result(0) = result(0) - &
                    (v(1)*visual_cfg_pol_vector(1)+  &
                     v(2)*visual_cfg_pol_vector(2)+  &
                     v(3)*visual_cfg_pol_vector(3))*cos(f)/visual_cfg_freq

            end if

!           full light-matter interaction: sine part
            if (visual_cfg_property(imat) == iq_tsin) then

               call get_j(v,            &
                          cval,         &
                          D(imat)%irep, &
                          D(imat)%nrow, &
                          D(imat)%elms, &
                          buffer,       &
                          ao)

               f = (visual_cfg_wave_vector(1)*px+   &
                    visual_cfg_wave_vector(2)*py+   &
                    visual_cfg_wave_vector(3)*pz)*visual_cfg_freq/cval
               result(0) = result(0) - &
                    (v(1)*visual_cfg_pol_vector(1)+  &
                     v(2)*visual_cfg_pol_vector(2)+  &
                     v(3)*visual_cfg_pol_vector(3))*sin(f)/visual_cfg_freq

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
                  v = -0.5d0*v
               else
                  if (present(jb_input)) then
                    v = jb_input
                  else
                  call get_j(v,            &
                             cval,         &
                             D(imat)%irep, &
                             D(imat)%nrow, &
                             D(imat)%elms, &
                             buffer,       &
                             ao)
               end if
               end if
               if (visual_cfg_london .and. .not. visual_cfg_london_skip_direct) then
                  if (levyle) then
                     call get_jblao_ll(visual_cfg_london_component, D_0%elms,                      &
                                       ntbas(0),                      &
                                       ao, &
                                       (/px, py, pz/),                &
                                       v2)
                     v = v2-v
                  else
                     call get_jblao(visual_cfg_london_component, cval, D_0%elms,                      &
                                    ntbas(0),                      &
                                    ao, &
                                    (/px, py, pz/),                &
                                    v2)
                     v = v2+v
                  end if
               end if

!              we either calculate the NMR shielding of a given nucleus
!              or a "NICS" value in a given point
!              (for a true NICS value multiply by "-1")
               if (visual_cfg_nics) then
                  rx = px-visual_cfg_nics_origin(1)
                  ry = py-visual_cfg_nics_origin(2)
                  rz = pz-visual_cfg_nics_origin(3)
               else
                  rx = px-cord(1, visual_cfg_ref_nucleus(imat))
                  ry = py-cord(2, visual_cfg_ref_nucleus(imat))
                  rz = pz-cord(3, visual_cfg_ref_nucleus(imat))
               end if
                  
               f = -1.0d0/((rx*rx+ry*ry+rz*rz)**1.5d0)
               
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
                          D(imat)%irep,            &
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

!              we either calculate the NMR shielding of a given nucleus
!              or a "NICS" value in a given point
!              (for a true NICS value multiply by "-1")
               if (visual_cfg_nics) then
                  r(1) = px-visual_cfg_nics_origin(1)
                  r(2) = py-visual_cfg_nics_origin(2)
                  r(3) = pz-visual_cfg_nics_origin(3)
               else
                  r(1) = px-cord(1, visual_cfg_ref_nucleus(imat))
                  r(2) = py-cord(2, visual_cfg_ref_nucleus(imat))
                  r(3) = pz-cord(3, visual_cfg_ref_nucleus(imat))
               end if

               v2(1) = r(2)*v(3) - r(3)*v(2)
               v2(2) = r(3)*v(1) - r(1)*v(3)
               v2(3) = r(1)*v(2) - r(2)*v(1)

               f = -1.0d0/((r(1)*r(1) + r(2)*r(2) + r(3)*r(3))**1.5d0)
               result(0) = result(0) + f*v2(i)

            end if

         end if
         if (gosia_debug_print) then
           write(*, *) 'imat, visual_cfg_property(imat), result_max_dim: ', &
                        imat, visual_cfg_property(imat), result_max_dim
           write(*, *) 'result1 after: ', result(1)
           if (result_max_dim > 3) then
             write(*, *) 'result4 after: ', result(4)
           end if
         end if
      end do

!     Multiply with Cartesian powers
      do i = 1, visual_cfg_cartesian_power(1)
        result = px*result    
      enddo
      do i = 1, visual_cfg_cartesian_power(2)
        result = py*result    
      enddo
      do i = 1, visual_cfg_cartesian_power(3)
        result = pz*result    
      enddo
      
!     Scale result 
      result = visual_cfg_scale*result

   end subroutine

end module
