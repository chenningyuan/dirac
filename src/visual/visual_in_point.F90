module visual_in_point

  use memory_allocator
  use electrostatic_potential
  use interface_ao
  use interface_mo
  use visual_cfg
  use dirac_ao_eval
  use density_eval
  use matrix_defop_old
  use xc_london_c1

  use visual_datatypes
  use dirac_dmat_factory
  use visual_dirac

  implicit none

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

  interface get_quantity_in_point_refactored
    module procedure get_quantity_in_point_single
    module procedure get_quantity_in_point_composite
  end interface

contains

! ------------------------------------------------------------------------------
  subroutine get_quantity_in_point_single(dmat, &
                                          ao,buffer,&
                                          grid_function,&
                                          grid_function_values,&
                                          px,py,pz,&
                                          mol,&
                                          ipoint)
    type(dmat_t), intent(in):: dmat
    type(mol_t), intent(in):: mol
    real(8), intent(in):: ao(*)
    real(8), intent(in):: buffer(*)
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)
    real(8), intent(in):: px, py, pz
    integer, intent(in), optional  :: ipoint !thiswill not be optional

    integer:: mat_dim, mat_nz, mat_irep
    real(8), pointer::mat_elements(:)
    integer:: imat,i,j,itemp(3)
    real(8):: a, b, c, f, t0, t1(3), t2_33(3,3)
    real(8):: temp1(1)
    real(8):: elf_chi, tau
    real(8):: eigenvalues(3),eigenvectors(3, 3)
    real(8):: temp(3), temp2(9), temp_33(3,3)
    integer:: ioff
    real(8):: rx, ry, rz


     !px = grid%points(1, ipoint)
     !py = grid%points(2, ipoint)
     !pz = grid%points(3, ipoint)
     !todo: call grid%points%get_cartesian_coords()

    select case (grid_function%name)

      case ('gamma5')
!       gamma5 density: \psi^\dag \gamma^5 \psi
        t0 = 0.0d0
        call get_gamma5(t0,        &
                        dmat%irep, &
                        dmat%mat_dim,  &
                        dmat%elements, &
                        buffer,    &
                        ao)
        grid_function_values(1:1,ipoint:ipoint) = t0

      case ('ed')
!       electron density
        t0 = 0.0d0
        call get_n(t0,        &
                   dmat%irep, &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   buffer,    &
                   ao)
        grid_function_values(1:1,ipoint:ipoint) = t0

      case ('ed_gradient')
!       electron density gradient
        t0 = 0.0d0
        t1 = 0.0d0
        call get_gn(t0,        &
                    t1,        &
                    dmat%irep, &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,    &
                    ao)               
        grid_function_values(1:1,ipoint:ipoint) = t1(1)
        grid_function_values(2:2,ipoint:ipoint) = t1(2)
        grid_function_values(3:3,ipoint:ipoint) = t1(3)

      case ('ed_hessian')
!       gosia fixme:
!       electron density Hessian is implemented for hermitian density matrices (such as dmat)
!       todo: add a warning
        t0 = 0.0d0
        t1 = 0.0d0
        t2_33 = 0.0d0
        call get_hn(t0,       &
                    t1,       &
                    t2_33,       &
                    dmat%irep,     &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,   &
                    ao)               
        grid_function_values(1:1,ipoint:ipoint) = t2_33(1, 1)   ! xx
        grid_function_values(2:2,ipoint:ipoint) = t2_33(1, 2)   ! xy
        grid_function_values(3:3,ipoint:ipoint) = t2_33(1, 3)   ! xz
        grid_function_values(4:4,ipoint:ipoint) = t2_33(2, 1)   ! yx
        grid_function_values(5:5,ipoint:ipoint) = t2_33(2, 2)   ! yy
        grid_function_values(6:6,ipoint:ipoint) = t2_33(2, 3)   ! yz
        grid_function_values(7:7,ipoint:ipoint) = t2_33(3, 1)   ! zx
        grid_function_values(8:8,ipoint:ipoint) = t2_33(3, 2)   ! zy
        grid_function_values(9:9,ipoint:ipoint) = t2_33(3, 3)   ! zz

      case ('rdg')
!       reduced density gradient
        t0 = 0.0d0
        t1 = 0.0d0
        call get_gn(t0,       &
                    t1,       &
                    dmat%irep,     &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,   &
                    ao)

        b = sqrt(t1(1)*t1(1) + t1(2)*t1(2) + t1(3)*t1(3))
        c = t0**(4.0d0/3.0d0)
!       0.161620459674d0 = 1/(2.0d0*((3.0d0*pi*pi)**(1.0d0/3.0d0))
        grid_function_values(1:1,ipoint:ipoint) = 0.161620459674d0*(b/c)

      case ('ed_sign_lambda2')
!       sign(lambda2), where lambda2=the second eigenvalue of the Hessian
!       of the electron density
!
!       gosia fixme:
!       electron density Hessian is implemented for hermitian density matrices (such as dmat)
!       todo: add a warning
        t0 = 0.0d0
        t1 = 0.0d0
        t2_33 = 0.0d0
        call get_hn(t0,       &
                    t1,       &
                    t2_33,       &
                    dmat%irep,     &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,   &
                    ao)               
        eigenvalues  = 0.0d0
        eigenvectors = 1.0d0
        temp         = 0.0d0
        itemp        = 0
        temp2 = reshape(t2_33, (/9/))

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
           grid_function_values(1:1,ipoint:ipoint) = t0
        else if (eigenvalues(2) .lt. 0.0d0) then
           grid_function_values(1:1,ipoint:ipoint) = -t0
        end if


      case ('ed_laplacian')
!
!       gosia fixme:
!       density Hessian is implemented for hermitian density matrices (such as the unperturbed density matrix)
!       todo: add a warning
        t0 = 0.0d0
        t1 = 0.0d0
        t2_33 = 0.0d0
        call get_hn(t0,        &
                    t1,        &
                    t2_33,        &
                    dmat%irep, &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,    &
                    ao)               
        eigenvectors = 1.0d0
        temp         = 0.0d0
        itemp        = 0
        temp2 = reshape(t2_33, (/9/))

        call jaco(temp2, eigenvectors, 3, 3, 3, temp, itemp)

        grid_function_values(1:1,ipoint:ipoint) = temp2(1) + temp2(3) + temp2(6)

      case ('elf')
!       elf: electron localization function
        tau = 0.0d0
        t0 = 0.0d0
        t1 = 0.0d0
        call get_kin_tau(tau,       &
                         dmat%irep, &
                         dmat%mat_dim,  &
                         dmat%elements, &
                         buffer,    &
                         ao)
        call get_gn(t0, &
                    t1, &
                    dmat%irep, &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,    &
                    ao)

        a = t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3)
        b = a/(8.0d0*t0)

!           2.87123400018819d0 = 0.3d0*((3.0d0*pi*pi)**(2.0d0/3.0d0))
        f = 2.87123400018819d0*(t0**(5.0d0/3.0d0))
        elf_chi = (tau-b)/f

        grid_function_values(1:1,ipoint:ipoint) = 1.0d0/(1.0d0+elf_chi*elf_chi)

      case ('kin')
!       kin: 
        t0 = 0.0d0
        call get_kin(t0,        &
                     rs_input%cvalue, &
                     dmat%irep, &
                     dmat%mat_dim,    &
                     dmat%elements,   &
                     buffer,    &
                     ao)

        grid_function_values(1:1,ipoint:ipoint) = t0


      case ('kinls')
!       kin: 
        t0 = 0.0d0
        call get_kin_ls(t0,        &
                        rs_input%cvalue, &
                        dmat%irep, &
                        dmat%mat_dim,    &
                        dmat%elements,   &
                        buffer,    &
                        ao)

        grid_function_values(1:1,ipoint:ipoint) = t0


      case ('kinsl')
!       kin: 
        t0 = 0.0d0
        call get_kin_sl(t0,       &
                        rs_input%cvalue,&
                        dmat%irep,      &
                        dmat%mat_dim,   &
                        dmat%elements,  &
                        buffer,   &
                        ao)

        grid_function_values(1:1,ipoint:ipoint) = t0


      case ('kintau')
!       kin: 
        t0 = 0.0d0
        call get_kin_tau(t0,        &
                         dmat%irep, &
                         dmat%mat_dim,  &
                         dmat%elements, &
                         buffer,    &
                         ao)

        grid_function_values(1:1,ipoint:ipoint) = t0


      case ('kinlap')
!       kin: 
        t0 = 0.0d0
        call get_kin_lap(t0,        &
                         dmat%irep, &
                         dmat%mat_dim,  &
                         dmat%elements, &
                         buffer,    &
                         ao)

        grid_function_values(1:1,ipoint:ipoint) = t0


      case ('kinnr')
!       kin: 
        t0 = 0.0d0
        call get_kin_tau(t0,            &
                         dmat%irep,     &
                         dmat%mat_dim,  &
                         dmat%elements, &
                         buffer,        &
                         ao)
 
        grid_function_values(1:1,ipoint:ipoint) = 0.5d0*t0

        call get_kin_lap(t0,         &
                         dmat%irep,  &
                         dmat%mat_dim,  &
                         dmat%elements, &
                         buffer,     &
                         ao)

        grid_function_values(1:1,ipoint:ipoint) = &
                 grid_function_values(1:1,ipoint:ipoint)+0.5d0*t0

      case ('s')
!       spin density: \psi^\dag \Sigma \psi
        t1 = 0.0d0
        call get_s(t1,        &
                   dmat%irep, &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   buffer,    &
                   ao)
        grid_function_values(1:1,ipoint:ipoint) = t1(1)
        grid_function_values(2:2,ipoint:ipoint) = t1(2)
        grid_function_values(3:3,ipoint:ipoint) = t1(3)

      case ('divs')
!       divergence of spin density:
        t1 = 0.0d0
        t2_33 = 0.0d0
        call get_gs(t1,        &
                    t2_33,        &
                    dmat%irep, &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,    &
                    ao)
        grid_function_values(1:1,ipoint:ipoint) = t2_33(1,1) + t2_33(2,2) + t2_33(3,3)

      case ('rots')
!       divergence of spin density:
        t1 = 0.0d0
        t2_33 = 0.0d0
        call get_gs(t1,        &
                    t2_33,        &
                    dmat%irep, &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,    &
                    ao)
        grid_function_values(1:1,ipoint:ipoint) = -t2_33(2,3) + t2_33(3,2)
        grid_function_values(2:2,ipoint:ipoint) = -t2_33(3,1) + t2_33(1,3)
        grid_function_values(3:3,ipoint:ipoint) = -t2_33(1,2) + t2_33(2,1)

      case ('j')
!       current density: -ec \psi^\dag \alpha \psi

        if (rs_wf%hamiltonian=='LEVYLE') then
          t0 = 0.0d0
          t1 = 0.0d0
          call get_gn(t0,        &
                      t1,        &
                      dmat%irep, &
                      dmat%mat_dim,  &
                      dmat%elements, &
                      buffer,    &
                      ao)               
          grid_function_values(1:1,ipoint:ipoint) = 0.5d0*t1(1)
          grid_function_values(2:2,ipoint:ipoint) = 0.5d0*t1(2)
          grid_function_values(3:3,ipoint:ipoint) = 0.5d0*t1(3)

        else
          t1 = 0.0d0
          call get_j(t1,        &
                     rs_input%cvalue,      &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements, &
                     buffer,    &
                     ao)
          grid_function_values(1:1,ipoint:ipoint) = t1(1)
          grid_function_values(2:2,ipoint:ipoint) = t1(2)
          grid_function_values(3:3,ipoint:ipoint) = t1(3)

        end if

      case ('divj', &
            'gradj',&
            'rotj')
!       divergence or curl or gradient of current density

        if (rs_wf%hamiltonian=='LEVYLE') then
           print *, 'error: Levy-Leblond not implemented'
           stop
        end if

        t1 = 0.0d0
        t2_33 = 0.0d0

        call get_gj(t1,            &
                    t2_33,            &
                    rs_input%cvalue,          &
                    dmat%irep,     &
                    dmat%mat_dim,  &
                    dmat%elements, &
                    buffer,        &
                    ao)
!
!       the order of elements on t2_33(1:3, 1:3) from get_gj* subroutines is as follows:
!      
!       (1, 1) (1, 2) (1, 3)        d(jx)/dx d(jy)/dx d(jz)/dx
!       (2, 1) (2, 2) (2, 3)  -->   d(jx)/dy d(jy)/dy d(jz)/dy
!       (3, 1) (3, 2) (3, 3)        d(jx)/dz d(jy)/dz d(jz)/dz
!      
        if (grid_function%name == 'divj') then
          grid_function_values(1:1,ipoint:ipoint) = t2_33(1,1) + t2_33(2,2) + t2_33(3,3)
        end if

        if (grid_function%name == 'rotj') then
!         gosia fixme todo: i think the signs should be reversed here (to test):
          grid_function_values(1:1,ipoint:ipoint) = -t2_33(2,3) + t2_33(3,2)
          grid_function_values(2:2,ipoint:ipoint) = -t2_33(3,1) + t2_33(1,3)
          grid_function_values(3:3,ipoint:ipoint) = -t2_33(1,2) + t2_33(2,1)
        end if

        if (grid_function%name == 'gradj') then
          grid_function_values(1:1,ipoint:ipoint) = t2_33(1, 1)  ! xx: d(jx)/dx
          grid_function_values(2:2,ipoint:ipoint) = t2_33(2, 1)  ! xy: d(jx)/dy
          grid_function_values(3:3,ipoint:ipoint) = t2_33(3, 1)  ! xz: d(jx)/dz
          grid_function_values(4:4,ipoint:ipoint) = t2_33(1, 2)  ! yx: d(jy)/dx
          grid_function_values(5:5,ipoint:ipoint) = t2_33(2, 2)  ! yy: d(jy)/dy
          grid_function_values(6:6,ipoint:ipoint) = t2_33(3, 2)  ! yz: d(jy)/dz
          grid_function_values(7:7,ipoint:ipoint) = t2_33(1, 3)  ! zx: d(jz)/dx
          grid_function_values(8:8,ipoint:ipoint) = t2_33(2, 3)  ! zy: d(jz)/dy
          grid_function_values(9:9,ipoint:ipoint) = t2_33(3, 3)  ! zz: d(jz)/dz
        end if

! gosia todo fixme - there was a comment that it is "better not to run it with .J,
! see debug tests in visual_lao_current", check what was this about

      case ('jdiaBz')
!       nonrel diamagnetic current density (hardcoded for B along z)
        t2_33 = 0.0d0
        call j_dia(t2_33,      &
                   (grid_function%lao > 0),&
                   rs_input%gauge_origin,&
                   dmat%irep,     &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   ao,            &
                   (/px, py, pz/))

        grid_function_values(1:1,ipoint:ipoint) = t2_33(1, 3)
        grid_function_values(2:2,ipoint:ipoint) = t2_33(2, 3)
        grid_function_values(3:3,ipoint:ipoint) = t2_33(3, 3)

      case ('bdip')
!       responding operator is magnetic dipole

        if (rs_wf%hamiltonian=='LEVYLE') then
           print *, 'error: Levy-Leblond not implemented'
           stop
        else
          t1 = 0.0d0
          call get_j(t1,        &
                     rs_input%cvalue,      &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements, &
                     buffer,    &
                     ao)
        end if

        rx = px-rs_input%gauge_origin(1)
        ry = py-rs_input%gauge_origin(2)
        rz = pz-rs_input%gauge_origin(3)

        f = -0.5d0

!       'bdipx'
        grid_function_values(1:1,ipoint:ipoint) = f*(ry*t1(3) - rz*t1(2))
!       'bdipy'
        grid_function_values(2:2,ipoint:ipoint) = f*(rz*t1(1) - rx*t1(3))
!       'bdipz'
        grid_function_values(3:3,ipoint:ipoint) = f*(rx*t1(2) - ry*t1(1))

      case ('ndip')
!       responding operator is nuclear magnetic dipole

        if (rs_wf%hamiltonian=='LEVYLE') then
          t0 = 0.0d0
          t1 = 0.0d0
          call get_gn(t0,        &
                      t1,        &
                      dmat%irep, &
                      dmat%mat_dim,  &
                      dmat%elements, &
                      buffer,    &
                      ao)               
          t1(1) = -0.5d0*t1(1)
          t1(2) = -0.5d0*t1(2)
          t1(3) = -0.5d0*t1(3)

        else

          t1 = 0.0d0
          call get_j(t1,        &
                     rs_input%cvalue,      &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements, &
                     buffer,    &
                     ao)
        end if

!       we either calculate the NMR shielding of a given nucleus
!       or a "NICS" value in a given point
!       (for a true NICS value multiply by "-1")
        !if (visual_cfg_nics) then
           !FIXME
           !rx = px-visual_cfg_nics_origin(1)
           !ry = py-visual_cfg_nics_origin(2)
           !rz = pz-visual_cfg_nics_origin(3)
        !else
           rx = px-mol%geometry(1, grid_function%nuc)
           ry = py-mol%geometry(2, grid_function%nuc)
           rz = pz-mol%geometry(3, grid_function%nuc)
        !end if

        f = -1.0d0/((rx*rx+ry*ry+rz*rz)**1.5d0)

!       'ndipx'
        grid_function_values(1:1,ipoint:ipoint) = f*(ry*t1(3) - rz*t1(2))
!       'ndipy'
        grid_function_values(2:2,ipoint:ipoint) = f*(rz*t1(1) - rx*t1(3))
!       'ndipz'
        grid_function_values(3:3,ipoint:ipoint) = f*(rx*t1(2) - ry*t1(1))


      case ('edip')
!       responding operator is electric dipole

        t0 = 0.0d0
        call get_n(t0,        &
                   dmat%irep, &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   buffer,    &
                   ao)

!       'edipx'
        grid_function_values(1:1,ipoint:ipoint) = px*t0
!       'edipy'
        grid_function_values(2:2,ipoint:ipoint) = py*t0
!       'edipz'
        grid_function_values(3:3,ipoint:ipoint) = pz*t0

      case ('esperho')
!       electrostatic potential-electronic part multiplied with density

        temp1 = 0.0d0
        call get_esp(1,     &
                     temp1, &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     (/px, py, pz/), &
                     include_nuc_part=.false., &
                     include_el_part=.true.)

        t0 = 0.0d0
        call get_n(t0,        &
                   dmat%irep, &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   buffer,    &
                   ao)
!       gosia: figure out why this should be scaled
        grid_function_values(1:1,ipoint:ipoint) = 0.5d0*temp1(1)*t0

      case ('espnrho')
!       electrostatic potential-nuclear part multiplied with density

        temp1 = 0.0d0
        call get_esp(1,     &
                     temp1, &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     (/px, py, pz/), &
                     include_nuc_part=.true., &
                     include_el_part=.false.)

        t0 = 0.0d0
        call get_n(t0,        &
                   dmat%irep, &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   buffer,    &
                   ao)
        grid_function_values(1:1,ipoint:ipoint) = temp1(1)*t0

      case ('esprho')
!       electrostatic potential-total multiplied with density

        temp1 = 0.0d0
        call get_esp(1,     &
                     temp1, &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     (/px, py, pz/), &
                     include_nuc_part=.false., &
                     include_el_part=.true.)
        grid_function_values(1:1,ipoint:ipoint) = 0.5d0*temp1(1)

        temp1 = 0.0d0
        call get_esp(1,     &
                     temp1, &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     (/px, py, pz/), &
                     include_nuc_part=.true., &
                     include_el_part=.false.)

        grid_function_values(1:1,ipoint:ipoint) = grid_function_values(1:1,ipoint:ipoint) + temp1(1)

        t0 = 0.0d0
        call get_n(t0,        &
                   dmat%irep, &
                   dmat%mat_dim,  &
                   dmat%elements, &
                   buffer,    &
                   ao)
        grid_function_values(1:1,ipoint:ipoint) = grid_function_values(1:1,ipoint:ipoint)*t0


 
      case ('bedcos')
!       full light-matter interaction: cosine part

!        t1 = 0.0d0
!!       hardcoded for testing:
!        CALL ANGLES(0.0d0,0.0d0,0.0d0,visual_cfg_wave_vector,visual_cfg_pol_vector)
!        visual_cfg_freq = 49.138
!
!        call get_j(t1,        &
!                   rs_input%cvalue,      &
!                   dmat%irep, &
!                   dmat%mat_dim,  &
!                   dmat%elements, &
!                   buffer,    &
!                   ao)
!
!        f = (visual_cfg_wave_vector(1)*px+   &
!             visual_cfg_wave_vector(2)*py+   &
!             visual_cfg_wave_vector(3)*pz)*visual_cfg_freq/rs_input%cvalue
!        grid_function_values(1:1,ipoint:ipoint) = & 
!            -(t1(1)*visual_cfg_pol_vector(1)+  &
!              t1(2)*visual_cfg_pol_vector(2)+  &
!              t1(3)*visual_cfg_pol_vector(3))*cos(f)/visual_cfg_freq
!

      case ('bedsin')
!       full light-matter interaction: sine part

!        t1 = 0.0d0
!!       hardcoded for testing:
!        CALL ANGLES(0.0d0,0.0d0,0.0d0,visual_cfg_wave_vector,visual_cfg_pol_vector)
!        visual_cfg_freq = 49.138
!
!
!        call get_j(t1,        &
!                   rs_input%cvalue,      &
!                   dmat%irep, &
!                   dmat%mat_dim,  &
!                   dmat%elements, &
!                   buffer,    &
!                   ao)
!
!        f = (visual_cfg_wave_vector(1)*px+   &
!             visual_cfg_wave_vector(2)*py+   &
!             visual_cfg_wave_vector(3)*pz)*visual_cfg_freq/rs_input%cvalue
!        grid_function_values(1:1,ipoint:ipoint) = &
!            -(t1(1)*visual_cfg_pol_vector(1)+  &
!              t1(2)*visual_cfg_pol_vector(2)+  &
!              t1(3)*visual_cfg_pol_vector(3))*sin(f)/visual_cfg_freq

    end select

    if (grid_function%modify == 1) then
      call modify_grid_function_pointwise(grid_function,&
                                          grid_function_values,&
                                          px,py,pz,&
                                          ipoint=ipoint)
    end if

  end subroutine



! ------------------------------------------------------------------------------
  subroutine get_quantity_in_point_composite(dmat_0, dmat_other, &
                                             ao,buffer,&
                                             grid_function,&
                                             grid_function_values,&
                                             px,py,pz,&
                                             mol,&
                                             ipoint)

    type(dmat_t), intent(in):: dmat_0
    type(dmat_t), intent(in):: dmat_other
    type(mol_t), intent(in):: mol
    real(8), intent(in):: ao(*)
    real(8), intent(in):: buffer(*)
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)
    real(8), intent(in):: px, py, pz
    integer, intent(in), optional  :: ipoint !thiswill not be optional

    integer:: mat_dim, mat_nz, mat_irep
    real(8), pointer::mat_elements(:)
    integer:: imat,i,j,itemp(3)
    real(8):: a, b, c, f, t0, t1(3), t2_33(3,3)
    real(8):: elf_chi, tau
    real(8):: eigenvalues(3),eigenvectors(3, 3)
    real(8):: temp(3), temp2(9), temp_33(3,3)
  integer:: ioff
  character(1):: xyzcomp
    real(8):: rx, ry, rz


     !px = grid%points(1, ipoint)
     !py = grid%points(2, ipoint)
     !pz = grid%points(3, ipoint)
     !todo: call grid%points%get_cartesian_coords()
    select case (grid_function%name)

      case ('j')
!       current density: -ec \psi^\dag \alpha \psi

        if (rs_wf%hamiltonian=='LEVYLE') then
          t0 = 0.0d0
          t1 = 0.0d0
          call get_gn(t0,        &
                      t1,        &
                      dmat_other%irep, &
                      dmat_other%mat_dim,  &
                      dmat_other%elements, &
                      buffer,    &
                      ao)               
          grid_function_values(1:1,ipoint:ipoint) = 0.5d0*t1(1)
          grid_function_values(2:2,ipoint:ipoint) = 0.5d0*t1(2)
          grid_function_values(3:3,ipoint:ipoint) = 0.5d0*t1(3)

        else

          t1 = 0.0d0
          call get_j(t1,       &
                     rs_input%cvalue,     &
                     dmat_other%irep, &
                     dmat_other%mat_dim,  &
                     dmat_other%elements, &
                     buffer,   &
                     ao)
          grid_function_values(1:1,ipoint:ipoint) = t1(1)
          grid_function_values(2:2,ipoint:ipoint) = t1(2)
          grid_function_values(3:3,ipoint:ipoint) = t1(3)

        end if

        if (grid_function%lao > 0) then
          if (grid_function%lao_skip(2) == 0) then
!           do not skip direct contribution ('nodirect')

            if (rs_wf%hamiltonian=='LEVYLE') then

              t1 = 0.0d0
              call get_jblao_ll(grid_function%lao_comp,&
                                dmat_0%elements, &
                                dmat_0%mat_dim,  &
                                ao, &
                                (/px, py, pz/),  &
                                t1)
            else
          
              t1 = 0.0d0
              call get_jblao(grid_function%lao_comp,&
                             rs_input%cvalue,&
                             dmat_0%elements, &
                             dmat_0%mat_dim,  &
                             ao, &
                             (/px, py, pz/),  &
                             t1)
            end if
                    
            grid_function_values(1:1,ipoint:ipoint) = &
              grid_function_values(1:1,ipoint:ipoint) + t1(1)
            grid_function_values(2:2,ipoint:ipoint) = &
              grid_function_values(2:2,ipoint:ipoint) + t1(2)
            grid_function_values(3:3,ipoint:ipoint) = &
              grid_function_values(3:3,ipoint:ipoint) + t1(3)

          end if
        end if

      case ('divj', &
            'gradj',&
            'rotj')
!       divergence or curl or gradient of current density

        if (rs_wf%hamiltonian=='LEVYLE') then
           print *, 'error: Levy-Leblond not implemented'
           stop
        end if

        t1 = 0.0d0
        t2_33 = 0.0d0

        call get_gj(t1,            &
                    t2_33,            &
                    rs_input%cvalue,          &
                    dmat_other%irep,     &
                    dmat_other%mat_dim,  &
                    dmat_other%elements, &
                    buffer,        &
                    ao)

        if (grid_function%lao > 0) then
          if (grid_function%lao_skip(2) == 0) then
!           do not skip direct contribution ('nodirect')
!
            if (rs_wf%hamiltonian=='LEVYLE') then
                print *, 'error: Levy-Leblond divj and rotj and gradj not implemented'
                stop
            else
              temp_33 = 0.0d0
              call get_gjblao(grid_function%lao_comp,&
                              rs_input%cvalue,            &
                              dmat_0%elements, &
                              dmat_0%mat_dim,  &
                              ao,              &
                              (/px, py, pz/),  &
                              temp_33)
            end if
            t2_33(1:3, 1:3) = t2_33(1:3, 1:3) + temp_33(1:3, 1:3)
          end if
        end if

!
!       the order of elements on t2_33(1:3, 1:3) from get_gj* subroutines is as follows:
!      
!       (1, 1) (1, 2) (1, 3)        d(jx)/dx d(jy)/dx d(jz)/dx
!       (2, 1) (2, 2) (2, 3)  -->   d(jx)/dy d(jy)/dy d(jz)/dy
!       (3, 1) (3, 2) (3, 3)        d(jx)/dz d(jy)/dz d(jz)/dz
!      
        if (grid_function%name == 'divj') then
          grid_function_values(1:1,ipoint:ipoint) = t2_33(1,1) + t2_33(2,2) + t2_33(3,3)
        end if

        if (grid_function%name == 'rotj') then
!         gosia fixme todo: i think the signs should be reversed here (to test):
          grid_function_values(1:1,ipoint:ipoint) = -t2_33(2,3) + t2_33(3,2)
          grid_function_values(2:2,ipoint:ipoint) = -t2_33(3,1) + t2_33(1,3)
          grid_function_values(3:3,ipoint:ipoint) = -t2_33(1,2) + t2_33(2,1)
        end if

        if (grid_function%name == 'gradj') then
          grid_function_values(1:1,ipoint:ipoint) = t2_33(1, 1)  ! xx: d(jx)/dx
          grid_function_values(2:2,ipoint:ipoint) = t2_33(2, 1)  ! xy: d(jx)/dy
          grid_function_values(3:3,ipoint:ipoint) = t2_33(3, 1)  ! xz: d(jx)/dz
          grid_function_values(4:4,ipoint:ipoint) = t2_33(1, 2)  ! yx: d(jy)/dx
          grid_function_values(5:5,ipoint:ipoint) = t2_33(2, 2)  ! yy: d(jy)/dy
          grid_function_values(6:6,ipoint:ipoint) = t2_33(3, 2)  ! yz: d(jy)/dz
          grid_function_values(7:7,ipoint:ipoint) = t2_33(1, 3)  ! zx: d(jz)/dx
          grid_function_values(8:8,ipoint:ipoint) = t2_33(2, 3)  ! zy: d(jz)/dy
          grid_function_values(9:9,ipoint:ipoint) = t2_33(3, 3)  ! zz: d(jz)/dz
        end if

! gosia todo fixme - there was a comment that it is "better not to run it with .J,
! see debug tests in visual_lao_current", check what was this about

      case ('jdiaBz')
!       nonrel diamagnetic current density (hardcoded for B along z)
        t2_33 = 0.0d0
        call j_dia(t2_33,      &
                   (grid_function%lao > 0),&
                   rs_input%gauge_origin,&
                   dmat_other%irep,     &
                   dmat_other%mat_dim,  &
                   dmat_other%elements, &
                   ao,            &
                   (/px, py, pz/))

        grid_function_values(1:1,ipoint:ipoint) = t2_33(1, 3)
        grid_function_values(2:2,ipoint:ipoint) = t2_33(2, 3)
        grid_function_values(3:3,ipoint:ipoint) = t2_33(3, 3)


      case ('ndip')
!       responding operator is nuclear magnetic dipole

        if (rs_wf%hamiltonian=='LEVYLE') then
          t0 = 0.0d0
          t1 = 0.0d0
          call get_gn(t0,        &
                      t1,        &
                      dmat_other%irep, &
                      dmat_other%mat_dim,  &
                      dmat_other%elements, &
                      buffer,    &
                      ao)               
          t1(1) = -0.5d0*t1(1)
          t1(2) = -0.5d0*t1(2)
          t1(3) = -0.5d0*t1(3)

        else

          t1 = 0.0d0
          
          call get_j(t1,        &
                     rs_input%cvalue,      &
                     dmat_other%irep, &
                     dmat_other%mat_dim,  &
                     dmat_other%elements, &
                     buffer,    &
                     ao)
        end if

        if (grid_function%lao > 0) then
          if (grid_function%lao_skip(2) == 0) then
!           do not skip direct contribution ('nodirect')

            if (rs_wf%hamiltonian=='LEVYLE') then

              temp = 0.0d0
              call get_jblao_ll(grid_function%lao_comp,&
                                dmat_0%elements, &
                                dmat_0%mat_dim,  &
                                ao, &
                                (/px, py, pz/),  &
                                temp)
              t1(1:3) = temp(1:3) - t1(1:3)
            else
          
              temp = 0.0d0
              call get_jblao(grid_function%lao_comp,&
                             rs_input%cvalue,            &
                             dmat_0%elements, &
                             dmat_0%mat_dim,  &
                             ao,              &
                             (/px, py, pz/),  &
                             temp)
              t1(1:3) = temp(1:3) + t1(1:3)
            end if
          end if
        end if


!       we either calculate the NMR shielding of a given nucleus
!       or a "NICS" value in a given point
!       (for a true NICS value multiply by "-1")
        !if (visual_cfg_nics) then
           !FIXME
           !rx = px-visual_cfg_nics_origin(1)
           !ry = py-visual_cfg_nics_origin(2)
           !rz = pz-visual_cfg_nics_origin(3)
        !else
           rx = px-mol%geometry(1, grid_function%nuc)
           ry = py-mol%geometry(2, grid_function%nuc)
           rz = pz-mol%geometry(3, grid_function%nuc)
        !end if

        f = -1.0d0/((rx*rx+ry*ry+rz*rz)**1.5d0)

!       'ndipx'
        grid_function_values(1:1,ipoint:ipoint) = f*(ry*t1(3) - rz*t1(2))
!       'ndipy'
        grid_function_values(2:2,ipoint:ipoint) = f*(rz*t1(1) - rx*t1(3))
!       'ndipz'
        grid_function_values(3:3,ipoint:ipoint) = f*(rx*t1(2) - ry*t1(1))


    end select

    if (grid_function%modify == 1) then
      call modify_grid_function_pointwise(grid_function,&
                                          grid_function_values,&
                                          px,py,pz,&
                                          ipoint=ipoint)
    end if

  end subroutine


! ------------------------------------------------------------------------------
  subroutine get_esp_functions_in_points(my_rank,&
                                         dmat,&
                                         ao,buffer,&
                                         ndim,nr_points,grid_points,&
                                         grid_function,grid_function_values)

    integer, intent(in) :: my_rank
    type(dmat_t), intent(in):: dmat
    real(8), intent(in):: ao(*)
    real(8), intent(in):: buffer(*)
    integer, intent(in):: ndim,nr_points
    real(8), intent(in):: grid_points(ndim,nr_points)
    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(out):: grid_function_values(:,:)
    real(8), allocatable:: temp(:,:)

    integer  :: mat_dim, mat_nz, mat_irep

    select case (grid_function%name)

      case ('espe')
!       electrostatic potential-electronic part

        call get_esp(nr_points,&
                     grid_function_values(1:1,1:nr_points), &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     grid_points(1:3,1:nr_points),&
                     include_nuc_part=.false., &
                     include_el_part=.true., my_rank=my_rank)
!       gosia: figure out why this should be scaled
        grid_function_values(1:1,1:nr_points) = 0.5d0*grid_function_values(1:1,1:nr_points)

      case ('espn')
!       electrostatic potential-nuclear part

        call get_esp(nr_points,&
                     grid_function_values(1:1,1:nr_points),&
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     grid_points(1:3,1:nr_points),&
                     include_nuc_part=.true., &
                     include_el_part=.false., my_rank=my_rank)

      case ('esp')
!       electrostatic potential-total

        allocate(temp(1,nr_points))
        temp = 0.0d0
        call get_esp(nr_points,&
                     temp, &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     grid_points(1:3,1:nr_points),&
                     include_nuc_part=.false., &
                     include_el_part=.true., my_rank=my_rank)
!       gosia: figure out why this should be scaled
        grid_function_values(1:1,1:nr_points) = 0.5d0*temp(1:1,1:nr_points)

        temp = 0.0d0
        call get_esp(nr_points,&
                     temp, &
                     dmat%irep, &
                     dmat%mat_dim,  &
                     dmat%elements,   &
                     grid_points(1:3,1:nr_points),&
                     include_nuc_part=.true., &
                     include_el_part=.false., my_rank=my_rank)
        grid_function_values(1:1,1:nr_points) = grid_function_values(1:1,1:nr_points)+temp(1:1,1:nr_points)

        deallocate(temp)

    end select

    if (grid_function%modify == 1) then
      call modify_grid_function(grid_function,&
                                grid_function_values,&
                                ndim, nr_points,&
                                grid_points)
    end if

  end subroutine

! ------------------------------------------------------------------------------
  subroutine modify_grid_function_pointwise(grid_function,&
                                            grid_function_values,&
                                            px,py,pz,&
                                            ipoint)

    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(inout):: grid_function_values(:,:)
    real(8), intent(in):: px, py, pz
    integer, intent(in):: ipoint
    real(8):: f, r2
    logical:: modified

    f = 1.0d0

!   Multiply with Cartesian powers
    if (grid_function%cartesian_power(1) .ne. 0) then
      f = f * (px**grid_function%cartesian_power(1))
    end if
    if (grid_function%cartesian_power(2) .ne. 0) then
      f = f * (py**grid_function%cartesian_power(2))
    end if
    if (grid_function%cartesian_power(3) .ne. 0) then
      f = f * (pz**grid_function%cartesian_power(3))
    end if

!   Multiply with radial powers
    if (grid_function%radial_power .ne. 0) then
      r2 = px**2 + py**2 + pz**2
      f  = f * r2**(grid_function%radial_power/2.0d0)
    end if

!   Scale result 
    if (grid_function%up_scale .ne. 1.0d0) then
      f = f * grid_function%up_scale 
    end if

    if (grid_function%down_scale .ne. 1.0d0 .and. &
        grid_function%down_scale .ne. 0.0d0) then
        f = f / grid_function%down_scale
    end if

    grid_function_values(:, ipoint:ipoint) = f*grid_function_values(:, ipoint:ipoint)

  end subroutine

! ------------------------------------------------------------------------------
  subroutine modify_grid_function(grid_function,&
                                  grid_function_values,&
                                  ndim, nr_points,&
                                  grid_points)

    type(grid_function_t), intent(inout):: grid_function
    real(8), intent(inout):: grid_function_values(:,:)
    integer, intent(in):: ndim,nr_points
    real(8), intent(in):: grid_points(ndim,nr_points)
    integer:: ipoint
    real(8):: f, px, py, pz

    f = 1.0d0

!   Multiply with Cartesian powers
    do ipoint = 1, nr_points
      px = grid_points(1,ipoint)
      py = grid_points(2,ipoint)
      pz = grid_points(3,ipoint)

      if (grid_function%cartesian_power(1) .ne. 0) then
        f = f * (px**grid_function%cartesian_power(1))
      end if
      if (grid_function%cartesian_power(2) .ne. 0) then
        f = f * (py**grid_function%cartesian_power(2))
      end if
      if (grid_function%cartesian_power(3) .ne. 0) then
        f = f * (pz**grid_function%cartesian_power(3))
      end if

      grid_function_values(:, ipoint:ipoint) = f*grid_function_values(:, ipoint:ipoint)
    end do

!   Scale result 
    if (grid_function%up_scale .ne. 1.0d0) then
      f = grid_function%up_scale 
      grid_function_values(:, 1:nr_points) = f*grid_function_values(:, 1:nr_points)
    end if
    if (grid_function%down_scale .ne. 1.0d0) then
      f = 1.0d0 / grid_function%down_scale
      grid_function_values(:, 1:nr_points) = f*grid_function_values(:, 1:nr_points)
    end if

  end subroutine

end module
