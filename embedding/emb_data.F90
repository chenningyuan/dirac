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

module fde_data
      
   use memory_allocator
   use fde_cfg

   implicit none

   private

   real(kind=8), save, allocatable :: fde_frozen_n(:)
   real(kind=8), save, allocatable :: fde_frozen_gn(:,:)
   real(kind=8), save, allocatable :: fde_frozen_elpot(:)

   real(kind=8), save, allocatable :: fde_active_n(:)
   real(kind=8), save, allocatable :: fde_active_gn(:,:)
   real(kind=8), save, allocatable :: fde_active_elpot(:)

   real(kind=8), save, allocatable :: fde_vemb(:)

!  public fde_importStaticVemb, fde_cleanupStaticVemb
!  public fde_importFrozenVars, fde_cleanupFrozenVars

!!!public fde_setFrozenN, fde_setFrozenGN, fde_setFrozenEP
!!!public fde_setActiveEP
!  public fde_setActiveN, fde_setActiveGN

   real(kind=8), pointer :: fde_grid(:,:), fde_frozen_prp(:,:) 

   public fde_frozen_n, fde_frozen_gn, fde_frozen_elpot
   public fde_active_n, fde_active_gn, fde_active_elpot
   public fde_vemb
      
   contains

      subroutine get_fde_grid_coords(i,x,y,z)
         integer, intent(in) :: i
         real(kind=8), intent (out) :: x, y, z

         x = fde_grid(1,i)
         y = fde_grid(2,i)
         z = fde_grid(3,i)

      end subroutine get_fde_grid_coords

      subroutine set_pointer_fde_grid(p) 
         real(kind=8), pointer :: p(:,:)
         p => fde_grid
      end subroutine set_pointer_fde_grid

!!!!!!!
!  subroutine fde_data_set_active_n(i,n)
!     integer, intent(in) :: i
!     real(kind=8), intent(in) : n, gnx,gny,gnz
!  end fde_data_set_active_n
!
!  subroutine fde_data_set_active_ngn(i,n,gnx,gny,gnz)
!     integer, intent(in) :: i
!     real(kind=8), intent(in) : n, gnx,gny,gnz
!
!     call fde_data_set_active_n(i,n)
!     call fde_data_set_active_gn(i,gnx,gny,gnz)
!
!     fde_active_n(i) = n
!     fde_active_gn(1,i) = gn(1)
!     fde_active_gn(2,i) = gn(2)
!     fde_active_gn(3,i) = gn(3)
!
!  end subroutine fde_data_set_active_ngn
!!!!!!!

   subroutine import_static
   end subroutine import_static

   subroutine import_frozen
   end subroutine import_frozen

   subroutine cleanup_static
   end subroutine cleanup_static

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_setupFDEgrid: opens the file containing the frozen density, reads in the density
!                   and whatever other quantities we have exportend (such as the 
!                   components of the density gradient etc) using read_grid (internally,
!                   it is read_grid which will allocate storage, and later we have to 
!                   free it using FDE_cleanupFDEgrid for instance). for convenience, 
!                   we also set fde_grid_np (number of grid points) and fde_nprop 
!                   (number of properties read in) 
!
         subroutine FDE_setupFDEgrid()
            integer :: frozdens_unit
            integer :: i
            real(kind=8), pointer :: fde_frozen_prp(:,:)

            call FDE_ImportFrozenDensity_FileOpen(frozdens_unit)
            call read_grid(frozdens_unit,fde_grid,fde_frozen_prp)
            call FDE_ImportFrozenDensity_FileClose()

            fde_grid_np = size(fde_grid,2)

            call alloc(fde_frozen_n,     fde_grid_np)
            call alloc(fde_frozen_gn,    fde_grid_np,3)
            call alloc(fde_frozen_elpot, fde_grid_np)

            do i=1, fde_grid_np
               fde_frozen_elpot(i) = fde_frozen_prp(i,1) + fde_frozen_prp(i,2) 
               fde_frozen_n(i)     = fde_frozen_prp(i,3)
               fde_frozen_gn(i,1)  = fde_frozen_prp(i,4)
               fde_frozen_gn(i,5)  = fde_frozen_prp(i,5)
               fde_frozen_gn(i,6)  = fde_frozen_prp(i,6)
            enddo
            
            call deallocp(fde_frozen_prp)

         end subroutine FDE_setupFDEgrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_testFrozenDensity: integrates the imported density, useful to check if the
!                        number of electrons does match that of the frozen subsystems
!
         subroutine FDE_testFrozenDensity()
            integer      :: point = 0
            real(kind=8) :: frozen_particle_number = 0.0

            do point = 1, fde_grid_np 
               frozen_particle_number = frozen_particle_number + fde_grid(4,point)*fde_frozen_prp(3,point) 
            end do
            
            write (*,*) '  FDE info: imported frozen density integrates to ',frozen_particle_number,' electrons'
            write (*,*) '            grid points processed ',fde_grid_np
            write (*,*) '            number of properties  ',fde_nprop
            frozen_particle_number = 0.0

         end subroutine FDE_testFrozenDensity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FDE_cleanupFDEgrid: dealocates the storage used for holding the grid, and whatever
!                     properties exported at the grid points
!
         subroutine FDE_cleanupFDEgrid()
            deallocate(fde_grid)
            deallocate(fde_frozen_prp) 
         end subroutine FDE_cleanupFDEgrid



end module fde_data
