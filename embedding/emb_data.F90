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
   use fde_types
   use fde_cfg
   use fde_io

   implicit none

   private
!
! routines
!
   public fde_data_cleanup
   public get_fde_grid_coord
   public get_fde_grid_weight
   public set_pointer_fde_grid
   public fde_test_frozen_density

   interface get_fde_grid_coord
      module procedure get_fde_grid_coords_single
      module procedure get_fde_grid_coords_block
   end interface get_fde_grid_coord

   interface get_fde_grid_weight
      module procedure get_fde_grid_weight_single
      module procedure get_fde_grid_weight_block
   end interface get_fde_grid_weight

!!!public fde_setFrozenN, fde_setFrozenGN, fde_setFrozenEP
!!!public fde_setActiveEP
!  public fde_setActiveN, fde_setActiveGN

!
! variables 
!
   public fde_frozen_n, fde_frozen_gn, fde_frozen_elpot
   public fde_active_n, fde_active_gn, fde_active_elpot
   public fde_vemb
   public fde_grid_np

   integer, save :: fde_grid_np = 0

   real(kind=8), save, allocatable :: fde_frozen_n(:)
   real(kind=8), save, allocatable :: fde_frozen_gn(:,:)
   real(kind=8), save, allocatable :: fde_frozen_elpot(:)

   real(kind=8), save, allocatable :: fde_active_n(:)
   real(kind=8), save, allocatable :: fde_active_gn(:,:)
   real(kind=8), save, allocatable :: fde_active_elpot(:)

   real(kind=8), save, allocatable :: fde_vemb(:)

   real(kind=8), pointer :: fde_grid(:,:) 

   contains

      subroutine get_fde_grid_coords_single(i,x,y,z)
         integer, intent(in) :: i
         real(kind=8), intent (out) :: x, y, z

         x = fde_grid(1,i)
         y = fde_grid(2,i)
         z = fde_grid(3,i)
      end subroutine get_fde_grid_coords_single

      subroutine get_fde_grid_coords_block(i,stride,x,y,z)
         integer, intent(in) :: i, stride
         integer :: j, k 
         real(kind=8), intent (out) :: x(:), y(:), z(:)

         k = 1
         do j = i, i+stride
            x(k) = fde_grid(1,j)
            y(k) = fde_grid(2,j)
            z(k) = fde_grid(3,j)
            k = k+1
         enddo
      end subroutine get_fde_grid_coords_block


      subroutine get_fde_grid_weight_single(i,w)
         integer, intent(in) :: i
         real(kind=8), intent (out) :: w

         w = fde_grid(4,i)
      end subroutine get_fde_grid_weight_single

      subroutine get_fde_grid_weight_block(i,stride,w)
         integer, intent(in) :: i, stride
         integer :: j, k
         real(kind=8), intent (out) :: w(:)

         k = 1
         do j = i, i+stride
            w(k) = fde_grid(4,j)
            k = k+1
         enddo
      end subroutine get_fde_grid_weight_block


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

   subroutine fde_import_static
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i, fde_grid_np
      real(kind=8), pointer :: fde_vtmp(:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%embpot%unit
      file_name = ftmp%embpot%name
!      call fde_get_embpot_unit(file_unit)
!      call fde_get_embpot_filename(file_name)

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,fde_grid,fde_vtmp)
      call fde_close_file(file_unit)

      fde_grid_np = size(fde_grid,2)

      call alloc(fde_vemb, fde_grid_np)

      do i=1, fde_grid_np
         fde_vemb(i) = fde_vtmp(i)
      enddo

      call deallocp(fde_vtmp)
   end subroutine fde_import_static


   subroutine fde_import_frozen
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i, fde_grid_np
      real(kind=8), pointer :: fde_frozen_prp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%frozen%unit
      file_name = ftmp%frozen%name

!      call fde_get_frozen_unit(file_unit)
!      call fde_get_frozen_filename(file_name)

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,fde_grid,fde_frozen_prp)
      call fde_close_file(file_unit)

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
   end subroutine fde_import_frozen


   subroutine fde_data_cleanup

      call dealloc(fde_frozen_n)
      call dealloc(fde_frozen_gn)
      call dealloc(fde_frozen_elpot)

      call dealloc(fde_active_n)
      call dealloc(fde_active_gn)
      call dealloc(fde_active_elpot)

      call dealloc(fde_vemb)

      call deallocp(fde_grid)
   end subroutine fde_data_cleanup


   subroutine fde_test_frozen_density
      integer      :: point = 0
      real(kind=8) :: frozen_particle_number = 0.0

      do point = 1, fde_grid_np 
         frozen_particle_number = frozen_particle_number + fde_grid(4,point)*fde_frozen_n(point) 
      end do
      
      write (*,*) '  FDE info: imported frozen density integrates to ',frozen_particle_number,' electrons'
      write (*,*) '            grid points processed ',fde_grid_np

      frozen_particle_number = 0.0

   end subroutine fde_test_frozen_density

end module fde_data
