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
   public fde_cleanup_data
   public get_fde_grid_coord
   public get_fde_grid_weight
   public fde_test_frozen_density
   public fde_test_active_density
   public fde_import_frozen
   public fde_import_static
   public fde_initialize_data

   interface get_fde_grid_coord
      module procedure get_fde_grid_coords_single
      module procedure get_fde_grid_coords_block
   end interface get_fde_grid_coord

   interface get_fde_grid_weight
      module procedure get_fde_grid_weight_single
      module procedure get_fde_grid_weight_block
   end interface get_fde_grid_weight


!
! variables 
!
   public fde_frozen_n, fde_frozen_gn, fde_frozen_elpot
   public fde_active_n, fde_active_gn, fde_active_elpot
   public fde_grid
   public fde_vemb
   public fde_grid_np

   integer, save :: fde_grid_np = 0
   logical, save :: data_is_initialized = .false.

   real(kind=8), save, allocatable, target :: fde_frozen_n(:)
   real(kind=8), save, allocatable, target :: fde_frozen_gn(:,:)
   real(kind=8), save, allocatable, target :: fde_frozen_elpot(:)

   real(kind=8), save, allocatable, target :: fde_active_n(:)
   real(kind=8), save, allocatable, target :: fde_active_gn(:,:)
   real(kind=8), save, allocatable, target :: fde_active_elpot(:)

   real(kind=8), save, allocatable, target :: fde_vemb(:)

   real(kind=8), save, allocatable, target :: fde_grid(:,:) 

   contains

!   ----------------------------------------------------------------------------
      subroutine get_fde_grid_coords_single(i,x,y,z)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: i
         real(kind=8), intent (out) :: x, y, z

         x = fde_grid(1,i)
         y = fde_grid(2,i)
         z = fde_grid(3,i)
      end subroutine get_fde_grid_coords_single


!   ----------------------------------------------------------------------------
      subroutine get_fde_grid_coords_block(first,last,x,y,z)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: first, last
         integer :: j, k 
         real(kind=8), intent (out) :: x(:), y(:), z(:)

         k = 1
         do j = first, first+last-1
            x(k) = fde_grid(1,j)
            y(k) = fde_grid(2,j)
            z(k) = fde_grid(3,j)
            k = k+1
         enddo
      end subroutine get_fde_grid_coords_block


!   ----------------------------------------------------------------------------
      subroutine get_fde_grid_weight_single(i,w)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: i
         real(kind=8), intent (out) :: w

         w = fde_grid(4,i)
      end subroutine get_fde_grid_weight_single


!   ----------------------------------------------------------------------------
      subroutine get_fde_grid_weight_block(first,last,w)
!   ----------------------------------------------------------------------------
         integer, intent(in) :: first, last
         integer :: j, k
         real(kind=8), intent (out) :: w(:)

         k = 1
         do j = first, first+last-1
            w(k) = fde_grid(4,j)
            k = k+1
         enddo
      end subroutine get_fde_grid_weight_block


!   ----------------------------------------------------------------------------
   subroutine fde_import_static
!   ----------------------------------------------------------------------------
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: fde_vtmp(:), gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%embpot%unit
      file_name = ftmp%embpot%name

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,gridp,fde_vtmp)
      call fde_close_file(file_unit)

      fde_grid_np = size(gridp,2)

      call fde_initialize_data

      do i=1, fde_grid_np
         fde_vemb(i)   = fde_vtmp(i)
         fde_grid(1,i) = gridp(1,i)
         fde_grid(2,i) = gridp(2,i)
         fde_grid(3,i) = gridp(3,i)
         fde_grid(4,i) = gridp(4,i)
      enddo

      call deallocp(fde_vtmp)
      call deallocp(gridp)
   end subroutine fde_import_static


!   ----------------------------------------------------------------------------
   subroutine fde_import_frozen
!   ----------------------------------------------------------------------------
      integer :: file_unit
      character(len=60) :: file_name
      integer :: i
      real(kind=8), pointer :: fde_frozen_prp(:,:), gridp(:,:)
      type(fde_files)  :: ftmp

      call fde_get_files_info(ftmp)
      
      file_unit = ftmp%frozen%unit
      file_name = ftmp%frozen%name

      call fde_open_file(file_name,file_unit)
      call read_grid(file_unit,gridp,fde_frozen_prp)
      call fde_close_file(file_unit)

      fde_grid_np = size(gridp,2)

      call fde_initialize_data
      
      do i=1, fde_grid_np
         fde_frozen_elpot(i) = fde_frozen_prp(1,i) + fde_frozen_prp(2,i)
         fde_frozen_n(i)     = fde_frozen_prp(3,i)
         fde_frozen_gn(i,1)  = fde_frozen_prp(4,i)
         fde_frozen_gn(i,2)  = fde_frozen_prp(5,i)
         fde_frozen_gn(i,3)  = fde_frozen_prp(6,i)
         
         fde_grid(1,i)       = gridp(1,i)
         fde_grid(2,i)       = gridp(2,i)
         fde_grid(3,i)       = gridp(3,i)
         fde_grid(4,i)       = gridp(4,i)
      enddo

      call deallocp(fde_frozen_prp)
      call deallocp(gridp)
   end subroutine fde_import_frozen


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_data
!   ----------------------------------------------------------------------------
      if (.not.data_is_initialized) then
         data_is_initialized = .true.
         if (fde_grid_np > 0) then
            
            if (.not.allocated(fde_grid)) then
               call alloc(fde_grid, 4,fde_grid_np, id="fde_grid")
            endif
            call fde_initialize_frozen_data
            call fde_initialize_active_data
         else
            call quit('there are no grid points for fde!')
         endif
      endif
      
   end subroutine fde_initialize_data


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_frozen_data
!   ----------------------------------------------------------------------------
         if (.not.allocated(fde_vemb)) then
            call alloc(fde_vemb, fde_grid_np, id="fde_vemb")
            fde_vemb         = 0.0d0
         endif
         
         if (.not.allocated(fde_frozen_n)) then
            call alloc(fde_frozen_n,     fde_grid_np, id="fde_frozen_n")
            fde_frozen_n     = 0.0d0
         endif
         
         if (.not.allocated(fde_frozen_gn)) then
            call alloc(fde_frozen_gn,    fde_grid_np,3, id="fde_frozen_gn")
            fde_frozen_gn    = 0.0d0
         endif
         
         if (.not.allocated(fde_frozen_elpot)) then
            call alloc(fde_frozen_elpot, fde_grid_np, id="fde_frozen_elpot")
            fde_frozen_elpot = 0.0d0
         endif
   end subroutine fde_initialize_frozen_data


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_active_data
!   ----------------------------------------------------------------------------
      if (.not.allocated(fde_active_n)) then
         write (*,*) 'debug, allocating fde_active_n(',fde_grid_np,')'
         call alloc(fde_active_n,     fde_grid_np, id="fde_active_n")
         fde_active_n     = 0.0d0
      endif

      if (.not.allocated(fde_active_gn)) then
         write (*,*) 'debug, allocating fde_active_gn(',fde_grid_np,',3)'
         call alloc(fde_active_gn,    fde_grid_np,3, id="fde_active_gn")
         fde_active_gn    = 0.0d0
      endif

      if (.not.allocated(fde_active_elpot)) then
         write (*,*) 'debug, allocating fde_active_elpot(',fde_grid_np,')'
         call alloc(fde_active_elpot, fde_grid_np, id="fde_active_elpot")
         fde_active_elpot = 0.0d0
      endif

   end subroutine fde_initialize_active_data
   

!   ----------------------------------------------------------------------------
   subroutine fde_cleanup_data
!   ----------------------------------------------------------------------------
      if (data_is_initialized) then
         write (*,*) 'debug, cleaning up internal fde data'
         call dealloc(fde_grid)
         call fde_cleanup_active_data
         call fde_cleanup_frozen_data
      endif
   end subroutine fde_cleanup_data


!   ----------------------------------------------------------------------------
   subroutine fde_cleanup_frozen_data
!   ----------------------------------------------------------------------------
      call dealloc(fde_vemb)
      call dealloc(fde_frozen_n)
      call dealloc(fde_frozen_gn)
      call dealloc(fde_frozen_elpot)
   end subroutine fde_cleanup_frozen_data

   
!   ----------------------------------------------------------------------------
   subroutine fde_cleanup_active_data
!   ----------------------------------------------------------------------------
      call dealloc(fde_active_n)
      call dealloc(fde_active_gn)
      call dealloc(fde_active_elpot)
   end subroutine fde_cleanup_active_data


!   ----------------------------------------------------------------------------
   subroutine fde_test_frozen_density
!   ----------------------------------------------------------------------------
      integer      :: point = 0
      real(kind=8) :: frozen_particle_number = 0.0

      do point = 1, fde_grid_np 
         frozen_particle_number = frozen_particle_number + fde_grid(4,point)*fde_frozen_n(point) 
      end do
      
      write (*,*) '  Imported (frozen) density integrates to ',frozen_particle_number,' electrons'
      write (*,*) '  Grid points processed ',fde_grid_np

      frozen_particle_number = 0.0

   end subroutine fde_test_frozen_density


!   ----------------------------------------------------------------------------
   subroutine fde_test_active_density
!   ----------------------------------------------------------------------------
      integer      :: point = 0
      real(kind=8) :: active_particle_number = 0.0

      do point = 1, fde_grid_np 
         active_particle_number = active_particle_number + fde_grid(4,point)*fde_active_n(point) 
      end do
      
      write (*,*) '  Exported (active) density integrates to ',active_particle_number,' electrons'
      write (*,*) '  Grid points processed ',fde_grid_np

      active_particle_number = 0.0

   end subroutine fde_test_active_density

   
end module fde_data
