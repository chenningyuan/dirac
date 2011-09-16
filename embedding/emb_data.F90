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
   use interface_ao
   use interface_mo
   use ao_eval
   use density_eval
   use xc_max_block_length


   implicit none

   private
!
! routines
!
   public fde_cleanup_data
   public get_fde_grid_coord
   public get_fde_grid_weight
   public fde_test_frozen_density
   public fde_import_frozen
   public fde_import_static
   public fde_get_active_density
   

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

      gridp => fde_grid
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
            call alloc(fde_grid, 4,fde_grid_np, id="fde_grid")
         
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
         call alloc(fde_vemb, fde_grid_np, id="fde_vemb")
         call alloc(fde_frozen_n,     fde_grid_np, id="fde_frozen_n")
         call alloc(fde_frozen_gn,    fde_grid_np,3, id="fde_frozen_gn")
         call alloc(fde_frozen_elpot, fde_grid_np, id="fde_frozen_elpot")
         
         fde_vemb         = 0.0d0  
         fde_frozen_n     = 0.0d0
         fde_frozen_gn    = 0.0d0
         fde_frozen_elpot = 0.0d0
   end subroutine fde_initialize_frozen_data


!   ----------------------------------------------------------------------------
   subroutine fde_initialize_active_data
!   ----------------------------------------------------------------------------
         call alloc(fde_active_n,     fde_grid_np, id="fde_active_n")
         call alloc(fde_active_gn,    fde_grid_np,3, id="fde_active_gn")
         call alloc(fde_active_elpot, fde_grid_np, id="fde_active_elpot")
         
         fde_active_n     = 0.0d0
         fde_active_gn    = 0.0d0
         fde_active_elpot = 0.0d0
   end subroutine fde_initialize_active_data
   

!   ----------------------------------------------------------------------------
   subroutine fde_cleanup_data
!   ----------------------------------------------------------------------------
      if (data_is_initialized) then
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
      
      write (*,*) '  Imported frozen density integrates to ',frozen_particle_number,' electrons'
      write (*,*) '  Grid points processed ',fde_grid_np

      frozen_particle_number = 0.0

   end subroutine fde_test_frozen_density


!   ----------------------------------------------------------------------------
   subroutine fde_get_active_density(dmat, mat_dim, get_gradient)
!   ----------------------------------------------------------------------------
      real(kind=8), intent(in) :: dmat(*)
      integer, intent(in)      :: mat_dim
      logical, intent(in) :: get_gradient
 
      real(kind=8) :: n_0(max_block_length)
      real(kind=8) :: gn_0_pointwise(3)
      real(8), allocatable :: ao(:, :)
      real(8), allocatable :: buffer(:, :)
      integer :: ipoint
      integer :: i
      integer :: nr_ao
      integer :: max_ao_g_order, max_ao_m_order


      max_ao_g_order = 0
      if (get_gradient) then
        max_ao_g_order = 1
      end if
     
      max_ao_m_order = 0
     
      call ao_eval_init(max_ao_g_order, &
                        max_ao_m_order, &
                        .true.)
                        
      nr_ao = nr_ao_slices*nr_ao_cartesian
      
      call alloc(ao,     1, nr_ao, id="ao")
      call alloc(buffer, 1, nr_ao, id="buffer")

      ao     = 0.0d0
      buffer = 0.0d0


      fde_active_n  = 0.0d0
      fde_active_gn = 0.0d0
      
      do ipoint = 1, fde_grid_np
         
!         time_ao_start = second()
         call get_ao(1,          &
                     fde_grid(ipoint,1), &
                     fde_grid(ipoint,2), &
                     fde_grid(ipoint,3), &
                     ao,         &
                     buffer)
                     
!         time_ao = time_ao + second() - time_ao_start

!         time_density_start = second()
         
         if (need_ao_order(1, 0)) then
!        density and density gradient
            call get_gn(n_0(1), gn_0_pointwise, 0, mat_dim, dmat, buffer, ao)
         else
!        density
            call get_n(n_0(1), 0, mat_dim, dmat, buffer, ao)
!           for mcscf/srdft also evaluate spin density here
!          call get_s(spin_density, 0, mat_dim, dmat_0, buffer, ao)
         end if

         if ( n_0(1) > fde_cfg_screening ) then
            fde_active_n(ipoint) = n_0(1)
            if (get_gradient) then
               fde_active_gn(1,ipoint) = gn_0_pointwise(1)
               fde_active_gn(2,ipoint) = gn_0_pointwise(2)
               fde_active_gn(3,ipoint) = gn_0_pointwise(3)
            endif                     
         endif

!         time_density = time_density + second() - time_density_start

      enddo
      call dealloc(ao)
      call dealloc(buffer)
   end subroutine fde_get_active_density
   
   
end module fde_data
