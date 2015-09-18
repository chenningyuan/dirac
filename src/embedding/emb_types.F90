

module fde_types

   use fde_mag_cfg

   implicit none

   public

   type fde_grid
      integer :: id
      integer :: npoints
      real(kind=8), pointer :: r(:,:)
      real(kind=8), pointer :: w(:) 
   end type

   type grid_function 
      integer :: id
      integer :: npoints
      integer :: nspin
      real(kind=8), pointer :: n(:)
      real(kind=8), pointer :: n_b(:,:)
      real(kind=8), pointer :: s_b(:,:,:)
      real(kind=8), pointer :: gn(:,:)
      real(kind=8), pointer :: gn_b(:,:,:)
      real(kind=8), pointer :: gs_b(:,:,:,:)
      real(kind=8), pointer :: elpot(:)
      real(kind=8), pointer :: nucpot(:)
      real(kind=8), pointer :: hn(:,:)
   end type

   public new_fde_grid
   public del_fde_grid
!gosia: unused   public set_fde_grid_size

   public new_grid_function
   public del_grid_function
   public set_grid_function_nspin

   contains


      subroutine new_grid_function(gf,np)
         type(grid_function) :: gf
         integer :: np

         gf%npoints = np
         gf%nspin   = 1

         allocate(gf%n(np))
         allocate(gf%gn(3,np))
         allocate(gf%hn(6,np))
         allocate(gf%elpot(np))
         allocate(gf%nucpot(np))

         if (fde_rsp_mag_lao_import .or. fde_rsp_mag_lao_export) then
            allocate(gf%n_b(3,np))
            allocate(gf%gn_b(3,3,np))
            if (.not. fde_rsp_mag_lao_no_s) then
               allocate(gf%s_b(3,3,np))
               allocate(gf%gs_b(3,3,3,np))
            end if
         end if
         
         gf%n  = 0.0d0
         gf%gn = 0.0d0
         gf%hn = 0.0d0
         gf%elpot = 0.0d0
         gf%nucpot = 0.0d0

         if (fde_rsp_mag_lao_import .or. fde_rsp_mag_lao_export) then
            gf%n_b  = 0.0d0
            gf%gn_b = 0.0d0
            if (.not. fde_rsp_mag_lao_no_s) then
               gf%s_b  = 0.0d0
               gf%gs_b = 0.0d0
            end if
         end if

      end subroutine new_grid_function


      subroutine del_grid_function(gf)
         type(grid_function) :: gf

         gf%npoints = 0
         deallocate(gf%n)
         deallocate(gf%gn)
         deallocate(gf%hn)
         deallocate(gf%elpot)
         deallocate(gf%nucpot)

!gosia: do rather "if allocated"
         if (fde_rsp_mag_lao_import .or. fde_rsp_mag_lao_export) then
            deallocate(gf%n_b)
            deallocate(gf%gn_b)
            if (.not. fde_rsp_mag_lao_no_s) then
               deallocate(gf%s_b)
               deallocate(gf%gs_b)
            end if
         end if

      end subroutine del_grid_function


      subroutine set_grid_function_nspin(gf,nspin)
         type(grid_function) :: gf
         integer :: nspin 

         gf%nspin = nspin
      end subroutine set_grid_function_nspin
      

     subroutine new_fde_grid(g,np)
       type(fde_grid) :: g
       integer :: np

       g%npoints = np
       allocate(g%r(3,np))
       allocate(g%w(np))
     end subroutine new_fde_grid


     subroutine del_fde_grid(g)
        type(fde_grid) :: g

        g%npoints = 0
        deallocate(g%r)
        deallocate(g%w)
     end subroutine del_fde_grid

end
