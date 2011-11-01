

module fde_types

   public

   type fde_export
      character(len=4) :: ex_level
      character(len=3) :: ex_format
      logical          :: ex_n
      logical          :: ex_gn
      logical          :: ex_coulomb
      logical          :: do_grid_out
   end type
   
   type fde_import
      character(len=4) :: im_level
      character(len=3) :: im_format
      logical          :: im_n
      logical          :: im_gn
      logical          :: im_coulomb
      logical          :: im_vemb
      logical          :: im_frozen
      logical          :: im_update_vemb
   end type

   type fde_file
      character(len=60) :: name
      integer           :: unit
   end type
   
   type fde_files
      type(fde_file) :: embpot
      type(fde_file) :: frozen
      type(fde_file) :: export
      type(fde_file) :: test
      type(fde_file) :: logfile
   end type

   type grid_function 
      integer :: npoints
      integer :: nspin
      real(kind=8), allocatable :: n(:)
      real(kind=8), allocatable :: gn(:,:)
      real(kind=8), allocatable :: elpot(:)
   end type

   public new_grid_function
   public del_grid_function
   public set_grid_function_nspin

   contains


      subroutine new_grid_function(gf,size)
         type(grid_function) :: gf
         integer :: size

         gf%npoints = size
         gf%nspin   = 1
         allocate(gf%n(size))
         allocate(gf%gn(3,size))
         allocate(gf%elpot(size))
         
         gf%n  = 0.0d0
         gf%gn = 0.0d0
         gf%elpot = 0.0d0
      end subroutine new_grid_function


      subroutine del_grid_function(gf)
         type(grid_function) :: gf
         integer :: size

         gf%npoints = 0
         deallocate(gf%n)
         deallocate(gf%gn)
         deallocate(gf%elpot)
      end subroutine del_grid_function


      subroutine set_grid_function_nspin(gf,nspin)
         type(grid_function) :: gf
         integer :: nspin 

         gf%nspin = nspin
      end subroutine set_grid_function_nspin
      
end module fde_types
