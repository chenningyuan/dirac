

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


end module fde_types
