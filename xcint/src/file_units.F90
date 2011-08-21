module file_units

   implicit none

   character(*), parameter, public :: interface_file_name_ao      = 'interface_ao'
   character(*), parameter, public :: interface_file_name_mo      = 'interface_mo'
   character(*), parameter, public :: interface_file_name_mo_coef = 'interface_mo_coef'
   integer,      parameter, public :: interface_file_unit = 137
   logical                , public :: interface_file_exists

   private

end module
