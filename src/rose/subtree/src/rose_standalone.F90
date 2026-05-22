program rose_standalone

   use rose_utils
   use rose_global

   implicit none

   logical :: input_found = .false.

   call get_free_fileunit(luinp)
   call read_menu_input('INPUT_GENIBO', luinp, 'ALL', input_found)
   if (.not. input_found) error stop 'INPUT_GENIBO not found'

   call rose

end program rose_standalone
