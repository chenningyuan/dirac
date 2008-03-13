! This module contains functions to communicate with the operating system.
! ulfek march 2008
module osutils
       contains
! Gets an environment variable and converts it to an integer.
! Return default if the variable is not defined or conversion fails
function get_environment_integer(name, default)
  implicit none
  integer :: default, get_environment_integer, length, status, i, ios
  character(len=*) :: name
  character(len=128) :: value
  call get_environment_variable('DIRWRK', value, length, status)
  if (status==0) then
     read (value,*,iostat=ios) i
     if (ios/=0) then
        print *,'Warning: Cannot parse environment variable ',name
        print *,'         Will use the default value:',default
        get_environment_integer = default
     else
        get_environment_integer = i
     endif
  else
        get_environment_integer = default
  endif
end function get_environment_integer
end module osutils