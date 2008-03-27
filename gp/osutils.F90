! This module contains functions to communicate with the operating system.
! ulfek march 2008
module os_utils
       contains
! Gets an environment variable and converts it to an integer.
! Return default if the variable is not defined or conversion fails
function get_environment_integer(name, default)
  implicit none
  integer :: default, get_environment_integer, length, status, i, ios
  character(len=*) :: name
  character(len=128) :: value
#ifndef FORTRAN2003  
  integer :: getenv
  status = getenv(name, value)
! status == 0: mismatch; status == 1: succes;
  if(status==0) then
     status=1
  else
     status=0
  end if
#else
  call get_environment_variable(name, value, length, status)
#endif  
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
end module os_utils
