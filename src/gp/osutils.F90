!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

! This module contains functions to communicate with the operating system.
! ulfek march 2008
module os_utils

       private

       public get_environment_integer
    
       interface get_environment_integer
          module procedure get_environment_integer_i4
          module procedure get_environment_integer_i8
       end interface get_environment_integer

       contains
! Gets an environment variable and converts it to an integer.
! Return default if the variable is not defined or conversion fails
function get_environment_integer_i4(name, default)
  implicit none
  integer(kind=4) :: default, get_environment_integer_i4, length, status, i, ios
  character(len=*) :: name
  character(len=128) :: value
#ifndef FORTRAN2003  
  call getenv(name, value)
  if (len_trim(value) == 0) then
     status = 1
  else  
     status = 0 
  endif
#else
  call get_environment_variable(name, value, length, status)
#endif  
  if (status==0) then
     read (value,*,iostat=ios) i
     if (ios/=0) then
        print *,'ERROR reading environment variable ',name
        print *,'ERROR cannot convert to 32 bit integer, value = ',value
        call quit('Invalid environment variable value')
        get_environment_integer_i4 = default
     else
        get_environment_integer_i4 = i
     endif
  else
        get_environment_integer_i4 = default
  endif
end function get_environment_integer_i4


function get_environment_integer_i8(name, default)
  implicit none
  integer(kind=8) :: default, get_environment_integer_i8, i
  integer(kind=4) :: length, status, ios
  character(len=*) :: name
  character(len=128) :: value
#ifndef FORTRAN2003  
  call getenv(name, value)
  if (len_trim(value) == 0) then
     status = 1
  else
     status = 0
  endif
#else
  call get_environment_variable(name, value, length, status)
#endif  
  if (status==0) then
     read (value,*,iostat=ios) i
     if (ios/=0) then
        print *,'ERROR reading environment variable ',name
        print *,'ERROR cannot convert to 64 bit integer, value = ',value
        call quit('Invalid environment variable value')
        get_environment_integer_i8 = default
     else
        get_environment_integer_i8 = i
     endif
  else
        get_environment_integer_i8 = default
  endif
end function get_environment_integer_i8


end module os_utils
