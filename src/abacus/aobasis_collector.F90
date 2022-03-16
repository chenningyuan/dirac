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

module aobasis_collector

implicit none


character(len=80), allocatable, public  :: aobasis_container(:)
logical,                        public  :: aobasis_container_reallocate = .true.
integer,                        private :: aobasis_counter = 1

    public add_aobasis2container
    public deallocate_aobasis_container


contains

    subroutine allocate_container()
#include "mxcent.h"

      if(.not.allocated(aobasis_container).and.aobasis_container_reallocate) &
      allocate(aobasis_container(MXATOM));aobasis_container = ""
      aobasis_container_reallocate = .false.
    end subroutine allocate_container

    subroutine deallocate_aobasis_container()
      if(allocated(aobasis_container)) deallocate(aobasis_container)
    end subroutine deallocate_aobasis_container

    subroutine add_aobasis2container(aobasis)
      character(len=80), intent(in) :: aobasis
      if(aobasis_container_reallocate)then
        call allocate_container()
      else
        return
      end if
      aobasis_container(aobasis_counter) = trim(aobasis)
      aobasis_counter = aobasis_counter + 1
    end subroutine add_aobasis2container

end module aobasis_collector
