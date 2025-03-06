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
!
! purpose: determine integer size in MPI library to 
!
module integer_model_in_mpi
#ifdef VAR_MPI
  use integer_kind_mpilib
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

  public determine_integer_model_in_mpi
  public get_integer_model_in_mpi
  
  private
  save

contains

  subroutine determine_integer_model_in_mpi

#ifdef VAR_MPI_32BIT_INT
  integer_kind_in_mpi = 4
#else
  integer_kind_in_mpi = sizeof(MPI_INTEGER)
#endif

  end subroutine determine_integer_model_in_mpi

  function get_integer_model_in_mpi() result(int_mpi)
    integer :: int_mpi
    int_mpi = integer_kind_in_mpi
  end function get_integer_model_in_mpi

#else
  implicit none
  integer         :: blabla_integer_model_in_mpi
#endif
end module integer_model_in_mpi
