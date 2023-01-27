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
! purpose: store default integer size in Dirac and MPI library
!
module integer_kind_mpilib

#ifdef VAR_MPI

#if defined (INT_STAR8)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

#endif

#ifdef INT_STAR8
  integer, parameter, public :: integer_kind        = 8
#else
  integer, parameter, public :: integer_kind        = 4
#endif
  integer, public            :: integer_kind_in_mpi = 4

!
! next variable is here so it is available both in f77 and f90 interface routines /hjaaj
!
  logical, public            :: MPI_INIT_called     = .false.

end module integer_kind_mpilib
