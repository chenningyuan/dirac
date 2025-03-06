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
!     MPI compatibility test
!
#ifdef VAR_MPI
      program selftest_mpi
!     Call MPI routines with some small tests to see if things look ok
!     implicit real*8 (A-H,O-Z)
      implicit none
#include "mpif.h"
      integer            :: mpierr, ireq
      integer            :: istat(mpi_status_size)
      integer            :: myrank, whatever
      integer            :: my_MPI_INTEGER
      integer, parameter :: master_id = 0

#ifdef INT_STAR8
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

      call mpi_init(mpierr)
      call mpi_comm_rank(mpi_comm_world, myrank, mpierr)
!
!     send dummy message (from MASTER to MASTER...)
!
      call mpi_isend(myrank,1,my_MPI_INTEGER,0,25,               &
                     mpi_comm_world,ireq,mpierr)
      call mpi_probe(mpi_any_source,mpi_any_tag,mpi_comm_world,  &
                     istat,mpierr)
      if( istat(mpi_source) .ne. master_id )then
         print *,' *** ERROR: MPI self test failed! *** '
         print *,' sender of message should have ID 0'
         print *,' but has ID ',istat(mpi_source)
         print *,' '
         print *,' Make sure that you are using the right mpirun for'
         print *,' your particular MPI library version.'
         print *,' '
#ifdef INT_STAR8
         print *,' DIRAC is configured to use integer*8 '        &
              //'(64 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*8, make sure that this is the case.'
#else
         print *,' DIRAC is configured to use integer*4 '        &
              //' (32 bit integers).'
         print *,' You need to link a MPI library that is compiled'
         print *,' using integer*4, make sure that this is the case.'
#endif
         print *,' '
         print *,' Cannot continue without a working MPI library.'
         print *,' Quitting.'
         stop ' *** ERROR: MPI self test failed ***'
      end if
!
!     ... it looks ok, receive dummy message
!
      call mpi_irecv(whatever,1,my_MPI_INTEGER,0,mpi_any_tag,    &
                     mpi_comm_world,ireq,mpierr)
      call mpi_wait(ireq,istat,mpierr)

      call mpi_finalize(mpierr)

      end 
#endif
