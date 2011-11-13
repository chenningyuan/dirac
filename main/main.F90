!
!
!      Copyright (c) 2008 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.sdu.dk
!
!
! this file by ulfek march 2008

program dirprg
  use memory_allocator
  use os_utils
  implicit none
#include "priunit.h"
#include "infpar.h"
#include "dcbgen.h"
#ifdef VAR_MPI 
#include "mpif.h"
#endif
  real(kind=8), allocatable :: memory(:)
  integer :: memsize, iparcal, status, length, ios, idummy, ierr, maxsize = -1
  integer :: iscr_mem
#ifdef VAR_MPI 
  call mpixinit()
  call stopwatch(1,numnod+1)
#else
  parcal  = .false.
  iparcal = 0
  mytid   = 0
  mparid  = 0
  numnod  = 0
  call stopwatch(1,numnod)
#endif
  if (parcal) iparcal = 1

  if (mytid == 0) then
     memsize = get_environment_integer('DIRWRK',INSTALL_WRKMEM)
     maxsize = get_environment_integer('DIRMAX',maxsize)
  else
!    memsize = get_environment_integer('DIRNOD',INSTALL_WRKMEM) 
!    SK - Sep 08: DIRNOD is passed to slaves in mpixinit and 
!                 stored on common block GNCPAM (dcbgen.h)
     if( CDIRNOD(1:LENDND) .eq. ' ' )then
       memsize = INSTALL_WRKMEM
     else
       read (CDIRNOD(1:LENDND),*,iostat=ios) iscr_mem
       if(ios/=0)then
         print *,'Warning: Cannot parse environment variable DIRNOD'
         print *,'         Will use the default value:',INSTALL_WRKMEM
       else
         memsize = iscr_mem
       endif
     endif
     maxsize = get_environment_integer('DIRMAX',maxsize)
  endif

#ifdef VAR_MPI
  if (parcal) then
    call mpi_barrier(mpi_comm_world,ierr)
    call allocator_init()
    call allocator_setmax(maxsize)
    call sync_nodes(memsize,0)
    call sync_nodes(maxsize,2)
  else
#else
    call allocator_init()
    call allocator_setmax(maxsize)
    print '(a,i0,a,i0,a)', 'DIRAC serial allocating ',memsize,' words (',memsize/2**17,' MB) of memory'
!    print '(a,i0,a,i0,a)', 'DIRAC serial will attempt to use at most ',  &
!         maxsize,' words (',memsize/2**17,' MB) of memory in dynamic allocations'
#endif
#ifdef VAR_MPI
  endif
#endif
  if (mytid == 0) then
     call alloc(memory,memsize)
     call printtitle
     call printlogo
     call printsubtitle
     call flush(lupri)
     call sysinf(iparcal)
     call flush(lupri)
     call dirac(iparcal,mytid,mparid,numnod,memory,memsize)
  else
#ifdef VAR_MPI
     call alloc(memory,memsize)
     call dirnod(iparcal,mytid,mparid,numnod,memory,memsize)
#else
     call quit('Internal error in main.F90. Quitting.')
#endif
  endif
  if (mytid == 0) then
     print '(a)', '*****************************************************'
  endif
#ifdef VAR_MPI
  if (parcal) then
    call sync_nodes(idummy,1)
  endif
#endif
  if (mytid == 0) then
     call stopwatch(0,mytid)
     print '(a)','********** E N D   of   D I R A C  output  **********'
     print '(a/)', '*****************************************************'
  endif
  call dealloc(memory)
#ifdef VAR_MPI
! if (parcal) call mpixfinalize()
!             SK - 21-07-2008 Since we have called mpixinit explicitly, a call to
!             mpixfinalize is mandatory.
              call mpixfinalize()
#endif
end program dirprg
