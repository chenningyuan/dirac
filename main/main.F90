!dirac_copyright_start
!      Copyright (c) 2011 by the authors of DIRAC.
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
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

! this file by ulfek march 2008

#if defined FORTRAN_MAIN_AS_SUBROUTINE
! ... needed to get another Fortran executable 
subroutine dirprg
#else
! default executable - dirac.x
program dirprg
#endif

  use memory_allocator
  use os_utils
  use dft_cfg
  use xmlout
  use interface_ao
  use interface_mo

  implicit none
#include "priunit.h"
#include "infpar.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#ifdef VAR_MPI
#include "mpif.h"
#endif
  real(kind=8), allocatable :: memory(:)
  integer :: memsize, iparcal, status, length, ios, idummy, ierr
  integer(kind=8) :: maxsize = -1
  integer :: iscr_mem
  integer :: high_water_mark
! initialize common variables (parallel and sequential run)
  mytid   = 0
  mparid  = 0
  numnod  = 0
  parcal  = .false.
  iparex  = 0
#ifdef VAR_MPI
! check for sequential run with parallel executable 
! iparex = 1 for parallel run mode with MPI executable (default)
  iparex = get_environment_integer('DIRPAR',1)
!
! initialize MPI environment and wake up slaves - set parcal == .true.
  if (iparex .gt. 0) call mpixinit()
  if (parcal) then 
    call stopwatch(1,numnod+1)
  else
    call stopwatch(1,numnod)
  endif
#else
  call stopwatch(1,numnod)
#endif
  iparcal = 0
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
    call allocator_set_max(maxsize)
    call sync_nodes(int(memsize,8),0)
    call sync_nodes(maxsize,2)
  else
!   iparex <= 0
    call allocator_init()
    call allocator_set_max(maxsize)
    print '(a,i20,a,i14,a)', 'DIRAC serial starts by allocating ',memsize,' words (',memsize/2**17,' MB)'
    if (maxsize.gt.0) then
       print '(a,i20,a,i14,a)', '             of memory out of the allowed',  &
          maxsize,' words (',maxsize/2**17,' MB)'
    else
       print '(a)', 'DIRAC serial has no limitations in place for the amount of dynamically allocated memory'
    endif
#else
    call allocator_init()
    call allocator_set_max(maxsize)
    print '(a,i0,a,i0,a)', 'DIRAC serial starts by allocating ',memsize,' words (',memsize/2**17,' MB) of memory'
    if (maxsize.gt.0) then
       print '(a,i20,a,i14,a)', '             out of the allowed',  &
          maxsize,' words (',maxsize/2**17,' MB)'
    else
       print '(a)', 'DIRAC serial has no limitations in place for the amount of dynamically allocated memory'
    endif
#endif
    print '(a)',' '
    print '(a)','Note: maximum allocatable memory for serial run can be set by pam --aw'
    print '(a)',' '
#ifdef VAR_MPI
  endif
#endif
  if (mytid == 0) then
     call alloc(memory,memsize, id="work array at main")
     call set_mem_start(memory)
     call printtitle
     call printlogo
     call printsubtitle
     call flush(lupri)

     write(lupri, *)
     write(lupri, *)
     write(lupri, '(a)') 'Compilation info'
     write(lupri, '(a)') '----------------'
     write(lupri, *)

     call print_compilation_info()
     call print_git_revision_info()

     call flush(lupri)
! miro: do the selftest after printing dirac system info
#ifndef NO_SELFTEST
     call selftest_lapack()
#ifdef VAR_MPI
!    MASTER checks MPI library compatibility (now done at compile time)
     if ( parcal ) call selftest_mpi()
#endif
#endif
     call dirac(iparcal,mytid,mparid,numnod,memory,memsize)
  else
#ifdef VAR_MPI
     call alloc(memory,memsize)
     call dirnod(iparcal,mytid,mparid,numnod,memory,memsize)
#else
     call quit('Internal error in main.F90. Quitting.')
#endif
  endif
#ifdef VAR_MPI
  if ( iparcal .eq. 1 ) then
    call sync_nodes(idummy,1)
  endif
#endif
  if (mytid == 0) then
     call stopwatch(0,mytid)
     print '(a)', '********** E N D   of   D I R A C  output  **********'
     print '(a/)','*****************************************************'
     if (doxml) call xml_close_file
  endif
  call dealloc(memory)
#ifdef VAR_MPI
  if (iparcal .eq. 1 ) call mpixfinalize()
#endif

 call allocator_report(mytid,0)
  if (mytid == 0) then
      print '(A,F8.2,A)',' MEMGET high-water mark: ',(high_water_mark()*8)*2.0D0**(-20),' MB'
      print '(/,a)','*****************************************************'
  endif

#ifdef MOD_HSCC 
 if (DOHSCC) then
   print '(/,2x,a)','Continuing with the external HSCC code.'
   call dirprg_hscc 
 endif
#endif

#if defined (MOD_HSFS) || defined (HSFS_DYN_ALLOC)
 if (DOHSFS) then
   print '(/,a)','*****************************************************'
   print '(2x,a)','Continuing with the external HighSect-FSCC code.'
   print '(/,a)','*****************************************************'
   call HIGH_SECTOR_FSCC_METHODS
 endif
#endif

   call allocator_cleanup()

   call interface_ao_clear()
   call interface_mo_clear()

#if defined FORTRAN_MAIN_AS_SUBROUTINE
end subroutine dirprg
#else
end program dirprg
#endif

