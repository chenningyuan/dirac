!
!...   Copyright (c) 2004 by the authors of Dirac (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   Dirac, a relativistic ab initio electronic structure program,
!...   release DIRAC04.0 (2004),
!...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
!...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
!...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
!...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
!...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
!...   (http://dirac.chem.sdu.dk).
!...
!...   For a suitable BibTEX entry, see:
!...      http://dirac.chem.sdu.dk/doc/reference.shtml
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the
!...   research group of the licence holder may not take any part of Dirac,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For questions concerning this copyright write to:
!...      dirac-admin@dirac.chem.sdu.dk
!...
!...   For information on how to get a licence see:
!...      http://dirac.chem.sdu.dk
!
! this file by ulfek march 2008

program dirprg
  use memory_allocator
  use os_utils
  implicit none
#include "infpar.h"
#include "dcbgen.h"
  real(kind=8), allocatable :: memory(:)
  integer :: memsize, iparcal, status, length, ios, idummy
#ifdef VAR_MPI 
  call mpixinit()
#else
  parcal  = .false.
  iparcal = 0
  mytid   = 0
  mparid  = 0
  numnod  = 0
#endif
  call stopwatch(1,numnod)
  if (parcal) iparcal = 1
  if (mytid == 0) then
     memsize = get_environment_integer('DIRWRK',INSTALL_WRKMEM)
  else
     memsize = get_environment_integer('DIRNOD',INSTALL_WRKMEM)
  endif
#ifdef VAR_MPI
  call sync_nodes(memsize,0)
#else
  print '(a,i0,a,i0,a)', 'DIRAC serial allocating ',memsize,' words (',memsize/2**17,' MB) of memory'
#endif
  if (mytid == 0) then
     call alloc(memory,memsize)
     call printtitle
     call printlogo
     call printsubtitle
     call sysinf(iparcal)
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
  call sync_nodes(idummy,1)
#endif
  if (mytid == 0) then
     call stopwatch(0,mytid)
     print '(a)','********** E N D   of   D I R A C  output  **********'
  endif
  call dealloc(memory)
#ifdef VAR_MPI
  call mpixfinalize()
#endif
end program dirprg
