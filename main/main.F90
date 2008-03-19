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
  use os_utils
  implicit none
#include "infpar.h"
#include "dcbgen.h"
  real(kind=8), allocatable :: memory(:)
  integer :: memsize, iparcal, status, length, ios
! initialize parcal (dcbgen.h) - SK - 19-03-2008
  parcal  = .false.
#ifdef VAR_MPI 
  call mpixinit()
#else
  iparcal = 0
  mytid   = 0
  mparid  = 0
  numnod  = 1
#endif
!    SK -19-03-2008 bugfix: call DIRAC / DIRNOD with iparcal (not parcal)
  if (parcal) iparcal = 1
  if (mytid == 0) then
     memsize = get_environment_integer('DIRWRK',INSTALL_WRKMEM)
     print '(a,i0,a,i0,a)', 'DIRAC master allocating ',memsize,' words (',memsize/2**16,' MB) of memory'     
     allocate(memory(memsize))
     call printlogo
     call dirac(iparcal,mytid,mparid,numnod,memory,memsize)
  else
     memsize = get_environment_integer('DIRNOD',INSTALL_WRKMEM)
     print '(a,i0,a,i0,a)', 'DIRAC slave allocating ',memsize,' words (',memsize/2**16,' MB) of memory'     
     allocate(memory(memsize))
     call dirnod(iparcal,mytid,mparid,numnod,memory,memsize)
  endif
#ifdef VAR_MPI
  call mpixfinalize()
#endif
end program dirprg
