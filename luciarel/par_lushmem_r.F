C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T.
CFleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V.
CLarsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T.
CWinther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the
Cresearch
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
C
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : par_lushmem_r.F
C
C
#if defined (VAR_MPI2)
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_MEM_T( MY_T_INFO )
C
C     Written by  S. Knecht         - November 17 2007
C
C**********************************************************************
C
C     routine for memory allocation of the integral/density block T 
C     to use in parallel shared-memory MPI LUCIAREL runs
C
C     pointer has to be defined in calling routine as 
C     POINTER (MY_T_PTR, T_BUFF(IDUMMY_BUFF_SZ))
C
C     LEN_T_BUFF stored on /IBUFF_LENGTH/   (krmc_shmem.h)
C     MY_T_PTR   stored on /MPI_MEMPOINTER/ (krmc_shmem.h)
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar_rel.inc"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      LOGICAL CINFO
      CINFO = .FALSE.
C
C     allocate T memory
C
      IF( MY_T_INFO .ne. 0 ) CINFO = .TRUE.
C
      CALL MPIXMEM_ALLOC( MY_T_PTR, 1, LEN_T_BUFF, MY_T_INFO, CINFO )
C
      END
***********************************************************************
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_T_WIN( T_BUFF, MY_T_COMM, MY_LEN_T )
C
C     Written by  S. Knecht         - November 17 2007
C
C**********************************************************************
C
C     routine for window creation of the integral/density block T to 
C     use in parallel shared-memory MPI LUCIAREL runs
C
C     T_BUFF is allocated using MPI_ALLOC_MEM.
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar_rel.inc"
#include "mxdim_mpi2.h"
#include "typesz_mpi2.h"
#include "krmc_shmem.h"
      INTEGER(KIND=MPI_ADDRESS_KIND) MY_LEN_T
C
C     open T window
C
      CALL MPIXWINDOW_CRT( T_BUFF, MY_LEN_T, 1, ISIZE_dp, 
     &                     MY_T_COMM, MY_T_WIN)
C
      END
***********************************************************************
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE REL_MEM_T_WIN( T_BUFF )
C
C     Written by  S. Knecht         - November 19 2007
C
C**********************************************************************
C
C     deallocate and and destroy memory window for integral/density 
C     block T. 
C
C     
C     Last revision:
C
***********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "krmc_shmem.h"
#include "lucipar_rel.inc"
C
C     destroy T memory window
C
      CALL MPIXWINDOW_FREE( MY_T_WIN )
C
C     release memory
C
      CALL MPIXMEM_FREE(T_BUFF)
C
      END
***********************************************************************
#else
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_LUSHMEM_REL
*     dummy routine for normal compilation
      END
#endif

