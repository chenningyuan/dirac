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
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE LUCI_TASK_NODES(CIRUN)
C
C     Introduce LUCIAREL action - Node part
C
C     Written by  S. Knecht         - April 13 2007 (Friday!)
C
C**********************************************************************
#include "implicit.h"
#include "lucipar_rel.inc"
C
      character*6 CIRUN
*
      NTEST = 0
*
      if (NTEST.ge.1) then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' ////////////////////////////////// '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' //  KR-MCSCF calling LUCIAREL   // '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' //                              // '
        write(MSLVOUT_REL,*) ' ////////////////////////////////// '
      end if
*
      if (CIRUN.eq.'NDET  ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' NDET   flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) ' The number of CI dets will be computed.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'CIINII') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' CIINII flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Determining configurational start 
     & guess'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'RSTRMC') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' RSTRMC flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Determining info for restart of MCSCF.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'SIGMA ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' SIGMA  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing a sigma vector.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'SIGMAD') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' SIGMAD flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing a sigma vector.'
        write(MSLVOUT_REL,*) '     Saving H diagonal also.'
        write(LUPRI,*)
      else if (CIRUN.eq.'DENS1 ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' DENS1  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing only 1-particle density 
     & matrices.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'DENS2 ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' DENS2  flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing 1- and 2-particle dens. 
     & matrices.'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'ANALYZ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' ANALYZ flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Analyzing the current CI vector(s).'
        write(MSLVOUT_REL,*)
      else if (CIRUN.eq.'DIAG  ') then
        write(MSLVOUT_REL,*)
        write(MSLVOUT_REL,*) ' DIAG   flag call '
        write(MSLVOUT_REL,*) ' ------           '
        write(MSLVOUT_REL,*) '   Computing CI diagonal.'
        call quit('NOT IMPLEMENTED YET!')
        write(MSLVOUT_REL,*)
      end if
*
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE SYNC_NODES_REL_DRV(F1,F2,ECORE_KRMC,CIRUN,LUWRT)
C
C     Synchronize MASTER and NODES - driver routine
C
C     Written by  S. Knecht         - April 13 2007 (Friday!)
C
C**********************************************************************
#include "implicit.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "lucipar_rel.inc"
      CHARACTER*6 CIRUN
      DIMENSION F1(*), F2(*)
C
C     introduce CIRUN for NODES
C
#if defined (VAR_MPI)
      CALL MPI_BCAST(CIRUN,6,MPI_CHARACTER,MASTER,MPI_COMM_WORLD,IERR)
#endif
C
C
C
C
      END
