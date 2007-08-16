C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
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
C...   This means also that persons (e.g. post-docs) leaving the research
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
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccmain.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      PROGRAM RELADC
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description-----------------------------------------------
C
C    Double group symmetry adapted one-particle propagator program
C
C---------------Routines called-------------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Markus Pernpointner
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      PARAMETER (MAXCOR = 120 000 000)
#include "../relccsd/files.inc"
#include "../relccsd/ccpar.inc"
#if defined (VAR_MPI)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC(MAXCOR)
      REAL*8 SEC,SEC0
      LOGICAL CA
      CHARACTER*10 INTFCE
C
C---------------Executable code--------------------------------------
C
      CALL CPUUSED(SEC0)
      IN = 5
      IW = 6
      INTFCE = 'DIRAC     '
      CA = .FALSE.
#if defined (VAR_MPI)
C
C     MPI initialization
C
      CALL MPI_INIT(IERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYPROC,IERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NMPROC,IERR)
      MASTER = 0
#else
      MASTER = 0
      MYPROC = 0
      NMPROC = 1
#endif
#if defined (ADC_ALONE)
      CALL ADCMAIN(MAXCOR,IN,IW,CA,INTFCE)
#else
      CALL ADCMAIN(CC,MAXCOR,IN,IW,CA,INTFCE)
#endif
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC-SEC0,SEC0,SEC
      WRITE (IW,'(//80A1)') ('#',I=1,80)
#if defined (VAR_MPI)
C
C     MPI finalization
C
      CALL MPI_FINALIZE(IERR)
#endif
C
 1000 FORMAT (//' CPU time (seconds) used in RELADC:',T50,F14.4/
     & ' CPU time (seconds) used before RELADC:',T50,F14.4/
     & ' CPU time (seconds) used in total sofar:',T50,F14.4//
     & '  --- Normal end of RELADC Run ---')
C
      END
