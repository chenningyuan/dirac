C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccmain.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      PROGRAM RELCCSD
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description-----------------------------------------------
C
C    Double group symmetry adapted unrestricted coupled cluster code
C    CRAY version uses dynamic memory allocation and variable file buffers
C
C---------------Routines called-------------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
C---------------Common Blocks--------------------------------------
C
      PARAMETER (MAXCOR = 13 000 000)
#if !defined (SYS_AIX) || defined (SYS_CRAY) || defined (SYS_IRIX) || defined (SYS_HPUX) || defined (SYS_T3D)
      COMMON // CC(MAXCOR)
#endif
      INCLUDE "files.inc"
C
C---------------Local variables--------------------------------------
C
      REAL*8 CC
      REAL*8 SEC
C
C---------------Executable code--------------------------------------
C
C     Define primary in and output units.
C     Other files units are defined in files include block.
C
      IN = 5
      IW = 6
C
      CALL CCHEADER
      MXCOR = MAXCOR
      CALL CCMAIN(CC,MXCOR)
C
      CALL CPUUSED(SEC)
      WRITE (IW,1000) SEC
      CALL QUIT("NORMAL END OF PROGRAM")
 1000 FORMAT (//' Total amount of CPU time (seconds) used :',F14.4,
     &//'  --- Normal end of Run ---')
      END
