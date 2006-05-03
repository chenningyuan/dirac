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
C$Id$
C
#if defined (VAR_SINGLE)
#include "single.h"
#endif
C
      PROGRAM NODPRG
C
#include "implicit.h"
#include "priunit.h"
#include "mtags.h"
      PARAMETER (MWORK = 10 000 000)
      DIMENSION WORK(MWORK)
C
C     Join PVM and return my name to master
C
      CALL PVMXMYTID(MYTID)
      CALL PVMXPARENT(MASTER)
      CALL SENDNAME(MASTER,MYTID)
C
C     Get calculation type (iprtyp)
C
      CALL PVMXRECV(MASTER,MTAG1)
      CALL PVMXUNPACK(IPRTYP,1,'INTEGER')
C
      DO WHILE(IPRTYP.GT.0)         
         CALL PVMXUNPACK(IPRINT,1,'INTEGER')
C
         IF (IPRTYP.EQ.1) THEN
            CALL HER_NODSTR(WORK,MWORK,IPRINT)
         ELSE IF (IPRTYP.EQ.2) THEN
            CALL ERI_NODSTR(WORK,MWORK,IPRINT)
C         ELSE IF (IPRTYP.EQ.3) THEN
C            CALL GEMNOD(WORK,MWORK,MASTER,MYTID,MYNUM,IPRINT)
         END IF
C
         CALL PVMXRECV(MASTER,MTAG1)
         CALL PVMXUNPACK(IPRTYP,1,'INTEGER')
      END DO
      CALL PVMXEXIT
      STOP
      END
