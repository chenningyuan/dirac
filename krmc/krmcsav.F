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
C...   group of the licence holder may not take any part of Dalton,
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
C FILE    : krmcsav.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck krmcsav */
      SUBROUTINE KRMCSAV(IOPT,WORK,LWORK)
C*****************************************************************************
C
C     Save info on LUKRMC
C
C     Input:
C        IOPT    - option (current only 0 supported)
C     
C     Written by J. Thyssen - Nov 26 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION WORK(LWORK)
      LOGICAL   FNDLAB, FND
C
#include <dcborb.h>
#include <dgroup.h>
#include <dcbopt.h>
#include <dcbgen.h>
C
#include <memint.h>
C
      CALL QENTER('KRMCSAV')
C
C 
C     Clean LUKRMC
C     ------------
C     
      REWIND(LUKRMC)
      FND = FNDLAB('SODLABEL',LUKRMC)
      BACKSPACE LUKRMC
      CALL NEWLAB('SODLABEL',LUKRMC,LUERR)
      REWIND(LUKRMC)
C 
C 
C     Save Kappa on LUKRMC
C     --------------------
C
      CALL MEMGET('REAL',KKAPPA,NZXOPTQ,WORK,KFREE,LFREE)
      REWIND(LUKRM1)
      write(6,*) 'krmcsav: nzxoptq:',nzxoptq
      CALL READT(LUKRM1,NZXOPTQ,WORK(KKAPPA))
      CALL WRTKRMC(LUKRMC,'KAPPA',WORK(KKAPPA),NZXOPTQ)
      CALL MEMREL('KRMCSAV.KAPPA',WORK,KWORK,KWORK,KFREE,LFREE)
C 
C 
C     Save old and new orbitals on LUKRMC
C     -----------------------------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
C
C     The new orbitals...
C
      REWIND(LUKRM2)
      IF (.NOT. FNDLAB('NEWORB  ',LUKRM2)) THEN
         CALL QUIT('KRMCSAV: label NEWORB not found on file "KRMC CMO"')
      END IF
      CALL READT(LUKRM2,NCMOTQ,WORK(KCMO))
      CALL WRTKRMC(LUKRMC,'NEWORB  ',WORK(KCMO),NCMOTQ)
      REWIND(LUKRMC)
C
C
C     The old orbitals...
C
      REWIND(LUKRM2)
      IF (.NOT. FNDLAB('OLDORB  ',LUKRM2)) THEN
         CALL QUIT('KRMCSAV: label OLDORB not found on file "KRMC CMO"')
      END IF
      CALL READT(LUKRM2,NCMOTQ,WORK(KCMO))
      CALL WRTKRMC(LUKRMC,'OLDORB  ',WORK(KCMO),NCMOTQ)
      REWIND(LUKRMC)
C
      CALL MEMREL('KRMCSAV.CMO',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('KRMCSAV')
      RETURN
      END
