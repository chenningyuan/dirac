!dirac_copyright_start
!      Copyright (c) 2017 by the authors of DIRAC.
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
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck effqed */
      SUBROUTINE EFFQED(FOCK,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Build contribution from all-electron pseudopotential representing
C     QED effects.
C
C     Written by Trond Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
C
#include "dcbgen.h"
#include "nuclei.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcb_effqed.h"      
      DIMENSION FOCK(*),WORK(*)
      CHARACTER TEXT*74
      KFRSAV = KFREE
      CALL MEMGET2('INTE','KINF',KINF,NTBAS(0),WORK,KFREE,LFREE)
C      
C.....Loop over atomic orbitals
C     ICENT runs over symmetry independent centers
C     INDEP runs over symmetry dependent centers
      INDEP=1
      ICENT=1
      DO IATOM = 1,NONTYP
C.......Find the total number of basis functions for this atom type
         CALL LABCOUNT(MBAS,WORK(KINF),NTBAS(0),1,-1,ICENT,1,-1)
C.......Read C1 coefficients and orbitalenergies
        CALL OPNFIL(LUCOEF,REFFIL(IATOM),'OLD','EFFQED')        
        READ (LUCOEF,END=10,ERR=20) TEXT,NSYM,NNEG,NPOS,NAO,DUM
        IF(NAO.NE.MBAS) THEN
           WRITE(LUPRI,'(A,A,I3,/,A,A6,A,I6)') 
     &     'Effective QED potential generation: ',
     &     '  Wrong number of basis function for atomic type ', IATOM,
     &     '  Number read from                     ',
     &     REFFIL(INDEP),' : '  ,MBAS
           WRITE(LUPRI,*) 'NAO, MBAS:',NAO, MBAS
           CALL QUIT
     &     ('EFFQED: Wrong number of basis functions !')
        ENDIF
        NTOT = NNEG+NPOS
        NQA  = NAO*NTOT*NZC1
        CALL MEMGET2('REAL','KQA',KQA,NQA ,WORK,KFREE,LFREE)
        CALL READT(LUCOEF,NQA,WORK(KQA))
        CALL MEMGET2('REAL','KEA',KEA,NTOT,WORK,KFREE,LFREE)
        CALL READT(LUCOEF,NTOT,WORK(KEA))
        CLOSE(LUCOEF,STATUS='KEEP')
C        NN = INDEP
        DO J = 0,NONT(IATOM)-1
          INDEP = INDEP + NUCDEG(ICENT+J)
        ENDDO
        ICENT = ICENT + NONT(IATOM)
      ENDDO
      CALL MEMREL('EFFQED',WORK,1,KFRSAV,KFREE,LFREE)
      STOP 'End of test'
C
      RETURN
 10   CONTINUE
      CALL QUIT('EFFQED: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('EFFQED: ERROR reading TEXT')
      END
      
