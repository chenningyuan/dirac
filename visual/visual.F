C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* VISUAL */
      SUBROUTINE VISUAL(WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "dcbgen.h"
#include "visual.h"
#include "maxorb.h"
#include "dcbham.h"
#include "dgroup.h"

      LOGICAL   FILE_EXISTS, LBIT
      DIMENSION WORK(LWORK)

      CALL QENTER('VISUAL')

#include "memint.h"

!     CALL VISUAL_CHECK_CONSISTENCY()
!     CALL VISUAL_REPORT()

!     allocate memory

      CALL MEMGET('REAL',KCMO         ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG         ,NORBT,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_MO     ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO     ,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO_TEMP,N2BBASXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING1     ,NESH(1),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSTRING2     ,NESH(2),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO         ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAO1        ,NTBAS(0)*3,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCNT        ,NTBAS(0),
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF         ,NTBAS(0)*4,
     &            WORK,KFREE,LFREE)

      CALL SETUP_ORBITAL_STRINGS(WORK(KSTRING1),WORK(KSTRING2))

!     coefficients

      CALL DZERO(WORK(KCMO),N2BBASXQ)
      CALL DZERO(WORK(KEIG),NORBT)

      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','VISUAL')
      CALL REACMO(LUCOEF,'DFCOEF',WORK(KCMO),WORK(KEIG),UNUSED,
     &            UNUSED,6,WORK,KFREE,LFREE)
      CLOSE(LUCOEF,STATUS = 'KEEP')

      CALL FNDOCC(WORK(KEIG),WORK,KFREE,LFREE)

!     density matrices

      IF(LBIT(LBIT_PLOT_GENERAL,4)) THEN !using the solution vector
        CALL SETUP_FO_DMAT_AO(WORK(KCMO),
     &                        WORK(KDMAT_AO),WORK(KDMAT_AO_TEMP),
     &                        WORK(KSTRING1),WORK(KSTRING2),
     &                        WORK(KFREE),LFREE)
        CALL SETUP_DMAT_AO(WORK(KCMO),
     &                     WORK(KDMAT_MO),
     &                     WORK(KDMAT_AO_TEMP),
     &                     WORK(KSTRING1),WORK(KSTRING2),
     &                     WORK(KFREE),LFREE)
      ELSE !"normal"
        CALL SETUP_DMAT_AO(WORK(KCMO),
     &                     WORK(KDMAT_MO),
     &                     WORK(KDMAT_AO),
     &                     WORK(KSTRING1),WORK(KSTRING2),
     &                     WORK(KFREE),LFREE)
      ENDIF

      CALL VISUAL_DRIVER(WORK(KDMAT_AO),
     &                   WORK(KDMAT_AO_TEMP),
     &                   WORK(KGAO),
     &                   WORK(KGAO1),
     &                   WORK(KNCNT),
     &                   WORK(KBUF),
     &                   WORK(KFREE),LFREE)

!     release memory

      CALL MEMREL('VISUAL',WORK,KCMO,KCMO,KFREE,LFREE)
      CALL QEXIT('VISUAL')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_ORBITAL_STRINGS */
      SUBROUTINE SETUP_ORBITAL_STRINGS(STRING1,STRING2)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"

      DIMENSION STRING1(*),
     &          STRING2(*)

      CALL DFILL(STRING1,D0,NESH(1))
      CALL DFILL(STRING1,D1,NISH(1))
      IF(NFSYM .GT. 1) THEN
        CALL DFILL(STRING2,D0,NESH(2))
        CALL DFILL(STRING2,D1,NISH(2))
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_DMAT_AO */
      SUBROUTINE SETUP_DMAT_AO(CMO,
     &                         DMAT_MO,
     &                         DMAT_AO,
     &                         STRING1,STRING2,
     &                         WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION CMO(*),
     &          DMAT_MO(*),
     &          DMAT_AO(*),
     &          STRING1(*),
     &          STRING2(*),
     &          WORK(LWORK)

#include "memint.h"

      CALL DZERO(DMAT_MO,N2ORBXQ)
      CALL DZERO(DMAT_AO,N2BBASXQ)

      K = 0
      DO I = 1,NFSYM
        IF(NISH(I) .GT. 0) THEN
          DO J = 1,NESH(I)
            K = K + 1
            JOFF = (NORB(1)*(I-1) + NPSH(I))*NORBT
            JOFF = JOFF + NORB(1)*(I-1) + NPSH(I)
            JOFF = JOFF + (J-1)*NORBT + J - 1

            IF(I .EQ. 1) DMAT_MO(1+JOFF) = STRING1(J)
            IF(I .EQ. 2) DMAT_MO(1+JOFF) = STRING2(J)
          ENDDO
        ENDIF
      ENDDO

      DO I = 1,NFSYM
        CALL QTRANS('MOAO','S',D0,
     &              NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &              DMAT_AO(1+I2BASX(I,I)),
     &              NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &              DMAT_MO(1+I2ORBX(I,I)),
     &              NORBT,NORBT,1,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,0)
      ENDDO

      CALL DSCAL(N2BBASXQ,D2,DMAT_AO,1)

      RETURN
      END
