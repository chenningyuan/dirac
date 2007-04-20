C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_ZO_DMAT_AO */
      SUBROUTINE SETUP_ZO_DMAT_AO(CMO,
     &                            DMAT_MO,
     &                            DMAT_AO,
     &                            STRING1,STRING2,
     &                            WORK,LWORK)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SCATTER_BVEC */
      SUBROUTINE SCATTER_BVEC(BVEC_EE,
     &                        BVEC_EP,
     &                        ICOUNTER_EE,
     &                        ICOUNTER_EP,
     &                        BMAT,
     &                        IZ,
     &                        IFSYM1,
     &                        IFSYM2)
C=======================================================================
C     scatter solution vector BVEC onto BMAT
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "dcborb.h"
#include "visual.h"
#include "dgroup.h"
#include "dcbbas.h"

      LOGICAL   LBIT
      DIMENSION BVEC_EE(*),
     &          BVEC_EP(*),
     &          BMAT(NORBT,NORBT,NZ),
     &          NEI(2),
     &          NES(2),
     &          NP(2)

      DO I = 1,2
        NEI(I) = NISH(I)
        NES(I) = NESH(I) - NISH(I)
        NP(I)  = NPSH(I)
      ENDDO

      II1 = 1
      II1 = II1 + (NESH(1) + NPSH(1))*(IFSYM1 - 1)
      II1 = II1 + NP(IFSYM1)
      II2 = II1 + NEI(IFSYM1) - 1
      IS1 = 1
      IS1 = IS1 + (NESH(1) + NPSH(1))*(IFSYM2 - 1)
      IS2 = IS1 + NP(IFSYM2)  - 1
      IS3 = IS2 + NEI(IFSYM2) + 1
      IS4 = IS3 + NES(IFSYM2) - 1

      DO II = II1,II2
        IF(.NOT. LBIT(LBIT_PLOT_GENERAL,5)) THEN !if not skip +- (SKIPEP)
          DO IS = IS1,IS2
            BMAT(IS,II,IZ) = BVEC_EP(ICOUNTER_EP)
            ICOUNTER_EP = ICOUNTER_EP + 1
          ENDDO
        ENDIF
        IF(.NOT. LBIT(LBIT_PLOT_GENERAL,6)) THEN !if not skip ++ (SKIPEE)
          DO IS = IS3,IS4
            BMAT(IS,II,IZ) = BVEC_EE(ICOUNTER_EE)
            ICOUNTER_EE = ICOUNTER_EE + 1
          ENDDO
        ENDIF
      ENDDO

      RETURN
      END
