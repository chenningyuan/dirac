C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_ZO_DMAT_AO */
      SUBROUTINE SETUP_ZO_DMAT_AO(CMO,
     &                            DMAT_AO,
     &                            STRING,
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
     &          DMAT_AO(*),
     &          STRING(*),
     &          WORK(LWORK)

#include "memint.h"

!     allocate memory

      CALL MEMGET('REAL',KDMAT_MO,N2ORBXQ,WORK,KFREE,LFREE)

      CALL DZERO(WORK(KDMAT_MO),N2ORBXQ)
      CALL DZERO(DMAT_AO,N2BBASXQ)

      K = 0
      DO I = 1,NFSYM
        IF(NISH(I) .GT. 0) THEN
          DO J = 1,NESH(I)
            K = K + 1
            JOFF = (NORB(1)*(I-1) + NPSH(I))*NORBT
            JOFF = JOFF + NORB(1)*(I-1) + NPSH(I)
            JOFF = JOFF + (J-1)*NORBT + J - 1

            WORK(KDMAT_MO+JOFF) = STRING(J+NESH(1)*(I-1))
          ENDDO
        ENDIF
      ENDDO

      DO I = 1,NFSYM
        CALL QTRANS('MOAO','S',D0,
     &              NFBAS(I,0),NFBAS(I,0),NORB(I),NORB(I),
     &              DMAT_AO(1+I2BASX(I,I)),
     &              NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &              WORK(KDMAT_MO+I2ORBX(I,I)),
     &              NORBT,NORBT,1,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              CMO(1+ICMOQ(I)),
     &              NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,0)
      ENDDO

      CALL DSCAL(N2BBASXQ,D2,DMAT_AO,1)

!     release memory

      CALL MEMREL('SETUP_ZO_DMAT_AO',WORK,KDMAT_MO,KDMAT_MO,KFREE,LFREE)

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
     &                        IFSYM2,
     &                        STRING)
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
     &          STRING(*),
     &          NEI(2),
     &          NES(2),
     &          NP(2)

      DO I = 1,2
        NEI(I) = NISH(I)
        NES(I) = NESH(I) - NISH(I)
        NP(I)  = NPSH(I)
      ENDDO

      IOFF = (NESH(1) + NPSH(1))*(IFSYM1 - 1) + NP(IFSYM1)

      II1 = IOFF + 1
      II2 = II1 + NEI(IFSYM1) - 1
      IS1 = 1
      IS1 = IS1 + (NESH(1) + NPSH(1))*(IFSYM2 - 1)
      IS2 = IS1 + NP(IFSYM2)  - 1
      IS3 = IS2 + NEI(IFSYM2) + 1
      IS4 = IS3 + NES(IFSYM2) - 1

      DO II = II1,II2

        FACTOR = STRING(II - IOFF + NESH(1)*(IFSYM1-1))

        IF(.NOT. LBIT(LBIT_VISUAL,5)) THEN !if not skip +- (SKIPEP)
          DO IS = IS1,IS2
            BMAT(IS,II,IZ) = BVEC_EP(ICOUNTER_EP)*FACTOR
            ICOUNTER_EP = ICOUNTER_EP + 1
          ENDDO
        ENDIF
        IF(.NOT. LBIT(LBIT_VISUAL,6)) THEN !if not skip ++ (SKIPEE)
          DO IS = IS3,IS4
            BMAT(IS,II,IZ) = BVEC_EE(ICOUNTER_EE)*FACTOR
            ICOUNTER_EE = ICOUNTER_EE + 1
          ENDDO
        ENDIF

      ENDDO

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SETUP_FO_DMAT_AO */
      SUBROUTINE SETUP_FO_DMAT_AO(CMO,
     &                            IXVECTOR,
     &                            DMAT_AO,
     &                            STRING,
     &                            WORK,LWORK)
C=======================================================================
C     description in here ...
C-----------------------------------------------------------------------
C                                                           radovan bast
C                                              last revision: april 2007
C=======================================================================

#include "implicit.h"
#include "priunit.h"
#include "dcbqal.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "visual.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dcbxrs.h"

      CHARACTER LAB1*16,
     &          LAB2*16,
     &          TYP1*2,
     &          TYP2*2,
     &          LBUF*16,
     &          PBUF*1
      LOGICAL   FILE_EXISTS,
     &          LBIT
      DIMENSION CMO(*),
     &          DMAT_AO(*),
     &          STRING(*),
     &          WORK(LWORK)

#include "memint.h"

      INQUIRE(FILE = XVECTOR_FILE(IXVECTOR),EXIST = FILE_EXISTS)
      IF(.NOT. FILE_EXISTS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'solution vector file not found'
        CALL QUIT('ERROR in VISUAL')
      ENDIF
      LUNIT_XVECTOR = 17
      OPEN(LUNIT_XVECTOR,
     &     FILE = XVECTOR_FILE(IXVECTOR),
     &     STATUS = 'OLD',
     &     ACCESS = 'SEQUENTIAL',
     &     FORM = 'UNFORMATTED')

!     get dimensions of record IWHICH_XVECTOR(IXVECTOR)

      NBDIM_EE    = 0
      NBDIM_EP    = 0
      NR_READS    = 0
      NR_READS_EE = 0
      NR_READS_EP = 0

      REWIND LUNIT_XVECTOR
  1   CONTINUE
      READ(LUNIT_XVECTOR) LAB1
      BACKSPACE(LUNIT_XVECTOR)
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 10
      ENDIF
      READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1
      BACKSPACE(LUNIT_XVECTOR)

      IF(TYP1(1:1) .EQ. 'E') THEN
        NR_READS    = NR_READS + 1
        NR_READS_EE = NR_READS_EE + 1
        IF(TYP1(2:2) .EQ. '+') THEN
          IH_EE = 1
        ELSEIF(TYP1(2:2) .EQ. '-') THEN
          IH_EE = -1
        ELSE
          GOTO 12
        ENDIF
        READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &                      JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &                      UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS
        IF(NR_READS_EE .EQ. IWHICH_XVECTOR(IXVECTOR)) THEN
          JSYMOP_EE = JSYMOP1
          NBDIM_EE  = LEN/NZ
        ENDIF
        READ(LUNIT_XVECTOR)

      ELSEIF(TYP1(1:1) .EQ. 'P') THEN
        NR_READS    = NR_READS + 1
        NR_READS_EP = NR_READS_EP + 1
        IF(TYP1(2:2) .EQ. '+') THEN
          IH_EP = 1
        ELSEIF(TYP1(2:2) .EQ. '-') THEN
          IH_EP = -1
        ELSE
          GOTO 12
        ENDIF
        READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &                      JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &                      UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS
        IF(NR_READS_EP .EQ. IWHICH_XVECTOR(IXVECTOR)) THEN
          JSYMOP_EP = JSYMOP1
          NBDIM_EP  = LEN/NZ
        ENDIF
        READ(LUNIT_XVECTOR)

      ELSE
        GOTO 12
      ENDIF
      GOTO 1
 10   CONTINUE

!     consistency

!     LBIT(LBIT_VISUAL,5) skip +- (SKIPEP)
!     LBIT(LBIT_VISUAL,6) skip ++ (SKIPEE)

      IF(NR_READS .EQ. 0)
     &   GOTO 11

      IF(.NOT. LBIT(LBIT_VISUAL,5) .AND.
     &   .NOT. LBIT(LBIT_VISUAL,6) .AND.
     &   .NOT. NR_READS_EE .EQ. NR_READS_EP)
     &   GOTO 11

      IF(.NOT. LBIT(LBIT_VISUAL,5) .AND.
     &   NR_READS_EP .EQ. 0)
     &   GOTO 13

      IF(.NOT. LBIT(LBIT_VISUAL,6) .AND.
     &   NR_READS_EE .EQ. 0)
     &   GOTO 14

      IF(.NOT. LBIT(LBIT_VISUAL,5) .AND.
     &   NR_READS_EP .LT. IWHICH_XVECTOR(IXVECTOR))
     &   GOTO 11

      IF(.NOT. LBIT(LBIT_VISUAL,6) .AND.
     &   NR_READS_EE .LT. IWHICH_XVECTOR(IXVECTOR))
     &   GOTO 11

      IF(LBIT(LBIT_VISUAL,5) .AND .
     &   LBIT(LBIT_VISUAL,6))
     &   GOTO 15

      GOTO 20
 11   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'invalid IWHICH_XVECTOR(IXVECTOR)'
        WRITE(LUPRI,'(A)') 'or solution vector file corrupt'
        CALL QUIT('ERROR in VISUAL')
 12   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'invalid operator label'
        WRITE(LUPRI,'(A)') 'solution vector file corrupt'
        CALL QUIT('ERROR in VISUAL')
 13   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'EP requested but no EP rotation'
        WRITE(LUPRI,'(A)') 'on solution vector file'
        CALL QUIT('ERROR in VISUAL')
 14   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'EE requested but no EE rotation'
        WRITE(LUPRI,'(A)') 'on solution vector file'
        CALL QUIT('ERROR in VISUAL')
 15   CONTINUE
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'SKIPEP = SKIPEE = .TRUE.'
        CALL QUIT('ERROR in VISUAL')
 20   CONTINUE

!     allocate memory

      CALL MEMGET('REAL',KBVEC_EE ,NBDIM_EE*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBVEC_EP ,NBDIM_EP*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT    ,N2ORBXQ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMAT_AO_TEMP,N2BBASXQ,
     &            WORK,KFREE,LFREE)

!     read vector

      REWIND LUNIT_XVECTOR
      DO I = 1,IWHICH_XVECTOR(IXVECTOR)-1
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)

        IF(LBIT(LBIT_VISUAL,9) .OR. LBIT(LBIT_VISUAL,10)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF
      ENDDO

      IF(.NOT. LBIT(LBIT_VISUAL,6)) THEN !"EE" part
        IF(LBIT(LBIT_VISUAL,10)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF

        READ(LUNIT_XVECTOR)
        CALL DZERO(WORK(KBVEC_EE),NBDIM_EE*NZ)
        CALL READT(LUNIT_XVECTOR,NBDIM_EE*NZ,WORK(KBVEC_EE))

        IF(LBIT(LBIT_VISUAL,9)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF
      ELSE
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        IF(LBIT(LBIT_VISUAL,9) .OR. LBIT(LBIT_VISUAL,10)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF
      ENDIF

      IF(.NOT. LBIT(LBIT_VISUAL,6)) THEN !"EP" part
        IF(LBIT(LBIT_VISUAL,10)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF

        READ(LUNIT_XVECTOR)
        CALL DZERO(WORK(KBVEC_EP),NBDIM_EP*NZ)
        CALL READT(LUNIT_XVECTOR,NBDIM_EP*NZ,WORK(KBVEC_EP))

        IF(LBIT(LBIT_VISUAL,9)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF
      ELSE
        READ(LUNIT_XVECTOR)
        READ(LUNIT_XVECTOR)
        IF(LBIT(LBIT_VISUAL,9) .OR. LBIT(LBIT_VISUAL,10)) THEN
          READ(LUNIT_XVECTOR)
          READ(LUNIT_XVECTOR)
        ENDIF
      ENDIF

      CLOSE(LUNIT_XVECTOR,STATUS = 'KEEP')

      IF(LBIT(LBIT_VISUAL,9)) THEN
        IH_EE = 1
        IH_EP = IH_EE
      ENDIF
      IF(LBIT(LBIT_VISUAL,10)) THEN
!       IH_EE = -1
!       dirty trick to multiply the j part with i
        IH_EE = 1
        IH_EP = IH_EE
      ENDIF

!     scatter vector onto MO density matrix

      IHERMITICITY = IH_EE
      ICOUNTER_EE  = 1
      ICOUNTER_EP  = 1
      CALL DZERO(WORK(KBMAT),N2ORBXQ)

      DO IZ = 1,NZ
        IF(NFSYM .EQ. 1) THEN
          CALL SCATTER_BVEC(WORK(KBVEC_EE),
     &                      WORK(KBVEC_EP),
     &                      ICOUNTER_EE,
     &                      ICOUNTER_EP,
     &                      WORK(KBMAT),
     &                      IZ,
     &                      1,1,
     &                      STRING)
        ELSE
          IF(IHERMITICITY .EQ. 1) THEN
            CALL SCATTER_BVEC(WORK(KBVEC_EE),
     &                        WORK(KBVEC_EP),
     &                        ICOUNTER_EE,
     &                        ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,2,
     &                        STRING)
            CALL SCATTER_BVEC(WORK(KBVEC_EE),
     &                        WORK(KBVEC_EP),
     &                        ICOUNTER_EE,
     &                        ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,1,
     &                        STRING)
          ELSE
            CALL SCATTER_BVEC(WORK(KBVEC_EE),
     &                        WORK(KBVEC_EP),
     &                        ICOUNTER_EE,
     &                        ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,1,
     &                        STRING)
            CALL SCATTER_BVEC(WORK(KBVEC_EE),
     &                        WORK(KBVEC_EP),
     &                        ICOUNTER_EE,
     &                        ICOUNTER_EP,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,2,
     &                        STRING)
          ENDIF
        ENDIF
      ENDDO
      CALL DSCAL(N2ORBXQ,-D1,WORK(KBMAT),1)

      JSYMOP                 = JSYMOP_EE
      IREP_XVECTOR(IXVECTOR) = JSYMOP - 1
      JOPSY                  = JBTOF(IREP_XVECTOR(IXVECTOR),1)

      DO IZ = 1,NZ
        IOFF = NORBT*NORBT*(IZ - 1)
        CALL FULMAT('S',NORBT,NORBT,WORK(KBMAT+IOFF)) !symmetrize
      ENDDO

!     get AO density matrix

      CALL DZERO(DMAT_AO,N2BBASXQ)
      IF(.NOT. LBIT(LBIT_VISUAL,6)) THEN !if not skip ++ (SKIPEE)
        CALL DMOLI3(0,IH_EE,DMAT_AO,CMO,WORK(KBMAT),
     &              WORK,KFREE,LFREE)
      ENDIF
      IF(.NOT. LBIT(LBIT_VISUAL,5)) THEN !if not skip +- (SKIPEP)
        CALL DZERO(WORK(KDMAT_AO_TEMP),N2BBASXQ)
        CALL DMOLI3(1,IH_EP,WORK(KDMAT_AO_TEMP),CMO,WORK(KBMAT),
     &              WORK,KFREE,LFREE)
        CALL DAXPY(N2BBASXQ,D1,WORK(KDMAT_AO_TEMP),1,DMAT_AO,1)
      ENDIF

      IF(NZ .LT. 4) THEN
        DO IZ = 1,NZ
          IQ = IPQTOQ(IZ,IREP_XVECTOR(IXVECTOR))
          IOFF = NTBAS(0)*NTBAS(0)*(IZ-1)
          CALL Q2BPHASE('D',IQ,1,DMAT_AO(1+IOFF))
        ENDDO
      ENDIF

      CALL DSCAL(N2BBASXQ,D2,DMAT_AO,1)

!     release memory

      CALL MEMREL('SETUP_FO_DMAT_AO',WORK,KBVEC_EE,KBVEC_EE,KFREE,LFREE)

      RETURN
      END
