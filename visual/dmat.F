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
#include "dcbqal.h"

      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)

#include "visual.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "pgroup.h"
#include "dcbxrs.h"

      DIMENSION CMO(*),
     &          DMAT_AO(*),
     &          STRING(*),
     &          WORK(LWORK)

#include "memint.h"

!     allocate memory

      CALL MEMGET('REAL',KDMAT_MO,N2ORBXQ,WORK,KFREE,LFREE)

!     get MO density matrix
     
      CALL MODENS(WORK(KDMAT_MO),STRING,IPRVIS)

!     get AO density matrix

      CALL DZERO(DMAT_AO,N2BBASXQ)

!#ifdef PRG_DIRAC
#ifdef allow_virtual_orbital_density
!     radovan: this code part is for Markus
!              to allow virtual orbital densities
!              in DMOLI3 this block gets filtered out
!              i will make the code nicer later, for
!              the moment this fix is ok ...
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
#else
      JSYMOP = 1
      IREPOP = 0
      JOPSY  = JBTOF(IREPOP,1)
      IH     = 1

      CALL DMOLI3(0,IH,DMAT_AO,CMO,WORK(KDMAT_MO),
     &            WORK,KFREE,LFREE)
#endif /* allow_virtual_orbital_density */

      IF(NZ .LT. 4) THEN
        DO IZ = 1,NZ
          IQ = IPQTOQ(IZ,IREPOP)
          IOFF = NTBAS(0)*NTBAS(0)*(IZ-1)
          CALL Q2BPHASE('D',IQ,1,DMAT_AO(1+IOFF))
        ENDDO
      ENDIF

!     release memory

      CALL MEMREL('SETUP_ZO_DMAT_AO',WORK,KDMAT_MO,KDMAT_MO,KFREE,LFREE)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* SCATTER_BVEC */
      SUBROUTINE SCATTER_BVEC(BVEC,
     &                        ICOUNTER,
     &                        TYP,
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

      CHARACTER TYP*2
      LOGICAL   LBIT
      DIMENSION BVEC(*),
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
        IF(TYP(1:1) .EQ. 'P') THEN
          DO IS = IS1,IS2
            BMAT(IS,II,IZ) = BVEC(ICOUNTER)*FACTOR
            ICOUNTER = ICOUNTER + 1
          ENDDO
        ENDIF
        IF(TYP(1:1) .EQ. 'E') THEN
          DO IS = IS3,IS4
            BMAT(IS,II,IZ) = BVEC(ICOUNTER)*FACTOR
            ICOUNTER = ICOUNTER + 1
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
        WRITE(LUPRI,'(A,A6,A)') 'solution vector file ',
     &    XVECTOR_FILE(IXVECTOR),' not found'
        CALL QUIT('ERROR in VISUAL')
      ENDIF
      LUNIT_XVECTOR = 17
      OPEN(LUNIT_XVECTOR,
     &     FILE = XVECTOR_FILE(IXVECTOR),
     &     STATUS = 'OLD',
     &     ACCESS = 'SEQUENTIAL',
     &     FORM = 'UNFORMATTED')

!     get number of vectors on file

      NR_VECTORS = 0
      REWIND LUNIT_XVECTOR
 1    CONTINUE
      READ(LUNIT_XVECTOR) LAB1
      IF(LAB1 .EQ. 'END_OF_THIS_FILE') THEN
        GOTO 2
      ELSE
        NR_VECTORS = NR_VECTORS + 1
        READ(LUNIT_XVECTOR)
        GOTO 1
      ENDIF
 2    CONTINUE

!     errors

      IF(NR_VECTORS .EQ. 0) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'NR_VECTORS .EQ. 0'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

      IF(IWHICH_XVECTOR(IXVECTOR) .GT. NR_VECTORS) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'IWHICH_XVECTOR(IXVECTOR) .GT. NR_VECTORS'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

      IF(IWHICH_XVECTOR(IXVECTOR) .LT. 1) THEN
        WRITE(LUPRI,'(A)') 'ERROR in VISUAL:'
        WRITE(LUPRI,'(A)') 'IWHICH_XVECTOR(IXVECTOR) .LT. 1'
        CALL QUIT('ERROR in VISUAL')
      ENDIF

!     move to the right record

      IVECTOR = 0
      REWIND LUNIT_XVECTOR
 3    CONTINUE
      READ(LUNIT_XVECTOR)
      IVECTOR = IVECTOR + 1
      IF(IVECTOR .EQ. IWHICH_XVECTOR(IXVECTOR)) THEN
        BACKSPACE(LUNIT_XVECTOR)
        GOTO 4
      ELSE
        READ(LUNIT_XVECTOR)
        GOTO 3
      ENDIF
 4    CONTINUE

!     get type and dimensions of record IWHICH_XVECTOR(IXVECTOR)

      READ(LUNIT_XVECTOR) LAB1, LAB2, TYP1, TYP2, FREQ1, FREQ2,
     &                    JSYMOP1, JSYMOP2, JTIMOP1, JTIMOP2,
     &                    UNUSED, LEN, INTFLG, ERGRSP, NBAS, NORBS

      IF(TYP1(2:2) .EQ. '+') IHERMITICITY =  1
      IF(TYP1(2:2) .EQ. '-') IHERMITICITY = -1

      JSYMOP                 = JSYMOP1
      IREP_XVECTOR(IXVECTOR) = JSYMOP - 1
      IREPOP                 = JSYMOP - 1
      JOPSY                  = JBTOF(IREP_XVECTOR(IXVECTOR),1)
      NBDIM                  = LEN/NZ

!     allocate memory

      CALL MEMGET('REAL',KBVEC    ,NBDIM*NZ,
     &            WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBMAT    ,N2ORBXQ,
     &            WORK,KFREE,LFREE)

!     read vector

      CALL DZERO(WORK(KBVEC),NBDIM*NZ)
      CALL READT(LUNIT_XVECTOR,NBDIM*NZ,WORK(KBVEC))

      CLOSE(LUNIT_XVECTOR,STATUS = 'KEEP')

!     scatter vector onto MO density matrix

      CALL DZERO(WORK(KBMAT),N2ORBXQ)
      ICOUNTER = 1

!     my static electric perturbations are ungerade
!               magnetic                     gerade
!
!     i don't know codewise how to extract gerade/ungerade
!     from the general operator irrep
!
!     with frequency-dependent magnetic perturbations
!     the code will fail when calculating with symmetry
!     but for the moment it does its job ...

      IGERADE = -IHERMITICITY
      IF(DABS(FREQ1) .GT. 1.0D-9) IGERADE = -1

      DO IZ = 1,NZ
        IF(NFSYM .EQ. 1) THEN
          CALL SCATTER_BVEC(WORK(KBVEC),
     &                      ICOUNTER,
     &                      TYP1,
     &                      WORK(KBMAT),
     &                      IZ,
     &                      1,1,
     &                      STRING)
        ELSE
          IF(IGERADE .EQ. -1) THEN
            CALL SCATTER_BVEC(WORK(KBVEC),
     &                        ICOUNTER,
     &                        TYP1,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,2,
     &                        STRING)
            CALL SCATTER_BVEC(WORK(KBVEC),
     &                        ICOUNTER,
     &                        TYP1,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,1,
     &                        STRING)
          ELSE
            CALL SCATTER_BVEC(WORK(KBVEC),
     &                        ICOUNTER,
     &                        TYP1,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        1,1,
     &                        STRING)
            CALL SCATTER_BVEC(WORK(KBVEC),
     &                        ICOUNTER,
     &                        TYP1,
     &                        WORK(KBMAT),
     &                        IZ,
     &                        2,2,
     &                        STRING)
          ENDIF
        ENDIF
      ENDDO
      CALL DSCAL(N2ORBXQ,-D1,WORK(KBMAT),1)

      DO IZ = 1,NZ
        IOFF = NORBT*NORBT*(IZ - 1)
        CALL FULMAT('S',NORBT,NORBT,WORK(KBMAT+IOFF)) !symmetrize
      ENDDO

!     get AO density matrix

      CALL DZERO(DMAT_AO,N2BBASXQ)
          IF(TYP1(1:1) .EQ. 'E') THEN
        CALL DMOLI3(0,IHERMITICITY,DMAT_AO,CMO,WORK(KBMAT),
     &              WORK,KFREE,LFREE)
      ELSEIF(TYP1(1:1) .EQ. 'P') THEN
        CALL DMOLI3(1,IHERMITICITY,DMAT_AO,CMO,WORK(KBMAT),
     &              WORK,KFREE,LFREE)
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

      CALL MEMREL('SETUP_FO_DMAT_AO',WORK,KBVEC,KBVEC,KFREE,LFREE)

      RETURN
      END
