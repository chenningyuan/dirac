C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fckden */
      SUBROUTINE FCKDEN(GETDC,GETDV,DCAO,DVAO,CMO,DV,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get backtransformed (contravariant) density matrices.
C     Based on DENSTY and DALTON/FCKDEN
C
C     Written by J. Thyssen
C     Last revision: 1997/09/22 - jth
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION DCAO(*),DVAO(*),CMO(*),DV(*)
      DIMENSION WORK(*)
      DIMENSION IFMUL(2,2), IASH(2)
      LOGICAL GETDC,GETDV
C
      CALL QENTER('FCKDEN')
      IF ( IPRINT .GE. 9 ) CALL HEADER('Output from FCKDEN',-1)
C
C     Get active density matrix
C
      IF (GETDV) THEN
         CALL DENSDV(DVAO,CMO,DV,IPRINT,WORK,LWORK)
      END IF
C
C     Get inactive density matrix
C   
      IF (GETDC) THEN
         CALL DENSTY(DCAO,CMO,0)
C         CALL DENSTY(DCAO,CMO,IPRINT)
      END IF
C
      CALL FLSHFO(LUPRI)
      CALL QEXIT('FCKDEN')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck densdv */
      SUBROUTINE DENSDV(DVAO,CMO,DV,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Get backtransformed (contravariant) active density.
C     Based on DALTON/FCKDEN
C
C     Written by J. Thyssen
C     Last revision: 1997/09/22 - jth
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION DVAO(*),CMO(*),DV(*)
      DIMENSION WORK(LWORK)
      DIMENSION IFMUL(2,2), IASH(2), IQBF(4)
      DATA IFMUL /1,2,2,1/
      PARAMETER ( D1 = 1.0 D00, D0 = 0.0 D00 )
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
C
#include <memint.h>
C
C
      CALL QENTER('DENSDV')
      IF ( IPRINT .GE. 10 ) CALL HEADER('Output from DENSDV',-1) 
      CALL DZERO(DVAO,N2BBASXQ)
C
C     Get active density matrix
C
      IF (NASHT .EQ. 0) THEN
         WRITE(LUPRI,*) 'WARNING!!! Active density matrix'//
     &                  ' requested but there are no active electrons.'
         WRITE(LUPRI,*) 'So DVAO = 0'
         GOTO 1000
      END IF
      WRITE(6,*) 'jth',DV(1),N2BBASXQ,NTBAS(0)
C 
C     Allocate memory for temp. DVMO matrix
C
      CALL MEMGET('REAL',KTEMP,NZ*NORBT*NORBT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTEMP),NZ*NORBT*NORBT)
C
      IASH(1) = 0
      IASH(2) = NASH(1) 
C
C
      DO IF = 1,NFSYM
         DO I = 1,NASH(IF)
            DO J = 1,NASH(IF)
               DO IZ = 1,NZ
                  IADR1 = I + (J-1) * NASH(IF) + 
     &                    (IZ-1) * NASH(IF) * NASH(IF)
                  II = NISH(IF)+NPSH(IF)+I-1
                  JJ = NISH(IF)+NPSH(IF)+J-1
                  WRITE(6,*) 'jth',I,J,II,JJ,NORBT
                  WRITE(6,*) 'jth',JJ*NORBT, (IZ-1) * NORBT * NORBT
                  IADR2 = II + JJ * NORBT + (IZ-1) * NORBT * NORBT
                  WRITE(6,*) 'jth IADR',IADR1,IADR2
                  CALL FLSHFO(6)
                  WORK(KTEMP+IADR2) = DV(IADR1)
               END DO
            END DO
         END DO
      END DO
      CALL HEADER('DENSDV "Expanded" DVMO',-1)
      CALL PRQMAT(WORK(KTEMP),NORBT,NORBT,NORBT,NORBT,
     &            NZ,IPQTOQ(1,0),LUPRI)
C
      CALL MEMGET('REAL',KTEMP2,NZ*NTBAS(0)*NORBT,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KTEMP2),NZ*NTBAS(0)*NORBT)
      NZA = NZ
      NZB = NZ
      NZC = NZ
      NZD = NZ
      DO I = 1,NFSYM
         IADR1 = KTEMP + I2ORBX(I,I)
         CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NORB(I),NORB(I),
     &               DVAO(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),
     &               NZA,IQBF,
     &               WORK(IADR1),NORBT,NORBT,
     &               NZB,IPQTOQ(1,0),
     &               CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZC,
     &               IPQTOQ(1,0),
     &               CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZD,
     &               IPQTOQ(1,0),
     &               WORK(KFREE),LFREE,IPRINT)
      END DO
C
      IF ( IPRINT .GE. 0 ) THEN
         CALL HEADER('DENSDV: Active density matrix',-1)
         CALL PRQMAT(DVAO,NTBAS(0),NTBAS(0),
     &                   NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
      END IF
C
      CALL MEMREL('DENSDV',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
 1000 CONTINUE
      CALL QEXIT('DENSDV')
      CALL FLSHFO(LUPRI)
C
      RETURN
      END


