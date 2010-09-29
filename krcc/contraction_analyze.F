      SUBROUTINE CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                               TOP,NTOP,
     &                               HOP,IHTP,RESOP,
     &                               RESOPS,RESOPF,
     &                               WORK,KFREE,LFREE)
*
* Will find which cluster operator TOP can contract NPART,NHOLE indices
* of a given Hamiltonian or intemediate operator (IHTP of HOP) and find
* a connection in resulting operator RESOP (aided by looking in
* RESOPS,RESOPF)
*
* Lasse 2010
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
*
      DIMENSION WORK(*)
*
      INTEGER TOP(4*NGAS,*),HOP(4*NGAS,*),RESOP(4*NGAS,*)
      INTEGER RESOPS(*),RESOPF(*) !for now
* Scratch
      INTEGER INDXP(2),INDXH(2) !up to 4 indices can be contracted 
      INTEGER INUMCON(4,4)
*
* First find which indices in HOP to be contracted
* See GAS_INDICES_CONTRACTED for assumptions on this
*
      CALL GAS_INDICES_CONTRACTED(NPARTTOT,NHOLETOT,
     &                            HOP(1,IHTP),INDXP,INDXH)
      stop
* Now find the number of ways NPART,NHOLE can be contracted among
* NPARTTOT,NHOLETOT
*
      CALL NUMBER_OF_CONTRACTIONS(NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                            INDXP,INDXH,INUMCON,NUMCON)
*
* Now find which T operator has the match
*
      DO INUM =1,NUMCON
        CALL CONTRACT_WITH_T(INI,INUMCON(INUM,1),
     &                       TOP,NTOP,NPART,NHOLE,
     &                       RESOP,RESOPS,RESOPF)
      END DO
*
      RETURN
      END
*
      SUBROUTINE GAS_INDICES_CONTRACTED(NPARTTOT,NHOLETOT,
     &                                  IOP,INDXP,INDXH)
*
* Find which NPART and NHOLE indices to be contracted 
* under assumption that The size of the particle space decrease with
* NGAS and hole space increase and have alpha before beta.
* Returns absolute index in operator CAAB form back.
*
* Lasse 2010
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
* Input
      INTEGER IOP(NGAS,4)
* Output
      INTEGER INDXP(2),INDXH(2)
*
      NTEST = 100
*
* Initialize
*
      DO I=1,2
        INDXP(I) = 0
        INDXH(I) = 0
      END DO
*
* Find particle and hole adresses
*
      IPINDX = 0
      IHINDX = 0
      DO ICA = 1, 2
        DO IAB = 1, 2
          ICAAB = (ICA-1)*2+ IAB
          DO IGAS = 1, NGAS
            IF(IOP(IGAS,ICAAB).EQ.1) THEN
              IF(ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
                IPINDX = IPINDX + 1
                INDXP(IPINDX) = (IGAS-1)*4+ICAAB
              ELSE IF(ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.1) THEN
                IHINDX = IHINDX + 1
                INDXH(IHINDX) = (IGAS-1)*4+ICAAB
              END IF
            ELSE IF(IOP(IGAS,ICAAB).EQ.2) THEN
              IF(ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
                IPINDX = IPINDX + 1
                INDXP(IPINDX) = (IGAS-1)*4+ICAAB
                IPINDX = IPINDX + 1
                INDXP(IPINDX) = (IGAS-1)*4+ICAAB
              ELSE IF(ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.1) THEN
                IHINDX = IHINDX + 1
                INDXH(IHINDX) = (IGAS-1)*4+ICAAB
                IHINDX = IHINDX + 1
                INDXH(IHINDX) = (IGAS-1)*4+ICAAB
              END IF
            END IF
          END DO 
        END DO
      END DO
*
      IF(IPINDX.NE.NPARTTOT.OR.IHINDX.NE.NHOLETOT) THEN
        WRITE(6,*) 'Lasse has made a bug in the code'
        STOP 'Lasse has made a bug in the code'
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The possible absolute adresses '
        WRITE(6,*) ' in H to be contracted '
        WRITE(6,*) ' Number of particle ',IPINDX
        DO I =1,IPINDX
          WRITE(6,*) INDXP(I)
        END DO 
        WRITE(6,*) ' Number of hole ',IHINDX
        DO I =1,IHINDX
          WRITE(6,*) INDXH(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE NUMBER_OF_CONTRACTIONS(NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                            INDXP,INDXH,INUMCON,NUMCON)
*
* Will find the different ways NPART and NHOLE can be contracted when
* having NPARTTOT and NHOLETOT in total
*
#include "implicit.inc"
*
      INTEGER INDXP(2),INDXH(2),INUMCON(4,4)
*
      NTEST = 100
*
* Initialize
*
      IZERO = 0
      CALL ISETVC(INUMCON,IZERO,16)
*
* Possible different particle contractions
*
      IF(NPART.EQ.NPARTTOT.OR.NPART.EQ.0) THEN
        NPOSPART = 1
      ELSE
        NPOSPART = 2
      END IF
*
      IF(NHOLE.EQ.NHOLETOT.OR.NHOLE.EQ.0) THEN
        NPOSHOLE = 1
      ELSE
        NPOSHOLE = 2
      END IF
*
      NPOSTOT = NPOSPART * NPOSHOLE
*
      IF(NPOSTOT.EQ.1) THEN
        DO I=1,2
          INUMCON(1,I) = INDXP(I)
          INUMCON(1,I+2) = INDXH(I)
        END DO
      ELSE IF(NPOSTOT.EQ.2) THEN
        IF(NPOSPART.EQ.2) THEN
          DO I=1,2
            DO J=1,2
              IF(I.EQ.J) THEN
                INUMCON(I,J) = INDXP(J)
              END IF
              INUMCON(I,J+2) = INDXH(J)
            END DO
          END DO
        ELSE
          DO I=1,2
            DO J=1,2
              INUMCON(I,J) = INDXP(J)
              IF(I.EQ.J) THEN
                INUMCON(I,J+2) = INDXH(J)
              END IF
            END DO
          END DO
        END IF
      ELSE IF(NPOSTOT.EQ.4) THEN
        DO I=1,4
          CALL EVENODD(I,ITEST)
          DO J=1,2
            IF(J.EQ.ITEST) THEN
              INUMCON(I,J) = INDXP(I)
            END IF
            IF(I.GE.3) THEN
              IF(J.EQ.2) THEN
                INUMCON(I,J+2) = INDXH(J)
              END IF
            ELSE
              IF(J.EQ.1) THEN
                INUMCON(I,J+2) = INDXH(J)
              END IF
            END IF
          END DO
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of possible ways to contract ',NPOSTOT
        DO I =1,NPOSTOT
          WRITE(6,*) ' Absolute entries of contraction ',I
          DO J =1,4
            WRITE(6,*) INUMCON(I,J)
          END DO
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE CONTRACT_WITH_T(INI,INUMCON,
     &                           TOP,NTOP,NPART,NHOLE,
     &                           RESOP,RESOPS,RESOPF)
*
* Routine to find the possible contractions of a given set of indices
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      INTEGER TOP(4*NGAS,*),INUMCON(4),ITINDX(4)
      INTEGER RESOP(4*NGAS,*),RESOPS(*),RESOPF(*)
*
      CALL FROM_O_TO_T_INDX(NPART,NHOLE,INUMCON,ITINDX)
*
      DO ITOP =1,NTOP
        ICONNECT = 0
        IF(TOP(ITINDX(1),ITOP).GE.1) THEN
          IF(NPART+NHOLE.GE.2) THEN
* First index match. Now check the rest.
            DO IPH = 2,NPART+NHOLE
              IF(TOP(ITINDX(IPH),ITOP).GE.1) THEN
* Check if two succesive indices are the same
                IF(ITINDX(IPH).EQ.ITINDX(IPH-1)) THEN
                  IF(TOP(ITINDX(IPH),ITOP).GE.2) THEN
* we have found a connection
                    ICONNECT = 1
                  ELSE
                    ICONNECT = 0
                    EXIT
                  END IF
                END IF
* we so far have found a connection
                ICONNECT = 1
              ELSE
* No connection
                ICONNECT = 0
                EXIT
              END IF
            END DO
          ELSE
* we have found a connection with one index contraction
            ICONNECT = 1
          END IF
        END IF
* store connection
        IF(ICONNECT.EQ.1) THEN

        END IF
      END DO
*
      RETURN
      END
*
      SUBROUTINE FROM_O_TO_T_INDX(NPART,NHOLE,INUMCON,ITINDX)
*
* Reforms from indices contracted to indices to be found
* and moves indices to front and joins common indices
*
#include "implicit.inc"
*
      INTEGER INUMCON(4),ITEMP(4),ITINDX(4)
*
* First remove zeroes 
*
      ICOUNT = 0
      DO I =1,4
        ITEMP(I) = 0
        IF(INUMCON(I).NE.0) THEN
          ICOUNT = ICOUNT + 1
          ITEMP(ICOUNT) = INUMCON(I)
        END IF
      END DO
*
* Copy back
*
      CALL ICOPVE(ITEMP,INUMCON,4)
*
* Calculate the T indices
*
      DO IPART =1,NPART
        ITINDX(IPART) = INUMCON(IPART) - 2
      END DO
*
      DO IHOLE =1,NHOLE
        ITINDX(IHOLE) = INUMCON(IHOLE) + 2
      END DO
*
      RETURN
      END
