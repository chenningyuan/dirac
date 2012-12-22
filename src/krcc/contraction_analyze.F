      SUBROUTINE CONTRACTION_ANALYZE(INI,NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                               TOP,NTOP,
     &                               HOP,IHTP,RESOP,
     &                               RESOPS,RESOPF,
     &                               NUMCON,IFACTOR,
     &                               NCONNECTTOT,
     &                               RESULTANTOP,
     &                               CONTRACTEDT,
     &                               INUMNUM,
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
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
*
      DIMENSION WORK(*)
*
      INTEGER TOP(4*NGAS,*),HOP(4*NGAS,*),RESOP(4*NGAS,*)
      INTEGER RESOPS(*),RESOPF(*) !for now
      INTEGER RESULTANTOP(*),CONTRACTEDT(*)
* Scratch
      INTEGER INDXP(2),INDXH(2) !up to 4 indices can be contracted 
      INTEGER INUMCON(4,4),IFACTOR(4),INUMNUM(4)
*
      NTEST = 00
*
* Initialize IFACTOR and INUMNUM
*
      CALL MEMCHK_KRCC(WORK)
      CALL ISETVC(IFACTOR,1,4)
      CALL ISETVC(INUMNUM,0,4)
*
* First find which indices in HOP to be contracted
* See GAS_INDICES_CONTRACTED for assumptions on this
*
      CALL MEMCHK_KRCC(WORK)
      CALL GAS_INDICES_CONTRACTED(NPARTTOT,NHOLETOT,
     &                            HOP(1,IHTP),INDXP,INDXH)
*
* Now find the number of ways NPART,NHOLE can be contracted among
* NPARTTOT,NHOLETOT
*
      CALL MEMCHK_KRCC(WORK)
      CALL NUMBER_OF_CONTRACTIONS(NPART,NHOLE,NPARTTOT,NHOLETOT,
     &                            INDXP,INDXH,INUMCON,NUMCON)
*
* Will now reduce the contraction by eliminating all that will give the
* same result and thereby give a factor in front. Thereby eleminating
* the zeroes in INUMCON and change how it is read since it will now
* really depend on NHOLE and NPART.
*
      CALL MEMCHK_KRCC(WORK)
      CALL REDUCE_CONTRACTIONS(INUMCON,NUMCON,IFACTOR)
*
* Now find which T operator has the match
*
*
* Remember to store where the diferent contractions start and the type
* of contraction
*
      CALL MEMCHK_KRCC(WORK)
      NCONNECTTOT = 0
      DO INUM =1,NUMCON
        CALL CONTRACT_WITH_T(INI,INUMCON(1,INUM),HOP(1,IHTP),
     &                       TOP,NTOP,NPART,NHOLE,
     &                       RESOP,RESOPS,RESOPF,NCONNECT,
     &                       RESULTANTOP(NCONNECTTOT+1),
     &                       CONTRACTEDT(NCONNECTTOT+1))
        INUMNUM(INUM) = NCONNECT
        NCONNECTTOT = NCONNECTTOT + NCONNECT
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' For operator '
        CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IHTP),1)
        WRITE(6,*) ' With ',NPART,' NPART ',
     &             ' and ',NHOLE,' NHOLE contracted'
        WRITE(6,*) ' from a total ',NPARTTOT,' particles ',
     &             ' and ',NHOLETOT,' holes'
        WRITE(6,*) ' A total of ',NCONNECTTOT,' connections was found'
        IF(NCONNECTTOT.GE.1) THEN
          ICOUNT = 0
          DO I =1,NUMCON
            WRITE(6,*) ' Number of connections ',INUMNUM(I),
     &                 ' for combination ',I,
     &                 ' with permutation factor ',IFACTOR(I)       
            IF(INI.EQ.1) THEN
              WRITE(6,*) '======Connections====='      
              DO J =1+ICOUNT,INUMNUM(I)+ICOUNT
                WRITE(6,*) ' Cluster op ',CONTRACTEDT(J),
     &                     ' Res op ',RESULTANTOP(J)
              END DO
              ICOUNT = ICOUNT + INUMNUM(I)
            END IF        
          END DO
        END IF
      END IF
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
      NTEST = 00
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
                INDXP(IPINDX) = (ICAAB-1)*NGAS+IGAS
              ELSE IF(ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.1) THEN
                IHINDX = IHINDX + 1
                INDXH(IHINDX) = (ICAAB-1)*NGAS+IGAS
              END IF
            ELSE IF(IOP(IGAS,ICAAB).EQ.2) THEN
              IF(ICA.EQ.2.AND.IHPVGAS_AB(IGAS,IAB).EQ.2) THEN
                IPINDX = IPINDX + 1
                INDXP(IPINDX) = (ICAAB-1)*NGAS+IGAS
                IPINDX = IPINDX + 1
                INDXP(IPINDX) = (ICAAB-1)*NGAS+IGAS
              ELSE IF(ICA.EQ.1.AND.IHPVGAS_AB(IGAS,IAB).EQ.1) THEN
                IHINDX = IHINDX + 1
                INDXH(IHINDX) = (ICAAB-1)*NGAS+IGAS
                IHINDX = IHINDX + 1
                INDXH(IHINDX) = (ICAAB-1)*NGAS+IGAS
              END IF
            END IF
          END DO 
        END DO
      END DO
*
      IF(IPINDX.NE.NPARTTOT.OR.IHINDX.NE.NHOLETOT) THEN
        WRITE(6,*) 'Lasse has made a bug in the code'
        WRITE(6,*) ' IPINDX ',IPINDX,' NPARTTOT ',NPARTTOT
        WRITE(6,*) ' IHINDX ',IHINDX,' NHOLETOT ',NHOLETOT
        CALL WRT_SPOX_TP_CC_KRCC(IOP,1)
        STOP 'Lasse has made a bug in the code'
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' For the Hamiltonian '
        CALL WRT_SPOX_TP_CC_KRCC(IOP,1)
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
      NTEST = 00
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
          DO J=1,2
            IF(NPART.EQ.0) THEN
              INUMCON(J,I) = 0
            ELSE IF(NPART.EQ.2) THEN
              INUMCON(J,I) = INDXP(J)
            ELSE IF(NPART.EQ.1) THEN
              IF(J.EQ.1) THEN
                INUMCON(J,I) = INDXP(J)
              ELSE
                INUMCON(J,I) = 0
              END IF
            END IF
            IF(NHOLE.EQ.0) THEN
              INUMCON(J+2,I) = 0
            ELSE IF(NHOLE.EQ.2) THEN
              INUMCON(J+2,I) = INDXH(J)
            ELSE IF(NHOLE.EQ.1) THEN
              IF(J.EQ.1) THEN
                INUMCON(J+2,I) = INDXH(J)
              ELSE
                INUMCON(J+2,I) = 0
              END IF
            END IF
          END DO
        END DO
      ELSE IF(NPOSTOT.EQ.2) THEN
        IF(NPOSPART.EQ.2) THEN
          DO I=1,2
            DO J=1,2
              IF(I.EQ.J) THEN
                INUMCON(J,I) = INDXP(J)
              END IF
              IF(NHOLE.EQ.0) THEN
                INUMCON(J+2,I) = 0
              ELSE IF(NHOLE.EQ.2) THEN
                INUMCON(J+2,I) = INDXH(J)
              ELSE IF(NHOLE.EQ.1) THEN
                IF(J.EQ.1) THEN
                  INUMCON(J+2,I) = INDXH(J)
                ELSE
                  INUMCON(J+2,I) = 0
                END IF
              END IF
            END DO
          END DO
        ELSE
          DO I=1,2
            DO J=1,2
              IF(NPART.EQ.0) THEN
                INUMCON(J,I) = 0
              ELSE IF(NPART.EQ.2) THEN
                INUMCON(J,I) = INDXP(J)
              ELSE IF(NPART.EQ.1) THEN
                IF(J.EQ.1) THEN
                  INUMCON(J,I) = INDXP(J)
                ELSE
                  INUMCON(J,I) = 0
                END IF
              END IF
              IF(I.EQ.J) THEN
                INUMCON(J+2,I) = INDXH(J)
              END IF
            END DO
          END DO
        END IF
      ELSE IF(NPOSTOT.EQ.4) THEN
        DO I=1,4
          CALL EVENODD(ITEST,I)
          DO J=1,2
            IF(J.EQ.ITEST) THEN
              INUMCON(J,I) = INDXP(J)
            END IF
            IF(I.GE.3) THEN
              IF(J.EQ.2) THEN
                INUMCON(J+2,I) = INDXH(J)
              END IF
            ELSE
              IF(J.EQ.1) THEN
                INUMCON(J+2,I) = INDXH(J)
              END IF
            END IF
          END DO
        END DO
      END IF
      NUMCON = NPOSTOT
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of possible ways to contract ',NUMCON
        DO I =1,NUMCON
          WRITE(6,*) ' Absolute entries of contraction ',I
          DO J =1,4
            WRITE(6,*) INUMCON(J,I)
          END DO
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE REDUCE_CONTRACTIONS(INUMCON,NUMCON,IFACTOR)
*
* Simple routine to eliminate all of the same type of contractions
*
#include "implicit.inc"
*
      INTEGER INUMCON(4,4),IFACTOR(4)
* Scratch
      INTEGER ITEMP(4,4),ISIM(4)
*
      NTEST = 00
*
* First eliminate the zeroes
*
      DO I =1,NUMCON
        ICOUNT = 0
        DO J =1,4
          ITEMP(J,I) = 0
          IF(INUMCON(J,I).NE.0) THEN
            ICOUNT = ICOUNT + 1
            ITEMP(ICOUNT,I) = INUMCON(J,I)
          END IF
        END DO
      END DO
* Copy back 
      CALL ICOPVE(ITEMP,INUMCON,16)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Zeroes removed'
        CALL IWRTMA(INUMCON,4,4,4,4)
      END IF
*
*  compare first against rest
*
      ISTART = 2
      IZERO = 0
      IREST = NUMCON
      IF(NUMCON.GE.2) THEN
      DO
        NSIM = 0
        NOTSIM = 0
        DO I =ISTART,NUMCON
          ISIM(I) = 0
        END DO
        ILENGTH = 16-(ISTART-1)*4
        CALL ISETVC(ITEMP,IZERO,ILENGTH)
        DO I =ISTART,NUMCON
          DO J =1,ICOUNT
            IF(INUMCON(J,ISTART-1).NE.INUMCON(J,I)) THEN
              ISIM(I)  = 1
C             EXIT
            END IF
          END DO  
        END DO
* Check for zeroes in ISIM since these will be the same contraction
        DO I =ISTART,NUMCON
          IF(ISIM(I).EQ.1) THEN
C            NSIM = NSIM - 1
            NOTSIM = NOTSIM + 1
* Copy rest
            CALL ICOPVE(INUMCON(1,I),ITEMP(1,NOTSIM),4)
          ELSE
            NSIM = NSIM + 1
            IFACTOR(ISTART-1) = IFACTOR(ISTART-1) + 1
          END IF
        END DO
C       CALL IWRTMA(INUMCON,4,4,4,4)
C       CALL IWRTMA(ITEMP,4,4,4,4)
        CALL ICOPVE(ITEMP,INUMCON(1,ISTART),ILENGTH)
        ISTART = ISTART + 1
        NUMCON = NUMCON - NSIM
        IREST = IREST - NSIM
        IF(ISTART.GT.NUMCON) EXIT    
      END DO
      END IF
*      
      NUMCON = IREST
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Remaining terms ',IREST
        WRITE(6,*) ' Actual contractions '
        CALL IWRTMA(INUMCON,4,4,4,4)
        WRITE(6,*) ' Prefactor for each term '
        DO I =1,4
          WRITE(6,*) IFACTOR(I)
        END DO
      END IF        
*
      RETURN
      END
*
      SUBROUTINE CONTRACT_WITH_T(INI,INUMCON,HOP,
     &                           TOP,NTOP,NPART,NHOLE,
     &                           RESOP,RESOPS,RESOPF,NCONNECT,
     &                           RESULTANTOP,
     &                           CONTRACTEDT)
*
* Routine to find the possible contractions of a given set of indices
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      INTEGER TOP(4*NGAS,*),INUMCON(4),ITINDX(4),HOP(4*NGAS)
      INTEGER RESOP(4*NGAS,*),RESOPS(*),RESOPF(*)
      INTEGER RESULTANTOP(*),CONTRACTEDT(*)
* Scratch
      INTEGER IOP(4*NGAS)
*
      NTEST = 00
C     print*,'enter CONTRACT_WITH_T'
*
      CALL FROM_O_TO_T_INDX(NPART,NHOLE,NGAS,INUMCON,ITINDX)
*
      NCONNECT = 0
      DO ITOP =1,NTOP
        ICONNECT = 0
C     CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
        IF(TOP(ITINDX(1),ITOP).GE.1) THEN
C         print*,'inside'
C         print*,'NPART,NHOLE',NPART,NHOLE
          IF(NPART+NHOLE.GE.2) THEN
* First index match. Now check the rest.
            DO IPH = 2,NPART+NHOLE
              IF(TOP(ITINDX(IPH),ITOP).GE.1) THEN
C               print*,IPH,'index matches'
* Check if two succesive indices are the same
                IF(ITINDX(IPH).EQ.ITINDX(IPH-1)) THEN
                  IF(TOP(ITINDX(IPH),ITOP).GE.2) THEN
* we have found a connection so far
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
C       print*,'ICONNECT',ICONNECT
* connection found. Is it also a part of the resulting operator
        IF(ICONNECT.EQ.1) THEN
C         CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
* Find operator form of connection
          CALL OP_FORM(IOP,HOP,TOP(1,ITOP),
     &                 INUMCON,ITINDX,
     &                 NPART,NHOLE,NGAS)
* Analyse operator
          CALL ANALYZE_OP(IOP,NGAS,N,MK,MUB)
* Check if this operator is part of resulting operator
          CALL IS_IT_IN(IOP,N,MK,MUB,RESOP,RESOPS,RESOPF,
     &                  ICONNECT,INUM)
          IF(ICONNECT.EQ.1) THEN
* Connection complete. Store it. Coming soon! Is here
            NCONNECT = NCONNECT + 1
            IF(INI.EQ.1) THEN
* Store connection
              RESULTANTOP(NCONNECT) = INUM
              CONTRACTEDT(NCONNECT) = ITOP
            END IF        
C           CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
          ELSE
C           CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
          END IF
        ELSE
C         CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
        END IF
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of connections ',NCONNECT,' for operator '
        CALL WRT_SPOX_TP_CC_KRCC(HOP,1)
        IF(INI.EQ.1) THEN
          DO I =1, NCONNECT
            WRITE(6,*) '==========================='
            WRITE(6,*) ' Connected to T op ',CONTRACTEDT(I)
            CALL WRT_SPOX_TP_CC_KRCC(TOP(1,CONTRACTEDT(I)),1)
            WRITE(6,*) ' Resulting in op ',RESULTANTOP(I)
            CALL WRT_SPOX_TP_CC_KRCC(RESOP(1,RESULTANTOP(I)),1)
          END DO
        END IF
      END IF   
*
      RETURN
      END
*
      SUBROUTINE FROM_O_TO_T_INDX(NPART,NHOLE,NGAS,INUMCON,ITINDX)
*
* Reforms from indices contracted to indices to be found
* and moves indices to front and joins common indices
*
#include "implicit.inc"
*
      INTEGER INUMCON(4),ITEMP(4),ITINDX(4)
*
      NTEST = 00
*
* Calculate the T indices
*
      DO IPART =1,NPART
        ITINDX(IPART) = INUMCON(IPART) - 2*NGAS
      END DO
*
      DO IHOLE =1+NPART,NHOLE+NPART
        ITINDX(IHOLE) = INUMCON(IHOLE) + 2*NGAS
      END DO
*
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Indices in the Hamiltonian to be contracted '
        DO I = 1,NPART+NHOLE
          WRITE(6,*) INUMCON(I)
        END DO
        WRITE(6,*) ' The corresponding indices in T '
        DO I = 1,NPART+NHOLE
          WRITE(6,*) ITINDX(I)
        END DO
      END IF
*      
      RETURN
      END
*
      SUBROUTINE OP_FORM(IOP,HOP,TOP,
     &                   INUMCON,ITINDX,
     &                   NPART,NHOLE,NGAS)
*
* Will calculate the resultant operator from a contraction between HOP
* and TOP defined from INUMCON
*
#include "implicit.inc"
*
      INTEGER IOP(4*NGAS),HOP(4*NGAS),TOP(4*NGAS)
      INTEGER INUMCON(4),ITINDX(4)
*
      NTEST = 00
*      
* First add HOP and TOP
*   
      DO I =1,4*NGAS
        IOP(I) = HOP(I) + TOP(I)
      END DO
*
* Apply contraction by removing contracted operator from IOP
*
      DO I =1,NPART+NHOLE
        IOP(INUMCON(I)) = IOP(INUMCON(I)) - 1
        IOP(ITINDX(I)) = IOP(ITINDX(I)) - 1
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Hamiltonian Operator '
        CALL WRT_SPOX_TP_CC_KRCC(HOP,1)
        WRITE(6,*) ' Cluster Operator '
        CALL WRT_SPOX_TP_CC_KRCC(TOP,1)
        WRITE(6,*) ' Was contracted to '
        CALL WRT_SPOX_TP_CC_KRCC(IOP,1)
        WRITE(6,*) ' Number of particles contracted ',NPART
        WRITE(6,*) ' Number of holes contracted ',NHOLE
        WRITE(6,*) ' Hamiltonian indices contracted '
        DO I =1,NPART+NHOLE
          WRITE(6,*) INUMCON(I)
        END DO
        WRITE(6,*) ' Cluster operator indices contracted '
        DO I =1,NPART+NHOLE
          WRITE(6,*) ITINDX(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE IS_IT_IN(IOP,N,MK,MUB,RESOP,RESOPS,RESOPF,
     &                    ICONNECT,INUM)
*
* See if IOP belongs to RESOP
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
*
      INTEGER IOP(4*NGAS),RESOPS(*),RESOPF(*),RESOP(4*NGAS,*)
*
      NTEST = 00
*
* Find range to look for operator
*
      CALL OP_RANGE(N,MK,MUB,RESOPS,RESOPF,IFROM,ITO,ICONNECT,
     &              MX_EXC_LEVEL,IMAXKRFLIP,IMUBMAX)
*
* We now know the start and the end of the possible intermediate
* operators which we have to search though. Therefore search!
*
      INUM = 0
      IF(ICONNECT.EQ.1) THEN
        DO I =IFROM,ITO
          IDIFF = 0
          DO J =1,4*NGAS
            IDIFF = IOP(J) - RESOP(J,I)
            IF(IDIFF.NE.0) EXIT
          END DO
          IF(IDIFF.EQ.0) THEN
* We have found the operator now get out of here                  
            INUM = I
            EXIT
          END IF
        END DO
      END IF
* In case it is not in T dagger
      IF(INUM.EQ.0) ICONNECT = 0      
*
      IF(NTEST.GE.100) THEN
        IF(ICONNECT.EQ.1) THEN
          WRITE(6,*) ' A connection to operator ',INUM,' was found '
          CALL WRT_SPOX_TP_CC_KRCC(IOP,1)
        ELSE
          WRITE(6,*) ' No connection found for operator '
          CALL WRT_SPOX_TP_CC_KRCC(IOP,1)
        END IF
        IF(IFROM.NE.0) THEN        
          WRITE(6,*) ' Searched from ',IFROM,' TO ',ITO
        ELSE
          WRITE(6,*) ' No search was needed '
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE OP_RANGE(N,MK,MUB,RESOPS,RESOPF,IFROM,ITO,ICONNECT,
     &                    IEXCITE,IKRFLIP,IMUB)
*
* Checks if operator would be accepted on N,MK,MUB basis
* and calculates the range where to search for an operator
*
#include "implicit.inc"
*
      INTEGER RESOPS(0:IEXCITE,-IKRFLIP:IKRFLIP,-IMUB:IMUB)
      INTEGER RESOPF(0:IEXCITE,-IKRFLIP:IKRFLIP,-IMUB:IMUB)
*
* Check N,MK,MUB
*
      NTEST = 00
*
      IF(N.GT.IEXCITE) THEN
        ICONNECT = 0
        IFROM = 0
        ITO =0
        RETURN
      END IF
*
      IF(ABS(MK).GT.IKRFLIP) THEN
        ICONNECT = 0
        IFROM = 0
        ITO =0
        RETURN
      END IF
*
      IF(ABS(MUB).GT.IMUB) THEN
        ICONNECT = 0
        IFROM = 0
        ITO =0
        RETURN
      END IF
*
* Now for the start and end
*
      IFROM = RESOPS(N,MK,MUB)
      ITO = RESOPF(N,MK,MUB)
*
      ILENGTH = ITO-IFROM
      IF(ILENGTH.LT.0.OR.IFROM.EQ.0.OR.ITO.EQ.0) THEN
        ICONNECT = 0
      END IF        
*
      IF(NTEST.GE.100) THEN
        IDIFF = 0
        WRITE(6,*) ' Start and finish for Operators in oprange'
        DO I =0,IEXCITE
          DO J = -IKRFLIP,IKRFLIP
            DO K = -IMUB,IMUB
              IF(RESOPF(I,J,K)-RESOPS(I,J,K).GE.0) THEN
                IF(RESOPF(I,J,K).GE.1) THEN
                WRITE(6,*) ' Operators with ',
     &                     'N ',I,' Mk',J,' Mub ',K,
     &                     ' begins ',RESOPS(I,J,K),
     &                     ' ends ',RESOPF(I,J,K)
                IDIFF = IDIFF + 1
                END IF
              END IF
            END DO
          END DO
        END DO
        WRITE(6,*) ' Number of different operators ',IDIFF
      END IF
*
      RETURN
      END
*
C       CALL CONTRACTION_TESTER(NGAS,NH22,WORK(KNH22),
C    &                          NH22TOM12,WORK(KNH22TOM12),WORK(KNH22T12),
C    &                          WORK(KCOMBNH22TOM12),
C    &                          WORK(KFACNH22TOM12),
C    &                          WORK(KSNH22TOM12),WORK(KFNH22TOM12))
      SUBROUTINE CONTRACTION_TESTER(NOP,IOP,
     &                              NCONT,IOPRES,
     &                              IOPT,
     &                              IPOSCON,
     &                              IPERM,
     &                              IPREFAC,
     &                              IPOSCONPERM,
     &                              ISTARTOP,IENDOP,
     &                              ITOP,NTOP,
     &                              IRESOP,NRESOP,
     &                              IHOP,IHM)
*
* Routine for testing and printing of all possible contractions any
* given operator can have including the operator contracted with and the
* resultant operator along with prefactors of contractions
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*
      INTEGER IOP(NOP),IPOSCON(NOP),IPREFAC(4,NOP),IPERM(NOP)
      INTEGER IPOSCONPERM(4,NOP)
      INTEGER IOPRES(NCONT),IOPT(NCONT)
      INTEGER ISTARTOP(NOP),IENDOP(NOP)
      INTEGER ITOP(4*NGAS,NTOP),IHOP(4*NGAS,NOP),IRESOP(4*NGAS,NRESOP)
* Scratch
      INTEGER ITEMP(NOP)
*
      ICOUNT = 0
*
      DO IH =1,NOP
        ITEMP(IH) = ICOUNT
        WRITE(6,*) ' Contractions for operator ',IH
        IF(IHM.EQ.1) THEN
* For the Hamiltonian operators
        IJ = IOP(IH) 
        WRITE(6,*) ' Absolute address for op ', IJ
        CALL WRT_SPOX_TP_CC_KRCC(IHOP(1,IJ),1)
        ELSE IF(IHM.EQ.2) THEN
* For the intermediate operators
        WRITE(6,*) ' Absolute address for op ', IH
        CALL WRT_SPOX_TP_CC_KRCC(IHOP(1,IH),1)
        END IF
        WRITE(6,*) ' Number of possible contractions ',IPOSCON(IH)
        WRITE(6,*) ' The number of different permutations ',IPERM(IH)
        IF(IPOSCON(IH).GE.1) THEN
          DO I =1,IPERM(IH)
            WRITE(6,*) 'Number of possible contractions ',
     &      IPOSCONPERM(I,IH),' for contraction number ',I
            WRITE(6,*) 'Permutation factor for contraction ',
     &      IPREFAC(I,IH)
            IF(IPOSCONPERM(I,IH).GE.1) THEN
              DO J=1,IPOSCONPERM(I,IH)
                ICOUNT = ICOUNT + 1
                WRITE(6,*) ' Operator contracted with ',IOPT(ICOUNT)
                CALL WRT_SPOX_TP_CC_KRCC(ITOP(1,IOPT(ICOUNT)),1)
                WRITE(6,*) ' Operators contracted to ',IOPRES(ICOUNT)
                CALL WRT_SPOX_TP_CC_KRCC(IRESOP(1,IOPRES(ICOUNT)),1)
              END DO
            END IF
          END DO
        END IF
        WRITE(6,*) ' =================================== '
      END DO
*
* For internal count check
*
      DO IH =1,NOP
        IF(ISTARTOP(IH).NE.(ITEMP(IH)+1)) STOP ' Wrong start count '
        IF(IH.LT.NOP) THEN
          IF(IENDOP(IH).NE.ITEMP(IH+1)) STOP 'Wrong end count '
        ELSE
          IF(IENDOP(IH).NE.ICOUNT) STOP 'Wrong end count for end'
        END IF
      END DO
*
      RETURN
      END
*
      SUBROUTINE OP_COMPARE(IWAY,IOP,NOP,IMOP,NMOP,
     &                      IVECS,IVECF,IEXCITE,IKRFLIP,ITMUB,IMTOH,HOP)
* 
* Compares the Hamiltonian types to the intermediates
* IWAY = 1 maps Hamiltonians to intermediates
* IWAY = 2 maps intermediates to Hamiltonians
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
*
      INTEGER IOP(NOP),IMOP(4*NGAS,NMOP)
      INTEGER IVECS(0:IEXCITE,-IKRFLIP:IKRFLIP,-ITMUB:ITMUB)
      INTEGER IVECF(0:IEXCITE,-IKRFLIP:IKRFLIP,-ITMUB:ITMUB)
      INTEGER IMTOH(NOP) ! Can also be H to M if IWAY = 2
      INTEGER HOP(4*NGAS,*)
* Scratch
      INTEGER INDX(4)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to OP_COMPARE '
        WRITE(6,*) ' Input Hamiltonian types '
        DO IHOP =1,NOP
          CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(IHOP)),1)
        END DO
        WRITE(6,*) ' Input Cluster operator types '
        CALL WRT_SPOX_TP_CC_KRCC(IMOP,NMOP)
      END IF
C     print*,'NOP,NMOP',NOP,NMOP
*
      CALL ISETVC(IMTOH,0,NOP) 
*
      DO IHOP =1,NOP
* Analyze H operator
        CALL ANALYZE_OP(HOP(1,IOP(IHOP)),NGAS,INH,IMK,IMUB)
C       print*,'INH,IMK,IMUB',INH,IMK,IMUB
* Find non-zero indices
        CALL NONZERO_INDX(HOP(1,IOP(IHOP)),INH,INDX,IDIFFINDX,NGAS)
* Compare indices to intermediate
        ISTART = IVECS(INH,IMK,IMUB)
        IEND = IVECF(INH,IMK,IMUB)
C       print*,'ISTART,IEND',ISTART,IEND
        DO I = ISTART,IEND
          ISIM = 0
          DO J = 1,IDIFFINDX
            IF(HOP(INDX(J),IOP(IHOP)).NE.IMOP(INDX(J),I)) THEN
* No match
C             print*,' no match '
              ISIM = 1
              EXIT
            END IF
          END DO
          IF(ISIM.EQ.0) THEN
* We have a match. Store the connection
C           print*,' we have a match'
            IF(IWAY.EQ.1) THEN
              IMTOH(I) = IOP(IHOP)
            ELSE IF(IWAY.EQ.2) THEN
              IMTOH(IHOP) = I
            END IF
            EXIT
          END IF
        END DO 
      END DO
*
      IF(NTEST.EQ.100) THEN
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' H operator match for T operator '
          DO I = 1,NOP
            WRITE(6,*) ' T op ',I,' match H op ',IMTOH(I)
            IF(IMTOH(I).NE.0) THEN
              WRITE(6,*) ' H op match '
              CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IMTOH(I)),1)
              WRITE(6,*) ' T op match '
C             CALL WRT_SPOX_TP_CC_KRCC(IMOP(1,IMTOH(I)),1)
              CALL WRT_SPOX_TP_CC_KRCC(IMOP(1,I),1)
            ELSE
              WRITE(6,*) ' No match for H op '
              CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(I)),1)
            END IF
          END DO
        ELSE IF(IWAY.EQ.2) THEN
          WRITE(6,*) ' T operator match for H operator '
          DO I = 1,NOP
            WRITE(6,*) ' H op ',IOP(I),' match T op ', IMTOH(I)
            IF(IMTOH(I).NE.0) THEN
              WRITE(6,*) ' H op match '
              CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(I)),1)
              print*,'I,IOP(I)',I,IOP(I)
              WRITE(6,*) ' T op match '
              CALL WRT_SPOX_TP_CC_KRCC(IMOP(1,IMTOH(I)),1)
              print*,'I,IMTOH(I)',I,IMTOH(I)
            ELSE
              WRITE(6,*) ' No match for H op '
              CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(I)),1)
            END IF
          END DO
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE NONZERO_INDX(IOP,NINDX,INDX,ICOUNT,NGAS)
*
* Will find non-zero entries in a given operator type with a max of
* NINDX non-zero indices.
*
#include "implicit.inc"
      INTEGER IOP(4*NGAS),INDX(NINDX)
*
      NTEST = 00
*
      ICOUNT = 0
      IEND = 0
* Zero of INDX
      CALL ISETVC(INDX,0,NINDX)
*
      DO I=1,4*NGAS
        IF(IOP(I).NE.0) THEN
          ICOUNT = ICOUNT + 1
          INDX(ICOUNT) = I
          IEND = IEND + IOP(I)
          IF(IEND.EQ.2*NINDX) EXIT
        END IF
      END DO
*
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' Number of non-zero indices ',ICOUNT
        DO I=1,ICOUNT
          WRITE(6,*) ' Index number ',I,' Absolute address ',INDX(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE IDIM_INT(WORK,KFREE,LFREE,IMULTFAC,NPART,NHOLE,IWAY)
*
* Subroutine that wil dimension and allocate arrays for the
* intermediates.
* First version all possible intermediates will be allocated even if
* they are not reached ==> wasting memory
* For the final version only intermediates reached should be allocated
* for and an intermediate should only be in memory when needed. Meaning
* this routine should be obsolete and the context merged into
* CC_VEC_FNC_KRCC.
*                               Lasse 2011
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "ctccp.inc"
*
      DIMENSION WORK(*)
*
* Totally symmetric
*
      ISYMINT = 1
*
      NTEST = 100
      IDUM = 0
      ZERO = 0.0D0
*
      MX_ST_TSOSO_INI = MX_ST_TSOSO
      MX_ST_TSOSO_BLK_INI = MX_ST_TSOSO_BLK
      MX_TBLK_INI = MX_TBLK
      MX_SBSTR_INI = MX_SBSTR
*
*  Have to check these does not exceed those for T by simply comparing against them 
*  MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK. Largest will be choosen!
*
      IF(NPART.EQ.0.AND.NHOLE.EQ.1) THEN
        IF(IWAY.EQ.1) THEN
*
* Setup and allocate for M01
*
      CALL MEMGET('INTE',KLLSOBEX_M01,NINTER01,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX_M01,NINTER01,WORK,KFREE,LFREE)
*
      CALL IDIM_TCC_KRCC(WORK(KINTM01),NINTER01,ISYMINT,
*                        actual operator, number,sym
     &     MX_ST_TSOSO_INT,MX_ST_TSOSO_BLK_INT,MX_TBLK_INT,
     &     WORK(KLLSOBEX_M01),WORK(KLIBSOBEX_M01),LEN_M01_VEC_CC,
* output   number of elements, offset, total length. Names needs to be
* found
     &     MX_SBSTR_INT,
* seems to be dummy/not used
     &     IDUM,IDUM,-1,IDUM,IDUM)
*
      CALL MEMGET('REAL',KM01,IMULTFAC*LEN_M01_VEC_CC,
     &            WORK,KFREE,LFREE)
      CALL SETVEC(WORK(KM01),IMULTFAC*ZERO,LEN_M01_VEC_CC)
*
* Make sure MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR are the largest possible
* (at least for now)
*
      IF(MX_ST_TSOSO_INT.GT.MX_ST_TSOSO) THEN
        MX_ST_TSOSO = MX_ST_TSOSO_INT
      END IF
      IF(MX_ST_TSOSO_BLK_INT.GT.MX_ST_TSOSO_BLK) THEN
        MX_ST_TSOSO_BLK = MX_ST_TSOSO_BLK_INT
      END IF
      IF(MX_TBLK_INT.GT.MX_TBLK) THEN
        MX_TBLK = MX_TBLK_INT
      END IF
      IF(MX_SBSTR_INT.GT.MX_SBSTR) THEN
        MX_SBSTR = MX_SBSTR_INT
      END IF
*
        ELSE
          CALL MEMREL('M01',WORK,KLLSOBEX_M01,KLLSOBEX_M01,KFREE,LFREE)
        END IF
      END IF
*
      IF(NPART.EQ.1.AND.NHOLE.EQ.0) THEN
        IF(IWAY.EQ.1) THEN
*
* Setup and allocate for M10
*
      CALL MEMGET('INTE',KLLSOBEX_M10,NINTER10,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX_M10,NINTER10,WORK,KFREE,LFREE)
*
      CALL IDIM_TCC_KRCC(WORK(KINTM10),NINTER10,ISYMINT,
     &     MX_ST_TSOSO_INT,MX_ST_TSOSO_BLK_INT,MX_TBLK_INT,
     &     WORK(KLLSOBEX_M10),WORK(KLIBSOBEX_M10),LEN_M10_VEC_CC,
     &     MX_SBSTR_INT,
     &     IDUM,IDUM,-1,IDUM,IDUM)
*
      CALL MEMGET('REAL',KM10,IMULTFAC*LEN_M10_VEC_CC,
     &            WORK,KFREE,LFREE)
      CALL SETVEC(WORK(KM10),IMULTFAC*ZERO,LEN_M10_VEC_CC)
*
* Make sure MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR are the largest possible
* (at least for now)
*
      IF(MX_ST_TSOSO_INT.GT.MX_ST_TSOSO) THEN
        MX_ST_TSOSO = MX_ST_TSOSO_INT
      END IF
      IF(MX_ST_TSOSO_BLK_INT.GT.MX_ST_TSOSO_BLK) THEN
        MX_ST_TSOSO_BLK = MX_ST_TSOSO_BLK_INT
      END IF
      IF(MX_TBLK_INT.GT.MX_TBLK) THEN
        MX_TBLK = MX_TBLK_INT
      END IF
      IF(MX_SBSTR_INT.GT.MX_SBSTR) THEN
        MX_SBSTR = MX_SBSTR_INT
      END IF
*
        ELSE
          CALL MEMREL('M10',WORK,KLLSOBEX_M10,KLLSOBEX_M10,KFREE,LFREE)
        END IF
      END IF
*
      IF(NPART.EQ.0.AND.NHOLE.EQ.2) THEN
        IF(IWAY.EQ.1) THEN
*
* Setup and allocate for M02
*
      CALL MEMGET('INTE',KLLSOBEX_M02,NINTER02,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX_M02,NINTER02,WORK,KFREE,LFREE)
*
      CALL IDIM_TCC_KRCC(WORK(KINTM02),NINTER02,ISYMINT,
     &     MX_ST_TSOSO_INT,MX_ST_TSOSO_BLK_INT,MX_TBLK_INT,
     &     WORK(KLLSOBEX_M02),WORK(KLIBSOBEX_M02),LEN_M02_VEC_CC,
     &     MX_SBSTR_INT,
     &     IDUM,IDUM,-1,IDUM,IDUM)
*
      CALL MEMGET('REAL',KM02,IMULTFAC*LEN_M02_VEC_CC,
     &            WORK,KFREE,LFREE)
      CALL SETVEC(WORK(KM02),IMULTFAC*ZERO,LEN_M02_VEC_CC)
*
* Make sure MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR are the largest possible
* (at least for now)
*
      IF(MX_ST_TSOSO_INT.GT.MX_ST_TSOSO) THEN
        MX_ST_TSOSO = MX_ST_TSOSO_INT
      END IF
      IF(MX_ST_TSOSO_BLK_INT.GT.MX_ST_TSOSO_BLK) THEN
        MX_ST_TSOSO_BLK = MX_ST_TSOSO_BLK_INT
      END IF
      IF(MX_TBLK_INT.GT.MX_TBLK) THEN
        MX_TBLK = MX_TBLK_INT
      END IF
      IF(MX_SBSTR_INT.GT.MX_SBSTR) THEN
        MX_SBSTR = MX_SBSTR_INT
      END IF
*
        ELSE
          CALL MEMREL('M02',WORK,KLLSOBEX_M02,KLLSOBEX_M02,KFREE,LFREE)
        END IF
      END IF
*
      IF(NPART.EQ.1.AND.NHOLE.EQ.1) THEN
        IF(IWAY.EQ.1) THEN
*
* Setup and allocate for M11
*
      CALL MEMGET('INTE',KLLSOBEX_M11,NINTER11,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX_M11,NINTER11,WORK,KFREE,LFREE)
*
      CALL IDIM_TCC_KRCC(WORK(KINTM11),NINTER11,ISYMINT,
     &     MX_ST_TSOSO_INT,MX_ST_TSOSO_BLK_INT,MX_TBLK_INT,
     &     WORK(KLLSOBEX_M11),WORK(KLIBSOBEX_M11),LEN_M11_VEC_CC,
     &     MX_SBSTR_INT,
     &     IDUM,IDUM,-1,IDUM,IDUM)
*
      CALL MEMGET('REAL',KM11,IMULTFAC*LEN_M11_VEC_CC,
     &            WORK,KFREE,LFREE)
      CALL SETVEC(WORK(KM11),IMULTFAC*ZERO,LEN_M11_VEC_CC)
*
* Make sure MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR are the largest possible
* (at least for now)
*
      IF(MX_ST_TSOSO_INT.GT.MX_ST_TSOSO) THEN
        MX_ST_TSOSO = MX_ST_TSOSO_INT
      END IF
      IF(MX_ST_TSOSO_BLK_INT.GT.MX_ST_TSOSO_BLK) THEN
        MX_ST_TSOSO_BLK = MX_ST_TSOSO_BLK_INT
      END IF
      IF(MX_TBLK_INT.GT.MX_TBLK) THEN
        MX_TBLK = MX_TBLK_INT
      END IF
      IF(MX_SBSTR_INT.GT.MX_SBSTR) THEN
        MX_SBSTR = MX_SBSTR_INT
      END IF
*
        ELSE
          CALL MEMREL('M11',WORK,KLLSOBEX_M11,KLLSOBEX_M11,KFREE,LFREE)
        END IF
      END IF
*
      IF(NPART.EQ.1.AND.NHOLE.EQ.2) THEN
        IF(IWAY.EQ.1) THEN
*
* Setup and allocate for M12
*
      CALL MEMGET('INTE',KLLSOBEX_M12,NINTER12,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KLIBSOBEX_M12,NINTER12,WORK,KFREE,LFREE)
*
      CALL IDIM_TCC_KRCC(WORK(KINTM12),NINTER12,ISYMINT,
     &     MX_ST_TSOSO_INT,MX_ST_TSOSO_BLK_INT,MX_TBLK_INT,
     &     WORK(KLLSOBEX_M12),WORK(KLIBSOBEX_M12),LEN_M12_VEC_CC,
     &     MX_SBSTR_INT,
     &     IDUM,IDUM,-1,IDUM,IDUM)
*
      CALL MEMGET('REAL',KM12,IMULTFAC*LEN_M12_VEC_CC,
     &            WORK,KFREE,LFREE)
      CALL SETVEC(WORK(KM12),IMULTFAC*ZERO,LEN_M12_VEC_CC)
*
* Make sure MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK,MX_SBSTR are the largest possible
* (at least for now)
*
      IF(MX_ST_TSOSO_INT.GT.MX_ST_TSOSO) THEN
        MX_ST_TSOSO = MX_ST_TSOSO_INT
      END IF
      IF(MX_ST_TSOSO_BLK_INT.GT.MX_ST_TSOSO_BLK) THEN
        MX_ST_TSOSO_BLK = MX_ST_TSOSO_BLK_INT
      END IF
      IF(MX_TBLK_INT.GT.MX_TBLK) THEN
        MX_TBLK = MX_TBLK_INT
      END IF
      IF(MX_SBSTR_INT.GT.MX_SBSTR) THEN
        MX_SBSTR = MX_SBSTR_INT
      END IF
*
        ELSE
          CALL MEMREL('M12',WORK,KLLSOBEX_M12,KLLSOBEX_M12,KFREE,LFREE)
        END IF
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'Dimension changes for intermediate with: NPART ',
     &              NPART,' and NHOLE ',NHOLE
        IF(MX_ST_TSOSO_INI.NE.MX_ST_TSOSO) THEN
          WRITE(6,*) 'Increased MX_ST_TSOSO from ',
     &                MX_ST_TSOSO_INI,' to ',MX_ST_TSOSO
        END IF
        IF(MX_ST_TSOSO_BLK_INI.NE.MX_ST_TSOSO_BLK) THEN
          WRITE(6,*) 'Increased MX_ST_TSOSO_BLK from ',
     &                MX_ST_TSOSO_BLK_INI,' to ',MX_ST_TSOSO_BLK
        END IF
        IF(MX_TBLK_INI.NE.MX_TBLK) THEN
          WRITE(6,*) 'Increased MX_TBLK from ',
     &                MX_TBLK_INI,' to ',MX_TBLK
        END IF
        IF(MX_SBSTR_INI.NE.MX_SBSTR) THEN
          WRITE(6,*) 'Increased  MX_SBSTR from ',
     &                MX_SBSTR_INI,' to ',MX_SBSTR
        END IF
      END IF
      CALL MEMCHK_KRCC(WORK)
*
      RETURN
      END
*
      SUBROUTINE M_TO_EX_DX_FORMAT(IWAY,IDIMEXDX,OUTT,NRANK,
     &                             IHSM,ITSM,
     &                             IFHM,IHINDEX,
     &                             MTYPE,TTYPE,RESTYPE,T,
     &                             WORK,KFREE,LFREE)
*
* Subroutine to map a given CAAB-type of cluster operator to a given
* EX,DX set of operators. Notice can be an injective mapping...hence we
* need to first find the dimension of the new array
*
* Strategy:
* 1. Find EX,DX
* 2. Find dimension of EX,DX
* 3. Perform mapping
* Lasse 2011
*
* Extended to also include the Hamiltonian. Lasse 2012
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "symm.inc"
*
      INTEGER MTYPE(4*NGAS),TTYPE(4*NGAS),RESTYPE(4*NGAS),IHINDEX(*)
      DIMENSION OUTT(IDIMEXDX),T(*),WORK(*)
* Local scratch
      INTEGER IEX(4*NGAS),IDX(4*NGAS)
      INTEGER ICONN(2*NRANK),IEXD1234(2*NRANK),IORD(4)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.0) THEN
          WRITE(6,*) ' Will dimension the needed  array '
        ELSE
          WRITE(6,*) ' Will perform the mapping '
        END IF
        IF(IFHM.EQ.1) THEN
          WRITE(6,*) ' Hamiltonian '
        ELSE
          WRITE(6,*) ' Intermediate '
        END IF
        CALL WRT_SPOX_TP_CC_KRCC(MTYPE,1)
        WRITE(6,*) ' Cluster operator '
        CALL WRT_SPOX_TP_CC_KRCC(TTYPE,1)
        WRITE(6,*) ' Resulting operator '
        CALL WRT_SPOX_TP_CC_KRCC(RESTYPE,1)
        WRITE(6,*) ' IFHM = ',IFHM
        WRITE(6,*) ' Intermediate size = ',IDIMEXDX
        WRITE(6,*) ' IWAY = ',IWAY
      END IF
      DO jj=1,100
       print*,'kkp'
      END DO
*
* First find EX,DX by calling old routines
*
      CALL CONTR_POS_REL_1T(MTYPE,TTYPE,RESTYPE,NGAS,ICONN,NRANK)
      NTOP = 1
      CALL HTYPE_TO_ED_KRCC(MTYPE,IHINDEX,NHOP,ICONN,NTOP,NRANK,
     &                   IDX,IEX,IEXD1234,IORD,ISIGN_DE)
      CALL MEMCHK_KRCC(WORK)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Excitation part '
        CALL WRT_SPOX_TP_CC_KRCC(IEX,1)
        WRITE(6,*) ' De-excitation part '
        CALL WRT_SPOX_TP_CC_KRCC(IDX,1)
      END IF
*
* Second find out which type of mapping is needed!
*
      IF(IFHM.EQ.2) THEN
        CALL MAP_TYPE_ORI_TO_EX_DX_SORT(IEX,IDX,NGAS,IMAP)
      END IF
*
* Finally find the dimension of the new array if IWAY = 0
* else perform the mapping
*
      IF(IFHM.EQ.2) THEN
        CALL IDIM_EX_DX(IWAY,IMAP,IEX,IDX,MTYPE,NRANK,
     &                  IHSM,OUTT,IDIMEXDX,T,
     &                  WORK,KFREE,LFREE)
      ELSE
        CALL OPCT1234M_REL_SETUP(IWAY,IEX,IDX,TTYPE,OUTT,IDIMEXDX,
     &                           ITSM,IHSM,IEXD1234,ISIGN_DE,IFHM,
     &                           T,NRANK,
     &                           WORK,KFREE,LFREE)
      END IF
*
      RETURN
      END
*
      SUBROUTINE OPCT1234M_REL_SETUP(IWAY,IOEX,IO1DX,IT1,OUTT,IDIM,
     &                               IT1SM,IOPSM,IOPINDX,ISIGN,IFHM,
     &                               INTEGRALS,NRANK,
     &                               WORK,KFREE,LFREE)
*
* Master routine for contracting operator with upto 4 operators
*
* Jeppe Olsen, April 2003 ( I hope it will be working before May ..)
*
* Reducing this to just one operator 
*
* Lasse 2011
*
* This routine is redundant and should be eliminated
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "cc_scr2.inc"
C
      DIMENSION WORK(*),INTEGRALS(*)
*
      CALL OPCT1234_REL_SETUP(IWAY,IOEX,IO1DX,IT1,IDIM,
     &    IT1SM,IOPSM,LCCB,
     &    WORK(KIX1_CA),WORK(KSX1_CA),WORK(KIX1_CB),WORK(KSX1_CB),
     &    WORK(KIX1_AA),WORK(KSX1_AA),WORK(KIX1_AB),WORK(KSX1_AB),
     &    WORK(KIX2_CA),WORK(KSX2_CA),WORK(KIX2_CB),WORK(KSX2_CB),
     &    WORK(KIX2_AA),WORK(KSX2_AA),WORK(KIX2_AB),WORK(KSX2_AB),
     &    WORK(KLOCK1),WORK(KLOCL1),
     &    WORK(KLSTOCC1),WORK(KLSTOCC2),WORK(KLSTOCC3),WORK(KLSTOCC4),
     &    WORK(KLTSCR1),WORK(KLTSCR2),
     &    OUTT,
     &    WORK(KLSMD1),
     &    WORK(KLSMEX),WORK(KLSMK1),WORK(KLSML1),
     &    WORK(KLNMD1),
     &    WORK(KLNMEX),WORK(KLNMK1),WORK(KLNML1),
     &    IOPINDX,NRANK,
     &    WORK(KLZ),WORK(KLZSCR),WORK(KLSTREO),
     &    ISIGN,IFHM,INTEGRALS,
     &    WORK,KFREE,LFREE)
*
      RETURN
      END
*
      SUBROUTINE OPCT1234_REL_SETUP(IWAY,IOEX,IO1DX,IT1,IDIM,
     &           IT1SM,IOPSM,LB,
     &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,
     &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,
     &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &           IOC_K1,IOC_L1,
     &           ISTR_D1,
     &           ISTR_EX,ISTR_K,ISTR_L,
     &           TSCR1,TSCR2,
     &           OPSCR,
     &           ISM_CAAB_D1,
     &           ISM_CAAB_EX,ISM_CAAB_K,ISM_CAAB_L,
     &           INM_CAAB_D1,
     &           INM_CAAB_EX,INM_CAAB_K,INM_CAAB_L,
     &           IEXD1234_INDX,NRANK,
     &           IZ,IZSCR,ISTREO,
     &           ISIGNG,IFHM,INTEGRALS,
     &           WORK,KFREE,LFREE)
*
* Will perform some of the setup and dimensioning otherwise done in
* OPCT1234_REL. Here the integrals will be batched and stored in batches
* for easy fetching in OPCT1234_REL. This will reduce the fetching time
* of integrals since same contractions are fetched in the same manner.
* Only the needed setup to find the integrals is in use.
*
* There should be some cleaning of this routine. Lasse 2012
*
*
* Contract indeces of Operator O with indeces of 
* excitation operator T1,T2, T3, T4
*
* Operator O is defined by an operator part IOEX, and 
* deexcitation parts IO1DX,IO2DX,IO3DX,IO4DX.
*
* LB   : Batch size for D1 and EX at the moment. 
*        May have to do this dynamic
*
* Jeppe Olsen, May 2000, Finished March 2003 ( Well it became July 2003)
*                        and even November 2004
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "ctcc.inc"
#include "symm.inc"
*. Input   
      INTEGER IOEX(NGAS,4),IO1DX(NGAS,4)
      INTEGER IT1(NGAS,4)
      DIMENSION INTEGRALS(*)
      DIMENSION WORK(*)
*. Map of index EXD1234 operator to original index of Hamiltonian
      INTEGER IEXD1234_INDX(2*NRANK)
*. Output
C     DIMENSION OT1234(*)
*
*. Local scratch
*
*. Occupation of conjugated operators
      INTEGER IO1DX_DAG(4*NGAS)
*
*. Occupation of gasspaces for various strings 
*. Number of strings per sym of the various string supergroups
      INTEGER NOEX(8,4),NO1DX(8,4)
      INTEGER NT1(8,4)
      INTEGER NK1(8,4)
      INTEGER NL1(8,4)
*. Offsets of strings with given sym for the strings of various CC ops
      INTEGER IBOEX(8,4),IBO1DX(8,4)
*. Offset in operators to strings with given sym
      INTEGER IBT1_TCC(8,8,8)
      INTEGER IBL1_TCC(8,8,8)
*
      INTEGER IB_D1K1(8,8,4) 
      INTEGER IB_EXK1(8,8,4)
*
*
*. Scratch through parameter list. 
*
*. IX1_* : Number of operators in operator * Largest C or A string
      DIMENSION IX1_CA(*),SX1_CA(*),IX1_CB(*),SX1_CB(*)
      DIMENSION IX1_AA(*),SX1_AA(*),IX1_AB(*),SX1_AB(*)
*. IX2, SX2
      DIMENSION IX2_CA(*),SX2_CA(*),IX2_CB(*),SX2_CB(*)
      DIMENSION IX2_AA(*),SX2_AA(*),IX2_AB(*),SX2_AB(*)
*
*. ISTR_D1, ISTR_K, ISTR_L : for occupations of strings  
      INTEGER ISTR_D1(*)
      INTEGER ISTR_EX(*)
      INTEGER ISTR_K(*), ISTR_L(*) 
*. For occupation of intermediate strings
      INTEGER IOC_K1(NGAS,4)
      INTEGER IOC_L1(NGAS,4)
*. For intermediates with both strings batched
*. Changed from squared to cubed by Lasse
      DIMENSION TSCR1(*),TSCR2(*)
*. For a batch of coefficients for Operator
*. Operator accessed as OP(EX,D1) so
      DIMENSION OPSCR(*)
*. For part of Hamiltonian
      INTEGER ISM_CAAB_D1(4,*), INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_EX(4,*), INM_CAAB_EX(4,*)
*. For CC operators
      INTEGER ISM_CAAB_K(4,*), INM_CAAB_K(4,*)
      INTEGER ISM_CAAB_L(4,*), INM_CAAB_L(4,*)
*
      INTEGER IZ(*), IZSCR(*), ISTREO(*)
*
* The story goes as
*
* Part 1 : O, T1 => OT1(l1) in TSCR4
*
*    Loop over symmetry of D1
*    Loop over batches  of D1
*
*    Loop over batches of Ex
*      Fetch O as O(d1,ex)
*      Loop over batches of K1
*        T1(I) => T1(d1,k1)
*        OT1(ex,k1) = sum(d1) O(d1,ex)T1(d1,k1)
*        OT1(ex,k1) => OT1(l1)
*      End of loop over batches of K1
*    End of loop over bathces of Ex
*    End of loop over batches of d1
*    End of loop over symmetry of D1
* ( We now have OT1(l1), all l1 in batch)    
*
*
* - That's all she wrote
*
*
      ZERO = 0.0D0
      ONE = 1.0D0
      IOFFM = 1
      IDIM = 0
*
      NTEST = 00 
      IF(NTEST.GE.100 )THEN
        WRITE(6,*) ' OPCT1234_REL_SETUP entered '
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' Ex Deex1 Deex2 Deex3 Deex4 ' 
        WRITE(6,*)
        CALL WRT_SPOX_TP_CC_KRCC(IOEX,1)
        CALL WRT_SPOX_TP_CC_KRCC(IO1DX,1)
        WRITE(6,*) 
        WRITE(6,*) ' Form of T1 T2 T3 T4'
        CALL WRT_SPOX_TP_CC_KRCC(IT1,1)
        CALL MEMCHK_KRCC(WORK)
*
      END IF
C     NTEST = 0
      SIGNG = DFLOAT(ISIGNG)
*. Symmetry of final type  
      IF(NTEST.GE.100) WRITE(6,*) ' IT1SM : ',IT1SM
C correct multiplication. Still to be named IL4SM? No IL1SM
      IL1SM = IDBGMULT(IOPSM,IT1SM) ! check sym 2011 correct
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Symmetry of OT1234 = ', IL1SM
      END IF
*
* =================================================
* Occupation of the various gas spaces for strings 
* =================================================
*
*. K1 = O1DX  I1 
      CALL NEW_CAAB_OC(IOC_K1,IT1,IO1DX,1,1,NGAS)
*. L1 = OEX K1 
      CALL NEW_CAAB_OC(IOC_L1,IOC_K1,IOEX,2,1,NGAS)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of K1, K2, K3, K4 '
        CALL WRT_SPOX_TP_CC_KRCC(IOC_K1,1)
        WRITE(6,*) ' Occupation of L1, L2, L3, L4 '
        CALL WRT_SPOX_TP_CC_KRCC(IOC_L1,1)
      END IF
*
* =============================================================
*. Obtain symmetry-dimensions and -offsets for various strings
* =============================================================
*
      DO I_CAAB = 1, 4
       call evenodd(IEO,I_CAAB)
       IF(IEO.EQ.1) THEN
         IUB = 1
       ELSE
         IUB = 2
       END IF
       IF(I_CAAB.GE.3) THEN
         ICA = 2
       ELSE
         ICA = 1
       END IF
C NST_SPGP_CC_REL should be checked!
C check with nst_spgrp_dbg from sigden.
       CALL NST_SPGP_CC_REL(IOEX(1,I_CAAB),NOEX(1,I_CAAB),IUB,ICA)
       CALL ZBASE(NOEX(1,I_CAAB),IBOEX(1,I_CAAB),NSMST)
*. D1
       CALL NST_SPGP_CC_REL(IO1DX(1,I_CAAB),NO1DX(1,I_CAAB),IUB,ICA)
       CALL ZBASE(NO1DX(1,I_CAAB),IBO1DX(1,I_CAAB),NSMST)
*. T1
       CALL NST_SPGP_CC_REL(IT1(1,I_CAAB),NT1(1,I_CAAB),IUB,ICA)
*. K1
       CALL NST_SPGP_CC_REL(IOC_K1(1,I_CAAB),NK1(1,I_CAAB),IUB,ICA)
*. L1   
       CALL NST_SPGP_CC_REL(IOC_L1(1,I_CAAB),NL1(1,I_CAAB),IUB,ICA)
*
      END DO
*. We now have the various dimensions, so we can write T1,T2,T3,T4 if 
*. required
*
* ============================
*. Offsets to T1, T2, T3, T4
* ============================
*
*. T1
      CALL Z_TCC_OFF2_REL(IBT1_TCC,LEN_T1,NT1(1,1),NT1(1,2),NT1(1,3),
     &                NT1(1,4),IT1SM,NSMST)
*
*
* ========================
*.  D1{\dagger} K1  => T1 mapping
* =========================
*
*. Obtain D1{\dagger}
      CALL CONJ_CAAB(IO1DX,IO1DX_DAG,NGAS,SP_D1K1)
      CALL T1T2_TO_T12_MAP_REL(IO1DX_DAG,IOC_K1,IT1,
     &     IX1_CA,SX1_CA,IX1_CB,SX1_CB,IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &     IB_D1K1,ISTR_D1,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D1K1,
     &     WORK,KFREE,LFREE)
*
* ======================
* Ex K1 => L1 mapping
* ======================
*
      CALL T1T2_TO_T12_MAP_REL(IOEX,IOC_K1,IOC_L1,
     &     IX2_CA,SX2_CA,IX2_CB,SX2_CB,IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &     IB_EXK1,ISTR_K,ISTR_L,ISTREO,IZ,IZSCR,SIGN_EXK1,
     &     WORK,KFREE,LFREE)
*
*. And the individual strings : all symmetries constructed
*. D1 strings 
*. Since ISTR_D1 is only for scratch space in T1T2_TO_T12_MAP_REL
      CALL STR_CAAB_REL(IO1DX,ISTR_D1,
     &                  WORK,KFREE,LFREE)
*. Ex strings 
      CALL STR_CAAB_REL(IOEX,ISTR_EX,
     &                  WORK,KFREE,LFREE)
*. Batch length : Should allow to vary 
      LD1_BAT = LB
      LEX_BAT = LB
*. Symmetry of Op T1, Op T1 T2, Op T1 T2 T3, Op T1 T2 T3 T4
*
* Check Op for annihilation symmetry
*
* IOPT1SM now final symmetry
      IOPT1SM    = IDBGMULT(IOPSM,IT1SM)
C see if I have to rename 2011
      IOPT1234SM = IOPT1SM
*
* Needs to be changed Lasse 2011
*
      L1SM = IOPT1SM ! I think
*. Number of L1 strings with sym L1SM
      NL1_TOT = LEN_TCCBLK_REL(NL1(1,1),NL1(1,2),NL1(1,3),
     &              NL1(1,4),L1SM,NSMST)
*. Offsets for L1 with sym L1SM
      CALL Z_TCC_OFF2_REL(IBL1_TCC,LEN_L1,NL1(1,1),
     &     NL1(1,2),NL1(1,3),NL1(1,4),L1SM,NSMST)
*
* Thinks this should be new start Lasse
*
* ==============================================================
* Part 1 O1(d1,ex) T1(i) => OT1(l1) in TSCR4
* ==============================================================
*
* Part 1 : Obtain in TSCR4 OT1(l1) for all l1 of correct sym
*
* Length of OT1(l1)
      LEN_OT1 = NL1_TOT!*ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
*
C     ZERO = 0.0D0
C     CALL SETVEC(WORK(KLTSCR4),ZERO,LEN_OT1)
*. Loop 1 can be realized in a number of different ways. 
*. The present approach minimizes the number of times 
*  integrals must be fetched 
*
      DO ID1SM = 1, NSMST
*. Obtain exsm from OP(ex,d4,d3,d2,d1) have symmetry IOPSM
        IF(NTEST.GE.100) WRITE(6,*) ' ID1SM = ', ID1SM
C
C Not sure about EXsym
C
C again symmetry 2011
C            ID1234SM   = IDBGMULT(ID1SM,ID234SM)
        ID1234SM = ID1SM
        IEXSM = IDBGMULT(IOPSM,IADJSYM(ID1234SM))
        K1SM = IDBGMULT(IT1SM,ID1SM)
*. Number of K1 strings with sym K1SM
        NK1_TOT = LEN_TCCBLK_REL(NK1(1,1),NK1(1,2),NK1(1,3), !!
     &            NK1(1,4),K1SM,NSMST)
        IF(NTEST.GE.100) WRITE(6,*) 'K1SM, NK1_TOT ,IEXSM = ', 
     &                               K1SM, NK1_TOT ,IEXSM
*. Number of excitations strings in O with this symmetry
        NEX_TOT = LEN_TCCBLK_REL(NOEX(1,1),NOEX(1,2),NOEX(1,3),
     &            NOEX(1,4),IEXSM,NSMST)
*. Number of D1 strings with given symmetry
        ND1_TOT = LEN_TCCBLK_REL(NO1DX(1,1),NO1DX(1,2),NO1DX(1,3),
     &         NO1DX(1,4),ID1SM,NSMST)
C     IF(ND1_TOT.GE.LB.OR.NEX_TOT.GE.LB) THEN
C       print*,'ND1_TOT,NEX_TOT,LB',ND1_TOT,NEX_TOT,LB
C       STOP 'batch'
C     END IF
        IF(ND1_TOT.GT.0.AND.NEX_TOT.GT.0.AND.NK1_TOT.GT.0) THEN
*. Number of batches of excitation part of O
C         print*,'ND1_TOT,NEX_TOT',ND1_TOT,NEX_TOT
          NEX_BAT = NEX_TOT/LEX_BAT
          IF(NEX_BAT*LEX_BAT.LT.NEX_TOT) NEX_BAT = NEX_BAT + 1
*. Number of D1 batches
          ND1_BAT = ND1_TOT/LD1_BAT
          IF(ND1_BAT*LD1_BAT.LT.ND1_TOT) ND1_BAT=ND1_BAT+1
*. Loop over D1 batches
          ID1SM_NEW = 1
          IF(NTEST.GE.100) WRITE(6,*) 'ND1_TOT, ND1_BAT ', 
     &                                 ND1_TOT, ND1_BAT
          DO ID1_BAT = 1, ND1_BAT
            IF(NTEST.GE.100) WRITE(6,*) ' ID1_BAT = ', ID1_BAT
            ID1_START = (ID1_BAT-1)*LD1_BAT + 1
            ID1_STOP  = MIN(ND1_TOT,ID1_START+LD1_BAT-1)
            ID1_BATLEN = ID1_STOP-ID1_START+1
            IF(NTEST.GE.100) 
     &      WRITE(6,*) ' ID1_STOP, ID1_START, ID1_BATLEN ',
     &                   ID1_STOP, ID1_START, ID1_BATLEN
*. Generate ID1 strings for given sym and batch
            CALL ISMNM_FOR_TCC_BAT_REL(NO1DX,ISM_CAAB_D1,INM_CAAB_D1,
     &           ID1SM,ID1_BATLEN,ID1SM_NEW,
     &           ISM_C1_D1,ISM_CA1_D1,ISM_AA1_D1,
     &           INM_AB1_D1,INM_AA1_D1,INM_CA1_D1,INM_CB1_D1,
     &           ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,      
     &           INM_ABINI_D1,INM_AAINI_D1,INM_CAINI_D1,INM_CBINI_D1,
     &           0)
            ID1SM_NEW = 0
*
* all excitation operators in one shot. 
* As we only have limited ISM, INM, the integrals (d1,d2,d3,d4,ex)
* are obtained in batches 
            IEXSM_NEW = 1
            DO IEX_BAT = 1, NEX_BAT
              IF(NTEST.GE.100) WRITE(6,*) ' IEX_BAT = ', IEX_BAT
              IEX_START = (IEX_BAT-1)*LB + 1
              IEX_STOP  = MIN(NEX_TOT,IEX_START+LB-1)
              IEX_BATLEN = IEX_STOP - IEX_START + 1
*. Generate IEX strings for given sym and batch
              CALL ISMNM_FOR_TCC_BAT_REL(NOEX,ISM_CAAB_EX,INM_CAAB_EX,
     &             IEXSM,IEX_BATLEN,IEXSM_NEW,
     &             ISM_C1_EX,ISM_CA1_EX,ISM_AA1_EX,
     &             INM_AB1_EX,INM_AA1_EX,INM_CA1_EX,INM_CB1_EX,
     &             ISM_CINI_EX,ISM_CAINI_EX,ISM_AAINI_EX,      
     &             INM_ABINI_EX,INM_AAINI_EX,INM_CAINI_EX,INM_CBINI_EX,
     &             0)
              IEXSM_NEW = 0
*. Length of EXD2D3D4
              LEN_EXD234 = IEX_BATLEN
*. Obtain integrals OP(D1,EX)
              IOFF = 1 + (IEX_START-1)*ID1_BATLEN
*
              IF(IWAY.EQ.1) THEN
              IOFF = IOFFM
              CALL GET_OPINT4_REL(OPSCR(IOFF),
     &           IO1DX,ID1_BATLEN,INM_CAAB_D1(1,1),ISM_CAAB_D1(1,1),
     &           ISTR_D1,IBO1DX,
     &           IOEX,IEX_BATLEN,INM_CAAB_EX(1,1),ISM_CAAB_EX(1,1),
     &           ISTR_EX,IBOEX,IEXD1234_INDX,IFHM,
     &           NRANK,INTEGRALS,IOFFM,
     &           WORK,KFREE,LFREE)
               END IF
*
               IDIM = IDIM + ID1_BATLEN*IEX_BATLEN
            END DO
*           ^ End of loop over batches of excitation operators OEX
*
            SIGNX = SP_D1K1*SIGN_D1K1
            SIGN_LK = SIGN_EXK1
            SIGNG = FLOAT(ISIGNG)
*
          END DO
*         ^ End of loop over D1 batch
        END IF
*       ^ End if there were nonvanishing dimensions
      END DO
*     ^ End of loop over symmetries of D1
*
      RETURN
      END
*
      SUBROUTINE IDIM_EX_DX_PRETTY(IWAY,IMAP,IEX,IDX,MTYPE,NRANK,
     &                      ISM,OUTT,IDIMEXDX,T,
     &                      WORK,KFREE,LFREE)
*
* Do as the routine name says. Dimension EX and DX.
* Based on ANA_GENCCS_CC_KRCC just with more loops.
* Original and pretty version (slow and not in use)
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "ctcc.inc"
*
      INTEGER IEX(4*NGAS),IDX(4*NGAS),MTYPE(4*NGAS)
      DIMENSION OUTT(IDIMEXDX),T(*),WORK(*)
* Local scratch 
      INTEGER IOCC_CAE(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_CBE(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AAE(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_ABE(NRANK*MX_ST_TSOSO_BLK_MX)
*
      INTEGER IOCC_CAD(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_CBD(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AAD(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_ABD(NRANK*MX_ST_TSOSO_BLK_MX)
*
      INTEGER IOCC_CA(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_CB(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AA(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AB(NRANK*MX_ST_TSOSO_BLK_MX)
*
      INTEGER IGRP_CAE(MXPNGAS),IGRP_CBE(MXPNGAS)
      INTEGER IGRP_AAE(MXPNGAS),IGRP_ABE(MXPNGAS)
*
      INTEGER IGRP_CAD(MXPNGAS),IGRP_CBD(MXPNGAS)
      INTEGER IGRP_AAD(MXPNGAS),IGRP_ABD(MXPNGAS)
*
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS)
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*
      INTEGER IFREEOP(4*NGAS),IGRP_FREE(MXPNGAS)
      INTEGER IOCC_FREE(MX_ST_TSOSO_BLK_MX)
*
      NTEST = 00
*
      IZERO = 0
      NORB = 0 !not needed
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Rank of EX+DX ',NRANK
        WRITE(6,*) ' Excitation operator '
        CALL WRT_SPOX_TP_CC_KRCC(IEX,1)
        WRITE(6,*) ' De-excitation operator '
        CALL WRT_SPOX_TP_CC_KRCC(IDX,1)
        WRITE(6,*) ' Original operator '
        CALL WRT_SPOX_TP_CC_KRCC(MTYPE,1)
        IF(IWAY.EQ.0) THEN
          WRITE(6,*) ' Will dimension the array '
        ELSE
          WRITE(6,*) ' Will copy from original to EX,DX format '
        END IF
      END IF
*
      CALL ISETVC(IFREEOP,IZERO,4*NGAS)
*. Transform from occupations to groups for EX 
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+0*NGAS),IGRP_CAE,1)
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+1*NGAS),IGRP_CBE,1)
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+2*NGAS),IGRP_AAE,1)
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+3*NGAS),IGRP_ABE,1)
*. Transform from occupations to groups for DX 
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+0*NGAS),IGRP_CAD,1)
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+1*NGAS),IGRP_CBD,1)
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+2*NGAS),IGRP_AAD,1)
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+3*NGAS),IGRP_ABD,1)
*. Transform from occupations to groups for MTYPE
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+0*NGAS),IGRP_CA,1)
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+1*NGAS),IGRP_CB,1)
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+2*NGAS),IGRP_AA,1)
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+3*NGAS),IGRP_AB,1)
*
      NEL_CAE = IELSUM(IEX(1+0*NGAS),NGAS)
      NEL_CBE = IELSUM(IEX(1+1*NGAS),NGAS)
      NEL_AAE = IELSUM(IEX(1+2*NGAS),NGAS)
      NEL_ABE = IELSUM(IEX(1+3*NGAS),NGAS)
*
      NEL_CAD = IELSUM(IDX(1+0*NGAS),NGAS)
      NEL_CBD = IELSUM(IDX(1+1*NGAS),NGAS)
      NEL_AAD = IELSUM(IDX(1+2*NGAS),NGAS)
      NEL_ABD = IELSUM(IDX(1+3*NGAS),NGAS)
*
      NEL_CA  = IELSUM(MTYPE(1+0*NGAS),NGAS)
      NEL_CB  = IELSUM(MTYPE(1+1*NGAS),NGAS)
      NEL_AA  = IELSUM(MTYPE(1+2*NGAS),NGAS)
      NEL_AB  = IELSUM(MTYPE(1+3*NGAS),NGAS)
*
      DO ISM_C = 1, NSMST
       ISM_A =  IDBGMULT(ISM,INVELM(ISM_C))
       ISM_A = IADJSYM(ISM_A)
* We now have the symmetry of both C and A for EX and DX since their
* total sym is ISM
       DO ISM_CA = 1, NSMST
        DO ISM_CAE = 1, NSMST
         ISM_CAD = IDBGMULT(ISM_CA,INVELM(ISM_CAE))
*
         ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
         DO ISM_CBE =1, NSMST
          ISM_CBD = IDBGMULT(ISM_CB,INVELM(ISM_CBE))
*
          DO ISM_AAA = 1, NSMST
           ISM_AA = IADJSYM(ISM_AAA)
           DO ISM_AAE = 1, NSMST
            ISM_AAD = IDBGMULT(ISM_AA,INVELM(ISM_AAE))
            ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
            DO ISM_ABE = 1, NSMST
             ISM_ABD = IDBGMULT(ISM_AB,INVELM(ISM_ABE))
*. obtain strings for EX 
             IUB = 1
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CAE,NGAS,
     &            ISM_CAE,NEL_CAE,NSTR_CAE,
     &            IOCC_CAE, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CBE,NGAS,
     &            ISM_CBE,NEL_CBE,NSTR_CBE,
     &            IOCC_CBE, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 1
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AAE,NGAS,
     &            ISM_AAE,NEL_AAE,NSTR_AAE,
     &            IOCC_AAE, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_ABE,NGAS,
     &            ISM_ABE,NEL_ABE,NSTR_ABE,
     &            IOCC_ABE, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
*. obtain strings for DX 
             IUB = 1
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CAD,NGAS,
     &            ISM_CAD,NEL_CAD,NSTR_CAD,
     &            IOCC_CAD, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CBD,NGAS,
     &            ISM_CBD,NEL_CBD,NSTR_CBD,
     &            IOCC_CBD, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 1
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AAD,NGAS,
     &            ISM_AAD,NEL_AAD,NSTR_AAD,
     &            IOCC_AAD, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_ABD,NGAS,
     &            ISM_ABD,NEL_ABD,NSTR_ABD,
     &            IOCC_ABD, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
*. obtain strings for MTYPE 
             IUB = 1
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CA,NGAS,
     &            ISM_CA,NEL_CA,NSTR_CA,
     &            IOCC_CA, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CB,NGAS,
     &            ISM_CB,NEL_CB,NSTR_CB,
     &            IOCC_CB, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 1
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AA,NGAS,
     &            ISM_AA,NEL_AA,NSTR_AA,
     &            IOCC_AA, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AB,NGAS,
     &            ISM_AB,NEL_AB,NSTR_AB,
     &            IOCC_AB, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
             NSTR_TOTE = NSTR_ABE*NSTR_AAE*NSTR_CBE*NSTR_CAE
             NSTR_TOTD = NSTR_ABD*NSTR_AAD*NSTR_CBD*NSTR_CAD
             NSTR_TOT  = NSTR_AB *NSTR_AA *NSTR_CB *NSTR_CA 
             NTOT = NSTR_TOTE*NSTR_TOTD*NSTR_TOT
             IF(NTOT.GE.1) THEN
* Find dimensions 
               IDIMEXDX = NSTR_TOTE * NSTR_TOTD
               I_ORI_DIM_M = NSTR_TOT
*
C             print*,'NSTR_TOTE,NSTR_TOTD,NSTR_TOT',
C    &                NSTR_TOTE,NSTR_TOTD,NSTR_TOT
C             print*,'NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE',
C    &                NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE
C             print*,'NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD',
C    &                NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD
C             print*,'NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA',
C    &                NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA
              IF(IWAY.EQ.1) THEN
              IF(IMAP.EQ.1) THEN
               CALL SORT_1(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.2) THEN
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,3)
               STOP ' Not yet implemented '
              ELSE IF(IMAP.EQ.3) THEN
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,4)
               STOP ' Not yet implemented '
              ELSE IF(IMAP.EQ.4) THEN
* Strategy: First reduce the EX operator to only the free index
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,1)
* Find the number of strings for IFREEOP with IDX symmetry
               CALL OCC_TO_GRP_CC_KRCC(IFREEOP(1+0*NGAS),IGRP_FREE,1)
               IUB = 1
               CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_FREE,NGAS,
     &              ISM_CAD,NEL_CAD,NSTR_FREE,
     &              IOCC_FREE, NORB,0,IDUM,IDUM,
     &              WORK,KFREE,LFREE)
               IF(NSTR_FREE.EQ.0) STOP ' something is wrong '
               print*,'NSTR_FREE',NSTR_FREE
               CALL SORT_4(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_FREE,NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.5) THEN
* Strategy: First reduce the EX operator to only the free index
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,2)
* Find the number of strings for IFREEOP with IDX symmetry
               CALL OCC_TO_GRP_CC_KRCC(IFREEOP(1+1*NGAS),IGRP_FREE,1)
               IUB = 2
               CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_FREE,NGAS,
     &              ISM_CBD,NEL_CBD,NSTR_FREE,
     &              IOCC_FREE, NORB,0,IDUM,IDUM,
     &              WORK,KFREE,LFREE)
               IF(NSTR_FREE.EQ.0) STOP ' something is wrong '
               print*,'NSTR_FREE',NSTR_FREE
               CALL SORT_5(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_FREE,NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.6) THEN
               CALL SORT_6(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.7) THEN
               CALL SORT_7(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.8) THEN
              ELSE IF(IMAP.EQ.9) THEN
              ELSE 
               STOP ' IMAP not set in IDIM_EX_DX'
              END IF
              END IF
             END IF
*
            END DO
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*     ^ End of loop over symmetry blocks
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of elements in EX,DX form: ',IDIMEXDX
        WRITE(6,*) ' Number of elements in original form: ',I_ORI_DIM_M
        IF(IDIMEXDX.EQ.I_ORI_DIM_M) THEN
          WRITE(6,*) ' We have a bijective mapping '
        ELSE
          WRITE(6,*) ' We have an injective mapping '
        END IF 
        IF(IWAY.NE.0) THEN
          WRITE(6,*) ' The original order (elements only) '
          CALL WRTMAT(T,1,I_ORI_DIM_M,1,I_ORI_DIM_M)
          WRITE(6,*) ' The new EX,DX order '
          CALL WRTMAT(OUTT,1,IDIMEXDX,1,IDIMEXDX)
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE IDIM_EX_DX(IWAY,IMAP,IEX,IDX,MTYPE,NRANK,
     &                      ISM,OUTT,IDIMEXDX,T,
     &                      WORK,KFREE,LFREE)
*
* Do as the routine name says. Dimension EX and DX.
* Based on ANA_GENCCS_CC_KRCC just with more loops.
* Moved GETSTR2_TOTSM_SPGP_KRCC into loop structure for speed up and
* added a lot of IF statements. Not so pretty but faster.
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "ctcc.inc"
*
      INTEGER IEX(4*NGAS),IDX(4*NGAS),MTYPE(4*NGAS)
      DIMENSION OUTT(IDIMEXDX),T(*),WORK(*)
* Local scratch 
      INTEGER IOCC_CAE(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_CBE(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AAE(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_ABE(NRANK*MX_ST_TSOSO_BLK_MX)
*
      INTEGER IOCC_CAD(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_CBD(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AAD(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_ABD(NRANK*MX_ST_TSOSO_BLK_MX)
*
      INTEGER IOCC_CA(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_CB(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AA(NRANK*MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AB(NRANK*MX_ST_TSOSO_BLK_MX)
*
      INTEGER IGRP_CAE(MXPNGAS),IGRP_CBE(MXPNGAS)
      INTEGER IGRP_AAE(MXPNGAS),IGRP_ABE(MXPNGAS)
*
      INTEGER IGRP_CAD(MXPNGAS),IGRP_CBD(MXPNGAS)
      INTEGER IGRP_AAD(MXPNGAS),IGRP_ABD(MXPNGAS)
*
      INTEGER IGRP_CA(MXPNGAS),IGRP_CB(MXPNGAS)
      INTEGER IGRP_AA(MXPNGAS),IGRP_AB(MXPNGAS)
*
      INTEGER IFREEOP(4*NGAS),IGRP_FREE(MXPNGAS)
      INTEGER IOCC_FREE(MX_ST_TSOSO_BLK_MX)
*
      NTEST = 00
*
      IZERO = 0
      NORB = 0 !not needed
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Rank of EX+DX ',NRANK
        WRITE(6,*) ' Excitation operator '
        CALL WRT_SPOX_TP_CC_KRCC(IEX,1)
        WRITE(6,*) ' De-excitation operator '
        CALL WRT_SPOX_TP_CC_KRCC(IDX,1)
        WRITE(6,*) ' Original operator '
        CALL WRT_SPOX_TP_CC_KRCC(MTYPE,1)
        IF(IWAY.EQ.0) THEN
          WRITE(6,*) ' Will dimension the array '
        ELSE
          WRITE(6,*) ' Will copy from original to EX,DX format '
        END IF
      END IF
*
      CALL ISETVC(IFREEOP,IZERO,4*NGAS)
*. Transform from occupations to groups for EX 
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+0*NGAS),IGRP_CAE,1)
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+1*NGAS),IGRP_CBE,1)
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+2*NGAS),IGRP_AAE,1)
      CALL OCC_TO_GRP_CC_KRCC(IEX(1+3*NGAS),IGRP_ABE,1)
*. Transform from occupations to groups for DX 
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+0*NGAS),IGRP_CAD,1)
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+1*NGAS),IGRP_CBD,1)
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+2*NGAS),IGRP_AAD,1)
      CALL OCC_TO_GRP_CC_KRCC(IDX(1+3*NGAS),IGRP_ABD,1)
*. Transform from occupations to groups for MTYPE
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+0*NGAS),IGRP_CA,1)
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+1*NGAS),IGRP_CB,1)
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+2*NGAS),IGRP_AA,1)
      CALL OCC_TO_GRP_CC_KRCC(MTYPE(1+3*NGAS),IGRP_AB,1)
*
      NEL_CAE = IELSUM(IEX(1+0*NGAS),NGAS)
      NEL_CBE = IELSUM(IEX(1+1*NGAS),NGAS)
      NEL_AAE = IELSUM(IEX(1+2*NGAS),NGAS)
      NEL_ABE = IELSUM(IEX(1+3*NGAS),NGAS)
*
      NEL_CAD = IELSUM(IDX(1+0*NGAS),NGAS)
      NEL_CBD = IELSUM(IDX(1+1*NGAS),NGAS)
      NEL_AAD = IELSUM(IDX(1+2*NGAS),NGAS)
      NEL_ABD = IELSUM(IDX(1+3*NGAS),NGAS)
*
      NEL_CA  = IELSUM(MTYPE(1+0*NGAS),NGAS)
      NEL_CB  = IELSUM(MTYPE(1+1*NGAS),NGAS)
      NEL_AA  = IELSUM(MTYPE(1+2*NGAS),NGAS)
      NEL_AB  = IELSUM(MTYPE(1+3*NGAS),NGAS)
*
      DO ISM_C = 1, NSMST
       ISM_A =  IDBGMULT(ISM,INVELM(ISM_C))
       ISM_A = IADJSYM(ISM_A)
* We now have the symmetry of both C and A for EX and DX since their
* total sym is ISM
       DO ISM_CA = 1, NSMST
*. obtain strings for MTYPE CA
        IUB = 1
        CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CA,NGAS,
     &       ISM_CA,NEL_CA,NSTR_CA,
     &       IOCC_CA, NORB,0,IDUM,IDUM,
     &       WORK,KFREE,LFREE)
        IF(NSTR_CA.EQ.0) CYCLE
        ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
*. obtain strings for MTYPE CB
        IUB = 2
        CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CB,NGAS,
     &       ISM_CB,NEL_CB,NSTR_CB,
     &       IOCC_CB, NORB,0,IDUM,IDUM,
     &       WORK,KFREE,LFREE)
        IF(NSTR_CB.EQ.0) CYCLE
        DO ISM_CAE = 1, NSMST
*. obtain strings for EX CAE
         IUB = 1
         CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CAE,NGAS,
     &         ISM_CAE,NEL_CAE,NSTR_CAE,
     &         IOCC_CAE, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
         IF(NSTR_CAE.EQ.0) CYCLE
         ISM_CAD = IDBGMULT(ISM_CA,INVELM(ISM_CAE))
*. obtain strings for DX CAD
         IUB = 1
         CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CAD,NGAS,
     &        ISM_CAD,NEL_CAD,NSTR_CAD,
     &        IOCC_CAD, NORB,0,IDUM,IDUM,
     &        WORK,KFREE,LFREE)
         IF(NSTR_CAD.EQ.0) CYCLE
*
         DO ISM_CBE =1, NSMST
*. obtain strings for EX CBE
          IUB = 2
          CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CBE,NGAS,
     &         ISM_CBE,NEL_CBE,NSTR_CBE,
     &         IOCC_CBE, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
          IF(NSTR_CBE.EQ.0) CYCLE
          ISM_CBD = IDBGMULT(ISM_CB,INVELM(ISM_CBE))
*. obtain strings for DX CBD
          IUB = 2
          CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_CBD,NGAS,
     &         ISM_CBD,NEL_CBD,NSTR_CBD,
     &         IOCC_CBD, NORB,0,IDUM,IDUM,
     &         WORK,KFREE,LFREE)
          IF(NSTR_CBD.EQ.0) CYCLE
*
          DO ISM_AAA = 1, NSMST
           ISM_AA = IADJSYM(ISM_AAA)
*. obtain strings for MTYPE AA
           IUB = 1
           CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AA,NGAS,
     &          ISM_AA,NEL_AA,NSTR_AA,
     &          IOCC_AA, NORB,0,IDUM,IDUM,
     &          WORK,KFREE,LFREE)
           IF(NSTR_AA.EQ.0) CYCLE
           ISM_AB = IDBGMULT(ISM_A,INVELM(ISM_AA))
*. obtain strings for MTYPE 
           IUB = 2
           CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AB,NGAS,
     &          ISM_AB,NEL_AB,NSTR_AB,
     &          IOCC_AB, NORB,0,IDUM,IDUM,
     &          WORK,KFREE,LFREE)
           IF(NSTR_AB.EQ.0) CYCLE
           DO ISM_AAE = 1, NSMST
*. obtain strings for EX AAE
            IUB = 1
            CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AAE,NGAS,
     &           ISM_AAE,NEL_AAE,NSTR_AAE,
     &           IOCC_AAE, NORB,0,IDUM,IDUM,
     &           WORK,KFREE,LFREE)
            IF(NSTR_AAE.EQ.0) CYCLE
            ISM_AAD = IDBGMULT(ISM_AA,INVELM(ISM_AAE))
*. obtain strings for DX 
            IUB = 1
            CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_AAD,NGAS,
     &           ISM_AAD,NEL_AAD,NSTR_AAD,
     &           IOCC_AAD, NORB,0,IDUM,IDUM,
     &           WORK,KFREE,LFREE)
            IF(NSTR_AAD.EQ.0) CYCLE
            DO ISM_ABE = 1, NSMST
*. obtain strings for EX 
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_ABE,NGAS,
     &            ISM_ABE,NEL_ABE,NSTR_ABE,
     &            IOCC_ABE, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IF(NSTR_ABE.EQ.0) CYCLE
             ISM_ABD = IDBGMULT(ISM_AB,INVELM(ISM_ABE))
*. obtain strings for DX 
             IUB = 2
             CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_ABD,NGAS,
     &            ISM_ABD,NEL_ABD,NSTR_ABD,
     &            IOCC_ABD, NORB,0,IDUM,IDUM,
     &            WORK,KFREE,LFREE)
             IF(NSTR_ABD.EQ.0) CYCLE
*. Loop over T elements as  matrix T(I_CA, I_CB, IAA, I_AB)
             NSTR_TOTE = NSTR_ABE*NSTR_AAE*NSTR_CBE*NSTR_CAE
             NSTR_TOTD = NSTR_ABD*NSTR_AAD*NSTR_CBD*NSTR_CAD
             NSTR_TOT  = NSTR_AB *NSTR_AA *NSTR_CB *NSTR_CA 
             NTOT = NSTR_TOTE*NSTR_TOTD*NSTR_TOT
             IF(NTOT.GE.1) THEN
* Find dimensions 
               IDIMEXDX = NSTR_TOTE * NSTR_TOTD
               I_ORI_DIM_M = NSTR_TOT
*
C             print*,'NSTR_TOTE,NSTR_TOTD,NSTR_TOT',
C    &                NSTR_TOTE,NSTR_TOTD,NSTR_TOT
C             print*,'NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE',
C    &                NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE
C             print*,'NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD',
C    &                NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD
C             print*,'NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA',
C    &                NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA
              IF(IWAY.EQ.1) THEN
              IF(IMAP.EQ.1) THEN
C              IF(IFHM.EQ.1) THEN
C              CALL SORT_H(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
C    &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
C    &                     NEL_ABE,NEL_AAE,NEL_CBE,NEL_CAE,
C    &                     NEL_ABD,NEL_AAD,NEL_CBD,NEL_CAD,
C    &                     IOCC_ABE,IOCC_AAE,IOCC_CBE,IOCC_CAE,
C    &                     IOCC_ABD,IOCC_AAD,IOCC_CBD,IOCC_CAD,
C    &                     NSTR_TOT,IEXD1234,NRANK,ISIGN_DE,T,OUTT)
C              ELSE
               CALL SORT_1(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_TOT,T,OUTT)
C              END IF
              ELSE IF(IMAP.EQ.2) THEN
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,3)
               STOP ' Not yet implemented '
              ELSE IF(IMAP.EQ.3) THEN
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,4)
               STOP ' Not yet implemented '
              ELSE IF(IMAP.EQ.4) THEN
* Strategy: First reduce the EX operator to only the free index
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,1)
* Find the number of strings for IFREEOP with IDX symmetry
               CALL OCC_TO_GRP_CC_KRCC(IFREEOP(1+0*NGAS),IGRP_FREE,1)
               IUB = 1
               CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_FREE,NGAS,
     &              ISM_CAD,NEL_CAD,NSTR_FREE,
     &              IOCC_FREE, NORB,0,IDUM,IDUM,
     &              WORK,KFREE,LFREE)
               IF(NSTR_FREE.EQ.0) STOP ' something is wrong '
               print*,'NSTR_FREE',NSTR_FREE
               CALL SORT_4(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_FREE,NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.5) THEN
* Strategy: First reduce the EX operator to only the free index
               CALL REDUCE_OP(IEX,IDX,IFREEOP,NGAS,2)
* Find the number of strings for IFREEOP with IDX symmetry
               CALL OCC_TO_GRP_CC_KRCC(IFREEOP(1+1*NGAS),IGRP_FREE,1)
               IUB = 2
               CALL GETSTR2_TOTSM_SPGP_KRCC(IUB,IGRP_FREE,NGAS,
     &              ISM_CBD,NEL_CBD,NSTR_FREE,
     &              IOCC_FREE, NORB,0,IDUM,IDUM,
     &              WORK,KFREE,LFREE)
               IF(NSTR_FREE.EQ.0) STOP ' something is wrong '
               print*,'NSTR_FREE',NSTR_FREE
               CALL SORT_5(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_FREE,NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.6) THEN
               CALL SORT_6(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.7) THEN
               CALL SORT_7(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                     NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                     NSTR_TOT,T,OUTT)
              ELSE IF(IMAP.EQ.8) THEN
               STOP ' Not yet implemented '
              ELSE IF(IMAP.EQ.9) THEN
               STOP ' Not yet implemented '
              ELSE 
               STOP ' IMAP not set in IDIM_EX_DX'
              END IF
              END IF
             END IF
*
            END DO
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*     ^ End of loop over symmetry blocks
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Number of elements in EX,DX form: ',IDIMEXDX
        WRITE(6,*) ' Number of elements in original form: ',I_ORI_DIM_M
        IF(IDIMEXDX.EQ.I_ORI_DIM_M) THEN
          WRITE(6,*) ' We have a bijective mapping '
        ELSE
          WRITE(6,*) ' We have an injective mapping '
        END IF 
        IF(IWAY.NE.0) THEN
          WRITE(6,*) ' The original order (elements only) '
          CALL WRTMAT(T,1,I_ORI_DIM_M,1,I_ORI_DIM_M)
          WRITE(6,*) ' The new EX,DX order '
          CALL WRTMAT(OUTT,1,IDIMEXDX,1,IDIMEXDX)
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_H(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NEL_ABE,NEL_AAE,NEL_CBE,NEL_CAE,
     &                  NEL_ABD,NEL_AAD,NEL_CBD,NEL_CAD,
     &                  IOCC_ABE,IOCC_AAE,IOCC_CBE,IOCC_CAE,
     &                  IOCC_ABD,IOCC_AAD,IOCC_CBD,IOCC_CAD,
     &                  NSTR_TOT,IEXD1234,NRANK,ISIGN_DE,T,OUTT)
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
      INTEGER IOCC_ABE(*),IOCC_AAE(*),IOCC_CBE(*),IOCC_CAE(*)
      INTEGER IOCC_ABD(*),IOCC_AAD(*),IOCC_CBD(*),IOCC_CAD(*)
      INTEGER IEXD1234(2*NRANK)
* Local scratch
      INTEGER IOFF(4)
*
      NTEST = 00
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine for the Hamiltonian '
      END IF
*
      ICOUNT = 0
* Find the offset
      IF(NEL_CAE.GE.1) THEN
        CALL ICOPVE(IOCC_CAE,IOFF,NDIM)
      END IF
*
      DO ISTR_CAE = 1,NSTR_CAE
      DO ISTR_ABE = 1,NSTR_ABE
       DO ISTR_AAE = 1,NSTR_AAE
        DO ISTR_CBE = 1,NSTR_CBE
          DO ISTR_ABD = 1,NSTR_ABD
           DO ISTR_AAD = 1,NSTR_AAD
            DO ISTR_CBD = 1,NSTR_CBD
             DO ISTR_CAD = 1,NSTR_CAD
              ICOUNT = ICOUNT + 1
C             print*,'IOFF_CA,ICOUNT',IOFF_CA,ICOUNT
C             print*,'---------------------------'
              OUTT(ICOUNT) = T(IOFF_CA)
             END DO
            END DO
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_1(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NSTR_TOT,T,OUTT)
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
*
      NTEST = 00
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine number 1 '
      END IF
*
      ICOUNT = 0
*
      I = NSTR_CAE*NSTR_CAD
      J = I*NSTR_CBE*NSTR_CBD
      K = J*NSTR_AAE*NSTR_AAD
      L = J*NSTR_AAD
      M = K*NSTR_ABD
C     print*,'I,J,K,L,M',I,J,K,L,M
*
      DO ISTR_ABE = 1,NSTR_ABE
C      print*,'ISTR_ABE',ISTR_ABE
        IOFF_ABE = (ISTR_ABE-1)*M
       DO ISTR_AAE = 1,NSTR_AAE
        IOFF_AAE = (ISTR_AAE-1)*L
C       print*,'ISTR_AAE',ISTR_AAE
        DO ISTR_CBE = 1,NSTR_CBE
C        print*,'ISTR_CBE',ISTR_CBE
         DO ISTR_CAE = 1,NSTR_CAE
C         print*,'ISTR_CAE',ISTR_CAE
          DO ISTR_ABD = 1,NSTR_ABD
           IOFF_AB = ISTR_ABD 
           IOFF_AB = (IOFF_AB-1)*K !IS THIS CORRECT!
C          print*,'IOFF_AB',IOFF_AB
           DO ISTR_AAD = 1,NSTR_AAD
            IOFF_AA = ISTR_AAD
            IOFF_AA = (IOFF_AA-1)*J
C           print*,'IOFF_AA',IOFF_AA
            DO ISTR_CBD = 1,NSTR_CBD
             IOFF_CB = (NSTR_CBE * (ISTR_CBD - 1) + ISTR_CBE)
             IOFF_CB = (IOFF_CB-1)*I 
C            print*,'IOFF_CB',IOFF_CB
             DO ISTR_CAD = 1,NSTR_CAD
              IOFF_CA = (NSTR_CAE * (ISTR_CAD- 1)+ISTR_CAE)
              IOFF_CA=IOFF_CA+IOFF_CB+IOFF_AA+IOFF_AB+IOFF_AAE+IOFF_ABE
              ICOUNT = ICOUNT + 1
C             print*,'IOFF_CA,ICOUNT',IOFF_CA,ICOUNT
C             print*,'---------------------------'
              OUTT(ICOUNT) = T(IOFF_CA)
             END DO
            END DO
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE MAP_TYPE_ORI_TO_EX_DX_SORT(IEX,IDX,NGAS,IMAP)
*
* Will identify the sort of mapping needed and the way the resorting
* from original to EX,DX format. Assumes that never more than two
* indices in an intermediate are contracted.
* There are three distinct cases resulting in 9 different sortings: 
* First case : The index (or indicies) contracted are the only one of
*              it's kind i.e. there are no more possible indices of this
*              CAAB kind left in the resulting operator => IMAP = 1
* Second case : Only one of two possible indices are contracted will
*               give two seperate cases:
*               If the two indices are in a different GAS then :
*               If the contracted index are the first index then:
*               For a hole => IMAP = 1
*               For a particle => IMAP = 2,3 (Not yet implemented)
*               If the contracted index are the second index then:
*               For a hole => IMAP = 4,5
*               For a particle => IMAP = 1 
*               If the indices belongs to the same GAS then we do not
*               have a bijective mapping but an injective mapping. 
*               => IMAP = 6-9 (8-9 Not yet implemented)
* Third case : Two indices are contracted => IMAP = 1
*
#include "implicit.inc"
*
      INTEGER IEX(NGAS,4),IDX(NGAS,4)
*
      NTEST = 00
      IMAP = 0
*
* First find the index (indices) to be contracted
*
      NCA = 0
      NCB = 0
      NAA = 0
      NAB = 0
      DO IGAS = 1,NGAS
        NCA = NCA + IDX(IGAS,1)
        NCB = NCB + IDX(IGAS,2)
        NAA = NAA + IDX(IGAS,3)
        NAB = NAB + IDX(IGAS,4)
      END DO
*
* Two indices then IMAP = 1
* 
      NINDEX = NCA + NCB + NAA + NAB
c     print*,'NINDEX,NCA,NCB,NAA,NAB',NINDEX,NCA,NCB,NAA,NAB
      IF(NINDEX.EQ.2) THEN
        IMAP = 1
        RETURN
      END IF
*
* Look for order of indices
*
      NC = NCA + NCB
      IF(NC.EQ.1) THEN
        IF(NCA.EQ.1) THEN
          IROW = 1
        ELSE
          IROW = 2
        END IF
        IEXCOUNT = 0
        DO IGAS = 1,NGAS
          IF(IEX(IGAS,IROW).GE.1) IEXCOUNT = IEXCOUNT + 1
          IF(IDX(IGAS,IROW).GE.1) THEN
            IF(IEX(IGAS,IROW).EQ.0) THEN
              IF(IEXCOUNT.EQ.0) THEN
                IMAP = 1
              ELSE
                IF(NCA.EQ.1) THEN
                  IMAP = 4
                ELSE
                  IMAP = 5
                END IF
              END IF
            ELSE
              IF(NCA.EQ.1) THEN      
                IMAP = 6
              ELSE
                IMAP = 7
              END IF
            END IF
            EXIT
          END IF
        END DO
      END IF
      NA = NAA + NAB
      IF(NA.EQ.1) THEN
        IF(NAA.EQ.1) THEN
          IROW = 3
        ELSE
          IROW = 4
        END IF
        IEXCOUNT = 0
        DO IGAS = NGAS,1,-1
          IF(IEX(IGAS,IROW).GE.1) IEXCOUNT = IEXCOUNT + 1
          IF(IDX(IGAS,IROW).GE.1) THEN
            IF(IEX(IGAS,IROW).EQ.0) THEN
              IF(IEXCOUNT.EQ.0) THEN
                IMAP = 1
              ELSE
                IF(NAA.EQ.1) THEN
                  IMAP = 2
                ELSE 
                  IMAP = 3
                END IF
              END IF
            ELSE
              IF(NAA.EQ.1) THEN      
                IMAP = 8
              ELSE
                IMAP = 9
              END IF
            END IF
            EXIT
          END IF
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Type of mapping ',IMAP
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_4(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NSTR_FREE,NSTR_TOT,T,OUTT)
*
* This is for a free CA index in EX that is lower than the one in DX
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
*
      NTEST = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine number 4 '
      END IF
*
* Number of strings after 
*
      NSTR_CAE_REMAIN = NSTR_CAE/NSTR_FREE ! Yes this is an integer
*
      ICOUNT = 0
*
      II = NSTR_CAE_REMAIN*NSTR_CAD
      I = NSTR_CAE*NSTR_CAD
      J = I*NSTR_CBE*NSTR_CBD
      K = J*NSTR_AAE*NSTR_AAD
      L = J*NSTR_AAD
      M = K*NSTR_ABD
C     print*,'I,J,K,L,M',I,J,K,L,M
*
      DO ISTR_ABE = 1,NSTR_ABE
C      print*,'ISTR_ABE',ISTR_ABE
        IOFF_ABE = (ISTR_ABE-1)*M
       DO ISTR_AAE = 1,NSTR_AAE
        IOFF_AAE = (ISTR_AAE-1)*L
C       print*,'ISTR_AAE',ISTR_AAE
        DO ISTR_CBE = 1,NSTR_CBE
         IOFF_CBE = (ISTR_CBE-1)*I 
C        print*,'ISTR_CBE',ISTR_CBE
!        DO ISTR_CAE = 1,NSTR_CAE
         DO ISTR_CAE = 1,NSTR_FREE
          IOFF_CAE = (ISTR_CAE-1)*II
C         print*,'ISTR_CAE',ISTR_CAE
          DO ISTR_CAR = 1,NSTR_CAE_REMAIN
* Dx start
           DO ISTR_CAD = 1,NSTR_CAD
!           IOFF_CA = (NSTR_CAE * (ISTR_CAD- 1)+ISTR_CAE)
            IOFF_CA = (NSTR_CAE_REMAIN * (ISTR_CAD- 1)+ISTR_CAR)
            IOFF_CA=IOFF_CA+IOFF_CAE+IOFF_CBE+IOFF_AAE+IOFF_ABE
            ICOUNT = ICOUNT + 1
C           print*,'IOFF_CA,ICOUNT',IOFF_CA,ICOUNT
C           print*,'---------------------------'
            OUTT(ICOUNT) = T(IOFF_CA)
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_5(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NSTR_FREE,NSTR_TOT,T,OUTT)
*
* This is for a free CA index in EX that is lower than the one in DX
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
*
      NTEST = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine number 5 '
      END IF
*
* Number of strings after 
*
      NSTR_CBE_REMAIN = NSTR_CBE/NSTR_FREE ! Yes this is an integer
*
      ICOUNT = 0
*
C     I = NSTR_CAE*NSTR_CAD
      I = NSTR_CAE*NSTR_CAD*NSTR_CBE_REMAIN*NSTR_CBD
      J = I*NSTR_CBE!*NSTR_CBD
      K = J*NSTR_AAE*NSTR_AAD
      L = J*NSTR_AAD
      M = K*NSTR_ABD
C     print*,'I,J,K,L,M',I,J,K,L,M
*
      DO ISTR_ABE = 1,NSTR_ABE
C      print*,'ISTR_ABE',ISTR_ABE
        IOFF_ABE = (ISTR_ABE-1)*M
       DO ISTR_AAE = 1,NSTR_AAE
        IOFF_AAE = (ISTR_AAE-1)*L
C       print*,'ISTR_AAE',ISTR_AAE
        DO ISTR_CBE = 1,NSTR_FREE
         IOFF_CBE = (ISTR_CBE-1)*I 
C        print*,'ISTR_CBE',ISTR_CBE
         DO ISTR_CBR = 1,NSTR_CBE_REMAIN
          DO ISTR_CAE = 1,NSTR_CAE
           IOFF_CAE = ISTR_CAE-1
C          print*,'ISTR_CAE',ISTR_CAE
* Dx start
           DO ISTR_CBD = 1,NSTR_CBD
!           IOFF_CA = (NSTR_CAE * (ISTR_CAD- 1)+ISTR_CAE)
            IOFF_CB = (NSTR_CBE_REMAIN * (ISTR_CBD- 1)+ISTR_CBR)
            IOFF_CB=IOFF_CB+IOFF_CAE+IOFF_CBE+IOFF_AAE+IOFF_ABE
            ICOUNT = ICOUNT + 1
C           print*,'IOFF_CB,ICOUNT',IOFF_CB,ICOUNT
C           print*,'---------------------------'
            OUTT(ICOUNT) = T(IOFF_CB)
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE REDUCE_OP(IEX,IDX,IFREEOP,NGAS,INDX)
*
* Reduces the EX operator to the free index
*
#include "implicit.inc"
*
      INTEGER IEX(NGAS,4),IDX(NGAS,4),IFREEOP(NGAS,4)
*
      NTEST = 00
*
      IF(INDX.EQ.1.OR.INDX.EQ.2) THEN
        DO IGAS =1,NGAS
          IF(IDX(IGAS,INDX).EQ.0) THEN
            IFREEOP(IGAS,INDX) = IEX(IGAS,INDX)
          ELSE
            EXIT
          END IF
        END DO
      ELSE IF(INDX.EQ.3.OR.INDX.EQ.4) THEN
        DO IGAS =NGAS,1,-1
          IF(IDX(IGAS,INDX).EQ.0) THEN
            IFREEOP(IGAS,INDX) = IEX(IGAS,INDX)
          ELSE
            EXIT
          END IF
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN
        IF(INDX.EQ.1) THEN
          WRITE(6,*) ' A creator alpha free index '
        ELSE IF(INDX.EQ.2) THEN
          WRITE(6,*) ' A creator beta free index '
        ELSE IF(INDX.EQ.3) THEN
          WRITE(6,*) ' An annihilator alpha free index '
        ELSE IF(INDX.EQ.4) THEN
          WRITE(6,*) ' An annihilator beta free index '
        END IF
        WRITE(6,*) ' The DX operator '
        CALL WRT_SPOX_TP_CC_KRCC(IDX,1)
        WRITE(6,*) ' The EX operator '
        CALL WRT_SPOX_TP_CC_KRCC(IEX,1)
        WRITE(6,*) ' The de-excitation part of EX '
        CALL WRT_SPOX_TP_CC_KRCC(IFREEOP,1)
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_62(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NSTR_TOT,T,OUTT)
* A save routine
*
* This is for a free CA index in EX that belongs to the same GAS as the
* free index in DX
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
* Local scratch
      INTEGER ISET(NSTR_CAD,NSTR_CAD)
*
      NTEST = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine number 6 '
      END IF
*
* Calculate the length of remaining indices of CA
*
      NSTR_CAE_REMAIN = NSTR_CAE/NSTR_CAD ! Yes this is an integer!
*
      CALL SET_MARK(ISET,NSTR_CAD,NSTR_CAE_REMAIN)
*
      ICOUNT = 0
*
* different I
*
      I = 0 
      DO ILOOP = 1,NSTR_CAD-1
        I = I + ILOOP 
      END DO
!     I = NSTR_CAE*NSTR_CAD
      I = I*NSTR_CAE_REMAIN
      J = I*NSTR_CBE*NSTR_CBD
      K = J*NSTR_AAE*NSTR_AAD
      L = J*NSTR_AAD
      M = K*NSTR_ABD
C     print*,'I,J,K,L,M',I,J,K,L,M
*
      DO ISTR_ABE = 1,NSTR_ABE
C      print*,'ISTR_ABE',ISTR_ABE
        IOFF_ABE = (ISTR_ABE-1)*M
       DO ISTR_AAE = 1,NSTR_AAE
        IOFF_AAE = (ISTR_AAE-1)*L
C       print*,'ISTR_AAE',ISTR_AAE
        DO ISTR_CBE = 1,NSTR_CBE
C        print*,'ISTR_CBE',ISTR_CBE
* Splitting loop over CA for EX
C        DO ISTR_CAE = 1,NSTR_CAE
         DO ISTR_CAE = 1,NSTR_CAD ! Not a typo
C        print*,'ISTR_CAE',ISTR_CAE
         DO ISTR_CAR = 1,NSTR_CAE_REMAIN
* DX start
          DO ISTR_ABD = 1,NSTR_ABD
           IOFF_AB = ISTR_ABD 
           IOFF_AB = (IOFF_AB-1)*K
C          print*,'IOFF_AB',IOFF_AB
           DO ISTR_AAD = 1,NSTR_AAD
            IOFF_AA = ISTR_AAD
            IOFF_AA = (IOFF_AA-1)*J
C           print*,'IOFF_AA',IOFF_AA
            DO ISTR_CBD = 1,NSTR_CBD
             IOFF_CB = (NSTR_CBE * (ISTR_CBD - 1) + ISTR_CBE)
             IOFF_CB = (IOFF_CB-1)*I 
C            print*,'IOFF_CB',IOFF_CB
             DO ISTR_CAD = 1,NSTR_CAD
C             IOFF_CA = (NSTR_CAE * (ISTR_CAD- 1)+ISTR_CAE)
              IF(ISET(ISTR_CAE,ISTR_CAD).NE.0) THEN
              IOFF_CA = ISET(ISTR_CAE,ISTR_CAD) + ISTR_CAR - 1 
!             IOFF_CA = NSTR_CAR * (ISTR_CAE- 1) + IOFF_CA
              IOFF_CA=IOFF_CA+IOFF_CB+IOFF_AA+IOFF_AB+IOFF_AAE+IOFF_ABE
              ICOUNT = ICOUNT + 1
               OUTT(ICOUNT) = T(IOFF_CA)
C              print*,'add'
              ELSE 
               ICOUNT = ICOUNT + 1
! Next line only for show
               IOFF_CA = 0
               OUTT(ICOUNT) = 0.0D0
              END IF
C             print*,'IOFF_CA,ICOUNT',IOFF_CA,ICOUNT
C             print*,'---------------------------'
             END DO
            END DO
           END DO
          END DO
         END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SET_MARK(ISET,NSTR_D,NSTR_E_REMAIN)
*
* To be used with sorting routines
*
      INTEGER ISET(NSTR_D,NSTR_D)
* Local scratch
      INTEGER IMARK(NSTR_D*(NSTR_D-1)/2),ISCR(NSTR_D)
*
      NTEST = 00
*
      NTOT = NSTR_D*(NSTR_D-1)/2
*
* Fill array with offsets
*
      IMARK(1) = 1
      DO ISTR = 2,NTOT
        IMARK(ISTR) = (ISTR-1)* NSTR_E_REMAIN + IMARK(ISTR-1)
      END DO 
*
* Fill the scratch array
*
      ISCR(1) = 0
      DO ISTR =2,NSTR_D
        ISCR(ISTR) = ISCR(ISTR-1) + NSTR_D - ISTR
      END DO
*
      DO ISTR = 1,NSTR_D
        DO JSTR = 1,NSTR_D
          KSTR = ABS(ISTR-JSTR)
          IF(KSTR.GE.1) THEN
            ISET(ISTR,JSTR) = IMARK(KSTR + ISCR(MIN(ISTR,JSTR)))
          ELSE
            ISET(ISTR,JSTR) = 0
          END IF
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' The offsets are '
        CALL IWRTMA(ISET,NSTR_D,NSTR_D,NSTR_D,NSTR_D)
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_6(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NSTR_TOT,T,OUTT)
*
* This is for a free CA index in EX that belongs to the same GAS as the
* free index in DX
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
* Local scratch
      INTEGER ISET(NSTR_CAD,NSTR_CAD)
*
      NTEST = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine number 6 '
      END IF
*
* Calculate the length of remaining indices of CA
*
      NSTR_CAE_REMAIN = NSTR_CAE/NSTR_CAD ! Yes this is an integer!
*
      CALL SET_MARK(ISET,NSTR_CAD,NSTR_CAE_REMAIN)
*
      ICOUNT = 0
*
* different I
*
      I = 0 
      DO ILOOP = 1,NSTR_CAD-1
        I = I + ILOOP 
      END DO
!     I = NSTR_CAE*NSTR_CAD
      I = I*NSTR_CAE_REMAIN
      J = I*NSTR_CBE*NSTR_CBD
      K = J*NSTR_AAE*NSTR_AAD
      L = J*NSTR_AAD
      M = K*NSTR_ABD
C     print*,'I,J,K,L,M',I,J,K,L,M
*
      DO ISTR_ABE = 1,NSTR_ABE
C      print*,'ISTR_ABE',ISTR_ABE
        IOFF_ABE = (ISTR_ABE-1)*M
       DO ISTR_AAE = 1,NSTR_AAE
        IOFF_AAE = (ISTR_AAE-1)*L
C       print*,'ISTR_AAE',ISTR_AAE
        DO ISTR_CBE = 1,NSTR_CBE
         IOFF_CBE = (ISTR_CBE-1)*I 
C        print*,'ISTR_CBE',ISTR_CBE
* Splitting loop over CA for EX
C        DO ISTR_CAE = 1,NSTR_CAE
         DO ISTR_CAE = 1,NSTR_CAD ! Not a typo
C        print*,'ISTR_CAE',ISTR_CAE
          DO ISTR_CAR = 1,NSTR_CAE_REMAIN
* DX start ! Since DX will only be one creation operator everyting else
* is eliminated
           DO ISTR_CAD = 1,NSTR_CAD
C           IOFF_CA = (NSTR_CAE * (ISTR_CAD- 1)+ISTR_CAE)
            IF(ISET(ISTR_CAE,ISTR_CAD).NE.0) THEN
             IOFF_CA = ISET(ISTR_CAE,ISTR_CAD) + ISTR_CAR - 1 
!            IOFF_CA = NSTR_CAR * (ISTR_CAE- 1) + IOFF_CA
             IOFF_CA=IOFF_CA+IOFF_CBE+IOFF_AAE+IOFF_ABE
             ICOUNT = ICOUNT + 1
             IF(ISTR_CAE.GT.ISTR_CAD) THEN
               SIGN = 1.0D0
             ELSE
               SIGN = -1.0D0
             END IF
             OUTT(ICOUNT) = SIGN*T(IOFF_CA)
            ELSE 
             ICOUNT = ICOUNT + 1
! Next line only for show
             IOFF_CA = 0
             OUTT(ICOUNT) = 0.0D0
            END IF
C           print*,'IOFF_CA,ICOUNT',IOFF_CA,ICOUNT
C           print*,'---------------------------'
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
*
      SUBROUTINE SORT_7(NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE,
     &                  NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD,
     &                  NSTR_TOT,T,OUTT)
*
* This is for a free CB index in EX that belongs to the same GAS as the
* free index in DX
*
#include "implicit.inc"
*
      DIMENSION T(*),OUTT(*)
* Local scratch
      INTEGER ISET(NSTR_CBD,NSTR_CBD)
*
      NTEST = 0
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' Welcome to sorting routine number 7 '
      END IF
*
* Calculate the length of remaining indices of CB
*
      NSTR_CBE_REMAIN = NSTR_CBE/NSTR_CBD ! Yes this is an integer!
*
      CALL SET_MARK(ISET,NSTR_CBD,NSTR_CBE_REMAIN)
*
      ICOUNT = 0
*
* different I
*
      J = 0 
      DO JLOOP = 1,NSTR_CBD-1
        J = J + JLOOP 
      END DO
!     I = NSTR_CAE*NSTR_CAD
      I = NSTR_CAE !I*NSTR_CBE_REMAIN
!     J = I*NSTR_CBE*NSTR_CBD
      J = J*I*NSTR_CBE_REMAIN!*NSTR_CBD
      K = J*NSTR_AAE*NSTR_AAD
      L = J*NSTR_AAD
      M = K*NSTR_ABD
C     print*,'I,J,K,L,M',I,J,K,L,M
*
      DO ISTR_ABE = 1,NSTR_ABE
C      print*,'ISTR_ABE',ISTR_ABE
        IOFF_ABE = (ISTR_ABE-1)*M
       DO ISTR_AAE = 1,NSTR_AAE
        IOFF_AAE = (ISTR_AAE-1)*L
C       print*,'ISTR_AAE',ISTR_AAE
* Splitting loop over CB for EX
        DO ISTR_CBE = 1,NSTR_CBD ! Not a typo
!        IOFF_CBE = (ISTR_CBE-1)*I 
C        print*,'ISTR_CBE',ISTR_CBE
         DO ISTR_CBR = 1,NSTR_CBE_REMAIN
C        print*,'ISTR_CBR',ISTR_CBR
          DO ISTR_CAE = 1,NSTR_CAE
           IOFF_CAE = ISTR_CAE-1 
C         print*,'ISTR_CAE',ISTR_CAE
* DX start ! Since DX will only be one creation operator everyting else
* is eliminated
           DO ISTR_CBD = 1,NSTR_CBD
C           IOFF_CA = (NSTR_CAE * (ISTR_CAD- 1)+ISTR_CAE)
            IF(ISET(ISTR_CBE,ISTR_CBD).NE.0) THEN
             IOFF_CB = ISET(ISTR_CBE,ISTR_CBD) + ISTR_CBR - 1 
!            IOFF_CA = NSTR_CAR * (ISTR_CAE- 1) + IOFF_CA
             IOFF_CB=IOFF_CB+IOFF_CAE+IOFF_AAE+IOFF_ABE
             ICOUNT = ICOUNT + 1
             IF(ISTR_CBE.GT.ISTR_CBD) THEN
               SIGN = 1.0D0
             ELSE
               SIGN = -1.0D0
             END IF
             OUTT(ICOUNT) = SIGN*T(IOFF_CB)
            ELSE 
             ICOUNT = ICOUNT + 1
! Next line only for show
             IOFF_CB = 0
             OUTT(ICOUNT) = 0.0D0
            END IF
C           print*,'IOFF_CB,ICOUNT',IOFF_CB,ICOUNT
C           print*,'---------------------------'
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original order of M '
        DO I =1,NSTR_TOT
          WRITE(6,*) ' M ',T(I)
        END DO
        WRITE(6,*) ' EX,DX order '
        DO I =1,ICOUNT
          WRITE(6,*) ' EX,DX ',OUTT(I)
        END DO
      END IF
*
      RETURN
      END
