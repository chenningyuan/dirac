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
      NTEST = 100
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
      NTEST = 100
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
            IF(INUMCON(J,I-1).NE.INUMCON(J,I)) THEN
              ISIM(I)  = 1
              EXIT
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
        CALL IWRTMA(INUMCON,4,4,4,4)
        CALL IWRTMA(ITEMP,4,4,4,4)
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
      NTEST = 100
*
      CALL FROM_O_TO_T_INDX(NPART,NHOLE,NGAS,INUMCON,ITINDX)
*
      NCONNECT = 0
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
* connection found. Is it also a part of the resulting operator
        IF(ICONNECT.EQ.1) THEN
          CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
          WRITE(6,*) ' will check'
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
            CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
          ELSE
            CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
          END IF
        ELSE
          CALL WRT_SPOX_TP_CC_KRCC(TOP(1,ITOP),1)
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
      NTEST = 100
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
      NTEST = 100
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
      NTEST = 100
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
        WRITE(6,*) ' Searched from ',IFROM,' TO ',ITO
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
      IF(N.GT.IEXCITE) THEN
        ICONNECT = 0
        RETURN
      END IF
*
      IF(ABS(MK).GT.IKRFLIP) THEN
        ICONNECT = 0
        RETURN
      END IF
*
      IF(ABS(MUB).GT.IMUB) THEN
        ICONNECT = 0
        RETURN
      END IF
*
* Now for the start and end
*
      IFROM = RESOPS(N,MK,MUB)
      ITO = RESOPF(N,MK,MUB)
*
      ILENGTH = ITO-IFROM
      IF(ILENGTH.LT.0) THEN
        ICONNECT = 0
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
      SUBROUTINE OP_COMPARE(IOP,NOP,IMOP,NMOP,
     &                      IVECS,IVECF,IEXCITE,IKRFLIP,IMTOH,HOP)
* 
* Compares the Hamiltonian types to the intermediates
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "crun.inc"
#include "cgas.inc"
*
      INTEGER IOP(NOP),IMOP(4*NGAS,NMOP)
      INTEGER IVECS(0:IEXCITE,-IKRFLIP:IKRFLIP,-IEXCITE:IEXCITE)
      INTEGER IVECF(0:IEXCITE,-IKRFLIP:IKRFLIP,-IEXCITE:IEXCITE)
      INTEGER IMTOH(NOP) !Misleading name since it is actually H TO M
      INTEGER HOP(4*NGAS,*)
* Scratch
      INTEGER INDX(4)
*
      NTEST = 00
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input Hamiltonian types '
        DO IHOP =1,NOP
          CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(IHOP)),1)
        END DO
        WRITE(6,*) ' Input Cluster operator types '
        CALL WRT_SPOX_TP_CC_KRCC(IMOP,NMOP)
      END IF
*
      CALL ISETVC(IMTOH,0,NOP) 
*
      DO IHOP =1,NOP
* Analyze H operator
        CALL ANALYZE_OP(HOP(1,IOP(IHOP)),NGAS,INH,IMK,IMUB)
* Find non-zero indices
        CALL NONZERO_INDX(HOP(1,IOP(IHOP)),INH,INDX,IDIFFINDX,NGAS)
* Compare indices to intermediate
        ISTART = IVECS(INH,IMK,IMUB)
        IEND = IVECF(INH,IMK,IMUB)
        DO I = ISTART,IEND
          ISIM = 0
          DO J = 1,IDIFFINDX
            IF(HOP(INDX(J),IOP(IHOP)).NE.IMOP(INDX(J),I)) THEN
* No match
              print*,' no match '
              ISIM = 1
              EXIT
            END IF
          END DO
          IF(ISIM.EQ.0) THEN
* We have a match. Store the connection
            print*,' we have a match'
            IMTOH(IHOP) = I
            EXIT
          END IF
        END DO 
      END DO
*
      IF(NTEST.EQ.100) THEN
        WRITE(6,*) ' T operator match for H operator '
        DO I = 1,NOP
          WRITE(6,*) ' H op ',IOP(I),' match T op ', IMTOH (I)
          IF(IMTOH(I).NE.0) THEN
            WRITE(6,*) ' H op match '
            CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(I)),1)
            WRITE(6,*) ' T op match '
            CALL WRT_SPOX_TP_CC_KRCC(IMOP(1,IMTOH(I)),1)
          ELSE
            WRITE(6,*) ' No match for H op '
            CALL WRT_SPOX_TP_CC_KRCC(HOP(1,IOP(I)),1)
          END IF
        END DO
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
      SUBROUTINE IDIM_INT(WORK,KFREE,LFREE,IMULTFAC)
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
      NTEST = 00
      IDUM = 0
*
      MX_ST_TSOSO_INI = MX_ST_TSOSO
      MX_ST_TSOSO_BLK_INI = MX_ST_TSOSO_BLK
      MX_TBLK_INI = MX_TBLK
      MX_SBSTR_INI = MX_SBSTR
*
*  Have to check these does not exceed those for T by simply comparing against them 
*  MX_ST_TSOSO,MX_ST_TSOSO_BLK,MX_TBLK. Largest will be choosen!
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
      IF(NTEST.GE.100) THEN
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
      SUBROUTINE M_TO_EX_DX_FORMAT(IWAY,IDIMEXDX,IOUT,NRANK,ISM,
     &                             MTYPE,TTYPE,RESTYPE,T,
     &                             WORK,KFREE,LFREE)
*
* Subroutine to map a given CAAB-type of cluster operator to a given
* EX,DX set of operators. Notice this is an injective mapping...hence we
* need to first find the dimension of the new array
*
* Strategy:
* 1. Find EX,DX
* 2. Find dimension of EX,DX
* 3. Perform mapping
* Lasse 2011
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
      INTEGER MTYPE(4*NGAS),TTYPE(4*NGAS),RESTYPE(4*NGAS)
      DIMENSION IOUT(IDIMEXDX),T(*),WORK(*)
* Local scratch
      INTEGER IEX(4*NGAS),IDX(4*NGAS)
      INTEGER ICONN(2*NRANK),IHINDEX(40),IEXD1234(2*NRANK),IORD(4)
*
      NTEST = 100
*
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.0) THEN
          WRITE(6,*) ' Will dimension the needed  array '
        ELSE
          WRITE(6,*) ' Will perform the mapping '
        END IF
        WRITE(6,*) ' Intermediate '
        CALL WRT_SPOX_TP_CC_KRCC(MTYPE,1)
        WRITE(6,*) ' Cluster operator '
        CALL WRT_SPOX_TP_CC_KRCC(TTYPE,1)
        WRITE(6,*) ' Resulting operator '
        CALL WRT_SPOX_TP_CC_KRCC(RESTYPE,1)
      END IF
*
* First find EX,DX by calling old routines
*
      CALL CONTR_POS_REL_1T(MTYPE,TTYPE,RESTYPE,NGAS,ICONN,NRANK)
      NTOP = 1
* Fill IHINDEX
      DO I=1,40
        IHINDEX(I) = I
      END DO
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
* Second find the dimension of the new array if IWAY = 0
* else perform the mapping
*
      CALL IDIM_EX_DX(IEX,IDX,MTYPE,NRANK,ISM,IOUT,IDIMEXDX,T,
     &            WORK,KFREE,LFREE)
*
      RETURN
      END
*
      SUBROUTINE IDIM_EX_DX(IEX,IDX,MTYPE,NRANK,ISM,IOUT,IDIMEXDX,T,
     &            WORK,KFREE,LFREE)
*
* Do as the routine name says. Dimension EX and DX.
* Based on ANA_GENCCS_CC_KRCC just with more loops.
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
      DIMENSION IOUT(IDIMEXDX),T(*),WORK(*)
* Local scratch
      INTEGER IOCC_CAE(MX_ST_TSOSO_BLK_MX),IOCC_CBE(MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AAE(MX_ST_TSOSO_BLK_MX),IOCC_ABE(MX_ST_TSOSO_BLK_MX)
*
      INTEGER IOCC_CAD(MX_ST_TSOSO_BLK_MX),IOCC_CBD(MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AAD(MX_ST_TSOSO_BLK_MX),IOCC_ABD(MX_ST_TSOSO_BLK_MX)
*
      INTEGER IOCC_CA(MX_ST_TSOSO_BLK_MX),IOCC_CB(MX_ST_TSOSO_BLK_MX)
      INTEGER IOCC_AA(MX_ST_TSOSO_BLK_MX),IOCC_AB(MX_ST_TSOSO_BLK_MX)
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
      NTEST = 100
      NSAVE = 0
      NSAVE2 = 0
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Rank of EX+DX ',NRANK
        WRITE(6,*) ' Excitation operator '
        CALL WRT_SPOX_TP_CC_KRCC(IEX,1)
        WRITE(6,*) ' De-excitation operator '
        CALL WRT_SPOX_TP_CC_KRCC(IDX,1)
        WRITE(6,*) ' Original operator '
        CALL WRT_SPOX_TP_CC_KRCC(MTYPE,1)
      END IF
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
C        print*,'NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE',
C    &           NSTR_ABE,NSTR_AAE,NSTR_CBE,NSTR_CAE
C        print*,'NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD',
C    &           NSTR_ABD,NSTR_AAD,NSTR_CBD,NSTR_CAD
C        print*,'NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA',
C    &           NSTR_AB,NSTR_AA,NSTR_CB,NSTR_CA
         NSTR_TOTE = NSTR_ABE*NSTR_AAE*NSTR_CBE*NSTR_CAE
         NSTR_TOTD = NSTR_ABD*NSTR_AAD*NSTR_CBD*NSTR_CAD
         NSTR_TOT  = NSTR_AB *NSTR_AA *NSTR_CB *NSTR_CA 
         NTOT =NSTR_TOTE*NSTR_TOTD*NSTR_TOT
         NTOT2 = NSTR_TOTE*NSTR_TOTD
         NSAVE = MAX(NTOT,NSAVE)
         NSAVE2 = MAX(NSTR_TOT,NSAVE2)
C        IF(NTOT.GE.1.OR.NSTR_TOTE.GE.1.OR.NSTR_TOTD.GE.1) THEN
         IF(NTOT.GE.1) THEN
         print*,'NSTR_TOTE,NSTR_TOTD,NSTR_TOT',
     &           NSTR_TOTE,NSTR_TOTD,NSTR_TOT
         END IF
         IF(NSTR_TOT.GE.10877) THEN
         DO I_AB = 1, NSTR_AB
            I_AA_MIN = 1
          DO I_AA = I_AA_MIN, NSTR_AA
           DO I_CB = 1, NSTR_CB
              I_CA_MIN = 1
            DO I_CA = I_CA_MIN, NSTR_CA
             IT = IT + 1
             WRITE(6,*) ' IT, T(IT) = ', IT,T(IT)
*
             IF(ILOOP .EQ. 1 ) THEN
               NTVAR = NTVAR + 1
               TNORM = TNORM + T(IT)**2
*. Classify element according to size
C              DO IPOT = 1, 10
C                IF(10.0D0 ** (-IPOT+1).GE.ABS(T(IT)).AND.
C    &           ABS(T(IT)).GT. 10.0D0 ** ( - IPOT )      )THEN
C                  IOEX = ISOX_TO_OX(ITSS)
C                  NCPMT_OX(IPOT,IOEX)= NCPMT_OX(IPOT,IOEX)+1
C                  WCPMT_OX(IPOT,IOEX)= WCPMT_OX(IPOT,IOEX)+T(IT)**2
C                END IF
C              END DO
*            ^ End of loop over powers of ten
             END IF
*            ^ end if we are in loop 1
C            IF( XMAX .GE. ABS(T(IT)) .AND.
C    &       ABS(T(IT)).GT. XMIN .AND. ITRM.LE.MAXTRM) THEN
               ITRM = ITRM + 1
               IIT = IIT + 1
               TTNORM = TTNORM + T(IT) ** 2
               WRITE(6,'(A)')
               WRITE(6,'(A)')
     &         '                 =================== '
               WRITE(6,*)

      print*,'ISM_CA,ISM_CB,ISM_AA,ISM_AB',ISM_CA,ISM_CB,ISM_AA,ISM_AB
               WRITE(6,'(A,2I8,2X,E14.8)')
     &         '  Type, number, size of amplitude : ',
     &         ITSS, IT, T(IT)
C              IF(IUSLAB.EQ.0) THEN
                 WRITE(6,*) ' Creator alpha '
                 WRITE(6,'(4X,10I4)')
     &           (IOCC_CA(IEL+(I_CA-1)*NEL_CA),IEL = 1, NEL_CA)
                 WRITE(6,*) ' Creator beta '
                 WRITE(6,'(4X,10I4)')
     &           (IOCC_CB(IEL+(I_CB-1)*NEL_CB),IEL = 1, NEL_CB)
                 WRITE(6,*) ' Annihilator alpha '
                 WRITE(6,'(4X,10I4)')
     &           (IOCC_AA(IEL+(I_AA-1)*NEL_AA),IEL = 1, NEL_AA)
                 WRITE(6,*) ' Annihilator beta '
                 WRITE(6,'(4X,10I4)')
     &           (IOCC_AB(IEL+(I_AB-1)*NEL_AB),IEL = 1, NEL_AB)
C              END IF
C            END IF
*            ^ End if could and should be printed
            END DO
*           ^ End of loop over alpha creation strings
           END DO
*          ^ End of loop over beta creation strings
          END DO
*         ^ End of loop over alpha annihilation
         END DO
*        ^ End of loop over beta annihilation
         END IF
            END DO
           END DO
          END DO
         END DO
        END DO
       END DO
      END DO
*     ^ End of loop over symmetry blocks
      IF(NSAVE.EQ.0) print*,'NSAVE is equal to zero',NRANK
*
      RETURN
      END
