      SUBROUTINE MEMMAR(KBASE,KADD,TASK,IR,IDENT)
*
* Memory manager routine for LUCIAREL "R"
*
* KBASE : New base address
*         If TASK = INI, KBASE is offset for memory to be controlled
*         by MEMMAR
* KADD  : Dimension of array to be added
*         If TASK = INI, KADD is total length of array
* TASK  : = INI  : Initialize                 Character*6
*         = ADDS : Add static memory
*         = ADDL : Add Local memory
*         = FLUSH : Flush local memory
*         = CHECK : Check memory paddings
*         = FREE  : Return first Free word in KBASE
*         = MARK  : Set a mark at current free adress
*         = FLUSM : Flush local memory to previous mark
*         = PRINT : Print memory allocations
* IR    : 1 => integer , 2 => real,
*         ratio between integer and real is IRAT
* IDENT : identifier of memory slice,Character*6
*
* Local Memory not flushed before allocation of additional static memory
* is tranferred to static memory
*
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER*6 TASK,IDENT,IIDENT,MARKC
      PARAMETER(NPAD = 1 )
      PARAMETER(MAXLVL = 4096 )
      PARAMETER(MAXMRK = 1024)
      COMMON/CMEMOR/NWORD,KFREES,KFREEL,NS,NL,NM,IBASE(MAXLVL),
     &              LENGTH(MAXLVL),IIDENT(MAXLVL),IMARK(MAXMRK),
     &              MARKL(MAXMRK),MARKS(MAXMRK),MARKC(MAXMRK)
*. Two real*8 words, one added NPAD times before each array, another
*. added NPAD times after each array
      DATA PAD1/0.123456789D0/
      DATA PAD2/0.987654321D0/
*
#include "mxpdim.inc"
#include "irat.inc"
#include "wrkspc.inc"
*
      ILPRNT = 0
      IF(ILPRNT .NE. 0 ) THEN
        WRITE(6,*) ' KBASE KADD, IR ', KBASE,KADD,IR
        WRITE(6,'(3A)') ' TASK , IDENT ', TASK,IDENT
        WRITE(6,*) '  NS NL NM ',NS,NL,NM
      END IF
*
      ISTOP = 0
      ITSOK  = 1
      IPRNTMP = 0
*
      IF(TASK(1:5).EQ.'PRINT') THEN
        IPRNTMP = 1
        NSNLI= NS + NL
*
      ELSE IF(TASK(1:3).EQ.'INI') THEN
*
**.Initialize
*
        NS = 0
        NL = 0
        NSNLI = 0
        NM = 0
        KFREES = KBASE
        KFREEL = KBASE
*       NWORD = MXPWRD
*
*       Compute pointer to end element
*
        NWORD = KBASE + MXPWRD-1
        IPRNTMP = 0
        ISTOP = 0
*
** First free word
*
      ELSE IF (TASK(1:4).EQ.'FREE') THEN
       KBASE = KFREEL
       NSNLI = NS+NL
*
**. Static memory
*
      ELSE IF(TASK(1:4).EQ.'ADDS') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREES = KFREEL +(KADD+1)/IRAT + 2*NPAD
        ELSE
          KFREES = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREES-1 .GT. NWORD ) THEN
          WRITE(6,*)
          WRITE(6,*) ' You can''t always get what you want'
          WRITE(6,*) ' No, you can''t always get what you want'
          WRITE(6,*) ' But if you try sometime, you may find '
          WRITE(6,*) ' you get what you need '
          WRITE(6,*) '                       Jagger/Richard '
*
          WRITE(6,*) ' MEMMAR : work array too short '
          WRITE(6,*) ' current and required length ',NWORD,KFREES-1
*
          WRITE(6,*) ' Trying to allocate : identifer,length'
          WRITE(6,'(20X,A,I9)')IDENT,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NS = NS + NL + 1
        NL = 0
        NSNLI = NS + NL
        IF(NS.GT.MAXLVL) THEN
          WRITE(6,*) ' Too many levels in MEMMAR '
          write(6,*) ' NS has value  ',NS
          WRITE(6,*) ' Increase MAXLVL from ', MAXLVL
          Call Abend1( 11 )
        END IF
        IIDENT(NS) = IDENT
        KFREEL = KFREES
        IBASE(NS) = KBASE
        DO 10 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   10   CONTINUE
*
**. Local memory
*
      ELSE IF(TASK(1:4).EQ.'ADDL') THEN
        KBASE = KFREEL+NPAD
        IF( IR .EQ. 1 ) THEN
          KFREEL = KFREEL +(KADD+1)/IRAT + 2*NPAD
        ELSE
          KFREEL = KFREEL + KADD + 2*NPAD
        END IF
        IF ( KFREEL-1 .GT. NWORD ) THEN
          WRITE(6,*)
          WRITE(6,*) ' You can''t always get what you want'
          WRITE(6,*) ' No, you can''t always get what you want'
          WRITE(6,*) ' But if you try sometime, you may find '
          WRITE(6,*) ' you get what you need '
          WRITE(6,*) '                       Jagger/Richard '

          WRITE(6,*) ' MEMMAR : work array too short '
          WRITE(6,*) ' current and required length ',NWORD,KFREEL-1
          WRITE(6,*) ' Trying to allocate : identifer,offset,  length'
          WRITE(6,'(24X,A,2I9)')IDENT,KBASE,KADD
          ISTOP = 1
          IPRNTMP = 1
          NSNLI = NS+NL
          GOTO 1001
        END IF
        NL =  NL + 1
        NSNLI = NS+NL
        IF(NS+NL.GT.MAXLVL) THEN
          WRITE(6,*) ' Too many levels in MEMMAR '
          WRITE(6,*) ' Increase MAXLVL from ', MAXLVL
          Call Abend1( 11 )
        END IF
        IIDENT(NS+NL) = IDENT
        IBASE(NS+NL) = KBASE
        DO 20 IPAD = 1, NPAD
          WORK(KBASE-NPAD-1+IPAD) = PAD1
          WORK(KFREEL-NPAD-1+IPAD) =  PAD2
   20   CONTINUE
*
** Flush local memory
*
      ELSE IF(TASK(1:5).EQ.'FLUSH') THEN
        NSNLI = NS+NL
        KFREEL = KFREES
        NL = 0
*. Flush output unit
        LU6 = 6
        CALL GFLUSH(LU6)
      ELSE IF(TASK(1:4).EQ.'MARK') THEN
*. Set a mark at current free address
        NM = NM + 1
        IF(NM.GT.MAXMRK) THEN
          WRITE(6,*) ' Too many marks  in MEMMAR '
          WRITE(6,*) ' Increase MAXMRK from ', MAXMRK
          Call Abend1( 11 )
        END IF
        MARKC(NM) = IDENT
        IMARK(NM) = KFREEL
        MARKL(NM) = NL
        MARKS(NM) = NS
        NSNLI = NS + NL
      ELSE IF (TASK(1:5).EQ.'FLUSM') THEN
        NSNLI = NS+NL
        KFREEL = IMARK(NM)
        IF(KFREES.GT.IMARK(NM)) KFREES = IMARK(NM)
        NL = MARKL(NM)
        NS = MARKS(NM)
        NM = NM - 1
*. Flush output unit
        LU6 = 6
        CALL GFLUSH(LU6)
      ELSE IF( TASK(1:5).EQ.'CHECK') THEN
        NSNLI = NS+ NL
      ELSE
          WRITE(6,*) ' MEMMAR : Unknown task parameter ',TASK
          WRITE(6,*) ' Too confused to continue  '
          Call Abend1( 11 )
      END IF
*
**. Check paddings
*
      ICHECK = 1
      ITSOK = 1
      IF(TASK(1:5).EQ.'CHECK'.OR.ICHECK.EQ.1) THEN
        DO 100 IL = 1, NSNLI
          JBASE = IBASE(IL)
          IF(IL.NE.NSNLI) THEN
           JBASEN = IBASE(IL+1)
          ELSE
           JBASEN = KFREEL + 1
          END IF
          L1OK = 1
          L2OK = 1
          DO IPAD = 1, NPAD
            IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1 .OR.
     &         WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) THEN
               ITSOK = 0
               WRITE(6,*) ' Memory problem for : '
               WRITE(6,*) '   Level (IL) ',IL
               IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
               IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
               IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
                WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &          IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
               ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
                WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &          IIDENT(IL),IBASE(IL),'    OKAY       -  '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
                WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &          IIDENT(IL),IBASE(IL),'     -       OKAY '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
                  WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &            IIDENT(IL),IBASE(IL),'     -       OKAY '
               ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
                WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &          IIDENT(IL),IBASE(IL),'     -       -    '
               END IF
            END IF
          END DO
  100   CONTINUE
      END IF
 1001 CONTINUE
*
        IF(ITSOK.EQ.0.OR.IPRNTMP.NE.0) THEN
          WRITE(6,'(A,A)') ' Current task : ', TASK
          WRITE(6,*) ' NS, NL, NSNLI',NS,NL,NSNLI
          WRITE(6,'(A,A)') ' Current Identifier : ', IDENT
          IF(ITSOK.EQ.0)
     &    WRITE(6,*) '  Sorry to say it , but memory is CORRUPTED '
          WRITE(6,*) '  Memory map : '
          WRITE(6,*) '  Identifier   Offset    Pad1 okay Pad2 okay '
          WRITE(6,*) '  ========== ==========  ========= ========= '
          DO 200 IL = 1, NSNLI
            JBASE = IBASE(IL)
            IF(IL.NE.NSNLI) THEN
             JBASEN = IBASE(IL+1)
            ELSE
             JBASEN = KFREEL + 1
            END IF
            L1OK = 1
            L2OK = 1
            DO 40 IPAD = 1, NPAD
              IF(WORK(JBASE-NPAD-1+IPAD).NE.PAD1) L1OK = 0
   40       CONTINUE
            DO 50 IPAD = 1, NPAD
              IF(WORK(JBASEN-2*NPAD-1+IPAD).NE.PAD2) L2OK = 0
   50       CONTINUE
            IF(L1OK.EQ.1.AND.L2OK.EQ.1) THEN
               WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY     OKAY '
            ELSE IF(L1OK.EQ.1.AND.L2OK.EQ.0) THEN
               WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &         IIDENT(IL),IBASE(IL),'    OKAY       -  '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.1) THEN
               WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       OKAY '
            ELSE IF(L1OK.EQ.0.AND.L2OK.EQ.0) THEN
               WRITE(6,'(1H ,4X,A,I10,4X,A)')
     &         IIDENT(IL),IBASE(IL),'     -       -    '
            END IF
  200     CONTINUE
*
* Marks
*
          WRITE(6,*)
          WRITE(6,*) '======='
          WRITE(6,*) ' Marks '
          WRITE(6,*) '======='
          WRITE(6,*)
*
          WRITE(6,*) ' Identifier  Start of free memory '
          WRITE(6,*) ' ================================='
          DO JMARK = 1, NM
            WRITE(6,'(3X,A6,10X,I10)') MARKC(JMARK),IMARK(JMARK)
          END DO
*
        IF(ITSOK.EQ.0) 
     &     Call Abend2( ' Error observed by  memory manager ' )
        END IF
*
      IF(ISTOP.NE.0) call Abend2
     &               ( ' Error observed by  memory manager ' )
      RETURN
      END
