#ifdef UNDEF
/* Comdeck qpack_log */
 1-Mar-1989 Hans Joergen Aa. Jensen
 Package for trace back, and maybe later call trees and timing
 statistics.

 Modules:
 QDUMP(LUPRI)
 QENTER(IDENT)
 QEXIT(IDENT)
 QINIT
 QTRACE(LUPRI)

 Common blocks:
 QSTACK
 QSTACC (for character variables)

 Ideas for new modules:
 QTREE(LUPRI)
 QSTAT(LUPRI)
 QDUMP(LUPRI) call QTRACE and dump QTREE

 This could maybe be implemented with:
 REAL QTIME
 COMMON /QSTACK/ QTIME(3,MXTREE),
                 LTIME(MXQLVL), LVLQTR(MXTREE),
                 LEVELQ,        LQTREE
 COMMON /QSTACC/ QNAMES(MXQLVL),QTREE(MXTREE)
 in QENTER:
   LEVELQ = LEVELQ + 1
   LQTREE = LQTREE + 1
   IF (LQTREE .LE. MXTREE) THEN
      LVLQTR(LQTREE) = LEVELQ
      QTREE (LQTREE) = IDENT
      QTIME (1,LQTREE) = system time
      QTIME (2,LQTREE) = cpu time
      QTIME (3,LQTREE) = wall time
   END IF
   LTIME(LEVELQ)  = LQTREE
 in QEXIT:
   JTREE  = LVLQTR(LEVELQ)
   IF (JTREE .LE. MXTREE)
  *   QTIME(i,JTREE) = i_time - QTIME(i,JTREE)
   LEVELQ = LEVELQ - 1
 in QINIT:
   LEVELQ = 0
   LQTREE = 0
 with this information, child time can be subtracted from module time:
      TIMCHL = 0.0D0
      DO 100 I = JTREE+1,MIN(LQTREE,MXTREE)
         IF (LVLQTR(I) .LE. LVLQTR(JTREE)) GO TO 200
         TIMCHL = TIMCHL + QTIME(I)
  100 CONTINUE
  200 CONTINUE
      WRITE (LUPRI,'()')
     *   QTREE(JTREE),LVLQTR(JTREE),
     *   QTIME(JTREE),TIMCHL,QTIME(JTREE)-TIMCHL
 a summary can also be made, with accumulated values for each called
 module.

#endif
C  /* Deck qpackage */
      SUBROUTINE QENTER(IDENT)
C
C  1-Mar-1989 Hans Joergen Aa. Jensen
C
C  Package for tracebacks and statistics.
C
      CHARACTER*(*) IDENT
C
#include <priunit.h>
C
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*8   QNAMES(MXQLVL), IDENT8, BLANK8
      PARAMETER (BLANK8 = '        ')
      COMMON /QSTACC/ QNAMES
C
      LOGICAL FIRST
      SAVE    FIRST
      DATA    FIRST /.TRUE./
C
      IF (FIRST) THEN
         FIRST = .FALSE.
         LEVELQ = 0
         LQTREE = 0
      END IF
      LEVELQ = LEVELQ + 1
      LQTREE = LQTREE + 1
      IF (LEVELQ .LE. MXQLVL) THEN
         IDENT8 = IDENT//BLANK8
         QNAMES(LEVELQ) = IDENT8
      END IF
      RETURN
C
      ENTRY QEXIT(IDENT)
      IDENT8 = IDENT//BLANK8
      IF (QNAMES(LEVELQ) .EQ. IDENT8) THEN
         LEVELQ = LEVELQ - 1
      ELSE
         WRITE (LUPRI,'(//A/A,I6,2(/2A),//A)')
     *      ' QEXIT error, exit ID does not match entry ID',
     *      ' ----- stack pt =',LEVELQ,
     *      ' ----- entry ID =',QNAMES(LEVELQ),
     *      ' ----- exit  ID =',IDENT8,
     *      ' Stack dump follows:'
         CALL QDUMP(LUPRI)
         CALL QUIT('QEXIT error, exit ID does not match entry ID')
      END IF
      RETURN
C
      ENTRY QINIT
      LEVELQ = 0
      LQTREE = 0
      FIRST  = .FALSE.
      RETURN
C
      END
C  /* Deck qdump */
      SUBROUTINE QDUMP(LUPRI)
C
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*8   QNAMES(MXQLVL)
      COMMON /QSTACC/ QNAMES
C
      IF (LEVELQ .GT. 0) THEN
         CALL QTRACE(LUPRI)
      END IF
      RETURN
      END
C  /* Deck qtrace */
      SUBROUTINE QTRACE(LUPRI)
C
C  1-Mar-1989 Hans Joergen Aa. Jensen
C
C  Printe traceback from Q stack
C
      COMMON /QSTACK/ LEVELQ, LQTREE
      PARAMETER (MXQLVL = 100)
      CHARACTER*8   QNAMES(MXQLVL)
      COMMON /QSTACC/ QNAMES
C
      IF (LEVELQ .EQ. 0) GO TO 9000
      WRITE (LUPRI,'(//A)') ' QTRACE dump of trace stack'
      WRITE (LUPRI,'(3(/A))')
     *   ' ====================',
     *   '      level    module',
     *   ' ===================='
      DO 100 I = LEVELQ,1,-1
         IF (I .LE. MXQLVL) THEN
            WRITE (LUPRI,'(I11,4X,A)') I,QNAMES(I)
         ELSE
            WRITE (LUPRI,'(I11,4X,A)') I,'not recorded'
         END IF
  100 CONTINUE
      WRITE (LUPRI,'(A/)') ' ===================='
C920522-hjaaj: do not call TRACE because it will abort on some computers
C9000 CALL TRACE
 9000 CONTINUE
      RETURN
      END
