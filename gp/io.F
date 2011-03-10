      SUBROUTINE PRQMAT(QMAT,NROW,NCOL,LRQ,LCQ,NZ,IQP,IUNIT)
#include "implicit.h"
      DIMENSION QMAT(LRQ,LCQ,NZ),IQP(*)
      CHARACTER M(4)*6
      M(1) = 'A real'
      M(2) = 'A imag'
      M(3) = 'B real'
      M(4) = 'B imag'
      DO 10 IZ = 1,NZ
        IQ = IQP(IZ)
        WRITE(IUNIT,'(/A,A6,A,4I4/)')
     &  '*** ',M(IQ),' part ***; NROW,NCOL,LRQ,LCQ:',NROW,NCOL,LRQ,LCQ
        CALL OUTPUT(QMAT(1,1,IZ),1,NROW,1,NCOL,LRQ,LCQ,-1,IUNIT)
   10 CONTINUE
      END
!C  /* Deck readi */
      SUBROUTINE READI (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
      DIMENSION INTX(N)
      READ (IT) INTX
      RETURN
      END
C  /* Deck readi2 */
      SUBROUTINE READI2 (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
#if defined (VAR_STAR2)
      INTEGER*2 INTX(N)
#else
      INTEGER   INTX(N)
#endif
      READ (IT) INTX
      RETURN
      END
C  /* Deck readdi */
      SUBROUTINE READDI(IT,IU,N,IX)
      DIMENSION IX(N)
#if defined (VAR_FASTIO)
      ERROR - READDI is not programmed for VAR_FASTIO
#else
      READ(IT,REC=IU) IX
#endif
      RETURN
      END
C  /* Deck readt */
      SUBROUTINE READT (IT,N,X)
#include "implicit.h"
      CHARACTER*30 STRING
      DIMENSION X(N)
      READ (IT,END=10,ERR=20) X
      RETURN
 10   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: END reading unit  ',IT
      CALL QUIT(STRING)
 20   CONTINUE
      WRITE (STRING,'(A25,I5)') 'READT: Error reading unit',IT
      CALL QUIT(STRING)
      END
!C  /* Deck writi */
      SUBROUTINE WRITI (IT,N,INTX)
!C
!C (30-Jan-1984) hjaaj
!C
      DIMENSION INTX(N)
      WRITE (IT) INTX
      RETURN
      END
C  /* Deck writi2 */
      SUBROUTINE WRITI2 (IT,N,INTX)
C
C (30-Jan-1984) hjaaj
C
#if defined (VAR_STAR2)
      INTEGER*2 INTX(N)
#else
      INTEGER   INTX(N)
#endif
      WRITE (IT) INTX
      RETURN
      END
C  /* Deck writdi */
      SUBROUTINE WRITDI(IT,IU,N,IX)
      DIMENSION IX(N)
#if !defined (VAR_FASTIO)
      WRITE(IT,REC=IU) IX
#endif
      RETURN
      END
C  /* Deck writt */
      SUBROUTINE WRITT(IT,N,X)
#include "implicit.h"
      DIMENSION X(N)
      WRITE (IT) X
      RETURN
C
      END
