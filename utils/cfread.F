      PROGRAM CFREAD
C***********************************************************************
C
C     Get basis information from a coefficient file
C
C     Written by T. Saue Oct 7 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (LUINP = 1)
      CHARACTER TEXT*74,FMT*6,MXFORM*6
      DIMENSION IDIM(3,2)      
      INTEGER   NARG,IARG,JARG
      CHARACTER*80 ARGC
      NARG = IARGC()
      IF (NARG .EQ. 0) THEN
         CALL GETARG(0,ARGC)
         JARG = LNBLNK(ARGC)
         WRITE(6,'(/3A/)') 'Usage: ',ARGC(1:JARG),' file1 file2 ...'
         STOP ' '
      END IF
      DO IARG = 1,NARG
        CALL GETARG(IARG,ARGC)
        JARG = LNBLNK(ARGC)
        WRITE (LUPRI, '(/2A/)') 
     &     ' CFREAD; Information read from file : ',ARGC(1:JARG)
C
        OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
        REWIND LUINP
        READ (LUINP,END=10,ERR=20) TEXT,NSYM,
     &       ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
        WRITE(LUPRI,'(2A)') ' - Heading :',TEXT
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(A,'//FMT//')') ' - Total energy: ',TOTERG
        WRITE(LUPRI,'(A,2I6)') 
     &       ' - Fermion ircops          : ',(I,I=1,NSYM)
        WRITE(LUPRI,'(A,2A )') 
     &       '-----------------------------',('------',I=1,NSYM)
        WRITE(LUPRI,'(A,2I6)') 
     &       ' - Negative-energy orbitals: ',(IDIM(1,I),I=1,NSYM)
        WRITE(LUPRI,'(A,2I6)') 
     &       ' - Positive-energy orbitals: ',(IDIM(2,I),I=1,NSYM)
        WRITE(LUPRI,'(A,2I6)') 
     &       ' - AO-basis functions      : ',(IDIM(3,I),I=1,NSYM)
        GOTO 30
 10     CONTINUE
        CALL QUIT('REACMO: END OF FILE reading TEXT')
 20     CONTINUE
        CALL QUIT('REACMO: ERROR reading TEXT')
 30     CONTINUE
        CLOSE(LUINP)
      ENDDO
      END

