      PROGRAM RSPREAD
C
C Written by Trond Saue July 4 2001
C based on program LABREAD
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (LUINP = 1)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16
      INTEGER   NARG,IARG
      CHARACTER*80 ARGC
C
      NARG = IARGC()
      IF (NARG .EQ. 0) THEN
         CALL GETARG(0,ARGC)
         WRITE(6,'(/A,A,A/)') 'Usage: ',ARGC,' file1 file2 ...'
         STOP ' '
      END IF
      DO IARG = 1,NARG
         CALL GETARG(IARG,ARGC)
         WRITE (LUPRI, '(/A,A/)') 
     &      ' Solution vectors found on file : ',ARGC
C
         OPEN(LUINP,FILE=ARGC,STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND LUINP
         IREC = 0
         IERR = 0
         IVEC = 0
         LBUF = '                '
    1    CONTINUE
         READ (LUINP,END=3,ERR=2) LAB1
         IF(LAB1.EQ.'END_OF_THIS_FILE') GOTO 10
         BACKSPACE LUINP
         READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,
     &          FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &          JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &          ERGRSP,NBAS,NORB
         IVEC = IVEC + 1
         IREC = IREC + 1
         IF(LAB1.NE.LBUF) THEN
           WRITE(LUPRI,'(A,A16,2(A,I3))') 
     &      '** Solution vector : ',LAB1,
     &           ' Irrep: ',JSYMOP1,' Trev: ',JTIMOP1
           LBUF = LAB1
         ENDIF
         WRITE(LUPRI,'(I3,A,A2,A,F9.4,A,E10.2,A,I10)')
     &     IVEC, ' Type : ',TYP1,' Freq.: ',FREQ1,
     &     ' Rnorm: ',RNORM,' Length:', LEN
         READ(LUINP,END=3,ERR=2)
         IREC = IREC + 1
         GO TO 1
C
    2    CONTINUE
         IREC = IREC + 1
         IERR = IERR + 1
         WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
         GOTO 10
C
    3    CONTINUE
         WRITE (LUPRI,'(/I10,A)') IREC,
     *      ' records read before EOF on file.'
         GOTO 10
 10      CONTINUE
         CLOSE(LUINP)
      ENDDO
C
      END
