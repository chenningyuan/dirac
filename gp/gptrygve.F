C
C...   Copyright (c) 1997 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dalton, an ab initio electronic structure program, Release 1.0
C...   (1997), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, T. Andersen, K. L. Bak, V. Bakken,
C...   O. Christiansen, P. Dahle, E. K. Dalskov, T. Enevoldsen,
C...   H. Heiberg, D. Jonsson, S. Kirpekar, R. Kobayashi, H. Koch,
C...   K. V. Mikkelsen, P. Norman, M. J. Packer, T.Saue,
C...   P. R. Taylor, and O. Vahtras"
C...
C...   This source code is provided under a written licence and may be 
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may 
C...   be distributed outside the research group of the licence holder. 
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence. 
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html 
C
C
C$Id$
C
C  /* Deck header */
      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
#include <priunit.h>
C
      LENGTH = LEN(HEAD)
      IF (IN .GE. 0) THEN
         INDENT = IN + 1
      ELSE
         INDENT = (72 - LENGTH)/2 + 1
      END IF
      WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
      WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), ('-',I=1,LENGTH)
      WRITE (LUPRI, '()')
      RETURN
      END
C  /* Deck timer */
      SUBROUTINE TIMER(TEXT,TIMSTR,TIMEND)
#include <implicit.h>
      CHARACTER TEXT*6
#include <priunit.h>
C
      IF (TEXT .EQ. 'START ') THEN
         TIMSTR = SECOND()
      ELSE
#if !defined (SYS_T3D)
         TIMEND = SECOND()
#else
         TIMEND = SECOND() + 1.0D0
#endif
         TIME   = TIMEND - TIMSTR
         TIMSTR = TIMEND
         WRITE (LUPRI,'()')
         CALL TIMTXT('>>> Time used in '//TEXT//' is',TIME,LUPRI)
         WRITE (LUPRI,'()')
         CALL FLSHFO(LUPRI)
      END IF
      RETURN
      END
C  /* Deck timdat */
      SUBROUTINE TIMDAT
#include <implicit.h>
#include <priunit.h>
#if defined (SYS_VAX)
      CHARACTER CDATE*9,CTIME*8
      CALL DATE(CDATE)
      CALL TIME(CTIME)
      WRITE(LUPRI,100) CDATE,CTIME
  100 FORMAT(/,'  Date and time (VAX) : ',A9,2X,A8)
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
      INTEGER CDATE, CTIME
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE(LUPRI,100) CDATE, CTIME
  100 FORMAT(/,'  Date and time (CRAY) : ',A8,2X,A8)
#endif
#if defined (SYS_ALLIANT)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Alliant) : ',A)
#endif
#if defined (SYS_CONVEX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Convex) : ',A)
#endif
#if defined (SYS_AIX)
C     AIX XL FORTRAN version 2.3+
      CHARACTER*(24) fdate_
      WRITE (LUPRI,100) fdate_()
  100 FORMAT(/,' Date and time (IBM-AIX): ',A)
#endif
#if defined (SYS_HPUX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (HP-UX)  : ',A)
#endif
#if defined (SYS_DEC)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (DEC)  : ',A)
#endif
#if defined (SYS_IRIX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (IRIX)   : ',A)
#endif
C#if defined (SYS_LINUX)
C      CHARACTER*(24) FTIME
C      WRITE (LUPRI,100) FTIME()
C  100 FORMAT(/,' Date and time (Linux) : ',A)
C#endif
#if defined (SYS_SUN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (SUN)   : ',A)
#endif
#if defined (SYS_PARAGON)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Paragon): ',A)
#endif
#if defined (SYS_IBM)
      CHARACTER DAOFWK*9, LDATE*8, LTIME*8
      CALL DATE(DAOFWK,LDATE)
      CALL TIME(LTIME)
      WRITE(LUPRI,100) DAOFWK, LDATE, LTIME
  100 FORMAT(/,'  Date and time (IBM) : ',A9,2X,A8,2X,A8)
#endif
      RETURN
      END
C  /* Deck timpri */
      SUBROUTINE TIMPRI(TEXT,TIME,TIMALL)
#include <implicit.h>
      CHARACTER TEXT*6
      PARAMETER (HUN = 100.0D00)
#include <priunit.h>
      SAVE TREST
      DATA TREST /0.0D0/
      IF (TEXT .EQ. 'REST  ') TIME = TREST
      ITIME = NINT(HUN*TIME/TIMALL)
      IF (ITIME .GT. 0 .AND. INT(TIME) .GT. 0) THEN
         MINUTE = INT(TIME)/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = NINT(TIME) - 3600*IHOURS - 60*MINUTE
         WRITE(LUPRI,100) TEXT, IHOURS, MINUTE, ISECND, ITIME
      ELSE
         TREST = TREST + TIME
      END IF
  100 FORMAT(1X,A6,'     ',I2.2,':',I2.2,':',I2.2,5X,I3,' %')
      RETURN
      END
C  /* Deck titler */
      SUBROUTINE TITLER(HEAD,A,IN)
      CHARACTER HEAD*(*), A
#include <priunit.h>
C
      LENGTH = LEN(HEAD)
      IF (IN .EQ. 200) THEN
         LENGTH = LENGTH + 2
      ELSE IF (IN .GE. 100) THEN
         MARG = IN - 100
         IF (MARG .GT. 0) MARG = MARG + 1
         LENGTH = LENGTH + 2*MARG
      END IF
      IF (IN .GE. 0 .AND. IN .LT. 100) THEN
         INDENT = IN + 1
      ELSE
         INDENT = MAX(1,(72 - LENGTH)/2 + 1)
      END IF
      IF (IN .EQ. 200) THEN
         WRITE (LUPRI, '(//,80A)')
     *      (' ',I=1,INDENT), '.', ('-',I=1,LENGTH), '.'
         WRITE (LUPRI,'(80A)') (' ',I=1,INDENT),'| ', HEAD, ' |'
         WRITE (LUPRI, '(80A)')
     *      (' ',I=1,INDENT), '`', ('-',I=1,LENGTH), ''''
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD, ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//,80A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(80A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck around */
      SUBROUTINE AROUND(HEAD)
      CHARACTER HEAD*(*)
#include <priunit.h>
      LNG = LEN(HEAD) + 2
      IND = (72 - LNG)/2 + 1
      WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '! ', HEAD, ' !'
      WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '+', ('-',I=1,LNG), '+'
C     WRITE (LUPRI, '(//,80A)') (' ',I=1,IND), '.', ('-',I=1,LNG), '.'
C     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '| ', HEAD, ' |'
C     WRITE (LUPRI, '(80A)')    (' ',I=1,IND), '`', ('-',I=1,LNG), ''''
      WRITE (LUPRI, '()')
      RETURN
      END
C  /* Deck stopit */
      SUBROUTINE STOPIT(SUB,PLACE,INT1,INT2)
#include <implicit.h>
#include <priunit.h>
      CHARACTER*(*) SUB, PLACE
      WRITE (LUPRI,'(//3A)') ' Work space exceeded in subroutine ',
     *                         SUB,'.'
      IF ((LEN(PLACE) .GT. 1) .OR. (PLACE .NE. ' ')) THEN
         WRITE (LUPRI,'(/2A)') ' Location: ',PLACE
      END IF
      JNT1  = ABS(INT1)
      JNT2  = ABS(INT2)
      LWRKR = MAX(JNT1,JNT2)
      LWRKA = MIN(JNT1,JNT2)
      IF (INT1 .LT. 0 .OR. INT2 .LT. 0) THEN
         WRITE (LUPRI,'(/A,I10)  ') ' Space required  >',LWRKR
      ELSE
         WRITE (LUPRI,'(/A,I10)  ') ' Space required  =',LWRKR
      END IF
      WRITE (LUPRI,'( A,I10)') ' Space available =',LWRKA
      CALL QUIT('Work space exceeded.')
      END
C  /* Deck quit */
      SUBROUTINE QUIT(TEXT)
#include <implicit.h>
#if defined (VAR_MPI)
      INCLUDE 'mpif.h'
#endif
#include <priunit.h>
      CHARACTER TEXT*(*)
      WRITE (LUPRI,'(/A)')
     *   '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
      CALL TIMDAT
      WRITE (LUPRI,'(1X,A)') TEXT
#if defined (SYS_CRAY) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_UNIX) || defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_IRIX) || defined (SYS_DEC) || defined (SYS_PARAGON) || defined (SYS_T3D) || defined (SYS_LINUX) || defined (SYS_SUN) || defined (SYS_HPUX)
C     Write to stderr
      WRITE (0,'(/A/1X,A)')
     *   '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---',TEXT
#endif
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI)
      CALL MPI_ABORT(MPI_COMM_WORLD,IERR,IERR)
#endif
#if defined (VAR_PVM)
      CALL PVMFEXIT(INFO)
#endif
      CALL TRACE
      STOP
      END
C  /* Deck opendx */
      SUBROUTINE OPENDX (LUDX,NAME,NELEM,NREC,STATUS,LRDX,NBDX,OLDDX)
C
C 15-Jun-1985 hjaaj
C
C Revisions :  9-Dec-1987 hjaaj (Alliant version)
C
C Purpose:
C   Open files for direct access through WRITDX and READDX routines.
C   The ....DX routines enables direct access, even when the number
C   of elements per record (the logical record length) is greater
C   than the maximum physical record length.
C   >>> THIS IS MACHINE DEPENDENT <<<
C
C Input:
C  LUDX     file unit number
C  NELEM    number of integer words per logical record
C  NREC     number of logical records
C  STATUS   file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output:
C  LRDX     physical record length (in integers)
C  NBDX     number of physical records per logical record
C  OLDDX    logical, true if old LUDX file was opened
C
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLDDX
#include <priunit.h>
#if defined (SYS_VAX)
C
C     VAX/VMS 3.5 OPEN STATEMENTS (APRIL 1985)
C     VAX/VMS CALCULATES RECORD LENGTH, WHEN UNFORMATTED, IN LONGWORDS
C             (4 BYTES) AND SIZE IN DISK BLOCKS (512 BYTES).
C             MAXIMUM RECORD LENGTH (UNFORMATTED) FOR DIRECT ACCESS
C             IS 8191 LONGWORDS, 4095 REAL*8 WORDS, 2047 REAL*16 WORDS.
C             A RECORD LENGTH OF 8188 IS USED BELOW, THIS ENABLES
C             THIS ROUTINE TO BE USED WITH QUAUDRUPLE PRECISION.
C
C             (THIS MEANS MAXIMUM 90 ORBITALS FOR LUSFDA,
C              AND MAXIMUM NELEM = 4095 FOR LUDX, FOR ONE PHYSICAL
C              RECORD PER LOGICAL RECORD.)
C
      NBDX   = (NELEM-1)/8188 + 1
      LRDX   = (NELEM-1)/NBDX + 1
      ILDX   = NREC*NBDX * (LRDX / 128 + 1)
C
      IF (STATUS .EQ. 'NEW') GO TO 300
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRDX)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRDX,INITIALSIZE=ILDX)
         OLDDX = .FALSE.
  600 CONTINUE
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
C
C     CRAY has 8 byte integers.
C
      NBDX   = 1
      LRDX   = NELEM
      LRECL  = 8*LRDX
C
      IF (STATUS .EQ. 'NEW') GO TO 300
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600 CONTINUE
#endif
#if defined (SYS_IBM)
      LOGICAL FEXIST
C
C     IBM version -- 24-Jun-1988 hjaaj (= Alliant version)
C
      NBDX   = 1
      LRDX   = NELEM
      LRECL  = 4*LRDX
C     ILDX   = NREC*NBDX * LRECL for "initialsize" in bytes
C
      IF (STATUS .EQ. 'OLD') GO TO 200
      IF (STATUS .EQ. 'NEW') GO TO 300
      IF (STATUS .EQ. 'UNKNOWN') THEN
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            GO TO 200
         ELSE
            GO TO 300
         END IF
      END IF
C     illegal STATUS keyword
      GO TO 9000
C
C     OPEN OLD FILE
C
  200    OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',
     *        ERR=9100,ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600  CONTINUE
#endif
#if defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_PARAGON) || defined (SYS_LINUX) || defined (SYS_SUN)
C
C     Alliant version -- 9-Dec-1987 hjaaj; recl in bytes
C     Assume same for CONVEX --  1-Nov-1989 hjaaj
C     Assume same for IBM-AIX--  1-Oct-1990 hjaaj
C     Assume same for HP-UX  -- 21-Aug-1991 hjaaj
C     Assume same for DEC    -- 21-May-1992 hjaaj
C     Assume same for IRIX   --  3-Feb-1994 hjaaj; except recl in *4 words
C     Assume same for PARAGON-- 13-Oct-1994 hjaaj
C
      NBDX   = 1
      LRDX   = NELEM
#if defined (SYS_IRIX)
      LRECL  = LRDX
#else
      LRECL  = 4*LRDX
#endif
C     ILDX   = NREC*NBDX * LRECL for "initialsize" in bytes
C
      IF (STATUS .EQ. 'NEW') GO TO 300
      IF (STATUS .NE. 'OLD' .AND. STATUS .NE. 'UNKNOWN') GO TO 9000
C
C     OPEN OLD FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='OLD',FORM='UNFORMATTED',ERR=300,
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .TRUE.
      GO TO 600
C
  300 CONTINUE
      IF (STATUS .EQ. 'OLD') GO TO 9100
C
C     OPEN NEW FILE
C
         OPEN(LUDX,FILE=NAME,STATUS='NEW',FORM='UNFORMATTED',
     *        ACCESS='DIRECT',RECL=LRECL)
         OLDDX = .FALSE.
  600  CONTINUE
#endif
#if !defined (SYS_VAX) && !defined (SYS_CRAY) && !defined (SYS_IBM) && !defined (SYS_ALLIANT) && !defined (SYS_CONVEX) && !defined (SYS_AIX) && !defined (SYS_HPUX) && !defined (SYS_DEC) && !defined (SYS_IRIX) && !defined (SYS_PARAGON) && !defined (SYS_T3D) && !defined (SYS_LINUX) && !defined (SYS_SUN)
 >>>>> insert appropriate OPEN statements in OPENDX.
#endif
      RETURN
C
C error branches
C
 9000 CONTINUE
      WRITE (LUPRI,'(//A,A/A,I5)')
     *   ' *** ERROR (OPENDX) INVALID STATUS KEYWORD: ',STATUS,
     *   '                    FILE NUMBER =',LUDX
      CALL QTRACE(LUPRI)
      CALL QUIT('*** ERROR (OPENDX) INVALID STATUS KEYWORD')
C
 9100 CONTINUE
      WRITE (LUPRI,'(//A/A,I5/A)')
     *   ' *** ERROR (OPENDX) OLD FILE NOT FOUND',
     *   '                    FILE NUMBER =',LUDX,
     *   ' --- or wrong record length on old file.'
      CALL QTRACE(LUPRI)
      CALL QUIT('*** ERROR (OPENDX) FILE NOT FOUND')
C
C end of OPENDX
C
      END
C  /* Deck finddx */
      LOGICAL FUNCTION FINDDX(LU,LRDX,I,LEN,IVEC)
C
C 27-Jun-1985 Hans Jorgen Aa. Jensen
C
C For direct access find record,
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      IF (LEN .LE. LRDX) THEN
         READ (LU, REC=I, IOSTAT=IOS) IVEC
         IF (IOS .NE. 0) GO TO 900
      ELSE
         NBUF = (LEN-1)/LRDX + 1
         IREC = 1 + NBUF*(I-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF-1
            READ (LU, REC=IREC, IOSTAT=IOS) (IVEC(JADD+J), J = 1,LRDX)
            IF (IOS .NE. 0) GO TO 900
            IREC = IREC + 1
            JADD = JADD + LRDX
  100    CONTINUE
         READ (LU, REC=IREC, IOSTAT=IOS) (IVEC(J), J = JADD+1,LEN)
         IF (IOS .NE. 0) GO TO 900
      END IF
      FINDDX = .TRUE.
      RETURN
C
  900 CONTINUE
      FINDDX = .FALSE.
      RETURN
      END
C  /* Deck readdx */
      SUBROUTINE READDX(LU,LRDX,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access read
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      IF (LEN .LE. LRDX) THEN
         READ (LU, REC = I) IVEC
      ELSE
         NBUF = (LEN-1)/LRDX + 1
         IREC = 1 + NBUF*(I-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF-1
            READ (LU, REC = IREC) (IVEC(JADD+J), J = 1,LRDX)
            IREC = IREC + 1
            JADD = JADD + LRDX
  100    CONTINUE
         READ (LU, REC = IREC) (IVEC(J), J = JADD+1,LEN)
      END IF
      RETURN
      END
C  /* Deck writdx */
      SUBROUTINE WRITDX(LU,LRDX,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access write
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      IF (LEN .LE. LRDX) THEN
         WRITE (LU, REC = I) IVEC
      ELSE
         NBUF = (LEN-1)/LRDX + 1
         IREC = 1 + NBUF*(I-1)
         JADD = 0
         DO 100 IBUF = 1,NBUF-1
            WRITE (LU, REC = IREC) (IVEC(JADD+J), J = 1,LRDX)
            IREC = IREC + 1
            JADD = JADD + LRDX
  100    CONTINUE
         WRITE (LU, REC = IREC) (IVEC(J), J = JADD+1,LEN)
      END IF
      RETURN
      END
C  /* Deck aolab4 */
      SUBROUTINE AOLAB4(IINDPK,NMAX,NIBUF,NBITS,IINDX4,N)
C
C     Written by Henrik Koch 22-Nov-1991 (as AOLABE)
C     Generalized 25-Oct-1993 hjaaj
C                 31-Jul-1996 hjaaj (new NIBUF parameter)
C
C     Unpack 4 integer indices packed in NIBUF *4 integer(s).
C
#include <implicit.h>
      DIMENSION IINDPK(*), IINDX4(4,*)
#include <ibtdef.h>
C
      N = IINDPK(NIBUF*NMAX+1)
C
      IF (NIBUF .EQ. 1) THEN
      IF (NBITS .EQ. 8) THEN
         DO 100 I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,24),IBT08)
            IINDX4(2,I) = IBTAND(IBTSHR(LABEL,16),IBT08)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL, 8),IBT08)
            IINDX4(4,I) = IBTAND(       LABEL,    IBT08)
  100    CONTINUE
      ELSE IF (NBITS .EQ. 16) THEN
         DO 200 I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,48),IBT16)
            IINDX4(2,I) = IBTAND(IBTSHR(LABEL,32),IBT16)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(4,I) = IBTAND(       LABEL,    IBT16)
  200    CONTINUE
      ELSE
	 CALL QUIT('AOLAB4 error')
      END IF
      ELSE
         DO 300 I = 1,N
            LABEL = IINDPK(I)
            IINDX4(1,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(2,I) = IBTAND(       LABEL    ,IBT16)
            LABEL = IINDPK(NMAX+I)
            IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
            IINDX4(4,I) = IBTAND(       LABEL    ,IBT16)
  300    CONTINUE
      END IF
C
      RETURN
      END
