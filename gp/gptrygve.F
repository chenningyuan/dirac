C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C
C
C VERSION : $Revision: 1.17.4.2 $
C DATE    : $Date: 2004/09/01 08:10:57 $
C FILE    : gptrygve.F
C
C
C  /* Deck header */
      SUBROUTINE HEADER(HEAD,IN)
      CHARACTER HEAD*(*)
#include "priunit.h"
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
#include "implicit.h"
      CHARACTER TEXT*6
#include "priunit.h"
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
#include "implicit.h"
#include "priunit.h"
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      INTEGER CDATE, CTIME
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE(LUPRI,100) CDATE, CTIME
  100 FORMAT(/,'  Date and time (CRAY) : ',A8,2X,A8)
#endif
#if defined (SYS_CONVEX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Convex) : ',A)
#endif
#if defined (SYS_HAL)
      CHARACTER*(24) CFDATE
      CALL FDATE(CFDATE)
      WRITE (LUPRI,100) CFDATE
  100 FORMAT(/,' Date and time (HAL)    : ',A)
#endif
#if defined (SYS_AIX)
C     AIX XL FORTRAN version 2.3+
      CHARACTER*(24) fdate
      WRITE (LUPRI,100) fdate()
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
#if defined (SYS_LINUX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Linux) : ',A)
#endif
#if defined (SYS_DARWIN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Darwin) : ',A)
#endif
#if defined (SYS_SX)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (Nec-SX) : ',A)
#endif
#if defined (SYS_SUN)
      CHARACTER*(24) FDATE
      WRITE (LUPRI,100) FDATE()
  100 FORMAT(/,' Date and time (SUN)   : ',A)
#endif
      RETURN
      END
C  /* Deck timpri */
      SUBROUTINE TIMPRI(TEXT,TIME,TIMALL)
#include "implicit.h"
      CHARACTER TEXT*6
      PARAMETER (HUN = 100.0D00)
#include "priunit.h"
      SAVE TREST
      DATA TREST /0.0D0/
      IF (TEXT .EQ. 'REST  ') THEN
         TIME = TREST
C 980826-hjaaj: initialize TREST for next abacus call
         TREST = 0.0D0
      END IF
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
#include "priunit.h"
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
         WRITE (LUPRI, '(//,120A)')
     *      (' ',I=1,INDENT), '.', ('-',I=1,LENGTH), '.'
         WRITE (LUPRI,'(120A)') (' ',I=1,INDENT),'| ', HEAD, ' |'
         WRITE (LUPRI, '(120A)')
     *      (' ',I=1,INDENT), '`', ('-',I=1,LENGTH), ''''
      ELSE IF (IN .EQ. 100) THEN
         WRITE (LUPRI, '(//,120A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(120A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(120A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE IF (IN .GT. 100) THEN
         WRITE (LUPRI, '(//,120A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
         WRITE (LUPRI, '(120A)') (' ',I=1,INDENT),
     *      (A,I=1,MARG-1), ' ', HEAD, ' ', (A,I=1,MARG-1)
         WRITE (LUPRI, '(120A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      ELSE
         WRITE (LUPRI, '(//,120A)') (' ',I=1,INDENT), HEAD
         WRITE (LUPRI, '(120A)') (' ',I=1,INDENT), (A,I=1,LENGTH)
      END IF
      WRITE (LUPRI, '()')
      CALL FLSHFO(LUPRI)
      RETURN
      END
C  /* Deck around */
      SUBROUTINE AROUND(HEAD)
      CHARACTER HEAD*(*)
#include "priunit.h"
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
#include "implicit.h"
#include "priunit.h"
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
#include "implicit.h"
#if defined (VAR_MPI)
#include "infpar.h"
      INCLUDE 'mpif.h'
#endif
#include "priunit.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#endif
      CHARACTER TEXT*(*), STRING*70
C
#if defined (VAR_MPI)
      IF (MYTID .GT. 0) THEN
         WRITE (STRING,'(A,I4,A)') ' Slave node',MYTID,
     &   ' :  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
      ELSE
         WRITE (STRING,'(A)')
     &   '  Master node : --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
      END IF
#else
      WRITE (STRING,'(A)')
     &   '  --- SEVERE ERROR, PROGRAM WILL BE ABORTED ---'
#endif
      WRITE (LUPRI,'(/A)') STRING
      CALL TIMDAT
      WRITE (LUPRI,'(1X,A)') TEXT
#if defined (SYS_CRAY) || defined (SYS_CONVEX)\
|| defined (SYS_UNIX) || defined (SYS_AIX)\
|| defined (SYS_HPUX) || defined (SYS_IRIX)\
|| defined (SYS_NEC) || defined (SYS_DEC)\
|| defined (SYS_T3D) || defined (SYS_LINUX)\
|| defined (SYS_SUN) || defined (SYS_HPUX)\
|| defined (SYS_HAL) || defined (SYS_T90)\
|| defined (SYS_SX) || defined (SYS_DARWIN)
C     Write to stderr
      WRITE (0,'(/A/1X,A)') STRING, TEXT
#endif
      CALL FLSHFO(LUPRI)
#if defined (VAR_MPI)
#   if defined (PRG_DIRAC)
      IF(PARCAL) CALL MPI_ABORT(MPI_COMM_WORLD,IERR,IERR)
#   else
      CALL MPI_ABORT(MPI_COMM_WORLD,IERR,IERR)
#   endif
#endif
#if defined (VAR_PVM)
      CALL PVMFEXIT(INFO)
#endif
      CALL TRACE
      STOP 1000
      END
C  /* Deck opendx */
      SUBROUTINE OPENDX (LUDX,NAME,NELEM,STATUS,OLDDX)
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
C  STATUS   file status: 'OLD', 'NEW', or 'UNKNOWN'
C
C Output:
C  OLDDX    logical, true if old LUDX file was opened
C
C
      CHARACTER*(*) NAME, STATUS
      LOGICAL OLDDX, FEXIST
#include "priunit.h"
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
C
C     CRAY has 8 byte integers.
C
      LRECL  = 8*NELEM
C
      IF (STATUS .EQ. 'NEW') THEN
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
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
     *        ACCESS='DIRECT',ERR=301,IOSTAT=IOS,RECL=LRECL)
         OLDDX = .FALSE.
C
C     ERROR HANDLING (JAA 1997-06-06)
C
 301     IF (IOS .GT. 0) THEN
            WRITE(102,*) '*** ERROR (OPENDX) NEW FILE OPEN FAILED'
            WRITE(102,*) '    FILE=',NAME,', UNIT=',LUDX,
     $           ', IOSTAT=',IOS
            IF (IOS .LT. 1000) THEN
               WRITE(102,'(A,I3.3,A)') '    (try explain sys-',IOS,')'
            ELSE IF (IOS .LT. 6000) THEN
               WRITE(102,'(A,I4,A)') '    (try explain lib-',IOS,')'
            ENDIF
            CALL FLUSH(102,IOS)
            CALL TRBK()
#if defined(SYS_T3D)
            CALL STOP_ALL()
#else
            CALL QUIT('OPENDX: new file open failed')
#endif
         ENDIF
  600 CONTINUE
#endif
#if defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_HPUX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_NEC) || defined (SYS_LINUX) || defined (SYS_SUN) || defined (SYS_HAL) || defined (SYS_DARWIN)
C
C     Assume same for CONVEX --  1-Nov-1989 hjaaj
C     Assume same for IBM-AIX--  1-Oct-1990 hjaaj
C     Assume same for HP-UX  -- 21-Aug-1991 hjaaj
C     Assume same for DEC    -- 21-May-1992 hjaaj
C     Assume same for IRIX   --  3-Feb-1994 hjaaj; except recl in *4 words
C     Assume same for PARAGON-- 13-Oct-1994 hjaaj
C
#if defined (SYS_IRIX) || defined (SYS_NEC)
      LRECL  = NELEM
#else
      LRECL  = 4*NELEM
#endif
C
      IF (STATUS .EQ. 'NEW') THEN
C        some systems abort if we try to open as NEW a file,
C        which already exist
         INQUIRE(FILE=NAME,EXIST=FEXIST)
         IF (FEXIST) THEN
            OPEN(UNIT=LUDX,STATUS='OLD',FILE=NAME)
            CLOSE(UNIT=LUDX,STATUS='DELETE')
         END IF
         GO TO 300
      END IF
C
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
#if !defined (SYS_CRAY) && !defined (SYS_CONVEX)\
&& !defined (SYS_AIX) && !defined (SYS_HPUX)\
&& !defined (SYS_DEC) && !defined (SYS_IRIX)\
&& !defined (SYS_NEC) && !defined (SYS_T3D)\
&& !defined (SYS_LINUX) && !defined (SYS_SUN)\
&& !defined (SYS_HAL) && !defined (SYS_T90)\
&& !defined(SYS_DARWIN)
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
      LOGICAL FUNCTION FINDDX(LU,I,LEN,IVEC)
C
C 27-Jun-1985 Hans Jorgen Aa. Jensen
C
C For direct access find record,
C when LEN may be greater than maximum record length.
C
      INTEGER IVEC(LEN)
      READ (LU, REC=I, IOSTAT=IOS) IVEC
      IF (IOS .NE. 0) GO TO 900
      FINDDX = .TRUE.
      RETURN
C
  900 CONTINUE
      FINDDX = .FALSE.
      RETURN
      END
C  /* Deck readdx */
      SUBROUTINE READDX(LU,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
C For direct access read
C
#include "implicit.h"
#include "priunit.h"
#include "2gbdef.h"
#include "dummy.h"
      LOGICAL OLDDX
      CHARACTER*16 FNNAME, FNNM2
      INTEGER IVEC(LEN)
#include "inftap.h"
#include "chrnos.h"
C
      IF (LEN .EQ. 0) RETURN
#if defined (VAR_SPLITFILES)
      INQUIRE(UNIT=LU,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      LRECL = LRECL/8
#else
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
      LRECL = LRECL/4
#endif
#endif
      IRECPF = I2GB/LRECL
C
C     Error branch for very large record lengths
C
      IF (IRECPF .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I14,A)') 
     &        ' ERROR: The requested record length of ', LRECL, 
     &        ' is larger than 2Gb'
         WRITE (LUPRI,'(/A)') ' The program cannot write even a '//
     &        'single record to file!'
         CALL QUIT('Record length of file too long for this computer')
      END IF
C
      IF (I .GT. IRECPF) THEN
C
C     Ooops, this record is not in this file
C
         LFILNM = I/IRECPF - 1
         IF (LFILNM .GT. 9) THEN
            WRITE (LUPRI,'(/A)') ' DALTON needs to read from a file '//
     &           ' split more than 11 times.',
     &           ' This is currently not supported'
            CALL QUIT('Too many splittings of a file')
         END IF
         IRECT = MOD(I,IRECPF)
         IF (IRECT .EQ. 0) IRECT = IRECPF
         INQUIRE(UNIT=LU,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(LFILNM)
         CALL GPOPEN(LU,FNNM2(1:(LN+2)),'UNKNOWN','DIRECT',' ',LRECL,
     &               OLDDX)
      ELSE
         IRECT = I
      END IF
#else
      IRECT = I
#endif
      READ (LU, REC = IRECT) IVEC
#if defined (VAR_SPLITFILES)
      IF (I .GT. IRECPF) THEN
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'UNKNOWN','DIRECT',' ',LRECL,OLDDX)
      END IF
#endif
      RETURN
      END
C  /* Deck writdx */
      SUBROUTINE WRITDX(LU,I,LEN,IVEC)
C
C 30-Apr-1985 Hans Jorgen Aa. Jensen
C
#include "implicit.h"
#include "dummy.h"
#include "2gbdef.h"
#include "priunit.h"
      LOGICAL OLDDX
      CHARACTER*16 FNNAME, FNNM2
      INTEGER IVEC(LEN)
#include "inftap.h"
#include "chrnos.h"
C
#if defined (VAR_SPLITFILES)
      INQUIRE(UNIT=LU,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
      LRECL = LRECL/8
#else
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
      LRECL = LRECL/4
#endif
#endif
      IRECPF = I2GB/LRECL
C
C     Error branch for very large record lengths
C
      IF (IRECPF .EQ. 0) THEN
         WRITE (LUPRI,'(/A,I14,A)') 
     &        ' ERROR: The requested record length of ', LRECL, 
     &        ' is larger than 2Gb'
         WRITE (LUPRI,'(/A)') ' The program cannot write even a '//
     &        'single record to file!'
         CALL QUIT('Record length of file too long for this computer')
      END IF
C
      IF (I .GT. IRECPF) THEN
C
C     Ooops, this record will not fit in this file
C
         LFILNM = I/IRECPF - 1
         IF (LFILNM .GT. 9) THEN
            WRITE (LUPRI,'(/A)') ' DALTON need to split a file more '//
     &           ' than 11 times.', ' This is currently not supported'
            CALL QUIT('Too many splittings of a file')
         END IF
         IRECT = MOD(I,IRECPF)
         IF (IRECT .EQ. 0) IRECT = IRECPF
         INQUIRE(UNIT=LU,NAME=FNNAME)
         LN = 1
 10      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 10
         END IF
         LN = LN - 1
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(LFILNM)
         LN2 = LN + 2
         CALL GPOPEN(LU,FNNM2(1:LN2),'UNKNOWN','DIRECT',' ',LRECL,
     &               OLDDX)
      ELSE
         IRECT = I   
      END IF
#else
      IRECT = I
#endif
      WRITE (LU, REC = IRECT) IVEC
#if defined (VAR_SPLITFILES)
      IF (I .GT. IRECPF) THEN
         LUBKP = LU
         CALL GPCLOSE(LU,'KEEP')
         LU = LUBKP
         CALL GPOPEN(LU,FNNAME(1:LN),'UNKNOWN','DIRECT',' ',LRECL,OLDDX)
      END IF
#endif
      RETURN
      END
C
C panor: The 48 bit shift done below:
C        IINDX4(1,I) = IBTAND(IBTSHR(LABEL,48),IBT16)
C        genereates a compilation error on SUN Solaris 32 bit
C        architectures. The subroutine is commented out because 
C        it is presently not used in Dirac. 2004-08-28
C
C  /* Deck aolab4 */
C       SUBROUTINE AOLAB4(IINDPK,NMAX,NIBUF,NBITS,IINDX4,N)
C C
C C     Written by Henrik Koch 22-Nov-1991 (as AOLABE)
C C     Generalized 25-Oct-1993 hjaaj
C C                 31-Jul-1996 hjaaj (new NIBUF parameter)
C C
C C     Unpack 4 integer indices packed in NIBUF *4 integer(s).
C C
C #include "implicit.h"
C       DIMENSION IINDPK(*), IINDX4(4,*)
C #include "ibtdef.h"
C C
C       N = IINDPK(NIBUF*NMAX+1)
C C
C       IF (NIBUF .EQ. 1) THEN
C       IF (NBITS .EQ. 8) THEN
C          DO 100 I = 1,N
C             LABEL = IINDPK(I)
C             IINDX4(1,I) = IBTAND(IBTSHR(LABEL,24),IBT08)
C             IINDX4(2,I) = IBTAND(IBTSHR(LABEL,16),IBT08)
C             IINDX4(3,I) = IBTAND(IBTSHR(LABEL, 8),IBT08)
C             IINDX4(4,I) = IBTAND(       LABEL,    IBT08)
C   100    CONTINUE
C       ELSE IF (NBITS .EQ. 16) THEN
C          DO 200 I = 1,N
C             LABEL = IINDPK(I)
C             IINDX4(1,I) = IBTAND(IBTSHR(LABEL,48),IBT16)
C             IINDX4(2,I) = IBTAND(IBTSHR(LABEL,32),IBT16)
C             IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
C             IINDX4(4,I) = IBTAND(       LABEL,    IBT16)
C   200    CONTINUE
C       ELSE
C 	 CALL QUIT('AOLAB4 error')
C       END IF
C       ELSE
C          DO 300 I = 1,N
C             LABEL = IINDPK(I)
C             IINDX4(1,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
C             IINDX4(2,I) = IBTAND(       LABEL    ,IBT16)
C             LABEL = IINDPK(NMAX+I)
C             IINDX4(3,I) = IBTAND(IBTSHR(LABEL,16),IBT16)
C             IINDX4(4,I) = IBTAND(       LABEL    ,IBT16)
C   300    CONTINUE
C       END IF
C C
C       RETURN
C       END
C  /* Deck GPOPEN */
      SUBROUTINE GPOPEN(LUNIT,FILEIN,STATIN,ACCEIN,FORMIN,NELEM,OLDDX)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for opening files in the Dalton program.
C   The routine will dynamically allocate unit numbers that will
C   become available again when the file is closed using GPCLSE
C
C   Direct access files are provided with a UNIT-number, but the rest
C   of the file opening process is taken care of by a call to OPENDX.
C
C   These files are strongly machine dependent, although care has been 
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton 
C   program
C
C Input:
C  LUNIT    Suggested unit number (OPTIONAL, but is mandatory if an unnamed
C           file is reopened after having been closed with STATUS='KEEP')
C           Otherwise it is recommended to not assign this.
C  FILEIN   Suggested name for the file (OPTIONAL, but strongly recommended)
C  STATIN   Suggested status of the file (OPTIONAL and maybe not recommended)
C  ACCEIN   Access method. 'DIRECT' or 'SEQUENTIAL'. Default is sequential.
C  FORMIN   Formatted or unformatted file format. Default is 'UNFORMATTED'
C
C Direct access-specific input:
C
C  NELEM    number of integer words per logical record
C  NREC     number of logical records
C
C
C Output:
C  LUNIT    Assigned file unit number
C
C  Direct access-specific output (through 
C  LRDX     physical record length (in integers)
C  NBDX     number of physical records per logical record
C  OLDDX    logical, true if old LUDX file was opened
C
C
#include "implicit.h"
#include "priunit.h"
#include "chrnos.h"
      CHARACTER*80 FILENM, STATUF, ACCESF, FORMF
      CHARACTER*(*) FILEIN, STATIN, ACCEIN, FORMIN
      INTEGER FILELN, STATLN, ACCELN, FORMLN
      LOGICAL OLDDX, FEXST, FIRST
      COMMON /UNITAR/ IUNTAB(1:99)
      DATA FIRST /.TRUE./
      SAVE FIRST
C
C     Initialization
C     
      IF (FIRST) THEN
         CALL IZERO(IUNTAB,99)
         FIRST = .FALSE.
      END IF
C
C     We first deal with the unit number
C
      FILELN = LEN(FILEIN)
      STATLN = LEN(STATIN)
      ACCELN = LEN(ACCEIN)
      FORMLN = LEN(FORMIN)
      FILENM(1:FILELN) = FILEIN(1:FILELN)
      STATUF(1:STATLN) = STATIN(1:STATLN)
      ACCESF(1:ACCELN) = ACCEIN(1:ACCELN)
      FORMF (1:FORMLN) = FORMIN(1:FORMLN)
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99)) THEN
C
C     Unit number left unassigned, we get to decide!
C
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
         IF (LUNIT .GT. 99) GOTO 9001
      ELSE
C
C     The user has requested a specific unit number. We don't 
C     quite trust the user, so we check that
C     1) It is not unit 5, or 6
C     2) The file either has been closed with status='KEEP' or is
C        not currently in use
C
         IF ((LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9002
         IF (IUNTAB(LUNIT) .EQ. 1) GOTO 9003
      END IF
C
C     We've got a file number now, get rid of the direct access case
C
      IF (ACCESF(1:6) .EQ. 'DIRECT') THEN
         IF (FILELN .EQ. 1) THEN
            FILELN = 6
            FILENM(1:6) = 'UNIT'//CHRNOS(LUNIT/10)//
     &                    CHRNOS(MOD(LUNIT,10))
         END IF
         IF (STATLN .EQ. 1) THEN
            STATLN = 7
            STATUF(1:7) = 'UNKNOWN'
         END IF
         CALL OPENDX(LUNIT,FILENM(1:FILELN),NELEM,STATUF(1:STATLN),
     &               OLDDX)
      ELSE
         IF (STATLN .EQ. 1) THEN
            STATLN = 7
            STATUF(1:7) = 'UNKNOWN'
         END IF
         IF (STATUF(1:3) .EQ. 'OLD' .AND. IUNTAB(LUNIT) .EQ. 0) THEN
C
C     This better be a file with a name, and it better exist
C
            IF (FILELN .EQ. 1) GOTO 9005
            INQUIRE(FILE=FILENM(1:FILELN),EXIST=FEXST,IOSTAT=IOS,
     &              ERR=9006)
            IF (.NOT. FEXST) GOTO 9007
            IF (FORMLN .EQ. 1) THEN
               FORMLN = 9
               FORMF(1:FORMLN)='FORMATTED'
            END IF
            OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),STATUS='OLD',
     &           ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),ERR=9004,
     &           IOSTAT=IOS)
         ELSE
            IF (STATUF(1:3) .EQ. 'NEW' .AND. IUNTAB(LUNIT) .NE. 0)
     &           GOTO 9006
            IF (FILELN .EQ. 1) THEN
               FILELN = 6
               FILENM(1:6) = 'UNIT'//CHRNOS(LUNIT/10)//
     &                               CHRNOS(MOD(LUNIT,10))
            END IF
            IF (FORMLN .EQ. 1) THEN
               FORMLN = 11
               FORMF(1:11)='UNFORMATTED'
            END IF
            IF (STATUF(1:7) .EQ. 'SCRATCH') THEN
               OPEN(UNIT=LUNIT,STATUS='SCRATCH',
     &              ACCESS='SEQUENTIAL',FORM=FORMF(1:FORMLN),
     &              ERR=9004,IOSTAT=IOS)
            ELSE
               OPEN(UNIT=LUNIT,FILE=FILENM(1:FILELN),
     &              STATUS=STATUF(1:STATLN),ACCESS='SEQUENTIAL',
     &              FORM=FORMF(1:FORMLN),ERR=9004,IOSTAT=IOS)
            END IF
         END IF
      END IF
      IUNTAB(LUNIT) = 1
      RETURN
C     
C error branches
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) NO MORE AVAILABLE FILENUMBERS! ',
     &     ' *** THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &     'FILES OR', ' *** SOMEBODY HAS FORGOT TO CLOSE FILES IN'//
     &     'THE SOURCE CODE'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) NO MORE FILE NUMBERS')
C
 9002 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN AN ILLEGAL FILE NUMBER',
     &     ' *** SOMEBODY HAS TRIED TO OPEN UNITS 0, 3, 5, OR 6',
     &     ' *** THE PROGRAM DOES NOT ALLOW THE USE OF THESE RESERVED'//
     &     ' UNIT NUMBERS'
      WRITE (6,'(//A,I5)') ' Offending UNIT number: ',LUNIT
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) ILLEGAL FILE NUMBER REQUESTED')
C
 9003 CONTINUE
      WRITE (6,'(//A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN A FILE ALREADY IN USE',
     &     ' *** SOMEBODY IS TRYING TO USE A FILENUMBER THAT IS '//
     &     'ALREADY IN USE'
      WRITE (6,'(//A,I5)') ' Offending UNIT number: ',LUNIT
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) TRYING TO OPEN A FILE ALREADY '//
     &          'IN USE')
C
 9004 CONTINUE
      WRITE (6,'(//A,I3/A,I7)')
     &     ' *** ERROR (GPOPEN) UPON TRYING TO OPEN FILE ON UNIT',
     &     LUNIT, ' *** ERROR CODE RETURNED ',IOS
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) UPON OPENING A FILE')
C
 9005 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &     ' *** A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE TABLE',
     &     ' *** ENTRY INDICATE THAT IT DOES NOT EXIST'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) TRYING TO OPEN A NON_EXISTING '//
     &          'FILE AS OLD')
C
 9006 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN AN EXISTING NEW FILE',
     &     ' *** A FILE HAS BEEN SPECIFIED TO BE NEW, BUT THE TABLE',
     &     ' *** ENTRY INDICATE THAT IT ALREADY EXISTS'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) TRYING TO OPEN AN EXISTING '//
     &          'FILE AS NEW')
C
 9007 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPOPEN) TRYING TO OPEN A NON-EXISTING OLD FILE',
     &     ' *** A FILE HAS BEEN SPECIFIED TO BE OLD, BUT THE FILE',
     &     ' *** DOES NOT EXIST'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPOPEN) TRYING TO OPEN A NON_EXISTING '//
     &          'FILE AS OLD')
C
C end of GPOPEN
C
      END
C  /* Deck gpclse */
      SUBROUTINE GPCLOSE(LUNIT,DISP)
C
C 23-Feb-2000 K.Ruud
C
C Purpose:
C   General purpose routine for closing files in the Dalton program.
C   The routine will ensure that files that are closed and that will not
C   be reopened in a later part of the program will make their unit 
C   number available for reuse.
C
C   These files are strongly machine dependent, although care has been 
C   taken to avoid using unit numbers that are illegal or reserved on a
C   particular architecture. However, the use of GPOPEN and GPCLSE should
C   remove much of the machine dependence in the rest of the Dalton 
C   program
C
C Input:
C  LUNIT    Mandatory
C  STATUS   Indicate whether the file should be removed or kept
C           ('KEEP' or 'DELETE'). Note that the default is that the file
C           is to be deleted
C
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) DISP
      CHARACTER*16 FNNAME, FNNM2
      CHARACTER*3 IODIR
      LOGICAL FEXIST, OPND, OLDDX
      COMMON /UNITAR/ IUNTAB(1:99)
#include "chrnos.h"
C
C     We first deal with the unit number
C
      IF ((LUNIT .LT. 1) .OR. (LUNIT .GT. 99) .OR.
     &    (LUNIT .EQ. 5) .OR. (LUNIT .EQ. 6)) GOTO 9001
C
C     Check that the file actually have been opened
C
      IF (IUNTAB(LUNIT) .EQ. 0) GOTO 9002
C
      INQUIRE(UNIT=LUNIT,EXIST=FEXIST,OPENED=OPND)
      IF (.NOT. FEXIST .OR. .NOT. OPND) WRITE (LUPRI,'(/A/A,I3)')
     &     ' WARNING: Tried to close a non-existent or already '//
     &     'closed file', ' Unit number was: ',LUNIT
      CALL FLSHFO(LUPRI)
      IF (DISP .EQ. 'KEEP') THEN
         CLOSE (LUNIT,STATUS='KEEP',ERR=9003)
         GOTO 30
      ELSE
#if defined (VAR_SPLITFILES)
         INQUIRE(UNIT=LUNIT,DIRECT=IODIR)
         INQUIRE(UNIT=LUNIT,NAME=FNNAME)
         IF (IODIR .EQ. 'YES') THEN
            INQUIRE(UNIT=LUNIT,RECL=LRECL)
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (SYS_T90)
            LRECL = LRECL/8
#endif
#if !defined (SYS_NEC) && !defined (SYS_IRIX)
            LRECL = LRECL/4
#endif
         END IF
         LN = 1
 15      CONTINUE
         IF (FNNAME(LN:LN) .NE. ' ') THEN
            LN = LN + 1
            GOTO 15
         END IF
         LN = LN - 1
         IF (IODIR .NE. 'YES') THEN
            IF (FNNAME((LN-1):(LN-1)) .EQ. '-') THEN
               CLOSE(LUNIT,STATUS='KEEP')
               LN = LN - 2
               OPEN (UNIT=LUNIT,FILE=FNNAME(1:LN),STATUS='OLD')
            END IF
         END IF
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
         IUNIT = 0
 16      CONTINUE
         FNNM2 = FNNAME(1:LN)//'-'//CHRNOS(IUNIT)
         LN2 = LN + 2
         INQUIRE(FILE=FNNM2(1:LN2),EXIST=FEXIST)
         IF (FEXIST) THEN
            IF (IODIR .EQ. 'YES') THEN
               CALL OPENDX(LUNIT,FNNM2(1:LN2),LRECL,'OLD',OLDDX)
            ELSE
               OPEN(UNIT=LUNIT,FILE=FNNM2(1:LN2),STATUS='OLD')
            END IF
            CLOSE(LUNIT,STATUS='DELETE')
            IUNIT = IUNIT + 1
            GOTO 16
         END IF
#else
         CLOSE (LUNIT,STATUS='DELETE',ERR=20)
#endif
         GOTO 30
      END IF
C
C     We might get an error for trying to close a SCRATCH file with STATUS
C     We just close it
C     
 20   CLOSE (LUNIT)
C
C     We release that saved unit number by resetting it to 0
C
 30   CONTINUE
      IUNTAB(LUNIT) = 0
      LUNIT = 0
      RETURN
C     
C error branches
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A)')
     &    ' *** ERROR (GPCLOSE) TRYING TO CLOSE AN ILLEGAL FILE NUMBER',
     &    ' *** SOMEBODY HAS TRIED TO CLOSE UNITS 0, 3, 5, OR 6',
     &    ' *** THE PROGRAM DOES NOT ALLOW THE USE OF THESE RESERVED'//
     &    ' UNIT NUMBERS'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPCLOSE) ILLEGAL FILE NUMBER REQUESTED')
C
 9002 CONTINUE
      WRITE (6,'(//A/A)')
     &     ' *** ERROR (GPCLOSE) TRYING TO CLOSE A FILE NOT IN USE',
     &     ' *** SOMEBODY IS TRYING TO USE A FILENUMBER THAT HAS '//
     &     'NOT BEEN USED'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPCLOSE) TRYING TO CLOSE A FILE NOT '//
     &          'IN USE')
C
 9003 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (GPCLOSE) TRYING TO KEEP A SCRATCH FILE',
     &     ' *** A FILE HAS BEEN INDICATED TO BE CLOSED AND KEPT,',
     &     ' *** BUT IT APPEARS THE FILE IS A SCRATCH FILE'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dirac-admin@dirac.chem.sdu.dk'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (GPCLOSE) TRYING TO CLOSE A SCRATCH FILE ')
C
C end of GPCLOSE
C
      END
#if ! defined (PRG_DIRAC)
C  /* Deck gpioin */
      SUBROUTINE GPIOIN
C
C 1-Mar-2000 K.Ruud
C
#include "priunit.h"
#include "inftap.h"
C
      LUCMD  = 0
      LUPRI  = 0
C
      LUCME  = 0
      LUMOL  = 0
      LUPROP = 0
      LUSOL  = 0
      LUINTA = 0
      LUONEL = 0
      LUSUPM = 0
      LUTLM  = 0
      LUDA1  = 0
      LUITMP = 0
      LU2DER = 0
      LUDASP = 0
      LUAAT  = 0
      LURDR  = 0
      LURDI  = 0
      LUGDR  = 0
      LUGDI  = 0
      LUGDT  = 0
      LURDT  = 0
      LUDFCK = 0
      LUSFDA = 0
      LUFDC  = 0
      LUWLK  = 0
      LUPAO  = 0
      LUPAS  = 0
      LUNR1  = 0
      LUNR3  = 0
      LUNR5  = 0
      LUCMD  = 0
      LUINTR = 0
      LUAORC = 0
      LUMOLDEN=0
C
C     /RSPTAP/
C
      LUAHSO = 0
      LUCRV1 = 0
      LUCRV2 = 0
      LUXYVE = 0
      LUCRVE = 0
      LUQRVE = 0
      LURSP3 = 0
      LURSP4 = 0
      LURSP5 = 0
      LUMHSO = 0
C
C     /SIRTAP/
C
      LUINTM = 0
      LUIT1  = 0
      LUINF = 0
      LUSIFC = 0
C
      FNSOL  = 'AOSOLINT'
      ABARDR = 'ABACUS.RD'
      ABARDI = 'ABACUS.RDI'
      ABAGDR = 'ABACUS.GD'
      ABAGDI = 'ABACUS.GDI'
      ABAGDT = 'ABACUS.GDT'
      ABARDT = 'ABACUS.RDT'
      ABADFK = 'ABACUS.DFK'
      ABASF  = 'ABACUS.SF'
      ABATLM = 'ABACUS.TLM'
      ABAWLK = 'ABACUS.WLK'
      ABANR1 = 'ABAENR.RST'
      ABANR3 = 'ABAENR.BVC'
      ABANR5 = 'ABAENR.SVC'
      FNINTM = 'MOTWOINT'
      FNSUPM = 'AOSUPINT'
      FNONEL = 'AOONEINT'
      FNSIFC = 'SIRIFC'
      LBSIFC = 'SIR IPH '
C
      RETURN
      END
C /* Deck daintf */
      SUBROUTINE DAINTF(LUNIT,NAME)
C
C     Interface for assigning a unit number to a unit before opening it in
C     DANAME. Needed to ensure compatibility with use of GPOPEN/GPCLOSE
C     K.Ruud, March 28 1999
C
#include "implicit.h"
      CHARACTER NAME*(*)
      COMMON /UNITAR/ IUNTAB(1:99)
C      
      IF (LUNIT .EQ. 0) THEN
         IUN = 0
 10      CONTINUE
         IUN = IUN + 1
         IF ((IUN .EQ. 5) .OR. (IUN .EQ. 6))
     &        GOTO 10
         IF (IUNTAB(IUN) .NE. 0) GOTO 10
         LUNIT = IUN
         IF (LUNIT .GT. 99) GOTO 9001
      END IF
      CALL DANAME(LUNIT,NAME)
      RETURN
C     
C error branch
C
 9001 CONTINUE
      WRITE (6,'(//A/A/A)')
     &     ' *** ERROR (DAINTF) NO MORE AVAILABLE FILENUMBERS! ',
     &     ' *** THIS CALCULATION EITHER NEEDS TOO MANY SIMULTANEOUS '//
     &     'FILES OR', ' *** SOMEBODY HAS FORGOT TO CLOSE FILES IN'//
     &     'THE SOURCE CODE'
      WRITE (6,'(/A)') ' ### Please report the problem to '//
     &     'dalton-admin@kjemi.uio.no'
      CALL QTRACE(6)
      CALL QUIT('*** ERROR (DAINTF) NO MORE FILE NUMBERS')
C
      END
#endif
