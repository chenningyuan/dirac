C  /* Deck getdat */
      SUBROUTINE GETDAT(CDATE,CTIME)
C
C     24-Jan-1988 Hans Joergen Aa. Jensen
C
C     Return date and time as character*8, for labels.
C
      CHARACTER*(8) CDATE, CTIME
C
C     The fortran90 variables
C
      character*8 date_f90,time_f90*10,zone_f90*5
      integer values(8)
CLV : Kept the old code, we may consider calling always the f90 implementation...
#if defined (SYS_CRAY) || defined (SYS_T3D)
      INTEGER IDATE, ITIME, DATE, CLOCK
      IDATE = DATE()
      ITIME = CLOCK()
      WRITE (CDATE,'(A8)') IDATE
      WRITE (CTIME,'(A8)') ITIME
#else
#if defined (SYS_AIX) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) \
   || defined (SYS_DEC) || defined (SYS_IRIX)  || defined (SYS_PARAGON) \
   || defined (SYS_SUN) || defined (SYS_LINUX) || defined(SYS_HPUX) \
   || defined (SYS_SX) || defined (SYS_DARWIN)
      CHARACTER*(24) FFDATE, FDATE
      FFDATE = FDATE()
      CDATE  = FFDATE(9:10)//FFDATE(5:7)//FFDATE(23:24)//' '
      CTIME  = FFDATE(12:19)
#else
#if defined (SYS_IBM)
      CHARACTER DAOFWK*9
      CALL DATE(DAOFWK,CDATE)
      CALL TIME(CTIME)
#else
C
C     LV : Use fortran90 implementation
      CALL DATE_AND_TIME(date_f90,time_f90,zone_f90,values)
      CDATE = DATE_F90(1:8)
      CTIME = TIME_F90(1:2)//':'//TIME_F90(3:4)//':'//TIME_F90(5:6)
#endif
#endif
#endif
      RETURN
      END
C  /* Deck gettim */
      SUBROUTINE GETTIM(TIMCPU,TIMWAL)
C
C Written 18-Dec-1984 Hans Joergen Aa. Jensen
C Rewrite Feb 6 2014 Trond Saue
C
C Return elapsed CPU time and elapsed real time.
C
      real(8) :: timcpu,timwal
      integer :: clock_count, clock_rate, clock_max
      call cpu_time (timcpu)
      call system_clock ( clock_count, clock_rate, clock_max )
      timwal = real ( clock_count ) / real ( clock_rate )
      return
      end
C  /* Deck settim */
#if defined (SYS_IBM)
      SUBROUTINE SETTIM
C
C IBM timing routine 24-Jun-1988 hjaaj
C (emulating Lund IBM timing routines)
C
      REAL*4 ECPU, ETIME, ETCPU
      REAL*4 SCPU, STIME, STCPU
      REAL*4 PCPU, PTIME, PTCPU
      REAL*8 CPUE, CPUA, TIOE, TIOA
      CHARACTER*23 DATTIM
      SAVE   SCPU, STIME, STCPU
      SAVE   PCPU, PTIME, PTCPU
C
      CALL DATETM(DATTIM,SCPU,STIME,STCPU)
      PCPU  = SCPU
      PTIME = STIME
      PTCPU = STCPU
      RETURN
C
C     ----- entry timing -----
C
      ENTRY TIMING(CPUA,CPUE,TIOA,TIOE)
      CALL DATETM(DATTIM,ECPU,ETIME,ETCPU)
      CPUA   = ETCPU - STCPU
      CPUE   = ETCPU - PTCPU
      PTCPU  = ETCPU
      TIOA   = ETIME - STIME
      TIOE   = ETIME - PTIME
      PTIME  = ETIME
C
      RETURN
      END
#endif
#if !defined (VAR_MOTECC) && !defined (SYS_IBM)
      SUBROUTINE SETTIM()
#include "implicit.h"
      PARAMETER ( D0 = 0.0D0 )
      PRINT *,'DUMMY SETTIM called'
      RETURN
      ENTRY TIMING(CPUA,CPUE,TIOA,TIOE)
      CPUA = D0
      CPUE = D0
      TIOA = D0
      TIOE = D0
      PRINT *,'DUMMY TIMING called'
      RETURN
      END
#endif
C  /* Deck timtxt */
      SUBROUTINE TIMTXT(TEXT,TIMUSD,LUPRIN)
C
C TIMTXT based on TIMER by TUH //900709-hjaaj
C
#include "implicit.h"
      CHARACTER*(*) TEXT
      CHARACTER AHOUR*6, ASEC*8, AMIN*8
C
      ISECND = NINT(TIMUSD)
      IF (ISECND .GE. 60) THEN
         MINUTE = ISECND/60
         IHOURS = MINUTE/60
         MINUTE = MINUTE - 60*IHOURS
         ISECND = ISECND - 3600*IHOURS - 60*MINUTE
         IF (IHOURS .EQ. 1) THEN
            AHOUR = ' hour '
         ELSE
            AHOUR = ' hours'
         END IF
         IF (MINUTE .EQ. 1) THEN
            AMIN = ' minute '
         ELSE
            AMIN = ' minutes'
         END IF
         IF (ISECND .EQ. 1) THEN
            ASEC = ' second '
         ELSE
            ASEC = ' seconds'
         END IF
         IF (IHOURS .GT. 0) THEN
            WRITE(LUPRIN,100)
     *            TEXT, IHOURS, AHOUR, MINUTE, AMIN, ISECND, ASEC
         ELSE
            WRITE(LUPRIN,200) TEXT, MINUTE, AMIN, ISECND, ASEC
         END IF
      ELSE
         WRITE(LUPRIN,300) TEXT,TIMUSD
      END IF
  100 FORMAT(1X,A,I4,A,I3,A,I3,A)
  200 FORMAT(1X,A,     I3,A,I3,A)
  300 FORMAT(1X,A,F7.2,' seconds')
      RETURN
      END
C  /* Deck tstamp */
      SUBROUTINE TSTAMP(TEXT,LUPRIN)
C
C Copyright Hans Joergen Aa. Jensen 9-Jul-1990
C
C Purpose: To stamp as many as possible of
C          text, date, time, computer, and hostname to LUPRIN
C
#include "implicit.h"
      CHARACTER*(*) TEXT
C
#if defined (SYS_IBM)
      CHARACTER CDATE*8,CTIME*8,DAOFWK*9
#endif
#if defined (SYS_UNIX) || defined (SYS_AIX) || defined (SYS_ALLIANT) ||defined (SYS_CONVEX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_PARAGON) || defined (SYS_SUN) || defined (SYS_SX) || defined (SYS_LINUX)
      CHARACTER*(24) FDATE
      CHARACTER*(40) HSTNAM
#endif
C
      LTEXT = LEN(TEXT)
      IF (LTEXT .GT. 0) THEN
         WRITE (LUPRIN,'(/A)') TEXT
      ELSE
         WRITE (LUPRIN,'()')
      END IF
#if defined (SYS_IBM)
      CALL DATE(DAOFWK,CDATE)
      CALL TIME(CTIME)
      WRITE (LUPRIN,'(T6,A,A8,1X,A8)')
     *   'Date and time (IBM) : ',CDATE,CTIME
#endif
#if defined (SYS_CRAY) || defined (SYS_T3D)
      CDATE = DATE()
      CTIME = CLOCK()
      WRITE (LUPRIN,'(T6,A,A8,1X,A8)')
     *   'Date and time (CRAY) : ',CDATE,CTIME
#endif
#if 0 && defined (SYS_HPUX)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (HP-UX)  : ',FDATE()
#endif
#if defined (SYS_DEC)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (DEC)    : ',FDATE()
#endif
#if defined (SYS_IRIX)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (IRIX)   : ',FDATE()
#endif
#if defined (SYS_SUN)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (SUN)    : ',FDATE()
#endif
#if defined (SYS_LINUX)
      WRITE (LUPRIN,'(T6,2A)') 'Date and time (Linux)  : ',FDATE()
#endif
#if defined (SYS_UNIX) || defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_HPUX) || defined (SYS_DEC) || defined (SYS_IRIX) || defined (SYS_SUN) || defined (SYS_SX) || defined (SYS_LINUX)
#ifndef VAR_EKO
C
C     CALL HOSTNM(HSTNAM)
C
      CALL GETENV('HOSTNAME',HSTNAM)
      if (HSTNAM=="") HSTNAM="not found"
      WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
#endif
#if defined (SYS_AIX)
!C 930414-hjaaj: apparent error IBM's xlf library routines:
!C    when 'T6' then column 1-5 not blanked but contains text
!C    from a previous print statement!
!C     AIX XL FORTRAN version 2.3+
      WRITE (LUPRIN,'(2A)') '     Date and time (IBM-AIX): ',fdate()
C     WRITE (LUPRIN,'(T6,2A)') 'Date and time (IBM-AIX): ',fdate_()
C     CALL hostnm_(HSTNAM)
      CALL hostnm(HSTNAM)
      if (HSTNAM=="") HSTNAM="not found"
      WRITE (LUPRIN,'(2A)') '     Host name              : ',HSTNAM
C     WRITE (LUPRIN,'(T6,2A)') 'Host name              : ',HSTNAM
#endif
C
      RETURN
      END
