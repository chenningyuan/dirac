C
C...   Copyright (c) 1998 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of 
C...   "Dirac, a relativistic ab initio electronic structure program, 
C...   Release 3.1 (1998), written by T. Saue, T. Enevoldsen, T. Helgaker,
C...   H. J. Aa. Jensen, J. Laerdahl, K. Ruud, J. Thyssen, and L. Visscher"
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
C...      dirac-admin@dirac.chem.ou.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.ou.dk/Dirac
C
C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : ccgenp.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTIMRD (TIMERN,ELTIM)
C
C  purpose:  read out a specific timer previously set by XTIME
C

      IMPLICIT INTEGER (A-Z)
      REAL*8 ELTIM

#include "timtbl.inc"
#include "param.inc"
#include "files.inc"
      
      IF(TIMERN.LT.1.OR.TIMERN.GT.MXTIMF) THEN
        WRITE(IW,10)
        ELTIM=D0
      ELSE
        ELTIM=CPUTFT(TIMERN)
      ENDIF

 10   FORMAT(/' **WARNING** Timer # out of range'/,
     &        ' ** zero returned')

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTIME (TIMERN,FUNCT,NAME)

C
C     Routine written by M. Pernpointner
C
C     Input parameters 
C         - TIMERN : ID of the timer (used with FUNCT 1 or 2)
C         - FUNCT  : TCLEAR=0 : Set ALL timers to zero
C                    TPRNTI=3 : Write contents of ALL timers to output (serial version, not used)
C                    TPRNTA=4 : Write contents of ALL timers to output (parallel version, always used)
C                    TSTART=1 : Add time until stop to the timer
C                    TSTOP=2 :  Stop the timer
C         - NAME   : Name of the timer (max 30 characters)
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*30 NAME

#include "timtbl.inc"
#include "param.inc"
#include "files.inc"
#include "ccpar.inc"
#if defined (VAR_MPI)
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
      PARAMETER(TITAG=12)
#endif

      REAL*8 A
      REAL*8 SUMT

C ** check requested function

      IF(TIMERN.LT.1.OR.TIMERN.GT.MXTIMF) THEN
        WRITE(IW,130) MXTIMF
        RETURN
      ENDIF

      IF(FUNCT.EQ.TCLEAR) THEN
        DO I=1,MXTIMF
          CPUTF(I)=D0
          CPUTFT(I)=D0
          CPUTAV(I)=D0
          CPUTMI(I)=D0
          CPUTMA(I)=D0
          CPUTT(I)='                    '
        ENDDO
      ELSE IF (FUNCT.EQ.TSTART) THEN
        CALL CPUUSED(A)
        CPUTF(TIMERN)=A
        CPUTT(TIMERN)=NAME
      ELSE IF (FUNCT.EQ.TSTOP) THEN
        CALL CPUUSED(A)
        CPUTF(TIMERN)=A-CPUTF(TIMERN)
        CPUTFT(TIMERN)=CPUTFT(TIMERN)+CPUTF(TIMERN)
      ELSE IF (FUNCT.EQ.TPRNTI) THEN
        IF(CPUTF(TIMERN).NE.D0) THEN
          WRITE(IW,222) CPUTT(TIMERN),CPUTF(TIMERN)
        ELSE
          WRITE(IW,224) TIMERN
        ENDIF
      ELSE IF (FUNCT.EQ.TPRNTA) THEN
#if defined (VAR_MPI)
        IF(NMPROC.EQ.1) THEN
          WRITE(IW,220)
          DO I=1,MXTIMF
            IF(CPUTFT(I).NE.D0) WRITE(IW,222) CPUTT(I),CPUTFT(I)
          ENDDO
        ELSE
          IF(MYPROC.EQ.MASTER) THEN
            DO I=1,MXTIMF
              CPUTMI(I)=CPUTFT(I)
              CPUTMA(I)=CPUTFT(I)
              CPUTAV(I)=CPUTFT(I)
            ENDDO
            WRITE(IW,300)
            WRITE(IW,310) MASTER
            DO I=1,MXTIMF
              IF(CPUTFT(I).NE.D0) WRITE(IW,222) CPUTT(I),CPUTFT(I)
            ENDDO
C
            DO ND=1,NMPROC-1
              CALL MPI_RECV(CPUTFT,MXTIMF,MPI_DOUBLE_PRECISION,
     &           ND,TITAG,MPI_COMM_WORLD,ISTAT,IERR)
              IF(IERR.NE.MPI_SUCCESS) CALL 
     &          QUIT('MPI ERROR IN XTIME')
Cluuk         Uncomment for detailed output, but you do not want to have this
Cluuk         on a machine with many active nodes....
cluuk         WRITE(IW,310) ND
              DO I=1,MXTIMF
                IF(CPUTFT(I).LT.CPUTMI(I)) CPUTMI(I)=CPUTFT(I)
                IF(CPUTFT(I).GT.CPUTMA(I)) CPUTMA(I)=CPUTFT(I)
                CPUTAV(I)=CPUTAV(I)+CPUTFT(I)
Cluuk           IF(CPUTFT(I).NE.D0) WRITE(IW,222) CPUTT(I),CPUTFT(I)
              ENDDO
            ENDDO

            WRITE(IW,330)
            DO I=1,MXTIMF
              IF(CPUTAV(I).NE.D0) WRITE(IW,340) CPUTT(I),CPUTMI(I),
     &          CPUTMA(I),CPUTAV(I)/DBLE(NMPROC)
            ENDDO
          ELSE
            CALL MPI_SEND(CPUTFT,MXTIMF,MPI_DOUBLE_PRECISION,
     &          MASTER,TITAG,MPI_COMM_WORLD,IERR)
            IF(IERR.NE.0) CALL QUIT('MPI ERROR in XTIME')
            WRITE(IW,310) MYPROC
            DO I=1,MXTIMF
               IF(CPUTFT(I).NE.D0) WRITE(IW,222) CPUTT(I),CPUTFT(I)
            ENDDO
          ENDIF
        ENDIF
        WRITE(IW,350)
#else
        WRITE(IW,220)
        DO I=1,MXTIMF
          IF(CPUTFT(I).NE.D0) WRITE(IW,222) CPUTT(I),CPUTFT(I)
        ENDDO
#endif
      ELSE
        WRITE(IW,135)
      ENDIF
C
 130  FORMAT(/' **WARNING** Timer # out of range'/,
     &        ' ** no action taken')
 135  FORMAT(/' **WARNING** Unknown timer function'/,
     &        ' ** no action taken')
 220  FORMAT(/'########################################'/,
     &        '###'/'### Timing results:'/'###'/,
     &        '########################################'/)
 222  FORMAT(' Time in ',A30,T33,F12.3,' seconds')
 223  FORMAT(T37,'-------------',/,T25,'Sum: ',F12.3,' seconds',
     &       /,T37,'-------------')
 224  FORMAT(' Timer ',I3,' not used for measurements')
 300  FORMAT(/'########################################'/,
     &        '###'/'### Timing results for all nodes'/'###'/,
     &        '########################################'/)
 310  FORMAT(/,T10,'Results node ',I3,/,T10,'----------------'/)
 330  FORMAT(//,T30,'Summary with min/max/average'/,
     &          T30,'----------------------------'//,
     &          T2,'Timer',T36,'Minimum',T48,'Maximum',
     &          T59,'Average'/)
 340  FORMAT(T2,A30,T33,F12.3,T43,F12.3,T53,F12.3)
 350  FORMAT(/'########################################'/,
     &        '###'/'### End of timing report'/'###'/,
     &        '########################################'/)
     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION SUMT(A,N)
      REAL*8 SUMT,A(*)
      INTEGER N

#include "param.inc"

      SUMT=D0
      DO I=1,N
        SUMT=SUMT+A(I)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FLIN(A,IDIM,IN,IM,DET)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      DIMENSION A(IDIM,*)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      N=IN
      NR=IM
      JMAX=N+NR
      SIGN=ONE
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
      TEMP=ZERO
      DO 41 I=M,N
      IF(M.GT.1)A(I,M)=A(I,M)-DOTT(A(I,1),IDIM,A(1,M),1,M-1)
      AVAL=A(I,M)
      IF(ABS(AVAL).LE.TEMP)GOTO 41
      TEMP=ABS(AVAL)
      IMAX=I
 41   CONTINUE
      IF(TEMP.LE.ZERO)GOTO 999
      IF(IMAX.EQ.M)GOTO 45
      SIGN=-SIGN
      DO 44 J=1,JMAX
      STOR=A(M,J)
      A(M,J)=A(IMAX,J)
      A(IMAX,J)=STOR
 44   CONTINUE
 45   CONTINUE
      JJ=M+1
      IF(JJ.GT.JMAX)GOTO 49
      IF(M.GT.1)GOTO 47
      DO 46 J=JJ,JMAX
      A(1,J)=A(1,J)/A(1,1)
 46   CONTINUE
      D=A(1,1)
      GOTO 49
 47   CONTINUE
      DO 48 J=JJ,JMAX
      A(M,J)=(A(M,J)-DOTT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48   CONTINUE
      D=D*A(M,M)
 49   CONTINUE
      IF(NR.EQ.0) RETURN
      DO 59 I=1,NR
      NPI=N+I
      DO 58 K=2,N
      J=N+1-K
      A(J,NPI)=A(J,NPI)-DOTT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58   CONTINUE
 59   CONTINUE
C***  IF(ABS(D).GE.1.0D+36) D=1.0D+36
C***  IF(ABS(D).LE.1.0D-36) D=1.0D-36
      DET=D*SIGN
      RETURN
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
 999  DET=ZERO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      FUNCTION DOTT(A,NA,B,NB,N)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION A(*),B(*)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
C
      IAPT=1
      IBPT=1
      D   =ZERO
      DO 10 I=1,N
        D   = D+A(IAPT)*B(IBPT)
        IAPT= IAPT+NA
        IBPT= IBPT+NB
 10   CONTINUE
      DOTT=D
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE QUITC
      CALL QUIT ("Error in one of the C I/O c-routines")
      RETURN
      END
