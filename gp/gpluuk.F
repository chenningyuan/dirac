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
C FILE    : gpluuk.F
C
C
C  /* Deck inpck */
      INTEGER FUNCTION INPCK (NI,INTS)
C
C     Packs small positive integers into a compound integer
C     L. Visscher
C
C     August 1997
C
#include "implicit.h"
#if defined (SYS_CRAY) || defined (SYS_T3D)
      PARAMETER (NBITS=64)
#else
      PARAMETER (NBITS=32)
#endif
C
      DIMENSION INTS(NI)
C
      NIBIT = NBITS / NI
C
      IWORD = 0
C
      NFACT = 1
      MFACT = 2**NIBIT
C
      DO I = 1, NI
         IWORD = IWORD + INTS(I) * NFACT
         NFACT = NFACT * MFACT
      ENDDO
C
      INPCK = IWORD
C
      RETURN
C
      END
C  /* Deck iunpck */
      SUBROUTINE IUNPCK (INT,NI,INTS)
C
C     Unpacks small integers from a compound integer
C     L. Visscher
C
C     August 1997
C
#include "implicit.h"
#if defined (SYS_CRAY) || defined (SYS_T3D)
      PARAMETER (NBITS=64)
#else
      PARAMETER (NBITS=32)
#endif
C
      DIMENSION INTS(NI)
C
      IWORD = INT
      NIBIT = NBITS / NI
C
      NFACT = 2**((NI-1)*NIBIT)
      MFACT = 2**NIBIT
C
      DO I = NI, 1, -1
         INTS(I) = IWORD / NFACT
         IWORD = MOD(IWORD,NFACT)
         NFACT = NFACT / MFACT
      ENDDO
C
      RETURN
C
      END
C  /* Deck indexx */
      SUBROUTINE INDEXX(N,ARR,INDX)
C***********************************************************************
C
C     This subroutine goes through array ARR and returns the ordering
C     index array INDX giving indices from smallest to largest
C
C***********************************************************************
      INTEGER N,INDX(N),M,NSTACK
      REAL*8 ARR(N)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,INDXT,IR,ITEMP,J,JSTACK,K,L,ISTACK(NSTACK)
      REAL*8 A
      DO J=1,N
         INDX(J)=J
      END DO
      JSTACK=0
      L=1
      IR=N
 1    IF(IR-L.GE.M)THEN
         K=(L+IR)/2
         ITEMP=INDX(K)
         INDX(K)=INDX(L+1)
         INDX(L+1)=ITEMP
         IF(ARR(INDX(L+1)).GT.ARR(INDX(IR)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ARR(INDX(L)).GT.ARR(INDX(IR)))THEN
            ITEMP=INDX(L)
            INDX(L)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ARR(INDX(L+1)).GT.ARR(INDX(L)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(L)
            INDX(L)=ITEMP
         ENDIF
         I=L+1
         J=IR
         INDXT=INDX(L)
         A=ARR(INDXT)
 3       CONTINUE
         I=I+1
         IF(ARR(INDX(I)).LT.A)GOTO 3
 4       CONTINUE
         J=J-1
         IF(ARR(INDX(J)).GT.A)GOTO 4
         IF(J.LT.I)GOTO 5
         ITEMP=INDX(I)
         INDX(I)=INDX(J)
         INDX(J)=ITEMP
         GOTO 3
 5       INDX(L)=INDX(J)
         INDX(J)=INDXT
         JSTACK=JSTACK+2
         IF(JSTACK.GT.NSTACK)
     $        CALL QUIT('NSTACK TOO SMALL IN INDEXX')
         IF(IR-I+1.GE.J-L)THEN
            ISTACK(JSTACK)=IR
            ISTACK(JSTACK-1)=I
            IR=J-1
         ELSE
            ISTACK(JSTACK)=J-1
            ISTACK(JSTACK-1)=L
            L=I
         ENDIF
      ELSE
         DO J=L+1,IR
            INDXT=INDX(J)
            A=ARR(INDXT)
            DO I=J-1,1,-1
               IF(ARR(INDX(I)).LE.A)GOTO 2
               INDX(I+1)=INDX(I)
            END DO
            I=0
 2          INDX(I+1)=INDXT
         END DO
         IF(JSTACK.EQ.0)RETURN
         IR=ISTACK(JSTACK)
         L=ISTACK(JSTACK-1)
         JSTACK=JSTACK-2
      ENDIF
      GOTO 1
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck IRECLEN*/
      FUNCTION IRECLEN (NREAL,NINTE,NLOGI)
C
C     Returns the record length for a record consisting of the specified
C     number of variables.
C
C     Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      INTEGER NREAL,NINTE,NLOGI,IRECLEN
C
#if defined (SYS_VAX)
      IRECL = 2 * NREAL + 1 * NINTE + 1 * NLOGI
#else
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IRECL = 8 * NREAL + 8 * NINTE + 8 * NLOGI
#else
#if defined (SYS_IRIX) || defined (SYS_SX) || defined (SYS_DEC)
C     SGI IRIX counts recl in *4 words
      IRECL = 2 * NREAL + 1 * NINTE + 1 * NLOGI
#else
#if defined (SYS_ALLIANT) || defined (SYS_CONVEX) || defined (SYS_AIX) || defined (SYS_PARAGON) || defined (SYS_DEC) || defined (SYS_HPUX) || defined (SYS_LINUX) || defined (SYS_DARWIN)
      IRECL = 8 * NREAL + 4 * NINTE + 4 * NLOGI
#else
#if defined (SYS_IBM) || defined (SYS_SUN)
      IRECL = 8 * NREAL + 4 * NINTE + 4 * NLOGI
#else
     ===== insert definition for specification of record length
#endif
#endif
#endif
#endif
#endif
C
      IRECLEN = IRECL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLANKR (N,A)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Circular shift of character array such that the blanks are to the
C     right side of the array.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N
      CHARACTER*1 A(N)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IS = 0
      DO I = N, 1, -1
         IF (A(I).NE.' ') GOTO 2
         IS = IS + 1
      ENDDO
C
    2 IF (IS.EQ.0) RETURN
C
      DO I = N-IS, 1, -1
            A(I+IS) = A(I)
      ENDDO
      DO I = 1, IS
         A(I) = ' '
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CONJUGA (N,X,INCX)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Take complex conjugate of vector X
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      COMPLEX*16 X(*)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
      INTEGER II
C
C---------------Executable code--------------------------------------
C
      II = 1
      DO I = 1, N
         X(II) = DCONJG(X(II))
C64B     X(II) =  CONJG(X(II))
         II = II + INCX
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine CpuUsed (TIMCPU)
C     Little wrapper for GETTIM
      Real*8 TIMCPU,TIMWAL
      CALL GETTIM (TIMCPU,TIMWAL)
      Return
      End
C
      Real*8 Function CpuDelta()
C
      Real*8 Sec1,Sec0
      Save Sec0
      Data Sec0 /0.d0/
C
      Call CPUused(Sec1)
      CPUdelta = Sec1 - Sec0
      Sec0 = Sec1
C
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine WallUsed (TIMWAL)
C     Another little wrapper for GETTIM
      Real*8 TIMCPU,TIMWAL
      CALL GETTIM (TIMCPU,TIMWAL)
      Return
      End
C
      Real*8 Function WallDelta()
C
      Real*8 Sec1,Sec0
      Save Sec0
      Data Sec0 /0.d0/
C
      Call WallUsed(Sec1)
      WallDelta = Sec1 - Sec0
      Sec0 = Sec1
C
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine DayTime (Datum,Tijd)
      Character Datum*10,Tijd*8

#if defined (SYS_CRAY) || defined (SYS_T3D)
      Call Date(Datum)
      Call Clock(Tijd)
#else
#if defined (SYS_HPUX) || defined (SYS_LINUX) || defined (SYS_DARWIN) \
   || defined (SYS_DEC) || defined (SYS_IRIX)  || defined (SYS_AIX) \
   || defined (SYS_SUN) || defined (SYS_SX)
      CHARACTER FDATE*24, String*24
C
      STRING = FDATE()
C     format from fdate(): 'Fri Sep 13 00:00:00 1986'
      DATUM = STRING(9:10) // STRING(4:8) // STRING(23:24) // ' '
C     format of DATUM: '13 Sep 86 '
      TIJD  = STRING(12:19)
C     format of TIJD: '00:00:00'
#else
      Datum = '          '
      Tijd  = '        '
#endif
#endif
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION INTOWP(N)
#if defined (SYS_CRAY) || defined (SYS_T3D)
      INTOWP = N
#else
      INTOWP=2*N
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRT_EV (LUPRI,N_IRPS,N_EV,REPNA,E_OFF,EV,
     &                   INDX,IRPS,IND)
C 
#include "implicit.h"
C
C---------------Description--------------------------------------------
C
C     Print list of eigenvalues from different symmetries, sorted on
C     increasing energy value.
C 
C     An offset is added to this value so that e.g. the constant
C     CCSD energy can be added to the Fockspace eigenvalues.
C
C     The routine is based on an older routine in GOSCIP but changed
C     to get rid of the common blocks.
C
C     Parameters :
C     - LUPRI   : Output unit
C     - N_IRPS  : Number of irreps
C     - N_EV    : Number of eigenvalues in each irrep
C     - REPNA   : Name of each irrep
C     - E_OFF   : Energy offset
C     - EV      : Eigenvalues
C     - INDX    : Work array of the same length as EV
C     - IRPS    : Work array of the same length as EV
C     - IND     : Work array of length N_IRPS
C
C---------------Routines called----------------------------------------
C
C     INDEXX
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher, august 2003.
C
C---------------Calling variables--------------------------------------
C
      DIMENSION N_EV(*),INDX(*)
      DIMENSION EV(*)
      DIMENSION IRPS(*),IND(*)
      CHARACTER*4 REPNA(*)

C
C---------------Common Blocks--------------------------------------
C
#include "codata.h"
C
C---------------Local variables--------------------------------------
C
      PARAMETER (DEGEN=1.D-9)
C
C---------------Executable code----------------
C     
C     Compute the total number of eigenvalues and fill
C     the array with the identication
C
      N_TOTAL = 0
      DO IRP = 1, N_IRPS
         CALL ICOPY(N_EV(IRP),IRP,0,IRPS(N_TOTAL+1),1)
         N_TOTAL = N_TOTAL + N_EV(IRP)
      ENDDO
C
      IF (N_TOTAL.EQ.0) RETURN
C     
C     Sort them from low to high
C
      CALL INDEXX(N_TOTAL,EV,INDX)
C
      WRITE (LUPRI,1005)
      ILEVEL=0
      EW=0.0D0
      EWF=EV(INDX(1))
      NEW=0
      AVER=0.0D0
      DO I=1,N_TOTAL
        E = EV(INDX(I))
        AVER=AVER+E
        IF (ABS(E-EW).GT.DEGEN) THEN
          IF (NEW.NE.0)
     &      WRITE(LUPRI,1010) ILEVEL,EW-EWF,EW+E_OFF,NEW
          ILEVEL=ILEVEL+1
          NEW=1
          EW=E
        ELSE
          NEW=NEW+1
          EW=((NEW-1)*EW+E)/DBLE(NEW)
        ENDIF
      ENDDO
      IF (NEW.NE.0)
     &WRITE(LUPRI,1010) ILEVEL,EW-EWF,EW+E_OFF,NEW
      WRITE(LUPRI,1020) AVER/DBLE(N_TOTAL)+E_OFF
C
      DO IRP1 = 1, N_IRPS, 8
      IRP2 = MIN0(N_IRPS,IRP1+7)
      WRITE(LUPRI,1025) (REPNA(I),I=IRP1,IRP2)
      ILEVEL=0
      EW=0.0D0
      EWF=EV(INDX(1))
      NEW=0
      DO JRP = IRP1, IRP2
         IND(JRP) = 0
      END DO
C
      DO 140 I=1,N_TOTAL
        E=EV(INDX(I))
        IF (ABS(E-EW).GT.DEGEN) THEN
          IF (NEW.NE.0)
     &      WRITE(LUPRI,1030) ILEVEL,
     &                     XTEV*(EW-EWF),XTKAYS*(EW-EWF),
     &                    (IND(JRP),JRP=IRP1,IRP2)
          DO JRP = IRP1, IRP2
            IND(JRP) = 0
          ENDDO
          ILEVEL = ILEVEL + 1
          IND(IRPS(INDX(I))) = 1
          NEW=1
          EW=E
        ELSE
          NEW=NEW+1
          IND(IRPS(INDX(I))) = IND(IRPS(INDX(I))) + 1
          EW=((NEW-1)*EW+E)/DBLE(NEW)
        ENDIF
 140    CONTINUE
      IF (NEW.NE.0) WRITE(LUPRI,1030) ILEVEL,
     &              XTEV*(EW-EWF),XTKAYS*(EW-EWF),
     &              (IND(JRP),JRP=IRP1,IRP2)
      ENDDO
 1000 FORMAT(F17.8)
 1008 FORMAT(17X,I8,2X,16(A1),F16.8)
 1005 FORMAT(//' Energy eigenvalues in atomic units'
     +            //' Level   Rel eigenvalue     ',
     +            'Abs eigenvalue  Degeneracy'/)
 1010 FORMAT(I5,F17.12,F20.12,' (',I4,' * )')
 1020 FORMAT(//' Total average: ',F20.10)
 1025 FORMAT(//' Relative eigenvalues in other units;'/' Symmetry',
     +       ' Classification in the Abelian subgroup',
     +       //' Level  eigenvalue (eV) ',
     +      ' Eigenvalue (cm-1)',2X,8A4/)
 1030 FORMAT(I5,F19.9,F18.6,8(I4))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTIME (TIMERT,FUNCT,NAME)
C
C     Initial routine written by M. Pernpointner
C     Documented and revised by L. Visscher
C
C     Input parameters 
C         - TIMERT : Type of timer (Used in MPE calls 
C                                   0:Do not log in MPE
C                                   1:Initialization
C                                   2:Wait,
C                                   3:Send/receive task,
C                                   4:Send/receive Data,
C                                   5:Perform computation)
C         - FUNCT  : TCLEAR=0 : Set the timer to zero
C                    TSTART=1 : Start adding time to this timer
C                    TSTOP=2  : Stop the timer
C                    TPRNTI=3 : Write contents of the timer to output
C                    TPRNTA=4 : Write contents of all timers to output
C                    TINIT=-1 : Initialize all timers
C         - NAME   : Name of the timer (max 30 characters)
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*30 NAME,T_REPORT
      PARAMETER ( D0 = 0.0D0 )
#include "timtbl.h"
#include "priunit.h"
#if defined (VAR_MPI)
#include "mpif.h"
#endif
#include "infpar.h"
      REAL*8 A,MY_T,CPU_AV,CPU_MI,CPU_MA
C
C     Check first for operations that involve all timers
C
      IF (FUNCT.EQ.TINIT) THEN
C     
C        Initialization of all timers
C
         DO I = 1, MXTIMF
          CPUTF(I)=D0
          CPUTFT(I)=D0
          CPUTT(I)=TIMERF_NAME
         ENDDO
         LAST_MPE_ACTION = 0
         LAST_MPE_STATE  = 0
         RETURN
      ELSEIF (FUNCT.EQ.TPRNTA) THEN
C     
C        Print out of all timers
C
        IF (MYTID.EQ.MPARID) WRITE(LUPRI,220) NAME
#if ! defined (VAR_MPI)
C       Serial code : no need to communicate with slaves
        DO I=1,MXTIMF
          IF(CPUTFT(I).NE.D0) WRITE(LUPRI,222) CPUTT(I),CPUTFT(I)
        ENDDO
#else
C
C       We may have the parallel version running in serial mode, 
C       handle this first (same code as serial)
C
        IF (NUMNOD.EQ.0) THEN
          DO I=1,MXTIMF
            IF(CPUTFT(I).NE.D0) WRITE(LUPRI,222) CPUTT(I),CPUTFT(I)
          ENDDO
        ELSE
C      
C        The parallel version runs in parallel. We need to communicate
C        to gather all data on the master.
C
C        Run over all master timers and collect data from slaves
C        NB : This also gets the inactive timers, but the loop structure
C        is easier this way
C
         IF (MYTID.EQ.MPARID)  WRITE(LUPRI,330)
C
         DO I = 1, MXTIMF
C          Tell the slaves to look for timer I
           T_REPORT = CPUTT(I)
           CALL MPI_BCAST (T_REPORT,30,MPI_CHARACTER,
     &                     MPARID,MPI_COMM_WORLD,IERR)
C          Find the timer (we do this code on master also for simplicity)
           TIMERN = -1
           DO J = 1, MXTIMF
              IF (T_REPORT.EQ.CPUTT(J)) TIMERN = J
           ENDDO
           IF (TIMERN.EQ.-1) THEN
C             Timer not used, give zero values
              MY_T = D0
           ELSE
              MY_T = CPUTFT(TIMERN)
           ENDIF
C          NB : The following computations are only correct on the master,
C          but since we report only there there is no need for an ALLREDUCE call.
C          Compute the average time spent by the slaves
           IF (MYTID.EQ.MPARID) THEN
              A = D0
           ELSE
              A = MY_T
           ENDIF
           CALL MPI_REDUCE (A,CPU_AV,1,MPI_DOUBLE_PRECISION,
     &                      MPI_SUM,MPARID,MPI_COMM_WORLD,IERR)
           CPU_AV = (CPU_AV)/DBLE(NUMNOD)
C          Compute the minimum time
           IF (MYTID.EQ.MPARID) A = 1.D8
           CALL MPI_REDUCE (A,CPU_MI,1,MPI_DOUBLE_PRECISION,
     &                      MPI_MIN,MPARID,MPI_COMM_WORLD,IERR)
C          Compute the maximum time
           IF (MYTID.EQ.MPARID) A = -1.D8
           CALL MPI_REDUCE (A,CPU_MA,1,MPI_DOUBLE_PRECISION,
     &                      MPI_MAX,MPARID,MPI_COMM_WORLD,IERR)
C
C          Write the result to the output
C
           IF (MYTID.EQ.MPARID.AND.(CPU_AV.NE.D0.OR.MY_T.NE.D0))
     &        WRITE(LUPRI,340) T_REPORT,MY_T,CPU_MI,CPU_MA,CPU_AV
C
         ENDDO
C
C        The slaves may also have timers that master does not have, let
C        the primus inter pares report.
C
         DO I = 1, MXTIMF
           SPARTACUS = MPARID+1
           T_REPORT = CPUTT(I)
           CALL MPI_BCAST (T_REPORT,30,MPI_CHARACTER,
     &                     SPARTACUS,MPI_COMM_WORLD,IERR)
           TIMERN = -1
           DO J = 1, MXTIMF
              IF (T_REPORT.EQ.CPUTT(J)) TIMERN = J
           ENDDO
           IF (TIMERN.EQ.-1) THEN
C             Timer not used, give zero values
              MY_T = D0
           ELSE
              MY_T = CPUTFT(TIMERN)
           ENDIF
           IF (MYTID.EQ.MPARID) THEN
              A = D0
           ELSE
              A = MY_T
           ENDIF
C          Compute the average time spent by the slaves
           CALL MPI_REDUCE (A,CPU_AV,1,MPI_DOUBLE_PRECISION,
     &                      MPI_SUM,MPARID,MPI_COMM_WORLD,IERR)
           CPU_AV = (CPU_AV)/DBLE(NUMNOD)
C          Compute the minimum time
           IF (MYTID.EQ.MPARID) A = 1.D8
           CALL MPI_REDUCE (A,CPU_MI,1,MPI_DOUBLE_PRECISION,
     &                      MPI_MIN,MPARID,MPI_COMM_WORLD,IERR)
C          Compute the maximum time
           IF (MYTID.EQ.MPARID) A = -1.D8
           CALL MPI_REDUCE (A,CPU_MA,1,MPI_DOUBLE_PRECISION,
     &                      MPI_MAX,MPARID,MPI_COMM_WORLD,IERR)
C
C          Write the result to the output
C          The master has to have TIMERN -1, otherwise this timer
C          has already been reported in the previous loop.
C
           IF (MYTID.EQ.MPARID.AND.TIMERN.EQ.-1.AND.CPU_AV.NE.D0)
     &        WRITE(LUPRI,340) T_REPORT,MY_T,CPU_MI,CPU_MA,CPU_AV
        ENDDO
      ENDIF ! end of 1 / multiple processor clause
      IF (MYTID.EQ.MPARID) WRITE(LUPRI,350)
#endif
      RETURN
      ENDIF ! end of all timers IF-THEN-ELSE clause
C
C     We come here only if we work on a single timer, find out which one
C
      TIMERN = -1
      DO I = 1, MXTIMF
         IF (NAME.EQ.CPUTT(I)) TIMERN = I
      ENDDO
C 
C     We may have a new timer, look for free slots
C
      IF (TIMERN.EQ.-1) THEN
         TIMERF = 0
         DO I = MXTIMF,1,-1
            IF (CPUTT(I).EQ.TIMERF_NAME) TIMERF = I
         ENDDO
C
         IF (TIMERF.NE.0) THEN
C           Use the first free slot and initialize
            CPUTT(TIMERF)=NAME
            CPUTFT(TIMERF)=D0
            TIMERN = TIMERF
         ELSE 
C           Too bad, we ran out of slots
            WRITE(LUPRI,130) MYTID,NAME
            RETURN
         ENDIF
      ENDIF
C
C     Perform the desired action
C
      IF(FUNCT.EQ.TCLEAR) THEN
          CPUTF(TIMERN)=D0
          CPUTFT(TIMERN)=D0
          CPUTT(TIMERN)=TIMERF_NAME
      ELSE IF (FUNCT.EQ.TSTART) THEN
        CALL CPUUSED(A)
        CPUTF(TIMERN)=A
        CPUTT(TIMERN)=NAME
#if defined (MPE)
        IF (TIMERT.NE.0) THEN
C          We may need to end the previous state
           IF (TIMERT.NE.LAST_MPE_STATE.AND.LAST_MPE_ACTION.EQ.1) THEN
              I_MPE = 2 * LAST_MPE_STATE
              CALL MPE_LOG_EVENT(I_MPE,0,NAME)
           ENDIF
C          Announce the new state if necessary
           IF (TIMERT.NE.LAST_MPE_STATE.OR.LAST_MPE_ACTION.NE.1) THEN
              I_MPE = 2 * TIMERT - 1
              CALL MPE_LOG_EVENT(I_MPE,0,NAME)
              LAST_MPE_STATE  = TIMERT
              LAST_MPE_ACTION = 1    
           ENDIF
        ENDIF
#endif
      ELSE IF (FUNCT.EQ.TSTOP) THEN
        CALL CPUUSED(A)
        CPUTF(TIMERN)=A-CPUTF(TIMERN)
        CPUTFT(TIMERN)=CPUTFT(TIMERN)+CPUTF(TIMERN)
#if defined (MPE)
        IF (TIMERT.NE.0) THEN
C          The state to be ended may be interrupted and ended by a start
C          Only do something if everything matches
           IF (TIMERT.EQ.LAST_MPE_STATE.AND.LAST_MPE_ACTION.EQ.1) THEN
              I_MPE = 2 * TIMERT
              CALL MPE_LOG_EVENT(I_MPE,0,NAME)
              LAST_MPE_STATE  = TIMERT
              LAST_MPE_ACTION = 2    
           ENDIF
        ENDIF
#endif
      ELSE IF (FUNCT.EQ.TPRNTI) THEN
        IF (MYTID.EQ.MPARID) THEN
           IF(CPUTF(TIMERN).NE.D0) THEN
             WRITE(LUPRI,222) CPUTT(TIMERN),CPUTFT(TIMERN)
           ELSE
             WRITE(LUPRI,224) TIMERN
           ENDIF
        ENDIF
      ELSE
        WRITE(LUPRI,135)
      ENDIF
C
 130  FORMAT(/' **WARNING** No free timer slots for proces',I4/,
     &        ' ** no action taken')
 135  FORMAT(/' **WARNING** Unknown timer function'/,
     &        ' ** no action taken')
 220  FORMAT(/'------ Timing report (in CPU seconds) of module ',A/)
 222  FORMAT(' Time in ',A30,T33,F12.3,' seconds')
 223  FORMAT(T37,'-------------',/,T25,'Sum: ',F12.3,' seconds',
     &       /,T37,'-------------')
 224  FORMAT(' Timer ',I3,' not used for measurements')
 330  FORMAT(//,T37,'Master',T46,'--------- Slaves ------------'//,
     &          T2,'Timer',T37,'Master',T46,'Minimum',T56,'Maximum',
     &          T68,'Average'/)
 340  FORMAT(T2,A30,T33,F10.1,T43,F10.1,T53,F10.1,T63,F10.1)
 350  FORMAT(/'------ End of timing report ------'/)
     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTIMRD (ELTIM,NAME)
C
C  purpose:  read out a specific timer previously set by XTIME
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*30 NAME
      REAL*8 ELTIM
      PARAMETER ( D0 = 0.0D0 )
#include "timtbl.h"
#include "priunit.h"
      TIMERN = -1
      DO I = 1, MXTIMF
         IF (NAME.EQ.CPUTT(I)) TIMERN = I
      ENDDO
      IF(TIMERN.LT.1.OR.TIMERN.GT.MXTIMF) THEN
        ELTIM=D0
      ELSE
        ELTIM=CPUTFT(TIMERN)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_INFO (NBUF,NSIZE)
C     
C     Report the maximum number of buffer and record size that 
C     will be used in the sort
C
C     Calling parameters
C     
C     NBUF   : Maximum number of buffers (hardwired in the common block)
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
C
      NBUF = N_BUFFER
      NSIZE = NGBFSZ
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_INIT (I_FILE,IROWS,NROWS)
C     
C     Initialize a buffer file to transpose a large matrix that is generated in parts
C     Currently used in MOLTRA for the halftransformed two-electron
C     integrals but it may be useful at other places as well.
C
C     Calling parameters
C     
C     I_FILE   : Bucket sort indentifier (we may have more than one active)
C     IROWS    : Pointer from rows to buffer and position therein
C     NROWS    : Number of rows
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
      DIMENSION IROWS(3,NROWS)
      CHARACTER*10 NAME*5,FNODE,LFNAM
C
C     The record size needs to be NGBFSZ real + (2*NGBFSZ + 2) integer
C     because we also need to store the labels
C
      IRECL = IRECLEN(NGBFSZ,2*NGBFSZ+2,0)
C
C     The opening of a file is a bit tricky in FORTRAN, pray that the unit
C     is not yet in use....
C
      WRITE (NAME,'(A4,I1)') 'BUFF',I_FILE 
      FNODE = LFNAM(NAME)
      LGFIL = I_BF_UNIT+I_FILE
      OPEN (UNIT=LGFIL,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
C
C     Initialize the data that is kept in the common block
C
      I_BF_REC(I_FILE) = 0
      BF_NAME(I_FILE) = FNODE
      CALL IZERO(LGREC(1,I_FILE),N_BUFFER)
      CALL IZERO(LGBUF(1,I_FILE),N_BUFFER)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_WRITE (I_FILE,RGBUF,IGBUF,IROWS,
     &                              A,NROWS,ICOL,TRESHOLD)                              
C     
C     Scatter row of matrix A over the buffers, flush them if necessary
C     
C     Calling parameters
C     
C     LGFIL    : File unit
C     RGBUF    : Real data in buffers
C     IGBUF    : Integer data data in buffers
C     IROWS    : Pointer from the rows to buffer
C     A        : Column of matrix that is to be buffered away
C     NROWS    : Number of rows in A
C     ICOL     : Column index of A
C     TRESHOLD : Treshold for writing the matrix element
C
C     Written by Luuk Visscher, August 2004
C     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
      REAL*8 RGBUF(NGBFSZ,*), A(NROWS)
      DIMENSION IGBUF(2,NGBFSZ,*)
      DIMENSION IROWS(3,NROWS)
C
C     Get the information stored in the common block 
C
      LGFIL    = I_BF_UNIT+I_FILE
      IREC     = I_BF_REC(I_FILE) 
C
C     Loop through the row and check if the data need be stored
C
      DO I = 1, NROWS
         IBUF = IROWS(1,I)
         IF (IBUF.GT.0.AND.ABS(A(I)).GT.TRESHOLD) THEN
            IROW = IROWS(2,I)
            N = LGBUF(IBUF,I_FILE) + 1
            RGBUF(N,IBUF) = A(I)
            IGBUF(1,N,IBUF) = ICOL
            IGBUF(2,N,IBUF) = IROW
            IF (N.GE.NGBFSZ) THEN 
               IREC = IREC + 1
               WRITE (LGFIL,REC=IREC) LGREC(IBUF,I_FILE),N,
     &               (RGBUF(JBUF,IBUF),JBUF=1,NGBFSZ),
     &               ((IGBUF(IX,JBUF,IBUF),IX=1,2),JBUF=1,NGBFSZ)
               LGREC(IBUF,I_FILE) = IREC
               N = 0
            ENDIF
            LGBUF(IBUF,I_FILE) = N
        ENDIF
      ENDDO
C
      I_BF_REC(I_FILE) = IREC
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_FLUSH (I_FILE,RGBUF,IGBUF,IROWS,NROWS)                              
C     
C     Flush the buffers
C     
C     Calling parameters
C     
C     LGFIL    : File unit
C     RBUF     : Real data in buffers
C     IBUF     : Integer data data in buffers
C     IROWS    : Pointer from rows to buffer and position therein
C     NROWS    : Number of rows
C     
C     Written by Luuk Visscher, August 2004
C     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
      REAL*8 RGBUF(NGBFSZ,*)
      DIMENSION IGBUF(2,NGBFSZ,*)
      DIMENSION IROWS(3,NROWS)
C
C     Get the information stored in the common block 
C
      LGFIL    = I_BF_UNIT+I_FILE
      IREC     = I_BF_REC(I_FILE) 
C
C     Loop through the rows and check if the data need be flushed
C
      DO I = 1, NROWS
         IBUF = IROWS(1,I)
         N = LGBUF(IBUF,I_FILE)
         IF (IBUF.GT.0.AND.N.GT.0) THEN 
            IREC = IREC + 1
            WRITE (LGFIL,REC=IREC) LGREC(IBUF,I_FILE),N,
     &            (RGBUF(JBUF,IBUF),JBUF=1,NGBFSZ),
     &            ((IGBUF(IX,JBUF,IBUF),IX=1,2),JBUF=1,NGBFSZ)
            LGREC(IBUF,I_FILE) = IREC
            N = 0
         ENDIF
         LGBUF(IBUF,I_FILE) = N
      ENDDO
C
      I_BF_REC(I_FILE) = IREC
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_READ (I_FILE,IBUF,A,NROW,ICOL)                              
C
C     Read back the sorted integrals and put them at the right position
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
      REAL*8 RGBUF(NGBFSZ),A(NROW)
      DIMENSION IGBUF(2,NGBFSZ)
C
      LGFIL    = I_BF_UNIT+I_FILE
      IREC     = LGREC(IBUF,I_FILE)
C
      CALL DZERO (A,NROW)
C
    1 IF (IREC.EQ.0) RETURN
C
C     Read a batch and the link to the previous record.
C
      READ (LGFIL,REC=IREC) IRECN,NABFSZ,RGBUF,IGBUF
      IREC = IRECN
C
C     Pick out the integrals that we need.
C
      DO JBUF = 1, NABFSZ
         IF (IGBUF(2,JBUF).EQ.ICOL) THEN
            A(IGBUF(1,JBUF)) = RGBUF(JBUF)
         ENDIF
      ENDDO
C
C     We read in reverse order, if we have the first block we're done
C
      GOTO 1
C
      RETURN
      END
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_DIST (I_FILE_1,IROWS,NROWS)    
C     
C     Send the buffers to the correct place
C     
C     Calling parameters
C     
C     LGFIL    : File unit
C     IROWS    : Pointer from rows to buffer and position therein
C     NROWS    : Number of rows
C     
C     Written by Luuk Visscher, August 2004
C     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
      INCLUDE 'mpif.h'
      DIMENSION ISTAT(MPI_STATUS_SIZE),ISTATS(MPI_STATUS_SIZE,3)
#include "infpar.h"
      REAL*8 RGBUF(NGBFSZ),RHBUF(NGBFSZ)
      DIMENSION IGBUF(2,NGBFSZ),IHBUF(2,NGBFSZ)
      DIMENSION IROWS(3,NROWS)
      INTEGER INFO(2),JNFO(2),IREQ(4)
      LOGICAL NEWMAIL, READY
C
C     Open a new file for the redistributed integrals
C
      I_FILE_2 = I_FILE_1 + 1
      CALL BUFFER_FILE_INIT (I_FILE_2,IROWS,NROWS)
C     
C     Get the file units of the old and new file
C
      LGFIL_1   = I_BF_UNIT + I_FILE_1
      LGFIL_2   = I_BF_UNIT + I_FILE_2
C
C     Get the start record for the new file
C
      IREC_2    = I_BF_REC(I_FILE_2) 
C
C     Put all local integrals that should stay on this node on the new file
C
      DO I = 1, NROWS
         IBUF = IROWS(1,I)
         IPOS = IROWS(2,I)
         IF (IPOS.EQ.1) THEN
            IBUF_PREVIOUS = IBUF
            JNOD = IROWS(3,I)
            IF (JNOD.EQ.MYTID) THEN
C
C              Read the integrals
C
               IREC_1 = LGREC(IBUF,I_FILE_1)
    1          IF (IREC_1.EQ.0) GOTO 10
               READ (LGFIL_1,REC=IREC_1) IRECN,N,RGBUF,IGBUF
C
C              Write to the new file
C
               IREC_2 = IREC_2 + 1
               WRITE (LGFIL_2,REC=IREC_2) LGREC(IBUF,I_FILE_2),N,
     &                                    RGBUF,IGBUF
               LGREC(IBUF,I_FILE_2) = IREC_2
               IREC_1 = IRECN
               GOTO 1
   10          CONTINUE
            ENDIF
         ENDIF
      ENDDO
C
C     Send all other local integrals to their destination and receive them from the others
C     
      DO IR = 1, 4
         IREQ(IR) = MPI_REQUEST_NULL
      ENDDO
C
      DO I = 1, NROWS 
         IBUF = IROWS(1,I)
         IPOS = IROWS(2,I)
         IF (IPOS.EQ.1) THEN
            IBUF_PREVIOUS = IBUF
            JNOD = IROWS(3,I)
            IF (JNOD.NE.MYTID) THEN
C           
C              Read the integrals
C              
               IREC_1 = LGREC(IBUF,I_FILE_1)
   11          IF (IREC_1.EQ.0) GOTO 20
C
C              Check whether the previous send has completed
C
   12          CALL MPI_TESTALL(3, IREQ, READY, ISTATS, IERR)
               IF (.NOT.READY) THEN
C
C                 We use the waiting time to receive integrals from the others
C
                  CALL MPI_IPROBE(MPI_ANY_SOURCE,95,
     &                            MPI_COMM_WORLD,NEWMAIL,ISTAT,IERR)
                  IF (NEWMAIL) THEN 
                     INOD = ISTAT(MPI_SOURCE)
                     CALL MPI_RECV(JNFO,2,MPI_INTEGER,INOD,95,
     &                             MPI_COMM_WORLD,ISTAT,IERR)
                     M    = JNFO(1)
                     JBUF = JNFO(2)
                     CALL MPI_RECV(IHBUF,2*M,MPI_INTEGER,INOD,96,
     &                             MPI_COMM_WORLD,ISTAT,IERR)
                     CALL MPI_RECV(RHBUF,M,MPI_DOUBLE_PRECISION,INOD,97,
     &                             MPI_COMM_WORLD,ISTAT,IERR)
                     IREC_2 = IREC_2 + 1
                     WRITE (LGFIL_2,REC=IREC_2) LGREC(JBUF,I_FILE_2),M,
     &                                          RHBUF,IHBUF
                     LGREC(JBUF,I_FILE_2) = IREC_2
                  ENDIF
                  GOTO 12
               ENDIF
C
C              Our read buffer is free, read integrals
C
               READ (LGFIL_1,REC=IREC_1) IRECN,N,RGBUF,IGBUF
C              
C              Send them to their destination
C
               INFO(1) = N
               INFO(2) = IBUF
               CALL MPI_ISEND(INFO,2,MPI_INTEGER,JNOD,95,
     &                        MPI_COMM_WORLD,IREQ(1),IERR)
               CALL MPI_ISEND(IGBUF,2*N,MPI_INTEGER,JNOD,96,
     &                        MPI_COMM_WORLD,IREQ(2),IERR)
               CALL MPI_ISEND(RGBUF,N,MPI_DOUBLE_PRECISION,JNOD,97,
     &                        MPI_COMM_WORLD,IREQ(3),IERR)
C
               IREC_1 = IRECN
               GOTO 11
   20          CONTINUE
            ENDIF
         ENDIF
      ENDDO
C
C     We have sent all our integrals, tell the others
C
      DO INOD = 0, NODES
         IF (INOD.NE.MYTID) CALL MPI_ISEND(1,1,MPI_INTEGER,INOD,98,
     &                        MPI_COMM_WORLD,IREQ(4),IERR)
      ENDDO
C
C     Go into loop to wait for integrals and completion signals
C
      IDONE = 0
  100 CONTINUE
      CALL MPI_IPROBE(MPI_ANY_SOURCE,95,
     &                MPI_COMM_WORLD,NEWMAIL,ISTAT,IERR)
      IF (NEWMAIL) THEN 
         INOD = ISTAT(MPI_SOURCE)
         CALL MPI_RECV(JNFO,2,MPI_INTEGER,INOD,95,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
         M    = JNFO(1)
         JBUF = JNFO(2)
         CALL MPI_RECV(IHBUF,2*M,MPI_INTEGER,INOD,96,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
         CALL MPI_RECV(RHBUF,M,MPI_DOUBLE_PRECISION,INOD,97,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
         IREC_2 = IREC_2 + 1
         WRITE (LGFIL_2,REC=IREC_2) LGREC(JBUF,I_FILE_2),M,
     &                              RHBUF,IHBUF
         LGREC(JBUF,I_FILE_2) = IREC_2
      ENDIF
      CALL MPI_IPROBE(MPI_ANY_SOURCE,98,
     &                MPI_COMM_WORLD,NEWMAIL,ISTAT,IERR)
      IF (NEWMAIL) THEN 
         IDONE = IDONE + 1
         INOD = ISTAT(MPI_SOURCE)
         CALL MPI_RECV(M,1,MPI_INTEGER,INOD,98,
     &                 MPI_COMM_WORLD,ISTAT,IERR)
      ENDIF
      IF (IDONE.LT.NODES) GOTO 100
C
      I_BF_REC(I_FILE_2) = IREC_2
C
      CLOSE (LGFIL_1,STATUS='DELETE')
C
      RETURN
      END
#endif
