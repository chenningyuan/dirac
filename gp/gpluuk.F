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

