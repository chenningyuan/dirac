C
C VERSION : $Revision$
C DATE    : $Date$
C FILE    : cigosc.F
C SOURCE  : $Source$
C
      SUBROUTINE PAMGOS
C     ===============================================================
C     SMALL CI PROGRAM TO PERFORM COMPLETE OPEN SHELL CI CALCULATIONS
C     METHOD DESCRIBED IN JCP 96, 2910, (1992).
C     FIRST VERSION WRITTEN NOVEMBER 1990 BY O. VISSER
C     THEORETICAL CHEMISTRY DEPARTMENT,
C     UNIVERSITY OF GRONINGEN (RUG), GRONINGEN,THE NETHERLANDS
C     ==========================================================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      CHARACTER*8 DATEX*10,TIMEX
C
      EXTERNAL INIT
      CALL DAYTIME(DATEX,TIMEX)
      WRITE(6,1000) DATEX,TIMEX
      CALL CPUUSED(SEC)
      CALL READGOS
      CALL RDTWO
      DO 10 IRP=1,NREP
        IRRP = IRP
        IF (MOD(NELEC,2).EQ.0) IRRP = IRP + NREP
        CALL CPUUSED(SEC0)
        CALL GENDET(IRP,IRRP)
        IF (NDET.EQ.0) GOTO 10
        CALL CPUUSED(SEC1)
        CALL MAKEH
        CALL CPUUSED(SEC2)
        CALL DIAGH(IRP)
        CALL CPUUSED(SEC3)
        CALL WRVEC(IRP,IRRP)
        CALL CPUUSED(SEC4)
        WRITE (6,1010) SEC1-SEC0,SEC2-SEC1,SEC3-SEC2,SEC4-SEC3
        CALL PRTOUT(IRP,IRRP)
 10     CONTINUE
      INTREP = MOD(NELEC,2)
      CALL PRTOUT(0,INTREP)
      CALL CPUUSED(SEC5)
      WRITE (6,1020) SEC5-SEC
 1000 FORMAT(///1X,'GOSCIP (Version 1.8)'//
     &' Today is :',T15,A10/' The time is :',T17,A8//)
 1010 FORMAT(/' CPU time for this representation '//
     &' Generation of determinants :',T30,F12.4/
     &' Building the CI matrix :',T30,F12.4/
     &' Diagonalization :',T30,F12.4/
     &' Writing CI vectors :',T30,F12.4)
 1020 FORMAT (//'Total CPU time :',T30,F12.4//8X,
     &          '(NORMAL END OF PROGRAM)'///)
C
C VERSION 0.1: SEEMS TO BE WORKING FOR P-SHELLS OF ATOM
C              POPULATIONS ARE GIVEN ONLY IF $POPAN IS PRESENT
C                                                   OV, 19-11-1990
C VERSION 0.2: COMPRESSED LIST OF EIGENVALUES GIVEN AT THE END OF THE OUTPUT
C                                                   OV, 17-12-1990
C VERSION 0.3: TAKE CARE OF BREIT INTERACTION IN CALCULATION OF ENERGY
C                                                   LV, 21-12-1990
C VERSION 0.4: ADAPTED TO  CRAY
C                                                   LV, 23-1-1991
C VERSION 0.5: GIVE POPULATION ANALYSIS IN DIFFERENT FORM IF DESIRED
C                                                   OV, 7-2-1991
C VERSION 0.6: BUG FIX IN MAKONE (OPEN-CLOSED AND CLOSED-OPEN INTEGRALS)
C                                                   LV, 17-7-1991
C
C VERSION 0.7: READ ONE ELECTRON INTEGRALS AND OTHER INFORMATION FROM
C              MRCONEE                              LV, 10-9-1991
C
C VERSION 0.8: BUGFIX WRITING CI VECTORS (CASE OF 0 DETS IN REPRESENTATION)
C                                                   LV, 25-11-1991
C VERSION 1.0: WORKS WITH ABELIAN SYMMETRY
C                                                   LV, 17-02-1992
C VERSION 1.1: MINOR BUGFIXES LV 2-3-1992
C VERSION 1.2: OPTION TO SPECIFY NUMBER OF ELECTRONS IN IREF (LV 8-9-1992)
C VERSION 1.3: REPLACE ' AND " BY 1 AND 2 IN NAMES OF CIVECFILES (LV SEPT.'92)
C VERSION 1.4: ADAPTED TO HP WORKSTATIONS   (BDJ MAY 1993)
C VERSION 1.5: SX2 CALLS REMOVED. PROGRAM ONLY WORKS FOR HP, CRAY, CONVEX AND
C              CYBER  (BDJ JULY '93)
C VERSION 1.6: Cleaned up, Ported to IBM RS6000
C              LV (April 8, 1994)
C VERSION 1.7: Uses time-reversal unique canonical list of integrals
C              LV (June 28, 1994)
C VERSION 1.8: Now also real groups can be handled al real
      END
C
C     =====================
      SUBROUTINE DIAGH(IRP)
C     =====================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      REAL*8 E(N2),EE(N2,2)
C
      CALL CRDIAG (N2,NDET,HR,HI,D(1,IRP),EE(1,1),EE(1,2),E,E,EE,IERR)
      RETURN
      END

C     ===========================
      SUBROUTINE GENDET(IRP,IRRP)
C     ===========================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      INTEGER N(N2)
C
#include <ibtfun.h>
      IF (IPRNT.GE.10) WRITE(6,1000) REPNA(IRRP)
C
      NDET = 0
      DO 10 I = 1, 2 ** NORB
c       IF ((NBITS(I).EQ.NELEC).AND.(NBITS(I.AND.IREF).EQ.IREFE)
        IF ((NBITS(I).EQ.NELEC).AND.(NBITS(IBTAND(I,IREF)).EQ.IREFE)
     +     .AND.(IRPDET(I).EQ.IRRP))  THEN
          NDET = NDET + 1
          IF (NDET .GT. N2) STOP ' N2 TOO SMALL'
          IDET (NDET) = I
          DO 20 J=1,NORB
            N(J)=0
 20         IF (IBTAND (I, IBTSHL (1, J - 1)) .NE. 0)  N(J)=1
          IF (IPRNT.GE.10)
     +        WRITE(6,1010) NDET,(N(J),J=1,NORB)
        ENDIF
 10     CONTINUE
      RETURN
 1000 FORMAT(/1X,'LIST OF DETERMINANTS WITH SYMMETRY ',A4,':'//)
 1010 FORMAT(1X,I4,32I2)
      END
C
C     ================================
      SUBROUTINE HMAT (ID, JD, HMR, HMI)
C     ================================
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'goscom.inc'
C
#include <ibtfun.h>
C
      HMR = 0.0D0
      HMI = 0.0D0
C
      IDV = IBTXOR (ID, JD)
      IT = NBITS (IDV)
C
      IF (IT .EQ. 0) THEN
C
C     ======================
C     DETERMINANTS ARE EQUAL
C     ======================
C
      DO 30 I = 1, NORB
      IF (IBTAND (ID, IBTSHL (1, I - 1)) .NE. 0) THEN
        HMR = HMR + ONER (I, I)
        DO 20 J = I + 1, NORB
          IF (IBTAND (ID, IBTSHL (1, J - 1)) .NE. 0) THEN
C---------------------
C H = H + (I, I, J, J)
C---------------------
            INDR = INDTWR(I,I,J,J)
            INDI = INDTWI(I,I,J,J)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
            HMR = HMR + SR * TWOINR (INDR)
            HMI = HMI + SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,I,J,J)
               INDI = INDGAI(I,I,J,J)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
               HMR = HMR + SR * TWOINR (INDR)
               HMI = HMI + SI * TWOINI (INDI)
            ENDIF
C---------------------
C H = H - (I, J, J, I)
C---------------------
            INDR = INDTWR(I,J,J,I)
            INDI = INDTWI(I,J,J,I)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
            HMR = HMR - SR * TWOINR (INDR)
            HMI = HMI - SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,J,J,I)
               INDI = INDGAI(I,J,J,I)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
               HMR = HMR - SR * TWOINR (INDR)
               HMI = HMI - SI * TWOINI (INDI)
            ENDIF
          END IF
 20       CONTINUE
        END IF
 30     CONTINUE
      RETURN
      END IF
C
      IF(IT.EQ.2) THEN
C
C     ===================================
C     DETERMINANTS DIFFER BY ONE FUNCTION
C     ===================================
C
      IA = IBTAND (IDV, ID)
      IB = IBTAND (IDV, JD)
      DO 35 I = 1, NORB
        IF (IBTAND (IA, IBTSHL (1, I - 1)) .NE. 0) K = I
        IF (IBTAND (IB, IBTSHL (1, I - 1)) .NE. 0) L = I
35      CONTINUE
C
      HMR = HMR + ONER (K, L)
      HMI = HMI + ONEI (K, L)
C
      DO 40 I = 1, NORB
        IF ((I .EQ. K) .OR. (I .EQ. L)) GO TO 40
        IF (IBTAND (ID, IBTSHL (1, I - 1)) .NE. 0) THEN
C---------------------
C H = H + (I, I, K, L)
C---------------------
            INDR = INDTWR(I,I,K,L)
            INDI = INDTWI(I,I,K,L)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
            HMR = HMR + SR * TWOINR (INDR)
            HMI = HMI + SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,I,K,L)
               INDI = INDGAI(I,I,K,L)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
               HMR = HMR + SR * TWOINR (INDR)
               HMI = HMI + SI * TWOINI (INDI)
            ENDIF
C---------------------
C H = H - (I, L, K, I)
C---------------------
            INDR = INDTWR(I,L,K,I)
            INDI = INDTWI(I,L,K,I)
            SR = SIGN(1,INDR)
            SI = SIGN(1,INDI)
            INDR = ABS(INDR)
            INDI = ABS(INDI)
            HMR = HMR - SR * TWOINR (INDR)
            HMI = HMI - SI * TWOINI (INDI)
            IF (BREIT) THEN
               INDR = INDGAR(I,L,K,I)
               INDI = INDGAI(I,L,K,I)
               SR = SIGN(1,INDR)
               SI = SIGN(1,INDI)
               INDR = ABS(INDR)
               INDI = ABS(INDI)
               HMR = HMR - SR * TWOINR (INDR)
               HMI = HMI - SI * TWOINI (INDI)
            ENDIF
        END IF
 40     CONTINUE
      PH = PHASE (K, ID) * PHASE (L, JD)
      HMR = PH * HMR
      HMI = PH * HMI
      RETURN
      END IF
C
      IF(IT.EQ.4) THEN
C
C     ====================================
C     DETERMINANTS DIFFER BY TWO FUNCTIONS
C     ====================================
C
      IA = IBTAND (IDV, ID)
      IB = IBTAND (IDV, JD)
      L = 0
      N = 0
      DO 50 I = 1, NORB
        IF (IBTAND (IA, IBTSHL (1, I - 1)) .NE. 0) THEN
          K = L
          L = I
        END IF
        IF (IBTAND (IB, IBTSHL (1, I - 1)) .NE. 0) THEN
          M = N
          N = I
        END IF
 50     CONTINUE
C
C---------------------
C H = H + (K, M, L, N)
C---------------------
      INDR = INDTWR(K,M,L,N)
      INDI = INDTWI(K,M,L,N)
      SR = SIGN(1,INDR)
      SI = SIGN(1,INDI)
      INDR = ABS(INDR)
      INDI = ABS(INDI)
      HMR = HMR + SR * TWOINR (INDR)
      HMI = HMI + SI * TWOINI (INDI)
      IF (BREIT) THEN
         INDR = INDGAR(K,M,L,N)
         INDI = INDGAI(K,M,L,N)
         SR = SIGN(1,INDR)
         SI = SIGN(1,INDI)
         INDR = ABS(INDR)
         INDI = ABS(INDI)
         HMR = HMR + SR * TWOINR (INDR)
         HMI = HMI + SI * TWOINI (INDI)
      ENDIF
C---------------------
C H = H - (K, N, L, M)
C---------------------
      INDR = INDTWR(K,N,L,M)
      INDI = INDTWI(K,N,L,M)
      SR = SIGN(1,INDR)
      SI = SIGN(1,INDI)
      INDR = ABS(INDR)
      INDI = ABS(INDI)
      HMR = HMR - SR * TWOINR (INDR)
      HMI = HMI - SI * TWOINI (INDI)
      IF (BREIT) THEN
         INDR = INDGAR(K,N,L,M)
         INDI = INDGAI(K,N,L,M)
         SR = SIGN(1,INDR)
         SI = SIGN(1,INDI)
         INDR = ABS(INDR)
         INDI = ABS(INDI)
         HMR = HMR - SR * TWOINR (INDR)
         HMI = HMI - SI * TWOINI (INDI)
      ENDIF
      PH = PHASE (K, ID) * PHASE (L, ID) *
     +     PHASE (M, JD) * PHASE (N, JD)
      HMR = HMR * PH
      HMI = HMI * PH
      RETURN
      END IF
C
C     =============================================
C     DETERMINANTS DIFFER BY MORE THEN TWO FUNCTIONS
C     ==============================================
C
      HMR = 0.0D0
      HMI = 0.0D0
      RETURN
      END
C
C     ===============
      BLOCK DATA INIT
C     ===============
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      DATA REPNA / N5A * '    ' /
      DATA MULTB / N5B * 0 /
      DATA ONER, ONEI / N3 * 0.0D0, N3 * 0.0D0 /
      DATA TWOINR, TWOINI / 0.0D0, N4A * 0.0D0, 0.0D0, N4A * 0.0D0 /
      DATA INDTWR, INDTWI / N4 * 0, N4 * 0 /
      DATA INDGAR, INDGAI / N4 * 0, N4 * 0 /
      DATA IDET / N2 * 0 /
      END
C
C     =====================
      FUNCTION IRPDET(JDET)
C     =====================
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'goscom.inc'
C
      JRRP = 0
      ID=JDET
      DO 10 I=1,NORB
        IF (MOD(ID,2).EQ.1) THEN
          JRRP = MULTB(IRPAMO(I),JRRP)
        ENDIF
        ID=ID/2
 10     CONTINUE
      IRPDET=JRRP
      RETURN
      END
C
C     ================
      SUBROUTINE MAKEH
C     ================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      DO 10 I = 1, NDET
        DO 10 J = 1, I
          CALL HMAT (IDET (I), IDET (J), HR (I, J), HI (I, J))
          HR (J, I) = HR (I, J)
          HI (J, I) = - HI (I, J)
 10       CONTINUE

C     PRINT*, ' ** HR **'
C     DO 100 I=1,NDET
C       WRITE(6,1010) (HR(I,J),J=1,I)
C100    CONTINUE
C     PRINT*, ' ** HI **'
C     DO 110 I=1,NDET
C       WRITE(6,1010) (HI(I,J),J=1,I)
C110    CONTINUE
C1010 FORMAT(15F8.4)
      RETURN
      END
C
C     =============================
      INTEGER FUNCTION NBITS (JDET)
C     =============================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      ID = JDET
      NB = 0
      DO 10 I = 1, NORB
        NB = NB + MOD (ID, 2)
 10     ID = ID / 2
      NBITS = NB
      RETURN
      END
C
C     =====================
      FUNCTION PHASE(I, ID)
C     =====================
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'goscom.inc'
C
#include <ibtfun.h>
C
      P = 1.0D0
      DO 10 K = I + 1, NORB
        IF (IBTAND (ID, IBTSHL (1, K - 1)) .NE. 0) P = - P
 10     CONTINUE
      PHASE = P
      RETURN
      END
C
C     ===========================
      SUBROUTINE PRTOUT(IRP,IRRP)
C     ===========================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      PARAMETER (EVOLT=27.21139590988D0,CM=219474.63068D0)
C
      CHARACTER*1 N(N1)
C
      INTEGER IRPS(N2*N5),IND(N5)
      REAL*8 DS(N2*N5)
C
      NAMELIST/POPAN/THRESH,DEGEN,SELPOP
C
#include <ibtfun.h>
      THRESH=1.0D-3
      SELPOP=1.0D2
      DEGEN=1.0D-10
      REWIND(5)
      READ(5,POPAN,END=100)
 100  IF (IRP.EQ.0) GOTO 110
      NDETI(IRP)=NDET
      IF (IPRNT.LT.5) RETURN
C
C     ------------------------------------------------
C     PRINT POPULATION ANALYSIS FOR REPRESENTATION IRP
C     ------------------------------------------------
C
      DO 10 I=1,NDET
        DO 20 J=1,NDET
          HR(J,I)=HR(J,I)*HR(J,I)+HI(J,I)*HI(J,I)
 20       CONTINUE
 10     CONTINUE
C
      DO 70 I=1,NDET
        IF (D(I,IRP)-D(1,IRP).LT.SELPOP) MX=I
 70     CONTINUE
C
      WRITE(6,1006) REPNA(IRRP),MX
 1006 FORMAT(/1X,'POPULATION ANALYSIS FOR REPRESENTATION ',A4/
     +        1X,'           THE FIRST ',I4,' VECTORS ARE ANALYZED'//)
C
      DO 30 I=1,MX
        WRITE(6,1000) D(I,IRP)+ECORE
        DO 40 J=1,NDET
          IF (HR(J,I).GT.THRESH) THEN
            DO 50 JJ=1,NORB
              N(JJ)='0'
              IF (IBTAND(IDET(J),IBTSHL(1,JJ-1)).NE.0) N(JJ)='1'
 50           CONTINUE
            DO 60 JJ=NORB+1,N1
              N(JJ)=' '
 60           CONTINUE
            WRITE(6,1008) J,(N(JJ),JJ=1,N1),HR(J,I)
          ENDIF
 40       CONTINUE
 30     CONTINUE
      RETURN
C
C     ------------------------------------
C     PRINT COMPRESSED LIST OF EIGENVALUES
C     ------------------------------------
C
 110  CALL SORTEW (DS,NDETT,IRPS)
      WRITE(6,1005)
      ILEVEL=0
      EW=0.0D0
      EWF=DS(1)
      NEW=0
      AVER=0.0D0
      DO 120 I=1,NDETT
        AVER=AVER+DS(I)
        IF (ABS(DS(I)-EW).GT.DEGEN) THEN
          IF (NEW.NE.0)
     &      WRITE(6,1010) ILEVEL,EW-EWF,EW+ECORE,NEW
          ILEVEL=ILEVEL+1
          NEW=1
          EW=DS(I)
        ELSE
          NEW=NEW+1
          EW=((NEW-1)*EW+DS(I))/DBLE(NEW)
        ENDIF
 120    CONTINUE
      IF (NEW.NE.0)
     &WRITE(6,1010) ILEVEL,EW-EWF,EW+ECORE,NEW
      WRITE(6,1020) AVER/DBLE(NDETT)+ECORE
C
      IF (IRRP.EQ.0) WRITE(6,1025) (REPNA(I),I=NREP+1,2*NREP)
      IF (IRRP.EQ.1) WRITE(6,1025) (REPNA(I),I=1,NREP)
      ILEVEL=0
      EW=0.0D0
      EWF=DS(1)
      NEW=0
      DO 140 I=1,NDETT
        IF (ABS(DS(I)-EW).GT.DEGEN) THEN
         IF (NEW.NE.0)
     &      WRITE(6,1030) ILEVEL,EVOLT*(EW-EWF),CM*(EW-EWF),
     &                    (IND(JRP),JRP=1,NREP)
          DO 130 JRP = 1, NREP
            IND(JRP) = 0
 130        CONTINUE
          ILEVEL = ILEVEL + 1
          IND(IRPS(I)) = 1
          NEW=1
          EW=DS(I)
        ELSE
          NEW=NEW+1
          IND(IRPS(I)) = IND(IRPS(I)) + 1
          EW=((NEW-1)*EW+DS(I))/DBLE(NEW)
        ENDIF
 140    CONTINUE
      IF (NEW.NE.0) WRITE(6,1030) ILEVEL,EVOLT*(EW-EWF),CM*(EW-EWF),
     &              (IND(JRP),JRP=1,NREP)
 1000 FORMAT(1X,F16.8)
 1008 FORMAT(17X,I8,2X,16(A1),F16.8)
 1005 FORMAT(//1X,'ENERGY EIGENVALUES IN ATOMIC UNITS'
     +            //1X,'LEVEL   REL EIGENVALUE     ',
     +            'ABS EIGENVALUE  DEGENERACY'/)
 1010 FORMAT(1X,I4,1X,F16.12,F20.12,1X,'(',I4,' * )')
 1020 FORMAT(//1X,'TOTAL AVERAGE: ',F20.10)
 1025 FORMAT(//1X,'RELATIVE EIGENVALUES IN OTHER UNITS;'/' SYMMETRY',
     +       ' CLASSIFICATION IN THE ABELIAN SUBGROUP',
     +       //1X,'LEVEL  EIGENVALUE (EV) ',
     +      ' EIGENVALUE (CM-1)',2X,8A4/)
 1030 FORMAT(1X,I4,1X,F18.9,F18.6,1X,8(I3,1X))
      RETURN
      END
C
C     ================
      SUBROUTINE RDTWO
C     ================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      PARAMETER(D0 = 0.0D0)
      CHARACTER*10 DATEX,TIMEX*8
      DIMENSION KR(-N1:N1)
      COMPLEX*16 RPHASE(N1),ONE
      INTEGER INDK(N3),INDL(N3)
C
      MDINT = 1
      NUNIQ = 1
      TWOINR(0) = D0
      TWOINI(0) = D0
      OPEN(MDINT,FILE='MDCINT',FORM='UNFORMATTED')
      Read (mdint,err=10000,end=10000) datex,timex,nkr,
     &                                 (KR(I),KR(-I),I=1,NKR) 
C--!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      nkrcoul = nkr
  10  read (mdint,ERR=10010,END=10010) ikr,jkr,nz,
     &             (indk(inz),indl(inz),inz=1,nz),
     &             (twoinr(inz),twoini(inz),inz=nuniq+1,nuniq+nz)
      If (nuniq+nz.GT.N4A) Call Q2Job (3,'P.RDTWO','N4A',nuniq+nz)
      If (ikr.EQ.0) GoTo 20 ! End of File
C-->  Distribute this set and generate time-reversed integrals.
C-->  K|ia> = |ib> ; K|ib> = - |ia>
      i = kr(ikr)
      itr = kr(-ikr)
      j = kr(jkr)
      jtr = kr(-jkr)
      SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
      do inz = 1, nz
         nuniq = nuniq + 1
         kkr = indk(inz)
         k = kr(kkr)
         ktr = kr(-kkr)
         lkr = indl(inz)
         l = kr(lkr)
         ltr = kr(-lkr)
         SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral plus time-reversed partners
         INDTWR(I,J,K,L) = NUNIQ
         INDTWR(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDTWR(I,J,LTR,KTR) = NUNIQ * SignKL
         INDTWR(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
         INDTWI(I,J,K,L) = NUNIQ
         INDTWI(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDTWI(I,J,LTR,KTR) = NUNIQ * SignKL
         INDTWI(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
C--> Complex conjugate plus time-reversed partners
         INDTWR(J,I,L,K) = NUNIQ
         INDTWR(ITR,JTR,L,K) = NUNIQ * SignIJ
         INDTWR(J,I,KTR,LTR) = NUNIQ * SignKL
         INDTWR(ITR,JTR,KTR,LTR) = NUNIQ * SignIJ * SignKL
         INDTWI(J,I,L,K) = - NUNIQ
         INDTWI(ITR,JTR,L,K) = - NUNIQ * SignIJ
         INDTWI(J,I,KTR,LTR) = - NUNIQ * SignKL
         INDTWI(ITR,JTR,KTR,LTR) = - NUNIQ * SignIJ * SignKL
C--> Particle interchanged plus time-reversed partners
         INDTWR(K,L,I,J) = NUNIQ
         INDTWR(LTR,KTR,I,J) = NUNIQ * SignKL
         INDTWR(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDTWR(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
         INDTWI(K,L,I,J) = NUNIQ
         INDTWI(LTR,KTR,I,J) = NUNIQ * SignKL
         INDTWI(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDTWI(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
C--> Particle interchanged and complex conjugated plus time-reversed partners
         INDTWR(L,K,J,I) = NUNIQ
         INDTWR(KTR,LTR,J,I) = NUNIQ * SignKL
         INDTWR(L,K,ITR,JTR) = NUNIQ * SignIJ
         INDTWR(KTR,LTR,ITR,JTR) = NUNIQ * SignIJ * SignKL
         INDTWI(L,K,J,I) = - NUNIQ
         INDTWI(KTR,LTR,J,I) = - NUNIQ * SignKL
         INDTWI(L,K,ITR,JTR) = - NUNIQ * SignIJ
         INDTWI(KTR,LTR,ITR,JTR) = - NUNIQ * SignIJ * SignKL
      Enddo
      Goto 10
 20   Continue
      CLOSE(1)
      Write (6,1000) "Coulomb",datex,timex,nuniq-1
      If (.NOT.BREIT) Return
      NUNIQC = NUNIQ
      OPEN(MDINT,FILE='MDBINT',FORM='UNFORMATTED')
      Read (mdint,err=10001,end=10001) datex,timex,nkr,
     &                                 (KR(I),KR(-I),I=1,NKR) 
      If (NKR.NE.NKRCOUL) Then
         Write (6,*) 'Basis set size on Coulomb file :',nkrcoul
         Write (6,*) 'Basis set size on Gaunt file :  ',nkr
         Stop 'Gaunt and Coulomb 2-el files incompatible'
      Endif
  30  read (mdint,ERR=10011,END=10011) ikr,jkr,nz,
     &             (indk(inz),indl(inz),inz=1,nz),
     &             (twoinr(inz),twoini(inz),inz=nuniq+1,nuniq+nz)
      If (nuniq+nz.GT.N4A) Call Q2Job (3,'P.RDTWO','N4A',nuniq+nz)
      If (ikr.EQ.0) GoTo 40 ! End of File
C-->  Distribute this set and generate time-reversed integrals.
C-->  K|ia> = |ib> ; K|ib> = - |ia>
      i = kr(ikr)
      itr = kr(-ikr)
      j = kr(jkr)
      jtr = kr(-jkr)
      SignIJ = - SIGN(1,ikr) * SIGN(1,jkr)
      do inz = 1, nz
         nuniq = nuniq + 1
         kkr = indk(inz)
         k = kr(kkr)
         ktr = kr(-kkr)
         lkr = indl(inz)
         l = kr(lkr)
         ltr = kr(-lkr)
         SignKL = - SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral plus time-reversed partners
         INDGAR(I,J,K,L) = NUNIQ
         INDGAR(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDGAR(I,J,LTR,KTR) = NUNIQ * SignKL
         INDGAR(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
         INDGAI(I,J,K,L) = NUNIQ
         INDGAI(JTR,ITR,K,L) = NUNIQ * SignIJ
         INDGAI(I,J,LTR,KTR) = NUNIQ * SignKL
         INDGAI(JTR,ITR,LTR,KTR) = NUNIQ * SignIJ * SignKL
C--> Complex conjugate plus time-reversed partners
         INDGAR(J,I,L,K) = NUNIQ
         INDGAR(ITR,JTR,L,K) = NUNIQ * SignIJ
         INDGAR(J,I,KTR,LTR) = NUNIQ * SignKL
         INDGAR(ITR,JTR,KTR,LTR) = NUNIQ * SignIJ * SignKL
         INDGAI(J,I,L,K) = - NUNIQ
         INDGAI(ITR,JTR,L,K) = - NUNIQ * SignIJ
         INDGAI(J,I,KTR,LTR) = - NUNIQ * SignKL
         INDGAI(ITR,JTR,KTR,LTR) = - NUNIQ * SignIJ * SignKL
C--> Particle interchanged plus time-reversed partners
         INDGAR(K,L,I,J) = NUNIQ
         INDGAR(LTR,KTR,I,J) = NUNIQ * SignKL
         INDGAR(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDGAR(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
         INDGAI(K,L,I,J) = NUNIQ
         INDGAI(LTR,KTR,I,J) = NUNIQ * SignKL
         INDGAI(K,L,JTR,ITR) = NUNIQ * SignIJ
         INDGAI(LTR,KTR,JTR,ITR) = NUNIQ * SignIJ * SignKL
C--> Particle interchanged and complex conjugated plus time-reversed partners
         INDGAR(L,K,J,I) = NUNIQ
         INDGAR(KTR,LTR,J,I) = NUNIQ * SignKL
         INDGAR(L,K,ITR,JTR) = NUNIQ * SignIJ
         INDGAR(KTR,LTR,ITR,JTR) = NUNIQ * SignIJ * SignKL
         INDGAI(L,K,J,I) = - NUNIQ
         INDGAI(KTR,LTR,J,I) = - NUNIQ * SignKL
         INDGAI(L,K,ITR,JTR) = - NUNIQ * SignIJ
         INDGAI(KTR,LTR,ITR,JTR) = - NUNIQ * SignIJ * SignKL
      Enddo
      Goto 30
 40   Continue
      CLOSE(1)
      Write (6,1000) "Gaunt",datex,timex,nuniq-nuniqc
      RETURN
 1000 Format (/1X,A7,' integral file was generated at ',A10,1x,A8,
     &/' Read',I8,' unique integrals')
10000 STOP 'ERROR READING HEADER OF MDCINT'
10001 STOP 'ERROR READING HEADER OF MDBINT'
10010 STOP 'ERROR READING INTEGRALS FROM MDCINT'
10011 STOP 'ERROR READING INTEGRALS FROM MDBINT'
      END
C
C     =================
      SUBROUTINE READGOS
C     =================
      IMPLICIT REAL*8 (A-H, O-Z)
C
      INCLUDE 'goscom.inc'
C
      CHARACTER*14 REPN(16)
C
      LOGICAL TWOC
      INTEGER REFDET,IRPMO(N1)
      REAL*8 SELPOP
C
      NAMELIST/GOSCIP/NELEC,REFDET,IREFE,IPRNT
      NAMELIST/GENERAL/TWOC
      NAMELIST/POPAN/THRESH,DEGEN,SELPOP
C
      MRCONEE=75
C
      REWIND(5)
      NELEC=0
      IPRNT=0
      REFDET=0
      IREFE=0
      READ (5,GOSCIP,END=200)
  200 IREF=REFDET
      IF (IREFE.LE.0) IREFE = NBITS(IREF)
      TWOC=.FALSE.
      REWIND(5)
      READ(5,GENERAL,END=201)
  201 WRITE(6,GENERAL)
C
      REWIND(5)
      THRESH=1.0D-3
      DEGEN=1.0D-10
      SELPOP=1.0D2
      REWIND(5)
      READ(5,POPAN,END=100)
  100 WRITE(6,POPAN)
C
      OPEN (MRCONEE,FILE='MRCONEE',FORM='UNFORMATTED')
      READ (MRCONEE) NORB,BREIT,ECORE
      IF (NORB.GT.N1) CALL Q2JOB(3,'P.READIN','N1',NORB)
      READ (MRCONEE) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
      READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
      READ (MRCONEE) ((MULTB(I,J),I=1,2*NREP),J=1,2*NREP)
      READ (MRCONEE) (IRPMO(IMO),IRPAMO(IMO),EPS(IMO),IMO=1,NORB)
      READ (MRCONEE)
     & ((ONER(IMO,JMO),ONEI(IMO,JMO),JMO=1,NORB),IMO=1,NORB)
      CLOSE (MRCONEE)
C
      DO 250 I = 1, 2 * NREP
      MULTB(I,0) = I
  250 CONTINUE
C
      WRITE(6,1200)
      WRITE(6,1210) NELEC,IREF
      WRITE(6,1220) BREIT,ECORE,NORB
      DO 300 I=1,NORB
      WRITE(6,1230) I,REPN(IRPMO(I)),EPS(I)
 300  CONTINUE
C
      RETURN
 1200 FORMAT(//1X,'ALL INFORMATION READ FROM MRCONEE (AND INPUT)')
 1210 FORMAT(//1X,'NUMBER OF ELECTRONS: ',I4/
     +         1X,'DECIMAL REPRESENTATION OF FROZEN ORBITALS:',I4)
 1220 FORMAT(//1X,'BREIT INTERACTION INCLUDED :',T50,L1
     +        /1X,'CORE ENERGY :',T30,G20.10
     +        /1X,'NUMBER OF ACTIVE ORBITALS :',T47,I4
     +       //1X,'ORBITAL',T10,'REPRESENTATION',T32,'ORBITAL ENERGY')
 1230 FORMAT(2X,I4,T10,A14,T30,G20.10)
      END
C
C     ==============================================
      SUBROUTINE SORTEW (DS,NDETT,IRPS)
C     ==============================================
      IMPLICIT REAL*8 (A-H,O-Z)
C
      INCLUDE 'goscom.inc'
C
      REAL*8 DS(N2*N5)
      INTEGER NDETT,IND(N5),IRPS(N2*N5)
C
      INDEX = 1
      IRPMIN = 1
      DO 1 IRP = 2, NREP
    1   IF (D(1,IRP).LT.D(1,IRPMIN)) IRPMIN = IRP
C
      DO 2 IRP = 1, NREP
    2   IND(IRP) = 1
C
      NDETT = 0
      DO 3 IRP = 1, NREP
        NDETT = NDETT + NDETI(IRP)
    3   CONTINUE
C
 10   DS(INDEX)=D(IND(IRPMIN),IRPMIN)
      IRPS(INDEX)=IRPMIN
      INDEX=INDEX+1
      IF (INDEX.GT.NDETT) GOTO 30
      IND(IRPMIN)=IND(IRPMIN)+1
      IF (IND(IRPMIN).GT.NDETI(IRPMIN)) D(IND(IRPMIN),IRPMIN)=1.D12
      DO 20 IRP = 1, NREP
      IF (D(IND(IRPMIN),IRPMIN).GT.D(IND(IRP),IRP)) THEN
        IRPMIN=IRP
      ENDIF
 20   CONTINUE
      GOTO 10
 30   CONTINUE
C
      RETURN
      END
C
C     ===========================
      SUBROUTINE WRVEC (IRP,IRRP)
C     ===========================
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     WRITE TRIAL VECTORS FOR DIRRCI TO FILE MRCTRIV
C     RECORD 1 : RECORD LENGTH,NUMBER OF VECTORS, ENERGIES
C     RECORD 2 : BIT REPRESENTATION DETERMINANTS
C     RECORD 3 : FIRST VECTOR
C     RECORD 4 : SECOND VECTOR, ETC.
C
      INCLUDE 'goscom.inc'
C
      CHARACTER*11 VECFIL
      DIMENSION DV (N2)
      NAMELIST/CIVEC/NVEC
C
      IF (NDET.EQ.0) RETURN
      DEGEN = 1.D-10
      NVEC = 5
      REWIND(5)
      READ (5,CIVEC,END=1,ERR=1)
    1 CONTINUE
      IF (NVEC.GT.NDET) NVEC = NDET
C
C     DETERMINE RECORD LENGTH
C
      LEN = MAX0(NDET,2)
      LENREC = 2 * LEN * 8
      VECFIL = 'MDTRIV_'//REPNA(IRRP)
      IF (VECFIL(8:8).EQ.'''') VECFIL(8:8) = '1'
      IF (VECFIL(8:8).EQ.'"') VECFIL(8:8) = '2'
      OPEN (10,FILE=VECFIL,ACCESS='DIRECT',RECL=LENREC)
      WRITE (10,REC=2) NDET,(IDET(I),I=1,NDET)
      IVEC = 0
      EPREV = D(1,IRP) - DEGEN - 1.D0
      DO 10 I = 1, NDET
      IF (D(I,IRP)-EPREV.GT.DEGEN) THEN
         IVEC = IVEC + 1
         DV (IVEC) = D(I,IRP) + ECORE
         WRITE (10,REC=IVEC+2) (HR(J,I),HI(J,I),J=1,NDET)
         EPREV = D(I,IRP)
      ENDIF
      IF (IVEC.EQ.NVEC) GOTO 100
   10 CONTINUE
  100 WRITE (10,REC=1) LENREC,NVEC,(DV(IVEC),IVEC=1,NVEC)
      CLOSE (10)
      WRITE (6,1000) REPNA(IRRP),NVEC,VECFIL,(DV(IVEC),IVEC=1,NVEC)
      RETURN
 1000 FORMAT (/' REPRESENTATION ',A4,'; ',I2,
     &' VECTORS WRITTEN TO FILE ',A12//' EIGENVALUES ',(T20,F20.12))
      END
