      SUBROUTINE RSBB1E(ISCSM,ISCTP,ISBCTP,ICCSM,ICCTP,ICBCTP,
     &                  IGRP,NROW,
     &                  NGAS,ISEL,ICEL,
     &                  SB,CB,
     &                  ADSXA,SXSTST,STSTSX,
     &                  MXPNGAS,NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,H,
     &                  NSMOB,NSMST,NSMSX,MXPOBS,MOC,
     &                  NSCOL,MXSXST,IEXSTR,FACSTR,IH2TRM,SCLFAC,NTESTG)
*
* One electron excitations on column strings
*. If IH2TRM .ne. 0 then the diagonal and one-electron
*  excitations arising from the two body operator is also
* included
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ISBCTP : Base for sigma column types
* ICCSM,ICCTP : Symmetry and type of C     columns
* ICBCTP : Base for C     column types
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active sets
* ISEL : Occupation in each active set for sigma block
* ICEL : Occupation in each active set for C     block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals,strings,
*       single excitations, double excitations
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* MOC  : Use MOC method ( instead of N-1 resolution method )
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* H : Space for one electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),SXSTST(NSMSX,NSMST),
     &        STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Specific Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),I1(MAXK,*),XI1S(MAXK,*),H(*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
      DIMENSION IEXSTR(MXSXST,*),FACSTR(MXSXST,*)
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
      CALL QENTER('RS1  ')
* Type of single excitations that connects the two column strings
C     MOC = 1
      NTESTL = 00
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500)THEN
        WRITE(6,*)
        WRITE(6,*) ' ======================= '
        WRITE(6,*) ' Information from RSBB1E '
        WRITE(6,*) ' ======================= '
        WRITE(6,*)
        WRITE(6,*) ' RSBB1E : MOC,IH2TRM ', MOC,IH2TRM
        WRITE(6,*) ' ISEL : '
        CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
        WRITE(6,*) ' ICEL : '
        CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
      END IF
*
C?    IF(SCLFAC.NE.1.0D0) THEN
C?      WRITE(6,*) ' Warning : RSBB1E, SCLFAC .ne 1 '
C?    END IF
      IFRST = 1
      JFRST = 1
*. Types of single excitations that connect ISEL and ICEL
      CALL SXTYP_GAS(NSXTP,ITP,JTP,NGAS,MXPNGAS,ISEL,ICEL)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = STSTSX(ISCSM,ICCSM)
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP=  1, NSXTP
        ITYP = ITP(IJTP)
        JTYP = JTP(IJTP)
        IF(NTEST.GE.2000)
     &  write(6,*) ' ITYP JTYP ', ITYP,JTYP
        DO 800 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
*. New intermediate strings will be accessed so
          KFRST = 1
          IF(JSM.EQ.0) GOTO 800
          IF(NTEST.GE.2000)
     &    write(6,*) ' ISM JSM ', ISM,JSM
          NIORB = NOBPTS(ITYP,ISM)
          NJORB = NOBPTS(JTYP,JSM)
          IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
          IF(MOC.EQ.0) THEN
*
*
* ======================================================================
*.                   Use N-1 resolution method
* ======================================================================
*
*
*. Loop over partitionings of the row strings
          NIPART = NROW/MAXI
          IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Loop over partitionings of N-1 strings
            KBOT = 1-MAXK
            KTOP = 0
  700       CONTINUE
              KBOT = KBOT + MAXK
              KTOP = KTOP + MAXK
*.set up I1(K,JORB) =  XI1S(K,JORB) a JORB !J STRING >
              JBORB = IOBPTS(JTYP,JSM)
              J12 = 1
              K12 = 1
              IF(JFRST.EQ.1) KFRST= 1
              CALL ADST_GAS(JBORB,NJORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &             I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,J12,K12,
     &             KKJACT,SCLFAC)
              JFRST = 0
              KFRST = 0
*.set up I2(K, IORB) = a IORB !I STRING >
              IBORB = IOBPTS(ITYP,ISM)
              I12 = 2
              IF(IFRST.EQ.1) KFRST = 1
              ONE = 1.0D0
              CALL ADST_GAS(IBORB,NIORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &        I2,XI2S,MAXK,NKBTC,KEND,IFRST,KFRST,
     &        I12,K12,KKIACT,ONE   )
              IFRST = 0
              KFRST = 0
*. This is the place to start over partitioning of I strings
          DO 701 IPART = 1, NIPART
            IBOT = (IPART-1)*MAXI+1
            ITOP = MIN(IBOT+MAXI-1,NROW)
            NIBTC = ITOP - IBOT + 1

* Obtain CSCR(I,K,JORB) = SUM(J)<K!A JORB!J>C(I,J)
              DO JJORB = 1,NJORB
                ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                     NKBTC,I1(1,JJORB),XI1S(1,JJORB))
              END DO
*.Obtain one electron integrals (JSM,JTP,ISM,ITP)
              CALL GETH1(H,JSM,JTYP,ISM,ITYP)
              IF(NTEST.GE.1000) THEN
                WRITE(6,*) ' RSBB1E H BLOCK '
                CALL WRTMAT(H,NJORB,NIORB,NJORB,NIORB)
              END IF
*.Sscr(I,K,i) = CSCR(I,K,j)*h(j,i)
              NIK = NIBTC*NKBTC
              FACTORC = 0.0D0
              FACTORAB = 1.0D0
              CALL MATML7(SSCR,CSCR,H,NIK,NIORB,NIK,NJORB,NJORB,NIORB,
     &                    FACTORC,FACTORAB,0)
*.S(I,a+ K) =  S(I, a+ K) + sgn*Sscr(I,K,i)
*.set up I1(IORB,K) = a IORB !I STRING >
              DO IIORB = 1,NIORB
                ISBOFF = 1+(IIORB-1)*NKBTC*NIBTC
                CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &               NKBTC,I2(1,IIORB),XI2S(1,IIORB))
              END DO
  701     CONTINUE
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
C 701     CONTINUE
*. End of loop over I partitioninigs
          ELSE
*
*
* ===================================================================
*                       Use MOC algorithm
* ===================================================================
*
*
COLD        CALL GETH1(H,ISM,ITYP,JSM,JTYP)
            CALL GETH1(H,JSM,JTYP,ISM,ITYP)
*. Construct list of string matrix elements => sum(ij) h(ji)* <J!E(ji)!I>
            CALL H1STR(H,ISCSM,ISCTP,IGRP,JSM,JTYP,
     &                ISM,ITYP,MXSXST,IEXSTR,FACSTR,
     &                I1,IH2TRM,NTEST)
            DO ICOL = 1, NSCOL
*. Msaxpy update
               CALL MSAXPYN6(SB((ICOL-1)*NROW+1),
     &              FACSTR(1,ICOL),CB,NROW,
     &              I1(ICOL,1),IEXSTR(1,ICOL))
            END DO
*
          END IF
*.(End of algorithm switch)
  800   CONTINUE
*.(end of loop over symmetries)
  900 CONTINUE
 1001 CONTINUE
*
      CALL QEXIT('RS1  ')
C!    WRITE(6,*) ' Enforced stop in RSBB1E '
C!    STOP' Enforced stop in RSBB1E '
      RETURN
      END
