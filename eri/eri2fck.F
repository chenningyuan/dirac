C  /* Deck erifck */
      SUBROUTINE HRFCK2(FMAT,DMAT,NDMT,ISYMDM,IFCTYP,PARFCK,IPRFCK,
     *                  WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <aovec.h>
#include <maxaqn.h>
#include <mxorb.h>
C
C	I don't know if IFCTYP and ISYMDM is used properly in this 
C	routine at the current stage, K.Ruud-jan96
C
      LOGICAL PARFCK
      DIMENSION FMAT(*), DMAT(*), WORK(LWORK), IFCTYP(NDMT),
     &          ISYMDM(NDMT)
#include <ccom.h>
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <erimem.h>
#include <aobtch.h>
#include <odbtch.h>
#include <symmet.h>
#include <iratdef.h>
C
      CALL TIMER('START ',TIMSTR,TIMEND)
C
C     Initialization in ER2INP
C
      CALL ER2INP('*ER2FCK')
C
      IPRINT = MAX(IPRERI,IPRFCK)
      THRSH  = MAX(THRS,1.00D-15)
      NDMAT  = NDMT
C
      IF (PARFCK) THEN
         WRITE (LUPRI,'(/,A,/A)') ' Parallel not implemented in ERIFCK',
     &                            ' Process aborted.'
         CALL QUIT('ERIFCK. not parallelized yet.')
      END IF
      IF (IFCTYP(1).NE.1) THEN
         WRITE (LUPRI,'(/,A,I5,/A)')
     &        ' IFCTYP different from 1 or 3 in ERIFCK',
     &          IFCTYP,
     &        ' Process aborted.'
         CALL QUIT('IFCTYP .ne. 1 in ERIFCK.')
      END IF
      IF (NDMAT .GT. 1) THEN
         WRITE (LUPRI,'(/A,I5,/A)')
     &        ' NDMAT .gt. 1 in ERIFCK.', NDMAT, ' Process aborted'
         CALL QUIT('NDMAT .gt. 1 in ERIFCK')
      END IF
      IF (ISYMDM(1) .NE. 0) THEN
         WRITE (LUPRI, '(/A,I5,/A)')
     &        ' ISYMDM .ne. 1 in ERIFCK', ISYMDM(1), ' Process aborted'
         CALL QUIT('ISYMDM .ne. 1 in ERIFCK')
      END IF
C
C     Memory
C
      MEMOK  = .TRUE.
      MEMADD = 0
      MODAB  = 0
      MODCD  = 0
C
C     Symmetry transformation
C
      SYMTRA = DOSTRA .OR. MAXREP .GT. 0
C
      WRTINT = .FALSE.
      FCKINT = .TRUE.
C
C     Pi
C
      CALL GETPI
C
C     AO batches
C     ==========
C
      CALL SETAOB(WORK,LWORK,IPRINT)
C
C     Density matrix (not needed yet)
C     ===============================
C
C     KDAO  = 1
C     KLAST = KDAO + N2BASX*NDMAT
C     IF(KLAST.GT.LWORK) CALL STOPIT('HERFCK','DSOTAO',KLAST,LWORK)
C     LWRK  = LWORK - KLAST + 1
C     DO 200 IMAT = 1,NDMAT
C        JDAO = KDAO + (IMAT-1)*N2BASX
C        CALL DSOTAO(DMAT(1,IMAT),WORK(JDAO),NBAST,ISYMDM,IPRINT)
C 200 CONTINUE
C
C     OD batches
C     ==========
C
C     This subroutine returns several arrays for each electron
C     starting at addresses K????1 and K????2. These are to be
C     transferred to ODCDRV.
C
      KWRK = 1
      LWRK = LWORK
      CALL ODBCHS(KODCL1,KODCL2,
     &            KODBC1,KODBC2,KRDBC1,KRDBC2,
     &            KODPP1,KODPP2,KRDPP1,KRDPP2,
     &            KFREE,LFREE,WORK(KWRK),LWRK,IPRINT)
C
      IF (IPRINT .GT. 2) THEN
         WRITE (LUPRI,'(2(/,2X,A,I10))')
     &      ' Memory requirements for ODBCHS:',LWRK - LFREE,
     &      ' Memory left for ODCDRV:        ',LFREE
      END IF
C
      ICALL = 0
      CALL GETDST(ICALL,IPRINT)
C
C     Select integrals to be calculated
C     =================================
C
      CALL PICKAO(IPRINT)
C
C     Information about distributions
C     ===============================
C
      CALL ERIDSI(IPRINT)
C
      KLAST = KFREE
      LWRK  = LFREE
C
C     Calculate integrals
C     ===================
C
      IF (.NOT.INTSKP) THEN
         CALL ODCDRV(WORK(KODCL1),WORK(KODCL2),
     &               WORK(KODBC1),WORK(KODBC2),
     &               WORK(KRDBC1),WORK(KRDBC2),
     &               WORK(KODPP1),WORK(KODPP2),
     &               WORK(KRDPP1),WORK(KRDPP2),
     &               FMAT,DMAT,WORK(KLAST),LWRK,IPRINT)
C
C        Error message in case of insufficient memory
C
         IF (.NOT.MEMOK) THEN
            WRITE (LUPRI,'(//,1X,A,3(/,1X,A,I10))')
     &         ' Not enough memory for this run of ERIFCK.',
     &         ' Available memory in ERIFCK:',LWORK,
     &         ' Required memory for ERIFCK:',LWORK + MEMADD,
     &         ' Increase memory (LWORK) by:',MEMADD
            WRITE (LUPRI,'(/,1X,A,2I5)')
     &         ' Memory requirements largest for OD classes :',
     &           MODAB,MODCD
            CALL QUIT('Insufficient memory in ERIFCK.')
         END IF
      END IF
C
C     Symmetrize Fock matrix
C     ======================
C
      CALL ERISFK(FMAT,NBASIS)
C
C     Print densities and Fock matrix
C     ===============================
C
      IF (IPRINT.GT.4) THEN
         CALL HEADER('Density and Fock matrices in ERIFCK',-1)
         KSTR = 1
         DO 300 I = 1, NDMAT
            WRITE (LUPRI,'(//,1X,A,I3)') ' Density matrix No.',I
            CALL OUTPUT(DMAT(KSTR),1,NBASIS,1,NBASIS,NBASIS,NBASIS,1,
     &                  LUPRI)
            WRITE (LUPRI,'(//,1X,A,I3)') ' Fock matrix No.',I
            CALL OUTPUT(FMAT(KSTR),1,NBASIS,1,NBASIS,NBASIS,NBASIS,1,
     &                  LUPRI)
            KSTR = NBASIS*NBASIS
  300    CONTINUE
      END IF
C
      CALL TIMER('ERIFCK',TIMSTR,TIMEND)
      CALL FLSHFO(LUPRI)
C
      RETURN
      END
C  /* Deck erifok */
      SUBROUTINE ERIFOK(SO,IPNTCR,IODDCC,IPNTUV,FMAT,DMAT,WORK,
     &                  LWORK,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
#include <eridst.h>
      DIMENSION SO(*), IPNTCR(MAXBCH,4),
     &          IPNTUV(KC2MAX,2,2), WORK(LWORK), IODDCC(NRTOP),
     &          FMAT(NBASIS,NBASIS), DMAT(NBASIS,NBASIS)
#include <cbieri.h>
#include <aobtch.h>
#include <ericom.h>
#include <eribuf.h>
#include <symmet.h>
#include <hertop.h>
C
C     Allocation for ERIFOK
C
      LBIN   = NCCX*KHKTA*KHKTB*KHKTC*KHKTD
      KBIN   = 1
      KIBIN  = KBIN   +  LBIN
      KINDEX = KIBIN  + (4*LBIN - 1)/IRAT + 1
      KLAST  = KINDEX + (4*LBIN - 1)/IRAT + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ERIFOK',' ',KLAST,LWORK)
      CALL ERIFO1(SO,WORK(KINDEX),IPNTCR,IODDCC,IPNTUV,
     &            WORK(KBIN),WORK(KIBIN),LBIN,FMAT,DMAT,IPRINT)
C
      RETURN
      END
C  /* Deck erifo1 */
      SUBROUTINE ERIFO1(SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,LBIN,FMAT,DMAT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
C
#include <eridst.h>
      DIMENSION SO(*), INDEX(LBIN,4),
     &          IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,2,2),
     &          BIN(LBIN), IBIN(LBIN,4),
     &          FMAT(NBASIS,NBASIS,NDMAT), DMAT(NBASIS,NBASIS,NDMAT)
#include <cbieri.h>
#include <ericom.h>
#include <eribuf.h>
#include <aobtch.h>
#include <hertop.h>
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERIFO1',-1)
C
C     Collect (non-zero) integrals and attach indices
C     ===============================================
C
      CALL ERINDF(SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &            BIN,IBIN,LBIN,INT,IPRINT)
C
C     Construct Fock matrix contribution
C     ==================================
C
      CALL FOKDI1(FMAT,DMAT,NDMAT,BIN,IBIN,LBIN,INT,IPRINT)
C
      RETURN
      END
C  /* Deck erindf */
      SUBROUTINE ERINDF(SO,INDEX,IPNTCR,IODDCC,IPNTUV,
     &                  BIN,IBIN,LBIN,INT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
#include <mxcent.h>
#include <mxorb.h>
#include <aovec.h>
      INTEGER A, B, C, D, I, A1, B1, C1, D1, AB, CD, R, S, T
      LOGICAL DOREP(0:7,4), CTRIAB, CTRICD, CTRIAC, CTRIBD, CTRIPQ
      DIMENSION SO(NCCS,MLTPR,MLTPS,MLTPT,KHKTAB,KHKTCD),
     &          INDEX(NCCS,4),
     &          IPNTCR(MAXBCH,4),
     &          IODDCC(NRTOP), IPNTUV(KC2MAX,2,2),
     &          BIN(LBIN), IBIN(LBIN,4),
     &          IPNRST(0:7,3),
     &          IADCMP(MXAQN,MXAQN,2)
#include <cbieri.h>
#include <ericom.h>
#include <erithr.h>
#include <aobtch.h>
#include <hertop.h>
#include <symmet.h>
C
#include <ibtfun.h>
      IBTEST(I,J,K,L) = IBTAND(I,IBTXOR(J,ISYMAO(K,L)))
      IBTREP(I,J,K,L) = IBTXOR(J,IBTAND(I,ISYMAO(K,L)))
C
      IF (IPRINT .GT. 6) CALL HEADER('Subroutine ERINDF',-1)
C
      IF (IPRINT .GT. 15) THEN
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCX,NPPBCX ',NPQBCX,NPPBCX
         WRITE (LUPRI,'(2X,A,2I5)') ' NPQBCS,NPPBCS ',NPQBCS,NPPBCS
      END IF
C
      CALL PRPREP(DOREP(0,1),NHKTA,KHKTA,ISTBLA)
      CALL PRPREP(DOREP(0,2),NHKTB,KHKTB,ISTBLB)
      CALL PRPREP(DOREP(0,3),NHKTC,KHKTC,ISTBLC)
      CALL PRPREP(DOREP(0,4),NHKTD,KHKTD,ISTBLD)
C
      CALL CMPADR(IADCMP(1,1,1),KHKTA,KHKTB,TKMPAB)
      CALL CMPADR(IADCMP(1,1,2),KHKTC,KHKTD,TKMPCD)
C
      CALL GETRST(IPNRST(0,1),ISTBLR)
      CALL GETRST(IPNRST(0,2),ISTBLS)
      CALL GETRST(IPNRST(0,3),ISTBLT)
C
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/,2X,A,8L2)')'DOREP A  ',(DOREP(I,1),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP B  ',(DOREP(I,2),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP C  ',(DOREP(I,3),I=0,MAXREP)
         WRITE (LUPRI,'(2X,A,8L2)')  'DOREP D  ',(DOREP(I,4),I=0,MAXREP)
      END IF
C
      NCCS1 = (NPQBCS - NPPBCS)*NCTFAB*NCTFCD
      NCCS2 = NPPBCS*NCTFAB*NCTFCD
C
      INT = 0
C
      IF (NCCS1 .GT. 0) THEN
C
         DO 100 A = 0, MAXREP
         IF (DOREP(A,1)) THEN
            DO 101 I = 1, NCCS1
               INDEX(I,1) = INDXBT(KNDXBT(IPNTCR(I,1)),A)
  101       CONTINUE
         DO 110 B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            DO 111 I = 1, NCCS1
               INDEX(I,2) = INDXBT(KNDXBT(IPNTCR(I,2)),B)
  111       CONTINUE
         DO 120 C = 0, MAXREP
         IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
            D = IBTXOR(IBTXOR(A,B),C)
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 190
            IF (DIAGCD .AND. D.GT.C) GO TO 190
C
            DO 121 I = 1, NCCS1
               INDEX(I,3) = INDXBT(KNDXBT(IPNTCR(I,3)),C)
               INDEX(I,4) = INDXBT(KNDXBT(IPNTCR(I,4)),D)
  121       CONTINUE
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            IA   = -1
            MAXB = KHKTB
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
            DO 300 ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO 310 ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,2,1))
                  DO 320 ICMPC = 1, KHKTC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRICD) MAXD = ICMPC
                     DO 330 ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,2,2))
                        IF (IODDAB .EQ. IODDCD) THEN
                           IF (WRTSCR) THEN
                              DO 400 I = 1, NCCS1
                                 SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD)
                                 IF (ABS(SOABCD) .GT. THRSH) THEN
                                    INT = INT + 1
                                    BIN (INT  ) = SOABCD
                                    IBIN(INT,1) = INDEX(I,1) + IA
                                    IBIN(INT,2) = INDEX(I,2) + IB
                                    IBIN(INT,3) = INDEX(I,3) + IC
                                    IBIN(INT,4) = INDEX(I,4) + ID
                                 END IF
  400                         CONTINUE
                           ELSE
                              DO 450 I = 1, NCCS1
                                 INT = INT + 1
                                 BIN(INT) = SO(I,R,S,T,ICMPAB,ICMPCD)
                                 IBIN(INT,1) = INDEX(I,1) + IA
                                 IBIN(INT,2) = INDEX(I,2) + IB
                                 IBIN(INT,3) = INDEX(I,3) + IC
                                 IBIN(INT,4) = INDEX(I,4) + ID
  450                         CONTINUE
                           END IF
                        END IF
                     END IF
  330                CONTINUE
                  END IF
  320             CONTINUE
               END IF
  310          CONTINUE
            END IF
  300       CONTINUE
C
  190       CONTINUE
C
         END IF
  120    CONTINUE
         END IF
  110    CONTINUE
         END IF
  100    CONTINUE
      END IF
C
      IF (NCCS2 .GT. 0) THEN
C
         DO 105 A = 0, MAXREP
         IF (DOREP(A,1)) THEN
            DO 106 I = NCCS1 + 1, NCCS
               INDEX(I,1) = INDXBT(KNDXBT(IPNTCR(I,1)),A)
  106       CONTINUE
         DO 115 B = 0, MAXREP
         IF (DOREP(B,2)) THEN
            DO 116 I = NCCS1 + 1, NCCS
               INDEX(I,2) = INDXBT(KNDXBT(IPNTCR(I,2)),B)
  116       CONTINUE
         DO 125 C = 0, MAXREP
         IF (DOREP(C,3) .AND. DOREP(IBTXOR(IBTXOR(A,B),C),4)) THEN
            D = IBTXOR(IBTXOR(A,B),C)
            CD = IBTXOR(C,D)
C
            IF (DIAGAB .AND. B.GT.A) GO TO 195
            IF (DIAGCD .AND. D.GT.C) GO TO 195
            IF (C.GT.A .OR. (C.EQ.A .AND. D.GT.B)) GO TO 195
C
            DO 126 I = NCCS1 + 1, NCCS
               INDEX(I,3) = INDXBT(KNDXBT(IPNTCR(I,3)),C)
               INDEX(I,4) = INDXBT(KNDXBT(IPNTCR(I,4)),D)
  126       CONTINUE
C
            R = IPNRST(B,1)
            S = IPNRST(D,2)
            T = IPNRST(CD,3)
C
            MAXB = KHKTB
            MAXC = KHKTC
            MAXD = KHKTD
            CTRIAB = DIAGAB .AND. A.EQ.B
            CTRICD = DIAGCD .AND. C.EQ.D
            CTRIAC = A.EQ.C
            CTRIBD = B.EQ.D
            CTRIPQ = A.EQ.C .AND. B.EQ.D
C
            IA = -1
            DO 305 ICMPA = 1, KHKTA
            IVARA = IBTEST(ISTBLA,A,NHKTA,ICMPA)
            IF (IVARA.EQ.0) THEN
               IA = IA + 1
               IB = -1
               IF (CTRIAB) MAXB = ICMPA
               DO 315 ICMPB = 1, MAXB
               IVARB = IBTEST(ISTBLB,B,NHKTB,ICMPB)
               IF (IVARB.EQ.0) THEN
                  IB = IB + 1
                  IC = -1
                  ICMPAB = IADCMP(ICMPA,ICMPB,1)
                  IODDAB = IODDCC(IPNTUV(ICMPAB,2,1))
                  IF (CTRIAC) MAXC = ICMPA
                  DO 325 ICMPC = 1, MAXC
                  IVARC = IBTEST(ISTBLC,C,NHKTC,ICMPC)
                  IF (IVARC.EQ.0) THEN
                     IC = IC + 1
                     ID = -1
                     IF (CTRIPQ .AND. ICMPA.EQ.ICMPC) THEN
                        MAXD = ICMPB
                     ELSE
                        MAXD = KHKTD
                        IF (CTRICD) MAXD = ICMPC
                     END IF
                     DO 335 ICMPD = 1, MAXD
                     IVARD = IBTEST(ISTBLD,D,NHKTD,ICMPD)
                     IF (IVARD.EQ.0) THEN
                        ID = ID + 1
                        ICMPCD = IADCMP(ICMPC,ICMPD,2)
                        IODDCD = IODDCC(IPNTUV(ICMPCD,2,2))
                        IF (IODDAB .EQ. IODDCD) THEN
                           IF (WRTSCR) THEN
                              DO 405 I = NCCS1 + 1, NCCS
                                 SOABCD = SO(I,R,S,T,ICMPAB,ICMPCD)
                                 IF (ABS(SOABCD) .GT. THRSH) THEN
                                    INT = INT + 1
                                    BIN (INT  ) = SOABCD
                                    IBIN(INT,1) = INDEX(I,1) + IA
                                    IBIN(INT,2) = INDEX(I,2) + IB
                                    IBIN(INT,3) = INDEX(I,3) + IC
                                    IBIN(INT,4) = INDEX(I,4) + ID
                                 END IF
  405                         CONTINUE
                           ELSE
                              DO 455 I = NCCS1 + 1, NCCS
                                 INT = INT + 1
                                 BIN(INT) = SO(I,R,S,T,ICMPAB,ICMPCD)
                                 IBIN(INT,1) = INDEX(I,1) + IA
                                 IBIN(INT,2) = INDEX(I,2) + IB
                                 IBIN(INT,3) = INDEX(I,3) + IC
                                 IBIN(INT,4) = INDEX(I,4) + ID
  455                         CONTINUE
                           END IF
                        END IF
                     END IF
  335                CONTINUE
                  END IF
  325             CONTINUE
               END IF
  315          CONTINUE
            END IF
  305       CONTINUE
C
  195       CONTINUE
C
         END IF
  125    CONTINUE
         END IF
  115    CONTINUE
         END IF
  105    CONTINUE
      END IF
      RETURN
      END
C  /* Deck fokdi1 */
      SUBROUTINE FOKDI1(FMAT,DMAT,NDMAT,BUF,IBUF,LBIN,LENGTH,
     &                  IPRINT)
C
C     Henrik Koch and Trygve Helgaker 18-NOV-1991.
C
C     This subroutine adds derivative two-electron integrals to
C     Fock matrices. The Fock matrices are assumed
C     to be square matrices in full dimension without symmetry reduction
C     in size. Remember to zero out the fock matrices before starting
C     to accumulate.
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D2 = 2.0D00, DP5 = 0.50D00, DP25 = 0.25D00)
      INTEGER P, Q, R, S
#include <inforb.h>
      DIMENSION FMAT(NBAST,NBAST), DMAT(NBAST,NBAST),
     &          BUF(LBIN), IBUF(LBIN,4)
#include <ibtfun.h>
C
      DO 100 INT = 1, LENGTH
         DINT = D2*BUF(INT)
         P    = IBUF(INT,1)
         Q    = IBUF(INT,2)
         R    = IBUF(INT,3)
         S    = IBUF(INT,4)
C
         IF  (P.EQ.Q)                    DINT = DP5*DINT
         IF  (R.EQ.S)                    DINT = DP5*DINT
         IF ((P.EQ.R .AND. Q.EQ.S) .OR.
     &       (P.EQ.S .AND. Q.EQ.R))      DINT = DP5*DINT
c         IF (IFCTYP .EQ. 1) THEN
C
            EINT = DP25*DINT
            FMAT(P,Q) = FMAT(P,Q) + DINT*DMAT(R,S)
            FMAT(R,S) = FMAT(R,S) + DINT*DMAT(P,Q)
            FMAT(P,R) = FMAT(P,R) - EINT*DMAT(S,Q)
            FMAT(P,S) = FMAT(P,S) - EINT*DMAT(R,Q)
            FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(S,P)
            FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(R,P)
c         ELSE IF (IFCTYP .EQ. 3) THEN
c            GRS = DINT*(DMAT(R,S) + DMAT(S,R))
c            FMAT(P,Q) = FMAT(P,Q) + GRS
c            FMAT(Q,P) = FMAT(Q,P) + GRS
c            GPQ = DINT*(DMAT(P,Q) + DMAT(Q,P))
c            FMAT(R,S) = FMAT(R,S) + GPQ
c            FMAT(S,R) = FMAT(S,R) + GPQ
c            EINT = DP5*DINT
c            FMAT(R,P) = FMAT(R,P) - EINT*DMAT(Q,S)
c            FMAT(S,P) = FMAT(S,P) - EINT*DMAT(Q,R)
c            FMAT(R,Q) = FMAT(R,Q) - EINT*DMAT(P,S)
c            FMAT(S,Q) = FMAT(S,Q) - EINT*DMAT(P,R)
c            FMAT(P,R) = FMAT(P,R) - EINT*DMAT(S,Q)
c            FMAT(P,S) = FMAT(P,S) - EINT*DMAT(R,Q)
c            FMAT(Q,R) = FMAT(Q,R) - EINT*DMAT(S,P)
c            FMAT(Q,S) = FMAT(Q,S) - EINT*DMAT(R,P)
c         END IF
  100 CONTINUE
      RETURN
      END
C  /* Deck erisfk */
      SUBROUTINE ERISFK(FMAT,NBASIS)
#include <implicit.h>
#include <priunit.h>
      DIMENSION FMAT(NBASIS,NBASIS)
      DO 100 I = 1, NBASIS
      DO 100 J = 1, I
         FMT = FMAT(I,J) + FMAT(J,I)
         FMAT(I,J) = FMT
         FMAT(J,I) = FMT
  100 CONTINUE
      RETURN
      END
