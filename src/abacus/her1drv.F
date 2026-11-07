!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

C
C
#ifdef UNDEF
/* Comdeck log */
tsaue - 941119 : array DOINT implemented
#endif
C  /* Deck oneinp */
      SUBROUTINE ONEINP(WORD)
C
C     Trygve Helgaker, University of Oslo, Norway
C
      use cbione_mod, only: CUT, IPRINT, NODC, NODV, SKIP
      use priunit_mod, only: LUCMD, LUPRI
!      use abainf_mod, only: IPRDEF
      use kinds_mod, only: dp
#include "implicit.h"
      PARAMETER (NTABLE = 6)
      LOGICAL SET, NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "abainf.h"
      SAVE SET
      DATA TABLE /'.SKIP  ', '.PRINT ','XXXXXXX',
     *            '.NODC  ', '.NODV  ','.STOP  '/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (WORD .NE. '*END OF') THEN
 969        READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .NE. '*') GO TO 969
         END IF
         RETURN
      END IF
C
      SET = .TRUE.
      CALL ONEINI
C
      NEWDEF = WORD .EQ. '*ONEINT'
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ONEINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ONEINP')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
                  READ (LUCMD, '(I5)') IPRINT
                  IF (IPRINT .EQ. IPRDEF) ICHANG = ICHANG - 1
               GO TO 100
    3             CONTINUE
               GO TO 100
    4             NODC = .TRUE.
               GO TO 100
    5             NODV = .TRUE.
               GO TO 100
    6             CUT  = .TRUE.
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ONEINP.'
            END IF
            CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
            CALL QUIT('Illegal prompt in ONEINP')
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for ONEINT:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' ONEINT skipped in this run.'
         ELSE
            IF (IPRINT .NE. IPRDEF) THEN
               WRITE (LUPRI,'(A,I5)') ' Print level in ONEINT:',IPRINT
            END IF
            IF (NODC) WRITE (LUPRI,'(/,2A)') ' Inactive one-electron',
     *      ' density matrix neglected in ONEINT.'
            IF (NODV) WRITE (LUPRI,'(/,2A)') ' Active one-electron',
     *      ' density matrix neglected in ONEINT.'
            IF (CUT) THEN
               WRITE (LUPRI,'(/,A)') ' Program is stopped after ONEINT.'
            END IF
         END IF
      END IF
      RETURN
      END
C  /* Deck oneini */
      SUBROUTINE ONEINI
C
C     Initialize /CBIONE/
C
!      use abainf_mod, only: DIPDER, IPRDEF, MOLGRD, MOLHES
      use cbione_mod, only: CUT, DIFDIP, IPRINT, MAXDIF, NODC, NODV,
     &  SKIP
#include "implicit.h"
#include "abainf.h"
C
      IPRINT = IPRDEF
      SKIP   = .FALSE.
      CUT    = .FALSE.
      IF (MOLHES) THEN
         MAXDIF = 2
      ELSE IF (DIPDER .OR. MOLGRD) THEN
         MAXDIF = 1
      ELSE
         SKIP = .TRUE.
      END IF
      NODC   = .FALSE.
      NODV   = .FALSE.
      DIFDIP = DIPDER
      RETURN
      END
C  /* Deck oneint */
      SUBROUTINE ONEINT(WORK,LWORK,PASS,PROPTY)
C
C     TUH
C
      use cbione_mod, only: CUT, DIFDIP, DIFINT, IPRINT, MAXDIF, NODC,
     &  NODV, SKIP
      use priunit_mod, only: LUPRI
#include "implicit.h"
      LOGICAL PASS, PROPTY
      DIMENSION WORK(LWORK)
      IF (SKIP) RETURN
      IF (IPRINT .GE. 2) CALL TIMER('START ',TIMSTR,TIMEND)
      IF (IPRINT .GT. 2) CALL TITLER('Output from ONEINT','*',103)
      PROPTY = .TRUE.
      CALL ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,NODV,
     &            DIFDIP)
      IF (IPRINT .GE. 2) CALL TIMER ('ONEINT',TIMSTR,TIMEND)
      PASS = .TRUE.
      IF (CUT) THEN
         WRITE (LUPRI,'(/,A)')
     &          ' Program stopped after ONEINT as required.'
         WRITE (LUPRI,'(A)') ' No restart file has been written.'
         CALL QUIT(' ***** End of ABACUS (in ONEINT) *****')
      END IF
      RETURN
      END
C  /* Deck onedrv */
      SUBROUTINE ONEDRV(WORK,LWORK,IPRINT,PROPTY,MAXDIF,DIFINT,NODC,
     &                  NODV,DIFDIP)
C
      use nuclei_mod, only: NUCDEP
      use symmet_mod, only: MAXREP
      use maxorb_mod, only: MXCORB, MXSHEL
      use shells_mod, only: ISTBAO, KHKT, KMAX
!      use cbisol_mod, only: LMNTOT, SOLVNT
      use maxaqn_mod, only: MXAQN, MXQN
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use iratdef_mod, only: IRAT
#include "implicit.h"
C
      LOGICAL DIFINT, NODC, NODV, DIFDIP, PROPTY, DOINT(2,2)
      CHARACTER*4 OMITVNUC(2)
      DIMENSION WORK(LWORK)
#include "cbisol.h"
C
C
      DOINT(1,1) = .TRUE.
      DOINT(2,1) = .TRUE.
      DOINT(1,2) = .TRUE.
      DOINT(2,2) = .TRUE.
      OMITVNUC(1)   = 'FFFF'
      OMITVNUC(2)   = 'FFFF'
C
C     ***** Number of basis functions *****
C
      NBAST  = 0
      NNBAST = 0
      DO 100 KB = 0, MAXREP
         NBASI = 0
         DO 200 ISHELL = 1, KMAX
            IF (IAND(KB,ISTBAO(ISHELL)).EQ.0) THEN
               NBASI = NBASI + KHKT(ISHELL)
            END IF
  200    CONTINUE
         NBAST  = NBAST  + NBASI
         NNBAST = NNBAST + NBASI*(NBASI + 1)/2
  100 CONTINUE
      NNBASX = NBAST*(NBAST + 1)/2
C
      KSTHMA = 1
      KDENMA = KSTHMA + 3*NNBASX
      KFOCMA = KDENMA +   NNBASX
      KFACIN = KFOCMA +   NNBASX
      KCOORC = KFACIN + 2*NUCDEP
C     allocate for nuclear charge AND an electronic core charge
C     for modified nuclear attraction (e.g. for the small charge
C     in relativistic calculations). /Mar2001 hjaaj
      KSIGNC = KCOORC + 3*NUCDEP
      KNCENT = KSIGNC + 3*NUCDEP
      KJSYMC = KNCENT +  (NUCDEP + 1)/IRAT
      KJCENT = KJSYMC +  (NUCDEP + 1)/IRAT
      KGEXP  = KJCENT +  (NUCDEP + 1)/IRAT
      KTLMD  = KGEXP  +  NUCDEP
      IF (SOLVNT .AND. MAXDIF .GE. 2) THEN
         KLAST  = KTLMD + LMNTOT*6*NUCDEP
      ELSE
         KLAST  = KTLMD
      END IF
      LWRK   = LWORK  - KLAST + 1
      IF (KLAST .GT. LWORK) CALL STOPIT('ONEDRV',' ',KLAST,LWORK)
      CALL ONEDR1(WORK(KSTHMA),WORK(KDENMA),WORK(KFOCMA),WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KTLMD),WORK(KLAST),
     &            LWRK,IPRINT,PROPTY,MAXDIF,NODC,NODV,
     &            DIFDIP,NBAST,NNBASX,NNBAST,DOINT,OMITVNUC,.TRUE.)
      RETURN
      END
C  /* Deck onedr1 */
      SUBROUTINE ONEDR1(STHMAT,DENMAT,FOCMAT,FACINT,COORC,SIGNC,NCENTC,
     &                  JSYMC,JCENTC,GEXP,TLMD,WORK,LWORK,IPRINT,
     &                  PROPTY,MAXDIF,NODC,NODV,DIFDIP,
     &                  NBAST,NNBASX,NNBAST,DOINT,OMITVNUC,TOFILE)
C
C     TUH
C
C     This program calculates overlap and one-electron Hamiltonian
C     matrix elements and the first and second derivatives of these
C     elements using the McMurchie/Davidson scheme.  See L. E. McMurchie
C     & E. R. Davidson, J. Comp. Phys. 26 (1978) 218, and also V. R.
C     Saunders in "Methods in Computational Molecular Physics", Reidel
C     1983.
C
C     Symmetry included  880406  TUH & PRT
C
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use symmet_mod, only: FMULT, ISYMAO, ISYMAX, MAXOPR, MAXREP, MULT,
     &  PT
!      use energy_mod, only: ENERKE, ENERNA, NA
      use lmns_mod, only: LVALUA, LVALUB, MVALUA, MVALUB, NVALUA, NVALUB
!      use ccom_mod, only: THRS
      use onecom_mod
!      use dcbgrd_mod, only: DODERM, GRADKN, GRADNU, GRADRO, I, INTCLASS
      use shells_mod, only: CENT, ISTBAO, JSTRT, KCKT, KHKT, KMAX,
     &  LCLASS, NCENT, NHKT, NUCO, NUMCF, SPHR
      use priunit_mod, only: LUPRI
      use nuclei_mod, only: CHARGE, CORD, GNUEXP, ISTBNU, NUCDEP,
     &  NUCIND, NUCNUM
      use symind_mod, only: INDFA, INDFB, ISOFRA, ISOFRB
      use dipole_mod, only: DDIPE, DDIPS
!      use cbisol_mod, only: LMNTOT, SOLVNT
      use csym1_mod, only: BUF, IBUF, INDMAX, LENGTH
      use maxaqn_mod, only: MXAQN, MXQN
!      use dorps_mod, only: DOREPS
      use maxorb_mod, only: MXCORB, MXSHEL
      use kinds_mod, only: dp
#include "implicit.h"
      PARAMETER (LUTEMP = 48)
      REAL(dp), PARAMETER :: D0 = 0.00D00
C
      LOGICAL NODC, NODV, FRSDER, SECDER, DIFDIP, PROPTY, DOINT(2,2),
     &        LOMITVNUC(2,5),TOFILE
      CHARACTER*4 OMITVNUC(2)
      DIMENSION DENMAT(*), FOCMAT(*), STHMAT(NNBASX,3),
     &          FACINT(NUCDEP,*), WORK(LWORK), COORC(3,NUCDEP),
     &          SIGNC(3,NUCDEP), NCENTC(NUCDEP), JSYMC(NUCDEP),
     &          JCENTC(NUCDEP), GEXP(NUCDEP),TLMD(3*NUCDEP*LMNTOT,2)
#ifdef PRG_DIRAC
#include "dcbgrd.h"
#else
#include "energy.h"
#endif
#include "ccom.h"
#include "dorps.h"
#include "cbisol.h"
C
      IF (IPRINT .GE. 5) CALL TITLER('Output from ONEDR1','*',103)
C
      TIMHER = D0
      TIMINT = D0
      FRSDER = MAXDIF .GE. 1
      SECDER = MAXDIF .EQ. 2
      TOLS   = THRS*THRS
      TOLOG  = - LOG(TOLS)
      IF (PROPTY .AND. (SECDER .OR. DIFDIP)) THEN
         REWIND LUTEMP
         INDMAX = 0
         LENGTH = 0
      END IF
      DO I = 1,2
         LOMITVNUC(I,1) = (OMITVNUC(I)(1:1) .EQ. 'T')
         LOMITVNUC(I,2) = (OMITVNUC(I)(2:2) .EQ. 'T')
         LOMITVNUC(I,3) = (OMITVNUC(I)(3:3) .EQ. 'T')
         LOMITVNUC(I,4) = (OMITVNUC(I)(4:4) .EQ. 'T')
         LOMITVNUC(I,5) = LOMITVNUC(I,1) .OR. LOMITVNUC(I,2)
     &             .OR. LOMITVNUC(I,3) .OR. LOMITVNUC(I,4)
      ENDDO
C
Cjth
#if !defined(PRG_DIRAC)
C
C     ***** Nuclear contributions *****
C
      IF (PROPTY) CALL DIPNUC(IPRINT,DIFDIP)
C
#endif
C
C     Gradient and Hessian elements equal to zero
C
      IF (PROPTY) THEN
         ENERKE = D0
         ENERNA = D0
#ifdef PRG_DIRAC
C
         CALL DZERO(GRADKN(1,0),3*NUCDEP)
         CALL DZERO(GRADKN(1,1),3*NUCDEP)
         CALL DZERO(GRADKN(1,2),3*NUCDEP)
C
         CALL DZERO(GRADNU(1,0),3*NUCDEP)
         CALL DZERO(GRADNU(1,1),3*NUCDEP)
         CALL DZERO(GRADNU(1,2),3*NUCDEP)
C
C
         CALL DZERO(GRADRO(1,0),3*NUCDEP)
         CALL DZERO(GRADRO(1,1),3*NUCDEP)
         CALL DZERO(GRADRO(1,2),3*NUCDEP)
C
#else 
         CALL DZERO(GRADKE,3*NUCDEP)
         CALL DZERO(GRADNA,3*NUCDEP)
         CALL DZERO(GRADFS,3*NUCDEP)
         IF (SECDER) THEN
            HESSKE(:,:) = 0.0D0
            HESSNA(:,:) = 0.0D0
            HESFS2(:,:) = 0.0D0
         END IF
         IF (SOLVNT) THEN
            ESOLTT = D0
            ESOLNN = D0
            CALL DZERO(GSOLTT,3*NUCDEP)
            CALL DZERO(GSOLNN,3*NUCDEP)
            IF (SECDER) THEN
               HSOLT2(:,:) = 0.0D0
               HSOLNN(:,:) = 0.0D0
            END IF
         END IF
#endif
         IF (DIFDIP) THEN
            CALL DZERO(DDIPE,9*NUCDEP)
            CALL DZERO(DDIPS,9*NUCDEP)
         END IF
      END IF
C
Cjth
C     The contravariant Fock matrix and the density matrix are
C     allready given from DIRAC, so skip the next.
C
#if !defined(PRG_DIRAC)
C
C     **************************************************************
C     ***** Set up total density and Fock matrices in AO basis *****
C     **************************************************************
C
      IF (PROPTY) THEN
         KDSO  = 1
         KFSO  = KDSO  + NNBAST
         KLAST = KFSO  + NNBAST
         LWRK  = LWORK - KLAST + 1
         IF (KLAST.GT.LWORK) CALL STOPIT('ONEDR1','DSOFSO',KLAST,LWORK)
         CALL DSOFSO(WORK(KDSO),WORK(KFSO),WORK(KLAST),LWRK,IPRINT,NODC,
     &               NODV)
         CALL DSYM1(DENMAT,FOCMAT,WORK(KDSO),WORK(KFSO),NBAST,IPRINT)
      END IF
#endif
      KLAST = 1
      KFCM  = KLAST
#if !defined(PRG_DIRAC)
C
C     ********************************************************
C     ***** Preparation for solvent contributions        *****
C     ***** Calculation of solvent nuclear contributions *****
C     ********************************************************
C
      IF (PROPTY .AND. SOLVNT) THEN
         KFRSAV = 1
         KFREE  = KFRSAV
         LFREE  = LWORK
         CALL MEMGET('REAL',KFCM,LMNTOT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KGLM, LMTOT,WORK,KFREE,LFREE)
         CALL MEMGET('REAL',KTLM, LMTOT,WORK,KFREE,LFREE)
         LWRK  = LWORK - KFREE + 1
         CALL SOLFL(WORK(KGLM),EPDIEL,RCAV,LCAVMX)
         CALL RDSIFC('SOLTLM',WORK(KTLM),WORK(KFREE),LWRK)
         CALL FCMFAC(LCAVMX,LMNTOT,LMTOT,WORK(KGLM),WORK(KTLM),
     &               WORK(KFCM),WORK(KFREE),LWRK,IPRINT)
         CALL MEMREL('ONEDR1.FCMFAC',WORK,KFRSAV,KTLM,KFREE,LFREE)
         IF (MAXDIF .GE. 2) CALL DZERO(TLMD,6*NUCDEP*LMNTOT)
         CALL SOLNUC(.TRUE.,MAXDIF,TLMD(1,2),
     &               WORK(KFCM),WORK,KFREE,LFREE,IPRINT)
C        CALL SOLNUC(PROPTY,MAXDIF,TLMND,FCM,WORK,KFRSAV,LFRSAV,IPRINT)
         KLAST = KFREE
      END IF
#endif
      LWRK  = LWORK - KLAST + 1
C
C     ************************************************************
C     ***** Triangular loop over symmetry independent shells *****
C     ************************************************************
C
      CALL DZERO(STHMAT,3*NNBASX)
      CALL IZERO(ISOFRA, 8)
      IDENA = 0
      DO 100 ISHELA = 1,KMAX
         NHKTA = NHKT(ISHELA)
         KHKTA = KHKT(ISHELA)
         KCKTA = KCKT(ISHELA)
         ICA   = LCLASS(ISHELA)
         SPHRA = SPHR(ISHELA)
         CALL LMNVAL(NHKTA,KCKTA,LVALUA,MVALUA,NVALUA)
         NCENTA = NCENT(ISHELA)
         ICENTA = NUCNUM(NCENTA,1)
         MULA   = ISTBAO(ISHELA)
         MULTA  = MULT(MULA)
         NUCA   = NUCO(ISHELA)
         NUMCFA = NUMCF(ISHELA)
         JSTA   = JSTRT(ISHELA)
         CORAX  = CENT(ISHELA,1,1)
         CORAY  = CENT(ISHELA,2,1)
         CORAZ  = CENT(ISHELA,3,1)
         IDENB0 = 0
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included.
C
         DO 600 I = 1, 8
            ISOFRB(I) = 0
            DO 610 J = 1, MXAQN
              INDFA(I,J) = -10 000 000
610         CONTINUE
600      CONTINUE
         DO 620 NA = 1, KHKTA
            DO 630 IREP = 0, MAXREP
            IF (IAND(MULA,IEOR(IREP,ISYMAO(NHKTA,NA))).EQ.0) THEN
               ISOFRA(IREP+1)    = ISOFRA(IREP+1) + 1
               INDFA (IREP+1,NA) = ISOFRA(IREP+1)
            END IF
630         CONTINUE
620      CONTINUE
         IF (IPRINT .GT. 20) THEN
            WRITE(LUPRI,'(A,I4)')' IA address offsets for shell ',ISHELA
            DO 640 NA = 1,KHKTA
               WRITE(LUPRI,'(8(1X,I5))') (INDFA(I,NA), I = 1,MAXREP+1)
640         CONTINUE
         END IF
      DO 110 ISHELB = 1,ISHELA
         LDIAG = ISHELA .EQ. ISHELB
         NHKTB = NHKT(ISHELB)
         KHKTB = KHKT(ISHELB)
         KCKTB = KCKT(ISHELB)
         ICB   = LCLASS(ISHELB)
         SPHRB = SPHR(ISHELB)
         CALL LMNVAL(NHKTB,KCKTB,LVALUB,MVALUB,NVALUB)
         NCENTB = NCENT(ISHELB)
         NHKTAB = NHKTA + NHKTB
         MULB   = ISTBAO(ISHELB)
         MULTB  = MULT(MULB)
         NUCB   = NUCO(ISHELB)
         NUMCFB = NUMCF(ISHELB)
         JSTB   = JSTRT(ISHELB)
         CORBX0 = CENT(ISHELB,1,1)
         CORBY0 = CENT(ISHELB,2,1)
         CORBZ0 = CENT(ISHELB,3,1)
         KHKTAB = KHKTA*KHKTB
         KCKTAB = KCKTA*KCKTB
         MAB    = IOR(MULA,MULB)
         KAB    = IAND(MULA,MULB)
         HKAB   = FMULT(KAB)
C
         SPHRAB = SPHRA .OR. SPHRB
C
C        Compute symmetry integral pointers for contributions
C        from this block.  Note that at present this assumes all
C        components from a shell are included
C
         DO 700 I = 1, 8
            DO 710 J = 1, MXAQN
              INDFB(I,J) = -10 000 000
710         CONTINUE
700      CONTINUE
         DO 720 NB = 1, KHKTB
            DO 730 IREP = 0, MAXREP
            IF (IAND(MULB,IEOR(IREP,ISYMAO(NHKTB,NB))).EQ.0) THEN
               ISOFRB(IREP+1)    = ISOFRB(IREP+1) + 1
               INDFB (IREP+1,NB) = ISOFRB(IREP+1)
            END IF
730         CONTINUE
720      CONTINUE
         IF (IPRINT .GT. 20) THEN
            WRITE(LUPRI,'(A,I4)')' IB address offsets for shell ',ISHELB
            DO 740 NB = 1, KHKTB
               WRITE(LUPRI,'(8(1X,I5))') (INDFB(I,NB), I = 1,MAXREP+1)
740         CONTINUE
         ENDIF
#ifdef PRG_DIRAC
C        jth - CLASS (which integral class: LL,LS,SL,SS)
C                                           0  1  1  2
         IF (ICA .EQ. 1) THEN
            IF (ICB .EQ. 1) THEN
               INTCLASS = 0
            ELSE
               INTCLASS = 1
            END IF
         ELSE 
            IF (ICB .EQ. 1) THEN
               INTCLASS = 1
            ELSE
               INTCLASS = 2
            END IF
         END IF           
         DODERM = .FALSE.
C        jth - if it is a LS integral and we have requested PROPTY 
C              (= derivatives) then calculate the derivatives of
C              the momentum integral. If DODERM then we do not do
C              (derivatives) of nuclear attraction and \beta matrix
C              etc.
         IF ( (ICA.NE.ICB) .AND. PROPTY ) THEN
            DODERM = .TRUE.
         ELSE IF (.NOT.DOINT(ICA,ICB)) THEN
            GOTO 110
         ENDIF
#else
         IF(.NOT.DOINT(ICA,ICB)) GOTO 110
#endif
         IF (IPRINT .GE. 05) WRITE (LUPRI, 1000) ISHELA, ISHELB
         IF (IPRINT .GE. 10) THEN
             WRITE (LUPRI,'(A,2I10)') ' NHKT   ', NHKTA, NHKTB
             WRITE (LUPRI,'(A,2I10)') ' KHKT   ', KHKTA, KHKTB
             WRITE (LUPRI,'(A,2I10)') ' KCKT   ', KCKTA, KCKTB
             WRITE (LUPRI,'(A,2I10)') ' NCENT  ', NCENTA, NCENTB
             WRITE (LUPRI,'(A,2I10)') ' ISTBAO ', MULA, MULB
             WRITE (LUPRI,'(A,2I10)') ' MULT   ', MULTA, MULTB
             WRITE (LUPRI,'(A,2I10)') ' NUC    ', NUCA, NUCB
             WRITE (LUPRI,'(A,2I10)') ' NUMCF  ', NUMCFA, NUMCFB
             WRITE (LUPRI,'(A,2I10)') ' JST    ', JSTA, JSTB
             WRITE (LUPRI,'(A,2F12.6)') ' CORAX    ', CORAX, CORBX0
             WRITE (LUPRI,'(A,2F12.6)') ' CORAY    ', CORAY, CORBY0
             WRITE (LUPRI,'(A,2F12.6)') ' CORAZ    ', CORAZ, CORBZ0
         END IF
C
C        Initialization for nuclear attraction integrals
C
         JMAX = NHKTAB - 2
         IF (PROPTY) JMAX = JMAX + MAXDIF
         ISTEPU = JMAX + 1
         ISTEPV = ISTEPU*ISTEPU
         NAHGTF = ISTEPU*ISTEPV
         NATOMC = 0
         DO 120 IATOMC = 1,NUCIND
            MULC   = ISTBNU(IATOMC)
            MABC   = IOR(MULC,KAB)
            CORCX0 = CORD(1,IATOMC)
            CORCY0 = CORD(2,IATOMC)
            CORCZ0 = CORD(3,IATOMC)
            CHARG1 = CHARGE(IATOMC)
            FACTOR = - FMULT(IAND(MULC,KAB))/HKAB
            DO 130 ISYMOP = 0, MAXOPR
               IF (IAND(ISYMOP,MABC) .EQ. 0) THEN
                  NATOMC = NATOMC + 1
                  JSYMC(NATOMC)   = ISYMOP
                  JCENTC(NATOMC)  = IATOMC
                  SIGNC(1,NATOMC) = PT(IAND(ISYMAX(1,1),ISYMOP))
                  SIGNC(2,NATOMC) = PT(IAND(ISYMAX(2,1),ISYMOP))
                  SIGNC(3,NATOMC) = PT(IAND(ISYMAX(3,1),ISYMOP))
                  COORC(1,NATOMC) = SIGNC(1,NATOMC)*CORCX0
                  COORC(2,NATOMC) = SIGNC(2,NATOMC)*CORCY0
                  COORC(3,NATOMC) = SIGNC(3,NATOMC)*CORCZ0
                  GEXP(NATOMC)    = GNUEXP(IATOMC)
                  FACINT(NATOMC,1)= FACTOR*CHARG1
                  FACINT(NATOMC,2)= D0
                  NCENTC(NATOMC)  = NUCNUM(IATOMC,ISYMOP+1)
               END IF
  130       CONTINUE
  120    CONTINUE
C         IF (NATOMC.NE.NUCDEP)
C     &      print *,'ONEDR1 NATOMC, NUCDEP:',NATOMC,NUCDEP
C
         CALL ONESOP(STHMAT,DENMAT,FOCMAT,FACINT,COORC,
     &               WORK(KLAST),LWRK,
     &               IPRINT,PROPTY,MAXDIF,IDENB0,CORBX0,CORBY0,CORBZ0,
     &               DIFDIP,SECDER,NATOMC,TOLOG,TOLS,JSYMC,JCENTC,
     &               NCENTC,SIGNC,GEXP,NNBASX,WORK(KFCM),TLMD,LOMITVNUC)
  110    IDENB0 = IDENB0 + KHKTB*MULTB
         IDENA = IDENA + KHKTA*MULTA
  100 CONTINUE
C
C     ***** End loop over symmetry independent orbitals *****
C
C     ****************************************************
C     ***** Write final buffers and dipole integrals *****
C     ***** ELSE write undifferentiated integrals    *****
C     ****************************************************
C
      IF (PROPTY) THEN
         IF (SECDER .OR. DIFDIP) THEN
C
C           Write final buffer of first derivative integrals
C
            IF (LENGTH .GT. 0) WRITE (LUTEMP) BUF, IBUF, LENGTH
            WRITE (LUTEMP) BUF, IBUF, -LENGTH
            IF (IPRINT .GT. 5) WRITE (LUPRI,'(/A,I4,A)')
     *       ' Last buffer of length',LENGTH,
     *       ' has been written on LUTEMP in ONEDR1.'
         END IF
      ELSE
C
C        Write undifferentiated one-electron integrals
C
        IF(TOFILE) THEN
          DO 800 I = 1, NNBASX
             STHMAT(I,2) = STHMAT(I,2) + STHMAT(I,3)
 800      CONTINUE
          CALL WRTUND(STHMAT,NBAST,NNBASX,IPRINT)
        ENDIF
      END IF
C
C     *******************************
C     ***** Symmetrize Hessians *****
C     *******************************
C
#ifndef PRG_DIRAC
      IF (PROPTY .AND. SECDER) THEN
         DO 500 I = 2,3*NUCDEP
            DO 510 J = 1,I-1
               HESSKE(I,J) = HESSKE(I,J) + HESSKE(J,I)
               HESSKE(J,I) = HESSKE(I,J)
               HESSNA(I,J) = HESSNA(I,J) + HESSNA(J,I)
               HESSNA(J,I) = HESSNA(I,J)
               HESFS2(I,J) = HESFS2(I,J) + HESFS2(J,I)
               HESFS2(J,I) = HESFS2(I,J)
               IF (SOLVNT) THEN
                  HSOLT2(I,J) = HSOLT2(I,J) + HSOLT2(J,I)
                  HSOLT2(J,I) = HSOLT2(I,J)
               END IF
  510       CONTINUE
  500    CONTINUE
      END IF
#endif
C
C     *************************
C     ***** Print Section *****
C     *************************
C
#ifdef PRG_DIRAC
      IF (PROPTY .AND. (IPRINT .GE. 5)) THEN
         KCSTRA = 1
         KSCTRA = KCSTRA + 9*NUCDEP*NUCDEP
         KLAST  = KSCTRA + 9*NUCDEP*NUCDEP
         CALL HEADER('Momentum integral gradient',-1)
         CALL PRIGRD(GRADKN(1,1),WORK(KCSTRA),WORK(KSCTRA))
         CALL HEADER('Nuclear attraction integral gradient',-1)
         CALL PRIGRD(GRADNU(1,1),WORK(KCSTRA),WORK(KSCTRA))
         CALL PRIGRD(GRADNU(1,2),WORK(KCSTRA),WORK(KSCTRA))
         CALL HEADER('Reorthonormalization gradient',-1)
         CALL PRIGRD(GRADRO(1,1),WORK(KCSTRA),WORK(KSCTRA))
         CALL PRIGRD(GRADRO(1,2),WORK(KCSTRA),WORK(KSCTRA))
         CALL HEADER('Beta matrix gradient',-1)
         CALL PRIGRD(GRADKN(1,2),WORK(KCSTRA),WORK(KSCTRA))
#else 
      IF (PROPTY .AND. (IPRINT .GT. 0)) THEN
         CALL HEADER('Kinetic energy integral gradient',-1)
         CALL PRIGRD(GRADKE)
         CALL HEADER('Nuclear attraction integral gradient',-1)
         CALL PRIGRD(GRADNA)
         CALL HEADER('Reorthonormalization gradient',-1)
         CALL PRIGRD(GRADFS)
         IF (SOLVNT) THEN
            CALL HEADER('Solvent energy contributions',-1)
            WRITE(LUPRI,'(3(/A,F15.10))')
     &         ' Nuclear solvent energy:',ESOLNN,
     &         ' Electr. solvent energy:',ESOLTT,
     &         ' Total   solvent energy:',ESOLNN+ESOLTT
            CALL HEADER('Electronic solvent part of gradient',-1)
            CALL PRIGRD(GSOLTT)
            CALL HEADER('Nuclear solvent part of gradient',-1)
            CALL PRIGRD(GSOLNN)
         END IF
#endif
         IF (SECDER) THEN
            CALL HEADER('Kinetic energy integral Hessian',-1)
            CALL PRIHES(HESSKE,'CENTERS')
            CALL HEADER('Nuclear attraction integral Hessian',-1)
            CALL PRIHES(HESSNA,'CENTERS')
            CALL HEADER('Highest order reorthonormalization Hessian',-1)
            CALL PRIHES(HESFS2,'CENTERS')
            IF (SOLVNT) THEN
               CALL HEADER('HSOLT2 part of solvent Hessian',-1)
               CALL PRIHES(HSOLT2,'CENTERS')
            END IF
         END IF
         IF (DIFDIP) THEN
            CALL HEADER('Electronic contributions to static part'//
     &                  ' of dipole gradient',-1)
            CALL FCPRI(DDIPE,'APT')
         END IF
      END IF
C
C     ********************************************************
C     ***** Sort gradient elements on direct access unit *****
C     ********************************************************
C
      IF (PROPTY) THEN
         IF (SECDER .OR. DIFDIP) THEN
            IF (SOLVNT) THEN
               CALL TLMWRT(TLMD,WORK(KFCM),WORK(KLAST),LWRK,IPRINT)
            END IF
            NMATS = 3*NUCDEP*(MAXREP+1)
            CALL SORONE(WORK,LWORK,NMATS,INDMAX,IPRINT)
            CLOSE(LUTEMP, STATUS='DELETE')
            IF (IPRINT .GE. 5) CALL SHDPRI(WORK,LWORK)
            IF (DOREPS(0)) CALL SHDCHK(WORK,LWORK,NODC,NODV,IPRINT)
         END IF
      END IF
      RETURN
 1000 FORMAT (//,2X,'***************************************',
     *         /,2X,'********** ISHELA/B =',I3,',',I3,' **********',
     *         /,2X,'***************************************',/)
      END
C  /* Deck dsofso */
      SUBROUTINE DSOFSO(DSO,FSO,WORK,LWORK,IPRINT,NODC,NODV)
#include "implicit.h"
#include "inforb.h"
!      use inforb_mod, only: N2ORBT, NCMOT, NNASHX, NNBAST
      LOGICAL NODC, NODV
      DIMENSION DSO(NNBAST), FSO(NNBAST), WORK(LWORK)
C
      KCMO  = 1
      KFT   = KCMO + NCMOT
      KDV   = KFT  + N2ORBT
      KLAST = KDV  + NNASHX
      IF (KLAST .GT. LWORK) CALL STOPIT('DSOFSO',' ',KLAST,LWORK)
      CALL DSOFS1(WORK(KCMO),WORK(KFT),WORK(KDV),DSO,FSO,IPRINT,NODC,
     &            NODV)
      RETURN
      END
C  /* Deck dsofs1 */
      SUBROUTINE DSOFS1(CMO,FT,DV,DSO,FSO,IPRINT,NODC,NODV)
C
C     This subroutine calculates the folded total one-electron density
C     and Fock matrices in AO basis (contravariant).  Input is
C     one-electron active density matrix and total Fock matrix in
C     MO basis.
C
C     Symmetry included        880418  PRT & TUH
C
      use priunit_mod, only: LUPRI
!      use abainf_mod, only: CCSD
!      use inforb_mod, only: N2ORBT, NASHT, NBAST, NCMOT, NISHT, NNASHX, NNASHY, NNORBT, NOCCT, NORBT, NSSHT, NSYM
      use iratdef_mod, only: IRAT
      use nuctap_mod, only: LUMC
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D0, TWO = 2.0D0, HALF = 0.5D0
      LOGICAL NODC, NODV
      DIMENSION CMO(*), FT(*), DV(*), DSO(NNBAST), FSO(NNBAST)
C
C Used from common blocks:
C   ABAINF : CCSD
C
#include "abainf.h"
#include "inforb.h"
      INTEGER R, S, RS, U, V, UV
C
C
C     ***** Read input from LUMC *****
C
      REWIND LUMC
      IF (CCSD) THEN
         CALL MOLLAB('ABACCIPH',LUMC,LUPRI)
         CALL MOLLAB('CANORB  ',LUMC,LUPRI)
         CALL READSQ(LUMC,IRAT*NCMOT,CMO)
         CALL MOLLAB('CCEFF   ',LUMC,LUPRI)
         CALL READSQ(LUMC,IRAT*NNASHX,DV)
         CALL READSQ(LUMC,IRAT*N2ORBT,FT)
      ELSE
         CALL MOLLAB('SIR IPH ',LUMC,LUPRI)
         READ (LUMC)
         READ (LUMC) NISHT,NASHT,NOCCT,NORBT,NBAST,NCONF,NWOPT,NWOPH,
     *               NCMOT,NNASHX,NNASHY,NNORBT,N2ORBT
         NSSHT  = NORBT - NOCCT
         CALL READSQ(LUMC,IRAT*NCMOT,CMO)
         IF (NASHT .GT. 0) THEN
            CALL READSQ(LUMC,IRAT*NNASHX,DV)
         ELSE
            READ (LUMC)
         END IF
         CALL READSQ(LUMC,IRAT*N2ORBT,FT)
      END IF
C
C     ***** Print Section *****
C
      IF (IPRINT .GT. 05) THEN
         WRITE (LUPRI, '(//A/)') ' <<<<< SUBROUTINE DSOFS1 >>>>>'
         WRITE (LUPRI, '(A,8I5)') ' NISH ', (NISH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NASH ', (NASH(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NOCC ', (NOCC(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NORB ', (NORB(I),I = 1,NSYM)
         WRITE (LUPRI, '(A,8I5)') ' NBAS ', (NBAS(I),I = 1,NSYM)
         IF (IPRINT .GE. 10) THEN
            CALL HEADER('Occupied molecular orbitals',0)
            IEND = 0
            DO 1000 ISYM = 1,NSYM
               IF (NBAS(ISYM) .EQ. 0) GOTO 1000
               IF (NOCC(ISYM) .EQ. 0) GOTO 1100
               WRITE (LUPRI, '(//,A,I5,/)') ' Symmetry ', ISYM
               IENDI = IEND
               DO 1200 I = 1, NOCC(ISYM)
                  WRITE (LUPRI,'(/,A,I5,/)') ' Molecular orbital ', I
                  WRITE (LUPRI,'(6F12.6)') (CMO(IENDI+J),J=1,NBAS(ISYM))
                  IENDI = IENDI + NBAS(ISYM)
1200           CONTINUE
               CALL HEADER('Total Fock matrix (MO basis)',-1)
               CALL OUTPUT(FT(I2ORB(ISYM)+1),1,NORB(ISYM),1,NORB(ISYM),
     &                     NORB(ISYM),NORB(ISYM),1,LUPRI)
1100           CONTINUE
               IEND = IEND + NORB(ISYM)*NBAS(ISYM)
1000        CONTINUE
            CALL HEADER('Active density matrix (MO basis)',-1)
            CALL OUTPAK(DV,NASHT,1,LUPRI)
         END IF
      END IF
C
C     ***** Construct contravariant SO matrices *****
C
      ISEND = 0
      ICEND = 0
      IFEND = 0
      DO 110 ISYM = 1,NSYM
         NORBI = NORB(ISYM)
         NISHI = NISH(ISYM)
         NASHI = NASH(ISYM)
         IASHI = IASH(ISYM)
         NBASI = NBAS(ISYM)
         IF (NBASI .EQ. 0) GOTO 120
         IF (NOCC(ISYM) .EQ. 0) THEN
            CALL DZERO(DSO(ISEND+1),NNBAS(ISYM))
            CALL DZERO(FSO(ISEND+1),NNBAS(ISYM))
            GO TO 120
         END IF
         RS = 0
         DO 100 R = 1, NBASI
            DO 200 S = 1,R
               RS = RS + 1
C
C              ***************************************
C              ***** One-electron density matrix *****
C              ***************************************
C
               DTRS = D0
C
C              (I) Inactive contribution
C
               IF (NISHI .GT. 0) THEN
                  ICENDI = ICEND
                  DO 300 I = 1, NISHI
                     DTRS = DTRS + CMO(ICENDI+R)*CMO(ICENDI+S)
                     ICENDI = ICENDI + NBASI
  300             CONTINUE
                  DTRS = DTRS + DTRS
               END IF
               IF (NODC) DTRS = D0
C
C              (II) Active contribution
C
               IF (.NOT. NODV) THEN
                  IF (NASHI .GT. 0) THEN
                     UV = ((IASHI + 1)*(IASHI + 2))/2
                     IDVEND = ICEND + NISHI*NBASI
                     ICENDU = IDVEND
                     DO 400 U = 1,NASHI
                        ICENDV = IDVEND
                        DO 410 V = 1, U
                           DUV = DV(UV)
                           IF (ABS(DUV) .GT. D0) THEN
                              TEMP = CMO(ICENDU+R)*CMO(ICENDV+S)
                              IF (U .NE. V) TEMP = TEMP
     *                             + CMO(ICENDU+S)*CMO(ICENDV+R)
                              DTRS = DTRS + DUV*TEMP
                           END IF
                           UV = UV + 1
                           ICENDV = ICENDV + NBASI
  410                   CONTINUE
                        UV = UV + IASHI
                        ICENDU = ICENDU + NBASI
  400                CONTINUE
                  END IF
               END IF
               IF (R .NE. S) DTRS = DTRS + DTRS
               DSO(ISEND+RS) = DTRS
C
C              ***********************
C              ***** Fock matrix *****
C              ***********************
C
               FTRS = D0
               IJ = 0
               ICENDI = 0
               DO 500 I = 1, NORBI
                  ICENDJ = 0
                  DO 510 J = 1, NOCC(ISYM)
                     IJ = IJ + 1
                     FTIJ = FT(IFEND+IJ)
                     IF (ABS(FTIJ) .GT. D0) THEN
                        TEMP = CMO(ICEND+ICENDI+R)*CMO(ICEND+ICENDJ+S)
                        IF (R .NE. S) TEMP = TEMP +
     *                         CMO(ICEND+ICENDI+S)*CMO(ICEND+ICENDJ+R)
                        FTRS = FTRS + FTIJ*TEMP
                     END IF
                     ICENDJ = ICENDJ + NBASI
  510             CONTINUE
                  IJ = IJ + NSSH(ISYM)
                  ICENDI = ICENDI + NBASI
  500          CONTINUE
               FSO(ISEND+RS) = FTRS
  200       CONTINUE
  100    CONTINUE
C
C        ***** Print Section *****
C
         IF (IPRINT .GE. 10) THEN
            WRITE (LUPRI,'(1X,A,I5)') ' Symmetry', ISYM
            CALL HEADER('Total density matrix (SO basis)',-1)
            CALL OUTPAK(DSO(ISEND+1),NBASI,1,LUPRI)
            CALL HEADER('Total Fock matrix (SO basis)',-1)
            CALL OUTPAK(FSO(ISEND+1),NBASI,1,LUPRI)
         END IF
120      CONTINUE
         IFEND = IFEND + NORBI*NORBI
         ISEND = ISEND + (NBASI*(NBASI + 1))/2
         ICEND = ICEND + NORBI*NBASI
110   CONTINUE
      RETURN
      END
C  /* Deck dipnuc */
      SUBROUTINE DIPNUC(IPRINT,DIFDIP)
C
C     Calculates nuclear contributions to electric dipole moments
C     and dipole gradients
C
C     1985 tuh
C     symmetry Dec 1988 tuh
C
      use maxaqn_mod, only: MXAQN, MXQN
      use symmet_mod, only: FMULT, IPTAX, IPTCNT, ISYMAX, MAXREP, NAXREP
      use nuclei_mod, only: CHARGE, CORD, ISTBNU, NUCDEP, NUCIND
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use maxorb_mod, only: MXCORB
      use dipole_mod, only: DDIPN, DIPMN
!      use dorps_mod, only: DOREPS
!      use orgcom_mod, only DIPORG
      use kinds_mod, only: dp
#include "implicit.h"
      REAL(dp), PARAMETER :: D0 = 0.0D0, D100 = 100.0D0
      LOGICAL DIFDIP
#include "orgcom.h"
#include "dorps.h"
C
C     Dipole moment
C
      CALL DZERO(DIPMN,3)
      DO 100 IATOM = 1, NUCIND
         CHA = CHARGE(IATOM)
CTROND         IF (ABS(CHA) .LT. D100) THEN
            FAC = FMULT(ISTBNU(IATOM))*CHA
            DO 200 ICOOR = 1, 3
               IF (ISYMAX(ICOOR,1) .EQ. 0) THEN
                  DIPMN(ICOOR) = DIPMN(ICOOR) 
     &                         + FAC*(CORD(ICOOR,IATOM)- DIPORG(ICOOR))
               END IF
  200       CONTINUE
CTROND         END IF
  100 CONTINUE
C
C     Dipole moment gradient
C
      IF (DIFDIP) THEN
         CALL DZERO(DDIPN,9*NUCDEP)
         DO 300 IREP = 0, MAXREP
            IF (DOREPS(IREP) .AND. (NAXREP(IREP,1).GT.0)) THEN
               DO 400 IATOM = 1, NUCIND
                  CHA = CHARGE(IATOM)
                  IF (ABS(CHA) .LT. D100) THEN
                     FAC = FMULT(ISTBNU(IATOM))*CHA
                     DO 500 ICOOR = 1, 3
                        IF (ISYMAX(ICOOR,1) .EQ. IREP) THEN
                           ISC = IPTCNT(3*(IATOM - 1) + ICOOR,IREP,1)
                           IF (ISC.GT.0) DDIPN(IPTAX(ICOOR,1),ISC) = FAC
                        END IF
  500                CONTINUE
                  END IF
  400          CONTINUE
            END IF
  300    CONTINUE
      END IF
      IF (IPRINT .GT. 0) THEN
         CALL HEADER('Nuclear contribution to dipole moments',-1)
         CALL DP0PRI(DIPMN)
         IF (DIFDIP) THEN
            CALL HEADER('Nuclear contribution to dipole gradient',-1)
            CALL FCPRI(DDIPN,'APT')
         END IF
      END IF
      RETURN
      END
C  /* Deck tlmwrt */
      SUBROUTINE TLMWRT(TLMD,FCM,WORK,LWORK,IPRINT)
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use nuclei_mod, only: NUCDEP
#include "implicit.h"
      DIMENSION TLMD(*), FCM(*), WORK(LWORK)
      KGRADE = 1
      KGRADN = KGRADE + 3*NUCDEP
      KGDIFE = KGRADN + 3*NUCDEP
      KGDIFN = KGDIFE + 3*NUCDEP
      KLAST  = KGDIFN + 3*NUCDEP
      IF (KLAST .GT. LWORK) CALL STOPIT('TLMWRT',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      CALL TLMWR1(TLMD,FCM,WORK(KGRADE),WORK(KGRADN),WORK(KGDIFE),
     &            WORK(KGDIFN),WORK(KLAST),LWRK,IPRINT)
      RETURN
      END
C  /* Deck tlmwr1 */
      SUBROUTINE TLMWR1(TLMD,FCM,GRADE,GRADN,GDIFE,GDIFN,WORK,LWORK,
     &                  IPRINT)
      use nuctap_mod, only: ABATLM, LRTLM, LUTLM, NBTLM
      use taysol_mod, only: GSOLNN, GSOLTT
      use nuclei_mod, only: NUCDEP
      use iratdef_mod, only: IRAT
      use mxcent_mod, only: MXATOM, MXBSETS, MXCENT, MXCOOR
      use priunit_mod, only: LUPRI
#include "implicit.h"
      PARAMETER (DM1 = -1.0D00, THRSH = 1.0D-07)
      LOGICAL ERRORE, ERRORN
      DIMENSION TLMD(LMNTOT,3*NUCDEP,2), FCM(LMNTOT), GRADE(*),
     &          GRADN(*), GDIFE(*), GDIFN(*), WORK(LWORK)
#include "cbisol.h"
C
C     Test TLMD
C
      IF (IPRINT .GT. 6) THEN
         CALL HEADER('Electronic TLMED in Cart. basis',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: lmn and ICOOR'
         CALL OUTPUT(TLMD(1,1,1),1,LMNTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
         CALL HEADER('Nuclear TLMND in Cart. basis',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: lmn and ICOOR'
         CALL OUTPUT(TLMD(1,1,2),1,LMNTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
      END IF
      DO 100 I = 1, 3*NUCDEP
         GRADE(I) =   DDOT(LMNTOT,FCM,1,TLMD(1,I,1),1)
         GRADN(I) = - DDOT(LMNTOT,FCM,1,TLMD(1,I,2),1)
         GDIFE(I) = GSOLTT(I) - GRADE(I)
         GDIFN(I) = GSOLNN(I) - GRADN(I)
  100 CONTINUE
      ERRORE = DNORM2(3*NUCDEP,GDIFE,1) .GT. THRSH
      ERRORN = DNORM2(3*NUCDEP,GDIFN,1) .GT. THRSH
      IF (IPRINT .GT. 2 .OR. ERRORE) THEN
         IF (ERRORE) WRITE (LUPRI,'(A)')
     &     ' Error in electronic solvent part of gradient '
         CALL HEADER('Difference GSOLTT and GRADE in TLMSTS',-1)
         WRITE (LUPRI,'(11X,A)')
     &      'GSOLTT             GRADE               difference'
         WRITE (LUPRI,'((2X,3F20.10))')
     &         (GSOLTT(I), GRADE(I), GDIFE(I),I=1,3*NUCDEP)
      END IF
      IF (IPRINT .GT. 2 .OR. ERRORN) THEN
         IF (ERRORN) WRITE (LUPRI,'(A)')
     &     ' Error in nuclear solvent part of gradient '
         CALL HEADER('Difference GSOLNN and GRADN in TLMSTS',-1)
         WRITE (LUPRI,'(11X,A)')
     &      'GSOLNN             GRADN               difference'
         WRITE (LUPRI,'((2X,3F20.10))')
     &         (GSOLNN(I), GRADN(I), GDIFN(I),I=1,3*NUCDEP)
      END IF
C
C     Sum electronic and nuclear contributions
C     ========================================
C
      DO 200 I = 1, 3*NUCDEP
         CALL DAXPY(LMNTOT,DM1,TLMD(1,I,1),1,TLMD(1,I,2),1)
  200 CONTINUE
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('TLMD in Cart. basis',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: lmn and ICOOR'
         CALL OUTPUT(TLMD(1,1,2),1,LMNTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
      END IF
C
C     Transform TLMD to spherical basis
C     =================================
C
      CALL TLMTRA(3*NUCDEP,TLMD(1,1,2),TLMD(1,1,1),WORK,LWORK,IPRINT)
C
C     Write TLMD to disk
C     ==================
C
      IF (IPRINT .GT. 5) THEN
         CALL HEADER('TLMD in spher. basis as written on LUTLM',-1)
         WRITE(LUPRI,'(A)') ' Rows and columns: l,m and ICOOR'
         CALL OUTPUT(TLMD(1,1,1),1,LMTOT,1,3*NUCDEP,
     &               LMNTOT,3*NUCDEP,1,LUPRI)
      END IF
      CALL OPENDX(LUTLM,ABATLM,IRAT*LMTOT,3*NUCDEP,'UNKNOWN',
     &            LRTLM,NBTLM,OLDDX)
      DO 300 I = 1, 3*NUCDEP
         CALL WRITDX(LUTLM,LRTLM,I,IRAT*LMTOT,TLMD(1,I,1))
  300 CONTINUE
      CLOSE(LUTLM,STATUS='KEEP')
      RETURN
      END
C  /* Deck tlmtra */
      SUBROUTINE TLMTRA(NMAT,CARMAT,SPHMAT,WORK,LWORK,IPRINT)
#include "implicit.h"
      DIMENSION CARMAT(*), SPHMAT(*), WORK(LWORK)
#include "cbisol.h"
C
      MXLM  = 2*LCAVMX + 1
      MXXYZ = (LCAVMX+1)*(LCAVMX+2)/2
      KTRAMA = 1
      KWRK1  = KTRAMA + MXLM*MXXYZ
      IF (KWRK1 .GT. LWORK) CALL STOPIT('TLMTRA',' ',LWORK,KWRK1)
      LWRK1  = LWORK + 1 - KWRK1
      CALL TLMTR1(NMAT,CARMAT,SPHMAT,WORK(KTRAMA),MXXYZ,MXLM,
     &            WORK(KWRK1),LWRK1,IPRINT)
      RETURN
      END
C  /* Deck tlmtr1 */
      SUBROUTINE TLMTR1(NMAT,CARMAT,SPHMAT,TRAMAT,MXXYZ,MXLM,
     &                  WORK,LWORK,IPRINT)
#include "implicit.h"
      PARAMETER (D0 = 0.0D0)
      PARAMETER (MORDER = 1, MINTEG = 0)
      DIMENSION CARMAT(LMNTOT,*), SPHMAT(LMNTOT,*), TRAMAT(MXXYZ,MXLM),
     &          WORK(LWORK)
#include "cbisol.h"
C
      JLM  = 0
      JXYZ = 0
      DO 100 IO = 0, LCAVMX
         NLM  = 2*IO + 1
         NXYZ = (IO+1)*(IO+2)/2
         CALL SPHCOM(IO,TRAMAT,MXLM,MXXYZ,MORDER,MINTEG,WORK,LWORK,
     &               IPRINT)
         DO 200 ICOOR = 1, NMAT
            DO 300 I = 1,NLM
               SPHIJ = D0
               DO 400 J = 1,NXYZ
                  SPHIJ = SPHIJ + CARMAT(JXYZ+J,ICOOR) * TRAMAT(J,I)
  400          CONTINUE
               SPHMAT(JLM+I,ICOOR) = SPHIJ
  300       CONTINUE
  200    CONTINUE
         JLM  = JLM  + NLM
         JXYZ = JXYZ + NXYZ
  100 CONTINUE
      RETURN
      END
