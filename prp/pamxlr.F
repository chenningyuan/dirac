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
C FILE    : pamxlr.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRPXLR */
      SUBROUTINE PRPXLR(WORK,KLRF,KCNV,KLRO,KFREE,LFREE)
C*****************************************************************************
C
C
C     Calculate linear response values defined in /XCBXLR/
C
C     Linear response functions are defined as:
C
C        <<A,;B>>_w = -(E[1]_A)^{+}(E[2] - wS[2])^{-1}(E[1]_B)
C
C        E[1]_A   is gradient of perturbation operator A  (h|A)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C     
C     where h is the operator manifold. In RPA h is limited to single 
C     excitations.
C
C     The linear response functions are found by iterative techniques.
C     We write
C
C        (E[2] - wS[2])^{-1}(E[1]_B) = X
C
C     which can be rearranged to the linear response equation
C
C        (E[2] - wS[2])X = E[1]_B
C       
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the reduced equation
C
C        (RE[2] - wRS[2])a = RE[1]_B
C
C     in which
C
C       RE[2] = B^{+}E[2]B
C       RS[2] = B^{+}S[2]B
C       RE[1] = B^{+}E[1}
C    
C
C     The evaluation of linear reponse functions thus proceeds in
C     three steps:
C
C       1. Generate property gradient      : E[1]_B
C       2. Solve reduced response equation : (RE[2] - wRS[2])a = RE[1]_B
C       3. Generate full solution vector   : X
C
C
C     Each operator A and B is assumed to span a given boson irrep
C     and be either symmetric(+) or antisymmetric(-) with respect
C     to time reversal.
C
C     Written by Trond Saue and Hans Joergen Aa. Jensen July 1996
C     Last revision July 9 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcborb.h>
#include <dcbxlr.h>
#include <dgroup.h>
      DIMENSION WORK(*)
C
      CALL QENTER('PRPXLR')
C
      NRSPF  = NLRBPT*NBFREQ
      IF(NRSPF.EQ.0) GOTO 10
      IF(ORBXLR) THEN
        CALL MEMGET('INTE',KLRO,2*NORBT,WORK,KFREE,LFREE)
        CALL DEFORB(WORK(KLRO),NCLS,3,XLR_INDSTR,
     &              WORK,KFREE,LFREE,IPRXLR)
        NI = 0
C       Number of inactive orbitals
CTROND  Check out active ones...
        DO I = 1,NFSYM
          NI = NI + NCLS(1,1,I) + NCLS(2,1,I)
        ENDDO        
      ELSE
        NI = 0
        KLRO = KFREE
      ENDIF
      IF(SKIPEE) THEN
        NEE = 0
      ELSE
        NEE = NI + 1
      ENDIF
      IF(SKIPEP) THEN
        NEP = 0
      ELSE
        NEP = NI + 1
      ENDIF
      NCI = 0
      N2LRF  = NLRAPT*NRSPF*(1+NEE+NEP+NCI)
      CALL MEMGET('REAL',KLRF,N2LRF,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KLRF),N2LRF)
      CALL MEMGET('REAL',KCNV,NRSPF,WORK,KFREE,LFREE)
C
      CALL PRPXL1(WORK(KLRF),WORK(KCNV),WORK(KLRO),WORK,KFREE,LFREE)
C
 10   CONTINUE
      IF(XLRANA) CALL ANAXLR(WORK,KFREE,LFREE)
C
      CALL QEXIT('PRPXLR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRPXL1 */
      SUBROUTINE PRPXL1(ABLRF,RCNV,IORBCL,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate linear response values defined in /XCBXLR/
C
C     Linear response functions are defined as:
C
C        <<A,;B>>_w = -(E[1]_A)^{+}(E[2] - wS[2])^{-1}(E[1]_B)
C
C        E[1]_A   is gradient of perturbation operator A  (h|A)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C     
C     where h is the operator manifold. In RPA h is limited to single 
C     excitations.
C
C     The linear response functions are found by iterative techniques.
C     We write
C
C        (E[2] - wS[2])^{-1}(E[1]_B) = X
C
C     which can be rearranged to the linear response equation
C
C        (E[2] - wS[2])X = E[1]_B
C       
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the reduced equation
C
C        (RE[2] - wRS[2])a = RE[1]_B
C
C     in which
C
C       RE[2] = B^{+}E[2]B
C       RS[2] = B^{+}S[2]B
C       RE[1] = B^{+}E[1}
C    
C
C     The evaluation of linear reponse functions thus proceeds in
C     three steps:
C
C       1. Generate property gradient      : E[1]_B
C       2. Solve reduced response equation : (RE[2] - wRS[2])a = RE[1]_B
C       3. Generate full solution vector   : X
C
C
C     Each operator A and B is assumed to span a given boson irrep
C     and be either symmetric(+) or antisymmetric(-) with respect
C     to time reversal.
C
C     Written by Trond Saue and Hans Joergen Aa. Jensen July 1996
C     Last revision Nov 26 1997 - tsaue - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(GPTHRS = 1.0D-9,D0 = 0.0D0,D2=2.0D0)
#if defined (VAR_MPI)
#include <infmpi.h>
      INCLUDE 'mpif.h'
#endif
#if defined (VAR_PVM)
#include <fpvm3h.h>
#endif
#include <aovec.h>
#include <mxorb.h>
#include <infpar.h>
C
#include <dcbgen.h>
#include <dcbxlr.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <pgroup.h>
#include <dcbxpr.h>
#include <dcbxrs.h>
C
#include <blocks.h>
      CHARACTER TTYP(-1:1)*2
      DIMENSION ABLRF(NFREQ,NLRAPT,NLRBPT,*),RCNV(NBFREQ,NLRBPT),
     &          IORBCL(2,*),WORK(*)
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
C
      CALL QENTER('PRPXL1')
      KFRSAV = KFREE
C
C     Get hold of the slaves if we run in parallel.
C        ( ITASK = 1 for Fock matrices )
C
      IF (PARCAL) CALL DIRAC_PARCTL( 1 )
C
C     ****************************
C     *** I N I T I A L I Z E  ***
C     ****************************
C
      CALL SETXLR(WORK,KFREE,LFREE)
C
C     If orbital analysis; extract relevant orbital energies
C     ======================================================
C
      IF(ORBXLR) THEN
        NI = 0
C       Number of inactive orbitals
CTROND  Check out active ones...
        DO I = 1,NFSYM
          NI = NI + NCLS(1,1,I) + NCLS(2,1,I)
        ENDDO        
        CALL MEMGET('REAL',KEIG,NI,WORK,KFREE,LFREE)
        CALL EIGEXT(WORK(KEIG),IORBCL,1,WORK,KFREE,LFREE)
      ELSE
        KEIG = KFREE
      ENDIF
C
C     **********************************************************
C     ***                                                    ***
C     *** L I N E A R   R E S P O N S E   F U N C T I O N S  ***
C     ***                                                    ***
C     **********************************************************
C
C     Memory allocation
C     =================
C
C     Dimensions of reduced system
C
      KXLR = KFREE
      NEVECR    = NREDM*NFREQ
      CALL MEMGET('INTE',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBCVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBPVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVALR,NFREQ  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVECR,NEVECR ,WORK,KFREE,LFREE)
C
C     Transfer frequencies to EVALR
C     =============================
C
      CALL DCOPY(NFREQ,BFREQ,1,WORK(KEVALR),1)
C
C       Loop over operators
C       ===================
C
      DO 10 JOPSY = 1,NFSYM
        NFC = NBSYM/NFSYM
        DO 20 IS = 1,NFC
          JSYMOP  = JFSYM(IS,JOPSY)
          NOPB    = NLRBP(JSYMOP)
          JOPB    = JLRBP(JSYMOP)
C
          DO 30 IOPB = 1,NOPB
            INDBP  = LLBPSU(JOPB+IOPB)
            INDPRB = LLRBPU(INDBP)
            JTIMOP = IPRPTIM(INDPRB)
            WRITE (LUPRI,'(///A/A,A,A3,2X,A2/)')
     &        ' <<<  SOLVING SETS OF LINEAR EQUATIONS '//
     &        'FOR LINEAR RESPONSE PROPERTY >>>',
     &        PRPNAM(INDPRB),
     &        ' - Symmetry: ',REP(JSYMOP-1),TTYP(JTIMOP)
C
C           Initialize configurational and orbital parameters
C           =================================================
C
            CALL XLRPAR(WORK,KFREE,LFREE,IPRXLR)
C      
C           ******************************************
C           *** P R O P E R T Y    G R A D I E N T ***
C           ******************************************
C
            KGRAD = KFREE
            CALL MEMGET('REAL',KGPCI,NZCONFQ,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KGPOE,NZXOPEQ,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KGPOP,NZXOPPQ,WORK,KFREE,LFREE)
            CALL PAMGRD(INDPRB,WORK(KGPCI),WORK(KGPOE),WORK(KGPOP),
     &                  JOPSY,WORK(KZCONF),WORK(KZXOPE),WORK(KZXOPP),
     &                  NZCONF,NZXOPE,NZXOPP,WORK,KFREE,LFREE,IPRXLR)
C
C           Check norm of property gradient; skip if below threshold
C           ========================================================
C        
            GPCNRM = D2*DNRM2(NZCONFQ,WORK(KGPCI),1)
            GPENRM = D2*DNRM2(NZXOPEQ,WORK(KGPOE),1)
            GPPNRM = D2*DNRM2(NZXOPPQ,WORK(KGPOP),1)
            GPTNRM = SQRT(GPCNRM**2+GPENRM**2+GPPNRM**2)
            IF(GPTNRM.LT.GPTHRS) THEN
              WRITE(LUPRI,'(A,A16)') 
     &          'WARNING: Skipping B property ',PRPNAM(INDPRB)
              WRITE(LUPRI,'(9X,A,1P,D8.1)') 
     &          'Norm of gradient  :',GPTNRM,
     &          'is below threshold:',GPTHRS
              GOTO 30
            ENDIF
C
C           Compress number of orbital rotations
C           ====================================
C
            IF(THRCOM.NE.D0) THEN
              CALL REDPAR(WORK(KGPCI),WORK(KGPOE),WORK(KGPOP),
     &                    WORK(KZCONF),WORK(KZXOPE),WORK(KZXOPP),
     &                    THRCOM,WORK,KFREE,LFREE,IPRXLR)
            ENDIF
C
C           ********************************************************
C           *** S O L V E   R E S P O N S E    E Q U A T I O N S ***
C           ********************************************************
C
            NCRED  = 0
            NERED  = 0
            NPRED  = 0
            NZRED  = NCRED + NERED + NPRED
            RSPLAB = PRPNAM(INDPRB)
            CALL XRSCTL(WORK(KGPCI) ,WORK(KGPOE),WORK(KGPOP),
     &                 WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),
     &                 WORK(KIBPVC),RCNV(1,INDBP),
     &                 WORK(KEVALR),WORK(KEVECR),
     &                 WORK,KFREE,LFREE)
C          ... reuse KGP space for solution vectors.
            CALL MEMREL('PRPXL1.XRSCTL',WORK,1,KGRAD,KFREE,LFREE)
C
C           ************************************************************
C           *** E V A L U A T E   R E S P O N S E   F U N C T I O N  ***
C           ************************************************************
C
            CALL GETLRF(ABLRF,WORK(KIBTYP),WORK(KIBCVC),
     &                  WORK(KIBEVC),WORK(KIBPVC),RCNV(1,INDBP),
     &                  WORK(KEVALR),WORK(KEVECR),IORBCL,WORK(KEIG),
     &                  WORK,KFREE,LFREE)
 30       CONTINUE
 20     CONTINUE
 10   CONTINUE
      CALL MEMREL('PRPXL1',WORK,1,KXLR,KFREE,LFREE)
C
C     Final print section
C
      CALL XLROUT(ABLRF,RCNV,WORK(KEIG),WORK,KFREE,LFREE)
      CLOSE(LUXVC,STATUS = 'KEEP')
      CALL MEMREL('PRPXL1',WORK,1,KFRSAV,KFREE,LFREE)
C
C
C     Release the slaves if we run in parallel.
C        ( ITASK = -1 )
C
      IF (PARCAL) CALL DIRAC_PARCTL( -1 )
      CALL QEXIT('PRPXL1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xlrinp */
      SUBROUTINE XLRINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for linear response  <<A;B>>
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C     Last revision: June 14 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include <dummy.h>
#include <mxorb.h>
#include <mxcent.h>
      PARAMETER (NTABLE = 32,DMIL = 1.0D3)
C
#include <dcbgen.h>
#include <dcbprp.h>
#include <dcbxlr.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbham.h>
#include <dgroup.h>
#include <pgroup.h>
#include <dcbxpr.h>
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7,
     &          MXFORM*6,FMT*6
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA TABLE /'.PRINT ','.A OPER','.B OPER','.OPERAT',
     &            '.B FREQ','.ALLCMB','.THRESH','.MAXITR',
     &            '.SKIPEP','.MAXRED','.QED   ','.XLRNRM',
     &            '.ONLYSF','.ONLYSG','.INTFLG','.NOPREC',
     &            '.OCCUP ','.VIRTUA','.ACTIVE','.COMPRE',
     &            '.ITRINT','.CNVINT','.REAXVC','.TRIAB ',
     &            '.RESFAC','.ANALYZ','.STERNH','.SKIPEE',
     &            '.MXLOAD','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
#include <ibtfun.h>
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C     Local initialization
C
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
      INPERR = 0
C
C
C     Initialize /XCBLR/
C     ===================
C
      IPRXLR     = 0
      THCLR      = 1.0D-5
      THRCOM     = D0
      RESXLR     = DMIL
      CNVXLR(1)  = DUMMY
      CNVXLR(2)  = DUMMY
      MAXRM      = 200
      MXLOAD     = 0
      ITRXLR     = 50
      INTXLR     = ILLINT + 2*ISLINT + 4*ISSINT
      IF (LEVYLE) INTXLR = ILLINT
      ITRILR(1)  = 1
      ITRILR(2)  = 1
      NLRAPT     = 0
      NLRBPT     = 0
      CALL IZERO(NLRAP,8)
      CALL IZERO(NLRBP,8)
      NBFREQ     = 1
      BFREQ(1)   = D0
      ALLCMB     = .FALSE.
      SKIPEP     = .FALSE.
      SKIPEE     = .FALSE.
      QEDRUN     = .FALSE.
      LSFG(1)    = .TRUE.
      LSFG(2)    = .TRUE.
      XLRNRM     = .FALSE.
      XLRDIH     = .TRUE.
      TRIAB      = .FALSE.
      XVCFIL     = 'NORXVC'
      XLRANA     = .FALSE.
      XSTERN     = .FALSE.
      DO I = 1,NFSYM
         XLR_INDSTR(2,I) = ' '
         XLR_INDSTR(1,I) = ' '
         XLR_INDSTR(3,I) = ' '
      ENDDO
C
C     Initialize user-defined linear response properties
C     ==================================================
C
C     Electric dipole polarizibilities
C     Magnetic susceptibilities
      IF(MSUSCP) THEN    
C        CALL DEF_MSUSCP(WORK,LWORK)                  
      ENDIF
C     NMR shielding constants
      IF(SHIELD) THEN    
C        CALL DEF_SHIELD(WORK,LWORK)                  
      ENDIF
C     NMR indirect spin-spin coupling constants
C
C     Process input from CBIXLR
C     =========================
C
      NEWDEF = (WORD .EQ. '*LINEAR')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1, 2, 3, 4, 5, 6, 7, 8, 9,10,
     &                     11,12,13,14,15,16,17,18,19,20,
     &                     21,22,23,24,25,26,27,28,29,30,31,32), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in XLRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in XLRINP.')
    1          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRXLR
               GO TO 100
    2          CONTINUE
C&&&& A OPERATOR: Left side operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPRXLR)
                  CALL OP1IND('NLRAPT',IND1OP,LLRAPU,NLRAPT,INDXPR,
     &                        MAXLLR)
               GO TO 100
    3          CONTINUE
C&&&& B OPERATOR: Right side operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPRXLR)
                  CALL OP1IND('NLRBPT',IND1OP,LLRBPU,NLRBPT,INDXPR,
     &                        MAXLLR)
               GO TO 100
    4          CONTINUE
C&&&& OPERATOR: Both left/right side operator
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPRXLR)
                  CALL OP1IND('NLRAPT',IND1OP,LLRAPU,NLRAPT,INDXPR,
     &                        MAXLLR)
                  CALL OP1IND('NLRBPT',IND1OP,LLRBPU,NLRBPT,INDXPR,
     &                        MAXLLR)
               GO TO 100
    5          CONTINUE
C&&&& B FREQuency
                  READ (LUCMD, *) NBFREQ
                  IF (NBFREQ.GT.0) THEN
                    IF (NBFREQ.LE.MAXFLR) THEN
                       READ (LUCMD,*) (BFREQ(J),J=1,NBFREQ)
                    ELSE
                       INPERR = INPERR + 1
                       WRITE (LUPRI,'(3(/,A,I5),/)')
     &                 '@Number of B frequencies specified :',
     &                 NBFREQ,
     &                 '@is greater than the allowed number :',
     &                 MAXFLR
                       READ (LUCMD,*) (BFREQ(J),J=1,MAXFLR),
     &                              (FFFF,J=MAXFLR+1,NBFREQ)
                       NBFREQ = MAXFLR
                    END IF
                  END IF
               GO TO 100
    6          CONTINUE
C&&&& ALLCMB: form all possible <<A;B>>_w even if imaginary
               ALLCMB = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& THRESH - Threshold for convergence:
                  READ (LUCMD,*) THCLR
               GO TO 100
    8          CONTINUE
C&&&& MAXITR - Maximum number of iterations
                  READ (LUCMD, '(I5)') ITRXLR
               GO TO 100
    9          CONTINUE
C&&&& SKIPEP -exclude all electron-positron rotations
                  SKIPEP = .TRUE.
               GOTO 100
 10            CONTINUE
C&&&& MAXRED: max dimension of matrix in reduced system
                  READ(LUCMD,*) MAXRM
               GOTO 100
 11            CONTINUE
C&&&& QED  : allow rotations between positrons and secondary orbitals
                  QEDRUN = .TRUE.
               GOTO 100
 12            CONTINUE
C&&&& XLRNRM: normalize trial vectors
                 XLRNRM = .TRUE.
               GOTO 100
C&&&& ONLYSF: debug keyword: only call FMOLI in sigmavector routine
 13            CONTINUE
                  LSFG(1) = .FALSE.
               GOTO 100
C&&&& ONLYSG: debug keyword: only call GMOLI in sigmavector routine
 14            CONTINUE
                  LSFG(2) = .FALSE.
               GOTO 100
C&&&& INTFLG: specify what two-electron integrals to include
 15            CONTINUE
                  READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  INTXLR = ILLINT + 2*ISLINT + 4*ISSINT
               GOTO 100
 16            CONTINUE
C&&&& NOPREC: no preconditioning of initial trial vectors
                  XLRDIH = .FALSE.
               GOTO 100
 17            CONTINUE
C&&&& OCCUP: String of inactive orbitals in XLR module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XLR_INDSTR(1,I)
                  ENDDO
               GO TO 100
 18            CONTINUE
C&&&& VIRTUA: String of secondary orbitals in XLR module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XLR_INDSTR(3,I)
                  ENDDO
               GO TO 100
 19            CONTINUE
C&&&& ACTIVE: String of active orbitals in XLR module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XLR_INDSTR(2,I)
                  ENDDO
               GO TO 100
 20            CONTINUE
C&&&& COMPRESS: Compress number of orbital variations by checking 
C               gradient elements againt a threshold
                  READ(LUCMD,*) THRCOM
               GO TO 100
 21            CONTINUE
C&&& ITRINT: Number of iterations before adding SL- and SS-integrals
                  READ(LUCMD,*) ITRILR(1),ITRILR(2)
               GO TO 100
 22            CONTINUE
C&&& CNVINT: Convergence thresholds for adding SL- and SS-integrals
                  READ(LUCMD,*) CNVXLR(1),CNVXLR(2)
               GO TO 100
 23            CONTINUE
C&&&& REAXVC: Read solution vectors from file XVCFIL
                  READ(LUCMD,'(A6)') XVCFIL
               GO TO 100
 24            CONTINUE
C&&&& TRIAB:   Enforce triangularity of <<A;B>>
                  TRIAB = .TRUE.
               GO TO 100
 25            CONTINUE
C&&&& RESFAC:  New trial vectors only for parameters whose residual
C              has a norm that is a fraction RESFAC of the max norm
                  READ(LUCMD,*) RESXLR
               GO TO 100
 26            CONTINUE
C&&&& ANALYZe: Analyze individual contributions to repsonse functions
                  XLRANA = .TRUE.
               GO TO 100
 27            CONTINUE
C&&&& STERNH: Diagonal part of ep,ep-Hessian approximated by -2mc^2
                  XSTERN = .TRUE.
               GO TO 100
 28            CONTINUE
C&&&& SKIPEE: Skip all e-e rotations
                  SKIPEE = .TRUE.
               GO TO 100
 29            CONTINUE
C&&&& MALOAD: Maximum number of trial vectors to be loaded in each microiteration
                  READ(LUCMD,*) MXLOAD
               GO TO 100
 30            CONTINUE
               GO TO 100
 31            CONTINUE
               GO TO 100
 32            CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in XLRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in XLRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C     ===============
C
      IF(CNVXLR(1).LT.DUMMY) ITRILR(1) = 1
      IF(CNVXLR(2).LT.DUMMY) ITRILR(2) = 1
      MXREDZ = MAXRM/2
      MAXRM  = 2*MXREDZ
      IF(QEDRUN) SKIPEP = .FALSE.      
C
C     Print section
C     =============
C      
      DOXLR = DOXLR.OR.(NLRAPT.GT.0.AND.NLRBPT.GT.0.AND.NBFREQ.GT.0)
      IF(.NOT.(DOXLR.OR.XLRANA)) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &    'XLRINP: Electronic linear reponse properties'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I5)') 
     &   '* Print level                 :',IPRXLR
      IF(XLRANA) THEN
        WRITE(LUPRI,'(1X,A)') '* Analysis of solution vectors.'
      ENDIF
      IF(SKIPEE) THEN
        WRITE(LUPRI,'(1X,A)') 
     &     '* NB!!! No e-e rotations'
      ENDIF
      IF(SKIPEP) THEN
        WRITE(LUPRI,'(1X,A)') 
     &     '* NB!!! No rotations between electrons and positrons.'
      ENDIF
      IF(QEDRUN) THEN
        WRITE(LUPRI,'(1X,A)') 
     &     '* Rotations between positrons and secondary orbitals.'
      ENDIF
      IF(.NOT.DOXLR) THEN
        DOXLR = DOXLR.OR.XLRANA
        GOTO 999
      ENDIF
      IF(XVCFIL.NE.'NORXVC') THEN
        WRITE(LUPRI,'(1X,A,A,A6)')
     &    '* The program will attemp restart on solution vectors ',
     &    'on file ',XVCFIL
      ENDIF
      WRITE(LUPRI,'(1X,A,I5)') 
     &   '* Maximum number of iterations:',ITRXLR
      WRITE(LUPRI,'(1X,A,D12.6)') 
     &   '* Threshold for convergence   :',THCLR
      WRITE(LUPRI,'(1X,A,I5)') 
     &   '* Max. size of reduced matrix :',MAXRM
      IF(TRIAB) THEN
        WRITE(LUPRI,'(4X,A)') '* Assuming <<A;B>> to be symmetric.'
      ENDIF
      IF(XSTERN) THEN
        WRITE(LUPRI,'(4X,A)') 
     &    '* Approximating diagonal Hessian for ep-rotations by -2mc^2'
      ENDIF
      IF(THRCOM.NE.D0) THEN
        WRITE(LUPRI,'(1X,A/1X,A/1X,A,D12.6)')
     &  '* Reducing orbital parameters by deleting rotations for which',
     &  '  for which the corresponding property gradient is below',
     &  '  threshold ',THRCOM
      ENDIF
C
      IF(ALLCMB) THEN
        WRITE(LUPRI,'(1X,A)') 
     &  '* All combinations of A and B operators will be formed,',
     &  '   even if not real.'
      ENDIF
      WRITE(LUPRI,'(1X,A)') '* Include contributions from '//
     +    'the following two-electron integrals:'
      IF(ILLINT.EQ.1) WRITE(LUPRI,'(6X,A)') '- LL-integrals'
      IF(ISLINT.EQ.1) THEN
        IF(CNVXLR(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(6X,A,E8.1)')
     +    '- SL-integrals  below convergence ',CNVXLR(1)
        ELSE
          WRITE(LUPRI,'(6X,A,I4)')
     +    '- SL-integrals from iteration ',ITRILR(1)
        ENDIF
      ENDIF
      IF(ISSINT.EQ.1) THEN
        IF(CNVXLR(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,E8.1)')
     +    '- SS-integrals  below convergence ',CNVXLR(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    '- SS-integrals from iteration ',ITRILR(2)
        ENDIF
      ENDIF
      IF(.NOT.LSFG(1)) WRITE(LUPRI,'(1X,A)')
     &     '* Only calling FMOLI in sigma-vector routine...'
      IF(.NOT.LSFG(2)) WRITE(LUPRI,'(1X,A)')
     &     '* Only calling GMOLI in sigma-vector routine...'
      IF(XLRNRM) WRITE(LUPRI,'(1X,A)')
     &     '* Trial vectors will be normalized.'
      IF(XLRDIH) THEN
        WRITE(LUPRI,'(1X,A)')
     &  '* Initial trial vectors will be preconditioned.'
      ELSE
        WRITE(LUPRI,'(1X,A)')
     &  '* No preconditioning of initial trial vectors.'
      ENDIF
      IF(MXLOAD.NE.0) THEN
        WRITE(LUPRI,'(1X,A,I5)') 
     &  '* Maximum number of trial loaded in each microtieration:',
     &  MXLOAD
      ENDIF
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)') '** A operators **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NLRAPT
        INDXPR = LLRAPU(I)
        CALL WRIXPR(I,INDXPR)
      ENDDO
C
C     B operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)') '** B operators **'
      CALL PRSYMB(LUPRI,'-',75,0)
      DO I = 1,NLRBPT
        INDXPR = LLRBPU(I)
        CALL WRIXPR(I,INDXPR)
      ENDDO
C
C     Frequencies
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(1X,A)') 
     &'* Linear response functions evaluated at the following '//
     &   'frequencies:'
      DO J = 1,NBFREQ
        FMT = MXFORM(BFREQ(J),12)
        WRITE(LUPRI,'(6X,'//FMT//',A)') BFREQ(J),' a.u.'
      ENDDO        
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck NMRINP */
      SUBROUTINE NMRINP(WORD,RESET)
C*****************************************************************************
C
C     Input section for NMR properties
C
C     Written by T.Enevoldsen - Nov 1997
C     Last revision: 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <orgcom.h>
      PARAMETER (NTABLE = 8)
      LOGICAL SET, NEWDEF, RESET, GORG
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.GAUGEO','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF
C
C
C     Initialize 
C     ===================
C
      GORG = .FALSE.
      CALL DZERO(GAGORG,3)
C
C     Process input 
C     =========================
C
      NEWDEF = (WORD .EQ. '*NMR')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
  110       CONTINUE
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4,5,6,7,8), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in NMRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in NMRINP.')
    1          CONTINUE
C Read user specified gauge origin
                  READ(LUCMD,*) (GAGORG(ICOOR), ICOOR = 1,3)
                  GORG = .TRUE.
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
               GO TO 100
    4          CONTINUE
               GO TO 100
    5          CONTINUE
               GO TO 100
    6          CONTINUE
               GO TO 100
    7          CONTINUE
               GO TO 100
    8          CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in NMRINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in NMRINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Print section
C     =============
C
      IF (GORG) THEN
         CALL PRSYMB(LUPRI,'=',75,0)
         WRITE(LUPRI,'(1X,A)') 'NMRINP: '
         CALL PRSYMB(LUPRI,'=',75,0)
         IF (GORG) THEN
            WRITE(LUPRI,'(/,A,3F12.6,/)') 'User defined gauge '//
     &           'origin used:',(GAGORG(I), I = 1, 3)
         END IF
         CALL DCOPY(3,GAGORG,1,DIPORG,1)
      END IF
  999 CONTINUE
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxlr */
      SUBROUTINE SETXLR(WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Linear response module: 
C       Initialize variables/Open necessary files etc.
C
C     Written by T.Saue and H.J.Aa.Jensen July 1996
C     Last revision: July 9 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0)
C
#include <dgroup.h>
#include <dcbxlr.h>
#include <dcbgen.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbxrs.h>
C
      LOGICAL TOBE,ORBDEF(3)
      DIMENSION WORK(*),IR(2,2),NR(2)
      CALL QENTER('SETXLR')
      KFRSAV = KFREE
C
C     ************************************
C     *** Set orbital strings          ***
C     ************************************
C
      CALL OSTRING(WORK,KFREE,LFREE)
C
C     ************************************
C     *** Sort operators on symmetries ***
C     ************************************
C
      CALL OP1SRT('LLAPSU',LLRAPU,LLAPSU,NLRAPT,NLRAP,JLRAP,
     &            WORK(KFREE),LFREE)
      CALL OP1SRT('LLBPSU',LLRBPU,LLBPSU,NLRBPT,NLRBP,JLRBP,
     &            WORK(KFREE),LFREE)
C
C     **************************************************
C     *** Transfer information from XCBXLR to WRKXRS ***
C     **************************************************
C
      LINEQ     = .TRUE.
      LSVCFG(1) = LSFG(1)
      LSVCFG(2) = LSFG(2)
      TKNORM    = XLRNRM
      DIAGHE    = XLRDIH
      IPRXRS    = IPRXLR
      THCXRS    = THCLR
      RESFAC    = RESXLR
      MAXITR    = ITRXLR
      NFREQ     = NBFREQ
      MAXSIM    = MXLOAD
      NREDM     = MAXRM
      N2REDM    = MAXRM*MAXRM
      LOFFTY    = 0
      STATIC    = (NFREQ.EQ.1.AND.BFREQ(1).EQ.D0)
      CNVINT(1) = CNVXLR(1)
      CNVINT(2) = CNVXLR(2)
      ITRINT(1) = ITRILR(1)
      ITRINT(2) = ITRILR(2)
      INTDEF    = INTXLR
      STERNH    = XSTERN
C
C     Indicate from vectors are read...
C
      IOPT  = 0
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','EIGEXT')
      CALL REACMO(LUBUF,'DFCOEF',DUM,DUM,DUM,DUM,IOPT,
     &            DUM,IDUM,IDUM)
      CLOSE(LUBUF,STATUS='KEEP')
C
C     Initialize file for solution vectors
C
      LUXVC = 29
      CALL OPNFIL(LUXVC,'PAMXVC','UNKNOWN','PRPXLR')      
      WRITE(LUXVC) 'END_OF_THIS_FILE'
C
C     Look for restart file for solution vectors
C
      RSREST = .FALSE.
      IF(XVCFIL.NE.'NORXVC') THEN
        INQUIRE(FILE=XVCFIL,EXIST=TOBE)
        RSREST = TOBE
      ENDIF
      IF(RSREST) THEN
        LURST = 28
        CALL OPNFIL(LURST,XVCFIL,'OLD','PRPXLR')      
      ENDIF
      CALL QEXIT('SETXLR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETLRF */
      SUBROUTINE GETLRF(ABLRF,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &                  EVALR,EVECR,IORBCL,EIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Evaluate linear response functions <<A;B>>.
C
C     Written by T.Saue Oct 22 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
      PARAMETER(D1 = 1.0D0)
C
#include <dcbxrs.h>
#include <dcbxlr.h>
#include <dcbibn.h>
#include <dcborb.h>
#include <dgroup.h>
      LOGICAL FIRST
      DIMENSION ABLRF(NFREQ,NLRAPT,NLRBPT,*),RCNV(*),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IORBCL(2,*),
     &          IBPVC(*),EVALR(*),EVECR(NREDM,*),EIG(*),WORK(*)
C
      CALL QENTER('GETLRF')
      KFRSAV = KFREE
      FIRST = .TRUE.
      IF(ORBXLR) THEN
        NI = 0
        DO I = 1,NFSYM
          NI = NI + NCLS(1,1,I) + NCLS(2,1,I)
        ENDDO        
      ELSE
        NI = 0
      ENDIF
      NRSP = NFREQ*NLRAPT*NLRBPT
C
C     First evaluate solution vectors for all frequencies
C     and write them to direct access file
C     ===================================================
C
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
      IFAC  = 1
      IF(.NOT.STATIC) IFAC = 2
C
C     Get coefficients
C     ----------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      IOPT  = 4
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','GETLRF')
      CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),DUM,DUM,TOTERG,IOPT,
     &            WORK,KFREE,LFREE)
      CLOSE(LUBUF,STATUS='KEEP')
C
      IOFF  = 1
C
C     Solution vector - orbital (e-e) part
C     ====================================   
C        
      IF(NZXOPE.GT.0) THEN
        NESIM = NFREQ * IFAC
        KEE   = KFREE
        CALL MEMGET('REAL',KGVEC,NZXOPEQ      ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KXVEC,NZXOPEQ*NFREQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBBUF,NZXOPEQ      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIVEC,NESIM        ,WORK,KFREE,LFREE)
        CALL GETLR1(ABLRF(1,1,1,IOFF+1),JBENDX,EVECR,EVALR,RCNV,
     &              IBTYP,IBEVC,
     &              WORK(KGVEC),WORK(KXVEC),WORK(KBBUF),NZXOPE,
     &              WORK(KCMO),IORBCL,EIG,WORK(KZXOPE),WORK(KIVEC),
     &              WORK,KFREE,LFREE)
        IF(FIRST) THEN
          CALL DCOPY(NRSP,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
          FIRST = .FALSE.
        ELSE
          CALL DAXPY(NRSP,D1,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
        ENDIF
        CALL MEMREL('GETLRF.e-e',WORK,1,KEE,KFREE,LFREE)
        IOFF = IOFF + NI + 1
      ENDIF
C
C     Solution vector - orbital (e-p) part
C     ====================================   
C        
      IF(NZXOPP.GT.0) THEN
        NPSIM = NFREQ * IFAC
        KEP   = KFREE
        CALL MEMGET('REAL',KXVEC,NZXOPPQ*NFREQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KGVEC,NZXOPPQ      ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBBUF,NZXOPPQ      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIVEC,NPSIM        ,WORK,KFREE,LFREE)
        CALL GETLR1(ABLRF(1,1,1,IOFF+1),JBPNDX,EVECR,EVALR,RCNV,
     &              IBTYP,IBPVC,
     &              WORK(KGVEC),WORK(KXVEC),WORK(KBBUF),NZXOPP,
     &              WORK(KCMO),IORBCL,EIG,WORK(KZXOPP),WORK(KIVEC),
     &              WORK,KFREE,LFREE)
        IF(FIRST) THEN
          CALL DCOPY(NRSP,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
          FIRST = .FALSE.
        ELSE
          CALL DAXPY(NRSP,D1,ABLRF(1,1,1,IOFF+1),1,ABLRF,1)
        ENDIF
        CALL MEMREL('GETLRF.e-p',WORK,1,KEP,KFREE,LFREE)
        IOFF = IOFF + NI + 1
      ENDIF
C
C     Solution vector - configurational part
C     ======================================
C        
      IF(NZCONF.GT.0) THEN
C     ..not yet written.....tsaue
      ENDIF
C
      CALL MEMREL('GETLRF.end',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GETLRF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETLR1 */
      SUBROUTINE GETLR1(ABLRF,NBTYP,EVECR,EVALR,RCNV,IBTYP,IBVEC,
     &                  GPO,XPO,BBUF,NXPAR,
     &                  CMO,IORBCL,EIG,JXOP,IVECS,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Get linear response function for orbital rotation parameter
C       contribution from given variational parameter type
C       specified by NBTYP
C
C     Written by T.Saue Oct 22 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D4 = 4.0D0,D0 = 0.0D0,DM4 = -4.0D0, D2 = 2.0D00)
C
#include <dcbibn.h>
#include <dcbxrs.h>
#include <dcbxlr.h>
#include <dgroup.h>
#include <dcbxpr.h>
      CHARACTER FILNAM*6,MEMTYP*6,BLANK16*16,BLANK2*2
      DIMENSION ABLRF(NFREQ,NLRAPT,NLRBPT,*),RCNV(*),
     &          GPO(NXPAR,NZ),XPO(NXPAR,NZ,NFREQ),BBUF(*),
     &          CMO(*),EIG(*),IORBCL(2,*),
     &          IBVEC(*),IBTYP(2,*),EVECR(NREDM,*),EVALR(*),
     &          JXOP(2,*),IVECS(*),WORK(*)
      DATA BLANK16/'                  '/
      DATA BLANK2/'  '/
C
#include <chrnos.h>
      CALL QENTER('GETLR1')
      KFRSAV = KFREE
      IF(IPRXLR.GE.2) CALL HEADER('Output from GETLR1',-1)
C
      NBDIM = NXPAR*NZ
      IF(NBTYP.EQ.JBENDX) THEN
        NBRED  = NERED
        LUBVEC = LUBOE
        FILNAM = 'PAMBOE'
        MEMTYP = 'GETLRE'
        ITYP   = IPEP
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBRED  = NPRED
        LUBVEC = LUBOP
        FILNAM = 'PAMBOP'
        MEMTYP = 'GETLRP'
        ITYP   = IPPP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'GETLR1 ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('GETLR1: Unknown NBTYP !')
      ENDIF
C
      INDPRB = LLRBPU(INDBP)
      IF(INDPRB.GT.0) THEN
        IF(INDPRB.LT.10) THEN
          NDIG = 1
        ELSEIF(INDPRB.LT.100) THEN
          NDIG = 2
        ELSE
          WRITE(LUPRI,'(A,I5)') 'GETLR1 ERROR: INDPRB = ',INDPRB
          CALL QUIT('GETLR1: Inappropriate INDPRB !')
        ENDIF
      ELSE
        WRITE(LUPRI,'(A,I5)') 'GETLR1 ERROR: INDPRB = ',INDPRB
        CALL QUIT('GETLR1: Inappropriate INDPRB !')
      ENDIF
C
C     First get solution vector
C     =========================
C      
      OPEN(LUBVEC,FILE=FILNAM,FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NBDIM,STATUS='UNKNOWN')
      ITIM = IPX
      NXDIM = NBDIM*NFREQ
C
C     If static we have one round with ITIM = IPX, if 
C     non-zero freq. then we have a second round with ITIM = IMX
C
      IROUND = 0
  100 IROUND = IROUND + 1
      CALL DZERO(XPO,NXDIM)
      CALL XRSXV1(ITIM,NBTYP,XPO,EVECR,NFREQ,IBTYP,IBVEC,IVECS,BBUF)
C
C     Print section
C     -------------
C
Cjth   IF(IPRXLR.GE.10) THEN
      IF(IPRXLR.GE.3) THEN
        WRITE(LUPRI,'(/1X,A2,A,I3,A)') BVTYP(ITYP),
     &   '- part of solution vectors for',NFREQ,' frequencies:'
        CALL PRBVEC(LUPRI,XPO,NFREQ,NXPAR)
        CALL FLSHFO(LUPRI)
      ENDIF
C
C     Next form linear response functions
C     ===================================
C     (for now we form only combinations between
C      operators with same point group symmetry)
C
      IF(ALLCMB) THEN
        WRITE(LUPRI,'(A)') 'GETLR1 WARNING: Sorry, not yet written !'
        ALLCMB = .FALSE.
      ENDIF
      IF(.NOT.ALLCMB) THEN
        NOPA   = NLRAP(JSYMOP)
        JOPA   = JLRAP(JSYMOP)
        DO IOPA = 1,NOPA
          INDAP  = LLAPSU(JOPA+IOPA)
          INDPRA = LLRAPU(INDAP)
          KTIMOP = IPRPTIM(INDPRA)
C
C         Get gradient
C         ------------
C
          IF (KTIMOP .EQ. ITIM) THEN
             CALL GPGET(INDPRA,GPO,JXOP,NXPAR,
     &                  CMO,NBTYP,WORK,KFREE,LFREE,IPRXLR)
          END IF
C
C         Form linear reponse function
C         ----------------------------
C
          DO IFREQ = 1,NFREQ
            IF (KTIMOP .EQ. ITIM) ABLRF(IFREQ,INDAP,INDBP,1)  =
     &        D2*DDOT(NBDIM,GPO,1,XPO(1,1,IFREQ),1)
            IF(IPRXLR.GE.2) THEN
              FREQ = EVALR(IFREQ)
              WRITE(LUPRI,1000) 
     &            '<<',INDAP,',',INDBP,'>>:',ABLRF(IFREQ,INDAP,INDBP,1),
     &            ' a.u.,  ',BVTYP(ITYP),
     &            '-contribution at frequency',FREQ,' a.u.'
            ENDIF            
            IF(ORBXLR .AND. KTIMOP.EQ.ITIM) THEN
              DO IZ = 1,NZ
                DO I = 1,NXPAR
                  I1 = JXOP(1,I)
                  IND = IORBCL(2,I1) + 1
                  ABLRF(IFREQ,INDAP,INDBP,IND)=
     &            ABLRF(IFREQ,INDAP,INDBP,IND)+D2*GPO(I,IZ)*
     &              XPO(I,IZ,IFREQ)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDDO        
      ENDIF
C
C
C     Write solution vectors to unformatted file for restart
C
      DO IFREQ = 1,NFREQ
        FREQ = EVALR(IFREQ)
        IF (IROUND.EQ.1 .OR. FREQ.NE.D0) THEN
          CALL WRTRSP(LUXVC,XPO(1,1,IFREQ),
     &         RSPLAB,BVTYP(ITYP),FREQ,JSYMOP,JTIMOP,
     &         BLANK16,BLANK2,D0,0,0,RCNV(IFREQ),INTXRS,NBDIM)
        END IF
      ENDDO
C
C     IF not static, then now calculate the part of the
C     B solution vector which is of -JTIMOP symmetry,
C     calculate contribution to properties and write them to file also.
C
      IF(.NOT.STATIC .AND. IROUND.EQ.1) THEN
        ITIM = IMX
        ITYP = ITYP + 1
        GO TO 100
      END IF
      CLOSE(LUBVEC,STATUS = 'KEEP')
C
      CALL QEXIT('GETLR1')      
      RETURN
 1000 FORMAT(/A,I2,A,I2,A,F16.10,A,3X,A2,A,F16.10,A)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XLROUT */
      SUBROUTINE XLROUT(ABLRF,RCNV,EIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Final output from linear response module
C
C     Written by T.Saue Oct 31 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER(D0 = 0.0D0)
#include <dgroup.h>
#include <pgroup.h>
#include <dcbxlr.h>
#include <dcbxrs.h>
#include <dcbxpr.h>
C
      CHARACTER MXFORM*6,FMT*6,TTYP(-1:1)*2
      DIMENSION ABLRF(NBFREQ,NLRAPT,NLRBPT,*),EIG(*),
     &          RCNV(NBFREQ,NLRBPT),IBUF(2),WORK(*)
      DATA TTYP /'T-','T0','T+'/
C
      IF(ORBXLR) THEN
        NI = 0
C       Number of inactive orbitals
CTROND  Check out active ones...
        DO I = 1,NFSYM
          NI = NI + NCLS(1,1,I) + NCLS(2,1,I)
        ENDDO        
        IOFF = 1
        NVAR = 0
        IF(.NOT.SKIPEE) THEN
          NVAR = NVAR + 1
          IBUF(NVAR)  = IOFF+1
          IOFF = IOFF + NI + 1
        ENDIF
        IF(.NOT.SKIPEP) THEN
          NVAR = NVAR + 1
          IBUF(NVAR)  = IOFF+1
          IOFF = IOFF + NI + 1
        ENDIF
      ENDIF
C
C     Loop over perturbation operators
C     ================================
C      
      CALL HEADER('Linear response; final output ',-1)
      IF(.NOT.ALLCMB) THEN
        DO IOPSY = 1,NFSYM
          NBF = NBSYM/NFSYM
          DO IBF = 1,NBF
            ISYM = JFSYM(IBF,IOPSY)
            NOPA = NLRAP(ISYM)
            JOPA = JLRAP(ISYM)
            NOPB = NLRBP(ISYM)
            JOPB = JLRBP(ISYM)
            DO IOPA = 1,NOPA
              INDAP  = LLAPSU(JOPA+IOPA)
              INDPRA = LLRAPU(INDAP)
              IREPA  = IPRPSYM(INDPRA)-1
              ITIMA  = IPRPTIM(INDPRA)
              DO IOPB = 1,NOPB
                INDBP  = LLBPSU(JOPB+IOPB)
                INDPRB = LLRBPU(INDBP)
                IREPB  = IPRPSYM(INDPRB)-1
                ITIMB  = IPRPTIM(INDPRB)
                ITIMAB = ITIMA*ITIMB
                IF    (ITIMAB.EQ. 1) THEN
                  WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &             '<<A(',INDAP,'),B(',INDBP,
     &             ')>> - linear response function (real):'
                ELSEIF(ITIMAB.EQ.-1) THEN
                  WRITE(LUPRI,'(/A,I2,A,I2,A)')
     &             '<<A(',INDAP,'),B(',INDBP,
     &             ')>> - linear response function (imaginary):'
                ENDIF
                CALL PRSYMB(LUPRI,'-',95,0)
                WRITE(LUPRI,'(3X,A,A16,5X,A3,2X,A2)')
     &             'A - ',PRPNAM(INDPRA),REP(IREPA),TTYP(ITIMA),
     &             'B - ',PRPNAM(INDPRB),REP(IREPB),TTYP(ITIMB)
                CALL PRSYMB(LUPRI,'-',95,0)
                DO I = 1,NBFREQ
                  FMT = MXFORM(ABLRF(I,INDAP,INDBP,1),15)
                  IF(RCNV(I,INDBP).LT.THCLR) THEN
                    WRITE(LUPRI,
     &                '(3X,'//FMT//',A,A,F15.8,A,3X,1P,E9.2,3X,A)')
     &                ABLRF(I,INDAP,INDBP,1),' a.u.',
     &                ' at frequency ',BFREQ(I),' a.u.',
     &                RCNV(I,INDBP),'(converged)'
                  ELSE
                    WRITE(LUPRI,
     &                '(3X,'//FMT//',A,A,F15.8,A,3X,1P,E9.2,3X,A)')
     &                ABLRF(I,INDAP,INDBP,1),' a.u.',
     &                ' at frequency ',BFREQ(I),' a.u.',
     &                RCNV(I,INDBP),'(**NOT** converged!)'
                  ENDIF
                  IF(ORBXLR) THEN
                    WRITE(LUPRI,'(2X,A)') '* Orbital analysis:'
                    IND = 0
                    DO IFRP = 1,NFSYM
C - positronic orbitals
                      DO J = 1,NCLS(1,1,IFRP)
                        IND = IND + 1
                        A = D0
                        DO K = 1,NVAR
                          A = A + ABLRF(I,INDAP,INDBP,IBUF(K)+IND)
                        ENDDO
                        WRITE(LUPRI,'(A3,A1,4(2X,F16.8))')
     &       FREP(IFRP),'p',EIG(IND),A,
     &       (ABLRF(I,INDAP,INDBP,IBUF(K)+IND),K=1,NVAR)
                      ENDDO
C - electronic orbitals
                      DO J = 1,NCLS(2,1,IFRP)
                        IND = IND + 1
                        A = D0
                        DO K = 1,NVAR
                          A = A + ABLRF(I,INDAP,INDBP,IBUF(K)+IND)
                        ENDDO
                        WRITE(LUPRI,'(A3,A1,4(2X,F16.8))')
     &       FREP(IFRP),'e',EIG(IND),A,
     &       (ABLRF(I,INDAP,INDBP,IBUF(K)+IND),K=1,NVAR)
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
                CALL PRSYMB(LUPRI,'-',95,0)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        CALL QUIT('XLROUT: No code for ALLCMB !')
      ENDIF
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck xlrpar */
      SUBROUTINE XLRPAR(WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C
C     Initialize configurational and orbital parameters for
C     RPA calculation
C
C     Written by T.Saue July 18 1997
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbxlr.h>
#include <dcbxrs.h>
#include <dcbbas.h>
#include <dcborb.h>
      DIMENSION WORK(*)
C
      CALL QENTER('XLRPAR')
C
C     Configurational parameters
C     ==========================
C
      KZCONF   = KFREE
      NZCONF   = 0
      NZCONFQ  = NZCONF *NZ
C
C     Orbital parameters
C     ==================
C
      CALL SETXOP(XLR_INDSTR,JOPSY,NZXOPE,KZXOPE,NZXOPP,KZXOPP,
     &            WORK,KFREE,LFREE,SKIPEP,SKIPEE,QEDRUN,IPRINT)      
      NZXOPEQ  = NZXOPE *NZ
      NZXOPPQ  = NZXOPP *NZ
      NZXOPT   = NZXOPE + NZXOPP
      NZXOPTQ  = NZXOPT *NZ
C
C     Total number of parameters
C     ==========================
C
      NZVAR    = NZCONF + NZXOPT
      NZVARQ   = NZVAR*NZ
C
      CALL QEXIT('XLRPAR')      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Anaxlr */
      SUBROUTINE ANAXLR(WORK,KFREE,LFREE)
C***********************************************************************
C
C     Analyze solution vectors on file PAMXVC
C
C     Written by Trond Saue July 4 2001
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbxrs.h>
#include <dcbxlr.h>
      DIMENSION WORK(*)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2,LBUF*16,PBUF*1
      LOGICAL TOBE   
C
      KFRSAV = KFREE   
      CALL HEADER('Analysis of LR solution vectors',-1)
C
C     Check if file of solution vectors is present
C
      INQUIRE(FILE='PAMXVC',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      LUINP = 1
      CALL OPNFIL(LUINP,'PAMXVC','OLD','XLRANA')
      IREC = 0
      IERR = 0
      IVEC = 0
      LBUF = '                '
      TOBE = .FALSE.
      JFRP = -1
      CALL OSTRING(WORK,KFREE,LFREE)
C
C     WHILE DO - loop
C
    1 CONTINUE
        READ (LUINP,END=3,ERR=2) LAB1
        IF(LAB1.EQ.'END_OF_THIS_FILE') THEN
          IF(TOBE) THEN
            NVEC = IVEC-IPAR
            NBDIM = NDIM/NZ
            IF(PBUF.EQ.'E') THEN
              IF(SKIPEE) GOTO 20
              CALL ANAXL1(NBDIM,NVEC,WORK(KZXOPE),LUINP,LAB1,LAB2,
     &                    TYP1,TYP2,WORK,KFREE,LFREE,IPRXLR)
            ELSEIF(PBUF.EQ.'P') THEN
              IF(SKIPEP) GOTO 20
              CALL ANAXL1(NBDIM,NVEC,WORK(KZXOPP),LUINP,LAB1,LAB2,
     &                    TYP1,TYP2,WORK,KFREE,LFREE,IPRXLR)
            ENDIF
            GOTO 20
          ENDIF
        ENDIF
        BACKSPACE LUINP
        READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,
     &            FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &            JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &            ERGRSP,NBAS,NORB
        IVEC = IVEC + 1
        IREC = IREC + 1
C       New parameter class
        IF(TYP1(1:1).NE.PBUF.OR.LAB1.NE.LBUF) THEN
          IF(TOBE) THEN
C           Analyze solution vectors
            NVEC = IVEC-IPAR
            NBDIM = NDIM/NZ
            IF(PBUF.EQ.'E') THEN
              IF(SKIPEE) GOTO 30
              CALL ANAXL1(NBDIM,NVEC,WORK(KZXOPE),LUINP,LAB1,LAB2,
     &                    TYP1,TYP2,WORK,KFREE,LFREE,IPRXLR)
            ELSEIF(PBUF.EQ.'P') THEN
              IF(SKIPEP) GOTO 30
              CALL ANAXL1(NBDIM,NVEC,WORK(KZXOPP),LUINP,LAB1,LAB2,
     &                    TYP1,TYP2,WORK,KFREE,LFREE,IPRXLR)
            ENDIF
            READ (LUINP,END=3,ERR=2) LAB1,LAB2,TYP1,TYP2,
     &              FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &              JTIMOP1,JTIMOP2,RNORM,LEN,INTFLG,
     &              ERGRSP,NBAS,NORB
          ENDIF
 30       CONTINUE
          PBUF = TYP1(1:1)
          IPAR = IVEC
          NDIM = LEN
          TOBE = .TRUE.
        ENDIF
C       New property
        IF(LAB1.NE.LBUF) THEN
          IREP = JSYMOP1-1
          IFRP = JBTOF(IREP,1)
C         Generate orbital rotation operator
          IF(IFRP.NE.JFRP) THEN
            IF(TOBE) THEN
              CALL MEMREL('ANAXLR',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
            ENDIF
            CALL SETXOP(XLR_INDSTR,IFRP,NZXOPE,KZXOPE,NZXOPP,KZXOPP,
     &                  WORK,KFREE,LFREE,SKIPEP,SKIPEE,QEDRUN,IPRXLR)
            JFRP = IFRP
          ENDIF
C
          LBUF = LAB1
        ENDIF
        READ(LUINP,END=3,ERR=2)
        IREC = IREC + 1
      GO TO 1
C
    2 CONTINUE
      IREC = IREC + 1
      IERR = IERR + 1
      WRITE (LUPRI, '(/A,I5/)') ' ERROR READING RECORD NO.',IREC
      GOTO 20
C
    3 CONTINUE
      WRITE (LUPRI,'(/I10,A)') IREC,
     *   ' records read before EOF on file.'
      GOTO 20
 20   CONTINUE
      CLOSE(LUINP)
      CALL MEMREL('ANAXLR',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'ANAXLR: No solution vectors found !'
      CALL QUIT('ANAXLR: Solution vectors not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck anaxl1 */
      SUBROUTINE ANAXL1(NBDIM,NVEC,JXOP,LUINP,LAB1,LAB2,TYP1,TYP2,
     &                  WORK,KFREE,LFREE,IPRINT)
C***********************************************************************
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
      DIMENSION JXOP(*),WORK(*)
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KVEC ,NBDIM*NZ*NVEC,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFREQ,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KNORM,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF1,NBDIM*NZ     ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF2,NBDIM        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KBTYP,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFLAG,NVEC         ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDX,NBDIM        ,WORK,KFREE,LFREE)
      CALL ANAXL2(WORK(KVEC),WORK(KFREQ),WORK(KNORM),
     &            WORK(KBTYP),WORK(KFLAG),WORK(KBUF1),
     &            WORK(KBUF2),WORK(KINDX),JXOP,
     &            NBDIM,NVEC,LUINP,LAB1,LAB2,TYP1,TYP2,IPRINT)
      CALL MEMREL('ANAXL1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Anaxl2 */
      SUBROUTINE ANAXL2(BVEC,FREQ,RNORM,IBTYP,LFLAG,BUF1,BUF2,INDX,
     &          JXOP,NBDIM,NVEC,LUINP,LAB1,LAB2,TYP1,TYP2,IPRINT)
C***********************************************************************
C
C     Written by T.Saue July 6 2001
C
C***********************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DTOL=1.0D-12,D0=0.0D0,D1=1.0D0)
C
#include <dgroup.h>
#include <pgroup.h>
      CHARACTER LAB1*16,LAB2*16,TYP1*2,TYP2*2
      LOGICAL LFLAG(NVEC)
      DIMENSION BVEC(NBDIM,NZ,NVEC),FREQ(NVEC),RNORM(NVEC),IBTYP(NVEC),
     &          BUF1(NBDIM,NZ),BUF2(NBDIM),INDX(NBDIM),JXOP(2,NBDIM)
C
      BACKSPACE LUINP
      DO I = 1,NVEC
        BACKSPACE LUINP
        BACKSPACE LUINP
      ENDDO
C     Read vectors
      READ (LUINP) LAB1,LAB2,TYP1,TYP2,
     &        FREQ1,FREQ2,JSYMOP1,JSYMOP2,
     &        JTIMOP1,JTIMOP2,DUM,LEN,INTFLG,
     &        ERGRSP,NBAS,NORB
      IREP = JSYMOP1 - 1
      WRITE(LUPRI,'(A,A1,A,A16,A,A3,A,I3,A,I10)') 
     &     '** ',TYP1(1:1),'-solution vectors : ',LAB1,
     &     ' Irrep: ',REP(IREP),' Trev: ',JTIMOP1,
     &     ' Length: ',NBDIM
      BACKSPACE LUINP
      DO I = 1,NVEC
        READ (LUINP) LAB1,LAB2,TYP1,TYP2,
     &        FREQ(I),FREQ2,JSYMOP1,JSYMOP2,
     &        JTIMOP1,JTIMOP2,RNORM(I),LEN,INTFLG,
     &        ERGRSP,NBAS,NORB
        IF(TYP1(2:2).EQ.'+') THEN
          IBTYP(I) = 1
        ELSEIF(TYP1(2:2).EQ.'-') THEN
          IBTYP(I) = -1
        ELSE
          CALL QUIT('XLRANA: Unknown IBTYP !')
        ENDIF
        CALL READT(LUINP,LEN,BVEC(1,1,I))
        IF(IPRINT.GE.4) THEN
          WRITE(LUPRI,'(3X,A2,A,F10.6)') 
     &      TYP1,' vector at frequency: ', FREQ(I)
          CALL PRQMAT(BVEC(1,1,I),NBDIM,1,NBDIM,1,
     &                NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDDO
C     Analyze major excitations
      CALL LSET(NVEC,.TRUE.,LFLAG)      
      DO I = 1,NVEC
      IF(LFLAG(I)) THEN
        LFLAG(I) = .FALSE.
        DO IZ = 1,NZ
          CALL DCOPY(NBDIM,BVEC(1,IZ,I),1,BUF1(1,IZ),1)
        ENDDO
        IF(DABS(FREQ(I)-D0).GT.DTOL) THEN
C       Frequency dependent case; search for partner
          DO J = 1,NVEC
          IF(LFLAG(J)) THEN
            IF(DABS(FREQ(I)-FREQ(J)).LT.DTOL) THEN
              LFLAG(J)=.FALSE.
              DO IZ = 1,NZ
                CALL DAXPY(NBDIM,D1,BVEC(1,IZ,J),1,BUF1(1,IZ),1)
              ENDDO
              GOTO 10
            ENDIF
          ENDIF
          ENDDO
          WRITE(LUPRI,'(A,I5)') 
     &      ' WARNING ! No partner found for vector. ',I
        ENDIF
 10     CONTINUE
C       Calculate norm
        DNTOT = D0
        DO J = 1,NBDIM
          BUF2(J) = BUF1(J,1)*BUF1(J,1)
        ENDDO
        DO IZ = 2,NZ
          DO J = 1,NBDIM
            BUF2(J) = BUF2(J) + BUF1(J,IZ)*BUF1(J,IZ)
          ENDDO
        ENDDO
        DO J = 1,NBDIM
          DNTOT = DNTOT + BUF2(J)
          BUF2(J) = SQRT(BUF2(J))
        ENDDO
        THR   = SQRT(DNTOT/NBDIM)
        DNTOT = SQRT(DNTOT)
        WRITE(LUPRI,'(A8,F10.6,A,E16.8,A,E10.2)')
     &        ' Freq.: ',FREQ(I),' Norm: ',DNTOT,
     &        ' Residual norm: ',RNORM(I)
C       Establish order
        CALL INDEXX (NBDIM,BUF2,INDX)
        WRITE(LUPRI,'(3X,A,E10.2/)') 
     &   'Excitation amplitudes larger than threshold ',THR
        DO J = NBDIM,1,-1
          JJ = INDX(J)
          IF(BUF2(JJ).GE.THR) THEN
            WRITE(LUPRI,'(3X,I6,A,I6,E16.8)') 
     &       JXOP(1,JJ),' ---> ',JXOP(2,JJ),BUF2(JJ)
          ENDIF
        ENDDO
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(3X,A/)') 'All excitation amplitudes :'
          DO J = 1,NBDIM
            WRITE(LUPRI,'(3X,I6,E16.8)') J,BUF2(J)
          ENDDO
        ENDIF
      ENDIF
      ENDDO
      RETURN
      END


