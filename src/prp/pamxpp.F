!dirac_copyright_start
!      Copyright (c) 2012 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org
!dirac_copyright_end

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRPXPP */
      SUBROUTINE PRPXPP(WORK,KPPF,KOMEGA,KCNV,KPPO,KSSYM,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate "polarization propagator" (PP) excitation energies
C     and transition moments defined in /XCBXPP/
C     (= poles and residues of principal linear response function)
C
C     Linear response functions are defined as:
C
C        <<A;B>>_w = -(E[1]_A)^{+} (E[2] - wS[2])^{-1} (E[1]_B)
C
C        E[1]_A   is gradient of perturbation operator A  (h|A)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C
C     where h is the operator manifold. In RPA h is limited to single
C     excitations.
C
C     The pole no. j of the principal linear response functions is found
C     by iterative techniques.  We write
C
C        (E[2] - w_j S[2]) X_j = 0
C
C     which can be rearranged to the generalized eigenvalue equation
C
C        E[2] X_j =  w_j S[2] X_j
C
C     Here X_j is the j'th excitation operator
C     and  w_j is the j'th excitation energy.
C
C     The eigenvalue equation is solved by projecting the equation onto
C     a reduced space of trial vectors:
C
C        B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the reduced equation
C
C        RE[2] a_j = Rw_j RS[2] a_j
C
C     in which
C
C       RE[2] = B^{+}E[2]B
C       RS[2] = B^{+}S[2]B
C       Rw_j  = approximate w_j
C
C
C     The evaluation of excitations thus proceeds in three steps:
C
C       1. Generate initial guess        : X_j^(0) = b_j
C       2. Solve reduced eigenvalue eq.  : RE[2] a_j = Rw_j RS[2] a_j
C       3. Generate full solution vector : X_j = B a_j
C
C
C     Each excitation operator X_j is assumed to span a given boson irrep
C     and is both symmetric(+) or antisymmetric(-) with respect
C     to time reversal.
C
C     Written by Hans Joergen Aa. Jensen, 2000
C     Based on PRPXLR
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dgroup.h"
      DIMENSION WORK(*)
C
      CALL QENTER('PRPXPP')
C
C     ****************************
C     *** I N I T I A L I Z E  ***
C     ****************************
C
      CALL SETXPP(WORK,KFREE,LFREE)
C
      IF (MAXEXC.EQ.0) THEN
        IF(XPPANA) GOTO 10
        CALL QUIT('PRPXPP called for zero excitations!')
      ENDIF
      NEXCF  = MAXEXC*NBSYM
      IF(ORBXPP) THEN
        CALL MEMGET('INTE',KPPO,4*NORBT,WORK,KFREE,LFREE)
        CALL DEFORB(WORK(KPPO),NCLS_XPP,3,XPP_INDSTR,
     &              WORK,KFREE,LFREE,IPRXPP)
        NI = 0
C       Number of inactive orbitals
CTROND  Check out active ones...
        DO I = 1,NFSYM
          NI = NI + NCLS_XPP(1,1,I) + NCLS_XPP(2,1,I)
        ENDDO
      ELSE
        NI = 0
        CALL MEMGET('INTE',KPPO,0,WORK,KFREE,LFREE)
      ENDIF
      IF(XPP_SKIPEE) THEN
        NEE = 0
      ELSE
        NEE = NI + 1
      ENDIF
      IF(XPP_SKIPEP) THEN
        NEP = 0
      ELSE
        NEP = NI + 1
      ENDIF
      NCI = 0
      N2PPF  = NPPAPT*NEXCF*(1+NEE+NEP+NCI)
      CALL MEMGET('REAL',KPPF,N2PPF,WORK,KFREE,LFREE)
C     ... ATMPPF: "A" transition moment for each excitation
      CALL MEMGET('REAL',KCNV,NEXCF,WORK,KFREE,LFREE)
C     ... RCNV: residual error for each excitation
      CALL MEMGET('REAL',KOMEGA,NEXCF,WORK,KFREE,LFREE)
C     ... OMEGA: excitation energies
      CALL MEMGET('INTE',KSSYM,NEXCF,WORK,KFREE,LFREE)
C     ... SYM: supersymmetry (if known)
C
C     Initialize:
      CALL DZERO(WORK(KPPF),N2PPF)
      CALL DZERO(WORK(KOMEGA),NEXCF)
      CALL PRPXP1(WORK(KPPF),WORK(KOMEGA),WORK(KCNV),WORK(KPPO),
     &            WORK(KSSYM),WORK,KFREE,LFREE)
C
 10   CONTINUE
      IF (XPPANA) CALL ANAXRS(XPP_INDSTR,XPP_SKIPEE,XPP_SKIPEP,
     &                        WORK,KFREE,LFREE,IPRXPP)
C
      CALL MEMCHK('PRPXPP memory check',WORK,1)
C     ... note: NOT MEMREL here as info is passed to calling
C         routine in WORK(KPPO), WORK(KPPF), WORK(KCNV), WORK(KOMEGA)
      CALL QEXIT('PRPXPP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRPXP1 */
      SUBROUTINE PRPXP1(ATMPPF,OMEGA,RCNV,IORBCL,ISSYM,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Calculate excitation energies and transition moments.
C     Called by PRPXPP (see PRPXPP for general description of task).
C
C     ATMPPF : "A" transition moment for each excitation
C     OMEGA  : calculated excitation energies
C     RCNV   : residual error for each excitation
C     IORBCL : for orbital analysis
C     ISSYM  : supersymmetry
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
#include "dcbxrs.h"
C
      DIMENSION ATMPPF(MAXEXC,NPPAPT,NBSYM,*),OMEGA(MAXEXC,NBSYM),
     &          RCNV(MAXEXC,NBSYM), IORBCL(4,*),
     &          ISSYM(MAXEXC,NBSYM),WORK(*)
      CHARACTER TTYP(-1:1)*2
      SAVE TTYP
      DATA TTYP /'T-','T0','T+'/
C
      CALL QENTER('PRPXP1')
      KFRSAV = KFREE
      CALL IZERO(ISSYM,MAXEXC*NBSYM)
C
C     If orbital analysis; extract relevant orbital energies
C     ======================================================
C
      IF(ORBXPP) THEN
        NI = 0
C       Number of inactive orbitals
CTROND  Check out active ones...
        DO I = 1,NFSYM
          NI = NI + NCLS_XPP(1,1,I) + NCLS_XPP(2,1,I)
        ENDDO
        
        CALL MEMGET('REAL',KEIG,NI,WORK,KFREE,LFREE)
        CALL EIGEXT(WORK(KEIG),IORBCL,1,WORK,KFREE,LFREE)
      ELSE
        NI = 0
        CALL MEMGET('REAL',KEIG,NI,WORK,KFREE,LFREE)
      ENDIF
C
C     ************************************************
C     ***                                          ***
C     *** E X C I T A T I O N   F U N C T I O N S  ***
C     ***                                          ***
C     ************************************************
C
C     Memory allocation
C     =================
C
C     Dimensions of reduced system
C
      KXPP   = KFREE
      NEVECR = NREDM*MAXEXC
      CALL MEMGET('INTE',KIBTYP,2*NREDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBCVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBPVC,NREDM  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVECR,NEVECR ,WORK,KFREE,LFREE)
C
C     ... no B property gradient in XPP
C
      CALL MEMGET('REAL',KGPCI,0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGPOE,0,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGPOP,0,WORK,KFREE,LFREE)
C
C     Loop over excitation symmetries
C     ===============================
C
      JTIMOP = 1
Chj: must be set to define IPX, IMX in pamrsp.F /Nov.2000
      DO 10 JOPSY = 1,NFSYM
        NFC = NBSYM/NFSYM
        DO 20 IS = 1,NFC
          JSYMOP = JFSYM(IS,JOPSY)
          NEXCNV = KEXCNV(JSYMOP)
          IF (NEXCNV .EQ. 0) GO TO 20
          WRITE (LUPRI,'(///A,A3/A,I3,A)')
     &      ' >>>  SOLVING PP EXCITATION ENERGIES - Symmetry: ',
     &        REP(JSYMOP-1),
     &      '      for the',NEXCNV,' lowest excitation energies.'
C
C         Initialize configurational and orbital parameters
C         =================================================
C
          CALL XRSPAR(XPP_INDSTR,XPP_SKIPEE,XPP_SKIPEP,
     &                WORK,KFREE,LFREE,IPRXPP)
C
C         Transfer number of excitations for this symmetry to dcbxrs.h
C         (NFREQ is used for allocations and loops)
C
          IF (NZCONF.EQ.0) THEN
             NEXMAX = NZXOPEQ
          ELSE
             NEXMAX = NZCONF
          END IF
C         ... NEXMAX is in agreement with the assignment
C             of start vectors in XRSST
          IF (NEXCNV .GT. NEXMAX) THEN
             WRITE (LUPRI,'(/A/A,I5)')
     &       ' WARNING: no. of excitation energies exceeds no. of var.',
     &       ' The number is decreased to',NEXMAX
             NEXCNV = NEXMAX
          END IF
          NEXSIM = MIN(KEXSIM(JSYMOP), NEXMAX)
          NFREQ  = NEXSIM
          NEXSTV = MIN(KEXSTV(JSYMOP), NEXMAX)
C
C         *******************************************************
C         *** S O L V E   E I G E N V A L U E   E Q U A T . S ***
C         *******************************************************
C
          NCRED  = 0
          NERED  = 0
          NPRED  = 0
          NZRED  = NCRED + NERED + NPRED
          RSPLAB = 'PP EXCITATION'//REP(JSYMOP-1)
          CALL XRSCTL(WORK(KGPCI) ,WORK(KGPOE),WORK(KGPOP),
     &                WORK(KIBTYP),WORK(KIBCVC),WORK(KIBEVC),
     &                WORK(KIBPVC),RCNV(1,JSYMOP),
     &                OMEGA(1,JSYMOP),WORK(KEVECR),
     &                WORK,KFREE,LFREE)
C
C         ************************************************************
C         *** E V A L U A T E   R E S P O N S E   F U N C T I O N  ***
C         ************************************************************
C
          CALL GETPPF(ATMPPF,WORK(KIBTYP),WORK(KIBCVC),
     &                WORK(KIBEVC),WORK(KIBPVC),RCNV(1,JSYMOP),
     &                OMEGA(1,JSYMOP),WORK(KEVECR),ISSYM(1,JSYMOP),
     &                IORBCL,WORK(KEIG),WORK,KFREE,LFREE)
    
          CALL XPPOUT_OLD(JSYMOP,ATMPPF,OMEGA,RCNV,WORK(KEIG),
     &                ISSYM,WORK,KFREE,LFREE)
 20     CONTINUE
 10   CONTINUE
      CALL MEMREL('PRPXP1',WORK,1,KXPP,KFREE,LFREE)
C
C     Final print section
C
      CALL XPPOUT(0,ATMPPF,OMEGA,RCNV,WORK(KEIG),ISSYM,WORK,KFREE,LFREE)
      CLOSE(LUXVC,STATUS = 'KEEP')
      CALL MEMREL('PRPXP1',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('PRPXP1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xppinp */
      SUBROUTINE XPPINP(WORD,RESET,WORK,LWORK)
C*****************************************************************************
C
C     Input section for PP excitation energies
C     (single residue of linear response function)
C
C     Written by H.J.Aa.Jensen 2000
C     Based on XLRINP
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
#include "dummy.h"
#include "maxorb.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 32,DRESFAC = 1.0D2)
C
#include "dcbgen.h"
#include "dcbprp.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
C
      LOGICAL SET, NEWDEF, RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
      DIMENSION WORK(LWORK)
C
      SAVE SET
      DATA SET/.FALSE./
      DATA TABLE /'.PRINT ','.XXXXXX','.ANALYZ','.OPERAT',
     &            '.EXCITA','.TRIPLE','.THRESH','.MAXITR',
     &            '.SKIPEP','.MAXRED','.E2CHEK','.NONORM',
     &            '.EPOLE ','.MPOLE ','.INTFLG','-------',
     &            '.OCCUP ','.VIRTUA','.ACTIVE','-------',
     &            '.ITRINT','.CNVINT','.REAXVC','-------',
     &            '.RESFAC','.ORBANA','-------','.SKIPEE',
     &            '.MXLOAD','.ONLYSF','.ONLYSG','-------'/
C
#include "ibtfun.h"
C
      NEWDEF = (WORD .EQ. '*EXCITA')
      IF (SET) THEN
         IF (NEWDEF)
     &      CALL QUIT('Only one "*EXCITA" input section allowed')
C        hjaaj: repeated input sections give infinite loop ...
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
chjMAERKE: ought to be dependent on wave function INTFLG ??
C
      ILLINT = IBTAND(INTGEN,1)
      ISLINT = IBTAND(INTGEN/2,1)
      ISSINT = IBTAND(INTGEN/4,1)
C
C
C     Initialize DCBXPP
C     =================
C
      IPRXPP     = 0
      THCPP      = 1.0D-5
      RESXPP     = DRESFAC
      CNVXPP(1)  = DUMMY
      CNVXPP(2)  = DUMMY
      ITRXPP     = 50
      INTXPP     = ILLINT + 2*ISLINT + 4*ISSINT
      IF (LEVYLE) INTXPP = ILLINT
      ITRIPP(1)  = 1
      ITRIPP(2)  = 1
      NPPAPT     = 0
      CALL IZERO(NPPAP,8)
      MAXEXC     = 0
      CALL IZERO(KEXCNV,8)
      CALL IZERO(KEXSIM,8)
      CALL IZERO(KEXSTV,8)
      XPPANA     = .FALSE.
      XPPNRM     = .TRUE.
      XPPFIL     = 'NORXVC'
      ORBXPP     = .FALSE.

Chj: initialize some parts of dcbxlr; (nov. 2000)
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
      MAXRM      = -1
      MXLOAD     = -1
      XPP_SKIPEP = .FALSE.
      XPP_SKIPEE = .FALSE.
      XPP_LSFG(1)= .TRUE.
      XPP_LSFG(2)= .TRUE.
      XPP_TRIPLET = .FALSE.
      XPP_E2CHEK = .FALSE.
      DO I = 1,NFSYM
         XPP_INDSTR(2,I) = ' '
         XPP_INDSTR(1,I) = ' '
         XPP_INDSTR(3,I) = ' '
      ENDDO
C
C     Initialize user-defined transition properties
C     =============================================
C
C     Electric dipole oscillator strengths
C
chjMAERKE: implement oscillator strengths here ? TODO
C
C     Process input
C     =============
C
      ICHANG = 0
      INPERR = 0
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
     *            '" not recognized in XPPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in XPPINP.')
    1          CONTINUE
C&&&& PRINT:  Print level
                  READ(LUCMD,*) IPRXPP
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
C&&&& ANALYZe: Analyze individual contributions to repsonse functions
                  XPPANA = .TRUE.
               GO TO 100
    4          CONTINUE
C&&&& OPERATOR: operator of transition moment
                  CALL XPRINP(LUCMD,WORD,INPERR,INDXPR,ISYXPR,ITRXPR,
     &                        WORK,LWORK,IPRXPP)
                  CALL OP1IND('NPPAPT',IND1OP,LPPAPU,NPPAPT,INDXPR,
     &                        MAXLPP)
               GO TO 100
    5          CONTINUE
C&&&& EXCITAtions - how many roots in symmetry I
                  READ(LUCMD,*) I, KEXCNV(I)
CHJMAERKE 2000: make input for KEXSIM and KEXSTV ?
                  KEXSIM(I) = KEXCNV(I)
                  KEXSTV(I) = KEXCNV(I)
               GO TO 100
    6          CONTINUE
C&&&& TRIPLET - allow triplet excitations in spinfree calculations
               XPP_TRIPLET = .TRUE.
               GO TO 100
    7          CONTINUE
C&&&& THRESH - Threshold for convergence:
                  READ (LUCMD,*) THCPP
               GO TO 100
    8          CONTINUE
C&&&& MAXITR - Maximum number of iterations
                  READ (LUCMD, '(I5)') ITRXPP
               GO TO 100
    9          CONTINUE
C&&&& SKIPEP -exclude all electron-positron rotations
                  XPP_SKIPEP = .TRUE.
               GOTO 100
 10            CONTINUE
C&&&& MAXRED: max dimension of matrix in reduced system
                  READ(LUCMD,*) MAXRM
               GOTO 100
 11            CONTINUE
C&&& E2CHEK: explicit contruction of the electronic Hessian E2
               XPP_E2CHEK = .TRUE.
               GOTO 100
 12            CONTINUE
C&&&& NONORM: normalize trial vectors
                 XPPNRM = .FALSE.
               GOTO 100
 13            CONTINUE
C&&&& EPOLE : Electric multipoles of order L
               READ(LUCMD,*) IORDER
               NDEG=(IORDER+2)*(IORDER+1)/2
               IF((NPPAPT+NDEG).GT.MAXLPP) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'XPPINP: Pointer array LPP',
     &           'APU out of bounds for electric multipole of order ',
     &           IORDER
                 WRITE(LUPRI,'(A,I5)') 
     &           'Increase MAXEXP to ',(NPPAPT+NDEG)
               ENDIF
               CALL DEF_EPOLE(IORDER,LPPAPU(NPPAPT+1),IPRXPP)
               NPPAPT = NPPAPT + NDEG
               GOTO 100
 14            CONTINUE
C&&&& MPOLE : Magnetic multipoles of order L 
               READ(LUCMD,*) IORDER
               NDEG=3*(IORDER+1)*IORDER/2
               IF((NPPAPT+NDEG).GT.MAXLPP) THEN
                 WRITE(LUPRI,'(A,A,I3)') 'XPPINP: Pointer array LPP',
     &           'APU out of bounds for magnetic multipole of order ',
     &           IORDER
                 WRITE(LUPRI,'(A,I5)') 
     &           'Increase MAXEXP to ',(NPPAPT+NDEG)
               ENDIF
               CALL DEF_MPOLE(IORDER,LPPAPU(NPPAPT+1),IPRXPP)
               NPPAPT = NPPAPT + NDEG
               GOTO 100
C&&&& INTFLG: specify what two-electron integrals to include
 15            CONTINUE
                  READ(LUCMD,*) ILLINT,ISLINT,ISSINT
                  INTXPP = ILLINT + 2*ISLINT + 4*ISSINT
               GOTO 100
 16            CONTINUE
               GOTO 100
 17            CONTINUE
C&&&& OCCUP: String of inactive orbitals in XPP module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XPP_INDSTR(1,I)
                  ENDDO
               GO TO 100
 18            CONTINUE
C&&&& VIRTUA: String of secondary orbitals in XPP module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XPP_INDSTR(3,I)
                  ENDDO
               GO TO 100
 19            CONTINUE
C&&&& ACTIVE: String of active orbitals in XPP module
                  DO I=1,NFSYM
                     READ(LUCMD,'(A)') XPP_INDSTR(2,I)
                  ENDDO
               GO TO 100
 20            CONTINUE
               GO TO 100
 21            CONTINUE
C&&& ITRINT: Number of iterations before adding SL- and SS-integrals
                  READ(LUCMD,*) ITRIPP(1),ITRIPP(2)
               GO TO 100
 22            CONTINUE
C&&& CNVINT: Convergence thresholds for adding SL- and SS-integrals
                  READ(LUCMD,*) CNVXPP(1),CNVXPP(2)
               GO TO 100
 23            CONTINUE
C&&&& REAXVC: Read solution vectors from file XPPFIL
                  READ(LUCMD,'(A6)') XPPFIL
               GO TO 100
 24            CONTINUE
               GO TO 100
 25            CONTINUE
C&&&& RESFAC:  New trial vectors only for parameters whose residual
C              has a norm that is a fraction RESFAC of the max norm
                  READ(LUCMD,*) RESXPP
               GO TO 100
 26            CONTINUE
C&&&& ORBANA: Analyze individual contributions to response functions
                  ORBXPP = .TRUE.
               GO TO 100
 27            CONTINUE
               GO TO 100
 28            CONTINUE
C&&&& SKIPEE: Skip all e-e rotations
                  XPP_SKIPEE = .TRUE.
               GO TO 100
 29            CONTINUE
C&&&& MXLOAD: Maximum number of trial vectors to be loaded
C             in each microiteration
                  READ(LUCMD,*) MXLOAD
               GO TO 100
 30            CONTINUE
C&&&& ONLYSF: Only generate one-index transformed Fock matrix
                  XPP_LSFG(1) = .FALSE.
               GO TO 100
 31            CONTINUE
C&&&& ONLYSG: Only generate G matrix (Fock-matrix of modified density)
                  XPP_LSFG(2) = .FALSE.
               GO TO 100
 32            CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in XPPINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in XPPINP.')
            END IF
      END IF
  300 CONTINUE
C
C     Process section
C     ===============
C
      IF(CNVXPP(1).LT.DUMMY) ITRIPP(1) = 1
      IF(CNVXPP(2).LT.DUMMY) ITRIPP(2) = 1
C
      MAXEXC = MAX( KEXCNV(1),KEXCNV(2),KEXCNV(3),KEXCNV(4),
     &              KEXCNV(5),KEXCNV(6),KEXCNV(7),KEXCNV(8) )
      IF (MAXRM.LT.0) MAXRM = MAX(200,100+MAXEXC*ITRXPP)
C        calculate default MAXRM
      MXREDZ = MAXRM/2
      MAXRM  = 2*MXREDZ
      IF(NPSHT .EQ. 0) XPP_SKIPEP = .TRUE.
C     ... no positronic orbitals
      IF(LEVYLE.OR.BSS.OR.X2C.OR.FREEPJ.OR.VEXTPJ) XPP_SKIPEP = .TRUE.
C     ... no positronic orbitals for Levy-Leblond or no-pair
      IF (XPP_SKIPEE .AND. XPP_SKIPEP) THEN
         WRITE(LUPRI,'(/A)') 'INPUT ERROR: both SKIPEE and SKIPEP'
         INPERR = INPERR + 1
      END IF
      IF(QED .AND. XPP_SKIPEP) THEN
         WRITE(LUPRI,'(/A)')
     &  ' INPUT ERROR: QED invalid for Levy-Leblond/no-pair/.SKIPEP'
         INPERR = INPERR + 1
      END IF
C
C     Print section
C     =============
C
      DOXPP = DOXPP.OR.MAXEXC.GT.0
      IF(.NOT.(DOXPP.OR.XPPANA)) GOTO 999
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(A)') ' XPPINP: Electronic PP excitation energies'//
     &    ' and transition moments'
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(1X,A,I5)')
     &   '* Print level                 :',IPRXPP
      IF(XPPANA) WRITE(LUPRI,'(A)')
     &   ' * Analysis of solution vectors.'
      IF(XPP_E2CHEK) WRITE(LUPRI,'(A)')
     &   ' * Explicit contruction of the electronic Hessian'
      IF(ORBXPP) WRITE(LUPRI,'(A)')
     &   ' * Analysis of individual orbital contributions.'
      IF(XPPFIL.NE.'NORXVC') WRITE(LUPRI,'(A,A,A6)')
     &    ' * The program will attemp restart on solution vectors ',
     &    'on file ',XPPFIL
      WRITE(LUPRI,'(A,I5)')
     &   ' * Maximum number of iterations:',ITRXPP
      WRITE(LUPRI,'(A,1P,D12.5)')
     &   ' * Threshold for convergence   :',THCPP
      WRITE(LUPRI,'(A,I5)')
     &   ' * Max. size of reduced matrix :',MAXRM
      IF(XPP_SKIPEE) WRITE(LUPRI,'(A)') ' * NB!!! No e-e rotations'
      IF(XPP_SKIPEP .AND. .NOT.(LEVYLE.OR.FREEPJ.OR.VEXTPJ))
     &   WRITE(LUPRI,'(A)')
     &   ' * NB!!! No rotations between electrons and positrons.'
      IF(QED) WRITE(LUPRI,'(A)')
     &  ' * QED: Rotations between positrons and all electron orbitals.'
      IF(.NOT.DOXPP) THEN
        DOXPP = DOXPP.OR.XPPANA
        GOTO 999
      ENDIF
      IF(SPINFR.AND.XPP_TRIPLET)WRITE(LUPRI,'(A)')
     &  ' * Triplet excitations allowed.'
C
      WRITE(LUPRI,'(A)') ' * Include contributions from '//
     +    'the following two-electron integrals:'
      IF(ILLINT.EQ.1) WRITE(LUPRI,'(3X,A)') '- LL-integrals'
      IF(ISLINT.EQ.1) THEN
        IF(CNVXPP(1).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,D8.1)')
     +    '- SL-integrals below convergence ',CNVXPP(1)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    '- SL-integrals from iteration ',ITRIPP(1)
        ENDIF
      ENDIF
      IF(ISSINT.EQ.1) THEN
        IF(CNVXPP(2).LT.DUMMY) THEN
          WRITE(LUPRI,'(3X,A,D8.1)')
     +    '- SS-integrals below convergence ',CNVXPP(2)
        ELSE
          WRITE(LUPRI,'(3X,A,I4)')
     +    '- SS-integrals from iteration ',ITRIPP(2)
        ENDIF
      ENDIF
      IF(XPPNRM) WRITE(LUPRI,'(A)')
     &     ' * Trial vectors will be normalized.'
      IF(MXLOAD.GT.0) WRITE(LUPRI,'(A,I5)')
     &  ' * Maximum number of trial vectors to be loaded'//
     &  ' in each microiteration:', MXLOAD
      IF(.NOT.XPP_LSFG(1)) WRITE(LUPRI,'(A)')
     &     ' * Debug: Only calling FMOLI in sigma-vector routine...'
      IF(.NOT.XPP_LSFG(2)) WRITE(LUPRI,'(A)')
     &     ' * Debug: Only calling GMOLI in sigma-vector routine...'
C
C     Number of excitations
C     =====================
C
      CALL PRSYMB(LUPRI,'-',75,0)
      WRITE(LUPRI,'(A,8(/10X,A,I2,3A,I5))')
     &   ' * No. of excitations evaluated in each boson symmetry:',
     &   (' boson symmetry',J,' (',REP(J-1),') :',KEXCNV(J),J=1,NBSYM)
C
C     A operators
C     ===========
C
      CALL PRSYMB(LUPRI,'-',75,0)
      IF (NPPAPT .EQ. 0) THEN
        WRITE(LUPRI,'(A)') ' ** No transition moment operators **'
      ELSE
        WRITE(LUPRI,'(A)') ' ** Transition moment operators **'
        CALL PRSYMB(LUPRI,'-',75,0)
        DO I = 1,NPPAPT
          INDXPR = LPPAPU(I)
          CALL WRIXPR(I,INDXPR)
        ENDDO
      END IF
C
  999 CONTINUE
      IF (INPERR.GT.0) CALL QUIT('Input error in *EXCITATIONS')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck setxpp */
      SUBROUTINE SETXPP(WORK,KFREE,LFREE)
C*****************************************************************************
C
C     PP excitation module:
C       Initialize variables/Open necessary files etc.
C
C     Written by H.J.Aa.Jensen 2000
C     Based on SETXLR
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0)
C
#include "dgroup.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbxrs.h"
C
      LOGICAL TOBE,ORBDEF(3)
      DIMENSION WORK(*),IR(2,2),NR(2)
      CALL QENTER('SETXPP')
      KFRSAV = KFREE
C
C     ************************************
C     *** Set orbital strings          ***
C     ************************************
C
      CALL OSTRING(XPP_INDSTR,XPP_SKIPEE,XPP_SKIPEP,WORK,KFREE,LFREE)
C
C     ************************************
C     *** Sort operators on symmetries ***
C     ************************************
C
      CALL OP1SRT('LPPAPS',LPPAPU,LPPAPS,NPPAPT,NPPAP,JPPAP,
     &            WORK(KFREE),LFREE)
C
C     **************************************************
C     *** Transfer information from XCBXPP to WRKXRS ***
C     **************************************************
C
      CALL SETRSP
      LINEQ     = .FALSE.
C     ... solve eigenvalue equation!
      FKRMC     = .FALSE.
C     ... not KRMC wave function optimization
      LSVCFG(1) = XPP_LSFG(1)
      LSVCFG(2) = XPP_LSFG(2)
      TKNORM    = XPPNRM
      DIAGHE    = .TRUE.
      IPRXRS    = IPRXPP
      THCXRS    = THCPP
      RESFAC    = RESXPP
      MAXITR    = ITRXPP
      MAXSIM    = MXLOAD
      NREDM     = MAXRM
      N2REDM    = MAXRM*MAXRM
      LOFFTY    = 0
      CNVINT(1) = CNVXPP(1)
      CNVINT(2) = CNVXPP(2)
      ITRINT(1) = ITRIPP(1)
      ITRINT(2) = ITRIPP(2)
      INTDEF    = INTXPP
C     special variables only for linear equations (LINEQ true):
      NFREQ     = 0
      STERNH    = .FALSE.
      STATIC    = .FALSE.
      TRIPLET   = XPP_TRIPLET
      E2CHEK    = XPP_E2CHEK
      IF(E2CHEK) DIAGHE=.FALSE.
C
C     Indicate from vectors are read...
C
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','EIGEXT')
      CALL REACMO(LUBUF,'DFCOEF',DUM,DUM,DUM,DUM,1,DUM,IDUM,IDUM)

      CLOSE(LUBUF,STATUS='KEEP')
C
C     Initialize file for solution vectors
C
      LUXVC = 29
      CALL OPNFIL(LUXVC,'PAMXVC','UNKNOWN','PRPXPP')
      WRITE(LUXVC) 'END_OF_THIS_FILE'
C
C     Look for restart file for solution vectors
C
      RSREST = .FALSE.
      IF(XPPFIL.NE.'NORXVC') THEN
        INQUIRE(FILE=XPPFIL,EXIST=TOBE)
        RSREST = TOBE
      ENDIF
      IF(RSREST) THEN
        LURST = 28
        CALL OPNFIL(LURST,XPPFIL,'OLD','PRPXPP')
      ENDIF
      IF (KFRSAV .NE. KFREE)
     &   CALL MEMREL('SETXPP',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SETXPP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETPPF */
      SUBROUTINE GETPPF(ATMPPF,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &                  EVALR,EVECR,ISSYM,IORBCL,EIG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Evaluate linear response excitation energies and transition moments
C
C     Written by H.J.Aa. Jensen 2000, based on GETLRF
C
C*****************************************************************************
      use orbital_rotation_indices

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,D0=0.0D0)
C
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dcbibn.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      LOGICAL FIRST
      DIMENSION ATMPPF(MAXEXC,NPPAPT,NBSYM,*),RCNV(*),ISSYM(*),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IORBCL(4,*),
     &          IBPVC(*),EVALR(*),EVECR(NREDM,*),EIG(*),WORK(*)
C
      CALL QENTER('GETPPF')
      KFRSAV = KFREE
      FIRST = .TRUE.
      IF(ORBXPP) THEN
        NI = 0
        DO I = 1,NFSYM
          NI = NI + NCLS_XPP(1,1,I) + NCLS_XPP(2,1,I)
        ENDDO
      ELSE
        NI = 0
      ENDIF
      NTRMOM = MAXEXC*NPPAPT
C
C     First evaluate solution vectors for all frequencies
C     and write them to direct access file
C     ===================================================
C
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
C     IFAC  = 1
C     IF(.NOT.STATIC) IFAC = 2
      IFAC  = 2
C     ... exc.en. always not static
C
C     Get coefficients
C     ----------------
C
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','GETPPF')
      CALL MEMGET('INTE',KIBRP,NTBAS(0),WORK,KFREE,LFREE)
      CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),DUM,WORK(KIBRP),
     &              TOTERG,10,WORK,KFREE,LFREE)
      CLOSE(LUBUF,STATUS='KEEP')
C
      IOFF  = 1
C
C     Solution vector - orbital (e-e) part
C     ====================================
C
      IF(NZXOPE.GT.0) THEN
        NESIM = NEXCNV * IFAC
        KEE   = KFREE
        CALL MEMGET('REAL',KGVEC,NZXOPEQ       ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KXVEC,NZXOPEQ*NEXCNV,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBBUF,NZXOPEQ       ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIVEC,NESIM         ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBMAX,NEXCNV        ,WORK,KFREE,LFREE)
!fixme ATMPPF can get out of bounds (run test dft_response)
        CALL GETPP1(ATMPPF(1,1,1,IOFF+1),JBENDX,EVECR,EVALR,RCNV,
     &              IBTYP,IBEVC,FIRST,ISSYM,
     &              WORK(KGVEC),WORK(KXVEC),WORK(KBBUF),NZXOPE,
     &              WORK(KCMO),WORK(KIBRP),IORBCL,EIG,
     &              get_orbital_rotation_indices_pp(),
     &              WORK(KIVEC),WORK(KBMAX),WORK,KFREE,LFREE)
        IF(FIRST) THEN
          CALL DCOPY(NTRMOM,ATMPPF(1,1,JSYMOP,IOFF+1),1,
     &               ATMPPF(1,1,JSYMOP,1),1)
          FIRST = .FALSE.
        ELSE
          CALL DAXPY(NTRMOM,D1,ATMPPF(1,1,JSYMOP,IOFF+1),1,
     &               ATMPPF(1,1,JSYMOP,1),1)
        ENDIF
        CALL MEMREL('GETPPF.e-e',WORK,1,KEE,KFREE,LFREE)
        IOFF = IOFF + NI + 1
      ENDIF
C
C     Solution vector - orbital (e-p) part
C     ====================================
C
      IF(NZXOPP.GT.0) THEN
        NPSIM = NEXCNV * IFAC
        KEP   = KFREE
        CALL MEMGET('REAL',KGVEC,NZXOPPQ       ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KXVEC,NZXOPPQ*NEXCNV,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBBUF,NZXOPPQ       ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIVEC,NPSIM         ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBMAX,NEXCNV        ,WORK,KFREE,LFREE)
        CALL GETPP1(ATMPPF(1,1,1,IOFF+1),JBPNDX,EVECR,EVALR,RCNV,
     &              IBTYP,IBPVC,FIRST,ISSYM,
     &              WORK(KGVEC),WORK(KXVEC),WORK(KBBUF),NZXOPP,
     &              WORK(KCMO),WORK(KIBRP),IORBCL,EIG,
     &              get_orbital_rotation_indices_pn(),
     &              WORK(KIVEC),WORK(KBMAX),WORK,KFREE,LFREE)
        IF(FIRST) THEN
          CALL DCOPY(NTRMOM,ATMPPF(1,1,JSYMOP,IOFF+1),1,
     &               ATMPPF(1,1,JSYMOP,1),1)
          FIRST = .FALSE.
        ELSE
          CALL DAXPY(NTRMOM,D1,ATMPPF(1,1,JSYMOP,IOFF+1),1,
     &               ATMPPF(1,1,JSYMOP,1),1)
        ENDIF
        CALL MEMREL('GETPPF.e-p',WORK,1,KEP,KFREE,LFREE)
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
      CALL MEMREL('GETPPF.end',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('GETPPF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GETPP1 */
      SUBROUTINE GETPP1(ATMPPF,NBTYP,EVECR,EVALR,RCNV,IBTYP,IBVEC,FIRST,
     &                  ISSYM,GPO,XPO,BBUF,NXPAR,CMO,IBEIG,
     &                  IORBCL,EIG,JXOP,IVECS,BMAX,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Get linear response excitation function for orbital rotation parameter
C       contribution from given variational parameter type
C       specified by NBTYP
C
C     Written by H.J.Aa. Jensen 2000, based on GETLR1
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D2 = 2.0D00)
C
#include "dcbibn.h"
#include "dcbxrs.h"
#include "dcbxpp.h"
#include "dcbham.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbxpr.h"
      LOGICAL FIRST,TEST
      CHARACTER FILNAM*6,MEMTYP*6,BLANK16*16,BLANK2*2
      DIMENSION ATMPPF(MAXEXC,NPPAPT,NBSYM,*),RCNV(*),
     &          GPO(NXPAR,NZ),XPO(NXPAR,NZ,MAXEXC),BBUF(*),
     &          CMO(*),EIG(*),IBEIG(*),IORBCL(4,*),ISSYM(*),
     &          IBVEC(*),IBTYP(2,*),EVECR(NREDM,*),EVALR(*),
     &          JXOP(2,*),IVECS(*),BMAX(NEXCNV),WORK(*)
      DATA BLANK16/'                '/
      DATA BLANK2/'  '/
C
#include "chrnos.h"
#include "ibtfun.h"
      CALL QENTER('GETPP1')
      KFRSAV = KFREE
      IF(IPRXPP.GE.3) CALL HEADER('Output from GETPP1',-1)
C
      NBDIM = NXPAR*NZ
      IF(NBTYP.EQ.JBENDX) THEN
        NBRED  = NERED
        LUBVEC = LUBOE
        FILNAM = 'PAMBOE'
        MEMTYP = 'GETPPE'
        ITYP   = IPEP
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBRED  = NPRED
        LUBVEC = LUBOP
        FILNAM = 'PAMBOP'
        MEMTYP = 'GETPPP'
        ITYP   = IPPP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'GETPP1 ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('GETPP1: Unknown NBTYP !')
      ENDIF
C
C
C     =========================
C
      OPEN(LUBVEC,FILE=FILNAM,FORM='UNFORMATTED',
     +         ACCESS='DIRECT',RECL=8*NBDIM,STATUS='UNKNOWN')
      ITIM = IPX
      NXDIM = NBDIM*NEXCNV
C
C     We have two rounds,
C     first with ITIM = IPX = JTIMOP and second with ITIM = IMX
C
      TEST=FIRST.AND.(LINEAR.OR.(SPINFR.AND.(NZ.EQ.1)))
      IF(TEST) CALL DZERO(BMAX,NEXCNV)
      XPOMAX = D0
      IROUND = 0
  100 IROUND = IROUND + 1
      CALL DZERO(XPO,NXDIM)
      CALL XRSXV1(ITIM,NBTYP,XPO,EVECR,NEXCNV,IBTYP,IBVEC,IVECS,BBUF)
C
C     Print section
C     -------------
C
      IF(IPRXPP.GE.5) THEN
        WRITE(LUPRI,'(/1X,A2,A,I3,A)') BVTYP(ITYP),
     &   ' part of solution vectors for',NEXCNV,' excitations:'
        CALL PRBVEC(LUPRI,XPO,NEXCNV,NXPAR)
        CALL FLSHFO(LUPRI)
      ENDIF
C
C     Identify largest element
C     ========================
C
      IF(TEST) THEN
        DO IEXC = 1,NEXCNV
          DO I = 1,NXPAR
            DTEMP = XPO(I,1,IEXC)*XPO(I,1,IEXC)
            DO IZ = 2,NZ
              DTEMP = DTEMP + XPO(I,IZ,IEXC)*XPO(I,IZ,IEXC)
            ENDDO
            IF(DTEMP.GT.BMAX(IEXC)) THEN
              ISSYM(IEXC) = I
              BMAX(IEXC) = DTEMP
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
C     Form transition moments
C     =======================
C
        NOPA   = NPPAP(JSYMOP)
        JOPA   = JPPAP(JSYMOP)
        DO IOPA = 1,NOPA
          INDAP  = LPPAPS(JOPA+IOPA)
          INDPRA = LPPAPU(INDAP)
          KTIMOP = IPRPTIM(INDPRA)
C
C         Get gradient
C         ------------
C
          IF (KTIMOP .EQ. ITIM) THEN
             CALL GPGET(INDPRA,GPO,JXOP,NXPAR,
     &                  CMO,IBEIG,NBTYP,WORK,KFREE,LFREE,IPRXPP)
          END IF
C
C         Form transition moment array
C         ----------------------------
C
          DO IEXC = 1,NEXCNV
            IF (KTIMOP .EQ. ITIM) THEN
              ATMPPF(IEXC,INDAP,JSYMOP,1)  =
     &        D2*DDOT(NBDIM,GPO,1,XPO(1,1,IEXC),1)
              IF(ORBXPP) THEN
                DO IZ = 1,NZ
                  DO I = 1,NXPAR
                    I1 = JXOP(1,I)
                    IND = IORBCL(2,I1) + 1
                    ATMPPF(IEXC,INDAP,JSYMOP,IND)=
     &              ATMPPF(IEXC,INDAP,JSYMOP,IND)+
     &                D2*GPO(I,IZ)*XPO(I,IZ,IEXC)
                  ENDDO
                ENDDO
              ENDIF
              IF(IPRXPP.GE.3) THEN
                FREQ = EVALR(IEXC)
Chj: MAERKE change to trans.mom. output
                WRITE(LUPRI,1000)
     &          '< 0 |',INDAP,' | ',IEXC,' > =',
     &          ATMPPF(IEXC,INDAP,JSYMOP,1),
     &          ' a.u.,',BVTYP(ITYP),
     &          ' contribution, excitation energy',FREQ,' a.u.'
              ENDIF
            ENDIF
          ENDDO
        ENDDO
C
C
C     Write solution vectors to unformatted file for restart
C
      DO IEXC = 1,NEXCNV
        FREQ = EVALR(IEXC)
        CALL WRTRSP(LUXVC,XPO(1,1,IEXC),JXOP,
     &       RSPLAB,BVTYP(ITYP),FREQ,JSYMOP,JTIMOP,
     &       BLANK16,BLANK2,D0,0,0,RCNV(IEXC),INTXRS,NBDIM)
      ENDDO
C
C     Now calculate the part of the excitation
C     solution vector which is of -JTIMOP symmetry,
C     calculate contribution to properties and write them to file also.
C
      IF(IROUND.EQ.1) THEN
        ITIM = IMX
        ITYP = ITYP + 1
        GO TO 100
      END IF
      CLOSE(LUBVEC,STATUS = 'KEEP')
C
      IF(TEST) THEN
        IREP = JSYMOP-1
        IFRP = JBTOF(IREP,1)
        IF(LINEAR) THEN
          IF(IREP.EQ.JSPINR(1,1,IFRP).OR.
     &       IREP.EQ.JSPINR(2,1,IFRP)) THEN
C.........excitation from unbarred to unbarred (or from barred to barred);
C.........Total mj of excited state is therefore mj(virtual)-mj(inactive)
            IFAC = -1
          ELSEIF(IREP.EQ.JSPINR(3,1,IFRP).OR.
     &           IREP.EQ.JSPINR(4,1,IFRP)) THEN
C.........excitation from unbarred to barred (or vice versa);
C.........Total mj of excited state is therefore mj(virtual)+mj(inactive)
            IFAC = +1
          ELSE
            WRITE(6,'(A,A3)') 
     &      'XPPsupersymmetry: something went wrong! Irrep: ',REP(IREP)
          ENDIF
          DO IEXC = 1,NEXCNV
            I = ISSYM(IEXC)
            INDI = JXOP(1,I)
            INDA = JXOP(2,I)
            ISSYM(IEXC) = IBEIG(INDA) + IFAC*IBEIG(INDI)
          ENDDO
        ELSEIF(SPINFR.AND.(NZ.EQ.1)) THEN
          DO IEXC = 1,NEXCNV
            I = ISSYM(IEXC)
            INDI = JXOP(1,I)
            INDA = JXOP(2,I)
            ISSYM(IEXC) = IBTXOR(IBEIG(INDA),IBEIG(INDI))
          ENDDO
        ENDIF
      ENDIF      
C
      CALL QEXIT('GETPP1')
      RETURN
 1000 FORMAT(/A,I2,A,I2,A,F16.10,A,/3X,A2,A,F16.10,A)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XPPOUT_old */
      SUBROUTINE XPPOUT_OLD(JSYMTEST,ATMPPF,OMEGA,RCNV,EIG,ISSYM,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Final output from linear response excitation energy module
C
C     Written by H.J.Aa. Jensen 2000, based on XLROUT
C
C*****************************************************************************
         use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      PARAMETER ( D0 = 0.0D0, D2R3 = (2.0D0/3.0D0), D2 = 2.0D0 )
#include "dgroup.h"
#include "pgroup.h"
#include "dcbham.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
#include "symmet.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dcbxrs.h"
#include "dcbxpr.h"
C
      CHARACTER MXFORM*6,PFMT*15,TTYP(-1:1)*2,SYMLAB*5
      DIMENSION ATMPPF(MAXEXC,NPPAPT,NBSYM,*),OMEGA(MAXEXC,NBSYM),
     &          EIG(*),RCNV(MAXEXC,NBSYM),IBUF(2),ISSYM(MAXEXC,NBSYM),
     &          WORK(*)
      DATA TTYP /'T-','T0','T+'/
      DATA SYMLAB/'     '/
#include "ibtfun.h"
C
      IF(ORBXPP) THEN
        NI = 0
C       Number of inactive orbitals
CTROND  Check out active ones...
        DO I = 1,NFSYM
          NI = NI + NCLS_XPP(1,1,I) + NCLS_XPP(2,1,I)
        ENDDO
        IOFF = 1
        NVAR = 0
        IF(.NOT.XPP_SKIPEE) THEN
          NVAR = NVAR + 1
          IBUF(NVAR)  = IOFF+1
          IOFF = IOFF + NI + 1
        ENDIF
        IF(.NOT.XPP_SKIPEP) THEN
          NVAR = NVAR + 1
          IBUF(NVAR)  = IOFF+1
          IOFF = IOFF + NI + 1
        ENDIF
      ENDIF
C
C     Loop over transition operators
C     ==============================
C
      CALL HEADER('Linear response excitations; final output ',-1)
      DO IOPSY = 1,NFSYM
        NBF = NBSYM/NFSYM
        IF(LINEAR.AND.NFSYM.EQ.2) THEN
          IF(IOPSY.EQ.1) THEN
            WRITE(SYMLAB(4:4),'(A1)') 'g'
          ELSE
            WRITE(SYMLAB(4:4),'(A1)') 'u'
          ENDIF
        ENDIF
        DO IBF = 1,NBF
          ISYM = JFSYM(IBF,IOPSY)
          IREP = ISYM-1
          IF (JSYMTEST .GT. 0) THEN
             IF (ISYM.NE.JSYMTEST) GO TO 500
          END IF
          NOPA = NPPAP(ISYM)
          JOPA = JPPAP(ISYM)
          IF (KEXCNV(ISYM).GT.0) THEN
             CALL PRSYMB(LUPRI,'=',80,0)
             WRITE(LUPRI,'(/A,I2,A,A3/)')
     &     ' *** Excitations of boson symmetry',ISYM,' : ',REP(IREP)
          END IF
          DO I = 1,KEXCNV(ISYM)
             CALL PRSYMB(LUPRI,'-',80,0)
             IF(RCNV(I,ISYM).LT.THCPP) THEN
               WRITE(LUPRI,'(A,I4,A,F15.8,A,3X,1P,E9.2,A)')
     &            ' Excitation no.',I,' excitation energy',
     &            OMEGA(I,ISYM),' a.u.',RCNV(I,ISYM),' (converged)'
             ELSE
               WRITE(LUPRI,'(A,I4,A,F15.8,A,3X,1P,E9.2,A)')
     &            ' Excitation no.',I,' excitation energy',
     &            OMEGA(I,ISYM),' a.u.',
     &            RCNV(I,ISYM),' (**NOT** converged!)'
             ENDIF
             IF(LINEAR) THEN
               WRITE(SYMLAB(1:3),'(I3)') ABS(ISSYM(I,ISYM))/2
               WRITE(LUPRI,'(0P,A17,A5,F30.4,A)')
     &           ' Linear symmetry:',SYMLAB,
     &            OMEGA(I,ISYM)*XTEV,' eV'
             ELSEIF(SPINFR.AND.(NZ.EQ.1)) THEN
               ISPACE = ISSYM(I,ISYM)
               ISPIN = IBTXOR(ISPACE,IREP)
               WRITE(SYMLAB(3:5),'(A3)') REP(ISPACE)
               IF(ISPIN.EQ.0) THEN
                 WRITE(SYMLAB(1:2),'(A2)') 'S0'
               ELSEIF(ISPIN.EQ.ISYMAX(3,2)) THEN
                 WRITE(SYMLAB(1:2),'(A2)') 'Tz'
               ELSEIF(ISPIN.EQ.ISYMAX(2,2)) THEN
                 WRITE(SYMLAB(1:2),'(A2)') 'Ty'
               ELSEIF(ISPIN.EQ.ISYMAX(1,2)) THEN
                 WRITE(SYMLAB(1:2),'(A2)') 'Tx'
               ELSE
                 WRITE(6,*)
     &            'ERROR: Spin-space symmetry failed',ISPIN,ISPACE
               ENDIF
               WRITE(LUPRI,'(0P,A17,A5,F30.4,A)')
     &           ' Nonrel. sym.:   ',SYMLAB,
     &            OMEGA(I,ISYM)*XTEV,' eV'
             ELSE
               WRITE(LUPRI,'(0P,22X,F30.4,A)')
     &            OMEGA(I,ISYM)*XTEV,' eV'
             ENDIF

          IF (NOPA .EQ. 0) GO TO 300
             WRITE(LUPRI,'(A/)')
     &           ' *** Transition moments <0|A|n> for this excitation:'
            OSCILL = D0
            OSCILV = D0
            DO IOPA = 1,NOPA
               INDAP  = LPPAPS(JOPA+IOPA)
               INDPRA = LPPAPU(INDAP)
               IREPA  = IPRPSYM(INDPRA)-1
               ITIMA  = IPRPTIM(INDPRA)
               PFMT = MXFORM(ATMPPF(I,INDAP,ISYM,1),15)
               WRITE(LUPRI,'(A,A16,5X,A3,2X,A2,8X,'//PFMT//',A)')
     &             ' A - ',PRPNAM(INDPRA),REP(IREPA),TTYP(ITIMA),
     &             ATMPPF(I,INDAP,ISYM,1),' a.u.'
               IF (INDEX(PRPNAM(INDPRA),'DIPLEN') .NE. 0) THEN
                   XMOM = ATMPPF(I,INDAP,ISYM,1)
                   OSCILL = OSCILL + D2R3*XMOM*XMOM*OMEGA(I,ISYM)
               END IF
               IF (INDEX(PRPNAM(INDPRA),'c_alpha') .NE. 0  .OR.
     &             INDEX(PRPNAM(INDPRA),'C_ALPHA') .NE. 0) THEN
                   XMOM = ATMPPF(I,INDAP,ISYM,1)
                   OSCILV = OSCILV + D2R3*XMOM*XMOM/OMEGA(I,ISYM)
               END IF
C
               IF (ORBXPP) THEN
                  WRITE(LUPRI,'(/3X,A)') '* Orbital analysis:'
                  IND = 0
                  DO IFRP = 1,NFSYM
C - positronic orbitals
                    DO J = 1,NCLS_XPP(1,1,IFRP)
                      IND = IND + 1
                      A = D0
                      DO K = 1,NVAR
                        A = A + ATMPPF(I,INDAP,ISYM,IBUF(K)+IND)
                      ENDDO
                      WRITE(LUPRI,'(3X,A3,A1,4(2X,F16.8))')
     &          FREP(IFRP),'p',EIG(IND),A,
     &          (ATMPPF(I,INDAP,ISYM,IBUF(K)+IND),K=1,NVAR)
                    ENDDO
C - electronic orbitals
                    DO J = 1,NCLS_XPP(2,1,IFRP)
                      IND = IND + 1
                      A = D0
                      DO K = 1,NVAR
                        A = A + ATMPPF(I,INDAP,ISYM,IBUF(K)+IND)
                      ENDDO
                      WRITE(LUPRI,'(3X,A3,A1,4(2X,F16.8))')
     &       FREP(IFRP),'e',EIG(IND),A,
     &       (ATMPPF(I,INDAP,ISYM,IBUF(K)+IND),K=1,NVAR)
                    ENDDO
                  ENDDO
               ENDIF
            ENDDO
C           end do iopa = 1,nopa
            OSCIL = MAX(OSCILL,OSCILV)
            IF (OSCIL .GT. 1.0D-10) THEN
C
C              Dipole oscillator strengths
C
               IF (OSCIL .GE. 1.0D-3) THEN
                   PFMT = '(5X,A,F12.8)'
                ELSE
                   PFMT = '(5X,A,1P,E16.8)'
                END IF
                WRITE(LUPRI,'()')
                WRITE(LUPRI,PFMT)
     &             '-> dipole oscillator strength (length)   :',OSCILL
                WRITE(LUPRI,PFMT)
     &             '-> dipole oscillator strength (velocity) :',OSCILV
c
c hj TODO : mixed require code with XMOML(3), XMOMV(3) to be sure the
c     right components are combined !!!!
c              IF (IDIPL .EQ. 1 .AND. IDIPV .EQ. 1) THEN
c                  OSCILM = OSCILM + D2R3*ABS(XMOML*XMOMV)
c              END IF
c               IF (OSCILM .GT. 1.0D-10) WRITE(LUPRI,PFMT)
c    &             '-> dipole oscillator strength (mixed)    :',OSCILM
C
C              Dipole radiation rates
C
c              WRITE(LUPRI,'(/A/A/2(A,I2)/)')
c    &            ' Phosphorescence transition rates ',
c    &            ' -------------------------------- ',
c    &            ' Symmetry ',MSYMC, ' Root ', I
c              RATE=0D0
c              DO I=1,3
c                 XDIP=0D0
c                 DO JHSO=1,3
c                    XDIP=XDIP+PHOSMAT(I,JHSO,I)
c                 END DO
c                 XRATE = 4*ALPHAC**3*omega**3/3*XDIPL**2/AUTIME
c                 WRITE(LUW6,'(2A,G12.5)')
c    &               ' Partial rates: ',
c    &               CHAR(ICHAR('X')+I-1)//'-polarization     ',
c    &               XRATE
c                 RATE=RATE+XRATE
c              END DO
c              WRITE(LUPRI,'(/A,G12.5,A)')
c    &            ' Total   rate                 ',  RATE, 's-1'
c              WRITE(LUPRI,'(A,G12.5,A)')
c    &            ' Total   lifetime             ',  1/RATE, 's'
               IF (OSCILL .GT. 1.0D-10) THEN
                  RATE = D2*ALPHAC**3*OSCILL*OMEGA(I,ISYM)**2/AUTIME
                  WRITE(LUPRI,'(2(/5X,A,1P,G12.5,A))')
     &            'Total dipole radiation rate (length)  ',RATE,' s-1',
     &            ' - corresponding radiation life time  ',1/RATE,' s'
               END IF
               IF (OSCILV .GT. 1.0D-10) THEN
                  RATE = D2*ALPHAC**3*OSCILV*OMEGA(I,ISYM)**2/AUTIME
                  WRITE(LUPRI,'(2(/5X,A,1P,G12.5,A))')
     &            'Total dipole radiation rate (velocity)',RATE,' s-1',
     &            ' - corresponding radiation life time  ',1/RATE,' s'
               END IF
c              IF (OSCILM .GT. 1.0D-10) THEN
c                 RATE = D2*ALPHAC**3*OSCILM*OMEGA(I,ISYM)**2/AUTIME
c                 WRITE(LUPRI,'(2(/5X,A,1P,G12.5,A))')
c    &            'Total dipole radiation rate (mixed)   ',RATE,' s-1',
c    &            ' - corresponding radiation life time  ',1/RATE,' s'
c              END IF
            END IF
C
  300     CONTINUE
          ENDDO
C         end do I = 1,KEXCNV(ISYM)
  500     CONTINUE
        ENDDO
      ENDDO
      CALL PRSYMB(LUPRI,'=',80,0)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Eigred */
      SUBROUTINE EIGRED(NEX,N,LDA,AE,AS,W,X,IERR,IPRINT,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     This subroutine solves the reduced eigenvalue equation
C
C       AE*X = W*AS*X
C
C     for the N lowest excitations
C
C     Written by T.Saue Oct 16 1996
C     Rewritten Dec. 2000 by H.J.Aa.Jensen
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D0 = 0.0D0)
C
C     DIMENSION AE(LDA,N),AS(LDA,N),W(N),X(LDA,N),
      DIMENSION AE(*),AS(*),W(N),X(*),
     &          WORK(*),ISNDX(3)
C
#include "dcbxrs.h"
C
      CALL QENTER('EIGRED')
C
C     Memory allocation
C
      KREL = KFREE
      N2A = LDA*N
      CALL MEMGET('REAL',KAE  ,N2A,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAS  ,N2A,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KALFR,N  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KALFI,N  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBETA,N  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIGV,N2A,WORK,KFREE,LFREE)
C
C     Transfer information
C
      CALL DCOPY(N2A,AE,1,WORK(KAE),1)
      CALL DCOPY(N2A,AS,1,WORK(KAS),1)
C
C     Solve generalized eigenvalue problem
C
      MATZ = 1
      CALL RGG(LDA,N,WORK(KAE),WORK(KAS),WORK(KALFR),WORK(KALFI),
     &         WORK(KBETA),MATZ,WORK(KEIGV),IERR)
C     IERR IS SET TO
C       ZERO       for normal return,
C       J          if the limit of 30*N iterations is exhausted
C                     while the J-th eigenvalue is being sought.
C
      IF (IERR.NE.0) THEN
         WRITE(LUPRI,*)
     &      'FATAL ERROR: EIGRED error after RGG, IERR =',IERR
         IPRINT = 999
      END IF
C
      IF (IPRINT .GT. 10) THEN
         WRITE(LUPRI,'(/A/A/,(3F20.10))')
     &   ' EIGRED: (ALFAR, ALFAI, BETA)',
     &   '         where (ALFAR + i ALFAI) / BETA are the eigenvalues;',
     &   (WORK(KALFR+I),WORK(KALFI+I),WORK(KBETA+I),I=0,N-1)
C
         WRITE(LUPRI,'(/A)') ' and the eigenvectors'
         CALL OUTPUT(WORK(KEIGV),1,N,1,N,LDA,N,1,LUPRI)
      END IF
      IF (IERR .GT. 0) CALL QUIT('EIGRED for excitation energies!')
C
C     Order eigenvalues and vectors, select those with positive norm,
C     and check for pairing.
C
      CALL XPPORD(LDA,N,AS,WORK(KEIGV),WORK(KALFR),WORK(KALFI),
     &            WORK(KBETA),IPRINT,WORK(KAE),WORK(KAS),ISNDX)
C
C     Transfer requested NEX electronic exitations eigenvalues and
C     vectors to W and X (skip the NPRED positronic "de-excitations).
C
      JPRED = NPRED/2
      JALFR = KALFR + JPRED
      JEIGV = KEIGV + JPRED*LDA
      CALL DCOPY(    NEX,WORK(JALFR),1,W,1)
      CALL DCOPY(LDA*NEX,WORK(JEIGV),1,X,1)
C
      CALL MEMREL('EIGRED',WORK,1,KREL,KFREE,LFREE)
      CALL QEXIT ('EIGRED')
      RETURN
      END
C  /* Deck xppord */
      SUBROUTINE XPPORD(LDA,KZYRED,SRED,EIVEC,ALFAR,ALFAI,BETA,
     &                  IPRINT,WRK1,WRK2,ISNDX)
C
C Analyze and order eigenvectors from EIGRED
C
C 13. Dec. 2000 HJAaJ
C Based on RSPORD (JO 1984-11-8)
C
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION SRED(LDA,LDA), EIVEC(LDA,LDA)
      DIMENSION WRK1(KZYRED,KZYRED), WRK2(KZYRED,KZYRED)
      DIMENSION ALFAR(KZYRED),ALFAI(KZYRED),BETA(KZYRED), ISNDX(3)
C
C eigenvalue k is (ALFAR(k)+i*ALFAI(k))/BETA(k)
C
C
      PARAMETER (COMPLX = 1.0D7 )
      PARAMETER (D0=0.0D0, D2=2.0D0 )
Chj jan 2001> accept all with non-zero norm;
C             we may get a very small norm
C             when XPPNRM false
      PARAMETER (ZEROT=1.0D-14, PAIRT = 1.0D-8)
C     This factor is connected to not normalizing excitation operators
      ROOT2 = SQRT( D2 )
C
      IF (IPRINT .GT. 5) THEN
         WRITE (LUPRI,'(//A/A)')
     *      '  (ALFAR + i ALFAI) / BETA are eigenvalues;',
     *      '     ALFAR           ALFAI          BETA'
         WRITE (LUPRI,'(1P,3D15.6)')
     *      (ALFAR(I),ALFAI(I),BETA(I),I=1,KZYRED)
      END IF
      ALFAMX = D0
      DO 10 I=1,KZYRED
         IF(ABS(BETA(I)).GT.ZEROT) THEN
            ALFAR(I)=ALFAR(I)/BETA(I)
            ALFAI(I)=ALFAI(I)/BETA(I)
            ALFAMX = MAX( ALFAMX, ABS(ALFAR(I)) )
         ELSE
C           singularities
            WRITE(LUPRI,1010)I, ALFAR(I),ALFAI(I),BETA(I)
         END IF
         IF(ALFAI(I).NE.D0) THEN
            WRITE(LUPRI,1020) I,ALFAR(I),ALFAI(I)
            WRITE(LUPRI,'(A,E7.1)')
     &   'The corresponding eigenvalue is set to :',COMPLX
C
C SET EIGENVALUE EQUAL TO a very large factor 
C COMPLX IN ORDER TO BE ABLE TO SKIP
C CONTRIBUTIONS FROM THIS ROOT WHEN SUMMING UP TERMS
C IN THE CALCULATION OF THE EFFECTIVE SPECTRUM IN C6 CALCULATIONS
C
            ALFAR(I) = COMPLX
         END IF
  10  CONTINUE
 1010 FORMAT(/' WARNING *** Singularity in reduced PP matrix:'
     &      ,/'     I,ALFAR(I),ALFAI(I),BETA(I):',I6,1P,3D13.6)
 1020 FORMAT(/' WARNING *** Complex eigenvalue in reduced PP matrix :'
     &      ,/'     Real and imaginary part :',I8,1P,2D20.8)
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI,'(/A)') ' Eigenvalues of E(2)  :'
         WRITE (LUPRI,'(I10,1P,D12.2)') (I,ALFAR(I),I=1,KZYRED)
      END IF
C
C     reduced S(2) in eigenvector basis
C
      CALL MPAB (SRED, KZYRED,KZYRED,LDA,LDA,
     +           EIVEC,KZYRED,KZYRED,LDA,LDA,
     +           WRK1, KZYRED,KZYRED)
      CALL MPATB(EIVEC,KZYRED,KZYRED,LDA,LDA,
     +           WRK1, KZYRED,KZYRED,KZYRED,KZYRED,
     +           WRK2, KZYRED,KZYRED)
      IF (IPRINT .GT. 10) THEN
         WRITE (LUPRI,'(/A)')
     &      ' Diagonal of reduced S(2) in eigenvector basis :'
         WRITE (LUPRI,'(I10,1P,D12.2)') (I,WRK2(I,I),I=1,KZYRED)
         CALL OUTPUT(WRK2,1,KZYRED,1,KZYRED,KZYRED,KZYRED,1,LUPRI)
      END IF
C
C     Select eigenvectors with positive norm and store them
C     as the first eigenvectors.
C
      IPOS = 0
      IZER = 0
      INEG = 0
      DO 20 I=1,KZYRED
         IF (ABS(BETA(I)) .LE. ZEROT) THEN
            IZER = IZER + 1
         ELSE IF (WRK2(I,I) .GT. D0) THEN
            IPOS = IPOS + 1
            SCALE= ROOT2/SQRT(WRK2(I,I))
            CALL DSCAL(KZYRED,SCALE,EIVEC(1,I),1)
            IF (IPOS.NE.I) THEN
               CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,IPOS),1)
               XSAVE = ALFAR(IPOS)
               ALFAR(IPOS)=ALFAR(I)
               ALFAR(I)   =XSAVE
               XSAVE = WRK2(IPOS,IPOS)
               WRK2(IPOS,IPOS) = WRK2(I,I)
               WRK2(I,I)       = XSAVE
               XSAVE = BETA(IPOS)
               BETA(IPOS) = BETA(I)
               BETA(I) = XSAVE
            END IF
         ELSE IF (WRK2(I,I) .LT. -D0) THEN
            INEG = INEG + 1
            SCALE= ROOT2/SQRT(ABS(WRK2(I,I)))
            CALL DSCAL(KZYRED,SCALE,EIVEC(1,I),1)
            IF (ALFAR(I) .EQ. COMPLX) ALFAR(I) = -COMPLX
         END IF
   20 CONTINUE
      NNEG = IPOS
      DO 22 I=IPOS+1,KZYRED
         IF (ABS(BETA(I)) .GT. ZEROT) THEN
            NNEG = NNEG + 1
            IF (NNEG.NE.I) THEN
               CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,NNEG),1)
               XSAVE = ALFAR(NNEG)
               ALFAR(NNEG)=ALFAR(I)
               ALFAR(I)   =XSAVE
               XSAVE = WRK2(NNEG,NNEG)
               WRK2(NNEG,NNEG) = WRK2(I,I)
               WRK2(I,I)       = XSAVE
            END IF
         END IF
   22 CONTINUE
      ISNDX(1) = IPOS
      ISNDX(2) = IZER
      ISNDX(3) = INEG
      IF (IPOS.NE.INEG .OR. IZER.NE.0)
     *   WRITE (LUPRI,2020) IPOS,IZER,INEG
 2020 FORMAT(/' *** EIGENVECTORS WITH POSITIVE METRIC:',I6,
     *       /'     EIGENVECTORS WITH ZERO METRIC:    ',I6,
     *       /'     EIGENVECTORS WITH NEGATIVE METRIC:',I6)
C
C Order eigensolutions in ascending order of eigenvalues.
C.....first the ones with positive metric
      DO 100 I=1,IPOS
         JMIN = I
         AMIN = ALFAR(I)
         DO 90 J=I+1,IPOS
            IF(ALFAR(J).LT.AMIN) THEN
               AMIN = ALFAR(J)
               JMIN = J
            ENDIF
   90    CONTINUE
         IF (JMIN.NE.I) THEN
            ALFAR(JMIN)=ALFAR(I)
            ALFAR(I)=AMIN
            CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,JMIN),1)
         ENDIF
  100 CONTINUE
      INOPAIR = 0
      TSTPAIR = ALFAMX*PAIRT
C.....next the ones with negative metric
      DO 110 I=IPOS+1,IPOS+INEG
         JMIN = I
         AMIN = ALFAR(I)
         DO 92 J=I+1,IPOS+INEG
            IF(ALFAR(J).GT.AMIN) THEN
               AMIN = ALFAR(J)
               JMIN = J
            ENDIF
   92    CONTINUE
         IF (JMIN.NE.I) THEN
            ALFAR(JMIN)=ALFAR(I)
            ALFAR(I)=AMIN
            CALL DSWAP(KZYRED,EIVEC(1,I),1,EIVEC(1,JMIN),1)
         ENDIF
         IF (ABS(ALFAR(I)+ALFAR(I-IPOS)).GT.TSTPAIR)
     &      INOPAIR = INOPAIR + 1
  110 CONTINUE
      IF (INEG.NE.IPOS) THEN
         WRITE(LUPRI,'(/3A)')' **WARNING** XPPORD:'
     *   ,' number of eigenvalues with negative metric differ from'
     *   ,' number with positive metric'
         WRITE(LUPRI,'(/A)') '   NUMBER    EIGENVALUE '
         DO 117 I=1,KZYRED
            WRITE(LUPRI,'(I10,1P,D20.8)') I,ALFAR(I)
 117     CONTINUE
      ELSE
         IF (INOPAIR.GT.0) WRITE(LUPRI,'(A/)')
     &      ' **WARNING** XPPORD: Eigenvalues not paired'
         IF (IPRINT.GT.5 .OR. INOPAIR .GT.0) WRITE(LUPRI,'(/A)')
     *      '      NUMBER    EIGENVALUE       PAIRED EIGENVALUE'
         DO 115 I=1,IPOS
            IF (IPRINT.GT.5 .OR. INOPAIR .GT.0)
     *         WRITE(LUPRI,'(I10,1P,2D20.8)')
     *                      I,ALFAR(I),ALFAR(IPOS+I)
 115     CONTINUE
         IF (IZER .GT. 0) THEN
            WRITE(LUPRI,'(/A/A)')
     *         ' **WARNING** XPPORD: Zero metric eigenvalue(s)',
     *         '      NUMBER    EIGENVALUE (real part)'
            DO 116 I=1,IZER
               WRITE(LUPRI,'(I10,1P,D20.8)')
     *                     I,ALFAR(IPOS+INEG+I)
 116        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XPPOUT */
      SUBROUTINE XPPOUT(JSYMTEST,ATMPPF,OMEGA,RCNV,EIG,ISSYM,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C
C***********************************************************************
         use codata
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
      PARAMETER ( D0 = 0.0D0, D2R3 = (2.0D0/3.0D0), D2 = 2.0D0 )
#include "dgroup.h"
#include "pgroup.h"
#include "dcbham.h"
#include "dcbxpp.h"
#include "dcbxlr.h"
#include "symmet.h"
Chj: dcbxlr includes some general information
Chj: which should be moved to dcbxrs (?) or another one:
Chj: MAXRM, MXLOAD
#include "dcbxrs.h"
#include "dcbxpr.h"
C
      CHARACTER MXFORM*6,PFMT*15,TTYP(-1:1)*2
      DIMENSION ATMPPF(MAXEXC,NPPAPT,NBSYM,*),OMEGA(MAXEXC,NBSYM),
     &          EIG(*),RCNV(MAXEXC,NBSYM),IBUF(2),ISSYM(MAXEXC,NBSYM),
     &          WORK(*)
      character (len=4),  allocatable   :: symlab(:)
      real(8), allocatable   :: excerg(:)
      integer, allocatable :: indx(:),irsp(:),ind(:)
      DATA TTYP /'T-','T0','T+'/
#include "ibtfun.h"
      CALL HEADER('Linear response excitations; final output ',-1)
      NEXC = ISUM(NBSYM,KEXCNV,1)
      allocate(symlab(nbsym))
      allocate(excerg(nexc))
      allocate(indx(nexc))
      allocate(irsp(nexc))
      allocate(ind(nbsym))
      IOFF = 2
      DO ISYM = 1,NBSYM
        SYMLAB(ISYM)=REP(ISYM-1)//' '
        MEXC=KEXCNV(ISYM)
        CALL DCOPY(MEXC,OMEGA(1,ISYM),1,EXCERG(IOFF),1)
        IOFF = IOFF + MEXC
      ENDDO
      EXCERG(1) = 0.0D0
      KEXCNV(1) = KEXCNV(1) + 1
      CALL PRT_EV(LUPRI,NBSYM,KEXCNV,SYMLAB,0.0D0,0.0D0,1,EXCERG,NEXC,
     &            INDX,IRPS,IND)
      deallocate(symlab)
      deallocate(excerg)
      deallocate(indx)
      deallocate(irsp)
      deallocate(ind)
C
      RETURN
      END

