C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsctl */
      SUBROUTINE XRSCTL(GPCI,GPOE,GPOP,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &                  EVALR,EVECR,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Module for solving the linear response equations.
C     Based on RSPCTL from RESPONS.
C
C     
C     LINEQ.EQ.FALSE  direct the solution the generalized response
C     eigenvalue problem
C
C        (E[2] - w(I)S[2]) X(I) = 0
C
C     LINEQ.EQ.TRUE   direct the solution the generalized response
C     linear equations
C
C        (E[2] - w(I)S[2]) X(I) = E[1]_B
C
C        w(I)     is frequency no. I
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C     
C     where h is the operator manifold. In RPA h is limited to single 
C     excitations.
C
C-----------------------------------------------------------------------------
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the corresponding reduced equations
C
C        (RE[2] - wRS[2])a = 0            (LINEQ.EQ.FALSE.)
C
C        (RE[2] - wRS[2])a = RE[1]_B      (LINEQ.EQ.TRUE.)
C
C     where
C
C        RE[2] = B^{+}E[2]B
C        RS[2] = B^{+}S[2]B
C        RE[1] = B^{+}E[1}
C    
C-----------------------------------------------------------------------------
C
C     Variables:
C        GPCI    - property gradient , configurational part
C        GPOE    - property gradient , orbital part(e-e rotations)
C        GPOP    - property gradient , orbital part(e-p rotations)
C        REDGP   - reduced property gradient
C        REDE    - reduced Hessian matrix E[2]
C        REDS    - reduced metric         S[2]
C        IBTYP   - trial vectors will be of three kinds:
C                  C-type: contains only configurational parameters
C                  E-type: contains only e-e rotations
C                  P-type: contains only e-p rotations
C        IBCVC   - pointer from C-type trial vectors to IBTYP
C        IBEVC   - pointer from E-type trial vectors to IBTYP
C        IBPVC   - pointer from P-type trial vectors to IBTYP
C        EVALR   - eigenvalues of reduced system
C        EVECR   - eigenvectors of reduced system
C        
C     Files:
C
C        LUBCI   - direct access file of C-type trial vectors
C        LUBOE   - direct access file of E-type trial vectors
C        LUBOP   - direct access file of B-tpye trial vectors
C
C        LUSCI   - direct access file of C-part of sigma vectors
C        LUSOE   - direct access file of E-part of sigma vectors
C        LUSOP   - direct access file of P-part of sigma vectors
C  
C        LURSP   - unformatted file with restart information
C-----------------------------------------------------------------------------
C
C     In the frequency-dependent case, the eigenvectors come in pairs
C
C     b_1^T = [Z  Y^*]        b_2^T = [Y  Z^*]
C
C     The paired structure is exploited in that only the upper half
C     (Z and Y respectively) of the pair of trial vectors is stored.
C
C-----------------------------------------------------------------------------
C
C NEXSTV: number of start vectors
C NEXSIM: desired number of simultaneous roots in a microiteration
C In each microiteration the reduced matrices are extended by NTSIM
C trial vectors.
C    NTSIM = NCSIM + NESIM + NPSIM
C            NCSIM - number of configurational trial vectors
C            NESIM - number of orbital trial vectors (e-e rotations)
C            NPSIM - number of orbital trial vectors (e-p rotations)
C
C IBTYP is an array running over ALL trial vectors:
C   IBTYP(1,*) - indicates type of this trial vector
C   IBTYP(2,*) - indicates Hermiticity of trial vector
C                 1 - [Z+Y  Y*+Z*] - anti-Hermitian combination
C                 2 - [Z-Y  Y*-Z*] - Hermitian combination
C                    
C  ( TRIAL VECTORS ARE ONLY ADDED FOR THE NON CONVERGED VECTORS AND
C    LINEAR DEPENDENCE IS REMOVED)
C
C MAXRIT: MAXIMUM NUMBER OF MICROITERATIONS
C
C     Written by Hans Joergen Aa. Jensen and T.Saue 1996
C     Last revision: July 30 1996 - tsaue
C  
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
C
      DIMENSION GPCI(*), GPOE(*),GPOP(*),
     &          IBTYP(2,*), IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(*),WORK(*)         
C
#include <dcbgen.h>
#include <dcbxrs.h>
#include <blocks.h>
C
      CALL QENTER('XRSCTL')
      KFRSAV = KFREE
      IF(LINEQ) THEN
        NREDGP = NREDM
      ELSE
        NREDGP = 0
      ENDIF
      CALL MEMGET('REAL',KREDE ,N2REDM ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDS ,N2REDM ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KREDGP,NREDGP,WORK,KFREE,LFREE)
      NPOS = 0
      IF(PARCAL) THEN
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NLRGBL*(NLRGBL+1)/2
        NPOS = NPOS + NSMLBL*(NSMLBL+1)/2
      ENDIF
      CALL MEMGET('INTE',KPOS,NPOS,WORK,KFREE,LFREE)
C
      CALL MEMGET('REAL',KRNORM,7*NFREQ,WORK,KFREE,LFREE)
C
      CALL XRSCT1(GPCI,GPOE,GPOP,RCNV,
     &            IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &            WORK(KREDE),WORK(KREDS),WORK(KREDGP),
     &            WORK(KRNORM),WORK(KPOS),
     &            WORK,KFREE,LFREE)
C
      CALL MEMREL('XRSCTL.XRSRED',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XRSCTL')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSCT1 */
      SUBROUTINE XRSCT1(GPCI,GPOE,GPOP,RCNV,
     &                  IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &                  REDE,REDS,REDGP,RNORM,
     &                  NPOS,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Module for solving the linear response equations.
C     Based on RSPCTL from RESPONS.
C
C     
C     LINEQ.EQ.FALSE  direct the solution the generalized response
C     eigenvalue problem
C
C        (E[2] - w(I)S[2]) X(I) = 0
C
C     LINEQ.EQ.TRUE   direct the solution the generalized response
C     linear equations
C
C        (E[2] - w(I)S[2]) X(I) = E[1]_B
C
C        w(I)     is frequency no. I
C        E[2]     is the molecular Hessian                (h|H|h)
C        S[2]     is a metric                             (h|h)
C        E[1}_B   is gradient of perturbation operator B  (h|B)
C     
C     where h is the operator manifold. In RPA h is limited to single 
C     excitations.
C
C-----------------------------------------------------------------------------
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the corresponding reduced equations
C
C        (RE[2] - wRS[2])a = 0            (LINEQ.EQ.FALSE.)
C
C        (RE[2] - wRS[2])a = RE[1]_B      (LINEQ.EQ.TRUE.)
C
C     where
C
C        RE[2] = B^{+}E[2]B
C        RS[2] = B^{+}S[2]B
C        RE[1] = B^{+}E[1}
C    
C-----------------------------------------------------------------------------
C
C     Variables:
C        GPCI    - property gradient , configurational part
C        GPOE    - property gradient , orbital part(e-e rotations)
C        GPOP    - property gradient , orbital part(e-p rotations)
C        REDGP   - reduced property gradient
C        REDE    - reduced Hessian matrix E[2]
C        REDS    - reduced metric         S[2]
C        IBTYP   - trial vectors will be of three kinds:
C                  C-type: contains only configurational parameters
C                  E-type: contains only e-e rotations
C                  P-type: contains only e-p rotations
C        IBCVC   - pointer from C-type trial vectors to IBTYP
C        IBEVC   - pointer from E-type trial vectors to IBTYP
C        IBPVC   - pointer from P-type trial vectors to IBTYP
C        EVALR   - eigenvalues of reduced system
C        EVECR   - eigenvectors of reduced system
C        
C     Files:
C
C        LUBCI   - direct access file of C-type trial vectors
C        LUBOE   - direct access file of E-type trial vectors
C        LUBOP   - direct access file of B-tpye trial vectors
C
C        LUSCI   - direct access file of C-part of sigma vectors
C        LUSOE   - direct access file of E-part of sigma vectors
C        LUSOP   - direct access file of P-part of sigma vectors
C  
C        LURSP   - unformatted file with restart information
C-----------------------------------------------------------------------------
C
C     In the frequency-dependent case, the eigenvectors come in pairs
C
C     b_1^T = [Z  Y^*]        b_2^T = [Y  Z^*]
C
C     The paired structure is exploited in that only the upper half
C     (Z and Y respectively) of the pair of trial vectors is stored.
C
C-----------------------------------------------------------------------------
C
C NEXSTV: number of start vectors
C NEXSIM: desired number of simultaneous roots in a microiteration
C In each microiteration the reduced matrices are extended by NTSIM
C trial vectors.
C    NTSIM = NCSIM + NESIM + NPSIM
C            NCSIM - number of configurational trial vectors
C            NESIM - number of orbital trial vectors (e-e rotations)
C            NPSIM - number of orbital trial vectors (e-p rotations)
C
C IBTYP is an array running over ALL trial vectors:
C   IBTYP(1,*) - indicates type of this trial vector
C   IBTYP(2,*) - indicates Hermiticity of trial vector
C                 1 - [Z+Y  Y*+Z*] - anti-Hermitian combination
C                 2 - [Z-Y  Y*-Z*] - Hermitian combination
C                    
C  ( TRIAL VECTORS ARE ONLY ADDED FOR THE NON CONVERGED VECTORS AND
C    LINEAR DEPENDENCE IS REMOVED)
C
C MAXRIT: MAXIMUM NUMBER OF MICROITERATIONS
C
C     Written by Hans Joergen Aa. Jensen and T.Saue 1996
C     Last revision: July 30 1996 - tsaue
C  
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER(D1 = 1.0D0)
C
      CHARACTER CPUTID*12,SECTID*12,TYP(6)*2
      DIMENSION GPCI(*), GPOE(*),GPOP(*),ISIM(6),
     &          IBTYP(2,*), IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(NREDM,NFREQ),
     &          RNORM(7,NFREQ),NPOS(*),WORK(*)         
C
C HJMAERKE: prepared for OPTORB, ORBSPC
      LOGICAL OPTORB,ORBSPC
      PARAMETER (ORBSPC=.FALSE.)
      DATA OPTORB/.FALSE./
C HJMAERKE: end
      DATA TYP/'E+','P+','C+','E-','P-','C-'/
C
#include <dcbibn.h>
C
C Used from common blocks:
C   XCBINF : MAXRM
C   XCBORB : N2ORBX
C   XCBDBG : ABCHK
C   XCBIO  : LUXXX*
C   XCBPRI : IPRINT, LUPRI, ?
C   WRKXRS : JSYMOP,LZYVAR,LZXOP*,?
C
C
#include <dcbxrs.h>
#include <pgroup.h>
#include <dgroup.h>
#include <dcborb.h>
C
      LOGICAL BINMEM
      CHARACTER LAB123(3)*8,RESTAR*8
C
      DATA LAB123/'********','********','XRSPONPP'/
      DATA RESTAR/'RESTART '/
C
      CALL QENTER('XRSCT1')
      KFRSAV = KFREE
C
C           WRITE (LUPRI,'(///A/A,A3,2X,A2/)')
C     &   ' <<< EXCITATION ENERGIES'//
C     &   ' AND TRANSITION MOMENT CALCULATION (PP) >>>',
C     &   ' Symmetry: ',REP(JSYMOP-1),TTYP(JTIMOP)
C
C     Determine maximum number of simultaneous trial vectors
C     and maximum number of vectors per batch in XRSNEX
C     ======================================================
C
CTROND Remember to activate XRSMCK !!!
CTROND      CALL XRSMCK(MAXSIM,MAXNEX,LFREE)
      MAXVEC = MAX(NEXSTV,NEXSIM)
C     ... balance load; i.e. set MAXSIM to optimal value .le. MAXVEC
CTROND      I = (MAXVEC-1)/MAXSIM + 1
CTROND      MAXSIM = (MAXVEC-1)/I + 1
CTROND      I = (MAXVEC-1)/MAXNEX + 1
CTROND      MAXNEX = (MAXVEC-1)/I + 1
CTROND ...these three lines should then be deleted....!
      MAXSIM = MAXVEC
      MAXNEX = MAXVEC
      NWRSP  = LFREE
CTROND
      IF (IPRXRS.GE.10) THEN
         WRITE(LUPRI,'(2(/A,I8))')
     *     'Maximum number of simultaneous trial vectors :',MAXSIM,
     *     'Maximum number of vectors per batch in XRSNEX:',MAXNEX
      END IF
C
C  Move outside this routine ???? HJMAERKE 941130
C  Test option ABCHK:
C    construct E2 and S2 explicitly using unit trial vectors
C  Test option ORBSPC:
C    solve orbital part of linear equations
C
C      IF ( ABCHK ) THEN
C         CALL XRSES2(GP,REDGP,REDE,REDS,EVALR,EVECR,
C     *               WORK,KFREE,LFREE)
C         WRITE(LUPRI,'(/A)')' .ABCHK E2 AND S2 CONSTRUCTED : STOP'
C         STOP 'XRSCT1: .ABCHK FINISHED'
C      END IF
C      IF (LINEQ .AND. ORBSPC) THEN
C#if defined (VAR_RESPONS)
C         DO 313 II = 1,LZXOPT
C            GP(II) = GP(NZCONFQ+II)
C 313     CONTINUE
CC        CALL DCOPY(LZXOPT,GP(1+NZCONFQ),1,GP,1)
C         DO 314 II = 1,LZXOPT
C            GP(LZXOPT+II) = GP(LZVAR+NZCONFQ+II)
C 314     CONTINUE
CC        CALL DCOPY(LZXOPT,GP(1+LZVAR+NZCONFQ),1,
CC    *              GP(1+LZXOPT),1)
C         CALL OROPTI(EVALR,IBTYP,GP,
C     *               WORK,kfree,lfree)
CC        CALL OROPTI(EVALR,IBTYP,A1,
CC    *               WORK,LWORK)
C#else
C         WRITE(LUPRI,*) 'XRSCT1 ERROR: .ORBSPC not implemented yet'
C#endif
C         WRITE(LUPRI,'(/A)')
C     * ' STOP: CHECK OF ORBITAL BLOCK OF E(2) AND S(2) FINISHED'
C         STOP 'XRSCT1: END OF .ORBSPC CHECK'
C      END IF
C
C     ***************************
C     *** I n i t i a l i z e ***
C     ***************************
C
      KFRITR = KFREE
C
      CALL XRSINI(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &            GPCI,GPOE,GPOP,REDGP,REDS,REDE,
     &            RCNV,EVALR,EVECR,RNORM,
     &            WORK,KFREE,LFREE)
      IF(NTSIM.LE.0) THEN
        KCONV = -3
        GOTO 40
      ENDIF
C
C     Go to equation solver if restart and no new
C     trial vectors
C     ===========================================
C
      IF (NZRED.GT.0.AND.NTSIM.LE.0) THEN
         IF (IPRXRS.GE.1) THEN
CTROND Restart warning from XRSST.....
           WRITE(LUPRI,'(/A,I5,A)') 
     &        ' ** XRSCTL Microiteration number',ITMIC,' **'
           CALL PRSYMB(LUPRI,'-',40,0)
         ENDIF
         TIMITR = SECOND()
         GO TO 30
      END IF
C
C     Set first coefficient for each frequency equal to one
C     for getting non-zero respons function with zero
C     microiterations
C
      IF(MAXITR.EQ.0) THEN
        DO I = 1,NFREQ
          EVECR(1,I) = D1
        ENDDO
        NZRED  = NTSIM
        NCRED  = NCSIM
        NERED  = NESIM
        NPRED  = NPSIM
      ENDIF
C
C     **********************************************
C     ***                                        ***
C     *** M I C R O I T E R A T I O N    L O O P ***
C     ***  - add new set of trial vectors to     ***
C     ***    reduced system                      ***
C     ***                                        ***
C     **********************************************
 10   CONTINUE
        IF(ITMIC.EQ.MAXITR) GOTO 40
        TIMITR = SECOND()
        CALL FLSHFO(LUPRI)
        ITMIC = ITMIC + 1
        IF (IPRXRS.GE.1) WRITE(LUPRI,'(/A,I5,A)') 
     &    ' ** XRSCTL microiteration number',ITMIC,' **'
        WRITE(LUPRI,'(A,A9)') '* Integral contributions: ',INTTYP
        CALL PRSYMB(LUPRI,'-',40,0)
C
CTROND       Check BINMEM !!!!!
C
        IF (IPRXRS.GE.5) WRITE(LUPRI,'(A,L1)')
     &      'BINMEM this microiteration: ',BINMEM
C
C       Fetch trial vectors
C       ===================       
C       If not in memory(BINMEM=F), read from file LUBVEC
C
        IF (.NOT.BINMEM) THEN
          NBTSIM = NCSIM + NESIM + NPSIM
          NBCSIM = NCSIM
          NBESIM = NESIM
          NBPSIM = NPSIM
        ENDIF
        JBCOFF = 0
        JBEOFF = 0
        JBPOFF = 0
        ITLOAD = 0
 20     CONTINUE
        ITLOAD = ITLOAD + 1
        IF (IPRXRS.GE.2 .AND. ITLOAD.GT.1) 
     &    WRITE(LUPRI,'(/A,I3)') '* Load no.',ITLOAD
C       Trial vectors not in memory; read from file
        IF(.NOT.BINMEM) THEN
          NWRST  = NWRSP
          CALL REABLD(NWRST,NBCSIM,NBESIM,NBPSIM,
     &                JBCOFF,JBEOFF,JBPOFF,
     &                WORK,KFREE,LFREE)
          JBCOFF = JBCOFF + NCSIM
          JBEOFF = JBEOFF + NESIM
          JBPOFF = JBPOFF + NPSIM
        ENDIF
        NOSIM = NESIM + NPSIM
        NTSIM = NCSIM + NOSIM
C
C       Print section
C       =============
C
        IF(IPRXRS.GE.10) THEN
          WRITE(LUPRI,'(A,I5)')
     &      'Load no.      :',ITLOAD,
     &      'C-type vectors:',NCSIM,
     &      'E-type vectors:',NESIM,
     &      'P-type vectors:',NPSIM
        ENDIF
C
        IF(IPRXRS.GE.5) THEN
          WRITE(LUPRI,'(/A)')' Before XRSSVC'
          IF(NCSIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NCSIM,' Config. trial vectors'
            CALL PRBVEC(LUPRI,WORK(KBCI),NCSIM,NZCONF)
          ENDIF
          IF(NESIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NESIM,' Electron trial vectors'
            CALL PRBVEC(LUPRI,WORK(KBOE),NESIM,NZXOPE)
          ENDIF
          IF(NPSIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NPSIM,' Positron trial vectors'
            CALL PRBVEC(LUPRI,WORK(KBOP),NPSIM,NZXOPP)
          ENDIF
        END IF
        CALL FLSHFO(LUPRI)
C
C       Calculate sigma-vectors (E[2]B)
C         Configurational part in WORK(KSCI)
C         e-e part in WORK(KSOE)
C         e-p part in WORK(KSOP)
C       ===============================================
C
        CALL MEMGET('REAL',KSCI,(NTSIM*NZCONFQ),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSOE,(NTSIM*NZXOPEQ),WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KSOP,(NTSIM*NZXOPPQ),WORK,KFREE,LFREE)
        CALL XRSSVC(WORK(KBCI),WORK(KBOE),WORK(KBOP),
     &              WORK(KSCI),WORK(KSOE),WORK(KSOP),
     &              IBTYP,IBCVC,IBEVC,IBPVC,NPOS,
     &              WORK,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
C
C       Write sigma-vectors to file 
C       ============================
C
        ICOUNT = NZRED
        CALL WRSVEC(ICOUNT,WORK(KSCI),WORK(KSOE),WORK(KSOP))
C
C       Update restart info file
C       ========================
C
        REWIND (LURSP)
        KTOT = NZRED + LOFFTY
        IF (IPRXRS.GE.7) THEN
          WRITE(LUPRI,'(/A,I5)')' Trial vectors for restart:',KTOT
        ENDIF
        WRITE(LURSP)LAB123,RESTAR
        WRITE(LURSP)JSYMOP,JTIMOP,KTOT,(IBTYP(1,I),IBTYP(2,I),I=1,KTOT)
        REWIND (LURSP)
C
C       Update matrices and vectors of reduced RSP equation
C       ===================================================
C
C
C       Check dimensions
C       ================
C
        IF (NZRED.GT.NREDM) THEN
           WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *       ' >>> ERROR IN XRSCTL >>>',
     *       ' Dimension of reduced space is  ',NZRED,
     *       ' which exceeds allowed dimension',NREDM
         CALL QUIT('XRSCTL: Too large dimension of reduced space')
        END IF
C
C       Update dimension of reduced equation
C       ====================================
C
        IZRED  = NZRED
        NZRED  = NZRED + NTSIM
        ICRED  = NCRED
        NCRED  = NCRED + NCSIM
        IERED  = NERED
        NERED  = NERED + NESIM
        IPRED  = NPRED
        NPRED  = NPRED + NPSIM
        CALL FLSHFO(LUPRI)
        CALL XRSRED (REDGP,REDE,REDS,GPCI,GPOE,GPOP,
     &               IBCVC,IBEVC,IBPVC,IBTYP,
     &               WORK(KBCI),WORK(KBOE),WORK(KBOP),
     &               WORK(KSCI),WORK(KSOE),WORK(KSOP),
     &               WORK,KFREE,LFREE)
        CALL MEMREL('XRSCT1.XRSRED',WORK,1,KFRITR,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
C
C       Perform a new load if some trial vectors still remains
C       ======================================================
        IF (.NOT.BINMEM) THEN
          NBCSIM = NBCSIM - NCSIM
          NBESIM = NBESIM - NESIM
          NBPSIM = NBPSIM - NPSIM
          NBTSIM = NBTSIM - NTSIM
          IF(NBTSIM.GT.0) GOTO 20
        ENDIF
C
 30     CONTINUE
C
C       Info on current trial vectors
C       -----------------------------
C
        CALL BVINFO(ISIM,IBTYP(1,IZRED+1),NTSIM)
C
C       Write info to PAMCYC
C       ====================
C
        WRITE(LUCYCL,'(I3,2X,A9,6(2X,A2,A1,I3))')
     &     ITMIC,INTTYP,(TYP(I),':',ISIM(I),I=1,6)
        WRITE(LUPRI,'(A,6(A2,A1,I3))') 
     &     '* New trial vectors: ',(TYP(I),':',ISIM(I),I=1,6)
C
C       Solve reduced RSP equation
C       ==========================
C
        CALL XRSSOL (REDGP,REDE,REDS,EVALR,EVECR,RCNV,
     &               IBTYP,WORK,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
C
C        CALCULATE TRANSITION MOMENT
C
C        CALL XRSMOM(IBTYP,EVALR,EVECR,BVECS,SVECS,PRVEC,NSIM)
C
C SET UP REDUCED MATRIX EXPLICLY AND CHECK BLOCK STRUCTURE IN E2 AND S2
C
C         IF (ABSYM) THEN
C            CALL XE2SYM(IBTYP,WORK,KFREE,LFREE)
C         END IF
C
C        Generate NEXSIM linear independent trial vectors for
C        next microiteration using solutions of reduced RSP equation
C        ===========================================================
C
         CALL XRSNEX(IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &               EVALR,EVECR,GPCI,GPOE,GPOP,
     &               RNORM,WORK,KFREE,LFREE)
C
         TIMITR = SECOND() - TIMITR
         CPUTID = SECTID(TIMITR)
C
C        Write info to PAMCYC
C        ====================
C
         WRITE(LUCYCL,'(6E10.2)') ((RNORM(I,J),I=1,6),J=1,NFREQ)
         WRITE(LUCYCL,'(A12)') CPUTID
         IF(IPRXRS.GE.1) THEN
           WRITE(LUPRI,'(/A,I5,A,A12/)')
     &       '* END OF MICROITERATION NO.',ITMIC,' IN ',CPUTID
         ENDIF
         CALL FLSHFO(LUPRI)
CTROND         OPTORB = OPTSAV
C
C        Check for state of the calculation
C        ==================================
C
         IF(KCONV.GT.0.AND.ITMIC.LT.MAXITR) THEN
           INTBUF = INTFLG
           GOTO 10
         ENDIF
C
C****************************************************
C****                                            ****
C**** E N D  O F   M I C R O I T E R A T I O N S ****
C****                                            ****
C****************************************************
C
 40   CONTINUE
      CALL XRSOUT(RCNV,RNORM,EVALR,IBTYP)
      IF(NZCONF.GT.0) THEN
        CLOSE(LUBCI,STATUS = 'KEEP')
        CLOSE(LUSCI,STATUS = 'DELETE')
      ENDIF
      IF(NZXOPE.GT.0) THEN
        CLOSE(LUBOE,STATUS = 'KEEP')
        CLOSE(LUSOE,STATUS = 'DELETE')
      ENDIF
      IF(NZXOPP.GT.0) THEN
        CLOSE(LUBOP,STATUS = 'KEEP')
        CLOSE(LUSOP,STATUS = 'DELETE')
      ENDIF
      CLOSE(LUCYCL,STATUS='DELETE')
      CALL MEMREL('XRSCT1',WORK,1,KFRSAV,KFREE,LFREE)
      CALL FLSHFO(LUPRI)
      CALL QEXIT('XRSCT1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSST */
      SUBROUTINE XRSST(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &                 IEVEC,IPVEC,ICVEC,GPCI,GPOE,GPOP,EVALR,
     &                 WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate start vectors in RSP
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision : July 30 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
      LOGICAL   BINMEM
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),
     &          GPCI(*),GPOE(*),GPOP(*),IEVEC(*),IPVEC(*),ICVEC(*),
     &          EVALR(*),WORK(*)
C
#include <dcbibn.h>
C
#include <dcbxrs.h>
#include <dgroup.h>
C
      CALL QENTER('XRSST')
      KFRSAV = KFREE
C
      IF(LINEQ) THEN
        IF(DIAGHE.OR.RSREST) THEN
          NVEC = NFREQ
        ELSE
          NVEC = 1
        ENDIF
        IF(STATIC) THEN
          NESIM = NVEC
          NPSIM = NVEC
          NCSIM = NVEC
        ELSE
          NESIM = NVEC*2
          NPSIM = NVEC*2
          NCSIM = NVEC*2
        ENDIF
        NTSIM = NESIM + NPSIM + NCSIM
        NBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
        CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
        KFIRST = KBTOT
C
        KBOE   = KFIRST
        KFIRST = KFIRST + NZXOPEQ*NESIM
C
        KBOP   = KFIRST
        KFIRST = KFIRST + NZXOPPQ*NPSIM
C
        KBCI   = KFIRST
        KFIRST = KFIRST + NZCONFQ*NCSIM
C
        KFRBUF = KFREE
        CALL MEMGET('INTE',KEVEC,NESIM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KPVEC,NPSIM,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KCVEC,NCSIM,WORK,KFREE,LFREE)
      ELSE
        CALL QUIT('XRSST: No code for LINEQ = F')
      ENDIF
C
C     If requested, allocate space for diagonal elements of Hessian (approximate)
C
      KEDIAH = KFREE
      KPDIAH = KFREE
      KCDIAH = KFREE
      IF(DIAGHE) THEN
        IF(NZXOPE.GT.0) 
     &    CALL MEMGET('REAL',KEDIAH,NZXOPE,WORK,KFREE,LFREE)
        IF(NZXOPP.GT.0) 
     &    CALL MEMGET('REAL',KPDIAH,NZXOPP,WORK,KFREE,LFREE)
        IF(NZCONF.GT.0) 
     &    CALL MEMGET('REAL',KCDIAH,NZCONF,WORK,KFREE,LFREE)
      ENDIF
      CALL  XRSST1(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &                 WORK(KEVEC),WORK(KPVEC),WORK(KCVEC),
     &                 GPCI,GPOE,GPOP,EVALR,
     &                 WORK(KEDIAH),WORK(KPDIAH),WORK(KCDIAH),
     &                 WORK,KFREE,LFREE)
      NOSIM = NESIM + NPSIM
      NTSIM = NCSIM + NOSIM
C
C     Compress trial vectors in memory
C     ================================
C
      MBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
      IF(MBTOT.NE.NBTOT) THEN
C     ...compress vectors in memory
        KFIRST = KBTOT
C
        JBOE   = KBOE
        KBOE   = KFIRST
        NEDIM  = NZXOPEQ*NESIM
        KFIRST = KFIRST + NEDIM
        IF(KBOE.NE.JBOE.AND.NEDIM.GT.0) 
     &   CALL DCOPY(NEDIM,WORK(JBOE),1,WORK(KBOE),1)
C
        JBOP   = KBOP
        KBOP   = KFIRST
        NPDIM  = NZXOPPQ*NPSIM
        KFIRST = KFIRST + NPDIM
        IF(KBOP.NE.JBOP.AND.NPDIM.GT.0) 
     &    CALL DCOPY(NPDIM,WORK(JBOP),1,WORK(KBOP),1)
C
        JBCI   = KBCI
        KBCI   = KFIRST
        NCDIM  = NZCONFQ*NCSIM
        KFIRST = KFIRST + NCDIM
        IF(KBCI.NE.JBCI.AND.NCDIM.GT.0) 
     &    CALL DCOPY(NCDIM,WORK(JBCI),1,WORK(KBCI),1)
C
        CALL MEMREL('XRSINI.CMP',WORK,1,KFRSAV,KFREE,LFREE)
        CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMREL('XRSINI.BUF',WORK,1,KFRBUF,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT ('XRSST')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSST1 */
      SUBROUTINE XRSST1(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &                 IEVEC,IPVEC,ICVEC,GPCI,GPOE,GPOP,EVALR,
     &                 HEDIAG,HPDIAG,HCDIAG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate start vectors in RSP
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision : July 30 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
C
#include <dcbibn.h>
C
#include <dcbxrs.h>
#include <dgroup.h>
      LOGICAL   BINMEM
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),EVALR(*),
     &          GPCI(*),GPOE(*),GPOP(*),IEVEC(*),IPVEC(*),ICVEC(*),
     &          HEDIAG(*),HPDIAG(*),HCDIAG(*),WORK(*)
C
      CALL QENTER('XRSST1')
      KFRSAV = KFREE
C
C     If requested, get diagonal elements of Hessian (approximate)
C
      IF(DIAGHE) CALL E2DIAG(HEDIAG,HPDIAG,HCDIAG,
     &                       WORK,KFREE,LFREE)
C
C
C     L I N E A R   R E S P O N S E   F U N C T I O N S
C     =================================================
C
C
      IND = NZRED
C
C     Orbtial (e-e) type trial vector
C     -------------------------------
C
      IF(NZXOPE.GT.0) THEN
        CALL XRSST2('E',WORK(KBOE),GPOE,HEDIAG,EVALR,IBTYP,IBEVC,IEVEC,
     &              JBENDX,NESIM,IND,NERED,NZXOPE,WORK,KFREE,LFREE)
        IND = IND + NESIM
      ELSE
        NESIM = 0
      END IF
C
C     Orbital (e-p) type trial vector
C     -------------------------------
C
      IF(NZXOPP.GT.0) THEN
        CALL XRSST2('P',WORK(KBOP),GPOP,HPDIAG,EVALR,IBTYP,IBPVC,IPVEC,
     &              JBPNDX,NPSIM,IND,NPRED,NZXOPP,WORK,KFREE,LFREE)
        IND = IND + NPSIM
      ELSE
        NPSIM = 0
      END IF
C
C     Configurational type trial vector
C     ---------------------------------
C
      IF(NZCONF.GT.0) THEN
        CALL XRSST2('C',WORK(KBCI),GPCI,HCDIAG,EVALR,IBTYP,IBCVC,ICVEC,
     &              JBCNDX,NCSIM,IND,NCRED,NZCONF,WORK,KFREE,LFREE)
        IND = IND + NCSIM
      ELSE
        NCSIM = 0
      END IF
C
      BINMEM = .TRUE.
C     HJMAERKE: insert appropriate test for BINMEM
C     CTROND: also reset KFREE for BINMEM = F !!!
      CALL QEXIT ('XRSST1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSORT */
      SUBROUTINE XRSORT(BVECS,NBTYP,NBNEW,NBPREV,
     &                  IBVEC,IBTYP,IVECS,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Purpose:
C     Orthogonalize new b-vectors against all previous b-vectors
C     and among themselves, and renormalize.
C
C     b-vectors generally have the structure
C
C              |  B  |
C       b(h) = | hB* |
C
C     where h = +/- 1 and only the upper part B is explicitly stored.
C    
C     b-vectors of opposite h-parameter are automatically 
C     orthogonal, provided that the have identical time reversal symmetry
C     structure, so only b-vectors of the same h need to be orthogonbalized
C     against each others.
C
C     Input:
C       BVECS,  non-orthogonal b-vectors
C       NBNEW,  number of b-vectors in BVECS
C       NBPREV, number of previous b-vectors on file
C       IBTYP(1,*) = NBTYP   b-vector type (CI, e-obital, p-orbital)
C       IBTYP(2,*) = H       h-parameter
C       IBVEC  - pointer to IBTYP
C
C     Output:
C       BVECS,  orthogonal b-vectors (also written to file)
C       NBNEW,  number of b-vectors written to file
C
C     Scratch:
C       BVBUF(*) - buffer for old b-vectors
C
C     Written by H.J.Aa.Jensen and T.Saue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbibn.h>
#include <thrldp.h>
      DIMENSION BVECS(*),IBVEC(*),IBTYP(2,*),IVECS(*),WORK(*)
C
C
C Used from common blocks:
C  WRKXRS: LZCONF,LZXOPT,JSYMOP
C  XCBPRI: LUPRI
C
#include <dgroup.h>
#include <dcbxrs.h>
C
      CALL QENTER('XRSORT')
C
C     Set dimensions; 
C     Set b vector input/output unit;
C     Set linear dependency limit =
C         sqrt( number of elements * machine # of digits squared ) =
C         sqrt( number of elements ) * machine # of digits
C
      KFRSAV = KFREE
CTROND : Rememebr LOFFTY !!!!!
      IF    (NBTYP.EQ.JBCNDX) THEN
        NBDIM  = NZCONFQ
        LUBVEC = LUBCI
        THRLDV = NZCONFQ*THRLDP
      ELSEIF(NBTYP.EQ.JBENDX) THEN
        NBDIM  = NZXOPEQ
        LUBVEC = LUBOE
        THRLDV = NZXOPEQ*THRLDP
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBDIM  = NZXOPPQ
        LUBVEC = LUBOP
        THRLDV = NZXOPPQ*THRLDP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'XRSORT ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('XRSORT: Unknown NBTYP !')
      ENDIF
      THRLDV = SQRT(THRLDV)
C
C     Memory allocation
C
      CALL MEMGET('REAL',KBVBUF,NBDIM,WORK,KFREE,LFREE)      
      CALL MEMGET('REAL',KBVNRM,NBNEW,WORK,KFREE,LFREE)      
      CALL XRSOR1(BVECS,NBDIM,NBNEW,NBPREV,THRLDV,IBVEC,IBTYP,
     &            WORK(KBVBUF),WORK(KBVNRM),LUBVEC,IVECS,NBTYP)
      CALL MEMREL('XRSORT',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XRSORT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XRSOR1 */
      SUBROUTINE XRSOR1(BVECS,NBDIM,NBNEW,NBPREV,THRLDP,IBVEC,IBTYP,
     &                  BVBUF,BVNRM,LUBVEC,IVECS,NBTYP)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (T1MIN = 1.D-8, THRRND = 1.D-4)
C
#include <dcbxrs.h>
      DIMENSION BVECS(NBDIM,NBNEW),IBTYP(2,NREDM),IBVEC(NREDM),
     &          BVBUF(NBDIM),IVECS(NBNEW),BVNRM(NBNEW)
      CHARACTER BTYP(-1:2)*1,XTYP(-1:1)*3
      DATA BTYP/'C',' ','E','P'/,XTYP/'(-)','   ','(+)'/
C
C
C     Take norm of new trial vectors
C     ==============================
C
      DO IBN = 1,NBNEW
        BVNRM(IBN) = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
      ENDDO
C
      IROUND = 0
      ITURN  = 0
 10   CONTINUE
      IREC   = 0
      ITURN  = ITURN + 1
C
C     Orthogonalize new b-vectors against previous b-vectors
C     ======================================================
C
      DO IB = 1,NBPREV
C  
C       Read old vector from LUBVEC
C
        IREC = IREC + 1
        CALL READAC(LUBVEC,NBDIM,BVBUF,IREC)
        IH   = IBTYP(2,IBVEC(IB))
        IF(.NOT.TKNORM) THEN
          BNBN = D4*DNRM2(NBDIM,BVBUF,1)**2
        ELSE
          BNBN = D1
        ENDIF
C
        DO 100 JBN = 1,NBNEW
          JB = JBN + NBPREV
          JH = IBTYP(2,IBVEC(JB))
          IF(IVECS(JBN).EQ.0.OR.(IH.NE.JH)) GOTO 100
C
C         Form -<BVECS|BVBUF>/<BVBUF|BVBUF> :
C
          BNBP = -D4*DDOT(NBDIM,BVECS(1,JBN),1,BVBUF,1)/BNBN
C
C         Form BVECS = BVECS -<BVECS|BVBUF>/<BVBUF|BVBUF>BVBUF
C
          CALL DAXPY(NBDIM,BNBP,BVBUF,1,BVECS(1,JBN),1)              
 100    CONTINUE
      ENDDO
C
C     Orthogonalize new vectors against each others
C     =============================================
C
      DO 200 IBN = 1,NBNEW        
        IF(IVECS(IBN).EQ.0) GOTO 200
        IB = IBN + NBPREV
        IH = IBTYP(2,IBVEC(IB))
C
        DO 300 JBN = 1,(IBN-1)
          JB = JBN + NBPREV
          JH = IBTYP(2,IBVEC(JB))
          IF(IVECS(JBN).EQ.0.OR.(IH.NE.JH)) GOTO 300
          IF(.NOT.TKNORM) THEN
            BNBN = D4*DNRM2(NBDIM,BVECS(1,JBN),1)**2
          ELSE
            BNBN = D1
          ENDIF
C
C         Form -<BVECS1|BVECS2>/<BVECS2|BVECS2> :
C
          BNBP = -D4*DDOT(NBDIM,BVECS(1,JBN),1,BVECS(1,IBN),1)/BNBN
C
C         Form BVECS1 = BVECS1 -<BVECS1|BVECS2>/<BVECS2|BVECS2>BVECS2
C
          CALL DAXPY(NBDIM,BNBP,BVECS(1,JBN),1,BVECS(1,IBN),1)              
 300    CONTINUE
C
C       Take norm of vector number IBN; 
C
        BNRM = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
C
C       Check size of norm
C       ==================
C
        IF    (BNRM.LT.THRLDP*BVNRM(IBN)) THEN
          WRITE(LUPRI,'(A,A1,A3,A,I4,A/8X,A,1P,D12.5)')
     &      'XRSORT: new ',BTYP(NBTYP),XTYP(IH),
     &      '-vector no. ',IBN,' is removed ',
     &      'because of linear dependence: ',BNRM
          IVECS(IBN) = 0
        ELSEIF(BNRM.LT.THRRND*BVNRM(IBN) .AND. ITURN.EQ.1) THEN
          IROUND = IROUND + 1
        ENDIF        
        IF(IPRXRS.GE.2) THEN
          WRITE (LUPRI,'(A,A1,A3,A,I6,5X,A,1P,D12.5)')
     &      'XRSORT: new ',BTYP(NBTYP),XTYP(IH),
     &      '-vector no.',IBN,
     &      'Norm : ',BNRM
        ENDIF
        IF(TKNORM) THEN
          IF(BNRM.LT.T1MIN) THEN
            BNRM = D1 / BNRM
            CALL DSCAL(NBDIM,BNRM,BVECS(1,IBN),1)
            BNRM = D2*DNRM2(NBDIM,BVECS(1,IBN),1)
          END IF
          BNRM = D1 / BNRM
          CALL DSCAL(NBDIM,BNRM,BVECS(1,IBN),1)
        ENDIF
 200  CONTINUE
C
C     Check whether a new round of orthonormalization is
C     necessary
C     ==================================================
C
      IF(IROUND.GT.0) THEN
         IROUND = 0
         IF (ITURN.EQ.1) GO TO 10
         CALL QUIT(
     &     'XRSORT error: IROUND .gt. 0 and ITURN .gt. 1')
      END IF
C
C     Save new vectors on file and update IBCVEC and IBTYP
C     ====================================================
C     (is is assumed that all new BVECS are consecutively 
C     ordered on IBTYP)
C
      IREC   = NBPREV
      DO 400 IBNB = 1,NBNEW
        IF(IVECS(IBNB).EQ.0) GOTO 400
        IB   = IBNB + NBPREV
        INDB = IBVEC(IB)
        IREC = IREC + 1
        CALL WRTDAC(LUBVEC,NBDIM,BVECS(1,IBNB),IREC)        
        IF(IB.NE.IREC) THEN
          IND          = IBVEC(IREC)
          IBTYP(1,IND) = IBTYP(1,INDB)
          IBTYP(2,IND) = IBTYP(2,INDB)
          IBN          = IREC - NBPREV
          CALL DCOPY(NBDIM,BVECS(1,IBNB),1,BVECS(1,IBN),1)
        ENDIF
 400  CONTINUE
      NBNEW = IREC - NBPREV
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Reabld */
      SUBROUTINE REABLD(NWRST,NBCSIM,NBESIM,NBPSIM,
     &                  JBCOFF,JBEOFF,JBPOFF,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Read a load of B-vectors, limited by available memory NWRST
C
C     Written by T.Saue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
      DIMENSION WORK(*)
C
      CALL QENTER('REABLD')
      NWREQ  = 0
      KFIRST = KBTOT
C
C     Electronic type trial vectors
C     -----------------------------
C
      KBOE = KFIRST
      IF(NBESIM.GT.0) THEN
        NESIM = MIN((NWRST/NZXOPEQ),NBESIM)
        NBUF   = NESIM*NZXOPEQ
        IF(NESIM.GT.0) THEN
          CALL MEMGET('REAL',KBOE,NBUF,WORK,KFREE,LFREE)
          IREC = NERED + JBEOFF
          JOFF = KBOE
          DO I = 1,NESIM
            IREC = IREC + 1
            CALL READAC(LUBOE,NZXOPEQ,WORK(JOFF),IREC)
            JOFF = JOFF + NZXOPEQ
          ENDDO
        ENDIF
        NWREQ  = NWREQ + NBUF
        NWRST  = NWRST - NBUF
        KFIRST = KFIRST + NZXOPEQ*NESIM
      ENDIF
C
C     Positronic type trial vectors
C     -----------------------------
C
      KBOP = KFIRST
      IF(NBPSIM.GT.0) THEN
        NPSIM = MIN((NWRST/NZXOPPQ),NBPSIM)
        NBUF   = NESIM*NZXOPPQ
        IF(NPSIM.GT.0) THEN
          CALL MEMGET('REAL',KBOP,NBUF,WORK,KFREE,LFREE)
          IREC = NPRED + JBPOFF
          JOFF = KBOP
          DO I = 1,NPSIM
            IREC = IREC + 1
            CALL READAC(LUBOP,NZXOPPQ,WORK(JOFF),IREC)
            JOFF = JOFF + NZXOPPQ
          ENDDO
        ENDIF
        NWREQ  = NWREQ + NBUF
        NWRST  = NWRST - NBUF
        KFIRST = KFIRST + NZXOPPQ*NPSIM
      ENDIF
C
C     Configurational type trial vectors
C     ----------------------------------
C
      KBCI = KFIRST
      IF(NBCSIM.GT.0) THEN
        NCSIM = MIN((NWRST/NZCONFQ),NBCSIM)
        NBUF  = NCSIM*NZCONFQ
        IF(NCSIM.GT.0) THEN
          CALL MEMGET('REAL',KBCI,NBUF,WORK,KFREE,LFREE)
          IREC = NCRED + JBCOFF
          JOFF = KBCI
          DO I = 1,NCSIM
            IREC = IREC + 1
            CALL READAC(LUBCI,NZCONFQ,WORK(JOFF),IREC)
            JOFF = JOFF + NZCONFQ
          ENDDO
        ENDIF
        NWREQ  = NWREQ + NBUF
        NWRST  = NWRST - NBUF
        KFIRST = KFIRST + NZCONFQ*NCSIM
      ENDIF
C
      CALL QEXIT('REABLD')      
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck wrsvec */
      SUBROUTINE WRSVEC(ICOUNT,SCVEC,SEVEC,SPVEC)
C*****************************************************************************
C
C     Write sigma-vectors to file
C
C     Written by T.Saue Sep 18 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
      DIMENSION SCVEC(*),SEVEC(*),SPVEC(*)
C
C       Configurational part
C       --------------------
C
      IF(NZCONF.GT.0) THEN
        IREC = ICOUNT
        JOFF = 1
        DO ISIM = 1,NTSIM
          IREC = IREC + 1
          CALL WRTDAC(LUSCI,NZCONFQ,SCVEC(JOFF),IREC)            
          JOFF = JOFF + NZCONFQ
        ENDDO
      ENDIF        
C
C     Electronic part
C     ---------------
C
      IF(NZXOPE.GT.0) THEN
        IREC = ICOUNT
        JOFF = 1
        DO ISIM = 1,NTSIM
          IREC = IREC + 1
          CALL WRTDAC(LUSOE,NZXOPEQ,SEVEC(JOFF),IREC)            
          JOFF = JOFF + NZXOPEQ
        ENDDO
      ENDIF        
C
C     Positronic part
C     ---------------
C
      IF(NZXOPP.GT.0) THEN
        IREC = ICOUNT
        JOFF = 1
        DO ISIM = 1,NTSIM
          IREC = IREC + 1
          CALL WRTDAC(LUSOP,NZXOPPQ,SPVEC(JOFF),IREC)            
          JOFF = JOFF + NZXOPPQ
        ENDDO
      ENDIF        
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsred */
      SUBROUTINE XRSRED(REDGP,REDE,REDS,GPCI,GPOE,GPOP,
     &                  IBCVC,IBEVC,IBPVC,IBTYP,BCVEC,BEVEC,BPVEC,
     &                  SCVEC,SEVEC,SPVEC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     PURPOSE:
C     Extend the reduced matrices with the NTSIM new trial vectors.
C
C     The reduced matrices are
C       Reduced Hessian:  RE[2] = B^D*E[2]*B     (D means dagger)
C       Reduced metric :  RS[2] = B^D*s[2]*B
C       E[2]B is already computed and given by the sigma-vector.
C     In the linear equation (LINEQ.EQ.TRUE) there is in addition the
C     reduced property gradient
C       Reduced gradient: RE[1] = B^D*E[1]
C
C     In the frequency-dependent case, trial b-vectors and sigma vectors
C     come in pairs
C
C        | Z  |    | Y  |
C        | Y* |  , | Z* |
C
C     Only the parts Z and Y are stored. Likewise the gradient has the
C     structure
C
C        | G* |
C        | G  |
C
C     Here only G* is stored.
C
C     In the frequency-independent case (STATIC) Z = Y, so the dimension of
C     the problem is reduced by a factor 2.
C
C   Input:
C     REDGP - reduced property gradient
C     REDE  - old reduced E[2]-matrix (dimension: NZRED-NTSIM)
C     REDS  - old reduced S[2]-matrix (dimension: NZRED-NTSIM)
C     BCVEC - new trial vectors corresponding to configurational parameters
C     BEVEC - new trial vectors corresponding to electron-electron rotations
C     BPVEC - new trial vectors corresponding to electron-positron rotations
C     SCVEC - configurational part of new sigma-vectors
C     SEVEC - electron-electron part of new sigma vectors
C     SPVEC - electron-postrion part of new sigma vectors
C     S times the N b-vectors are created in XRSSLI and stored in sigma-vectors
C     after extending the reduced Hessian
C
C   Output:
C     REDE, the new, extended reduced E[2]-matrix (dimension: NZRED)
C     REDS, rhe new, extended reduced S[2]-matrix (dimension: NZRED)
C
C   Scratch:
C
C     WORK, real scratch array dimension KZYVAR
C
C     Written by H.J.Aa.Jensen and T.Saue
C     Last revision Sep 19 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <ibndxdef.h>
      PARAMETER(DM1 = -1.0D0)
C
#include <dcbxrs.h>
#include <dgroup.h>
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),REDS(NREDM,NREDM),
     &          GPCI(*),GPOE(*),GPOP(*),
     &          IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          WORK(*)
C
      CALL QENTER('XRSRED')
C
C     Exit routine if no new trial vectors
C     ====================================
C
      IF (NTSIM.LT.1) GO TO 9999
C
C     Memory check
C     ============
C
      KFRSAV = KFREE
CTROND: Skipping this for now
CTROND#if defined (VAR_HJMAERKE)
CTROND      LNEED  = NCREF
CTRONDC     ... needed in XRSSLI (+ more)
CTRONDCTSMAERKE: skal det ikke vaere LZYCONF ???
CTROND      LNEED = MAX(LNEED,LZYXOPE,LZYXOPP,LZCONF)
CTROND      IF (LNEED .GT. LFREE) CALL ERRWORK('XRSRED',-LNEED,LFREE)
CTROND#endif
CTROND: End of skip
C
C     **********************************************************************
C     ***** If linear equations (LINEQ), then extend reduced gradients *****
C     ***** (we make the negative gradient)                            *****
C     **********************************************************************
C
C     The gradient is a +-type vector and only gets contributions for
C     trial vector of (+)-type, that is for which IBTYP(2,*) equals 1
C
      IF(LINEQ)THEN
        CALL EXTRVC(IPX,REDGP,GPCI,GPOE,GPOP,BCVEC,BEVEC,BPVEC,
     &              IBTYP,IBCVC,IBEVC,IBPVC)
        CALL DSCAL(NTSIM,DM1,REDGP(IZRED+1),1)
        IF(IPRXRS.GE.4) THEN
          WRITE(LUPRI,'(A)') '* XRSRED: Minus reduced gradient:'
          CALL PRIRVC(REDGP,NZRED,IBTYP)
C          WRITE(LUPRI,'(3X,F16.8)') (REDGP(I),I=1,NZRED)
        ENDIF
      ENDIF
C
C     **************************************
C     ****  Extend reduced E[2]-matrix  ****
C     **************************************
C
C     We need only make the upper triangle explicitly
C     since the matrices are real symmetric
C
C     The E[2}-matrix conserves the IH-parameter,
C     thus it only couples trial vectors of the same kind
C       (IH1 = IH2)
C
C     Allocate work space for previous trial vectors
      NBVEC = MAX(NZCONFQ,NZXOPEQ,NZXOPPQ)
      CALL MEMGET('REAL',KB1,NBVEC,WORK,KFREE,LFREE)
      CALL EXTRDE(REDE,IBTYP,IBCVC,IBEVC,IBPVC,
     &            BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
CTROND      CALL EXTRMT(1,REDE,IBTYP,IBCVC,IBEVC,IBPVC,
CTROND     &            BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
      CALL MEMREL('EXTRE2',WORK,1,KB1,KFREE,LFREE)
C
C     Print section
C     =============
C
      IF(IPRXRS.GE.4) THEN
        WRITE(LUPRI,'(A,I5)') 
     &    'XRSRED: Reduced Hessian. Dimension:',NZRED
        CALL PRIRMT(REDE,NREDM,NZRED,IBTYP)
      ENDIF
C
C     **************************************
C     ****  Extend reduced S[2]-matrix  ****
C     **************************************
C
C     We need only make the upper triangle explicitly
C     since the matrices are real symmetric
C
C     The S[2}-matrix switches the IH-parameter,
C     thus it only couples trial vectors of the opposite kind
C       (IH1 = - IH2)
C
C     In the static case, we only have the (+)-combination,
C     thus there can be no contribution from the
C     generalized metric.
C     
C     Metric times trial vectors is stored in SVECS...
C
      IF(.NOT.STATIC) THEN
C
C       Generate tau-vectors: S[2]b
C       ---------------------------
C      
        CALL XRSTVC(BEVEC,IBEVC,NESIM,IERED,
     &              BPVEC,IBPVC,NPSIM,IPRED,
     &              BCVEC,IBCVC,NCSIM,ICRED,
     &              SEVEC,SPVEC,SCVEC,NTSIM,IZRED,
     &              WORK,KFREE,LFREE)
C
C       Extend reduced S-matrix
C       -----------------------
C
        NBVEC = MAX(NZCONFQ,NZXOPEQ,NZXOPPQ)
        CALL MEMGET('REAL',KB1,NBVEC,WORK,KFREE,LFREE)
        CALL EXTRMT(-1,REDS,IBTYP,IBCVC,IBEVC,IBPVC,
     &              BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
        CALL MEMREL('EXTRS2',WORK,1,KB1,KFREE,LFREE)
C
C       Print section
C       =============
C
        IF(IPRXRS.GE.4) THEN
          WRITE(LUPRI,'(A,I5)') 
     &      'XRSRED: Reduced metric. Dimension:',NZRED
C            CALL OUTPUT(REDS,1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
        CALL PRIRMT(REDS,NREDM,NZRED,IBTYP)
        ENDIF
      ENDIF
C
C *** End of subroutine XRSRED
C
 9999 CONTINUE
      CALL QEXIT('XRSRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck xrsini */
      SUBROUTINE XRSINI(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &            GPCI,GPOE,GPOP,REDGP,REDS,REDE,
     &            RCNV,EVALR,EVECR,RNORM,
     &            WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Initialize XRS-module
C
C     Written by T.Saue Oct 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <dcbibn.h>
C
#include <dcbfir.h>
#include <dcbxrs.h>
#include <dgroup.h>
      LOGICAL   BINMEM
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),
     &          GPCI(*),GPOE(*),GPOP(*), WORK(*),
     &          REDGP(*), REDS(*), REDE(*), EVALR(*), 
     &          EVECR(NREDM,NFREQ),RNORM(*),RCNV(*)
C
      CALL QENTER('XRSINI')
      KFRSAV = KFREE
C
C     File control
C     ============
C
      LURSP  = 30
      LUBCI  = 31
      LUBOE  = 32
      LUBOP  = 33
      LUSCI  = 34
      LUSOE  = 35
      LUSOP  = 36
      LUCYCL = 37
C
C     LURSP - file of restart info
C
      CALL OPNFIL(LURSP,'PAMRSP','UNKNOWN','PAMRSP')      
C
C     LUB* - files of b-vectors
C     LUS* - files of sigma-vectors
C
      IF(NZCONF.GT.0) THEN
        OPEN(LUBCI,FILE='PAMBCI',FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*NZCONFQ,STATUS='UNKNOWN')
        OPEN(LUSCI,FILE='PAMSCI',FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*NZCONFQ,STATUS='UNKNOWN')
      ENDIF
      IF(NZXOPE.GT.0) THEN
        OPEN(LUBOE,FILE='PAMBOE',FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
        OPEN(LUSOE,FILE='PAMSOE',FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*NZXOPEQ,STATUS='UNKNOWN')
      ENDIF
      IF(NZXOPP.GT.0) THEN
        OPEN(LUBOP,FILE='PAMBOP',FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
        OPEN(LUSOP,FILE='PAMSOP',FORM='UNFORMATTED',
     +      ACCESS='DIRECT',RECL=8*NZXOPPQ,STATUS='UNKNOWN')
      ENDIF
      OPEN(LUCYCL,FILE ='PAMCYC',STATUS='UNKNOWN',
     &            ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
C
C     Initialize
C     ==========
C
      IPX  =  JTIMOP
      IMX  = -JTIMOP
      NTYP = 2
      IF(STATIC) NTYP = 1
      ITMIC = 0
      CALL DCOPY(NFREQ,DUMMY,0,RCNV,1)
      CALL DZERO(EVECR,NREDM*NFREQ)
      CALL DZERO(REDE,N2REDM)
      CALL DZERO(REDS,N2REDM)
      IF (LINEQ) CALL DZERO(REDGP,NREDM)
      NRNORM = 7*NFREQ
      CALL DZERO(RNORM,NRNORM)
      FIRST1 = .TRUE.
      FIRST2 = .TRUE.
      FIRST3 = .TRUE.
C
C     Determine integral contribution
C     ===============================
C
      INTBUF = 0
      ITNEX  = ITMIC + 1
      CALL INTCON(INTFLG,INTBUF,INTDEF,
     &            DUMMY,CNVINT,ITNEX,ITRINT,INTTYP)
C
C     Set up initial trial vectors
C     ============================     
C
      CALL XRSST(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &                 IEVEC,IPVEC,ICVEC,GPCI,GPOE,GPOP,EVALR,
     &                 WORK,KFREE,LFREE)
      CALL QEXIT ('XRSINI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrmt */
      SUBROUTINE EXTRMT(IPAR,REDMAT,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,BVBUF)
C*****************************************************************************
C
C     Extend reduced matrices in response equations
C     IPAR =  1	Hessian-type :
C               Coupling of trial vectors of same type (IH = JH)
C     IPAR = -1 Metric type:
C               Coupling of trial vectors of opposite type (IH = -JH)
C
C     Written by T.Saue Oct 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4=4.0D0)
C
#include <dcbxrs.h>
C
      DIMENSION IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          BVBUF(*),REDMAT(NREDM,NREDM)
C
C     C-type trial vectors
C     ====================
      DO I = 1,NCSIM
CTROND If LOFFTY, remember to skip reference vector
CTROND    Nothing written yet.....
      ENDDO
C
C     E-type trial vectors
C     ====================
C
C     Read in old vectors and multiply with all new sigma-vectors
C     -----------------------------------------------------------
C
      DO IREC = 1,IERED
        CALL READAC(LUBOE,NZXOPEQ,BVBUF,IREC)
        I  = IBEVC(IREC)
        IH = IBTYP(2,I)
        DO JN = 1,NTSIM
          J  = JN + IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPEQ,BVBUF,1,SEVEC(1,JN),1)
            REDMAT(I,J) = REA  
            REDMAT(J,I) = REA
          ENDIF
        ENDDO          
      ENDDO
C
C     For new vectors only upper triangular loop is needed
C     ----------------------------------------------------
C
      DO IN = 1,NESIM
        I  = IBEVC(IN+IERED)
        IH = IBTYP(2,I)
        DO J = I,NZRED
          JN = J - IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPEQ,BEVEC(1,IN),1,SEVEC(1,JN),1)
            REDMAT(I,J) = REA  
            REDMAT(J,I) = REA
          ENDIF
        ENDDO
      ENDDO
C
C     P-type trial vectors
C     ====================
C
C
C     Read in old vectors and multiply with all new sigma-vectors
C     -----------------------------------------------------------
C
      DO IREC = 1,IPRED
        CALL READAC(LUBOP,NZXOPPQ,BVBUF,IREC)
        I  = IBPVC(IREC)
        IH = IBTYP(2,I)
        DO JN = 1,NTSIM
          J  = JN + IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPPQ,BVBUF,1,SPVEC(1,JN),1)
            REDMAT(I,J) = REA  
            REDMAT(J,I) = REA
          ENDIF
        ENDDO          
      ENDDO
C
C     For new vectors only upper triangular loop is needed
C     ----------------------------------------------------
C
      DO IN = 1,NPSIM
        I  = IBPVC(IN+IPRED)
        IH = IBTYP(2,I)
        DO J = I,NZRED
          JN = J - IZRED
          JH = IBTYP(2,J)
          IF(IH*JH.EQ.IPAR) THEN
            REA = D4*DDOT(NZXOPPQ,BPVEC(1,IN),1,SPVEC(1,JN),1)
            REDMAT(I,J) = REA  
            REDMAT(J,I) = REA
          ENDIF
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrde */
      SUBROUTINE EXTRDE(REDMAT,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,BVBUF)
C*****************************************************************************
C
C     Extend reduced Hessian in response equations
C
C     Written by T.Saue and H.J.Aa.Jensen Sep 29 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4=4.0D0)
C
#include <dcbxrs.h>
C
      DIMENSION IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          BVBUF(*),REDMAT(NREDM,NREDM)
C
      IPAR = 1
C
C     C-type trial vectors
C     ====================
      DO I = 1,NCSIM
CTROND If LOFFTY, remember to skip reference vector
CTROND    Nothing written yet.....
      ENDDO
C
C     E-type trial vectors
C     ====================
C
      IF(NESIM.GT.0) THEN
C
C     Read in all old sigma-vectors and multiply with new E-type trial-vectors
C     ------------------------------------------------------------------------
C
        DO I = 1,IZRED
          CALL READAC(LUSOE,NZXOPEQ,BVBUF,I)
          IH = IBTYP(2,I)
          DO JN = 1,NESIM
            J  = IBEVC(JN+IERED)
            JH = IBTYP(2,J)
            IF(IH*JH.EQ.IPAR) THEN
              REA = D4*DDOT(NZXOPEQ,BVBUF,1,BEVEC(1,JN),1)
              REDMAT(I,J) = REA  
              REDMAT(J,I) = REA
            ENDIF
          ENDDO          
        ENDDO
C
C     For new vectors only upper triangular loop is needed
C     ----------------------------------------------------
C
        DO IN = 1,NESIM
          I  = IBEVC(IN+IERED)
          IH = IBTYP(2,I)
          DO J = I,NZRED
            JN = J - IZRED
            JH = IBTYP(2,J)
            IF(IH*JH.EQ.IPAR) THEN
              REA = D4*DDOT(NZXOPEQ,BEVEC(1,IN),1,SEVEC(1,JN),1)
              REDMAT(I,J) = REA  
              REDMAT(J,I) = REA
            ENDIF
          ENDDO
        ENDDO
      ENDIF     
C
C     P-type trial vectors
C     ====================
C
      IF(NPSIM.GT.0) THEN
C
C     Read in old sigma-vectors and multiply with all new trial-vectors
C     -----------------------------------------------------------
C
        DO I = 1,IZRED
          CALL READAC(LUSOP,NZXOPPQ,BVBUF,I)
          IH = IBTYP(2,I)
          DO JN = 1,NPSIM
            J  = IBPVC(JN+IPRED)
            JH = IBTYP(2,J)
            IF(IH*JH.EQ.IPAR) THEN
              REA = D4*DDOT(NZXOPPQ,BVBUF,1,BPVEC(1,JN),1)
              REDMAT(I,J) = REA  
              REDMAT(J,I) = REA
            ENDIF
          ENDDO          
        ENDDO
C
C     For new vectors only upper triangular loop is needed
C     ----------------------------------------------------
C
        DO IN = 1,NPSIM
          I  = IBPVC(IN+IPRED)
          IH = IBTYP(2,I)
          DO J = I,NZRED
            JN = J - IZRED
            JH = IBTYP(2,J)
            IF(IH*JH.EQ.IPAR) THEN
              REA = D4*DDOT(NZXOPPQ,BPVEC(1,IN),1,SPVEC(1,JN),1)
              REDMAT(I,J) = REA  
              REDMAT(J,I) = REA
            ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrvc */
      SUBROUTINE EXTRVC(IPAR,REDVC,FCVEC,FEVEC,FPVEC,
     &                  BCVEC,BEVEC,BPVEC,IBTYP,IBCVC,IBEVC,IBPVC)
C*****************************************************************************
C
C     Extend reduced vectors in response equations
C     IPAR =  1	(+)-vector (gradient-type)
C     IPAR = -1 (-)-vector
C
C     Written by T.Saue Oct 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D4=4.0D0)
C
#include <dcbxrs.h>
#include <dgroup.h>
C
      DIMENSION IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          FCVEC(NZCONFQ,*),FEVEC(NZXOPEQ,*),FPVEC(NZXOPPQ,*),
     &          REDVC(NREDM)
C
C     C-type trial vectors
C     ====================
      DO I = 1,NCSIM
CTROND    Nothing written yet.....
      ENDDO
C
C     E-type trial vectors
C     ====================
C
      DO I = 1,NESIM
        IND        = IBEVC(I+IERED)
        IF(IBTYP(2,IND).EQ.IPAR) THEN
          REDVC(IND) = D4*DDOT(NZXOPEQ,BEVEC(1,I),1,FEVEC,1)
        ENDIF
      ENDDO
C
C     P-type trial vectors
C     ====================
C
      DO I = 1,NPSIM
        IND        = IBPVC(I+IPRED)
        IF(IBTYP(2,IND).EQ.IPAR) THEN
          REDVC(IND) =  D4*DDOT(NZXOPPQ,BPVEC(1,I),1,FPVEC,1)
        ENDIF
      ENDDO
C
C     Print section
C     =============
C
      IF(IPRXRS.GE.200) THEN
        IREP = JSYMOP-1
        IF(NZXOPE.GT.0) THEN
          WRITE(LUPRI,'(A)')'EXTRVC : Full vector (e-e) part:'
          CALL PRQMAT(FEVEC,NZXOPE,1,NZXOPE,1,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
        IF(NZXOPP.GT.0) THEN
          WRITE(LUPRI,'(A)')'EXTRVC : Full vector (e-p) part'
          CALL PRQMAT(FPVEC,NZXOPP,1,NZXOPP,1,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDIF
      IF(IPRXRS.GE.50) THEN
        WRITE(LUPRI,'(A,I3)')
     &    'EXTRVC: Reduced vector Dim.: NZRED = ',NZRED
        CALL OUTPUT(REDVC,1,NZRED,1,1,NREDM,1,1,LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssol */
      SUBROUTINE XRSSOL (REDGP,REDE,REDS,EVALR,EVECR,
     *                   RCNV,IBTYP,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     This routine solves the response equation in reduced subspace
C     of dimension NZRED.
C
C     LINEQ = T:       (REDE - wREDS) X = REDG  (linear system of equations)
C
C     LINEQ = F:       (REDE - wREDS) X = 0     (eigenvalue equation)
C
C     Input variables:
C       REDE(NREDM,NZRED) - reduced Hessian
C       REDS(NREDM,NZRED) - reduced metric
C       REDG (NFREQ)      - reduced property gradient (if LINEQ)
C       EVALR(NFREQ)      - frequencies
C
C     Output variables:
C       EVALR(NZRED)         - real part of eigenvalues(excitation energies) 
C                              (if. NOT.LINEQ)
C       EVECR(NREDM,NZRED,*) - solutions vectors in reduced space
C
C
C     Written by T.Saue and Hans Joergen Aa. Jensen. 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
C
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),REDS(NREDM,NREDM),
     &          EVALR(NFREQ),EVECR(NREDM,NFREQ),IBTYP(2,*),RCNV(*),
     &          WORK(*)
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0 )
      PARAMETER ( DRTEST = 1.0D-5 )
C
      CALL QENTER('XRSSOL')
C
C     Check dimension of reduced problem
C     ==================================
C
      IF (NZRED.GT.NREDM) THEN
         WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *   ' >>> ERROR IN XRSSOL >>>',
     *   ' Dimension of reduced space is  ',NZRED,
     *   ' which exceeds allowed dimension',NREDM
         CALL QUIT('XRSSOL: Too large dimension of reduced space')
      END IF
C
C     Print section
C     =============
C
      IF(IPRXRS.GE.50) THEN
        WRITE(LUPRI,'(A,I5)') 
     &    'XRSSOL: Reduced Hessian. Dimension:',NZRED
C          CALL OUTPUT(REDE,1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
        CALL PRIRMT(REDE,NREDM,NZRED,IBTYP)
        IF(.NOT.STATIC) THEN
          WRITE(LUPRI,'(A,I5)') 
     &      'XRSSOL: Reduced metric. Dimension:',NZRED
C            CALL OUTPUT(REDS,1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
        CALL PRIRMT(REDS,NREDM,NZRED,IBTYP)
        ENDIF
      ENDIF
C
      IF (LINEQ) THEN
C
C     Solve reduced linear response problem in subspace
C     =================================================
C
         DO IFREQ = 1,NFREQ
           AFREQ = EVALR(IFREQ)
           CALL LINRED(NZRED,NREDM,REDE,REDS,AFREQ,REDGP,
     &                 EVECR(1,IFREQ),IBTYP,IPRXRS,WORK,KFREE,LFREE)
         ENDDO
C
C        Print section
C        -------------
C
         IF(IPRXRS.GE.1) THEN
           WRITE(LUPRI,'(/A,I4)') '<<B;B>> at iteration: ',ITMIC
           CALL PRSYMB(LUPRI,'-',44,0)
           DO IFREQ = 1,NFREQ
           VLRF = DDOT(NZRED,REDGP,1,EVECR(1,IFREQ),1)
           IF(RCNV(IFREQ).GT.THCXRS) THEN
             WRITE(LUPRI,'(3X,A,F15.8,5X,F15.8)')
     &         'Freq.:',EVALR(IFREQ),VLRF
           ELSE
              WRITE(LUPRI,'(3X,A,F15.8,5X,F15.8,3X,A)')
     &         'Freq.:',EVALR(IFREQ),VLRF,'(converged)'
           ENDIF
           ENDDO
        ENDIF
      ELSE
C
C
C     Solve reduced general problem in subspace
C     =========================================
CTROND: This needs to be finished !
        CALL EIGRED(NZRED,NREDM,REDE,REDS,EVALR,EVECR,
     &          IERR,IPRXRS,WORK,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT('XRSSOL')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Linred */
      SUBROUTINE LINRED(N,LDA,AE,AS,W,G,X,IBTYP,IPRINT,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Solve the reduced linear reponse equation
C
C       [AE - W*AS]X = G
C
C     Written by T.Saue Oct 16 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0)
C
      DIMENSION AE(LDA,N),AS(LDA,N),G(N),X(N),
     &          DET(2),INERT(3),IBTYP(2,*),WORK(*)
C
      CALL QENTER('LINRED')
      KFRSAV = KFREE
C
C     Memory allocation
C
      N2A = LDA*N
      NNA = (N*(N+1))/2
      CALL MEMGET('REAL',KA2,N2A,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAN,NNA,WORK,KFREE,LFREE)
C
C     Form matrix [AE - W*AS] in row-packed lower triangle form
C
      CALL DCOPY(N2A,AE,1,WORK(KA2),1)
      IF (W.NE.D0) THEN
        CALL DAXPY(N2A,-W,AS,1,WORK(KA2),1)
      END IF
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A)') 'LINRED: Reduced E[2] - WS[2} matrix'
C        CALL OUTPUT(WORK(KA2),1,N,1,N,LDA,N,1,LUPRI)
        CALL PRIRMT(WORK(KA2),LDA,N,IBTYP)
      ENDIF
      CALL DAMATR(N,WORK(KA2),LDA,WORK(KAN))
C
C     Place G in buffer X
C
      CALL DCOPY(N,G,1,X,1)
C
C     Solve linear system
C
      CALL DSPSLI(N,1,WORK(KAN),X,WORK(KA2),INFO,DET,INERT)
C
C     Print section
C
      IPAR = INERT(2) + INERT(3)
      IF(IPRINT.GE.2.OR.IPAR.GT.0) THEN
        WRITE(LUPRI,'(A,F12.5,A)') 
     &    'LINRED: Solution vector at frequency:',W,' a.u.'
C        WRITE(LUPRI,'(5X,5F12.5)') (X(I),I=1,N)
        CALL PRIRVC(X,N,IBTYP)
        WRITE(LUPRI,'(/3X,A)') '* Reduced matrix:'
        WRITE(LUPRI,'(5X,A,I5)')
     &    '* Number of positive eigenvalues:',INERT(1),
     &    '* Number of zero     eigenvalues:',INERT(2),
     &    '* Number of negative eigenvalues:',INERT(3)
      ENDIF
C
      CALL MEMREL('LINRED',WORK,1,KFREE,KFREE,LFREE)
      CALL QEXIT('LINRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Eigred */
      SUBROUTINE EIGRED(N,LDA,AE,AS,W,X,IERR,IPRINT,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     This subroutine solves the reduced eigenvalue equation
C
C       AE*X = W*AS*X
C
C     Written by T.Saue Oct 16 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION AE(LDA,N),AS(LDA,N),W(N),X(LDA,N),
     &          WORK(*)
C
      CALL QENTER('EIGRED')
C
C     Memory allocation
C
      N2A = LDA*N
      CALL MEMGET('REAL',KAE  ,N2A,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KAS  ,N2A,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KALFR,N  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KALFI,N  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBETA,N  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEVEC,N2A,WORK,KFREE,LFREE)
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
     &         WORK(KBETA),MATZ,WORK(KEVEC),IERR)
C     IERR IS SET TO
C       ZERO       for normal return,
C       J          if the limit of 30*N iterations is exhausted
C                     while the J-th eigenvalue is being sought.
C
C     Transfer information back
C
      CALL QUIT('EIGRED not finished !')
      CALL MEMREL('EIGRED',WORK,1,KBCI,KFREE,LFREE)
      CALL QEXIT ('EIGRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrsnex */
      SUBROUTINE XRSNEX(IBTYP,IBCVC,IBEVC,IBPVC,RCNV,EVALR,EVECR,
     &                  GPCI,GPOE,GPOP,RNORM,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Purpose:
C       1) Construct residual:
C
C             LINEQ = T:  R(I) = (E[2] - W(I)*S[2])*X(I) + G
C             LINEQ = F:  R(I) = (E[2] - W(I)*S[2])*X(I)
C
C       2) Test for convergence of NEXCNV eigenvectors
C          Convergence criterium: |R(I)| .LE. THCXRS * |X(I)|
C
C       3) Generate new trial vectors
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(NREDM,*),GPCI(*),GPOE(*),GPOP(*),
     &          RNORM(7,NFREQ),WORK(*)
      CALL QENTER('XRSNEX')
      KFRSAV = KFREE
C
C       Memory allocation 
C
      IFAC = 2
      IF(STATIC) THEN
        NESIM = NFREQ
        NPSIM = NFREQ
        NCSIM = NFREQ
      ELSE
        NESIM = NFREQ*2
        NPSIM = NFREQ*2
        NCSIM = NFREQ*2
      ENDIF
      NTSIM = NESIM + NPSIM + NCSIM
      NBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
      IFREQ = 1
      MFREQ = NFREQ
      CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
      KFIRST = KBTOT
C
      KBOE   = KFIRST
      KFIRST = KFIRST + NZXOPEQ*NESIM

C
      KBOP   = KFIRST
      KFIRST = KFIRST + NZXOPPQ*NPSIM
C
      KBCI   = KFIRST
      KFIRST = KFIRST + NZCONFQ*NCSIM
C
C     Allocate buffers
C
      KFRBUF = KFREE
      CALL MEMGET('REAL',KBEBF,NZXOPEQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBPBF,NZXOPPQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBCBF,NZCONFQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KECNV,NESIM  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KPCNV,NPSIM  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCCNV,NCSIM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KEVEC,NESIM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KPVEC,NPSIM  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KCVEC,NCSIM  ,WORK,KFREE,LFREE)
C
C     Generate new trial vectors
C     ==========================
C
      CALL XRSNE1(WORK(KBOE),WORK(KBOP),WORK(KBCI),
     &            MFREQ,EVALR(IFREQ),EVECR(1,IFREQ),
     &            GPOE,GPOP,GPCI,IBTYP,IBCVC,IBEVC,IBPVC,
     &            WORK(KBEBF),WORK(KBPBF),WORK(KBCBF),RCNV,
     &            WORK(KECNV),WORK(KPCNV),WORK(KCCNV),
     &            WORK(KEVEC),WORK(KPVEC),WORK(KCVEC),
     &            RNORM(1,IFREQ),WORK,KFREE,LFREE)
C
C     Compress trial vectors in memory
C     ================================
C
      MBTOT = NBTOT
      NBTOT = NESIM*NZXOPEQ + NPSIM*NZXOPPQ + NCSIM*NZCONFQ
      IF(MBTOT.NE.NBTOT) THEN
C     ...compress vectors in memory
        KFIRST = KBTOT
C
        JBOE   = KBOE
        KBOE   = KFIRST
        NEDIM  = NZXOPEQ*NESIM
        KFIRST = KFIRST + NEDIM
        IF(KBOE.NE.JBOE.AND.NEDIM.GT.0) 
     &   CALL DCOPY(NEDIM,WORK(JBOE),1,WORK(KBOE),1)
C
        JBOP   = KBOP
        KBOP   = KFIRST
        NPDIM  = NZXOPPQ*NPSIM
        KFIRST = KFIRST + NPDIM
        IF(KBOP.NE.JBOP.AND.NPDIM.GT.0) 
     &    CALL DCOPY(NPDIM,WORK(JBOP),1,WORK(KBOP),1)
C
        JBCI   = KBCI
        KBCI   = KFIRST
        NCDIM  = NZCONFQ*NCSIM
        KFIRST = KFIRST + NCDIM
        IF(KBCI.NE.JBCI.AND.NCDIM.GT.0) 
     &    CALL DCOPY(NCDIM,WORK(JBCI),1,WORK(KBCI),1)
C
        CALL MEMREL('XRSNEX.CMP',WORK,1,KFRSAV,KFREE,LFREE)
        CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMREL('XRSNEX.BUF',WORK,1,KFRBUF,KFREE,LFREE)
      ENDIF
      CALL QEXIT('XRSNEX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRBVEC */
      SUBROUTINE PRBVEC(IUNIT,BVEC,NBSIM,NBDIM)
C*****************************************************************************
C     Print trial vector
C
C     Written by T.Saue Oct 24 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dgroup.h>
#include <dcbxrs.h>
      DIMENSION BVEC(NBDIM,NZ,NBSIM)
C
      IREP = JSYMOP - 1
      DO I = 1,NBSIM
        CALL PRQMAT(BVEC(1,1,I),NBDIM,1,NBDIM,1,NZ,
     &              IPQTOQ(1,IREP),IUNIT)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck BVINFO */
      SUBROUTINE BVINFO(ISIM,IBTYP,NBVEC)
C*****************************************************************************
C
C     Give information on number of trial vectors
C
C     Written by T.Saue Nov 8 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
#include <dcbibn.h>
      DIMENSION IBTYP(2,NBVEC),ISIM(6)
C
C     Initialize
C
      CALL IZERO(ISIM,6)
C
C     Count vectors
C
      DO I = 1,NBVEC
        IBV = IBTYP(1,I)
        IBH = IBTYP(2,I)
        IF(IBH.EQ.1) THEN
C       X(+)-vectors
          IBH = 0
        ELSE
C       X(-)-vectors
          IBH = 3
        ENDIF
        IF    (IBV.EQ.JBENDX) THEN
          IBV = 1
        ELSEIF(IBV.EQ.JBPNDX) THEN
          IBV = 2
        ELSEIF(IBV.EQ.JBCNDX) THEN
          IBV = 3
        ENDIF
        IBV = IBV + IBH
        ISIM(IBV) = ISIM(IBV) + 1
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRIRMT */
      SUBROUTINE PRIRMT (REDMAT,LDA,NRED,IBTYP)
C*****************************************************************************
C
C     Print reduced matrix with information
C
C     Written by T.Saue Nov 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER BEGIN
      DIMENSION REDMAT(LDA,NRED),IBTYP(2,NRED)
      CHARACTER BTP(-1:2)*1,HTP(-1:1)*1
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA BTP/'C',' ','E','P'/,HTP/'-',' ','+'/
C
      NCTL = 1
C
C     First check if zero matrix
C
      AMAX = D0
      DO 10 J = 1,NRED
        DO 10 I = 1,NRED
          AMAX = MAX(AMAX,ABS(REDMAT(I,J)))
   10 CONTINUE
      IF (AMAX .EQ. D0) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero matrix.'
         GO TO 3
      END IF
      KCOL = 4
C
      LAST = MIN(NRED,KCOL)
      DO 2 BEGIN = 1,NRED,KCOL
         WRITE (LUPRI,'(/18X,4(1X,A3,I4,1X,4A1,2X))')
     &       ('Col',I,'(',BTP(IBTYP(1,I)),HTP(IBTYP(2,I)),')',
     &       I = BEGIN,LAST)
         DO 1 K = 1,NRED
            DO 4 I = BEGIN,LAST
               IF (REDMAT(K,I).NE.D0) GO TO 5
    4       CONTINUE
         GO TO 1
    5       WRITE (LUPRI,'(4X,A3,I4,1X,4A1,2X,1P,8D15.7)') 
     &        'Row',K,'(',BTP(IBTYP(1,K)),HTP(IBTYP(2,K)),')',
     &        (REDMAT(K,I), I = BEGIN,LAST)
    1    CONTINUE
    2 LAST = MIN(LAST+KCOL,NRED)
    3 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRIRVC */
      SUBROUTINE PRIRVC (REDVEC,NRED,IBTYP)
C*****************************************************************************
C
C     Print reduced vector with information
C
C     Written by T.Saue Nov 11 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION REDVEC(NRED),IBTYP(2,NRED)
      CHARACTER BTP(-1:2)*1,HTP(-1:1)*1,PFMT*26
      PARAMETER (FFMIN=1.D-3, FFMAX = 1.D3)
      DATA BTP/'C',' ','E','P'/,HTP/'-',' ','+'/
C
C
C     First check if zero matrix
C
      AMAX = D0
      DO 10 I = 1,NRED
        AMAX = MAX(AMAX,ABS(REDVEC(I)))
   10 CONTINUE
      IF (AMAX .EQ. D0) THEN
         WRITE (LUPRI,'(/T6,A)') 'Zero vector.'
         GO TO 3
      END IF
C
C     Set print format
C
      IF (FFMIN.LE.AMAX .AND. AMAX.LE.FFMAX) THEN
C        use F output format
         PFMT = '(4X,A3,I4,1X,4A1,2X,F15.8)'
      ELSE
C        use 1PD output format
         PFMT = '(4X,A3,I4,1X,4A1,2X,D15.6)'
      END IF
C
      DO 2 I = 1,NRED
        IF(REDVEC(I).EQ.D0) GOTO 2
        WRITE (LUPRI,PFMT) 
     &        'Row',I,'(',BTP(IBTYP(1,I)),HTP(IBTYP(2,I)),')',
     &        REDVEC(I)
    2 CONTINUE
    3 RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSOUT */
      SUBROUTINE XRSOUT(RCNV,RNORM,EVALR,IBTYP)
C*****************************************************************************
C
C     Final output from XRSCTL
C       RCNV  - total norm of residual relative to norm of solution vector
C       RNORM - residual norm
C             1: norm of E(+)-part
C             2: norm of P(+)-part
C             3: norm of C(+)-part
C             4: norm of E(-)-part
C             5: norm of P(-)-part
C             6: norm of C(-)-part
C
C     Written by T.Saue Nov 19 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
C
      CHARACTER CPUTID*12,WORD*3
      DIMENSION RNORM(7,NFREQ),EVALR(NFREQ),
     &          RCNV(NFREQ),IBTYP(2,NREDM),ISIM(6)
C
      CALL HEADER('Final output from XRSCTL',-1)
C
C     Check for state of the calculation
C     ==================================
C
      IF (KCONV.EQ.-3) THEN
C
C     Start vectors rejected
C     ----------------------
C
        WRITE (LUPRI,'(/A/4X,A)')
     &      '*** WARNING - Start vectors rejected !',
     &      'No microiterations.'
        RETURN
      ELSE IF (KCONV.EQ.-2) THEN
C
C     Maximum dimension of reduced space exceeded
C     -------------------------------------------
C
        WRITE (LUPRI,'(/A,I3,A/4X,A)')
     &      '*** WARNING - Microiterations stopped after',
     &      ITMIC,' iterations',
     &      'because maximum dimension of reduced space exceeded.'
      ELSE IF (KCONV.EQ.-1) THEN
C
C     Linear dependence between new trial vectors
C     -------------------------------------------
C
        WRITE (LUPRI,'(/A,I3,A/4X,A)')
     &      '*** WARNING - Microiterations stopped after',
     &      ITMIC,' iterations',
     &      'because of linear dependence between new trial vectors'
      ELSE IF(KCONV.EQ.0)THEN
C
C     Convergence of iterative solution
C     ---------------------------------
C
        WRITE(LUPRI,'(/A,I3,A)')
     &     '*** Microiterations converged after',ITMIC,
     &     ' microiterations'
      ELSE
C
C           No convergence
C           --------------
C
        IF (ITMIC.GE.MAXITR) THEN
C     
C         Max no of microiterations reached
C         ---------------------------------
C
c NWARN is not used anywhere./Feb 9 1998 - jth
c         NWARN = NWARN + 1
          WRITE(LUPRI,'(/A,I4,A)')
     *      ' *** WARNING-Maximum number of microiterations,',
     *      ITMIC,', reached'
        END IF
      END IF
C
C     Give summary of calculations
C
      REWIND LUCYCL
      WRITE(LUPRI,'(/12X,A9,3X,A10,2X,6(1X,A9,4X))')
     &  'Freq.(au)','Total norm',' E+ norm ',' P+ norm ',' C+ norm ',
     &  ' E- norm ',' P- norm ',' C- norm '
      DO ITR = 1,ITMIC
        CALL PRSYMB(LUPRI,'=',120,0)
        READ(LUCYCL,'(I3,2X,A9,6(2X,A3,I3))')
     &     JTR,INTTYP,(WORD,ISIM(I),I=1,6)
        READ(LUCYCL,'(6E10.2)') ((RNORM(I,J),I=1,6),J=1,NFREQ)
        READ(LUCYCL,'(A12)') CPUTID
        WRITE(LUPRI,'(A5,I3,3X,A12,A,A9,A)') 
     &        '* It.',ITR,CPUTID,' (',INTTYP,')'
        CALL PRSYMB(LUPRI,'-',120,0)
        DO IFREQ = 1,NFREQ
           WRITE(LUPRI,'(10X,F10.6,A3,E10.2,A3,6(E10.2,A1,I2,A1))') 
     &     EVALR(IFREQ),' : ',RNORM(1,IFREQ),' | ',
     &     (RNORM(I,IFREQ),'(',ISIM(I-1),')',I=2,7)
        ENDDO
      ENDDO
      CALL PRSYMB(LUPRI,'=',120,0)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2DIAG */
      SUBROUTINE E2DIAG(HEDIAG,HPDIAG,HCDIAG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Get diagonal elements of Hessian(approximate)
C
C     Written by T.Saue 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
      DIMENSION HEDIAG(*),HPDIAG(*),HCDIAG(*),WORK(*)
C
C     Orbital part
C
      IF((NZXOPE.GT.0).OR.(NZXOPP.GT.0)) THEN
        CALL E2ODIA(HEDIAG,HPDIAG,WORK(KZXOPE),WORK(KZXOPP),
     &              WORK,KFREE,LFREE)
      ENDIF
C
C     Configurational part
C
      IF(NZCONF.GT.0) THEN
C        CALL E2CDIA(HCDIAG,WORK,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2ODIA */
      SUBROUTINE E2ODIA(HEDIAG,HPDIAG,JXOPE,JXOPP,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Get orbital part of diagonal hessian elements
C
C     In RPA the diagonal elements are given by:
C
C       E[2}_{si,si} = F_{tt} - F_{ii} - [g_{tt,ii} - g_{ti,it}] 
C     
C     which is approximated by
C
C       E[2}_{si,si} ~ F_{tt} - F_{ii}
C
C     Written by T.Saue Dec 1 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DM2 = -2.0D0)
C
#include <dcbgen.h>
#include <dcbxrs.h>
#include <dcbbas.h>
#include <dcbibn.h>
C
      DIMENSION HEDIAG(*),HPDIAG(*),JXOPE(2,*),JXOPP(2,*),WORK(*)
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KFMO,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Get total Fock matrix in MO-basis
C
      CALL GETFCK(WORK(KFMO),IPRXRS,WORK,KFREE,LFREE)
C
C     Electronic part
C
      IF(NZXOPE.GT.0) THEN
        CALL E2ODI1(HEDIAG,WORK(KFMO),JXOPE,NZXOPE)
        IF(IPRXRS.GE.5) THEN
          WRITE(LUPRI,'(/A)') '** E2ODIA: Diagonal Hessian (e-e):'
          WRITE(LUPRI,'(2I5,3X,1P,E18.10)')
     &      (JXOPE(1,I),JXOPE(2,I),HEDIAG(I), I = 1,NZXOPE)
        ENDIF
      ENDIF
C
C     Positronic part
C
      IF(NZXOPP.GT.0) THEN
        IF(STERNH) THEN
          FAC = DM2*CVAL*CVAL
          CALL DCOPY(NZXOPP,FAC,0,HPDIAG,1)
        ELSE
          CALL E2ODI1(HPDIAG,WORK(KFMO),JXOPP,NZXOPP)
        ENDIF
        IF(IPRXRS.GE.5) THEN
          WRITE(LUPRI,'(/A)') '** E2ODIA: Diagonal Hessian (e-p):'
          WRITE(LUPRI,'(2I5,3X,1P,E18.10)')
     &      (JXOPP(1,I),JXOPP(2,I),HPDIAG(I), I = 1,NZXOPP)
        ENDIF
      ENDIF
C
C
      CALL MEMREL('XRSSEP',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2ODI1 */
      SUBROUTINE E2ODI1(HODIAG,FMO,JXOP,NXOP)
C*****************************************************************************
C
C     Get orbital part of diagonal hessian elements
C
C     In RPA the diagonal elements are given by:
C
C       E[2}_{si,si} = F_{tt} - F_{ii} - [g_{tt,ii} - g_{ti,it}] 
C     
C     which is approximated by
C
C       E[2}_{si,si} ~ F_{tt} - F_{ii}
C
C     Written by T.Saue Dec 1 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C     
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION HODIAG(NXOP),FMO(NORBT,NORBT,NZ),
     &          JXOP(2,NXOP)
C
      DO IXOP = 1,NXOP
        IS = JXOP(2,IXOP)
        II = JXOP(1,IXOP)
        HODIAG(IXOP) = FMO(IS,IS,1) - FMO(II,II,1)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSST2 */
      SUBROUTINE XRSST2(TYP,BVEC,GP,HDIAG,EVALR,IBTYP,IBVEC,IVEC,
     &                  JBTYP,NBVEC,IBOFF,NBPREV,NVAR,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate start vectors
C
C     Written by T.Saue Dec 1 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1 = -1.0D0)
C
#include <dcbxrs.h>
#include <dgroup.h>
C
      LOGICAL FNDXLB
      CHARACTER TYP*1
      DIMENSION BVEC(NVAR,NZ,*),GP(NVAR,NZ),HDIAG(NVAR),
     &          EVALR(NFREQ),IBTYP(2,*),IBVEC(*),IVEC(*),
     &          WORK(*)
#include <chrnos.h>
C
      NVARQ = NVAR*NZ
      IF(LINEQ) THEN
C
C     L I N E A R   R E S P O N S E   F U N C T I O N S
C     =================================================
C
        IF(RSREST) THEN
          CALL GETDAT(XRSLBL(1),XRSLBL(2))
          WRITE(XRSLBL(2)(1:1),'(A1)') TYP
          WRITE(XRSLBL(2)(2:2),'(A1)') '+'
          DO IP = 1,NFREQ
            XRSLBL(2)(2:2) = '+'
            XRSLBL(2)(3:3) = CHRNOS(IP/10)
            XRSLBL(2)(4:4) = CHRNOS(MOD(IP,10))
            IM   = IP + NFREQ
            FREQ = EVALR(IP)
            IF(FREQ.EQ.D0) THEN
              IF(FNDXLB(XRSLBL,LURST)) THEN
                WRITE(LUPRI,'(A)') '*** WARNING ***'
                WRITE(LUPRI,'(A,A2,A/3(3X,A8))')
     &            'Restart on ',XRSLBL(2)(1:2),'-vector:',
     &            (XRSLBL(I),I=1,3)
                CALL READT(LURST,NVARQ,BVEC(1,1,IP))
              ELSE
                CALL DCOPY(NVARQ,GP,1,BVEC(1,1,IP),1)
                IF(DIAGHE) THEN
                  DO IZ = 1,NZ
                    DO IVAR = 1,NVAR
                      BVEC(IVAR,IZ,IP) = (DM1/HDIAG(IVAR))
     &                                  *BVEC(IVAR,IZ,IP)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
              IBTYP(1,IBOFF+IP) = JBTYP
              IBTYP(2,IBOFF+IP) = IPX
              IBVEC(NBPREV+IP)  = IBOFF + IP
              IVEC(IP)         = 1        
              IF(.NOT.STATIC) THEN
                IBTYP(1,IBOFF+IM) = JBTYP
                IBTYP(2,IBOFF+IM) = IMX
                IBVEC(NBPREV+IM)  = IBOFF + IM
                IVEC(IM)         = 0
              ENDIF
            ELSE
              IF(FNDXLB(XRSLBL,LURST)) THEN
                WRITE(LUPRI,'(A)') '*** WARNING ***'
                WRITE(LUPRI,'(A,A2,A/3(3X,A8))')
     &            'Restart on ',XRSLBL(2)(1:2),'-vector:',
     &            (XRSLBL(I),I=1,3)
                CALL READT(LURST,NVARQ,BVEC(1,1,IP))
              ELSE
                CALL DCOPY(NVARQ,GP,1,BVEC(1,1,IP),1)
                IF(DIAGHE) THEN
                  DO IZ = 1,NZ
                    DO IVAR = 1,NVAR
                      FAC = (FREQ**2-HDIAG(IVAR)**2)
                      BVEC(IVAR,IZ,IP) = (HDIAG(IVAR)/FAC)
     &                                 *BVEC(IVAR,IZ,IP)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
              IBTYP(1,IBOFF+IM)     = JBTYP
              IBTYP(2,IBOFF+IM)     = IMX
              IBVEC(NBPREV+IM)      = IBOFF + IM
              IVEC(IM)              = 1
C
C             X(-) - vector
C
              XRSLBL(2)(2:2) = '-'
              IF(FNDXLB(XRSLBL,LURST)) THEN
                WRITE(LUPRI,'(A)') '*** WARNING ***'
                WRITE(LUPRI,'(A,A2,A/3(3X,A8))')
     &            'Restart on ',XRSLBL(2)(1:2),'-vector:',
     &            (XRSLBL(I),I=1,3)
                CALL READT(LURST,NVARQ,BVEC(1,1,IM))
              ELSE
                CALL DCOPY(NVARQ,GP,1,BVEC(1,1,IP),1)
                IF(DIAGHE) THEN
                  DO IZ = 1,NZ
                    DO IVAR = 1,NVAR
                      FAC = (FREQ**2-HDIAG(IVAR)**2)
                      BVEC(IVAR,IZ,IM) = (FREQ/FAC)
     &                                 *BVEC(IVAR,IZ,IM)
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF
              IBTYP(1,IBOFF+IM)     = JBTYP
              IBTYP(2,IBOFF+IM)     = IMX
              IBVEC(NBPREV+IM)      = IBOFF + IM
              IVEC(IM)              = 1
            ENDIF
          ENDDO
        ELSEIF(DIAGHE) THEN
          DO IP = 1,NFREQ
            IM = IP + NFREQ
            FREQ = EVALR(IP)
            IF(FREQ.EQ.D0) THEN
              CALL DCOPY(NVARQ,GP,1,BVEC(1,1,IP),1)
              DO IZ = 1,NZ
                DO IVAR = 1,NVAR
                  BVEC(IVAR,IZ,IP) = (DM1/HDIAG(IVAR))*BVEC(IVAR,IZ,IP)
                ENDDO
              ENDDO
              IBTYP(1,IBOFF+IP) = JBTYP
              IBTYP(2,IBOFF+IP) = IPX
              IBVEC(NBPREV+IP)  = IBOFF + IP
              IVEC(IP)         = 1
              IBTYP(1,IBOFF+IM) = JBTYP
              IBTYP(2,IBOFF+IM) = IMX
              IBVEC(NBPREV+IM)  = IBOFF + IM
              IVEC(IM)         = 0
            ELSE
              CALL DCOPY(NVARQ,GP,1,BVEC(1,1,IP),1)
              CALL DCOPY(NVARQ,GP,1,BVEC(1,1,IM),1)
              DO IZ = 1,NZ
                DO IVAR = 1,NVAR
                  FAC = (FREQ**2-HDIAG(IVAR)**2)
                  BVEC(IVAR,IZ,IP) = (HDIAG(IVAR)/FAC)*BVEC(IVAR,IZ,IP)
                  BVEC(IVAR,IZ,IM) = (FREQ/FAC)*BVEC(IVAR,IZ,IM)
                ENDDO
              ENDDO
              IBTYP(1,IBOFF+IP)     = JBTYP
              IBTYP(2,IBOFF+IP)     = IPX
              IBVEC(NBPREV+IP)       = IBOFF + IP
              IVEC(IP)             = 1
              IBTYP(1,IBOFF+IM)     = JBTYP
              IBTYP(2,IBOFF+IM)     = IMX
              IBVEC(NBPREV+IM)      = IBOFF + IM
              IVEC(IM)              = 1
            ENDIF
          ENDDO          
        ELSE
          IF(STATIC) THEN
            CALL DCOPY(NVARQ,GP,1,BVEC(1,1,1),1)
            IBTYP(1,IBOFF+1) = JBTYP
            IBTYP(2,IBOFF+1) = IPX
            IBVEC(NBPREV+1)  = IBOFF + 1
            IVEC(1)         = 1
          ELSE
            CALL DCOPY(NVARQ,GP,1,BVEC(1,1,1),1)
            IBTYP(1,IBOFF+1) = JBTYP
            IBTYP(2,IBOFF+1) = IPX
            IBVEC(NBPREV+1)  = IBOFF + 1
            IVEC(1)         = 1
            CALL DCOPY(NVARQ,GP,1,BVEC(1,1,2),1)
            IBTYP(1,IBOFF+2) = JBTYP
            IBTYP(2,IBOFF+2) = IMX
            IBVEC(NBPREV+2)  = IBOFF + 2
            IVEC(2)         = 1
          ENDIF
        ENDIF
        CALL XRSORT(BVEC,JBTYP,NBVEC,NBPREV,
     &            IBVEC,IBTYP,IVEC,WORK,KFREE,LFREE)
      ELSE
C
C     E X C I T A T I O N  E N E R G I E S etc.
C     ========================================
C
         CALL QUIT(
     &    ' Eigenvalue equation (PP) not implemented yet in XRSST.')
      END IF
C
      RETURN
      END
C  /* Deck fndlb3 */
      LOGICAL FUNCTION FNDXLB(RTNLBL,LU)
C
C  5-Aug-1986 hjaaj
C  Feb 1997 - tsaue
C
      CHARACTER*8 RTNLBL(3), B(4), STAR8
      PARAMETER (STAR8 = '********')
      REWIND LU
      IRDERR = 0
    1 READ (LU,END=3,ERR=6) B
      IRDERR = 0
      IF (B(1).NE.STAR8)     GO TO 1
      IF (B(4).NE.RTNLBL(3)) GO TO 1
      IF (B(3)(1:4).NE.RTNLBL(2)(1:4)) GO TO 1
      FNDXLB    = .TRUE.
      RTNLBL(1) = B(2)
      RTNLBL(2) = B(3)
      GO TO 10
C
    6 IRDERR = IRDERR + 1
      IF (IRDERR .LT. 5) GO TO 1
#if defined (SYS_IBM)
C 880703-hjaaj -- attempt to fix an IBM problem
C IBM gave permanent io error on write after calling FNDLAB for
C     a non-existing file
      CLOSE(LU,ERR=3)
#endif
      GO TO 8
    3 CONTINUE
#if defined (SYS_IBM) || defined (SYS_CRAY) || defined (VAR_MFDS) || defined (SYS_T3D)
C 880916-hjaaj -- attempt to fix an IBM problem
C IBM shifts to new file after END= branch (e.g. FTxxF002 instead
C     of FTxxF001), backspace makes LU ready for append.
C 940510-hjaaj: same change for Cray's multifile datasets
      BACKSPACE LU
#endif
    8 FNDXLB = .FALSE.
C
   10 RETURN
      END
