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
C...   group of the licence holder may not take any part of Dalton,
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
C FILE    : krmcxnr.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xnrctl */
      SUBROUTINE XNRCTL(GPCI,GPOE,GPOP,IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &                  EVALR,EVECR,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Module for solving the NR equations from KR-MCSCF second order
C     optimization. Build on XRSCTL from prp/rsp
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
C     where h is the operator manifold. In a KR-MCSCF optimziation
C     h is limited to single excitations.
C
C-----------------------------------------------------------------------------
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the corresponding reduced equations
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
#include <dcbxnr.h>
#include <blocks.h>
C
      CALL QENTER('XNRCTL')
      KFRSAV = KFREE
      NREDGP = NREDM
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
C
C     Fix the number of variational parameters
C
      NVPAR = 0
      IF(NZXOPE.GT.0) THEN
        NVPAR        = NVPAR + 1
        IPEP         = NVPAR
        BVTYP(NVPAR) = 'E+'
      ENDIF
      IF(NZXOPP.GT.0) THEN
        NVPAR        = NVPAR + 1
        IPPP         = NVPAR
        BVTYP(NVPAR) = 'P+'
      ENDIF
      IF(NZCONF.GT.0) THEN
        NVPAR        = NVPAR + 1
        IPCP         = NVPAR
        BVTYP(NVPAR) = 'C+'
      ENDIF
      CALL MEMGET('REAL',KRNORM,NVPAR,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KISIM ,NVPAR,WORK,KFREE,LFREE)
C
      CALL XNRCT1(GPCI,GPOE,GPOP,RCNV,
     &            IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &            WORK(KREDE),WORK(KREDS),WORK(KREDGP),
     &            WORK(KRNORM),WORK(KISIM),WORK(KPOS),
     &            WORK,KFREE,LFREE)
C
      CALL MEMREL('XNRCTL.XNRRED',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XNRCTL')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XNRCT1 */
      SUBROUTINE XNRCT1(GPCI,GPOE,GPOP,RCNV,
     &                  IBTYP,IBCVC,IBEVC,IBPVC,EVALR,EVECR,
     &                  REDE,REDS,REDGP,RNORM,ISIM,
     &                  NPOS,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Module for solving the NR equations from KR-MCSCF second order
C     optimization. Build on XRSCTL from prp/rsp
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
C     where h is the operator manifold. In a KR-MCSCF optimziation
C     h is limited to single excitations.
C
C-----------------------------------------------------------------------------
C     The linear reponse equation is solved by projecting X onto
C     a reduced space of trial vectors:
C
C        X = (b_i)a_i           B = [(b_1)  (b_2)  ...  (b_n)]
C
C     and solving the corresponding reduced equations
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
#include <dcbxnr.h>
#include <pgroup.h>
#include <dgroup.h>
#include <dcborb.h>
#include <dcbkrmc_itinfo.h>
C
      CHARACTER CPUTID*12,SECTID*12,DAYTID*24
      DIMENSION GPCI(*), GPOE(*),GPOP(*),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(NREDM),
     &          RNORM(NVPAR,NFREQ),ISIM(NVPAR),NPOS(*),WORK(*)
C
#include <dcbibn.h>
      LOGICAL BINMEM
      CHARACTER LAB123(3)*8,RESTAR*8
C
      DATA LAB123/'********','********','XRSPONPP'/
      DATA RESTAR/'RESTART '/
C
      CALL QENTER('XNRCT1')
      CALL GETTIM(CPUTOT1,CPUWALL1)
      KFRSAV = KFREE
C
C     Determine maximum number of simultaneous trial vectors
C     and maximum number of vectors per batch in XNRNEX
C     ======================================================
C
      MAXVEC = MAX(NEXSTV,NEXSIM)
      MAXSIM = MAXVEC
      MAXNEX = MAXVEC
      NWRSP  = LFREE
C
      IF (IPRXNR.GE.10) THEN
         WRITE(LUPRI,'(2(/A,I8))')
     *     'Maximum number of simultaneous trial vectors :',MAXSIM,
     *     'Maximum number of vectors per batch in XNRNEX:',MAXNEX
      END IF
C
C
C     ***************************
C     *** I n i t i a l i z e ***
C     ***************************
C
      KFRITR = KFREE
C
      CALL XNRINI(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
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
         CALL GTINFO(DAYTID)
         WRITE(LUPRI,'(A,I4,A,3X,A24/)')
     +   '########## START MICROITERATION NO.',ITMIC,' ##########',
     +   DAYTID
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
          EVECR(1) = D1
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
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(A,I4,A,3X,A24/)')
     +  '########## START MICROITERATION NO.',ITMIC,' ##########',
     +  DAYTID
C
CTROND       Check BINMEM !!!!!
C
        IF (IPRNRS.GE.5) WRITE(LUPRI,'(A,L1)')
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
        IF (IPRXNR.GE.2 .AND. ITLOAD.GT.1) 
     &    WRITE(LUPRI,'(/A,I3)') '* Load no.',ITLOAD
C       Trial vectors not in memory; read from file
        IF(.NOT.BINMEM) THEN
          NWRST  = NWRSP
          CALL NREABLD(NWRST,NBCSIM,NBESIM,NBPSIM,
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
        IF(IPRNRS.GE.10) THEN
          WRITE(LUPRI,'(A,I5)')
     &      'Load no.      :',ITLOAD,
     &      'C-type vectors:',NCSIM,
     &      'E-type vectors:',NESIM,
     &      'P-type vectors:',NPSIM
        ENDIF
C
        IF(IPRNRS.GE.5) THEN
          WRITE(LUPRI,'(/A)')' Before XNRSVC'
          IF(NCSIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NCSIM,' Config. trial vectors'
            CALL NRPRBVEC(LUPRI,WORK(KBCI),NCSIM,NZCONF)
          ENDIF
          IF(NESIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NESIM,' Electron trial vectors'
            CALL NRPRBVEC(LUPRI,WORK(KBOE),NESIM,NZXOPE)
          ENDIF
          IF(NPSIM.GT.0) THEN
            WRITE(LUPRI,'(/I5,A)')NPSIM,' Positron trial vectors'
            CALL NRPRBVEC(LUPRI,WORK(KBOP),NPSIM,NZXOPP)
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
        CALL XNRSVC(WORK(KBCI),WORK(KBOE),WORK(KBOP),
     &              WORK(KSCI),WORK(KSOE),WORK(KSOP),
     &              IBTYP,IBCVC,IBEVC,IBPVC,NPOS,
     &              WORK,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
C
C       Write sigma-vectors to file 
C       ============================
C
        ICOUNT = NZRED
        CALL NRWRSVEC(ICOUNT,WORK(KSCI),WORK(KSOE),WORK(KSOP))
C
C       Update restart info file
C       ========================
C
        REWIND (LURSP)
        KTOT = NZRED + LOFFTY
        IF (IPRNRS.GE.7) THEN
          WRITE(LUPRI,'(/A,I5)')' Trial vectors for restart:',KTOT
        ENDIF
        WRITE(LURSP) LAB123,RESTAR
        WRITE(LURSP) JSYMOP,JTIMOP,KTOT,
     &               (IBTYP(1,I),IBTYP(2,I),I=1,KTOT)
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
     *       ' >>> ERROR IN XNRCTL >>>',
     *       ' Dimension of reduced space is  ',NZRED,
     *       ' which exceeds allowed dimension',NREDM
         CALL QUIT('XNRCTL: Too large dimension of reduced space')
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
        CALL XNRRED (REDGP,REDE,REDS,GPCI,GPOE,GPOP,
     &               IBCVC,IBEVC,IBPVC,IBTYP,
     &               WORK(KBCI),WORK(KBOE),WORK(KBOP),
     &               WORK(KSCI),WORK(KSOE),WORK(KSOP),
     &               WORK,KFREE,LFREE)
        CALL MEMREL('XNRCT1.XNRRED',WORK,1,KFRITR,KFREE,LFREE)
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
        CALL NRBVINFO(ISIM,IBTYP(1,IZRED+1),NTSIM)
C
C       Write info to PAMCYC
C       ====================
C
        WRITE(LUCYCL,'(I3,2X,A9,6(2X,A2,A1,I3))')
     &     ITMIC,INTTYP,(BVTYP(I),':',ISIM(I),I=1,NVPAR)
        WRITE(LUPRI,'(A,6(2X,A2,A1,I3))') 
     &     '* New trial vectors: ',(BVTYP(I),':',ISIM(I),I=1,NVPAR)
C
C       Solve reduced RSP equation
C       ==========================
C
        CALL XNRSOL (REDGP,REDE,REDS,EVALR,EVECR,RCNV,
     &               IBTYP,WORK,KFREE,LFREE)
        CALL FLSHFO(LUPRI)
C
C        Generate NEXSIM linear independent trial vectors for
C        next microiteration using solutions of reduced RSP equation
C        ===========================================================
C
         CALL XNRNEX(IBTYP,IBCVC,IBEVC,IBPVC,RCNV,
     &               EVALR,EVECR,GPCI,GPOE,GPOP,
     &               RNORM,WORK,KFREE,LFREE)
C
         TIMITR = SECOND() - TIMITR
         CPUTID = SECTID(TIMITR)
C
C        Write info to PAMCYC
C        ====================
C
         WRITE(LUCYCL,'(6E10.2)') ((RNORM(I,J),I=1,NVPAR),J=1,NFREQ)
         WRITE(LUCYCL,'(A12)') CPUTID
         IF(IPRXNR.GE.1) THEN
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
      CALL XNROUT(RCNV,RNORM,ISIM,EVALR,IBTYP)
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
      CALL MEMREL('XNRCT1',WORK,1,KFRSAV,KFREE,LFREE)
C
C
C     For KR-MCSCF iterations, calculate predicted energy change
C
      CALL GETTIM(CPUTOT2,CPUWALL2)
      CALL RSNRPRD(REDGP,REDE,EVECR,NREDM)
      DINFO(17) = CPUTOT2-CPUTOT1
      IINFO(2) = ITMIC
      IINFO(3) = NZRED
      CALL FLSHFO(LUPRI)
      CALL QEXIT('XNRCT1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rsnrprd */
      SUBROUTINE RSNRPRD(REDGP,REDE,EVECR,NREDM)
C*****************************************************************************
C
C     Calculate second order change in energy
C
C     Written by J. Thyssen - Nov 27 1998
C     Last revision : 
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (DP5 = 0.50D00, D0 = 0.00D00, DP25 = 0.25D00)
C
#include <dcbopt.h>
C
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),EVECR(NREDM)
C
      CALL HEADER('Output from RSNRPRD',-1)
      CALL HEADER('solution vector',-1)
      CALL OUTPUT(EVECR,1,NREDM,1,1,NREDM,1,1,LUPRI)
      CALL HEADER('gradient',-1)
      CALL OUTPUT(REDGP,1,NREDM,1,1,NREDM,1,1,LUPRI)
      CALL HEADER('hessian',-1)
      CALL OUTPUT(REDE,1,NREDM,1,NREDM,NREDM,NREDM,1,LUPRI)
      SUM = D0
      DO I = 1,NREDM
         SUM1 = D0
         DO J = 1,NREDM
            SUM1 = SUM1 + REDE(I,J)*EVECR(J)
         END DO
         SUM = SUM + EVECR(I)*(-REDGP(I) + DP5 * SUM1)
      END DO
      write(6,*) 'rsnrprd: predicted energy change: ',sum
      DEPRED = DP5*SUM
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XNRST */
      SUBROUTINE XNRST(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
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
#include <dcbxnr.h>
#include <dgroup.h>
C
      CALL QENTER('XNRST')
      KFRSAV = KFREE
C
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
      IF(DIAGHE.OR.RSREST) THEN
         NVEC = NFREQ
      ELSE
         NVEC = 1
      ENDIF
      IFAC = 1
      IF(NZXOPE.GT.0) NESIM = NVEC*IFAC
      IF(NZXOPP.GT.0) NPSIM = NVEC*IFAC
      IF(NZCONF.GT.0) NCSIM = NVEC*IFAC
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
C
C     If requested, allocate space for diagonal 
C     elements of Hessian (approximate)
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
      CALL  XNRST1(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
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
        CALL MEMREL('XNRINI.CMP',WORK,1,KFRSAV,KFREE,LFREE)
        CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMREL('XNRINI.BUF',WORK,1,KFRBUF,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT ('XNRST')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XNRST1 */
      SUBROUTINE XNRST1(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &                 IEVEC,IPVEC,ICVEC,GPCI,GPOE,GPOP,EVALR,
     &                 HEDIAG,HPDIAG,HCDIAG,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate start vectors
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision : July 30 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
C
#include <dcbibn.h>
C
#include <dcbxnr.h>
#include <dgroup.h>
      LOGICAL   BINMEM
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),EVALR(*),
     &          GPCI(*),GPOE(*),GPOP(*),IEVEC(*),IPVEC(*),ICVEC(*),
     &          HEDIAG(*),HPDIAG(*),HCDIAG(*),WORK(*)
C
      CALL QENTER('XNRST1')
      KFRSAV = KFREE
C
C     If requested, get diagonal elements of Hessian (approximate)
C
      IF(DIAGHE) CALL NRE2DIAG(HEDIAG,HPDIAG,HCDIAG,
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
        CALL XNRST2('E',WORK(KBOE),GPOE,HEDIAG,EVALR,IBTYP,IBEVC,IEVEC,
     &              JBENDX,NESIM,IND,NERED,NZXOPE,WORK,KFREE,LFREE)
        IND = IND + NESIM
      END IF
C
C     Orbital (e-p) type trial vector
C     -------------------------------
C
      IF(NZXOPP.GT.0) THEN
        CALL XNRST2('P',WORK(KBOP),GPOP,HPDIAG,EVALR,IBTYP,IBPVC,IPVEC,
     &              JBPNDX,NPSIM,IND,NPRED,NZXOPP,WORK,KFREE,LFREE)
        IND = IND + NPSIM
      END IF
C
C     Configurational type trial vector
C     ---------------------------------
C
      IF(NZCONF.GT.0) THEN
        CALL XNRST2('C',WORK(KBCI),GPCI,HCDIAG,EVALR,IBTYP,IBCVC,ICVEC,
     &              JBCNDX,NCSIM,IND,NCRED,NZCONF,WORK,KFREE,LFREE)
        IND = IND + NCSIM
      END IF
C
      BINMEM = .TRUE.
C     HJMAERKE: insert appropriate test for BINMEM
C     CTROND: also reset KFREE for BINMEM = F !!!
      CALL QEXIT ('XNRST1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XNRORT */
      SUBROUTINE XNRORT(BVECS,NBTYP,NBNEW,NBPREV,
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
C  WRKXNR: LZCONF,LZXOPT,JSYMOP
C  XCBPRI: LUPRI
C
#include <dgroup.h>
#include <dcbxnr.h>
C
      CALL QENTER('XNRORT')
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
        WRITE(LUPRI,'(A,A,I5)') 'XNRORT ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('XNRORT: Unknown NBTYP !')
      ENDIF
      THRLDV = SQRT(THRLDV)
C
C     Memory allocation
C
      CALL MEMGET('REAL',KBVBUF,NBDIM,WORK,KFREE,LFREE)      
      CALL MEMGET('REAL',KBVNRM,NBNEW,WORK,KFREE,LFREE)      
      CALL XNROR1(BVECS,NBDIM,NBNEW,NBPREV,THRLDV,IBVEC,IBTYP,
     &            WORK(KBVBUF),WORK(KBVNRM),LUBVEC,IVECS,NBTYP)
      CALL MEMREL('XNRORT',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('XNRORT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck XNROR1 */
      SUBROUTINE XNROR1(BVECS,NBDIM,NBNEW,NBPREV,THRLDP,IBVEC,IBTYP,
     &                  BVBUF,BVNRM,LUBVEC,IVECS,NBTYP)
C*****************************************************************************
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D2 = 2.0D0,D4 = 4.0D0)
      PARAMETER (T1MIN = 1.D-8, THRRND = 1.D-4)
C
#include <dcbxnr.h>
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
     &      'XNRORT: new ',BTYP(NBTYP),XTYP(IH),
     &      '-vector no. ',IBN,' is removed ',
     &      'because of linear dependence: ',BNRM
          IVECS(IBN) = 0
        ELSEIF(BNRM.LT.THRRND*BVNRM(IBN) .AND. ITURN.EQ.1) THEN
          IROUND = IROUND + 1
        ENDIF        
        IF(IPRXNR.GE.2) THEN
          WRITE (LUPRI,'(A,A1,A3,A,I6,5X,A,1P,D12.5)')
     &      'XNRORT: new ',BTYP(NBTYP),XTYP(IH),
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
     &     'XNRORT error: IROUND .gt. 0 and ITURN .gt. 1')
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
C/* Deck NReabld */
      SUBROUTINE NREABLD(NWRST,NBCSIM,NBESIM,NBPSIM,
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
#include <dcbxnr.h>
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
      SUBROUTINE NRWRSVEC(ICOUNT,SCVEC,SEVEC,SPVEC)
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
#include <dcbxnr.h>
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
C  /* Deck xnrred */
      SUBROUTINE XNRRED(REDGP,REDE,REDS,GPCI,GPOE,GPOP,
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
C     S times the N b-vectors are created in XNRSLI and stored in sigma-vectors
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
#include <dcbxnr.h>
#include <dgroup.h>
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),REDS(NREDM,NREDM),
     &          GPCI(*),GPOE(*),GPOP(*),
     &          IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),IBTYP(2,NREDM),
     &          BCVEC(NZCONFQ,*),BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          SCVEC(NZCONFQ,*),SEVEC(NZXOPEQ,*),SPVEC(NZXOPPQ,*),
     &          WORK(*)
C
      CALL QENTER('XNRRED')
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
CTRONDC     ... needed in XNRSLI (+ more)
CTRONDCTSMAERKE: skal det ikke vaere LZYCONF ???
CTROND      LNEED = MAX(LNEED,LZYXOPE,LZYXOPP,LZCONF)
CTROND      IF (LNEED .GT. LFREE) CALL ERRWORK('XNRRED',-LNEED,LFREE)
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
        CALL NREXTRVC(IPX,REDGP,GPCI,GPOE,GPOP,BCVEC,BEVEC,BPVEC,
     &              IBTYP,IBCVC,IBEVC,IBPVC)
        CALL DSCAL(NTSIM,DM1,REDGP(IZRED+1),1)
        IF(IPRXNR.GE.4) THEN
          WRITE(LUPRI,'(A)') '* XNRRED: Minus reduced gradient:'
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
      CALL NREXTRDE(REDE,IBTYP,IBCVC,IBEVC,IBPVC,
     &            BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
CTROND      CALL NREXTRMT(1,REDE,IBTYP,IBCVC,IBEVC,IBPVC,
CTROND     &            BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
      CALL MEMREL('EXTRE2',WORK,1,KB1,KFREE,LFREE)
C
C     Print section
C     =============
C
      IF(IPRXNR.GE.4) THEN
        WRITE(LUPRI,'(A,I5)') 
     &    'XNRRED: Reduced Hessian. Dimension:',NZRED
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
        CALL XNRTVC(BEVEC,IBEVC,NESIM,IERED,
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
        CALL NREXTRMT(-1,REDS,IBTYP,IBCVC,IBEVC,IBPVC,
     &              BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,WORK(KB1))
        CALL MEMREL('EXTRS2',WORK,1,KB1,KFREE,LFREE)
C
C       Print section
C       =============
C
        IF(IPRXNR.GE.4) THEN
          WRITE(LUPRI,'(A,I5)') 
     &      'XNRRED: Reduced metric. Dimension:',NZRED
C            CALL OUTPUT(REDS,1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
        CALL PRIRMT(REDS,NREDM,NZRED,IBTYP)
        ENDIF
      ENDIF
C
C *** End of subroutine XNRRED
C
 9999 CONTINUE
      CALL QEXIT('XNRRED')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck xnrini */
      SUBROUTINE XNRINI(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &            GPCI,GPOE,GPOP,REDGP,REDS,REDE,
     &            RCNV,EVALR,EVECR,RNORM,
     &            WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Initialize XNR-module
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
#include <dcbxnr.h>
#include <dgroup.h>
      LOGICAL   BINMEM
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),
     &          GPCI(*),GPOE(*),GPOP(*), WORK(*),
     &          REDGP(*), REDS(*), REDE(*), EVALR(*), 
     &          EVECR(NREDM,NFREQ),RNORM(*),RCNV(*)
C
      CALL QENTER('XNRINI')
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
      NRNORM = NVPAR*NFREQ
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
      CALL XNRST(BINMEM,MAXSIM,IBTYP,IBCVC,IBEVC,IBPVC,
     &                 IEVEC,IPVEC,ICVEC,GPCI,GPOE,GPOP,EVALR,
     &                 WORK,KFREE,LFREE)
      CALL QEXIT ('XNRINI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extrmt */
      SUBROUTINE NREXTRMT(IPAR,REDMAT,IBTYP,IBCVC,IBEVC,IBPVC,
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
#include <dcbxnr.h>
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
      SUBROUTINE NREXTRDE(REDMAT,IBTYP,IBCVC,IBEVC,IBPVC,
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
#include <dcbxnr.h>
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
      SUBROUTINE NREXTRVC(IPAR,REDVC,FCVEC,FEVEC,FPVEC,
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
#include <dcbxnr.h>
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
      IF(IPRXNR.GE.200) THEN
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
      IF(IPRXNR.GE.50) THEN
        WRITE(LUPRI,'(A,I3)')
     &    'EXTRVC: Reduced vector Dim.: NZRED = ',NZRED
        CALL OUTPUT(REDVC,1,NZRED,1,1,NREDM,1,1,LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xnrsol */
      SUBROUTINE XNRSOL (REDGP,REDE,REDS,EVALR,EVECR,
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
#include <dcbxnr.h>
C
      DIMENSION REDGP(NREDM),REDE(NREDM,NREDM),REDS(NREDM,NREDM),
     &          EVALR(NFREQ),EVECR(NREDM,NFREQ),IBTYP(2,*),RCNV(*),
     &          WORK(*)
C
      PARAMETER ( D0 = 0.0D0 , D1 = 1.0D0 , D2 = 2.0D0 )
      PARAMETER ( DRTEST = 1.0D-5 )
C
      CALL QENTER('XNRSOL')
C
C     Check dimension of reduced problem
C     ==================================
C
      IF (NZRED.GT.NREDM) THEN
         WRITE(LUPRI,'(//A/A,I5,/A,I5)')
     *   ' >>> ERROR IN XNRSOL >>>',
     *   ' Dimension of reduced space is  ',NZRED,
     *   ' which exceeds allowed dimension',NREDM
         CALL QUIT('XNRSOL: Too large dimension of reduced space')
      END IF
C
C     Print section
C     =============
C
      IF(IPRXNR.GE.50) THEN
        WRITE(LUPRI,'(A,I5)') 
     &    'XNRSOL: Reduced Hessian. Dimension:',NZRED
C          CALL OUTPUT(REDE,1,NZRED,1,NZRED,NREDM,NREDM,1,LUPRI)
        CALL PRIRMT(REDE,NREDM,NZRED,IBTYP)
        IF(.NOT.STATIC) THEN
          WRITE(LUPRI,'(A,I5)') 
     &      'XNRSOL: Reduced metric. Dimension:',NZRED
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
     &                 EVECR(1,IFREQ),IBTYP,IPRXNR,WORK,KFREE,LFREE)
         ENDDO
C
C        Print section
C        -------------
C
         IF(IPRXNR.GE.1) THEN
           WRITE(LUPRI,'(/A,I4)') '<<B;B>> at iteration: ',ITMIC
           CALL PRSYMB(LUPRI,'-',44,0)
           DO IFREQ = 1,NFREQ
           VLRF = DDOT(NZRED,REDGP,1,EVECR(1,IFREQ),1)
           IF(RCNV(IFREQ).GT.THCXNR) THEN
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
     &          IERR,IPRXNR,WORK,KFREE,LFREE)
      ENDIF
C
      CALL QEXIT('XNRSOL')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xnrnex */
      SUBROUTINE XNRNEX(IBTYP,IBCVC,IBEVC,IBPVC,RCNV,EVALR,EVECR,
     &                  GPCI,GPOE,GPOP,RNORM,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Purpose:
C       1) Construct residual:
C
C             LINEQ = T:  R(I) = (E[2] - W(I)*S[2])*X(I) + G
C
C       2) Test for convergence of NEXCNV eigenvectors
C          Convergence criterium: |R(I)| .LE. THCXNR * |X(I)|
C
C       3) Generate new trial vectors
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxnr.h>
      DIMENSION IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),RCNV(*),
     &          EVALR(*),EVECR(NREDM,*),GPCI(*),GPOE(*),GPOP(*),
     &          RNORM(NVPAR,NFREQ),WORK(*)
      CALL QENTER('XNRNEX')
      KFRSAV = KFREE
      NESIM = 0
      NPSIM = 0
      NCSIM = 0
C
C       Memory allocation 
C
      IFREQ = 1
      MFREQ =NFREQ
      IF(STATIC) THEN
        NSIM = MFREQ
      ELSE
        NSIM = MFREQ + MFREQ
      ENDIF
      NBTOT = NSIM*(NZXOPEQ+NZXOPPQ+NZCONFQ)
      CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
      KFIRST = KBTOT
      KBOE   = KFIRST
      KFIRST = KFIRST + NZXOPEQ*NSIM
      KBOP   = KFIRST
      KFIRST = KFIRST + NZXOPPQ*NSIM
      KBCI   = KFIRST
      KFIRST = KFIRST + NZCONFQ*NSIM
C
C     Allocate buffers
C
      KFRBUF = KFREE
      CALL MEMGET('REAL',KBEBUF,NZXOPEQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBPBUF,NZXOPPQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBCBUF,NZCONFQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXNORM,NVPAR*MFREQ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBVEC,NVPAR*MFREQ,WORK,KFREE,LFREE)
C
C     Generate new trial vectors
C     ==========================
C
      CALL XNRNE1(WORK(KBOE),WORK(KBOP),WORK(KBCI),
     &            MFREQ,EVALR(IFREQ),EVECR(1,IFREQ),
     &            GPOE,GPOP,GPCI,IBTYP,IBCVC,IBEVC,IBPVC,
     &            WORK(KBEBUF),WORK(KBPBUF),WORK(KBCBUF),RCNV,
     &            RNORM(1,IFREQ),WORK(KXNORM),WORK(KIBVEC),
     &            WORK,KFREE,LFREE)
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
        CALL MEMREL('XNRNEX.CMP',WORK,1,KFRSAV,KFREE,LFREE)
        CALL MEMGET('REAL',KBTOT,NBTOT,WORK,KFREE,LFREE)
      ELSE
        CALL MEMREL('XNRNEX.BUF',WORK,1,KFRBUF,KFREE,LFREE)
      ENDIF
      CALL QEXIT('XNRNEX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PRBVEC */
      SUBROUTINE NRPRBVEC(IUNIT,BVEC,NBSIM,NBDIM)
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
#include <dcbxnr.h>
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
      SUBROUTINE NRBVINFO(ISIM,IBTYP,NBVEC)
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
#include <dcbxnr.h>
#include <dcbibn.h>
      DIMENSION IBTYP(2,NBVEC),ISIM(*)
C
C     Initialize
C
      CALL IZERO(ISIM,NVPAR)
C
C     Count vectors
C
      DO I = 1,NBVEC
        IBV = IBTYP(1,I)
        IBH = IBTYP(2,I)
        IF    (IBV.EQ.JBENDX) THEN
          IF(IBH.EQ.IPX) THEN
            ISIM(IPEP) = ISIM(IPEP) + 1
          ELSE
            ISIM(IPEM) = ISIM(IPEM) + 1
          ENDIF
        ELSEIF(IBV.EQ.JBPNDX) THEN
          IF(IBH.EQ.IPX) THEN
            ISIM(IPPP) = ISIM(IPPP) + 1
          ELSE
            ISIM(IPPM) = ISIM(IPPM) + 1
          ENDIF
        ELSEIF(IBV.EQ.JBCNDX) THEN
          IF(IBH.EQ.IPX) THEN
            ISIM(IPCP) = ISIM(IPCP) + 1
          ELSE
            ISIM(IPCM) = ISIM(IPCM) + 1
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XNROUT */
      SUBROUTINE XNROUT(RCNV,RNORM,ISIM,EVALR,IBTYP)
C*****************************************************************************
C
C     Final output from XNRCTL
C       RCNV  - total norm of residual relative to norm of solution vector
C       RNORM - residual norm
C
C     Written by T.Saue Nov 19 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxnr.h>
C
      CHARACTER CPUTID*12,WORD*3
      DIMENSION RNORM(NVPAR,NFREQ),EVALR(NFREQ),ISIM(*),
     &          RCNV(NFREQ),IBTYP(2,NREDM)
C
      CALL HEADER('Final output from XNRCTL',-1)
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
          WRITE(LUPRI,'(/A,I4,A)')
     *      ' *** WARNING-Maximum number of microiterations,',
     *      ITMIC,', reached'
        END IF
      END IF
C
C     Give summary of calculations
C
      WRITE(LUPRI,'(/12X,A9,3X,A10,2X,6(1X,A2,A6,5X))')
     &  'Freq.(au)','Total norm',(BVTYP(I),' norm ',I=1,NVPAR)
      NCOL = 36 + NVPAR*14
      REWIND LUCYCL
      DO ITR = 1,ITMIC
        CALL PRSYMB(LUPRI,'=',NCOL,0)
        READ(LUCYCL,'(I3,2X,A9,6(2X,A3,I3))')
     &     JTR,INTTYP,(WORD,ISIM(I),I=1,NVPAR)
        READ(LUCYCL,'(6E10.2)') ((RNORM(I,J),I=1,NVPAR),J=1,NFREQ)
        READ(LUCYCL,'(A12)') CPUTID
        WRITE(LUPRI,'(A5,I3,3X,A12,A,A9,A)') 
     &        '* It.',ITR,CPUTID,' (',INTTYP,')'
        CALL PRSYMB(LUPRI,'-',NCOL,0)
        DO IFREQ = 1,NFREQ
           RTNORM = DNRM2(NVPAR,RNORM(1,IFREQ),1)
           WRITE(LUPRI,'(10X,F10.6,A3,E10.2,A3,6(E10.2,A1,I2,A1))') 
     &     EVALR(IFREQ),' : ',RTNORM,' | ',
     &     (RNORM(I,IFREQ),'(',ISIM(I),')',I=1,NVPAR)
        ENDDO
      ENDDO
      CALL PRSYMB(LUPRI,'=',NCOL,0)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck E2DIAG */
      SUBROUTINE NRE2DIAG(HEDIAG,HPDIAG,HCDIAG,WORK,KFREE,LFREE)
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
#include <dcbxnr.h>
      DIMENSION HEDIAG(*),HPDIAG(*),HCDIAG(*),WORK(*)
C
C     Orbital part
C
      IF((NZXOPE.GT.0).OR.(NZXOPP.GT.0)) THEN
        CALL NRE2ODIA(HEDIAG,HPDIAG,WORK(KZXOPE),WORK(KZXOPP),
     &              WORK,KFREE,LFREE)
      ENDIF
C
C     Configurational part
C
      IF(NZCONF.GT.0) THEN
C        CALL NRE2CDIA(HCDIAG,WORK,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck NRE2ODIA */
      SUBROUTINE NRE2ODIA(HEDIAG,HPDIAG,JXOPE,JXOPP,WORK,KFREE,LFREE)
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
#include <dcbxnr.h>
#include <dcbbas.h>
#include <dcbibn.h>
C
      DIMENSION HEDIAG(*),HPDIAG(*),JXOPE(2,*),JXOPP(2,*),WORK(*)
C
      KFRSAV = KFREE
      CALL MEMGET('REAL',KFMO,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Read from file KRMC FOCK
C
      CALL KRMC_GETFCK(WORK(KFMO),IPRXRS,WORK(KFREE),LFREE)
C
C     Electronic part
C
      IF(NZXOPE.GT.0) THEN
        IF(STERNH) THEN
          FAC = DM2*CVAL*CVAL
          CALL DCOPY(NZXOPE,FAC,0,HEDIAG,1)
        ELSE
          CALL E2ODI1(HEDIAG,WORK(KFMO),JXOPE,NZXOPE)
        ENDIF
        IF(IPRXNR.GE.5) THEN
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
        IF(IPRXNR.GE.5) THEN
          WRITE(LUPRI,'(/A)') '** E2ODIA: Diagonal Hessian (e-p):'
          WRITE(LUPRI,'(2I5,3X,1P,E18.10)')
     &      (JXOPP(1,I),JXOPP(2,I),HPDIAG(I), I = 1,NZXOPP)
        ENDIF
      ENDIF
C
C
      CALL MEMREL('XNRSEP',WORK,1,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XNRST2 */
      SUBROUTINE XNRST2(TYP,BVEC,GP,HDIAG,EVALR,IBTYP,IBVEC,IVEC,
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
#include <dcbxnr.h>
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
          CALL GETDAT(XNRLBL(1),XNRLBL(2))
          WRITE(XNRLBL(2)(1:1),'(A1)') TYP
          WRITE(XNRLBL(2)(2:2),'(A1)') '+'
          DO IP = 1,NFREQ
            XNRLBL(2)(2:2) = '+'
            XNRLBL(2)(3:3) = CHRNOS(IP/10)
            XNRLBL(2)(4:4) = CHRNOS(MOD(IP,10))
            IM   = IP + NFREQ
            FREQ = EVALR(IP)
            IF(FREQ.EQ.D0) THEN
              IF(FNDXLB(XNRLBL,LURST)) THEN
                WRITE(LUPRI,'(A)') '*** WARNING ***'
                WRITE(LUPRI,'(A,A2,A/3(3X,A8))')
     &            'Restart on ',XNRLBL(2)(1:2),'-vector:',
     &            (XNRLBL(I),I=1,3)
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
              IF(FNDXLB(XNRLBL,LURST)) THEN
                WRITE(LUPRI,'(A)') '*** WARNING ***'
                WRITE(LUPRI,'(A,A2,A/3(3X,A8))')
     &            'Restart on ',XNRLBL(2)(1:2),'-vector:',
     &            (XNRLBL(I),I=1,3)
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
              XNRLBL(2)(2:2) = '-'
              IF(FNDXLB(XNRLBL,LURST)) THEN
                WRITE(LUPRI,'(A)') '*** WARNING ***'
                WRITE(LUPRI,'(A,A2,A/3(3X,A8))')
     &            'Restart on ',XNRLBL(2)(1:2),'-vector:',
     &            (XNRLBL(I),I=1,3)
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
        CALL XNRORT(BVEC,JBTYP,NBVEC,NBPREV,
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
