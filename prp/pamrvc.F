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
C FILE    : pamrvc.F
C SOURCE  : $Source$
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssvc */
      SUBROUTINE XRSSVC(BCVEC,BEVEC,BPVEC,
     &                  SCVEC,SEVEC,SPVEC,IBTYP,
     &                  IBCVC,IBEVC,IBPVC,NPOS,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate sigma vector SVEC:
C     molecular Hessian E[2] times trial vector B(I)
C
C        SVEC(I) = E[2]*B(I)
C
C     The trial vectors come in three types:
C        C-type: contains only configurational parameters
C        E-type: contains only electron-eelectron (e-e) rotations
C        P-type: contains only electron-positron  (e-p) rotations
C
C     The sigma vectors are split up into a C, E and P -part.
C
C     Written by H.J.Aa.Jensen and Trond Saue 1996
C     Last revision July 31 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION BCVEC(*),BEVEC(*),BPVEC(*),
     &          SCVEC(*),SEVEC(*),SPVEC(*),
     &          IBCVC(*),IBEVC(*),IBPVC(*),
     &          IBTYP(2,*),NPOS(*),WORK(*)
C
C Used from common blocks:
C   DCBXRS : LZYVAR,NCSIM,NESIM,NPSIM,NOSIM
C
#include <dcbxrs.h>
C
C
      CALL QENTER('XRSSVC')
C
C     Initialize sigma vectors
C     ========================
C
      IF(NZCONF.GT.0) CALL DZERO(SCVEC,NTSIM*NZCONFQ)
      IF(NZXOPE.GT.0) CALL DZERO(SEVEC,NTSIM*NZXOPEQ)
      IF(NZXOPP.GT.0) CALL DZERO(SPVEC,NTSIM*NZXOPPQ)
C
C     Construct orbital part of sigma vectors
C     ====================================================
C
chj   IF (.NOT.RSPCI) THEN
         CALL XRSSEP(BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,
     &               IBEVC,IBPVC,IBTYP,NPOS,
     &               WORK,KFREE,LFREE)
chj   END IF
C
C     Construct configurational part of sigma vectors
C     for the configurational trial vectors
C     ===============================================
C
      IF (NCSIM.GT.0) THEN
         CALL XRSSCI(BCVEC,SCVEC,SEVEC,SPVEC,WORK,KFREE,LFREE)
      END IF
C
      CALL QEXIT('XRSSVC')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrscli */
      SUBROUTINE XRSSCI(BCVEC,SCVEC,SEVEC,SPVEC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate cc block of sigmavector
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision July 30 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION BCVEC(*),SCVEC(*),SEVEC(*),SPVEC(*),WORK(*)
C
C Used from common blocks:
C     XCBPRI: LUPRI
C     WRKXRS: NCSIM
C
#include <dcbxrs.h>
      CALL QENTER('XRSSCI')
      WRITE(LUPRI,'(A)') 
     &     '*** Sorry, XRSSCI called, but is not implemented yet!'
      CALL QUIT('XRSCLI called, but is not implemented yet!')
      CALL QEXIT ('XRSSCI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrssep */
      SUBROUTINE XRSSEP(BCVEC,BEVEC,BPVEC,SCVEC,SEVEC,SPVEC,
     &                  IBEVC,IBPVC,IBTYP,NPOS,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Construct orbital part of sigma vectors
C
C     SA_(pq) = E[2](o*o)_(pq,rs)Z_(rs) + E[2](o*o*)_(pq,rs)Y^*_(rs)
C             = -<0|[Qp,[Rs,H]]|0>Z_(rs) + <0|[Qp,[Sr,H]]|0>Y^*_(rs)
C             = -<0|[Qp,[Rs,H]]|0>Z_(rs) - <0|[Qp,[Sr,H]]|0>{-Y+}_(sr)
C             = -<0|[Qp,[Rs,H]]|0>WA_(r,s)
C
C     where the matrix WA has the structure:
C
C           P  I   A   S             P  I  A  S
C          --------------         ---------------
C        P | 0  0   0   0       P |  0  Z  Z  0
C        I | 0  0  -Y+ -Y+      I | -Y+ 0  0  0
C        A | 0  Z   X  -Y+      A | -Y+ 0  0  0
C        S | 0  Z   Z   0       S |  0  0  0  0
C
C           e-e rotations           e-p rotations
C
C     with the convention:
C        P - positronic spinors
C        I - inactive spinors
C        A - active spinors
C        S - secondary spinors
C
C     (X has -Y+ in upper diagonal and Z in lower diagonal)
C
C     Similarly
C
C     SB_(pq) = E[2](o*o)_(pq,rs)Y_(rs) + E[2](o*o*)_(pq,rs)Z^*_(rs)
C             = -<0|[Qp,[Rs,H]]|0>WB_(r,s)
C
C     where WB is identical to WA with Y and Z interchanged.
C
C-----------------------------------------------------------------------------
C
C     The key to efficient evaluation of SVEC is to reformulate it
C     as a gradient element:
C
C     SA_(pq) = -<0|[Qp,HA]]|0>
C
C     where HA is a one-index transformed Hamiltonian
C
C     HA = [Tu,H] WA_(t,u) = ha_(pq)Pq + ga_(pq,rs)PRsq
C
C     with
C
C     ha_(pq)    = WA_(pt)h_(tq)    - h_(pt)WA_(tq)
C
C     ga_(pq,rs) = WA_(pt)g_(tq,rs) - g_(pt,rs)WA_(tq)
C                + WA_(rt)g_(pq,ts) - g_(pq,rt)WA_(ts)
C
C     In RPA we need to evaluate
C
C     SA_(si) = -FA_(si) - GA_(si)
C     SB_(si) = -FB_(si) - GB_(si)
C
C     with
C      
C     FA_(si) = WA_(st)h_(ti) - h_(st)WA_(ti)
C
C     GA_(si) = [g_(si,xy)-g_(sy,xi)] DA_(yx)
C
C     where x and y are AO-indices and
C
C     DA_(yx) = DA1_(yx) - DA2_(yx)
C
C     DA1_(yx) = C_(yj)(-Y+)_(js)(C+)_(sx)
C     DA2_(yx) = C_(ys)( Z )_(sj)(C+)_(jx)
C
C-----------------------------------------------------------------------------
C     
C     For efficient evaluation we shall exploit the relations
C
C     Time - symmetric operators:      DA1 = -[DA2^+]      
C     Time - anti symmetric operators: DA1 = -[-i(DA2^+)i] (i-transform)
C
C     
C-----------------------------------------------------------------------------
C     
C     The algorithm proceeds in the following steps:
C
C     1. Generate standard and modified two-electron Fock matrices.
C     2. Accumulate GA/GB into SVEC
C     3. Accumulate FA/FB into SVEC
C
C     Step 1 is the most memory and CPU-intensive step, but after
C     step 2 most of the memory can be released.
C
C     The b-vectors are scatterd into W matrices twice to 
C     minimize memory.
C
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C     Last revision Aug 2 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DM1 = -1.0D0)
      DIMENSION BCVEC(*),BEVEC(*),BPVEC(*),
     &          SCVEC(*),SEVEC(*),SPVEC(*),
     &          IBEVC(*),IBPVC(*),IBTYP(2,*),NPOS(*),
     &          WORK(*)
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      CALL QENTER('XRSSEP')
C
      KFRSAV = KFREE
C
C     1. Generate two-electron Fock matrices
C     ======================================
C
      CALL MEMGET('REAL',KFMO,NOSIM*N2BBASXQ*NFOCKMAT,WORK,KFREE,LFREE)
C
C     1. Generate G-matrix
C     ====================
C
      CALL DZERO(WORK(KFMO),NOSIM*N2BBASXQ*NFOCKMAT)
      IF(LSVCFG(1)) THEN
        CALL GMOLI(WORK(KFMO),BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &           NPOS,WORK,KFREE,LFREE)
      ELSE
c       CALL DZERO(WORK(KFMO),NOSIM*N2BBASXQ*NFOCKMAT)
      ENDIF
C
C     2) Add one-index transformed Fock-matrices
C     ==========================================
C
      IF (LSVCFG(2)) CALL FMOLI(WORK(KFMO),BEVEC,BPVEC,
     &                          IBEVC,IBPVC,IBTYP,WORK,KFREE,LFREE)
C
C     3) Gather contributions to sigma-vectors
C     ===================================================
C
      IF(NZXOPE+NZXOPP.GT.0) THEN
        CALL SMOLI(SEVEC,SPVEC,WORK(KFMO),
     &             WORK(KZXOPE),NZXOPE,WORK(KZXOPP),NZXOPP,
     &             IBEVC,IBPVC,IBTYP,WORK(KFREE),LFREE)
        IF(IPRXRS.GE.5) THEN
          WRITE(LUPRI,'(/I5,A)') NTSIM,
     &          ' orbital (e-e) sigma vectors'
          CALL PRBVEC(LUPRI,SEVEC,NTSIM,NZXOPE)
          WRITE(LUPRI,'(/I5,A)') NTSIM,
     &          ' orbital (e-p) sigma vectors'
          CALL PRBVEC(LUPRI,SPVEC,NTSIM,NZXOPP)
        ENDIF
      ENDIF
C
      CALL MEMREL('XRSSEP',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT ('XRSSEP')
      RETURN
      END
C
C     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gmoli */
      SUBROUTINE GMOLI(FMO,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                 NPOS,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate G-matrix: 2-electron Fock matrix 
C                        with modified density matrices
C
C     Written by T.Saue Sep 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <aovec.h>
#include <mxorb.h>
C
#include <dcbxrs.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
C
      DIMENSION FMO(N2ORBXQ,*),BEVEC(*),BPVEC(*),
     &          IBEVC(*),IBPVC(*),IBTYP(2,*),
     &          NPOS(*),WORK(*)
C
      CALL QENTER('GMOLI')
      KFRSAV = KFREE
      NDMAT  = NOSIM*NZ
      CALL MEMGET('REAL',KCMO ,NCMOTQ       ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAO ,NDMAT*N2BBASX*NFOCKMAT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KSYM ,NOSIM*NFOCKMAT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHRM ,NOSIM*NFOCKMAT,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KFCK ,NOSIM*NFOCKMAT,WORK,KFREE,LFREE)
C
      CALL GMOLI1(FMO,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &            WORK(KCMO),WORK(KFAO),WORK(KSYM),
     &            WORK(KHRM),WORK(KFCK),
     &            NPOS,WORK,KFREE,LFREE)
C
      CALL MEMREL('GMOLI',WORK,1,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('GMOLI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gmoli1 */
      SUBROUTINE GMOLI1(FMO,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                  CMO,FAO,ISYMOP,IHRMOP,IFCKOP,
     &                  NPOS,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate G-matrix: 2-electron Fock matrix 
C                        with modified density matrices
C
C     Written by T.Saue Sep 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <cbihr2.h>
#include <dcbxrs.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbham.h>
      LOGICAL TOBE
      DIMENSION FMO(N2ORBXQ,*),FAO(N2BBASXQ,*),
     &          BEVEC(*),BPVEC(*),IBEVC(*),IBPVC(*),IBTYP(2,*),
     &          CMO(*),ISYMOP(*),IHRMOP(*),IFCKOP(*),
     &          NPOS(*),WORK(*)
C
C     a) Get coefficients
C     ===================
C
      LUBUF = 1
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','GMOLI1')
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUM,DUM,TOTERG,4)
      IF (IPRXRS .GE. 10) THEN
         CALL HEADER('MO Coefficients',-1)
         DO I = 1, NFSYM
            WRITE(LUPRI,'(/1X,A,I1/)') '*** Symmetry ',I
            CALL PRQMAT(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),
     &         NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         END DO
      END IF
      CLOSE(LUBUF,STATUS='KEEP')
C
C     b) Set arrays of IHRMOP,IFCKOP,ISYMOP
C     =====================================
C
      CALL SETGMT(ISYMOP,IHRMOP,IFCKOP,
     &            IBEVC,IBPVC,IBTYP,IPRXRS)
C
C     c) Accumulate density matrices
C        (temporary storage in FMO)
C     =================================
C
      CALL DMOLI(FMO,BEVEC,BPVEC,CMO,
     &           IHRMOP,WORK,KFREE,LFREE)
C
C     b) Accumulate two-electron Fock matrices
C     ----------------------------------------
C
      CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,
     &            FAO,FMO,NOSIM*NFOCKMAT,NPOS,INTFLG,
     &            IPRTWO,WORK(KFREE),LFREE)
C
C     If solvent; get contribution
C     ============================
C
      IF(SOLVEN) CALL SOLFCK(FAO,FMO,NOSIM*NFOCKMAT,ESOLVE,ESOLVN,
     &                       WORK(KFREE),LFREE,IPRSOL)
C
      IM = 0
C
C     c) Transform to MO-basis; e-e rotations
C     ---------------------------------------
C
      DO I = 1,NESIM
C
C         Inactive part...
C
          IF (NISHT .GT. 0) THEN
             ITYP  = 0*1 + 0*2
             IM = IM + 1
             CALL GMOLI3(ITYP,FMO(1,IM),CMO,FAO(1,IM),
     &            WORK,KFREE,LFREE)
C            Print section
             IF(IPRXRS.GE.8) THEN
               IREP = JSYMOP - 1
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: inactive (e-e) G-matrix (AO-basis) no.:',I
               CALL PRQMAT(FAO(1,IM),NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: inactive (e-e) G-matrix (MO-basis) no.:',I
               CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
             ENDIF
          END IF
C
C         Active part...
C
          IF (NASHT .GT. 0) THEN
             ITYP  = 0*1 + 1*2
             IM = IM + 1
             CALL GMOLI3(ITYP,FMO(1,IM),CMO,FAO(1,IM),
     &            WORK,KFREE,LFREE)
C            Print section
             IF(IPRXRS.GE.8) THEN
               IREP = JSYMOP - 1
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: active (e-e) G-matrix (AO-basis) no.:',I
               CALL PRQMAT(FAO(1,IM),NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: active (e-e) G-matrix (MO-basis) no.:',I
               CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
             ENDIF
          ENDIF
      ENDDO
C
C     d) Transform to MO-basis; e-p rotations
C     ---------------------------------------
C
      DO I = 1,NPSIM
C
        IF (NISHT .GT. 0) THEN
           IM = IM + 1
           ITYP  = 1*1 + 0*2
           CALL GMOLI3(ITYP,FMO(1,IM),CMO,FAO(1,IM),
     &            WORK,KFREE,LFREE)
C            Print section
             IF(IPRXRS.GE.8) THEN
               IREP = JSYMOP - 1
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: (e-p) G-matrix (AO-basis) no.:',I
               CALL PRQMAT(FAO(1,IM),NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: (e-p) G-matrix (MO-basis) no.:',I
               CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
             ENDIF
        ENDIF
        IF (NASHT .GT. 0) THEN
           IM = IM + 1
           ITYP  = 1*1 + 1*2
           CALL GMOLI3(ITYP,FMO(1,IM),CMO,FAO(1,IM),
     &            WORK,KFREE,LFREE)
C            Print section
             IF(IPRXRS.GE.8) THEN
               IREP = JSYMOP - 1
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: (e-p) G-matrix (AO-basis) no.:',I
               CALL PRQMAT(FAO(1,IM),NTBAS(0),NTBAS(0),
     &                     NTBAS(0),NTBAS(0),NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
               WRITE(LUPRI,'(A,I5)') 
     &           '* GMOLI2: (e-p) G-matrix (MO-basis) no.:',I
               CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &                     IPQTOQ(1,IREP),LUPRI)
             ENDIF
         ENDIF
      ENDDO
      RETURN
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'GMOLI1: Coefficient file not found !'
      CALL QUIT('GMOLI1: Coefficients not found !')
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmoli */
      SUBROUTINE DMOLI(DMAT,BEVEC,BPVEC,CMO,IHRMOP,
     &                 WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Accumulate density matrices for use in generation
C     of orbital part of sigma vectors
C
C     Written by T.Saue August 2 1996
C     Last revision Aug 2 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION DMAT(*),BEVEC(*),BPVEC(*),CMO(*),
     &          IHRMOP(*),WORK(*)
C
#include <dcbxrs.h>
#include <dcborb.h>
#include <dcbbas.h>
C
      CALL QENTER('DMOLI')
      KFRSAV = KFREE
C
      CALL MEMGET('REAL',KBMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL DMOLI2(DMAT,BEVEC,BPVEC,CMO,WORK(KBMAT),
     &            IHRMOP,WORK,KFREE,LFREE)
C
      CALL MEMREL('DMOLI',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('DMOLI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmoli2 */
      SUBROUTINE DMOLI2(DMAT,BEVEC,BPVEC,CMO,BMAT,
     &                  IHRMOP,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Accumulate density matrices for use in generation
C     of orbital part of sigma vectors
C
C     Written by T.Saue August 2 1996
C     Last revision Aug 2 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
C
#include <dcbxrs.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
C
      DIMENSION DMAT(N2BBASXQ,*),CMO(*),BEVEC(*),BPVEC(*),
     &          BMAT(*),IHRMOP(*),WORK(*)
C
      IM = 0
      IREP = JSYMOP - 1
C
C     Generate modified density matrices ; e-e rotations
C     ===================================================
C
      IBOFF = 1
      DO I = 1,NESIM
C
C        Inactive...
C
         IF (NISHT .GT. 0) THEN
         ITYP  = 0*1+0*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BEVEC(IBOFF),WORK(KZXOPE),NZXOPE,BMAT,NORBT,NZ)
C
C        Take symmetric/antisymmetric combinations to form BMAT
C        ======================================================
C
C        Dette skal vaere modsat af DMAT (i DMOLI3)
C
         IH = IHRMOP(IM)
         IZOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
           ENDIF
           IZOFF = IZOFF + N2ORBX
         ENDDO
C
C        Construct modified density matrix
         CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT,
     &                 WORK,KFREE,LFREE)
C        Print section
         IF(IPRXRS.GE.8) THEN
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Scattered inactive (e-e) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Modified inactive (e-e) matrix no.:',I
            CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         END IF
         END IF
C
         IF (NASHT .GT. 0) THEN
C
C        Active...
C
         ITYP  = 0*1+1*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BEVEC(IBOFF),WORK(KZXOPE),NZXOPE,BMAT,NORBT,NZ)
C
C        Take symmetric/antisymmetric combinations to form BMAT
C        ======================================================
C
C        Dette skal vaere modsat af DMAT (i DMOLI3)
C
         IH = IHRMOP(IM)
         IZOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
           ENDIF
           IZOFF = IZOFF + N2ORBX
         ENDDO
C        Construct modified density matrix
         CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT,
     &                 WORK,KFREE,LFREE)
C        Print section
         IF(IPRXRS.GE.8) THEN
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Scattered active (e-e) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Modified active (e-e) matrix no.:',I
            CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         END IF
         END IF
         IBOFF = IBOFF + NZXOPEQ
C
      ENDDO
C
C     Generate modified density matrices ; e-p rotations
C     ===================================================
C
      IBOFF = 1
      DO I = 1,NPSIM
C
         IF (NISHT .GT. 0) THEN
C
C        Inactive...
C
         ITYP = 1*1 + 0*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BPVEC(IBOFF),WORK(KZXOPP),NZXOPP,BMAT,NORBT,NZ)
C
C        Take symmetric/antisymmetric combinations to form BMAT
C        ======================================================
C
C        Dette skal vaere modsat af DMAT (i DMOLI3)
C
         IH = IHRMOP(IM)
         IZOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
           ENDIF
           IZOFF = IZOFF + N2ORBX
         ENDDO
C
C        Construct modified density matrix
         CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT,
     &               WORK,KFREE,LFREE)
C        Print section
         IF(IPRXRS.GE.8) THEN
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Scattered inactive (e-p) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Modified active (e-p) matrix no.:',I
            CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         END IF
         END IF
C
         IF (NASHT .GT. 0) THEN
C
C        Active...
C
         ITYP = 1*1 + 1*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BPVEC(IBOFF),WORK(KZXOPP),NZXOPP,BMAT,NORBT,NZ)
C
C        Take symmetric/antisymmetric combinations to form BMAT
C        ======================================================
C
C        Dette skal vaere modsat af DMAT (i DMOLI3)
C
         IH = IHRMOP(IM)
         IZOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IZOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IZOFF))
           ENDIF
           IZOFF = IZOFF + N2ORBX
         ENDDO
C
C        Construct modified density matrix
         CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT,
     &               WORK,KFREE,LFREE)
C        Print section
         IF(IPRXRS.GE.8) THEN
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Scattered inactive (e-p) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(A,I5)') 
     &         '* DMOLI: Modified active (e-p) matrix no.:',I
            CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         END IF
         END IF
         IBOFF = IBOFF + NZXOPPQ
C
C
      ENDDO
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dmolie */
      SUBROUTINE DMOLI3(ITYP,IH,DMAT,CMO,BMAT,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate modified density matrix for use in generation of
C     sigma vectors
C
C     Bit 1 of ITYP = 0   e-e  rotations in BMAT
C                     1   e-p rotations in BMAT
C     Bit 2 of ITYP = 0   inactive density matrix
C                     1   active density matrix
C
C     We define: DB_{ab} = c_{as}B_{si}c_{bi}^*
C                a,b - AO-indices
C                s   - all indeces (electronic or positronic)
C                i   - inactive index
C
C     IH   = 1        Hermitian density matrix:  DMAT = -[DB + DB^+]
C                       based on anti-Hermitian BMAT
C     IH   = 2   anti-Hermitian density matrix:  DMAT = -[DB - DB^+]
C                       based on anti-Hermitian BMAT
C
C     For time-symmetric operators DB^+ refers to the Hermite conjugate
C     For time-antisymmetric operators DB^+ refers to the Hermite conjugate
C     of the i-transform
C
C     Written by T.Saue Sep 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
      PARAMETER(D0 = 0.0D0, DP5 = 0.5D00)
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),BMAT(NORBT,NORBT,NZ),WORK(*)
      LOGICAL   LBIT
C
      EXTERNAL BDQALGEBRA
C
C     Initialization
C     ==============
C
      IREP = JSYMOP-1
      CALL DZERO(DMAT,N2BBASXQ)
C
C     Print section
C
      IF(IPRXRS.GE.10) THEN
        WRITE(LUPRI,'(A)') 'DMOLI3: Input B-matrix:'
        CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
      ENDIF
C
C     Generate DB 
C     ============
C
      DO 10 I1 = 1,NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for first index
C
        IF (.NOT. LBIT(ITYP,1)) THEN
C
C         electronic indeces...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
          JS = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ELSE
C
C         or positronic indeces...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ENDIF
C
C       Parameters for inactive index
C
         IF (.NOT. LBIT(ITYP,2)) THEN
C
C           inactive density matrix...
C
            NI  = NISH(I2) - NFRO(I2)
            II  = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
            JI  = ICMOQ(I2) + (NPSH(I2)+ NFRO(I2))*NFBAS(I2,0) + 1
C
C          Generate modified density matrix as transformation MO to AO
C
           IF(NI.EQ.0.OR.NS.EQ.0) GOTO 10
            CALL QTRANS('MOAO','S',D0,
     &                  NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                  DMAT(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,IREP),
     &                  BMAT(IS,II,1),NORBT,NORBT,
     &                  NZ,IPQTOQ(1,IREP),
     &                  CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                  CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                  WORK(KFREE),LFREE,IPRXRS)
C
         ELSE
C
C           or active density matrix...
C
            NI  = NASH(I2)
            II  = IORB(I2) + NPSH(I2) + NISH(I2) + 1
            JI  = ICMOQ(I2) + 
     &               (NPSH(I2)+ NISH(I2))*NFBAS(I2,0) + 1
C
C
C          Generate modified density matrix as transformation MO to AO
C
           IF(NI.EQ.0.OR.NS.EQ.0) GOTO 10
            CALL QTRANS('MOAO','S',D0,
     &                  NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                  DMAT(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),
     &                  NZ,IPQTOQ(1,IREP),
     &                  BMAT(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &                  CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,IREP),
     &                  CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,IREP),
     &                  WORK(KFREE),LFREE,IPRXRS)
C
C
         END IF
 10   CONTINUE
      IF(IPRXRS.GE.5) THEN
        CALL HEADER('DMOLI3: Non-symmetrized density matrix',-1)
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
C
      IF (LBIT(ITYP,2)) THEN
C
C        For active density matrices scale the density matrix by
C        1/2 because all Fock matrices are multiplied with 2 in
C        TWOFCK.
C
         CALL DSCAL(N2BBASXQ,DP5,DMAT,1)
C
      END IF
C
C
C     Take symmetric/antisymmetric combinations to form DMAT
C     ======================================================
C
      DO IZ = 1,NZ
        IQ = IPQTOQ(IZ,IREP)
        IM = IHQMAT(IQ,IH)
        IF(IM.EQ.1) THEN
          CALL FULMAT('S',NTBAS(0),NTBAS(0),DMAT(1,IZ))
        ELSEIF(IM.EQ.2) THEN
          CALL FULMAT('A',NTBAS(0),NTBAS(0),DMAT(1,IZ))
        ENDIF
      ENDDO
C
C     Print section
C
      IF(IPRXRS.GE.5) THEN
        CALL HEADER('DMOLI3: Modified density matrix',-1)
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fmoli */
      SUBROUTINE FMOLI(FMO,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                 WORK,KFREE,LFREE)
C*****************************************************************************
C
C     FMO contains G-matrix; add one-index transformed Fock-matrix
C
C     Written by T.Saue Sep 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbxrs.h>
      DIMENSION FMO(*),BEVEC(*),BPVEC(*),IBTYP(*),
     &          IBEVC(*),IBPVC(*),WORK(*)
C
      CALL QENTER('FMOLI')
      KFRSAV = KFREE
      CALL MEMGET('REAL',KBMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFBUF,N2BBASXQ*NFOCKMAT,WORK,KFREE,LFREE)
C
C     Get Fock matrices in MO-basis (FCmo and FVmo)
C     =============================================
C
      IF (NASHT .GT. 0) THEN
         CALL GETMAT('FCMO+FVMO',WORK(KFBUF),IPRXRS,WORK(KFREE),LFREE)
      ELSE
         CALL GETMAT('FCMO',WORK(KFBUF),IPRXRS,WORK(KFREE),LFREE)
      END IF
c     CALL GETFCK(WORK(KFBUF),IPRXRS,WORK,KFREE,LFREE)
C
C     Do one-index transformation
C
      CALL FMOLI2(FMO,WORK(KFBUF),WORK(KBMAT),BEVEC,BPVEC,
     &            IBEVC,IBPVC,IBTYP,WORK,KFREE,LFREE)
C
      CALL MEMREL('FMOLI',WORK,1,KFRSAV,KFREE,LFREE)      
      CALL QEXIT('FMOLI')      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fmoli2 */
      SUBROUTINE FMOLI2(FMO,FBUF,BMAT,BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     FMO contains G-matrix; add one-index transformed Fock-matrix
C     FBUF contains original Fock-matrix
C
C     Written by T.Saue Sep 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION FMO(N2ORBXQ,*),FBUF(*),BMAT(*),BEVEC(*),
     &          BPVEC(*),IBEVC(*),IBPVC(*),IBTYP(2,*),WORK(*)
C
      IM = 0
      IVEC = 0
      IREP = JSYMOP-1
C
cjth  IF (IPRXRS .GE. 10) THEN
cjth     KFRSAV = KFREE
cjth     CALL MEMGET('REAL',KTEMP,N2ORBXQ,WORK,KFREE,LFREE)
cjth  END IF
C
C     Add one-index tranformed Fock-matrices; e-e rotations
C     =====================================================
C
      IBOFF = 1
      DO I = 1,NESIM
C
         IVEC = IVEC + 1
         IH   = IBEVC(NERED+IVEC)
         IH = IBTYP(2,IH)
         IF (NISHT .GT. 0) THEN
C
C        Inactive...
C
         ITYP  = 0*1 + 0*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BEVEC(IBOFF),WORK(KZXOPE),NZXOPE,BMAT,NORBT,NZ)
C
C        Construct full BMAT
C        ===================
C
         IOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
           ENDIF
           IOFF = IOFF + N2ORBX
         ENDDO
C
C        Print section
         IF(IPRXRS.GE.6) THEN
           WRITE(LUPRI,'(A,I5)') 
     &       '* FMOLI: Scattered inactive (e-e) trial vector no.:',I
           CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                 IPQTOQ(1,IREP),LUPRI)
         ENDIF
C
cjth     IF (IPRXRS .GE. 10) THEN
cjth        CALL DCOPY(N2ORBXQ,FMO(1,IM),1,WORK(KTEMP),1)
cjth     END IF
C        Add one-index transformed Fock-matrix
         CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF,
     &                 WORK,KFREE,LFREE)
cjth     IF (IPRXRS .GE. 10) THEN
cjth        CALL DAXPY(N2ORBXQ,-1.0D00,FMO(1,IM),1,WORK(KTEMP),1)
cjth        CALL DSCAL(N2ORBXQ,-1.0D00,WORK(KTEMP),1)
cjth        WRITE(LUPRI,'(A,I5)') 
cjth &        '* FMOLI: Inactive 1-index transformed Fock matrix:',I
cjth        CALL PRQMAT(WORK(KTEMP),NORBT,NORBT,NORBT,NORBT,NZ,
cjth &                  IPQTOQ(1,IREP),LUPRI)
cjth     END IF
         END IF
C
         IF (NASHT .GT. 0) THEN
C
C        active...
C
         ITYP  = 0*1 + 1*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BEVEC(IBOFF),WORK(KZXOPE),NZXOPE,BMAT,NORBT,NZ)
C
C        Construct full BMAT
C        ===================
C
         IOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
           ENDIF
           IOFF = IOFF + N2ORBX
         ENDDO
C        Print section
         IF(IPRXRS.GE.6) THEN
           WRITE(LUPRI,'(A,I5)') 
     &       '* FMOLI: Scattered active (e-e) trial vector no.:',I
           CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                 IPQTOQ(1,IREP),LUPRI)
         ENDIF
C        Add one-index transformed Fock-matrix
         CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF(1+N2ORBXQ),
     &                 WORK,KFREE,LFREE)
C
         END IF
         IBOFF = IBOFF + NZXOPEQ
      ENDDO
C
C     Add one-index tranformed Fock-matrices; e-p rotations
C     =====================================================
C
      IBOFF = 1
      IVEC = 0
      DO I = 1,NPSIM
C
         IVEC = IVEC + 1
         IH   = IBPVC(NPRED+IVEC)
         IH = IBTYP(2,IH)
         IF (NISHT .GT. 0) THEN
C
C        Inactive...
C
         ITYP  = 1*1 + 0*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BPVEC(IBOFF),WORK(KZXOPP),NZXOPP,BMAT,NORBT,NZ)
C
C        Construct full BMAT
C        ===================
C
         IOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
           ENDIF
           IOFF = IOFF + N2ORBX
         ENDDO
C
C        Print section
         IF(IPRXRS.GE.6) THEN
            WRITE(LUPRI,'(A,I5)') 
     &         '* FMOLI: Scattered inactive (e-p) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         ENDIF
C        Add one-index transformed Fock-matrix
         CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF,
     &               WORK,KFREE,LFREE)
         END IF
C
         IF (NASHT .GT. 0) THEN
C
C        Active...
C
         ITYP  = 1*1 + 1*2
         IM = IM + 1
C        Scatter BVEC into BMAT
         CALL DZERO(BMAT,N2ORBXQ)
         CALL XOPSCT(BPVEC(IBOFF),WORK(KZXOPP),NZXOPP,BMAT,NORBT,NZ)
C
C        Construct full BMAT
C        ===================
C
         IOFF = 1
         DO IZ = 1,NZ
           IQ = IPQTOQ(IZ,IREP)
           IS = IHQMAT(IQ,IH)
           IF(IS.EQ.1) THEN
             CALL FULMAT('A',NORBT,NORBT,BMAT(IOFF))
           ELSEIF(IS.EQ.2) THEN
             CALL FULMAT('S',NORBT,NORBT,BMAT(IOFF))
           ENDIF
           IOFF = IOFF + N2ORBX
         ENDDO
C
C        Print section
         IF(IPRXRS.GE.6) THEN
            WRITE(LUPRI,'(A,I5)') 
     &         '* FMOLI: Scattered active (e-p) trial vector no.:',I
            CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
         ENDIF
C        Add one-index transformed Fock-matrix
         CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF(1+N2ORBXQ),
     &               WORK,KFREE,LFREE)
C
         END IF
         IBOFF = IBOFF + NZXOPPQ
      ENDDO
cjth  IF (IPRXRS .GE. 10) THEN
cjth     CALL MEMREL('GMOLI',WORK,1,KFRSAV,KFREE,LFREE)
cjth  END IF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Fmoli3 */
      SUBROUTINE FMOLI3(ITYP,FMO,BMAT,FBUF,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate one-index transformed Fock-matrices
C     for use in generation of sigma vectors
C
C
C     Bit 1 of ITYP = 0   e-e  rotations in BMAT
C                     1   e-p rotations in BMAT
C     Bit 2 of ITYP = 0   inactive density matrix
C                     1   active density matrix
C
C     FMO_{si} = BMAT_{sj}FBUF_{ji} - FBUF^t_{sj}Z_{ji}
C
C     For time-antisymmetric operators FBUF^t is i-transform
C
C     Written by T.Saue Sep 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbqal.h>
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0,DM1 = -1.0D0)
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      LOGICAL EPMIX,LBIT
      DIMENSION FMO(NORBT,NORBT,NZ),FBUF(NORBT,NORBT,NZ),
     &          BMAT(NORBT,NORBT,NZ),WORK(*)
C
      EXTERNAL BDQALGEBRA
C
C     Initialization
C     ==============
C
      IF(.NOT. LBIT(ITYP,1))     THEN
        EPMIX = NZXOPP.GT.0
      ELSE
        EPMIX = NZXOPE.GT.0
      ENDIF
      IREP = JSYMOP-1
C
C     Generate one-index transformed Fock matrix
C     ==========================================
C
      DO 10 I1 = 1,NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for first index
C
        NI = NOCC(I2) - NFRO(I2)
        II = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
C
C       Parameters for second index
C
        IF    (.NOT. LBIT(ITYP,1)) THEN
C
C         Electronic...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
C
        ELSE
C
C         Positronic...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
C
        ENDIF
C
        IF (NI.NE.0 .AND. NS.NE.0) THEN
C
C         First part: FMO_{si} = -FMO_{si} + BMAT_{sj}FBUF_{ji}
C
          CALL QGEMM(NS,NI,NI,D1,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(II,II,1),NORBT,NORBT,NZ,
     &          DM1,IPQTOQ(1,IREP),FMO(IS,II,1),NORBT,NORBT,NZ)
C
C
C         Second part: FMO_{si} = FMO_{si} - FBUF^t_{st}BMAT_{ti}
C
          CALL QGEMM(NS,NI,NS,DM1,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(IS,IS,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          D1,IPQTOQ(1,IREP),FMO(IS,II,1),NORBT,NORBT,NZ)
          IF(EPMIX) THEN
C
C           Parameters for secondary index(note the switch...)
C
            IF    (.NOT. LBIT(ITYP,1)) THEN
              NT = NPSH(I1)
              IT = IORB(I1) + 1
            ELSE
              NT = NESH(I1)
              IT = IORB(I1) + NPSH(I1) + 1
            ENDIF
            IF(NT.NE.0) CALL QGEMM(NT,NI,NS,DM1,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(IT,IS,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          DM1,IPQTOQ(1,IREP),FMO(IT,II,1),NORBT,NORBT,NZ)
          ENDIF
        ENDIF
C
 10   CONTINUE
C
C     Print section
C
      IF(IPRXRS.GE.5) THEN
        CALL HEADER('FMOLI3: F-G-matrix',-1)
        CALL PRQMAT(FMO,NORBT,NORBT,NORBT,
     &              NORBT,NZ,IPQTOQ(1,IREP),LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Gmoli3 */
      SUBROUTINE GMOLI3(ITYP,FMO,CMO,FAO,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Transform the G-matrix to MO-basis, but only the blocks needed
C
C     Bit 1 of ITYP = 0   e-e  rotations in BMAT
C                     1   e-p rotations in BMAT
C     Bit 2 of ITYP = 0   inactive density matrix
C                     1   active density matrix
C
C     Written by T.Saue Sep 18 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0)
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dcbxrs.h>
#include <dgroup.h>
      LOGICAL EPMIX,LBIT
      DIMENSION FMO(NORBT,NORBT,NZ),CMO(*),
     &          FAO(N2BBASXQ),WORK(*)
C
C
C     Transform only blocks of right symmetry combination
C     ===================================================
C
      IF(.NOT. LBIT(ITYP,1))     THEN
        EPMIX = NZXOPP.GT.0
      ELSE
        EPMIX = NZXOPE.GT.0
      ENDIF
      IREP  = JSYMOP-1
      ITIM  = JTIMOP
      DO 10 I1 = 1,NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for first index
C
        IF (.NOT. LBIT(ITYP,1)) THEN
C
C         electronic indeces...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
          JS = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ELSE
C
C         or positronic indeces...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ENDIF
C
C       Parameters for second index.
C
         NI = NOCC(I2)
         II = IORB(I2) + NPSH(I2) + 1
         JI = ICMOQ(I2) + NPSH(I2) * NFBAS(I2,0) + 1
C
        IF(NS.NE.0.AND.NI.NE.0) THEN
          CALL QTRANS('AOMO','S',D0,
     &                NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),
     &                FMO(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &                CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRXRS)
        ENDIF
C
C       Transform other block as well, 
C       if both electronic and postitronic rotations exist
C     
        IF(EPMIX) THEN
C
C         Parameters for secondary index(note the switch...)
C
C
C       Parameters for first index
C
        IF (.NOT. LBIT(ITYP,1)) THEN
C
C         or positronic indeces...
C
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ELSE
C
C         electronic indeces...
C
          NS = NESH(I1)
          IS = IORB(I1) + NPSH(I1) + 1
          JS = ICMOQ(I1) + NPSH(I1)*NFBAS(I1,0) + 1
        ENDIF
          IF(NS.NE.0.AND.NI.NE.0) THEN
            CALL QTRANS('AOMO','S',D0,
     &                NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),
     &                FMO(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &                CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRXRS)
          ENDIF
        ENDIF
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Smoli */
      SUBROUTINE SMOLI(SEVEC,SPVEC,FMO,
     &                 JXOPE,NXOPE,JXOPP,NXOPP,
     &                 IBEVC,IBPVC,IBTYP,WORK,LWORK)
C*****************************************************************************
C
C     Gather contribution into sigma-vector
C
C     Written by T.Saue Sep 18 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(DM1 = -1.0D0)
C
#include <dgroup.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dcbxrs.h>
      DIMENSION SEVEC(NXOPE,NZ,*),SPVEC(NXOPP,NZ,*),
     &          FMO(N2ORBXQ,*),JXOPE(*),JXOPP(*)
      DIMENSION WORK(*), DFAC(2)
      DIMENSION IBEVC(*),IBPVC(*),IBTYP(2,*)
Cold code
c     NXOPQ = NXOP*NZ
CTROND We assume first all NESIM and then all NPSIM.
c     DO I = 1,NFMAT
c       CALL XOPGAT(SVEC(1,1,I),JXOP,NXOP,FMO(1,I),NORBT,NZ)
c       CALL DSCAL(NXOPQ,DM1,SVEC(1,1,I),1)
c     ENDDO
Cnew code
#include <memint.h>
      LFPV = NZ*NASHT*NORBT
      CALL MEMGET('REAL',KFPV,LFPV,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KFPV),LFPV)
C
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHT*NZ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRXRS,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
C
      IREP = JSYMOP - 1
      IM = 0
C
C
C     Calculate trial vectors
C     =======================
C
      DO I = 1,NTSIM
C
C
C        Hermicity (used in RORBGRD)
C        ---------------------------
C
         IF ( I .LE. NESIM ) THEN
            IS = IBEVC(NERED+I)
            IS = IBTYP(2,IS)
            DFAC(1) = -1.0D00
            IF (IS .EQ. -1) THEN
               DFAC(2) = -1.0D00
            ELSE
               DFAC(2) = +1.0D00
            END IF
         ELSE
            IS = IBPVC(NPRED+I-NESIM)
            IS = IBTYP(2,IS)
            DFAC(1) = -1.0D00
            IF (IS .EQ. -1) THEN
               DFAC(2) = -1.0D00
            ELSE
               DFAC(2) = +1.0D00
            END IF
         END IF

C
C
C        Calculate orbital sigma vector
C        ------------------------------
C
         IF (NISHT .GT. 0) IM = IM + 1
C
C        Electronic
C
         IF (NXOPE .GT. 0) THEN
            CALL RORBGRD(WORK(KDV),FMO(1,IM),FMO(1,IM+1),DUM,
     &                   SEVEC(1,1,I),JXOPE,NXOPE,WORK(KFPV),DFAC,
     &                   IPQTOQ(1,IREP),IPRXRS)
         END IF
C
C        Positronic
C
         IF (NXOPP .GT. 0) THEN
            CALL RORBGRD(WORK(KDV),FMO(1,IM),FMO(1,IM+1),DUM,
     &                   SPVEC(1,1,I),JXOPP,NXOPP,WORK(KFPV),DFAC,
     &                   IPQTOQ(1,IREP),IPRXRS)
         END IF
         IF (NASHT .GT. 0) IM = IM + 1
      END DO
C
      CALL MEMREL('RORBGRD.KFPV',WORK,KFPV,KFPV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck setgmt */
      SUBROUTINE SETGMT(ISYMOP,IHRMOP,IFCKOP,IBEVC,IBPVC,IBTYP,
     &                  IPRINT)
C*****************************************************************************
C
C     Set arrays IHRMOP, ITIMOP, ISYMOP in preparation
C     for construction of G-matrices
C
C     Written by T.Saue Oct 15 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
#include <dgroup.h>
      DIMENSION IHRMOP(NOSIM),IFCKOP(NOSIM),ISYMOP(NOSIM),
     &          IBEVC(*),IBPVC(*),IBTYP(2,*)
C
      IM = 0      
C
C
C
C     e-e rotations
C     =============
C      
      DO I = 1,NESIM
         DO IMAT = 1,NFOCKMAT
           IM = IM + 1
           IS = IBEVC(NERED+I)
           IHRMOP(IM) = IBTYP(2,IS)
           IFCKOP(IM) = 1
           ISYMOP(IM) = JSYMOP
         ENDDO
      END DO
C
C
C
C     e-p rotations
C     =============
C
      DO I = 1,NPSIM
         DO IMAT = 1,NFOCKMAT
           IM = IM + 1
           IS = IBPVC(NPRED+I)
           IHRMOP(IM) = IBTYP(2,IS)
           IFCKOP(IM) = 1
           ISYMOP(IM) = JSYMOP
         ENDDO
      END DO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstvc */
      SUBROUTINE XRSTVC(BEVEC,IBEVC,NBEVC,IEOFF,
     &                  BPVEC,IBPVC,NBPVC,IPOFF,
     &                  BCVEC,IBCVC,NBCVC,ICOFF,
     &                  TEVEC,TPVEC,TCVEC,NBTOT,ITOFF,
     &                  IBTYP,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate tau vectors : metric S[2] times trial vector B(I)
C
C        TVEC(I) = S[2]*B(I)
C
C     The trial vectors come in three types:
C        C-type: contains only configurational parameters
C        E-type: contains only electron-eelectron (e-e) rotations
C        P-type: contains only electron-positron  (e-p) rotations
C
C     The tau vectors are correspondingly split up into a C, E and P -part.
C     The tau vectors are accumulated, so that one tau vector may
C     receive contributions from more than one trial vector.
C
C     Written by H.J.Aa.Jensen and Trond Saue 1996
C     Last revision Sep 19 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
C
      DIMENSION BCVEC(*),BEVEC(*),BPVEC(*),
     &          TCVEC(*),TEVEC(*),TPVEC(*),
     &          IBCVC(*),IBEVC(*),IBPVC(*),
     &          IBTYP(*),
     &          WORK(*)
C
#include <dcbxrs.h>
C
C
      CALL QENTER('XRSTVC')
C
C     Initialize tau vectors
C     ========================
C
      IF(NZCONF.GT.0) CALL DZERO(TCVEC,NBTOT*NZCONFQ)
      IF(NZXOPE.GT.0) CALL DZERO(TEVEC,NBTOT*NZXOPEQ)
      IF(NZXOPP.GT.0) CALL DZERO(TPVEC,NBTOT*NZXOPPQ)
C
C     Construct orbital part of tau vectors
C     =====================================
C
chj   IF (.NOT.RSPCI) THEN
          CALL XRSTEP(BEVEC,IBEVC,NBEVC,IEOFF,
     &                BPVEC,IBPVC,NBPVC,IPOFF,
     &                BCVEC,IBCVC,NBCVC,ICOFF,
     &                TEVEC,TPVEC,NBTOT,ITOFF,
     &                IBTYP,WORK,KFREE,LFREE)
chj   END IF
C
C     Construct configurational part of linear transformed vectors
C     for the configurational trial vectors
C     ============================================================
C
      IF (NCSIM.GT.0) THEN
         CALL XRSTCI(BCVEC,TCVEC,TEVEC,TPVEC,WORK,KFREE,LFREE)
      END IF
C
      CALL QEXIT('XRSTVC')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstci */
      SUBROUTINE XRSTCI(BCVEC,TCVEC,TEVEC,TPVEC,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate cc block of tau vector
C
C     Written by H.J.Aa.Jensen and T.Saue 1996
C     Last revision Oct 15 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION BCVEC(*),TCVEC(*),TEVEC(*), TPVEC(*),WORK(*)
C
C Used from common blocks:
C     XCBPRI: LUPRI
C     WRKXRS: NCSIM
C
#include <dcbxrs.h>
      CALL QENTER('XRSTCI')
      WRITE(LUPRI,'(A)') 
     &     '*** Sorry, XRSTCI called, but is not implemented yet!'
      CALL QUIT('XRSCLI called, but is not implemented yet!')
      CALL QEXIT ('XRSTCI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstov */
      SUBROUTINE XRSTEP(BEVEC,IBEVC,NBEVC,IEOFF,
     &                  BPVEC,IBPVC,NBPVC,IPOFF,
     &                  BCVEC,IBCVC,NBCVC,ICOFF,
     &                  TEVEC,TPVEC,NBTOT,ITOFF,
     &                  IBTYP,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Construct upper orbital part of tau vectors
C
C     TVEC_(mn) = S[2](o*o)_(mn,pq)U_(pq) + hS[2](o*o*)_(mn,pq)U^*_(pq)
C               = <0|[Nm,Pq]|0>U_(pq) - h<0|[Nm,Qp]|0>U^*_(pq)
C               = <0|[Nm,Pq]|0>W_(p,q)
C
C     where the matrix WA has the structure:
C
C           P  I   A   S             P  I  A  S
C          --------------         ---------------
C        P | 0  0   0   0       P |  0  U  U  0
C        I | 0  0 -hU+-hU+      I |-hU+ 0  0  0
C        A | 0  U   X -hU+      A |-hU+ 0  0  0
C        S | 0  U   U   0       S |  0  0  0  0
C
C           e-e rotations           e-p rotations
C
C     with the convention:
C        P - positronic spinors
C        I - inactive spinors
C        A - active spinors
C        S - secondary spinors
C
C     (X has -Y+ in upper diagonal and Z in lower diagonal)
C
C-----------------------------------------------------------------------------
C     In RPA we need to evaluate
C
C     TVEC_(si) = W_(si)
C
C     Written by T.Saue and H.J.Aa.Jensen 1996
C     Last revision Oct 15 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0, DM1 = -1.0D00)
      DIMENSION BEVEC(NZXOPEQ,*),IBEVC(*),
     &          BPVEC(NZXOPPQ,*),IBPVC(*),
     &          BCVEC(NZCONFQ,*),IBCVC(*),
     &          TEVEC(NZXOPEQ,*),TPVEC(NZXOPPQ,*),
     &          IBTYP(2,*),WORK(*)
      CHARACTER SSTR(2)*1
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      DATA SSTR /'S','A'/
C
C
C     Hermicity of operator
C     =====================
C
C
C     Allocate memory and get active density matrix
C     =============================================
C
      IREP = JSYMOP - 1
      KFRSAV = KFREE
      MNZXOP = MAX(NZXOPEQ,NZXOPPQ)
      CALL MEMGET('REAL',KUBVEC,N2ORBXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KOP   ,MNZXOP ,WORK,KFREE,LFREE)
      IF (NASHT .GT. 0) THEN
         CALL MEMGET('REAL',KDV,N2ASHT*NZ,WORK,KFREE,LFREE)
         CALL GETMAT('DVMO',WORK(KDV),IPRXRS,WORK(KFREE),LFREE)
      ELSE
         CALL MEMGET('REAL',KDV,0,WORK,KFREE,LFREE)
      END IF
C
C     Calculate tau vectors
C     =====================
C
C     Electronic...
C
      DO I = 1,NBEVC
         IH = IBEVC(I+IEOFF)
         IH = IBTYP(2,IH)
         IF (IH .EQ. -1) THEN
            IH  = 1
            IHH = 1
         ELSE IF (IH .EQ. 1) THEN
            IH  = -1
            IHH =  2
         ELSE
            CALL QUIT('argh ralle ralle')
         END IF
C
C        Construct full matrix
C        ---------------------
C
         CALL DZERO(WORK(KUBVEC),N2ORBXQ)
         CALL XOPSCT(BEVEC(1,I),WORK(KZXOPE),NZXOPE,
     &               WORK(KUBVEC),NORBT,NZ)

         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)') 
     &         'XRSTEP: Scattered trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IZOFF = 0
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            IS = IHQMAT(IQ,IH)
            CALL FULMAT(SSTR(IS),NORBT,NORBT,WORK(KUBVEC+IZOFF))
            IZOFF = IZOFF + N2ORBX
         END DO
         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)') 
     &         'XRSTEP: Full trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Calculate contribution to tau vector
C        ------------------------------------
C
         CALL DZERO(WORK(KOP),NZXOPEQ)
         CALL PRPORB(WORK(KOP),WORK(KZXOPE),NZXOPE,WORK(KUBVEC),
     &               NORBT,WORK(KDV),NZ,IHH)
C
C
         J  = IBEVC(I+IEOFF) - ITOFF
         IF (IPRXRS .GE. 20) THEN
             WRITE(LUPRI,'(A,I2)') 'XRSTEP: contribution to J = ',J
             WRITE(LUPRI,'(I5,3X,1P,E18.10)')
     &         (K,WORK(KOP+K-1), K = 1,NZXOPE)
         END IF
C
C
         CALL DAXPY(NZXOPEQ,D1,WORK(KOP),1,TEVEC(1,J),1)
C
      END DO
C
C
C     Positronic...
C
      DO I = 1,NBPVC
         IH = IBPVC(I+IPOFF)
         IH = IBTYP(2,IH)
         IF (IH .EQ. -1) THEN
            IH  = 1
            IHH = 1
         ELSE IF (IH .EQ. 1) THEN
            IH  = -1
            IHH =  2
         ELSE
            CALL QUIT('argh ralle ralle')
         END IF
C
C        Construct full matrix
C        ---------------------
C
         CALL DZERO(WORK(KUBVEC),N2ORBXQ)
         CALL XOPSCT(BPVEC(1,I),WORK(KZXOPP),NZXOPP,
     &               WORK(KUBVEC),NORBT,NZ)

         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)') 
     &         'XRSTEP: Scattered trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
         IZOFF = 0
         DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            IS = IHQMAT(IQ,IH)
            CALL FULMAT(SSTR(IS),NORBT,NORBT,WORK(KUBVEC+IZOFF))
            IZOFF = IZOFF + N2ORBX
         END DO
         IF (IPRXRS .GE. 100) THEN
            WRITE(LUPRI,'(1X,A,I3)') 
     &         'XRSTEP: Full trial vector no. ',I
            CALL PRQMAT(WORK(KUBVEC),NORBT,NORBT,NORBT,NORBT,
     &                  NZ,IPQTOQ(1,0),LUPRI)
         END IF
C
C        Calculate contribution to tau vector
C        ------------------------------------
C
         CALL DZERO(WORK(KOP),NZXOPPQ)
         CALL PRPORB(WORK(KOP),WORK(KZXOPP),NZXOPP,WORK(KUBVEC),
     &               NORBT,WORK(KDV),NZ,IHH)
C
C
         J  = IBPVC(I+IPOFF) - ITOFF
         IF (IPRXRS .GE. 20) THEN
             WRITE(LUPRI,'(A,I2)') 'XRSTEP: contribution to J = ',J
             WRITE(LUPRI,'(I5,3X,1P,E18.10)')
     &         (K,WORK(KOP+K-1), K = 1,NZXOPP)
         END IF
C
C
         CALL DAXPY(NZXOPPQ,D1,WORK(KOP),1,TPVEC(1,J),1)
C
      END DO
      CALL MEMREL('XRSTEP',WORK,1,KFRSAV,KFREE,LFREE)

Cnew code
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Xrsxvc */
      SUBROUTINE XRSXVC(ITYP,XEVEC,XPVEC,XCVEC,EVECR,MFREQ,
     &           IBTYP,IBCVC,IBEVC,IBPVC,BEBUF,BPBUF,BCBUF,
     &           WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Form solution vector in response equations.
C     The solution vector is split into a (+) and a (-) part.
C     In the static case the (-) part is automatically zero.
C
C     Written by T.Saue and H.J.Aa.Jensen Oct 16 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
#include <dcbibn.h>
      DIMENSION XEVEC(*),XPVEC(*),XCVEC(*),IBTYP(2,*),
     &          IBCVC(*),IBEVC(*),IBPVC(*),EVECR(*),
     &          BEBUF(*),BPBUF(*),BCBUF(*),WORK(*)
C
      CALL QENTER('XRSXVC')
      IF(NZXOPE.GT.0) CALL DZERO(XEVEC,MFREQ*NZXOPEQ)
      IF(NZXOPP.GT.0) CALL DZERO(XPVEC,MFREQ*NZXOPPQ)
      IF(NZCONF.GT.0) CALL DZERO(XCVEC,MFREQ*NZCONFQ)
C
C     e-e rotations
C      
      IF(NERED.GT.0) THEN
        CALL XRSXV1(ITYP,JBENDX,XEVEC,EVECR,MFREQ,
     &              IBTYP,IBEVC,IVECS,BEBUF)
      ENDIF
C
C     e-p rotations
C      
      IF(NPRED.GT.0) THEN
        CALL XRSXV1(ITYP,JBPNDX,XPVEC,EVECR,MFREQ,
     &              IBTYP,IBPVC,IVECS,BPBUF)
      ENDIF
C
C     Configurational parameters
C      
      IF(NCRED.GT.0) THEN
        CALL XRSXV1(ITYP,JBCNDX,XCVEC,EVECR,MFREQ,
     &              IBTYP,IBCVC,IVECS,BCBUF)
      ENDIF
C
      CALL QEXIT('XRSXVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Xrsxv1 */
      SUBROUTINE XRSXV1(ITYP,NBTYP,XVEC,EVECR,MFREQ,IBTYP,IBVEC,
     &                  IVECS,BVBUF)
C*****************************************************************************
C
C     Form electronic/positronic/configurational
C     part of solution vector in linear response
C
C     Written by T.Saue Oct 16 1996
C
C*****************************************************************************
#include <priunit.h>
#include <implicit.h>
C
#include <dcbxrs.h>      
#include <dcbibn.h>
      PARAMETER(DTOL=1.0D-10)
      DIMENSION XVEC(*),EVECR(NREDM,MFREQ),
     &               IBTYP(2,*),IBVEC(*),IVECS(*),BVBUF(*)
C
C
      IF    (NBTYP.EQ.JBCNDX) THEN
        NBDIM  = NZCONFQ
        NBRED  = NCRED
        LUBVEC = LUBCI
      ELSEIF(NBTYP.EQ.JBENDX) THEN
        NBDIM  = NZXOPEQ
        NBRED  = NERED
        LUBVEC = LUBOE
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBDIM  = NZXOPPQ
        NBRED  = NPRED
        LUBVEC = LUBOP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'XRSXVC ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('XRSXVC: Unknown NBTYP !')
      ENDIF
C
      CALL IZERO(IVECS,MFREQ)
      CALL DZERO(XVEC,NBDIM*MFREQ)
      DO IREC = 1,NBRED
        IND = IBVEC(IREC)
        IH  = IBTYP(2,IND)
        IF(IH.EQ.ITYP) THEN
          CALL READAC(LUBVEC,NBDIM,BVBUF,IREC)
          IBOFF = 1
          DO IFREQ = 1,MFREQ
            IF(ABS(EVECR(IND,IFREQ)).GT.DTOL) THEN
              IVECS(IFREQ) = IVECS(IFREQ) + 1
              CALL DAXPY(NBDIM,EVECR(IND,IFREQ),
     &                     BVBUF,1,XVEC(IBOFF),1)
            ENDIF
            IBOFF = IBOFF + NBDIM
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ACCS2X */
      SUBROUTINE ACCS2X(BEPVC,BEMVC,BPPVC,BPMVC,BCPVC,BCMVC,
     &                  IBTYP,IBCVC,IBEVC,IBPVC,
     &                  EVECR,EVALR,MFREQ,BEBUF,BPBUF,BCBUF,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     In the construction of the residual accumulate -W*S[2]X
C     for all frequencies W.
C
C     Input:
C       B(+)-vec contains current X(-) solution vectors
C       B(-)-vec contains current X(+) solution vectors
C
C     Output:
C       B(+)-vec contains T(+) = -W*S[2]X(-)
C       B(-)-vec contains T(+) = -W*S[2]X(+)
C
C     Written by T.Saue Oct 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0)
C
#include <dcbxrs.h>
#include <dcbibn.h>
      DIMENSION BEPVC(NZXOPEQ,*),BEMVC(NZXOPEQ,*),BEBUF(NZXOPEQ),
     &          BPPVC(NZXOPPQ,*),BPMVC(NZXOPPQ,*),BPBUF(NZXOPPQ),
     &          BCPVC(NZCONFQ,*),BCMVC(NZCONFQ,*),BCBUF(NZCONFQ),
     &          IBTYP(2,*),IBCVC(*),IBEVC(*),IBPVC(*),
     &          EVECR(NREDM,MFREQ),EVALR(MFREQ),WORK(*)
C
      CALL QENTER('ACCS2X')
C
      NBEVC = MIN(NZXOPE,1)
      NBPVC = MIN(NZXOPP,1)
      NBCVC = MIN(NZCONF,1)
C
      DO I = 1,MFREQ
        FAC = -EVALR(I)
        IF(FAC.NE.D0) THEN
C
C         Accumulate T(-) = S[2}X(+) in buffer
C         ------------------------------------
C          
          CALL XRSTVC(BEPVC(1,I),1,NBEVC,0,
     &                BPPVC(1,I),1,NBPVC,0,
     &                BCPVC(1,I),1,NBCVC,0,
     &                BEBUF,BPBUF,BCBUF,1,0,
     &                IBTYP,WORK,KFREE,LFREE)
C
C         Accumulate T(+) = S[2}X(-) in b(+)
C         ------------------------------------
C
          CALL XRSTVC(BEMVC(1,I),1,NBEVC,0,
     &                BPMVC(1,I),1,NBPVC,0,
     &                BCMVC(1,I),1,NBCVC,0,
     &                BEPVC(1,I),BPPVC(1,I),BCPVC(1,I),1,0,
     &                IBTYP,WORK,KFREE,LFREE)
C
C         Scale T = S[2]X with negative frequency
C         ---------------------------------------
C
          IF(NZXOPE.GT.0) THEN
            CALL DCOPY(NZXOPEQ,BEBUF,1,BEMVC(1,I),1)              
            CALL DSCAL(NZXOPEQ,FAC,BEMVC(1,I),1)
            CALL DSCAL(NZXOPEQ,FAC,BEPVC(1,I),1)
          ENDIF          
          IF(NZXOPP.GT.0) THEN
            CALL DCOPY(NZXOPPQ,BPBUF,1,BPMVC(1,I),1)              
            CALL DSCAL(NZXOPPQ,FAC,BPMVC(1,I),1)
            CALL DSCAL(NZXOPPQ,FAC,BPPVC(1,I),1)
          ENDIF          
          IF(NZCONF.GT.0) THEN
            CALL DCOPY(NZCONFQ,BCBUF,1,BCMVC(1,I),1)              
            CALL DSCAL(NZCONFQ,FAC,BCMVC(1,I),1)
            CALL DSCAL(NZCONFQ,FAC,BCPVC(1,I),1)
          ENDIF          
        ELSE
C
C       IF FREQ = 0, then R(-) is automatically zero
C       --------------------------------------------
C
          IF(NZXOPE.GT.0) THEN
            CALL DZERO(BEPVC(1,I),NZXOPEQ)
            CALL DZERO(BEMVC(1,I),NZXOPEQ)
          ENDIF
          IF(NZXOPP.GT.0) THEN
            CALL DZERO(BPPVC(1,I),NZXOPPQ)
            CALL DZERO(BPMVC(1,I),NZXOPPQ)
          ENDIF
          IF(NZCONF.GT.0) THEN
            CALL DZERO(BCPVC(1,I),NZCONFQ)
            CALL DZERO(BCMVC(1,I),NZCONFQ)
          ENDIF
        ENDIF
      ENDDO
C
      CALL QEXIT('ACCS2X')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ADDE2X */
      SUBROUTINE ADDE2X(ITYP,BEVEC,BPVEC,BCVEC,EVECR,MFREQ,
     &            IBTYP,BEBUF,BPBUF,BCBUF,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     In the construction of the residual vector, add E[2}X.
C
C     Written by T.Saue Oct 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
#include <dcbibn.h>
      DIMENSION BEVEC(*),BPVEC(*),BCVEC(*),IBTYP(2,*),
     &          EVECR(*),BEBUF(*),BPBUF(*),BCBUF(*),WORK(*)
C
      CALL QENTER('ADDSIG')
C
C     e-e rotations
C      
      IF(NERED.GT.0) THEN
        CALL ADDSV1(ITYP,JBENDX,BEVEC,EVECR,MFREQ,
     &              IBTYP,BEBUF)
      ENDIF
C
C     e-p rotations
C      
      IF(NPRED.GT.0) THEN
        CALL ADDSV1(ITYP,JBPNDX,BPVEC,EVECR,MFREQ,
     &              IBTYP,BPBUF)
      ENDIF
C
C     Configurational parameters
C      
      IF(NCRED.GT.0) THEN
        CALL ADDSV1(ITYP,JBCNDX,BCVEC,EVECR,MFREQ,
     &              IBTYP,BCBUF)
      ENDIF
C
      CALL QEXIT('ADDSIG')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck ADDSV1 */
      SUBROUTINE ADDSV1(ITYP,NBTYP,BVEC,EVECR,MFREQ,IBTYP,BVBUF)
C*****************************************************************************
C
C     Form electronic/positronic/configurational
C     part of solution vector in linear response
C
C     Written by T.Saue Oct 16 1996
C
C*****************************************************************************
#include <priunit.h>
#include <implicit.h>
C
#include <dcbxrs.h>      
#include <dcbibn.h>
      DIMENSION BVEC(*),EVECR(NREDM,MFREQ),
     &               IBTYP(2,*),BVBUF(*)
C
C
      IF    (NBTYP.EQ.JBCNDX) THEN
        NBDIM  = NZCONFQ
        LUBVEC = LUSCI
      ELSEIF(NBTYP.EQ.JBENDX) THEN
        NBDIM  = NZXOPEQ
        LUBVEC = LUSOE
      ELSEIF(NBTYP.EQ.JBPNDX) THEN
        NBDIM  = NZXOPPQ
        LUBVEC = LUSOP
      ELSE
        WRITE(LUPRI,'(A,A,I5)') 'ADDSV1 ERROR: ',
     &    'unknown NBTYP = ',NBTYP
        CALL QUIT('ADDSV1: Unknown NBTYP !')
      ENDIF
C
      DO IREC = 1,NZRED
        IH  = IBTYP(2,IREC)
        IF(IH.EQ.ITYP) THEN
          CALL READAC(LUBVEC,NBDIM,BVBUF,IREC)
          IBOFF = 1
          DO IFREQ = 1,MFREQ
            CALL DAXPY(NBDIM,EVECR(IREC,IFREQ),BVBUF,1,BVEC(IBOFF),1)
            IBOFF = IBOFF + NBDIM
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XVCNRM */
      SUBROUTINE XVCNRM(XNRM,XVEC,NDIM,MFREQ)
C*****************************************************************************
C
C     Take norm of resultant vector X
C
C     Written by T.Saue Oct 17 1996
C
C     Last revision Oct 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION XVEC(NDIM,*),XNRM(*)
C
      DO I = 1,MFREQ
        XNRM(I) = DNRM2(NDIM,XVEC(1,I),1)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/*Deck xrsne1 */
      SUBROUTINE XRSNE1(BEVEC,BPVEC,BCVEC,MFREQ,EVALR,EVECR,
     &                  GPOE,GPOP,GPCI,IBTYP,IBCVC,IBEVC,IBPVC,
     &                  BEBUF,BPBUF,BCBUF,RCNV,RNORM,XNORM,IVECS,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     This routine will form new trial vectors for the
C     next microiteration. 
C
C     Written by T.Saue Oct 18 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbxrs.h>
#include <dcbibn.h>
      DIMENSION BEVEC(*),BPVEC(*),BCVEC(*),EVALR(*),
     &          EVECR(NREDM,*),GPOE(*),GPOP(*),GPCI(*),
     &          IBTYP(2,NREDM),IBCVC(NREDM),IBEVC(NREDM),IBPVC(NREDM),
     &          BEBUF(*),BPBUF(*),BCBUF(*),IVECS(MFREQ,NVPAR),
     &          RNORM(NVPAR,MFREQ),XNORM(MFREQ,NVPAR),WORK(*)
C
      IF(IPRXRS.GE.3) CALL HEADER('Output from XRSNE1',-1)
C
C     ***********************************
C     ***** Form solution vectors X *****
C     ***********************************
C      
      IND = NZRED
C
C     Static case
C     ===========
C
      IF(STATIC) THEN
C
C       Orbital (e-e) part
C       ------------------
C
        IF(NZXOPE.GT.0) THEN
          NESIM = MFREQ
          IF((IND+NESIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          DO I = 1,NESIM
            IND = IND + 1
            IBTYP(1,IND)   = JBENDX
            IBTYP(2,IND)   = IPX
            IBEVC(NERED+I) = IND
          ENDDO
          CALL DZERO(BEVEC,NESIM*NZXOPEQ)
          IF(NERED.GT.0) THEN
            CALL XRSXV1(IPX,JBENDX,BEVEC,EVECR,MFREQ,
     &                  IBTYP,IBEVC,IVECS(1,IPEP),BEBUF)
            CALL XVCNRM(XNORM(1,IPEP),BEVEC,NZXOPEQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-e) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPEP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-e) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPEP),MFREQ)
          ENDIF
        ENDIF
C
C       Orbital (e-p) part
C       ------------------
C
        IF(NZXOPP.GT.0) THEN
          NPSIM = MFREQ
          IF((IND+NPSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          DO I = 1,NPSIM
            IND = IND + 1
            IBTYP(1,IND)   = JBPNDX
            IBTYP(2,IND)   = IPX
            IBPVC(NPRED+I) = IND
          ENDDO
          CALL DZERO(BPVEC,NPSIM*NZXOPPQ)
          IF(NPRED.GT.0) THEN
            CALL XRSXV1(IPX,JBPNDX,BPVEC,EVECR,MFREQ,
     &                  IBTYP,IBPVC,IVECS(1,IPPP),BPBUF)
            CALL XVCNRM(XNORM(1,IPPP),BPVEC,NZXOPPQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-p) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPPP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-p) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPPP),MFREQ)
          ENDIF
        ENDIF
C
C       Configurational part
C       --------------------
C
        IF(NZCONF.GT.0) THEN
          NCSIM = MFREQ         
          IF((IND+NCSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          DO I = 1,NCSIM
            IND = IND + 1
            IBTYP(1,IND)   = JBCNDX
            IBTYP(2,IND)   = IPX
            IBCVC(NCRED+I) = IND
          ENDDO
          CALL DZERO(BCVEC,NCSIM*NZCONFQ)
          IF(NCRED.GT.0) THEN
            CALL XRSXV1(IPX,JBCNDX,BCVEC,EVECR,MFREQ,
     &                  IBTYP,IBCVC,IVECS(1,IPCP),BCBUF)
            CALL XVCNRM(XNORM(1,IPCP),BCVEC,NZCONFQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of configurational X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPCP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' configurational X(+) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPCP),NCSIM)
          ENDIF
        ENDIF
      ELSE
C
C       In frequency-dependent case, the resultant X-vectors
C       split into X(+) and X(-). 
C
        KBE2 = 1 + MFREQ*NZXOPEQ
        KBP2 = 1 + MFREQ*NZXOPPQ
        KBC2 = 1 + MFREQ*NZCONFQ
        IXE2 = MFREQ
        IXP2 = MFREQ
        IXC2 = MFREQ
C
C       Orbital (e-e) part
C       ------------------
C
        IF(NZXOPE.GT.0) THEN
          NESIM = MFREQ + MFREQ
          IF((IND+NESIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          IOFF = NERED
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)   = JBENDX
            IBTYP(2,IND)   = IPX
            IBEVC(IOFF+I) = IND
          ENDDO
          IOFF = NERED + MFREQ
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)   = JBENDX
            IBTYP(2,IND)   = IMX
            IBEVC(IOFF+I) = IND
          ENDDO
          CALL DZERO(BEVEC,NESIM*NZXOPEQ)
          IF(NERED.GT.0) THEN
            CALL XRSXV1(IPX,JBENDX,BEVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBEVC,IVECS(1,IPEP),BEBUF)
            CALL XRSXV1(IMX,JBENDX,BEVEC(KBE2),EVECR,MFREQ,
     &                  IBTYP,IBEVC,IVECS(1,IPEM),BEBUF)
            CALL XVCNRM(XNORM(1,IPEP),BEVEC      ,NZXOPEQ,MFREQ)
            CALL XVCNRM(XNORM(1,IPEM),BEVEC(KBE2),NZXOPEQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-e) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPEP),I=1,MFREQ)
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-e) X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPEM),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-e) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BEVEC      ,MFREQ,NZXOPE)
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-e) X(-) solution vectors'
              CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPEP),MFREQ)
            CALL DZERO(XNORM(1,IPEM),MFREQ)
          ENDIF
        ENDIF
C
C       Orbital (e-p) part
C       ------------------
C
        IF(NZXOPP.GT.0) THEN
          NPSIM = MFREQ + MFREQ
          IF((IND+NPSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          IOFF = NPRED
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)   = JBPNDX
            IBTYP(2,IND)   = IPX
            IBPVC(IOFF+I) = IND
          ENDDO
          IOFF = NPRED + MFREQ
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)   = JBPNDX
            IBTYP(2,IND)   = IMX
            IBPVC(IOFF)    = IND
          ENDDO
          CALL DZERO(BPVEC,NPSIM*NZXOPPQ)
          IF(NPRED.GT.0) THEN
            CALL XRSXV1(IPX,JBPNDX,BPVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBPVC,IVECS(1,IPPP),BPBUF)
            CALL XRSXV1(IMX,JBPNDX,BPVEC(KBP2),EVECR,MFREQ,
     &                  IBTYP,IBPVC,IVECS(1,IPPM),BPBUF)
            CALL XVCNRM(XNORM(1,IPPP),BPVEC      ,NZXOPPQ,MFREQ)
            CALL XVCNRM(XNORM(1,IPPM),BPVEC(KBP2),NZXOPPQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-p) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPPP),I=1,MFREQ)
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-p) X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPPM),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-p) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-p) X(-) solution vectors'
              CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPPP),MFREQ)
            CALL DZERO(XNORM(1,IPPM),MFREQ)
          ENDIF
        ENDIF
C
C       Configurational part
C       --------------------
C
        IF(NZCONF.GT.0) THEN
          NCSIM = MFREQ + MFREQ
          IF((IND+NCSIM).GT.NREDM) THEN
            WRITE(LUPRI,'(A,I4,A,I4)')
     &        'XRSNE1 ERROR: IBTYP beyond NREDM(',NREDM,
     &        ') in microiteration no.',ITMIC
            CALL QUIT('XRSNE1: Exceeded reduced space')
          ENDIF
          IOFF = NCRED
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)   = JBCNDX
            IBTYP(2,IND)   = IPX
            IBCVC(IOFF+I) = IND
          ENDDO
          IOFF = NCRED + MFREQ
          DO I = 1,MFREQ
            IND = IND + 1
            IBTYP(1,IND)   = JBCNDX
            IBTYP(2,IND)   = IMX
            IBCVC(IOFF) = IND
          ENDDO
          CALL DZERO(BCVEC,NCSIM*NZCONFQ)
          IF(NCRED.GT.0) THEN
            CALL XRSXV1(IPX,JBCNDX,BCVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBCVC,IVECS(1,IPCP),BCBUF)
            CALL XRSXV1(IMX,JBCNDX,BCVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBCVC,IVECS(1,IPCM),BCBUF)
            CALL XVCNRM(XNORM(1,IPCP),BCVEC      ,NZCONFQ,MFREQ)
            CALL XVCNRM(XNORM(1,IPCM),BCVEC(KBC2),NZCONFQ,MFREQ)
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of configurational X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPCP),I=1,MFREQ)
              WRITE(LUPRI,'(A)') 
     &          'Norms of configurational X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XNORM(I,IPCP),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' configurational X(+) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' configurational X(-) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
            ENDIF
          ELSE
            CALL DZERO(XNORM(1,IPCP),MFREQ)
            CALL DZERO(XNORM(1,IPCM),MFREQ)
          ENDIF
        ENDIF
      ENDIF
C
      NTSIM = NESIM + NPSIM + NCSIM
C
C     Print section
C
      IF(IPRXRS.GT.0) THEN
        WRITE(LUPRI,'(/A)') '* XRSNEX: Norms of resultant vectors:'
        WRITE(LUPRI,'(A15,2X,A8,9(2X,A2,A6))')  
     &      '* Frequency    ','Total   ',(BVTYP(I),' part ',I=1,NVPAR)
        DO IFREQ = 1,MFREQ
          XTNORM = DNRM2(NVPAR,XNORM(IFREQ,1),MFREQ)
          WRITE(LUPRI,'(F15.8,1P,10(2X,E8.2))')
     &        EVALR(IFREQ),XTNORM,(XNORM(IFREQ,I),I=1,NVPAR)
        ENDDO
      ENDIF
      CALL FLSHFO(LUPRI)
C
C     ****************************
C     ***** Form residuals R *****
C     ****************************
C      
C
C     In frequency-dependent case, form -W*S[2]X
C     ==========================================
C
      IF(STATIC) THEN
        IF(NESIM.GT.0) CALL DZERO(BEVEC,NZXOPEQ*NESIM)
        IF(NPSIM.GT.0) CALL DZERO(BPVEC,NZXOPPQ*NPSIM)
        IF(NCSIM.GT.0) CALL DZERO(BCVEC,NZCONFQ*NCSIM)
      ELSE
        CALL ACCS2X(BEVEC,BEVEC(KBE2),BPVEC,BPVEC(KBP2),
     &              BCVEC,BCVEC(KBC2),IBTYP,IBCVC,IBEVC,IBPVC,
     &              EVECR,EVALR,MFREQ,BEBUF,BPBUF,BCBUF,
     &              WORK,KFREE,LFREE)
      ENDIF
C
C     Print section
C
      IF(IPRXRS.GE.5) THEN
        IF(NZXOPE.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' orbital (e-e) -WS[2}X(-) residual vectors'
          CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' orbital (e-e) -WS[2}X(+) residual vectors'
            CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
          ENDIF
        ENDIF
        IF(NZXOPP.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' orbital (e-p) -WS[2}X(-) residual vectors'
          CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' orbital (e-p) -WS[2}X(+) residual vectors'
            CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
          ENDIF
        ENDIF
        IF(NZCONF.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' configurational -WS[2}X(-) residual vectors'
          CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' configurational -WS[2}X(+) residual vectors'
            CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
          ENDIF
        ENDIF
      ENDIF
C
C     Add E[2]X
C     =========
C
      CALL ADDE2X(IPX,BEVEC,BPVEC,BCVEC,EVECR,MFREQ,
     &            IBTYP,BEBUF,BPBUF,BCBUF,WORK,KFREE,LFREE)
      IF(.NOT.STATIC) THEN
        CALL ADDE2X(IMX,BEVEC(KBE2),BPVEC(KBP2),BCVEC(KBC2),EVECR,MFREQ,
     &              IBTYP,BEBUF,BPBUF,BCBUF,WORK,KFREE,LFREE)

      ENDIF
C
C     Print section
C
      IF(IPRXRS.GE.5) THEN
        IF(NZXOPE.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' orbital (e-e) (E[2}X(+)-WS[2]X(-)) residual vectors'
          CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' orbital (e-e) (E[2}X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
          ENDIF
        ENDIF
        IF(NZXOPP.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' orbital (e-p) (E[2}X(+)-WS[2]X(-)) residual vectors'
          CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' orbital (e-p) (E[2}X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
          ENDIF
        ENDIF
        IF(NZCONF.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' configurational (E[2}X(+)-WS[2]X(-)) residual vectors'
          CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       'configurational (E[2}X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
          ENDIF
        ENDIF
      ENDIF
C
C     If (LINEQ) add gradient
C     =======================
C
      IF(LINEQ) THEN
        IF(NZXOPE.GT.0) THEN
          CALL GRDADD(BEVEC,GPOE,NZXOPEQ,MFREQ)
        ENDIF
        IF(NZXOPP.GT.0) THEN
          CALL GRDADD(BPVEC,GPOP,NZXOPPQ,MFREQ)
        ENDIF
        IF(NZCONF.GT.0) THEN
          CALL GRDADD(BCVEC,GPCI,NZCONFQ,MFREQ)
        ENDIF
C
C       Print section
C
      IF(IPRXRS.GE.5) THEN
        IF(NZXOPE.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' orbital (e-e) (E[2}X(+)-WS[2]X(-)+G) residual vectors'
          CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' orbital (e-e) (E[2}X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BEVEC(KBE2),MFREQ,NZXOPE)
          ENDIF
        ENDIF
        IF(NZXOPP.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' orbital (e-p) (E[2}X(+)-WS[2]X(-)+G) residual vectors'
          CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' orbital (e-p) (E[2}X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BPVEC(KBP2),MFREQ,NZXOPP)
          ENDIF
        ENDIF
        IF(NZCONF.GT.0) THEN
          WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &     ' configurational (E[2}X(+)-WS[2]X(-)+G) residual vectors'
          CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
          IF(.NOT.STATIC) THEN
            WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &       ' configurational (E[2}X(-)-WS[2]X(+)) residual vectors'
            CALL PRBVEC(LUPRI,BCVEC(KBC2),MFREQ,NZCONF)
          ENDIF
        ENDIF
      ENDIF
      ENDIF
C
C     **********************************
C     ***** Form new trial vectors *****
C     **********************************
C
      CALL XRSBVC(BEVEC,BEBUF,IBEVC,
     &            BPVEC,BPBUF,IBPVC,
     &            BCVEC,BCBUF,IBCVC,
     &            MFREQ,IBTYP,IVECS,EVALR,EVECR,RCNV,RNORM,XNORM,
     &            WORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GRDADD */
      SUBROUTINE GRDADD(BVEC,GP,NBDIM,MFREQ)
C*****************************************************************************
C
C     Add property gradient to residual
C
C     Written by T.Saue Oct 17 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>   
      PARAMETER(D1 = 1.0D0)
C
      DIMENSION BVEC(NBDIM,MFREQ),GP(NBDIM)
C
      DO I = 1,MFREQ
        CALL DAXPY(NBDIM,D1,GP,1,BVEC(1,I),1)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XRSBVC */
      SUBROUTINE XRSBVC(BEVEC,BEBUF,IBEVC,
     &                  BPVEC,BPBUF,IBPVC,
     &                  BCVEC,BCBUF,IBCVC,
     &                  MFREQ,IBTYP,IVECS,EVALR,EVECR,RCNV,RNORM,XNORM,
     &                  WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Construct new trial vectors for next microiteration
C
C     The selection of trial vectors proceed in three steps:
C       1) If |R| .LE. (THRS*|X|) 
C             where |R| is norm of residual
C             and   |X| is norm of current solution vector
C          then no new trial vectors are generated from this root.
C       2) If |R*| .LE. (RESFAC*|RMAX|)
C             where |R*| is the norm of a partial residual (E-,P-,C-type)
C             and |RMAX| is the largest partial residual
C          then no new trial vectors of this type is generated from this root.
C          If RESFAC = 1, then only one trial vector is generated from each root.
C          If RESFAC = 0, then all possible trial vectors are generated from 
C                         the root.
C
C     Input:
C       BEVEC - residual : orbital (e-e) part
C       BPVEC - residual : orbital (e-p) part
C       BCVEC - residual : configurational part
C
C       XECNV - norm of resultant vector: orbital (e-e) part
C       XPCNV - norm of resultant vector: orbital (e-p) part
C       XCCNV - norm of resultant vector: configurational part
C
C     Output:
C
C       BEVEC - NESIM trial vectors: orbital (e-e) type
C       BPVEC - NPSIM trial vectors: orbital (e-p) type
C       BCVEC - NCSIM trial vectors : configurational part
C
C       RNORM - residual norm
C
C     Written by T.Saue, Oct 18 1996
C     Last modification Nov 26 1997
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D5=5.0D0)
C
#include <dcbxrs.h>
#include <dcbibn.h>
      LOGICAL   USEXVC
      DIMENSION BEVEC(NZXOPEQ,*),BEBUF(*),IBEVC(NREDM),
     &          BPVEC(NZXOPPQ,*),BPBUF(*),IBPVC(NREDM),
     &          BCVEC(NZCONFQ,*),BCBUF(*),IBCVC(NREDM),
     &          IBTYP(2,NREDM),IVECS(MFREQ,NVPAR),EVALR(*),EVECR(*),
     &          RCNV(*),RNORM(NVPAR,MFREQ),XNORM(MFREQ,NVPAR),
     &          WORK(*)
C
      CNVMAX = D0
      WRITE(LUPRI,'(/A)') '* Norm of residual vectors:'
      WRITE(LUPRI,'(A15,2X,A8,9(2X,A2,A6))')  
     &      '  Frequency    ','Total   ',(BVTYP(I),' part ',I=1,NVPAR)
      DO IFREQ = 1,MFREQ
C
C      ************************************
C      ***** Get norms of residual vectors 
C      ************************************
C
C       Orbital (e-e) part
C       ------------------
C
        IF(NZXOPE.GT.0) THEN
          IVECS(IFREQ,IPEP) = 1
          RNORM(IPEP,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,IFREQ),1)
          IF(.NOT.STATIC) THEN
            IVECS(IFREQ,IPEM) = 1
            IOFF = IFREQ + MFREQ
            RNORM(IPEM,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,IOFF),1)
          ENDIF           
        ENDIF
C
C       Orbital (e-p) part
C       ------------------
C
        IF(NZXOPP.GT.0) THEN
          IVECS(IFREQ,IPPP) = 1
          RNORM(IPPP,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,IFREQ),1)
          IF(.NOT.STATIC) THEN
            IOFF = IFREQ + MFREQ
            IVECS(IFREQ,IPPM) = 1
            RNORM(IPPM,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,IOFF),1)
          ENDIF           
        ENDIF
C
C       Configurational part
C       --------------------
C
        IF(NZCONF.GT.0) THEN
          IVECS(IFREQ,IPCP) = 1
          RNORM(IPCP,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,IFREQ),1)
          IF(.NOT.STATIC) THEN
            IOFF = IFREQ + MFREQ
            IVECS(IFREQ,IPCM) = 1
            RNORM(IPCM,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,IOFF),1)
          ENDIF           
        ENDIF
C
C       *****************************
C       ***** Check convergence *****
C       *****************************
C       
        RNRMAX = D0
        DO I = 1,NVPAR
          RNRMAX = MAX(RNRMAX,RNORM(I,IFREQ))
        ENDDO
        RTNORM = DNRM2(NVPAR,RNORM(1,IFREQ),1)
        XTNORM = DNRM2(NVPAR,XNORM(IFREQ,1),MFREQ)
        IF(XTNORM.EQ.D0) THEN
          WRITE(LUPRI,'(A,A)') '*** WARNING! ***',
     &        'Zero norm of solution vector !'
        ELSE
          RCNV(IFREQ) = RTNORM/GPTNRM
          CNVMAX = MAX(CNVMAX,RCNV(IFREQ))
        ENDIF
        IF(RTNORM.LE.(GPTNRM*THCXRS)) THEN
C         ... this root converged
          DO I = 1,NVPAR
            IVECS(IFREQ,I) = 0
          ENDDO
        ELSE
C       hjj+patrick:
C       Check each type separately whether residual norm is less than 
C       a fraction of the maximun residual norm, or if less than
C       a fraction (0.2) of convergence threshold.
          DO I = 1,NVPAR
            IF(RESFAC*RNORM(I,IFREQ).LT.RNRMAX.OR.
     &         D5*RNORM(I,IFREQ).LE.GPTNRM*THCXRS) IVECS(IFREQ,I) = 0
          ENDDO
        ENDIF
        WRITE(LUPRI,'(F15.8,1P,10(2X,E8.2))')
     &        EVALR(IFREQ),RTNORM,(RNORM(I,IFREQ),I=1,NVPAR)
      ENDDO
C
C     Check what integral classes contribute
C     ======================================
C
      ITNEX  = ITMIC + 1
      INTBUF = INTFLG
      CALL INTCON(INTFLG,INTBUF,INTDEF,
     &            CNVMAX,CNVINT,ITNEX,ITRINT,INTTYP)
      USEXVC = INTFLG.NE.INTBUF
C
C     Check whether calculation has converged
C     =======================================
C
      IF(.NOT.USEXVC) THEN
        KCONV = ISUM(NTSIM,IVECS,1)
        IF(KCONV.EQ.0) RETURN
      ENDIF
C
C     If maximum number of iterations passed, return
C     ==============================================
C
      IF(ITMIC.GE.MAXITR) RETURN
C
C     ***************************************
C     ***** Generate new trial vectors  *****
C     ***** (update IBTYP and IB*VC)    *****
C     ***************************************
C      
C
C     If requested, allocate space for diagonal elements of Hessian (approximate)
C
      KFRSAV = KFREE
      KEDIAH = KFREE
      KPDIAH = KFREE
      KCDIAH = KFREE
      IF(DIAGHE.AND.(.NOT.USEXVC)) THEN
         IF (NZXOPE.GT.0) THEN
            CALL MEMGET('REAL',KEDIAH,NZXOPE,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KEDIAS,NZXOPE,WORK,KFREE,LFREE)
         END IF
         IF (NZXOPP.GT.0) THEN
            CALL MEMGET('REAL',KPDIAH,NZXOPP,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KPDIAS,NZXOPP,WORK,KFREE,LFREE)
         END IF
         IF (NZCONF.GT.0) THEN
            CALL MEMGET('REAL',KCDIAH,NZCONF,WORK,KFREE,LFREE)
            CALL MEMGET('REAL',KCDIAS,NZCONF,WORK,KFREE,LFREE)
         END IF
         CALL E2DIAG(WORK(KEDIAH),WORK(KPDIAH),WORK(KCDIAH),
     &               WORK(KEDIAS),WORK(KPDIAS),WORK(KCDIAS),
     &               WORK,KFREE,LFREE)
      ENDIF
      IF(USEXVC) THEN
        IBOFF  = 0
        NEPREV = 0
        NPPREV = 0
        NCPREV = 0
      ELSE
        IBOFF = NZRED
        NEPREV = NERED
        NPPREV = NPRED
        NCPREV = NCRED
      ENDIF
C
C     Orbital (e-e) type trial vectors
C
      IF(NZXOPE.GT.0) THEN
        CALL XBVGEN(BEVEC,USEXVC,WORK(KEDIAH),EVALR,EVECR,
     &              IBTYP,IBEVC,IVECS(1,IPEP),
     &              MFREQ,JBENDX,NESIM,IBOFF,NEPREV,NZXOPE,
     &              BEBUF,WORK,KFREE,LFREE)
        IBOFF = IBOFF + NESIM
      ENDIF
C
C     Orbital (e-p) type trial vectors
C
      IF(NZXOPP.GT.0) THEN
        CALL XBVGEN(BPVEC,USEXVC,WORK(KPDIAH),EVALR,EVECR,
     &              IBTYP,IBPVC,IVECS(1,IPPP),
     &              MFREQ,JBPNDX,NPSIM,IBOFF,NPPREV,NZXOPP,
     &              BPBUF,WORK,KFREE,LFREE)
        IBOFF = IBOFF + NPSIM
      ENDIF
C
C     Configuration type trial vectors
C
      IF(NZCONF.GT.0) THEN
        CALL XBVGEN(BCVEC,USEXVC,WORK(KCDIAH),EVALR,EVECR,
     &              IBTYP,IBCVC,IVECS(1,IPCP),
     &              MFREQ,JBCNDX,NCSIM,IBOFF,NCPREV,NZCONF,
     &              BCBUF,WORK,KFREE,LFREE)
        IBOFF = IBOFF + NCSIM
      ENDIF
C
C     Check final number of new trial vectors
C     =======================================
C
      NTSIM = NESIM + NPSIM + NCSIM
      IF(NTSIM.EQ.0) THEN
        KCONV = -1
        RETURN
      ENDIF
C
C     Check dimension of reduced space
C     ================================
C
      IF(USEXVC) RETURN
      MREDM = NZRED + NTSIM
      IF(MREDM.GT.NREDM) THEN
        KCONV = -2
        RETURN
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XBVGEN */
      SUBROUTINE XBVGEN(BVEC,USEXVC,HDIAG,EVALR,EVECR,IBTYP,IBVEC,IVECS,
     &            MFREQ,JBTYP,NBNEW,IBOFF,NBPREV,NVAR,BBUF,
     &            WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Generate new vectors
C
C     Written by T.saue Dec 1 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1 = -1.0D0)
C
#include <dcbxrs.h>
#include <dgroup.h>
C
      LOGICAL USEXVC,TEST
      DIMENSION BVEC(NVAR,NZ,*),HDIAG(*),
     &          EVALR(MFREQ),IBTYP(2,*),IBVEC(*),IVECS(*),
     &          BBUF(*),WORK(*)
C
      TEST = DIAGHE.AND.(.NOT.USEXVC)
      NVARQ = NVAR*NZ
      IF(LINEQ) THEN
C
C     L I N E A R   R E S P O N S E   F U N C T I O N S
C     =================================================
C
        IF(USEXVC) THEN
C         New integral class added; 
C         start from present solution vector
          CALL XRSXV1(IPX,JBTYP,BVEC,EVECR,MFREQ,
     &                    IBTYP,IBVEC,IVECS,BBUF)
          DO IP = 1,MFREQ
            IBTYP(1,IBOFF+IP) = JBTYP
            IBTYP(2,IBOFF+IP) = IPX
            IBVEC(NBPREV +IP) = IBOFF + IP
            IF(IVECS(IP).GT.0) IVECS(IP) = 1
          ENDDO
          IF(.NOT.STATIC) THEN
            CALL XRSXV1(IMX,JBTYP,BVEC(1,1,(MFREQ+1)),EVECR,MFREQ,
     &                    IBTYP,IBVEC,IVECS(MFREQ+1),BBUF)
            DO IP = 1,MFREQ
              IM = IP + MFREQ
              FREQ = EVALR(IP)
              IBTYP(1,IBOFF+IM) = JBTYP
              IBTYP(2,IBOFF+IM) = IMX
              IBVEC(NBPREV +IM) = IBOFF + IM
              IF(FREQ.EQ.D0) THEN
                IVECS(IM) = 0
              ELSE
                IF(IVECS(IM).GT.0) IVECS(IM) = 1
              ENDIF
            ENDDO
          ENDIF
          GOTO 10
        ENDIF
        DO IP = 1,MFREQ
          IM = IP + MFREQ
          FREQ = EVALR(IP)
          IF(FREQ.EQ.D0) THEN
C           X(+) - vector
            IBTYP(1,IBOFF+IP) = JBTYP
            IBTYP(2,IBOFF+IP) = IPX
            IBVEC(NBPREV +IP) = IBOFF + IP
            IF(TEST.AND.IVECS(IP).EQ.1) THEN
              BNRMP = DNRM2(NVARQ,BVEC(1,1,IP),1)
              DO IZ = 1,NZ
                DO IVAR = 1,NVAR
                  BVEC(IVAR,IZ,IP) = 
     &              (DM1/HDIAG(IVAR))*BVEC(IVAR,IZ,IP)
                ENDDO
              ENDDO
              BNRMP = BNRMP/DNRM2(NVARQ,BVEC(1,1,IP),1)
C             retain norm for numerical stability !
              CALL DSCAL(NVARQ,BNRMP,BVEC(1,1,IP),1)
            ENDIF
C           No X(-) - vector              
            IF(.NOT.STATIC) THEN
              IBTYP(1,IBOFF+IM) = JBTYP
              IBTYP(2,IBOFF+IM) = IMX
              IBVEC(NBPREV+IM)  = IBOFF + IM
              IVECS(IM)          = 0
            ENDIF
          ELSE
            IBTYP(1,IBOFF+IP)     = JBTYP
            IBTYP(2,IBOFF+IP)     = IPX
            IBVEC(NBPREV+IP)      = IBOFF + IP
            IBTYP(1,IBOFF+IM)     = JBTYP
            IBTYP(2,IBOFF+IM)     = IMX
            IBVEC(NBPREV+IM)      = IBOFF + IM
            IF(TEST) THEN
C
Chj-nov99 
C     Preconditioning: multiply with minus the elements of
C
C      / hdiag  -freq \ -1      / hdiag*fac   freq*fac \
C     |                |     = |                        |
C      \ -freq  hdiag /         \  freq*fac  hdiag*fac /
C
C     where fac = (hdiag**2 - freq**2)**(-1)
C     for the res(ip), res(im) pair
C
              FREQ2 = FREQ**2
              IF(IVECS(IP).EQ.1.AND.IVECS(IM).EQ.1) THEN
                BNRMP = DNRM2(NVARQ,BVEC(1,1,IP),1)
                BNRMM = DNRM2(NVARQ,BVEC(1,1,IM),1)
                DO IZ = 1,NZ
                  DO IVAR = 1,NVAR
                    FAC = D1/(FREQ2-HDIAG(IVAR)**2)
                    BP = HDIAG(IVAR)*FAC * BVEC(IVAR,IZ,IP)
     &                 + FREQ * FAC * BVEC(IVAR,IZ,IM)
                    BM = HDIAG(IVAR)*FAC * BVEC(IVAR,IZ,IM)
     &                 + FREQ * FAC * BVEC(IVAR,IZ,IP)
                    BVEC(IVAR,IZ,IP) = BP
                    BVEC(IVAR,IZ,IM) = BM
                  ENDDO
                ENDDO
                BNRMP = BNRMP/DNRM2(NVARQ,BVEC(1,1,IP),1)
                BNRMM = BNRMM/DNRM2(NVARQ,BVEC(1,1,IM),1)
C               retain norm for numerical stability !
                CALL DSCAL(NVARQ,BNRMP,BVEC(1,1,IP),1)
                CALL DSCAL(NVARQ,BNRMM,BVEC(1,1,IM),1)
              ELSEIF(IVECS(IP).EQ.1) THEN
                BNRMP = DNRM2(NVARQ,BVEC(1,1,IP),1)
                DO IZ = 1,NZ
                  DO IVAR = 1,NVAR
                    FAC = HDIAG(IVAR)/(FREQ2-HDIAG(IVAR)**2)
                    BVEC(IVAR,IZ,IP) = FAC * BVEC(IVAR,IZ,IP)
Chj:    Here bvec(ivar,iz,im) should be negligible
                  ENDDO
                ENDDO
                BNRMP = BNRMP/DNRM2(NVARQ,BVEC(1,1,IP),1)
C               retain norm for numerical stability !
                CALL DSCAL(NVARQ,BNRMP,BVEC(1,1,IP),1)
              ELSEIF(IVECS(IM).EQ.1) THEN
                BNRMM = DNRM2(NVARQ,BVEC(1,1,IM),1)
                DO IZ = 1,NZ
                  DO IVAR = 1,NVAR
                    FAC = HDIAG(IVAR)/(FREQ2-HDIAG(IVAR)**2)
                    BVEC(IVAR,IZ,IM) = FAC * BVEC(IVAR,IZ,IM)
Chj:    Here bvec(ivar,iz,ip) should be negligible
                  ENDDO
                BNRMM = BNRMM/DNRM2(NVARQ,BVEC(1,1,IM),1)
C               retain norm for numerical stability !
                CALL DSCAL(NVARQ,BNRMM,BVEC(1,1,IM),1)
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO          
C
 10     CONTINUE
        CALL XRSORT(BVEC,JBTYP,NBNEW,NBPREV,
     &            IBVEC,IBTYP,IVECS,WORK,KFREE,LFREE)
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
