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
      CALL MEMGET('REAL',KFMO,NOSIM*N2BBASXQ,WORK,KFREE,LFREE)
C
C     1. Generate G-matrix
C     ====================
C
      IF(LSVCFG(1)) THEN
        CALL GMOLI(WORK(KFMO),BEVEC,BPVEC,IBEVC,IBPVC,IBTYP,
     &           NPOS,WORK,KFREE,LFREE)
      ELSE
        CALL DZERO(WORK(KFMO),NOSIM*N2BBASXQ)
      ENDIF
C
C     2) Add one-index transformed Fock-matrices
C     ==========================================
C
      IF (LSVCFG(2)) CALL FMOLI(WORK(KFMO),BEVEC,BPVEC,
     &                          WORK,KFREE,LFREE)
C
C     3) Gather contributions to electronic sigma-vectors
C     ===================================================
C
      IF(NZXOPE.GT.0) THEN
        CALL SMOLI(SEVEC,WORK(KFMO),WORK(KZXOPE),NZXOPE,NTSIM)
        IF(IPRXRS.GE.5) THEN
          WRITE(LUPRI,'(/I5,A)') NTSIM,
     &          ' orbital (e-e) sigma vectors'
          CALL PRBVEC(LUPRI,SEVEC,NTSIM,NZXOPE)
        ENDIF
      ENDIF
C
C     4) Gather contributions to positronic sigma-vectors
C     ===================================================
C
      IF(NZXOPP.GT.0) THEN
        CALL SMOLI(SPVEC,WORK(KFMO),WORK(KZXOPP),NZXOPP,NTSIM)
        IF(IPRXRS.GE.5) THEN
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
      CALL MEMGET('REAL',KFAO ,NDMAT*N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KSYM ,NOSIM        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KHRM ,NOSIM        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KFCK ,NOSIM        ,WORK,KFREE,LFREE)
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
#include <dummy.h>
C
#include <cbihr2.h>
#include <dcbxrs.h>
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
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
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','DMOLI2')
      CALL REACMO(LUBUF,'DFCOEF',CMO,DUMMY,TOTERG,3)
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
     &            FAO,FMO,NOSIM,NPOS,INTFLG,
     &            IPRTWO,WORK(KFREE),LFREE)
C
      IM = 0
C
C     c) Transform to MO-basis; e-e rotations
C     ---------------------------------------
C
      ITYP  = 1
      DO I = 1,NESIM
          IM = IM + 1
          CALL GMOLI3(ITYP,FMO(1,IM),CMO,FAO(1,IM),
     &         WORK,KFREE,LFREE)
C         Print section
          IF(IPRXRS.GE.8) THEN
            IREP = JSYMOP - 1
            WRITE(LUPRI,'(A,I5)') 
     &        '* GMOLI2: (e-e) G-matrix (AO-basis) no.:',I
            CALL PRQMAT(FAO(1,IM),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(A,I5)') 
     &        '* GMOLI2: (e-e) G-matrix (MO-basis) no.:',I
            CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
          ENDIF
      ENDDO
C
C     d) Transform to MO-basis; e-p rotations
C     ---------------------------------------
C
      ITYP  = 2
      DO I = 1,NPSIM
        IM = IM + 1
        CALL GMOLI3(ITYP,FMO(1,IM),CMO,FAO(1,IM),
     &         WORK,KFREE,LFREE)
C         Print section
          IF(IPRXRS.GE.8) THEN
            IREP = JSYMOP - 1
            WRITE(LUPRI,'(A,I5)') 
     &        '* GMOLI2: (e-p) G-matrix (AO-basis) no.:',I
            CALL PRQMAT(FAO(1,IM),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
            WRITE(LUPRI,'(A,I5)') 
     &        '* GMOLI2: (e-p) G-matrix (MO-basis) no.:',I
            CALL PRQMAT(FMO(1,IM),NORBT,NORBT,NORBT,NORBT,NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
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
C
C     Generate modified density matrices ; e-e rotations
C     ===================================================
C
      ITYP  = 1
      IBOFF = 1
      DO I = 1,NESIM
        IM = IM + 1
C       Scatter BVEC into BMAT
        CALL DZERO(BMAT,N2ORBXQ)
        CALL XOPSCT(BEVEC(IBOFF),WORK(KZXOPE),NZXOPE,BMAT,NORBT,NZ)
        IBOFF = IBOFF + NZXOPEQ
C       Construct modified density matrix
        CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT,
     &                WORK,KFREE,LFREE)
C       Print section
        IF(IPRXRS.GE.8) THEN
          IREP = JSYMOP - 1
          WRITE(LUPRI,'(A,I5)') 
     &      '* DMOLI: Scattered (e-e) trial vector no.:',I
          CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
          WRITE(LUPRI,'(A,I5)') 
     &      '* DMOLI: Modified (e-e) matrix no.:',I
          CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDDO
C
C     Generate modified density matrices ; e-p rotations
C     ===================================================
C
      ITYP  = 2
      IBOFF = 1
      DO I = 1,NPSIM
        IM = IM + 1
C       Scatter BVEC into BMAT
        CALL DZERO(BMAT,N2ORBXQ)
        CALL XOPSCT(BPVEC(IBOFF),WORK(KZXOPP),NZXOPP,BMAT,NORBT,NZ)
        IBOFF = IBOFF + NZXOPPQ
C       Construct modified density matrix
        CALL DMOLI3(ITYP,IHRMOP(IM),DMAT(1,IM),CMO,BMAT,
     &              WORK,KFREE,LFREE)
C       Print section
        IF(IPRXRS.GE.8) THEN
          IREP = JSYMOP - 1
          WRITE(LUPRI,'(A,I5)') 
     &      '* DMOLI: Scattered (e-p) trial vector no.:',I
          CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
          WRITE(LUPRI,'(A,I5)') 
     &      '* DMOLI: Modified (e-p) matrix no.:',I
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(DMAT(1,IM),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
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
C     ITYP = 1   e-e  rotations in BMAT
C     ITYP = 2   e-p  rotations in BMAT
C
C     We define: DB_{ab} = c_{as}B_{si}c_{bi}^*
C                a,b - AO-indices
C                s   - secondary index
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
      PARAMETER(D0 = 0.0D0)
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      DIMENSION DMAT(N2BBASX,NZ),CMO(*),BMAT(NORBT,NORBT,NZ),WORK(*)
      DIMENSION IQBF(4)
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
C       Parameters for secondary index
C
        IF    (ITYP.EQ.1) THEN
          NS = NSSH(I1)
          IS = IORB(I1) + NPSH(I1) + NOCC(I1) + 1
          JS = ICMOQ(I1) + (NPSH(I1)+NOCC(I1))*NFBAS(I1,0) + 1
        ELSEIF(ITYP.EQ.2) THEN
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ENDIF
C
C       Parameters for inactive index
C
        NI  = NISH(I2) - NFRO(I2)
        II  = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
        JI  = ICMOQ(I2) + (NPSH(I2)+ NFRO(I2))*NFBAS(I2,0) + 1
C
C       Generate modified density matrix as transformation MO to AO
C
        IF(NI.EQ.0.OR.NS.EQ.0) GOTO 10
        CALL QTRANS('MOAO','S',D0,NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &              DMAT(I2BASX(I1,I2)+1,1),NTBAS(0),NTBAS(0),NBF,IQBF,
     &              BMAT(IS,II,1),NORBT,NORBT,NZ,IPQTOQ(1,IREP),
     &              CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &              CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,IPRXRS)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** DFSOLV: DMOLI3.NE.NZ !'
 10   CONTINUE
      IF(IPRXRS.GE.5) THEN
        CALL HEADER('DMOLI3: Non-symmetrized density matrix',-1)
        CALL PRQMAT(DMAT,NTBAS(0),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF
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
      SUBROUTINE FMOLI(FMO,BEVEC,BPVEC,WORK,KFREE,LFREE)
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
      DIMENSION FMO(*),BEVEC(*),BPVEC(*),WORK(*)
C
      CALL QENTER('FMOLI')
      KFRSAV = KFREE
      CALL MEMGET('REAL',KBMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFBUF,N2BBASXQ,WORK,KFREE,LFREE)
C
C     Get total Fock matrix in MO-basis
C
      CALL GETFCK(WORK(KFBUF),IPRXRS,WORK,KFREE,LFREE)
C
C     Do one-index transformation
C
      CALL FMOLI2(FMO,WORK(KFBUF),WORK(KBMAT),BEVEC,BPVEC,
     &            WORK,KFREE,LFREE)
C
      CALL MEMREL('FMOLI',WORK,1,KFRSAV,KFREE,LFREE)      
      CALL QEXIT('FMOLI')      
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fmoli2 */
      SUBROUTINE FMOLI2(FMO,FBUF,BMAT,BEVEC,BPVEC,WORK,KFREE,LFREE)
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
      DIMENSION FMO(N2ORBXQ,*),FBUF(*),BMAT(*),BEVEC(*),BPVEC(*),
     &          WORK(*)
C
      IM = 0
C
C     Add one-index tranformed Fock-matrices; e-e rotations
C     =====================================================
C
      ITYP  = 1
      IBOFF = 1
      DO I = 1,NESIM
        IM = IM + 1
C       Scatter BVEC into BMAT
        CALL DZERO(BMAT,N2ORBXQ)
        CALL XOPSCT(BEVEC(IBOFF),WORK(KZXOPE),NZXOPE,BMAT,NORBT,NZ)
C       Print section
        IF(IPRXRS.GE.6) THEN
          IREP = JSYMOP - 1
          WRITE(LUPRI,'(A,I5)') 
     &      '* FMOLI: Scattered (e-e) trial vector no.:',I
          CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
        IBOFF = IBOFF + NZXOPEQ
C       Add one-index transformed Fock-matrix
        CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF,
     &                WORK,KFREE,LFREE)
      ENDDO
C
C     Add one-index tranformed Fock-matrices; e-p rotations
C     =====================================================
C
      ITYP  = 2
      IBOFF = 1
      DO I = 1,NPSIM
        IM = IM + 1
C       Scatter BVEC into BMAT
        CALL DZERO(BMAT,N2ORBXQ)
        CALL XOPSCT(BPVEC(IBOFF),WORK(KZXOPP),NZXOPP,BMAT,NORBT,NZ)
C       Print section
        IF(IPRXRS.GE.6) THEN
          IREP = JSYMOP - 1
          WRITE(LUPRI,'(A,I5)') 
     &      '* FMOLI: Scattered (e-p) trial vector no.:',I
          CALL PRQMAT(BMAT,NORBT,NORBT,NORBT,NORBT,NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
        IBOFF = IBOFF + NZXOPPQ
C       Add one-index transformed Fock-matrix
        CALL FMOLI3(ITYP,FMO(1,IM),BMAT,FBUF,
     &              WORK,KFREE,LFREE)
      ENDDO
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
C     ITYP = 1   e-e  rotations in BMAT
C     ITYP = 2   e-p  rotations in BMAT
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
      LOGICAL EPMIX
      DIMENSION FMO(NORBT,NORBT,NZ),FBUF(NORBT,NORBT,NZ),
     &          BMAT(NORBT,NORBT,NZ),WORK(*)
C
      EXTERNAL BDQALGEBRA
C
C     Initialization
C     ==============
C
      IF(ITYP.EQ.1)     THEN
        EPMIX = NZXOPP.GT.0
      ELSEIF(ITYP.EQ.2) THEN
        EPMIX = NZXOPE.GT.0
      ELSE
        WRITE(LUPRI,'(A,I5)') 
     &    '*** ERROR *** FMOLI3: Unknown ITYP :',ITYP
        CALL QUIT('FMOLI3: Unknown ITYP !')
      ENDIF
      IREP = JSYMOP-1
C
C     Generate one-index transformed Fock matrix
C     ==========================================
C
      DO 10 I1 = 1,NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for secondary index
C
        IF    (ITYP.EQ.1) THEN
          NS = NSSH(I1)
          IS = IORB(I1) + NPSH(I1) + NOCC(I1) + 1
        ELSEIF(ITYP.EQ.2) THEN
          NS = NPSH(I1)
          IS = IORB(I1) + 1
        ENDIF
        IF(NS.EQ.0) GOTO 10
C
C       Parameters for inactive index
C
        NI  = NISH(I2) - NFRO(I2)
        II  = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
C
        IF(NI.NE.0) THEN
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
            IF    (ITYP.EQ.1) THEN
              NT = NPSH(I1)
              IT = IORB(I1) + 1
            ELSEIF(ITYP.EQ.2) THEN
              NT = NSSH(I1)
              IT = IORB(I1) + NPSH(I1) + NOCC(I1) + 1
            ENDIF
            IF(NT.NE.0) CALL QGEMM(NT,NI,NS,DM1,
     &          'N','N',IPQTOQ(1,0),
     &          FBUF(IT,IS,1),NORBT,NORBT,NZ,
     &          'N','N',IPQTOQ(1,IREP),
     &          BMAT(IS,II,1),NORBT,NORBT,NZ,
     &          DM1,IPQTOQ(1,IREP),FMO(IT,II,1),NORBT,NORBT,NZ)
          ENDIF
        ENDIF
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
      LOGICAL EPMIX
      DIMENSION FMO(NORBT,NORBT,NZ),CMO(*),
     &          FAO(N2BBASXQ),WORK(*),IQBF(4)
C
C
C     Transform only blocks of right symmetry combination
C     ===================================================
C
      IF(ITYP.EQ.1)     THEN
        EPMIX = NZXOPP.GT.0
      ELSEIF(ITYP.EQ.2) THEN
        EPMIX = NZXOPE.GT.0
      ELSE
        WRITE(LUPRI,'(A,I5)') 
     &    '*** ERROR *** GMOLI: Unknown ITYP :',ITYP
        CALL QUIT('GMOLI3: Unknown ITYP !')
      ENDIF
      IREP  = JSYMOP-1
      ITIM  = JTIMOP
      DO 10 I1 = 1,NFSYM
        I2  = MOD(I1+JOPSY,2) + 1
C
C       Parameters for secondary index
C
        IF    (ITYP.EQ.1) THEN
          NS = NSSH(I1)
          IS = IORB(I1) + NPSH(I1) + NOCC(I1) + 1
          JS = ICMOQ(I1) + (NPSH(I1)+NOCC(I1))*NFBAS(I1,0) + 1
        ELSEIF(ITYP.EQ.2) THEN
          NS = NPSH(I1)
          IS = IORB(I1) + 1
          JS = ICMOQ(I1) + 1
        ENDIF
C
C       Parameters for inactive index
C
        NI  = NISH(I2)-NFRO(I2)
        II  = IORB(I2) + NPSH(I2) + NFRO(I2) + 1
        JI  = ICMOQ(I2) + (NPSH(I2)+NFRO(I2))*NFBAS(I2,0) + 1
C
        IF(NS.NE.0.AND.NI.NE.0) THEN
          CALL QTRANS('AOMO','S',D0,NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),
     &                FMO(IS,II,1),NORBT,NORBT,NBF,IQBF,
     &                CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRXRS)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** GMOLI3: NBF.NE.NZ !'
        ENDIF
C
C       Transform other block as well, 
C       if both electronic and postitronic rotations exist
C     
        IF(EPMIX) THEN
C
C         Parameters for secondary index(note the switch...)
C
          IF    (ITYP.EQ.1) THEN
            NS = NPSH(I1)
            IS = IORB(I1) + 1
            JS = ICMOQ(I1) + 1
          ELSEIF(ITYP.EQ.2) THEN
            NS = NSSH(I1)
            IS = IORB(I1) + NPSH(I1) + NOCC(I1) + 1
            JS = ICMOQ(I1) + (NPSH(I1)+NOCC(I1))*NFBAS(I1,0) + 1
          ENDIF
          IF(NS.NE.0.AND.NI.NE.0) THEN
            CALL QTRANS('AOMO','S',D0,NFBAS(I1,0),NFBAS(I2,0),NS,NI,
     &                FAO(I2BASX(I1,I2)+1),NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),
     &                FMO(IS,II,1),NORBT,NORBT,NBF,IQBF,
     &                CMO(JS),NFBAS(I1,0),NORB(I1),NZ,IPQTOQ(1,0),
     &                CMO(JI),NFBAS(I2,0),NORB(I2),NZ,IPQTOQ(1,0),
     &                WORK(KFREE),LFREE,IPRXRS)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** GMOLI3: NBF.NE.NZ !'
          ENDIF
        ENDIF
 10   CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Smoli */
      SUBROUTINE SMOLI(SVEC,FMO,JXOP,NXOP,NFMAT)
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
      DIMENSION SVEC(NXOP,NZ,*),FMO(N2ORBXQ,*),JXOP(*)
      NXOPQ = NXOP*NZ
CTROND We assume first all NESIM and then all NPSIM.
      DO I = 1,NFMAT
        CALL XOPGAT(SVEC(1,1,I),JXOP,NXOP,FMO(1,I),NORBT,NZ)
        CALL DSCAL(NXOPQ,DM1,SVEC(1,1,I),1)
      ENDDO
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
C     e-e rotations
C     =============
C      
      DO I = 1,NESIM
        IS = IBEVC(NERED+I)
        IHRMOP(IM+I) = IBTYP(2,IS)
        IFCKOP(IM+I) = 1
        ISYMOP(IM+I) = JSYMOP
      ENDDO
C
C     e-p rotations
C     =============
C
      IM = NESIM
      DO I = 1,NPSIM
        IS = IBPVC(NPRED+I)
        IHRMOP(IM+I) = IBTYP(2,IS)
        IFCKOP(IM+I) = 1
        ISYMOP(IM+I) = JSYMOP
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xrstvc */
      SUBROUTINE XRSTVC(BEVEC,IBEVC,NBEVC,IEOFF,
     &                  BPVEC,IBPVC,NBPVC,IPOFF,
     &                  BCVEC,IBCVC,NBCVC,ICOFF,
     &                  TEVEC,TPVEC,TCVEC,NBTOT,ITOFF,
     &                  WORK,KFREE,LFREE)
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
     &                TEVEC,TPVEC,NBTOT,ITOFF,WORK,KFREE,LFREE)
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
     &                  TEVEC,TPVEC,NBTOT,ITOFF,WORK,KFREE,LFREE)
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
      PARAMETER(D1 = 1.0D0)
      DIMENSION BEVEC(NZXOPEQ,*),IBEVC(*),
     &          BPVEC(NZXOPPQ,*),IBPVC(*),
     &          BCVEC(NZCONFQ,*),IBCVC(*),
     &          TEVEC(NZXOPEQ,*),TPVEC(NZXOPPQ,*),
     &          WORK(*)
C
#include <dcborb.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbxrs.h>
C
      IF(NASHT.EQ.0) THEN
        DO I = 1,NBEVC
          J  = IBEVC(I+IEOFF) - ITOFF
          CALL DAXPY(NZXOPEQ,D1,BEVEC(1,I),1,TEVEC(1,J),1)
        ENDDO
        DO I = 1,NBPVC
          J  = IBPVC(I+IPOFF) - ITOFF
          CALL DAXPY(NZXOPPQ,D1,BPVEC(1,I),1,TPVEC(1,J),1)
        ENDDO
      ELSE
        WRITE(LUPRI,'(A)') 'XRSTOV: Open shell code not written !'
        CALL QUIT('XRSTEP: Open shell code not written !')
      ENDIF
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
     &              IBTYP,IBEVC,BEBUF)
      ENDIF
C
C     e-p rotations
C      
      IF(NPRED.GT.0) THEN
        CALL XRSXV1(ITYP,JBPNDX,XPVEC,EVECR,MFREQ,
     &              IBTYP,IBPVC,BPBUF)
      ENDIF
C
C     Configurational parameters
C      
      IF(NCRED.GT.0) THEN
        CALL XRSXV1(ITYP,JBCNDX,XCVEC,EVECR,MFREQ,
     &              IBTYP,IBCVC,BCBUF)
      ENDIF
C
      CALL QEXIT('XRSXVC')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Xrsxv1 */
      SUBROUTINE XRSXV1(ITYP,NBTYP,XVEC,EVECR,MFREQ,IBTYP,IBVEC,BVBUF)
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
      DIMENSION XVEC(*),EVECR(NREDM,MFREQ),
     &               IBTYP(2,*),IBVEC(*),BVBUF(*)
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
      CALL DZERO(XVEC,NBDIM*MFREQ)
      DO IREC = 1,NBRED
        IND = IBVEC(IREC)
        IH  = IBTYP(2,IND)
        IF(IH.EQ.ITYP) THEN
          CALL READAC(LUBVEC,NBDIM,BVBUF,IREC)
          IBOFF = 1
          DO IFREQ = 1,MFREQ
            CALL DAXPY(NBDIM,EVECR(IND,IFREQ),
     &                   BVBUF,1,XVEC(IBOFF),1)
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
     &                WORK,KFREE,LFREE)
C
C         Accumulate T(+) = S[2}X(-) in b(+)
C         ------------------------------------
C
          CALL XRSTVC(BEMVC(1,I),1,NBEVC,0,
     &                BPMVC(1,I),1,NBPVC,0,
     &                BCMVC(1,I),1,NBCVC,0,
     &                BEPVC(1,I),BPPVC(1,I),BCPVC(1,I),1,0,
     &                WORK,KFREE,LFREE)
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
     &                  BEBUF,BPBUF,BCBUF,RCNV,XECNV,XPCNV,XCCNV,
     &                  IEVEC,IPVEC,ICVEC,RNORM,WORK,KFREE,LFREE)
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
     &          BEBUF(*),BPBUF(*),BCBUF(*),
     &          XECNV(*),XPCNV(*),XCCNV(*),RCNV(*),
     &          IEVEC(*),IPVEC(*),ICVEC(*),RNORM(7,MFREQ),WORK(*)
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
     &                  IBTYP,IBEVC,BEBUF)
            CALL XVCNRM(XECNV,BEVEC,NZXOPEQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-e) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XECNV(I),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-e) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BEVEC,MFREQ,NZXOPE)
            ENDIF
          ELSE
            CALL DZERO(XECNV,NESIM)
          ENDIF
        ELSE
          NESIM = 0
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
     &                  IBTYP,IBPVC,BPBUF)
            CALL XVCNRM(XPCNV,BPVEC,NZXOPPQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-p) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XPCNV(I),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' orbital (e-p) X(+) solution vectors'
              CALL PRBVEC(LUPRI,BPVEC,MFREQ,NZXOPP)
            ENDIF
          ELSE
            CALL DZERO(XPCNV,NPSIM)
          ENDIF
        ELSE
          NPSIM = 0
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
     &                  IBTYP,IBCVC,BCBUF)
            CALL XVCNRM(XCCNV,BCVEC,NZCONFQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of configurational X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XCCNV(I),I=1,MFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/I5,A)') MFREQ, 
     &          ' configurational X(+) solution vectors'
              CALL PRBVEC(LUPRI,BCVEC,MFREQ,NZCONF)
            ENDIF
          ELSE
            CALL DZERO(XPCNV,NCSIM)
          ENDIF
        ELSE
          NCSIM = 0
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
     &                  IBTYP,IBEVC,BEBUF)
            CALL XVCNRM(XECNV      ,BEVEC      ,NZXOPEQ,MFREQ)
            CALL XRSXV1(IMX,JBENDX,BEVEC(KBE2),EVECR,MFREQ,
     &                  IBTYP,IBEVC,BEBUF)
            CALL XVCNRM(XECNV(IXE2+1),BEVEC(KBE2),NZXOPEQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-e) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XECNV(I),I=1,MFREQ)
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-e) X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XECNV(IXE2+I),I=1,MFREQ)
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
            CALL DZERO(XECNV,NESIM)
          ENDIF
        ELSE
          NESIM = 0
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
     &                  IBTYP,IBPVC,BPBUF)
            CALL XVCNRM(XPCNV      ,BPVEC      ,NZXOPPQ,MFREQ)
            CALL XRSXV1(IMX,JBPNDX,BPVEC(KBP2),EVECR,MFREQ,
     &                  IBTYP,IBPVC,BPBUF)
            CALL XVCNRM(XPCNV(IXP2+1),BPVEC(KBP2),NZXOPPQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-p) X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XPCNV(I),I=1,MFREQ)
              WRITE(LUPRI,'(A)') 
     &          'Norms of orbital(e-p) X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XPCNV(IXP2+I),I=1,MFREQ)
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
            CALL DZERO(XPCNV,NPSIM)
          ENDIF
        ELSE
          NPSIM = 0
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
     &                  IBTYP,IBCVC,BCBUF)
            CALL XVCNRM(XCCNV      ,BCVEC      ,NZCONFQ,MFREQ)
            CALL XRSXV1(IMX,JBCNDX,BCVEC      ,EVECR,MFREQ,
     &                  IBTYP,IBCVC,BCBUF)
            CALL XVCNRM(XCCNV(IXC2+1),BCVEC(KBC2),NZCONFQ,MFREQ)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A)') 
     &          'Norms of configurational X(+)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XCCNV(I),I=1,MFREQ)
              WRITE(LUPRI,'(A)') 
     &          'Norms of configurational X(-)-solution vectors:'
              WRITE(LUPRI,'(1P,8(3X,E8.2))') (XCCNV(IXC2+I),I=1,MFREQ)
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
            CALL DZERO(XPCNV,NCSIM)
          ENDIF
        ELSE
          NCSIM = 0
        ENDIF
      ENDIF
      NTSIM = NESIM + NPSIM + NCSIM
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
C     Add E[2}X
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
      CALL XRSBVC(BEVEC,BPVEC,BCVEC,XECNV,XPCNV,XCCNV,
     &            MFREQ,IBTYP,IBCVC,IBEVC,IBPVC,IEVEC,IPVEC,ICVEC,
     &            RCNV,EVALR,EVECR,RNORM,BEBUF,BPBUF,BCBUF,
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
      SUBROUTINE XRSBVC(BEVEC,BPVEC,BCVEC,XECNV,XPCNV,XCCNV,
     &            MFREQ,IBTYP,IBCVC,IBEVC,IBPVC,IEVEC,IPVEC,ICVEC,
     &            RCNV,EVALR,EVECR,RNORM,BEBUF,BPBUF,BCBUF,
     &            WORK,KFREE,LFREE)
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
C             1: norm of E(+)-part
C             2: norm of P(+)-part
C             3: norm of C(+)-part
C             4: norm of E(-)-part
C             5: norm of P(-)-part
C             6: norm of C(-)-part
C
C     Written by T.Saue, Oct 18 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include <dcbxrs.h>
#include <dcbibn.h>
      LOGICAL   USEXVC
      DIMENSION BEVEC(NZXOPEQ,*),BPVEC(NZXOPPQ,*),
     &          BCVEC(NZCONFQ,*),EVALR(*),EVECR(*),
     &          XECNV(*),XPCNV(*),XCCNV(*),RCNV(*),
     &          IBTYP(2,NREDM),IBEVC(NREDM),IBPVC(NREDM),IBCVC(NREDM),
     &          IEVEC(*),IPVEC(*),ICVEC(*),RNORM(7,MFREQ),
     &          BEBUF(*),BPBUF(*),BCBUF(*),WORK(*)
C
      IF(IPRXRS.GE.3) CALL HEADER('Output from XRSBVC',-1)
C
C     ********************************
C     ***** Select trial vectors *****
C     ********************************
C
C
      ISIM   = 0
      ISIM2  = 0
      CNVMAX = D0
C
C     ===================
C     *** Static code ***
C     ===================
C
      IF(STATIC) THEN
        DO IFREQ = 1,MFREQ
          ISIM  = ISIM + 1
C
C         Check norm of each residual
C         ===========================
C
C         Orbital (e-e) part
C         ------------------
C
          IF(NZXOPE.GT.0) THEN
            IEVEC(ISIM) = 1
            RNORM(2,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,ISIM),1)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A,F8.5,A)') 
     &          'Norm of orbital(e-e) X(+)-residual vector (freq = ',
     &          EVALR(IFREQ),')'
              WRITE(LUPRI,'(1P,3X,E8.2)') RNORM(2,IFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/A)')
     &          ' Orbital (e-e) X(+) residual vector:'
              CALL PRBVEC(LUPRI,BEVEC(1,ISIM),1,NZXOPE)
            ENDIF
          ELSE
            IEVEC(ISIM) = 0
          ENDIF
C
C         Orbital (e-p) part
C         ------------------
C
          IF(NZXOPP.GT.0) THEN
            IPVEC(ISIM) = 1
            RNORM(3,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,ISIM),1)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A,F8.5,A)') 
     &          'Norm of orbital(e-p) X(+)-residual vector (freq = ',
     &          EVALR(IFREQ),')'
              WRITE(LUPRI,'(1P,3X,E8.2)') RNORM(3,IFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/A)') 
     &          'Orbital (e-e) X(+) residual vector:'
              CALL PRBVEC(LUPRI,BPVEC(1,ISIM),1,NZXOPP)
            ENDIF
          ELSE
            IPVEC(ISIM) = 0
          ENDIF
C
C         Configurational part
C         --------------------
C
          IF(NZCONF.GT.0) THEN
            ICVEC(ISIM) = 1
            RNORM(4,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,ISIM),1)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(A,F8.5,A)') 
     &          'Norm of configurational X(+)-residual vector (freq = ',
     &          EVALR(IFREQ),')'
              WRITE(LUPRI,'(1P,3X,E8.2)') RNORM(4,IFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/A)') 
     &          'Configurational X(+) residual vector:'
              CALL PRBVEC(LUPRI,BCVEC(1,ISIM),1,NZCONF)
            ENDIF
          ELSE
            ICVEC(ISIM) = 0
          ENDIF
          RNORM(1,IFREQ) = SQRT(RNORM(2,IFREQ)**2 + RNORM(3,IFREQ)**2 
     &                        + RNORM(4,IFREQ)**2)
C
          XTNORM      = SQRT(XECNV(ISIM)**2    + XPCNV(ISIM)**2 
     &                        + XCCNV(ISIM)**2)
          IF(XTNORM.EQ.D0) THEN
            WRITE(LUPRI,'(A,A)') '*** WARNING! ***',
     &        'Zero norm of soution vector !'
          ELSE
            RCNV(IFREQ) = RNORM(1,IFREQ)/XTNORM
            CNVMAX = MAX(CNVMAX,RCNV(IFREQ))
          ENDIF
C
C         Convergence check
C
          IF(RNORM(1,IFREQ).LE.(XTNORM*THCXRS)) THEN
C         ... this root converged
            IEVEC(ISIM)  = 0
            IPVEC(ISIM)  = 0
            ICVEC(ISIM)  = 0
          ELSE
C         ... check each type separately
            RNRMAX = MAX(RNORM(2,IFREQ),RNORM(3,IFREQ),RNORM(4,IFREQ))
            IF(RESFAC*RNORM(2,IFREQ).LT.RNRMAX) IEVEC(ISIM) = 0
            IF(RESFAC*RNORM(3,IFREQ).LT.RNRMAX) IPVEC(ISIM) = 0
            IF(RESFAC*RNORM(4,IFREQ).LT.RNRMAX) ICVEC(ISIM) = 0
          ENDIF
          IF(IPRXRS.GE.1) THEN
            WRITE(LUPRI,'(A,A,F15.8/)') 'XRSBVC: Convergence check',
     &        ' at frequency:',EVALR(ISIM)
            WRITE(LUPRI,'(3X,A,1P,E8.2)') 
     &         '* Norm of total resultant vector :',XTNORM
            WRITE(LUPRI,'(3X,A,1P,E8.2)') 
     &         '* Norm of total residual vector  :',RNORM(1,IFREQ)
            IF(NZXOPE.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &         '- (e-e) R(+) - part           :',RNORM(2,IFREQ)
            IF(NZXOPP.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &         '- (e-p) R(+) - part           :',RNORM(3,IFREQ)
            IF(NZCONF.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &         '- conf. R(+) - part           :',RNORM(4,IFREQ)
          ENDIF
        ENDDO
        IF(IPRXRS.GE.3) THEN
          CALL PRSYMB(LUPRI,'-',48,3)
          WRITE(LUPRI,'(3X,A,3X,3(A3,I4,1X))')
     &      'Static convergence:', 
     &      'E+:',IEVEC(1),'P+:',IPVEC(1),'C+:',ICVEC(1)
        ENDIF
      ELSE
C
C     ================================
C     *** Frequency-dependent code ***
C     ================================
C
        DO IFREQ = 1,MFREQ
          ISIM  = ISIM + 1
          ISIM2 = ISIM + MFREQ
C
C         Check norm of each residual
C         ===========================
C
C
C         Orbital (e-e) part
C         ------------------
C
          IF(NZXOPE.GT.0) THEN
            IEVEC(ISIM)  = 1
            IEVEC(ISIM2) = 1
            RNORM(2,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,ISIM) ,1)
            RNORM(5,IFREQ) = DNRM2(NZXOPEQ,BEVEC(1,ISIM2),1)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(/3X,A,3X,F15.8)') 
     &          '*Frequency:',EVALR(IFREQ)
              WRITE(LUPRI,'(A,1P,E8.2)') 
     &          'Norm of orbital(e-e) X(+)-residual vector :',
     &          RNORM(2,IFREQ)
              WRITE(LUPRI,'(A,1P,E8.2)') 
     &          'Norm of orbital(e-e) X(-)-residual vector :',
     &          RNORM(5,IFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/A)')
     &          ' Orbital (e-e) X(+) residual vector:'
              CALL PRBVEC(LUPRI,BEVEC(1,ISIM),1,NZXOPE)
              WRITE(LUPRI,'(/A)')
     &          ' Orbital (e-e) X(-) residual vector:'
              CALL PRBVEC(LUPRI,BEVEC(1,ISIM2),1,NZXOPE)
            ENDIF
          ELSE
            IEVEC(ISIM)  = 0
            IEVEC(ISIM2) = 0
            RNORM(2,IFREQ) = D0
            RNORM(5,IFREQ) = D0
          ENDIF
C
          IF(NZXOPP.GT.0) THEN
            IPVEC(ISIM)  = 1
            IPVEC(ISIM2) = 1
            RNORM(3,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,ISIM) ,1)
            RNORM(6,IFREQ) = DNRM2(NZXOPPQ,BPVEC(1,ISIM2),1)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(/3X,A,3X,F15.8)') 
     &          '*Frequency:',EVALR(IFREQ)
              WRITE(LUPRI,'(A,1P,E8.2)') 
     &          'Norm of orbital(e-p) X(+)-residual vector :',
     &          RNORM(3,IFREQ)
              WRITE(LUPRI,'(A,1P,E8.2)') 
     &          'Norm of orbital(e-p) X(-)-residual vector :',
     &          RNORM(6,IFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/A)')
     &          ' Orbital (e-p) X(+) residual vector:'
              CALL PRBVEC(LUPRI,BPVEC(1,ISIM),1,NZXOPP)
              WRITE(LUPRI,'(/A)')
     &          ' Orbital (e-p) X(-) residual vector:'
              CALL PRBVEC(LUPRI,BPVEC(1,ISIM2),1,NZXOPP)
            ENDIF
          ELSE
            IPVEC(ISIM)  = 0
            IPVEC(ISIM2) = 0
            RNORM(3,IFREQ) = D0
            RNORM(6,IFREQ) = D0
          ENDIF
C
          IF(NZCONF.GT.0) THEN
            ICVEC(ISIM)  = 1
            ICVEC(ISIM2) = 1
            RNORM(4,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,ISIM) ,1)
            RNORM(7,IFREQ) = DNRM2(NZCONFQ,BCVEC(1,ISIM2),1)
C           Print section
            IF(IPRXRS.GE.3) THEN
              WRITE(LUPRI,'(/3X,A,3X,F15.8)') 
     &          '*Frequency:',EVALR(IFREQ)
              WRITE(LUPRI,'(A,1P,E8.2)') 
     &          'Norm of configurational X(+)-residual vector :',
     &          RNORM(4,IFREQ)
              WRITE(LUPRI,'(A,1P,E8.2)') 
     &          'Norm of configurational X(-)-residual vector :',
     &          RNORM(7,IFREQ)
            ENDIF
            IF(IPRXRS.GE.5) THEN
              WRITE(LUPRI,'(/A)')
     &          ' Configurational X(+) residual vector:'
              CALL PRBVEC(LUPRI,BCVEC(1,ISIM),1,NZCONF)
              WRITE(LUPRI,'(/A)')
     &          ' Configurational X(-) residual vector:'
              CALL PRBVEC(LUPRI,BCVEC(1,ISIM2),1,NZCONF)
            ENDIF
          ELSE
            ICVEC(ISIM)  = 0
            ICVEC(ISIM2) = 0
            RNORM(4,IFREQ) = D0
            RNORM(7,IFREQ) = D0
          ENDIF
C
          RTPNRM = SQRT(RNORM(2,IFREQ)**2 + RNORM(3,IFREQ)**2 
     &                + RNORM(4,IFREQ)**2)
          RTMNRM = SQRT(RNORM(5,IFREQ)**2 + RNORM(6,IFREQ)**2 
     &                + RNORM(7,IFREQ)**2)
C
          XTPNRM = SQRT(XECNV(ISIM)**2    + XPCNV(ISIM)**2  
     &                + XCCNV(ISIM)**2)
          XTMNRM = SQRT(XECNV(ISIM2)**2   + XPCNV(ISIM2)**2 
     &                + XCCNV(ISIM2)**2)
C
C         Convergence check
C         =================
C
          RNORM(1,IFREQ) = SQRT(RTPNRM**2+RTMNRM**2)
          XTNORM = SQRT(XTPNRM**2+XTMNRM**2)
          IF(XTNORM.EQ.D0) THEN
            WRITE(LUPRI,'(A,A,F9.6)') '*** WARNING! ***',
     &        'Zero norm of solution vector for frequency',EVALR(IFREQ)
          ELSE
            RCNV(IFREQ) = RNORM(1,IFREQ)/XTNORM
            CNVMAX = MAX(CNVMAX,RCNV(IFREQ))
          ENDIF
          IF(RNORM(1,IFREQ).LE.(XTNORM*THCXRS)) THEN
C          IF(RCNV(IFREQ).LE.THCXRS) THEN
C         ... this root converged
            IEVEC(ISIM)  = 0
            IPVEC(ISIM)  = 0
            ICVEC(ISIM)  = 0
            IEVEC(ISIM2) = 0
            IPVEC(ISIM2) = 0
            ICVEC(ISIM2) = 0
          ELSE
C         ... check each type separately
C           Check (+) convergence
            IF(RTPNRM.LE.(XTPNRM*THCXRS)) THEN
              IEVEC(ISIM) = 0
              IPVEC(ISIM) = 0
              ICVEC(ISIM) = 0
            ELSE
              RNPMAX=MAX(RNORM(2,IFREQ),RNORM(3,IFREQ),RNORM(4,IFREQ))
              IF(RESFAC*RNORM(2,IFREQ).LT.RNPMAX) IEVEC(ISIM) = 0
              IF(RESFAC*RNORM(3,IFREQ).LT.RNPMAX) IPVEC(ISIM) = 0
              IF(RESFAC*RNORM(4,IFREQ).LT.RNPMAX) ICVEC(ISIM) = 0
            ENDIF
C           Check (-) convergence
            IF(RTMNRM.LE.(XTMNRM*THCXRS)) THEN
              IEVEC(ISIM2) = 0
              IPVEC(ISIM2) = 0
              ICVEC(ISIM2) = 0
            ELSE
              RNMMAX=MAX(RNORM(5,IFREQ),RNORM(6,IFREQ),RNORM(7,IFREQ))
              IF(RESFAC*RNORM(5,IFREQ).LT.RNMMAX) IEVEC(ISIM2) = 0
              IF(RESFAC*RNORM(6,IFREQ).LT.RNMMAX) IPVEC(ISIM2) = 0
              IF(RESFAC*RNORM(7,IFREQ).LT.RNMMAX) ICVEC(ISIM2) = 0
            ENDIF
          ENDIF
C
C     Print section
C
          IF(IPRXRS.GE.1) THEN
            WRITE(LUPRI,'(A,A,F15.8/)') 'XRSBVC: Convergence check',
     &        ' at frequency:',EVALR(IFREQ)
            WRITE(LUPRI,'(3X,A,1P,E8.2)') 
     &         '* Norm of total resultant vector :',XTNORM
            WRITE(LUPRI,'(3X,A,1P,E8.2)') 
     &         '* Norm of total residual vector  :',RNORM(1,IFREQ)
            IF(NZXOPE.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &         '- (e-e) R(+) - part           :',RNORM(2,IFREQ)
            IF(NZXOPP.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &         '- (e-p) R(+) - part           :',RNORM(3,IFREQ)
            IF(NZCONF.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &         '- conf. R(+) - part           :',RNORM(4,IFREQ)
            IF(NZXOPE.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &        '- (e-e) R(-) - part           :',RNORM(5,IFREQ)
            IF(NZXOPP.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &        '- (e-p) R(-) - part           :',RNORM(6,IFREQ)
            IF(NZCONF.GT.0) WRITE(LUPRI,'(6X,A,1P,E8.2)') 
     &        '- conf. R(-) - part           :',RNORM(7,IFREQ)
          ENDIF
        ENDDO
        IF(IPRXRS.GE.3) THEN
          IOFF = MFREQ
          WRITE(LUPRI,'(/3X,A)')
     &      'Convergence of individual frequencies:'
          CALL PRSYMB(LUPRI,'-',48,3)
          WRITE(LUPRI,'(30X,6(A2,2X))')
     &      'E+','P+','C+','E-','P-','C-'
          WRITE(LUPRI,'(3X,A,F15.8,3X,6I4)')
     &      ('Freq.:',EVALR(I),
     &      IEVEC(I),IPVEC(I),ICVEC(I),
     &      IEVEC(IOFF+I),IPVEC(IOFF+I),ICVEC(IOFF+I),
     &       I = 1,MFREQ)
        ENDIF
      ENDIF
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
        NECNV = ISUM(NESIM,IEVEC,1)
        IF(NECNV.EQ.0) NESIM = 0
        NPCNV = ISUM(NPSIM,IPVEC,1)
        IF(NPCNV.EQ.0) NPSIM = 0
        NCCNV = ISUM(NCSIM,ICVEC,1)
        IF(NCCNV.EQ.0) NCSIM = 0
        KCONV = NECNV + NPCNV + NCCNV
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
        IF(NZXOPE.GT.0) 
     &    CALL MEMGET('REAL',KEDIAH,NZXOPE,WORK,KFREE,LFREE)
        IF(NZXOPP.GT.0) 
     &    CALL MEMGET('REAL',KPDIAH,NZXOPP,WORK,KFREE,LFREE)
        IF(NZCONF.GT.0) 
     &    CALL MEMGET('REAL',KCDIAH,NZCONF,WORK,KFREE,LFREE)
        CALL E2DIAG(WORK(KEDIAH),WORK(KPDIAH),WORK(KCDIAH),
     &              WORK,KFREE,LFREE)
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
      IF(NESIM.GT.0) THEN
        CALL XBVGEN(BEVEC,USEXVC,WORK(KEDIAH),EVALR,EVECR,
     &              IBTYP,IBEVC,IEVEC,
     &              MFREQ,JBENDX,NESIM,IBOFF,NEPREV,NZXOPE,
     &              BEBUF,WORK,KFREE,LFREE)
        IBOFF = IBOFF + NESIM
      ENDIF
C
C     Orbital (e-p) type trial vectors
C
      IF(NPSIM.GT.0) THEN
        CALL XBVGEN(BPVEC,USEXVC,WORK(KPDIAH),EVALR,EVECR,
     &              IBTYP,IBPVC,IPVEC,
     &              MFREQ,JBPNDX,NPSIM,IBOFF,NPPREV,NZXOPP,
     &              BPBUF,WORK,KFREE,LFREE)
        IBOFF = IBOFF + NPSIM
      ENDIF
C
C     Configuration type trial vectors
C
      IF(NCSIM.GT.0) THEN
        CALL XBVGEN(BCVEC,USEXVC,WORK(KCDIAH),EVALR,EVECR,
     &              IBTYP,IBCVC,ICVEC,
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
C     Zero reduced space if USEXVC
C     ============================
C   
      IF(USEXVC) THEN
        NERED = 0
        NPRED = 0
        NCRED = 0
        NZRED = 0
      ENDIF
C
C     Check dimension of reduced space
C     ================================
C
      MREDM = NZRED + NTSIM
      IF(MREDM.GT.NREDM) THEN
        KCONV = -2
        RETURN
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck XBVGEN */
      SUBROUTINE XBVGEN(BVEC,USEXVC,HDIAG,EVALR,EVECR,IBTYP,IBVEC,IVEC,
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
     &          EVALR(MFREQ),IBTYP(2,*),IBVEC(*),IVEC(*),
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
          CALL XRSXV1(IPX,JBTYP,BVEC,EVECR,MFREQ,
     &                    IBTYP,IBVEC,BBUF)
          IF(.NOT.STATIC) THEN
            CALL XRSXV1(IMX,JBTYP,BEVEC,EVECR,MFREQ,
     &                    IBTYP,IBVEC,BBUF)
          ENDIF
          CALL ICOPY(NBNEW,1,0,IVEC,1)
        ENDIF
        DO IP = 1,MFREQ
          IM = IP + MFREQ
          FREQ = EVALR(IP)
          IF(FREQ.EQ.D0) THEN
C           X(+) - vector
            IBTYP(1,IBOFF+IP) = JBTYP
            IBTYP(2,IBOFF+IP) = IPX
            IBVEC(NBPREV +IP) = IBOFF + IP
            IF(TEST.AND.IVEC(IP).EQ.1) THEN
              DO IZ = 1,NZ
                DO IVAR = 1,NVAR
                  BVEC(IVAR,IZ,IP) = 
     &              (DM1/HDIAG(IVAR))*BVEC(IVAR,IZ,IP)
                ENDDO
              ENDDO
            ENDIF
C           No X(-) - vector              
            IF(.NOT.STATIC) THEN
              IBTYP(1,IBOFF+IM) = JBTYP
              IBTYP(2,IBOFF+IM) = IMX
              IBVEC(NBPREV+IM)  = IBOFF + IM
              IVEC(IM)          = 0
            ENDIF
          ELSE
            IBTYP(1,IBOFF+IP)     = JBTYP
            IBTYP(2,IBOFF+IP)     = IPX
            IBVEC(NBPREV+IP)      = IBOFF + IP
            IBTYP(1,IBOFF+IM)     = JBTYP
            IBTYP(2,IBOFF+IM)     = IMX
            IBVEC(NBPREV+IM)      = IBOFF + IM
            IF(TEST) THEN
              IF(IVEC(IP).EQ.1.AND.IVEC(IM).EQ.1) THEN
                DO IZ = 1,NZ
                  DO IVAR = 1,NVAR
                    FAC = (FREQ**2-HDIAG(IVAR)**2)
                    BVEC(IVAR,IZ,IP) = (HDIAG(IVAR)/FAC)
     &                             *BVEC(IVAR,IZ,IP)
                    BVEC(IVAR,IZ,IM) = (FREQ/FAC)*BVEC(IVAR,IZ,IM)
                  ENDDO
                ENDDO
              ELSEIF(IVEC(IP).EQ.1) THEN
                DO IZ = 1,NZ
                  DO IVAR = 1,NVAR
                    FAC = (FREQ**2-HDIAG(IVAR)**2)
                    BVEC(IVAR,IZ,IP) = (HDIAG(IVAR)/FAC)
     &                               *BVEC(IVAR,IZ,IP)
                  ENDDO
                ENDDO
              ELSEIF(IVEC(IM).EQ.1) THEN
                DO IZ = 1,NZ
                  DO IVAR = 1,NVAR
                    FAC = (FREQ**2-HDIAG(IVAR)**2)
                    BVEC(IVAR,IZ,IM) = (FREQ/FAC)*BVEC(IVAR,IZ,IM)
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDIF
        ENDDO          
C
        CALL XRSORT(BVEC,JBTYP,NBNEW,NBPREV,
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
