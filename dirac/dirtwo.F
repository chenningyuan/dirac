C 
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
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
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C 
C
C
C VERSION : $Revision: 1.42.4.2 $
C DATE    : $Date: 2004/09/01 08:08:56 $
C FILE    : dirtwo.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twofck */
      SUBROUTINE TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMAT,DMAT,NFMAT,
     &                  NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for the construction of the two-electron Fock matrix in SO-basis
C
C     Written by T.Saue - January 1995
C     tsaue - Sep 3 1996    Extended to include parameter NFMAT, number of
C                           Fock matrices
C     
C     The reduced density matrices are generated within this routine,
C     whereas Cauchy-Schwarz integrals must be provided externally
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "dummy.h"
      PARAMETER(D0=0.0D0, D1=1.0D0)
C
#include "dcbind.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "dcbbas.h"
#include "ccom.h"
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),
     &          DMAT(N2BBASX,NZ,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),
     &          IFCKOP(NFMAT),NPOS(*),WORK(LWORK)
C
      CALL QENTER('TWOFCK')
#include "memint.h"
C
      IF(IPRINT.GE.5) THEN
        CALL TITLER('Output from TWOFCK','*',103)
      ENDIF
C 
C     In approximate density ZORA only the LL part of the density matrix
C     is to be used
C
      IF (ZORA.AND..NOT.ZORA4) THEN
         DO I = 1,NFMAT
            CALL KEEPLL(DMAT(1,1,I),NZ,'SRT')
         ENDDO
      ENDIF
C
C     Transform density matrix to unsorted basis
C     and insert half-phases
C     ==========================================
C
      DO I = 1,NFMAT
        IREP = ISYMOP(I)-1
        IF(IPRINT.GE.5) THEN
          CALL HEADER('Unsorted DMAT',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(DMAT(1,1,I),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
C
C       Insert half-phases
C       ==================
C
        IF(NZ.LT.4) THEN
          DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            CALL Q2BPHASE('D',IQ,1,DMAT(1,IZ,I))
          ENDDO
          IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted DMAT with half-phases inserted',-1)
            WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
            CALL PRQMAT(DMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
        CALL BSTOBU(DMAT(1,1,I),NZ,WORK,LWORK)
      ENDDO
C
C     ===========================
      IF (ONECAP.AND.SMLV1C) THEN
C     ===========================
C
C        Adapt density matrix in the case that we are going to neglect
C        multicenter integrals in the VSS part. Save original density
C        matrix to file as it should be unchanged once we leave the
C        routine. Use Fock matrix to store projector.
C        TODO : reduce memory use for DMAT in AO basis
C
         CALL MEMGET ('REAL',KBUF1,4*N2BBASX,WORK,KFREE,LFREE)
         CALL OPNFIL(LUPMAT,'DFPMAT','OLD','TWOFCK')
         CALL READT (LUPMAT,N2BBASX,FMAT)
         CALL WRITT (LUPMAT,N2BBASXQ*NFMAT,DMAT)
         CLOSE(LUPMAT,STATUS='KEEP')
         DO I = 1,NFMAT
            IREP = ISYMOP(I)-1
            CALL PR1CEX3 (NTBAS(0),IREP,DMAT(1,1,I),
     &                    FMAT,WORK(KBUF1),IPRINT)
         ENDDO
         CALL MEMREL('TWOFCK',WORK,KBUF1,KBUF1,KFREE,LFREE)
C     =====
      ENDIF
C     =====
C
      IF (DENSFIT) THEN
         NIND = 2 * NUMBER_OF_NUCLEI(2) + 2
         CALL MEMGET ('INTE',KINDAO,NTBAS(0),WORK,KFREE,LFREE)
         CALL MEMGET ('INTE',KINDAOS,NTBAS(0),WORK,KFREE,LFREE)
         CALL MEMGET ('INTE',KINDNAO,NIND,WORK,KFREE,LFREE)
         CALL MEMGET ('REAL',KSDMAT,N2BBASX,WORK,KFREE,LFREE)
         CALL MEMGET ('REAL',KADMAT,N2BBASX,WORK,KFREE,LFREE)
         CALL INDAOC (WORK(KINDAO),WORK(KINDAOS),WORK(KINDNAO))
C        We allocate the fit-coefficients in work for the time being.
C        When we have implemented the fit-basis read routine and the Coulomb build routine, they should
C        probably be allocated at a higher level and passed on to TWOFCK.
         nfit_coef = 1000 ! something to start with as long as we have not implemented the correct structure
         CALL MEMGET ('REAL',KCFIT,nfit_coef,WORK,KFREE,LFREE)
         DO I = 1,NFMAT
            IREP = ISYMOP(I)-1
            IF (IREP.NE.0) 
     &      CALL QUIT ("Only totally symmetric densities can be fitted")
            CALL DTSOUAO(DMAT(1,1,I),WORK(KSDMAT),WORK,KFREE,LFREE)
c           CALL DRIVE_FIT(WORK(KSDMAT),WORK(KADMAT),NTBAS(0),
c    &                     WORK(KINDAO),WORK(KINDAOS),WORK(KINDNAO),
c    &                     IPRINT,WORK(KCFIT))
         ENDDO
         CALL MEMREL('TWOFCK',WORK,KINDAO,KINDAO,KFREE,LFREE)
      ENDIF
C
C     Initialize Fock matrix
C     ======================
C
      CALL DZERO(FMAT,N2BBASXQ*NFMAT)
C
C     *************************************************
C     *****  D i r e c t  c o n s t r u c t i o n *****
C     *****  o f   F o c k   m a t r i c e s      *****
C     *************************************************
C
C     This module has to come before the conventional SCF-module;
C     otherwise FOCK must be copied into buffer; this stems from
C     the symmetrization procedure
C
      IF(IDFLAG.NE.0) THEN
        IF(SOFOCK) THEN
          CALL DIRSOF(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &                IFCKOP,NPOS,INTFLG,
     &                IPRINT,WORK,LWORK)
        ELSE
          CALL DIRSCF(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &                IFCKOP,NPOS,INTFLG,
     &                IPRINT,WORK,LWORK)
        ENDIF
      ENDIF
C
C     **********************************************
C     *****  C o n v e n t i o n a l    S C F  *****
C     **********************************************
C
      IF(IDFLAG.NE.7) THEN
        IF(HFXFAC.NE.D1 .AND. HFXFAC.NE.D0) 
     &       CALL QUIT('HFXFAC not implemented in CONSCF !')
C       hjaaj July 2004: HFXFAC .eq. D0 is defined through SETFCK
        CALL CONSCF(FMAT,DMAT,NFMAT,IFCKOP,ISYMOP,IHRMOP,INTFLG,
     &              IPRINT,WORK,LWORK)
      ENDIF
C
C     ===========================
      IF (ONECAP.AND.SMLV1C) THEN
C     ===========================
C
C        Correct for neglected multi-center contributions via a 1-center
C        expansion that approximates the missing integrals. Restore 
C        density matrix after using it as scratch space.
C        TODO : Reduce memory usage 
C
         CALL MEMGET ('REAL',KBUF1,4*N2BBASX,WORK,KFREE,LFREE)
         CALL OPNFIL(LUPMAT,'DFPMAT','OLD','TWOFCK')
C
C        CAUTION: Here we read the projector in the scratch array DMAT !
C
         CALL READT(LUPMAT,N2BBASX,DMAT)
         DO I = 1,NFMAT
            IREP = ISYMOP(I)-1
            CALL PR1CEX2 (NTBAS(0),IREP,FMAT,DMAT,WORK(KBUF1),IPRINT)
         ENDDO
C
C        Now we do read the original (full) DMAT !
C
         CALL READT(LUPMAT,N2BBASXQ*NFMAT,DMAT)
         CLOSE(LUPMAT,STATUS='KEEP')
         CALL MEMREL('TWOFCK',WORK,KBUF1,KBUF1,KFREE,LFREE)
C     =====
      ENDIF
C     =====
C
C     Insert half-phases and symmetrize
C     =================================
C
      DO I = 1,NFMAT
C
C       Reindex to sorted basis
C       =======================
C
        CALL BUTOBS(FMAT(1,1,I),NZ,WORK,LWORK)
        CALL BUTOBS(DMAT(1,1,I),NZ,WORK,LWORK)
        IREP = ISYMOP(I)-1
        IF(NZ.LT.4) THEN
          DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            CALL Q2BPHASE('F',IQ,1,FMAT(1,IZ,I))
            CALL Q2BPHASE('D',IQ,1,DMAT(1,IZ,I))
          ENDDO
          IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted FMAT with half-phases inserted',-1)
            WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
            CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
C
C       Symmetrize FMAT matrix
C       ======================
C       IH = 0   Non-symmetric
C       IH = 1   Symmetric
C       IH = 2   Anti-symmetric
C
        DO IZ = 1,NZ
          IQ = IPQTOQ(IZ,IREP)
          IH = IHQMAT(IQ,IHRMOP(I)) 
          IF    (IH.EQ.1) THEN
            CALL FULMAT('S',NTBAS(0),NTBAS(0),FMAT(1,IZ,I))
          ELSEIF(IH.EQ.2) THEN
            CALL FULMAT('A',NTBAS(0),NTBAS(0),FMAT(1,IZ,I))
          ENDIF
        ENDDO
        IF(IPRINT.GE.5) THEN
          CALL HEADER('Sorted FMAT matrix',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDDO
      CALL QEXIT('TWOFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pretwo */
      SUBROUTINE TWOPRP(WORK,LWORK)
C*****************************************************************************
C
C     Prepare for two-electron integrals
C
C     Written by T.Saue Sep 30 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "cbihr2.h"
#include "dcbind.h"
#include "dcbbas.h"
#include "hrunit.h"
      DIMENSION WORK(LWORK)
C
C     Initialize COMMON BLOCKS
C     ========================
C
      CALL PAOVEC(WORK,LWORK,0,IPRTWO)
C
C     Initialize HRUNIT
C     =================
C
      LUINTA = 11
      LUSUPM = 17
      LUDASP = 24
      LUSOL  = 49
C
C     Define bit packing
C     ==================
C
#if defined (SYS_CRAY) || defined (SYS_T3D)
      NIBUF = 1
      IF (NTBAS(0) .LE. 255) THEN
         NBITS = 8
      ELSE
         NBITS = 16
      END IF
#else
      IF (NTBAS(0) .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      NBIT1   = 2**NBITS
      NBIT2   = 2**(2*NBITS)
      IBIT1   = NBIT1 - 1
      IBIT2   = NBIT2 - 1
C
C     Calculate integrals to be stored on disk
C     ========================================
C
      CALL TWOGEN(WORK,LWORK)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twosrt */
      SUBROUTINE TWOSRT(ITYP,NCHUNK,IREPA,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C	This files provides sorting of XX-integrals based on index relations
C		Class D,Ca and Cb sorted in buffers of size N2BBASX
C		Class Ba and Bb placed in matrices of dimension NDIM2
C		Class A placed in array of length NDIM
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
#include "dcbbas.h"
#include "dcbind.h"
#include "dgroup.h"
      CHARACTER FILNAM*6
      DIMENSION	WORK(LWORK),IREPA(N2BBASX),LUINT(0:8)
      CALL QENTER('TWOSRT')
#include "memint.h"
      LU2HER    = 13
      DO I = 0,8
        LUINT(I) = LU2HER + 1 + I
      ENDDO
      GO TO (11,12,13,14), ITYP
   11 CONTINUE
C
C     *****************************
C     ***** LLLL - integrals ******
C     *****************************
C
C
C       **************************
C	***** Coarse sorting *****
C       **************************
C
C       File assignment:
C         LUINTA - file of sorted LL-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted LL-integrals (HERMIT-format, deleted in DO1SRT)
C         LU2BUF - buffer file of UNSORTED LL-integrals, other classes
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from LL1SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER  ,'DFTWLL','OLD'    ,'LL1SRT')
        CALL OPNFIL(LUINT(0),'DFLLSA','UNKNOWN','LL1SRT')
        DO I = 1,NBSYM
          FILNAM = 'DFLLBF'
          WRITE(FILNAM(6:6),'(I1)') I-1
          CALL OPNFIL(LUINT(I),FILNAM,'UNKNOWN','LL1SRT')          
        ENDDO
C
        MBUF = NBSYM*2
        CALL MEMGET('REAL',KTWOMT,N2BBASX*MBUF      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDMT,N2BBASX*MBUF*NIBUF,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBUF  ,NCHUNK            ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBUF ,NCHUNK*NIBUF      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIND  ,NCHUNK*2          ,WORK,KFREE,LFREE)
C
        CALL XX1SRT(LUINT,LU2HER,WORK(KTWOMT),WORK(KINDMT),N2BBASX,
     &              WORK(KBUF),WORK(KIBUF),WORK(KIND),IREPA,
     &              NCHUNK,IPRINT)
C
        CALL MEMREL('LL1SRT',WORK,KWORK,KWORK,KFREE,LFREE)
        CLOSE(LU2HER,STATUS='DELETE')
C
C       *************************
C       ***** Fine sorting ******
C       *************************
C
C
        NNDIM = NTBAS(1)*(NTBAS(1)-1)/2
        CALL MEMGET('REAL',KTWOMT,N2BBASX        ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDMT,N2BBASX*NIBUF  ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIND  ,N2BBASX*4      ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KCMAT ,2*N2BBASX      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KICMAT,2*NIBUF*N2BBASX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBMAT ,4*NNDIM        ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBMAT,4*NNDIM        ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KAVEC ,NTBAS(1)       ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIAVEC,NTBAS(1)*NIBUF ,WORK,KFREE,LFREE)
C
        CALL OPNFIL(LU2HER,'DFLLSB','UNKNOWN','LL2SRT')
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from LL2SRT',-1)
        ENDIF
        CALL XX2SRT(LUINT(0),LU2HER,LUINT(1),N2BBASX,
     &           NTBAS(1),NNDIM,WORK(KTWOMT),WORK(KINDMT),WORK(KIND),
     &           WORK(KCMAT),WORK(KICMAT),WORK(KBMAT),WORK(KIBMAT),
     &           WORK(KAVEC),WORK(KIAVEC),IPRINT)
        CLOSE(LUINT(0),STATUS='KEEP')
        CLOSE(LU2HER  ,STATUS='KEEP')
        DO I = 1,NBSYM
          CLOSE(LUINT(I),STATUS = 'DELETE')
        ENDDO
C
        CALL MEMREL('LL2SRT',WORK,KWORK,KWORK,KFREE,LFREE)
      GOTO 10
   12 CONTINUE
C
C     *****************************
C     ***** LLSS - integrals ******
C     *****************************
C
C
C       **************************
C	***** Coarse sorting *****
C       **************************
C
C       File assignment:
C         LUINTA - file of sorted LL-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted LL-integrals (HERMIT-format, deleted in DO1SRT)
C         LU2BUF - buffer file of UNSORTED LL-integrals, other classes
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SL1SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER  ,'DFTWSL','OLD'    ,'SL1SRT')
        CALL OPNFIL(LUINT(0),'DFSLSA','UNKNOWN','SL1SRT')
        DO I = 1,NBSYM
          FILNAM = 'DFSLBF'
          WRITE(FILNAM(6:6),'(I1)') I-1
          CALL OPNFIL(LUINT(I),FILNAM,'UNKNOWN','SL1SRT')          
        ENDDO
C
        MBUF = NBSYM*2
        CALL MEMGET('REAL',KTWOMT,N2BBASX*MBUF      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDMT,N2BBASX*MBUF*NIBUF,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBUF  ,NCHUNK            ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBUF ,NCHUNK*NIBUF      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIND  ,NCHUNK*2          ,WORK,KFREE,LFREE)
C
        CALL XY1SRT(LUINT,LU2HER,WORK(KTWOMT),WORK(KINDMT),N2BBASX,
     &              WORK(KBUF),WORK(KIBUF),WORK(KIND),IREPA,
     &              NCHUNK,IPRINT)
C
        CALL MEMREL('SL1SRT',WORK,KWORK,KWORK,KFREE,LFREE)
        CLOSE(LU2HER,STATUS='DELETE')
C
C       *************************
C       ***** Fine sorting ******
C       *************************
C
C       File assignment:
C         LUINTA - file of sorted LL-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted LL-integrals (HERMIT-format)
C         LU2BUF - buffer file of SORTED integrals, other classes
C
C	Memory assignment
C
        NSLDIM = NTBAS(1)*NTBAS(2)
        CALL MEMGET('REAL',KBUF  ,N2BBASX        ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBUF ,N2BBASX*NIBUF  ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIND  ,N2BBASX*4      ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KAMAT ,NSLDIM         ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIAMAT,NSLDIM         ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBCMAT,2*N2BBASX      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBCMT,2*N2BBASX*NIBUF,WORK,KFREE,LFREE)
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SL2SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER,'DFSLSB','UNKNOWN','SL2SRT')
        CALL XY2SRT(LUINT(0),LU2HER,LUINT(1),N2BBASX,
     &              NSLDIM,WORK(KBUF),WORK(KIBUF),WORK(KIND),
     &              WORK(KAMAT),WORK(KIAMAT),
     &              WORK(KBCMAT),WORK(KIBCMT),IPRINT)
        CLOSE(LUINT(0),STATUS='KEEP')
        CLOSE(LU2HER  ,STATUS='KEEP')
        DO I = 1,NBSYM
          CLOSE(LUINT(I),STATUS = 'DELETE')
        ENDDO
        CALL MEMREL('SL2SRT',WORK,KWORK,KWORK,KFREE,LFREE)
      GOTO 10
   13 CONTINUE
C
C     *****************************
C     ***** SSSS - integrals ******
C     *****************************
C
C       **************************
C	***** Coarse sorting *****
C       **************************
C
C       File assignment:
C         LUINTA - file of sorted SS-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted SS-integrals (HERMIT-format, deleted in DO1SRT)
C         LU2BUF - buffer file of UNSORTED SS-integrals, other classes
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SS1SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER  ,'DFTWSS','OLD'    ,'SS1SRT')
        CALL OPNFIL(LUINT(0),'DFSSSA','UNKNOWN','SS1SRT')
        DO I = 1,NBSYM
          FILNAM = 'DFSSBF'
          WRITE(FILNAM(6:6),'(I1)') I-1
          CALL OPNFIL(LUINT(I),FILNAM,'UNKNOWN','SS1SRT')          
        ENDDO
C
        MBUF = NBSYM*2
        CALL MEMGET('REAL',KTWOMT,N2BBASX*MBUF      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDMT,N2BBASX*MBUF*NIBUF,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBUF  ,NCHUNK            ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBUF ,NCHUNK*NIBUF      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIND  ,NCHUNK*2          ,WORK,KFREE,LFREE)
        CALL XX1SRT(LUINT,LU2HER,WORK(KTWOMT),WORK(KINDMT),N2BBASX,
     &              WORK(KBUF),WORK(KIBUF),WORK(KIND),IREPA,
     &              NCHUNK,IPRINT)
        CALL MEMREL('SS1SRT',WORK,KWORK,KBUF,KFREE,LFREE)
        CLOSE(LU2HER,STATUS='DELETE')
C
C       *************************
C       ***** Fine sorting ******
C       *************************
C
C
        NNDIM = NTBAS(2)*(NTBAS(2)-1)/2
        CALL MEMGET('REAL',KTWOMT,N2BBASX        ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KINDMT,N2BBASX*NIBUF  ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIND  ,N2BBASX*4      ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KCMAT ,2*N2BBASX      ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KICMAT,2*NIBUF*N2BBASX,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KBMAT ,4*NNDIM        ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBMAT,4*NNDIM        ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KAVEC ,NTBAS(2)       ,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIAVEC,NTBAS(2)*NIBUF ,WORK,KFREE,LFREE)
C
        CALL OPNFIL(LU2HER,'DFSSSB','UNKNOWN','SS2SRT')
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SS2SRT',-1)
        ENDIF
        CALL XX2SRT(LUINT(0),LU2HER,LUINT(1),N2BBASX,
     &           NTBAS(2),NNDIM,WORK(KTWOMT),WORK(KINDMT),WORK(KIND),
     &           WORK(KCMAT),WORK(KICMAT),WORK(KBMAT),WORK(KIBMAT),
     &           WORK(KAVEC),WORK(KIAVEC),IPRINT)
        CLOSE(LUINT(0),STATUS='KEEP')
        CLOSE(LU2HER  ,STATUS='KEEP')
        DO I = 1,NBSYM
          CLOSE(LUINT(I),STATUS = 'DELETE')
        ENDDO
C
        CALL MEMREL('SS2SRT',WORK,KWORK,KWORK,KFREE,LFREE)
      GOTO 10
   14 CONTINUE
C
C     Gaunt integrals
C
      GOTO 10
   10 CONTINUE
C
      CALL QEXIT('TWOSRT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xx1srt */
      SUBROUTINE XX1SRT(LTFIL,LUTWXX,TWOMAT,INDMAT,MAXBUF,
     &                  BUF,IBUF,IND,IREPA,NCHUNK,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Coarse sorting of 2-electron XX - integrals
C
C     METHOD:
C
C     FILE ASSIGNMENT:
C      LUTWXX   - 'AOTWOXX'    - file of unsorted XX-integrals
C      LTFIL(0) - 'DF2XXA'     - file of sorted XX-integrals (D,Cb)
C      LTFIL(1) - 'DF2XXC'     - buffer file of XX-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size N2BBASX) to be read
C       1 - XX, class D
C       2 - XX-buffer
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C
C************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION LTFIL(0:8),ICOUNT(0:7,0:1),NINTS(0:7,0:1)
      DIMENSION IBUF(NCHUNK,NIBUF),BUF(NCHUNK),IND(2,NCHUNK)
      DIMENSION TWOMAT(MAXBUF,0:MAXREP,0:1),
     &          INDMAT(MAXBUF,NIBUF,0:MAXREP,0:1),
     &          IREPA(*)
#include "ibtfun.h"
C
      CALL IZERO(NINTS  ,16)
      CALL IZERO(ICOUNT,16)
      REWIND LUTWXX
      CALL MOLLAB('BASTWOEL',LUTWXX,LUPRI)
C*** DO WHILE - loop: read a chunk of integrals, perform  sorting
      MINT = 0
   10 CONTINUE
      READ(LUTWXX,IOSTAT = IOS) BUF,IBUF,NUT
      IF(NUT.LE.0) GOTO 30
      MINT = MINT + NUT
      IF (NIBUF .EQ. 1) THEN
C
C       Unpack indices and extract information about
C         1. Combined symmetry of IJ
C         2. Whether any indices are equal
C
        DO INT = 1,NUT
          I     = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
          J     = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
          K     = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
          L     = IBTAND(       IBUF(INT,1),         IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - K)*(I - J)*(J - K)*(J - L)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      ELSE
C
C       Unpack indices and extract information about
C         1. Whether any indices are equal
C         2. Combined symmetry of IJ
C
        DO INT = 1,NUT
          I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
          J = IBTAND(       IBUF(INT,1),       IBIT1)
          K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
          L = IBTAND(       IBUF(INT,2),       IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - K)*(I - J)*(J - K)*(J - L)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
          INDMAT(INDEX,2,IREP,IT) = IBUF(INT,2)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      END IF
      GOTO 10
C
C*** END DO WHILE - loop: reading of integrals from LUTWXX
C
   30 CONTINUE
C*** Empty matrices:
      DO IT = 0,1
        DO IREP = 0,MAXREP
          IF(NINTS(IREP,IT).GT.0) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),NINTS(IREP,IT),IREP,MAXBUF,NIBUF,
     &               INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            NINTS(IREP,IT) = NINTS(IREP,IT) + ICOUNT(IREP,IT)*MAXBUF
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
          ELSE
            NINTS(IREP,IT) = ICOUNT(1,1)*MAXBUF
          ENDIF
        ENDDO
      ENDDO       
      DO IUNIT = 0,MAXREP+1
        WRITE(LTFIL(IUNIT)) -1,0
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - coarse sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class D integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Other classes    :'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,0),IREP=0,MAXREP)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xy1srt */
      SUBROUTINE XY1SRT(LTFIL,LUTWXX,TWOMAT,INDMAT,MAXBUF,
     &                  BUF,IBUF,IND,IREPA,NCHUNK,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Coarse sorting of 2-electron XY - integrals
C
C     METHOD:
C
C     FILE ASSIGNMENT:
C      LUTWXX   - 'AOTWOXX'    - file of unsorted XX-integrals
C      LTFIL(0) - 'DF2XXA'     - file of sorted XX-integrals (D,Cb)
C      LTFIL(1) - 'DF2XXC'     - buffer file of XX-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size N2BBASX) to be read
C       1 - XX, class D
C       2 - XX-buffer
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C
C************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION LTFIL(0:8),ICOUNT(0:7,0:1),NINTS(0:7,0:1)
      DIMENSION IBUF(NCHUNK,NIBUF),BUF(NCHUNK),IND(2,NCHUNK)
      DIMENSION TWOMAT(MAXBUF,0:MAXREP,0:1),
     &          INDMAT(MAXBUF,NIBUF,0:MAXREP,0:1),
     &          IREPA(*)
#include "ibtfun.h"
C
      CALL IZERO(NINTS  ,16)
      CALL IZERO(ICOUNT,16)
      REWIND LUTWXX
      CALL MOLLAB('BASTWOEL',LUTWXX,LUPRI)
C*** DO WHILE - loop: read a chunk of integrals, perform  sorting
      MINT = 0
   10 CONTINUE
      READ(LUTWXX,IOSTAT = IOS) BUF,IBUF,NUT
      IF(NUT.LE.0) GOTO 30
      MINT = MINT + NUT
      IF (NIBUF .EQ. 1) THEN
C
C       Unpack indices and extract information about
C         1. Combined symmetry of IJ
C         2. Whether any indices are equal
C
        DO INT = 1,NUT
          I     = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
          J     = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
          K     = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
          L     = IBTAND(       IBUF(INT,1),         IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - J)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      ELSE
C
C       Unpack indices and extract information about
C         1. Whether any indices are equal
C         2. Combined symmetry of IJ
C
        DO INT = 1,NUT
          I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
          J = IBTAND(       IBUF(INT,1),       IBIT1)
          K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
          L = IBTAND(       IBUF(INT,2),       IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - J)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
          INDMAT(INDEX,2,IREP,IT) = IBUF(INT,2)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      END IF
      GOTO 10
C
C*** END DO WHILE - loop: reading of integrals from LUTWXX
C
   30 CONTINUE
C*** Empty matrices:
      DO IT = 0,1
        DO IREP = 0,MAXREP
          IF(NINTS(IREP,IT).GT.0) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),NINTS(IREP,IT),IREP,MAXBUF,NIBUF,
     &               INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            NINTS(IREP,IT) = NINTS(IREP,IT) + ICOUNT(IREP,IT)*MAXBUF
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
          ELSE
            NINTS(IREP,IT) = ICOUNT(1,1)*MAXBUF
          ENDIF
        ENDDO
      ENDDO       
      DO IUNIT = 0,MAXREP+1
        WRITE(LTFIL(IUNIT)) -1,0
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - coarse sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class D integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Other classes    :'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,0),IREP=0,MAXREP)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xx2srt */
      SUBROUTINE XX2SRT(LUINTA,LUINTB,LUBUF,MAXBUF,
     &                  NDIM,NNDIM,BUF,IBUF,IND,CMAT,ICMAT,
     &                  BMAT,IBMAT,AVEC,IAVEC,IPRINT)
C*****************************************************************************
C
C     PURPOSE:
C        Fine sorting of 2-electron XX-integrals:
C
C     METHOD:
C        Twoelectron XX-integrals IJKL can be classified as
C        follows according to indices:
C
C        *Group A: All indices equal
C            Class  1: I = J = K = L
C        *Group B: Three indices equal
C            Class a1: I > J = K = L
C            Class a2: I = J = K > L
C            Class b1: I = J > K = L
C            Class b2: I = K > J = L
C        *Group C: two indices equal
C            Class a1: I = K > J > L
C            Class a2: I > J = K > L
C            Class a3: I > K > J = L
C            Class b1: I > K = L > J
C            Class b2: I > J > K = L
C            Class b3: I = J > K > L
C        *Group D : all indices different(sorted in XX1SRT):
C            Class  1: I > K > J > L
C            Class  2: I > K > L > J
C            Class  3: I > J > K > L
C
C        In this routine the integrals are sorted, that is
C        distributed on different files, according to how they
C        contribute to the twoelectron Fock matrix
C
C     FILE ASSIGNMENT:
C      LUINTA - 'DF2XXA'     - file of sorted XX-integrals (D,Cb)
C      LUINTB - 'DF2XXB'     - file of sorted XX-integrals (Ca,Bb,Ba,A)
C      LU2BUF - 'DF2XXC'     - buffer file of XX-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size N2BBASX) to be read
C       1 - Input: XX-buffer.  Output: XX, class Ca
C       2 - XX, class Cb
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C     Written by T.Saue Feb 27 1995
C     Last revision: Feb 27 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION BUF(MAXBUF),IBUF(MAXBUF,NIBUF),IND(4,MAXBUF),
     &          CMAT(MAXBUF,2),ICMAT(MAXBUF,NIBUF,2),
     &          BMAT(NNDIM,4)  ,IBMAT(NNDIM,4),
     &          AVEC(NDIM),IAVEC(NDIM),LUBUF(0:7)
      DIMENSION NC(0:7,2),ICOUNT(0:7,2),NB(0:7,4),NBT(4)
#include "ibtfun.h"
C
C     Initialize
C
      REWIND LUINTB
      CALL IZERO(ICOUNT,16)
      CALL IZERO(NC,16)
      CALL IZERO(NB,32)
      CALL IZERO(NBT,4)
      NA = 0
C
C     Loop over irreps
C      
      DO IREP = 0,MAXREP
        LUINTC = LUBUF(IREP)
        REWIND LUINTC
C
C       Loop over batches in buffer
C
   10   CONTINUE
        CALL REAINT(LUINTC,NINTS,IRP,MAXBUF,NIBUF,IBUF,BUF)
        IF(NINTS.EQ.-1) GOTO 30
C
C       Unpack indices
C
        IF (NIBUF .EQ. 1) THEN
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
            IND(2,INT) = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,1),         IBIT1)
          ENDDO
        ELSE
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
            IND(2,INT) = IBTAND(       IBUF(INT,1),       IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,2),       IBIT1)
          ENDDO
        END IF
C
C       Begin sorting
C
        DO 50 INT = 1,NINTS
          I = IND(1,INT)
          J = IND(2,INT)
          K = IND(3,INT)
          L = IND(4,INT)
          IJ = I - J
          JL = J - L
          IF((IJ*JL).NE.0) THEN
C
C***  Group C: Class a1,a2,b1,b2
C
            KL = K - L
            IF(KL.EQ.0) THEN
C
C*** Group C:  Class b1,b2 need not be distinguished
C
              IND1 = NBIT1*I + J
              IND2 = K
              NC(IREP,2) = NC(IREP,2) + 1
              CMAT(NC(IREP,2),2) = BUF(INT)
              IF(NIBUF.EQ.1) THEN
                ICMAT(NC(IREP,2),1,2) = IND1*NBIT2+IND2
              ELSE
                ICMAT(NC(IREP,2),1,2) = IND1
                ICMAT(NC(IREP,2),2,2) = IND2
              ENDIF
              IF(NC(IREP,2).EQ.MAXBUF) THEN
                CALL WRINT(LUINTB,NC(IREP,2),IREP,MAXBUF,NIBUF,
     &                     ICMAT(1,1,2),CMAT(1,2))
                ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
                NC(IREP,2) = 0
              ENDIF
            ELSE
C
C***  Group C: Class a1,a2
C
              IF(I.EQ.K) THEN
C**** Group C: Class a1:
                IND1 = NBIT1*J + L
                IND2 = I
              ELSE
C**** Group C: Class a2:
                IND1 = NBIT1*I + L
                IND2 = J
              ENDIF
              NC(IREP,1) = NC(IREP,1) + 1
              CMAT(NC(IREP,1),1) = BUF(INT)
              IF(NIBUF.EQ.1) THEN
                ICMAT(NC(IREP,1),1,1) = IND1*NBIT2+IND2
              ELSE
                ICMAT(NC(IREP,1),1,1) = IND1
                ICMAT(NC(IREP,1),2,1) = IND2
              ENDIF
              IF(NC(IREP,1).EQ.MAXBUF) THEN
                CALL WRINT(LUINTA,NC(IREP,1),IREP,MAXBUF,NIBUF,
     &                     ICMAT(1,1,1),CMAT(1,1))
                ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
                NC(IREP,1) = 0
              ENDIF
            ENDIF
          ELSE
            IK = I - K
            KL = K - L
            IF((IK*KL).GT.0) THEN
C
C*** Group C : Class a3,b3
C
              IF(IJ.EQ.0) THEN
C
C*** Group C: Class b3
C
                NC(IREP,2)         = NC(IREP,2) + 1
                CMAT(NC(IREP,2),2) = BUF(INT)
                IND1               = NBIT1*K + L
                IND2               = I
                IF(NIBUF.EQ.1) THEN
                  ICMAT(NC(IREP,2),1,2) = IND1*NBIT2+IND2
                ELSE
                  ICMAT(NC(IREP,2),1,2) = IND1
                  ICMAT(NC(IREP,2),2,2) = IND2
                ENDIF
                IF(NC(IREP,2).EQ.MAXBUF) THEN
                  CALL WRINT(LUINTB,NC(IREP,2),IREP,MAXBUF,NIBUF,
     &                       ICMAT(1,1,2),CMAT(1,2))
                  ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
                  NC(IREP,2) = 0
                ENDIF
              ELSE
C
C*** Group C: Class a3
C
                NC(IREP,1)         = NC(IREP,1) + 1
                CMAT(NC(IREP,1),1) = BUF(INT)
                IND1          = NBIT1*I + K
                IND2          = J
                IF(NIBUF.EQ.1) THEN
                  ICMAT(NC(IREP,1),1,1) = IND1*NBIT2+IND2
                ELSE
                  ICMAT(NC(IREP,1),1,1) = IND1
                  ICMAT(NC(IREP,1),2,1) = IND2
                ENDIF
                IF(NC(IREP,1).EQ.MAXBUF) THEN
                    CALL WRINT(LUINTA,NC(IREP,1),IREP,MAXBUF,NIBUF,
     &                       ICMAT(1,1,1),CMAT(1,1))
                  ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
                  NC(IREP,1) = 0
                ENDIF
              ENDIF
            ELSEIF(K.GT.J) THEN
C
C*** Group B - b2
C
              NB(IREP,1)      = NB(IREP,1) + 1
              NBT(1)          = NBT(1) + 1
              BMAT(NBT(1),1)  = BUF(INT)
              IBMAT(NBT(1),1) = NBIT1*J+K
            ELSEIF(K.LT.J) THEN
C
C*** Group B - b1
C
              NB(IREP,4)      = NB(IREP,4) + 1
              NBT(4)          = NBT(4) + 1
              BMAT(NBT(4), 4) = BUF(INT)
              IBMAT(NBT(4),4) = NBIT1*J+K
            ELSEIF(I.GT.J) THEN
C
C*** Group B - a1
C
              NB(IREP,2)      = NB(IREP,2)+1
              NBT(2)          = NBT(2)+1
              BMAT(NBT(2), 2) = BUF(INT)
              IBMAT(NBT(2),2) = NBIT1*I+J
            ELSEIF(I.GT.L) THEN
C
C*** Group B - a2
C
              NB(IREP,3)      = NB(IREP,3) + 1
              NBT(3)          = NBT(3) + 1
              BMAT(NBT(3), 3) = BUF(INT)
              IBMAT(NBT(3),3) = NBIT1*L + I
            ELSE
C
C*** Group A:
C
              NA        = NA + 1
              AVEC(NA)  = BUF(INT)
              IAVEC(NA) = I
            ENDIF
          ENDIF
 50       CONTINUE
        GOTO 10
   30   CONTINUE
C*** Empty matrices:
        IF(NC(IREP,1).GT.0) THEN
          CALL WRINT(LUINTA,NC(IREP,1),IREP,MAXBUF,NIBUF,
     &               ICMAT(1,1,1),CMAT(1,1))
          NC(IREP,1)  = NC(IREP,1) + ICOUNT(IREP,1)*MAXBUF
          ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
        ELSE
          NC(IREP,1) = ICOUNT(IREP,1)*MAXBUF
        ENDIF
        IF(NC(IREP,2).GT.0) THEN
          CALL WRINT(LUINTB,NC(IREP,2),IREP,MAXBUF,NIBUF,
     &               ICMAT(1,1,2),CMAT(1,2))
          NC(IREP,2) = NC(IREP,2) + ICOUNT(IREP,2)*MAXBUF
          ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
        ELSE
          NC(IREP,2) = ICOUNT(IREP,2)*MAXBUF
        ENDIF
      ENDDO
C
      WRITE(LUINTA) -1,0
      WRITE(LUINTB) -1,0
      IOFF = 1
      DO IREP = 0,MAXREP
        CALL WRINT(LUINTB,NB(IREP,1),IREP,NNDIM,1,
     &                      IBMAT(IOFF,1),BMAT(IOFF,1))
        IOFF = IOFF + NB(IREP,1)
      ENDDO
C     Class Ba1, Ba2 and Bb1 has totally symmetric densities !
      DO I = 2,4
        CALL WRINT(LUINTB,NB(0,I),0,NNDIM,1,
     &                      IBMAT(1,I),BMAT(1,I))
      ENDDO
      CALL WRINT(LUINTB,NA,0,NDIM,1,IAVEC,AVEC)
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - fine sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class Ca integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NC(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Cb integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NC(IREP,2),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Ba1 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,2),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Ba2 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,3),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Bb1 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,4),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Bb2 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class A   integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       'Irrep ',REP(0),NA
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrint */
      SUBROUTINE WRINT(IUNIT,NINTS,IREP,NDIM,NIBUF,IBUF,BUF)
C*****************************************************************************
C
C     Write integral and indices to unformatted file
C
C     Written by T.Saue, Feb 27 1995
C     Last revision: Feb 27 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BUF(NDIM),IBUF(NDIM,NIBUF)
C
      WRITE(IUNIT) NINTS,IREP
      IF(NINTS.LE.0) RETURN
      CALL WRITI(IUNIT,NINTS,IBUF(1,1))
      IF(NIBUF.EQ.2) CALL WRITI(IUNIT,NINTS,IBUF(1,2))
      CALL WRITT(IUNIT,NINTS,BUF)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reaint */
      SUBROUTINE REAINT(IUNIT,NINTS,IREP,NDIM,NIBUF,IBUF,BUF)
C*****************************************************************************
C
C     Read integral and indices from unformatted file
C
C     Written by T.Saue, Feb 27 1995
C     Last revision: Feb 27 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BUF(NDIM),IBUF(NDIM,NIBUF)
C
      READ(IUNIT) NINTS,IREP
      IF(NINTS.LE.0) RETURN
      CALL READI(IUNIT,NINTS,IBUF(1,1))
      IF(NIBUF.EQ.2) CALL READI(IUNIT,NINTS,IBUF(1,2))
      CALL READT(IUNIT,NINTS,BUF)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xy2srt */
      SUBROUTINE XY2SRT(LUINTA,LUINTB,LUBUF,MAXBUF,
     &                  NSLDIM,BUF,IBUF,IND,AMAT,IAMAT,
     &                  BCMAT,IBCMAT,IPRINT)
C*****************************************************************************
C
C     PURPOSE:
C        Fine sorting of 2-electron XY-integrals:
C
C     METHOD:
C        Two-electron XY-integrals IJKL can be classified as
C        follows according to indices:
C
C        *Group A: I = J,  K = L
C        *Group B: I > J,  K = L
C        *Group C: I = J,  K > L
C        *Group D: I > J,  K > L
C
C        In this routine the integrals are sorted, that is
C        distributed on different files, according to how they
C        contribute to the twoelectron Fock matrix
C
C     FILE ASSIGNMENT:
C      LUINTA - 'DF2XYA'     - file of sorted XY-integrals (D,B)
C      LUINTB - 'DF2XYB'     - file of sorted XY-integrals (C,A)
C      LU2BUF - 'DF2XYC'     - buffer file of XY-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size MAXBUF) to be read
C       1 - Input: XY-buffer.  Output: XY, class Ca
C       2 - XY, class Cb
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C     Written by T.Saue Feb 27 1995
C     Last revision : Feb 27 1995
C
C************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION BUF(MAXBUF),IBUF(MAXBUF,NIBUF),IND(4,MAXBUF),
     &          AMAT(NSLDIM),IAMAT(NSLDIM),
     &          BCMAT(MAXBUF,2),IBCMAT(MAXBUF,NIBUF,2),
     &          ICOUNT(0:7,2),LUBUF(0:7),NB(0:7),NC(0:7)
#include "ibtfun.h"
      REWIND LUINTB
      CALL IZERO(ICOUNT,16)
      CALL IZERO(NC,8)
      CALL IZERO(NB,8)
      NA = 0
C
C     Loop over irreps
C      
      DO IREP = 0,MAXREP
        LUINTC = LUBUF(IREP)
        REWIND LUINTC
C
C       Loop over batches in buffer
C
   10   CONTINUE
        CALL REAINT(LUINTC,NINTS,IRP,MAXBUF,NIBUF,IBUF,BUF)
        IF(NINTS.EQ.-1) GOTO 30
C
C       Unpack indices
C
        IF (NIBUF .EQ. 1) THEN
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
            IND(2,INT) = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,1),         IBIT1)
          ENDDO
        ELSE
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
            IND(2,INT) = IBTAND(       IBUF(INT,1),       IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,2),       IBIT1)
          ENDDO
        END IF
C
C       Begin sorting
C
        DO INT = 1,NINTS
          I = IND(1,INT)
          J = IND(2,INT)
          K = IND(3,INT)
          L = IND(4,INT)
          IF(I.GT.J) THEN
C
C*** Group B : IJKK
C
            NB(IREP)           = NB(IREP) + 1
            BCMAT(NB(IREP),1)  = BUF(INT)
            IND1         = NBIT1*I + J
            IND2         = K
            IF(NIBUF.EQ.1) THEN
              IBCMAT(NB(IREP),1,1) = IND1*NBIT2+IND2
            ELSE
              IBCMAT(NB(IREP),1,1) = IND1
              IBCMAT(NB(IREP),2,1) = IND2
            ENDIF
            IF(NB(IREP).EQ.MAXBUF) THEN
              CALL WRINT(LUINTA,NB(IREP),IREP,MAXBUF,NIBUF,
     &                   IBCMAT(1,1,1),BCMAT(1,1))
              ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
              NB(IREP) = 0
            ENDIF
          ELSEIF(K.GT.L) THEN
C
C*** Group C : IIKL
C
            IND1         = NBIT1*K + L
            IND2         = I
            NC(IREP)           = NC(IREP) + 1
            BCMAT(NC(IREP),2)  = BUF(INT)
            IF(NIBUF.EQ.1) THEN
              IBCMAT(NC(IREP),1,2) = IND1*NBIT2+IND2
            ELSE
              IBCMAT(NC(IREP),1,2) = IND1
              IBCMAT(NC(IREP),2,2) = IND2
            ENDIF
            IF(NC(IREP).EQ.MAXBUF) THEN
              CALL WRINT(LUINTB,NC(IREP),IREP,MAXBUF,NIBUF,
     &                   IBCMAT(1,1,2),BCMAT(1,2))
              ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
              NC(IREP) = 0
            ENDIF
          ELSE
C
C*** Group A : IIKK
C
            NA        = NA + 1
            AMAT(NA)  = BUF(INT)
            IAMAT(NA) = NBIT1*I + K
          ENDIF
        ENDDO
        GOTO 10
   30   CONTINUE
C*** Empty matrices:
        IF(NB(IREP).GT.0) THEN
          CALL WRINT(LUINTA,NB(IREP),IREP,MAXBUF,NIBUF,
     &               IBCMAT(1,1,1),BCMAT(1,1))
          NB(IREP) = NB(IREP) + ICOUNT(IREP,1)*MAXBUF
          ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
        ELSE
          NB(IREP) = ICOUNT(IREP,1)*MAXBUF
        ENDIF
        IF(NC(IREP).GT.0) THEN
          CALL WRINT(LUINTB,NC(IREP),IREP,MAXBUF,NIBUF,
     &              IBCMAT(1,1,2),BCMAT(1,2))
          NC(IREP) = NC(IREP) + ICOUNT(IREP,2)*MAXBUF
          ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
        ELSE
          NC(IREP) = ICOUNT(IREP,2)*MAXBUF
        ENDIF
      ENDDO
C
      WRITE(LUINTA) -1,0
      WRITE(LUINTB) -1,0
      CALL WRINT(LUINTB,NA,0,NSLDIM,1,IAMAT,AMAT)
C
C     Output section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - fine sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class C integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NC(IREP),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class B integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)')  ' - Class A integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       'Irrep ',REP(0),NA
      ENDIF
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twogen */
      SUBROUTINE TWOGEN(WORK,LWORK)
C*****************************************************************************
C
C     Do all necessary preparations for two-electron integrals
C
C     Written by T.Saue May 18 1995
C     Last revision : Sep 30 1996 - tsaue
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "cbihr2.h"
#include "dcbbas.h"
#include "dcbgen.h"
      DIMENSION WORK(LWORK)
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('INTE',KIREPA,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL TWOGE1(WORK(KIREPA),WORK(KFREE),LFREE)
C
C     Memory deallocation
C
      CALL MEMREL('TWOGEN',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twoge1 */
      SUBROUTINE TWOGE1(IREPA,WORK,LWORK)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr2.h"
#include "hrunit.h"
      LOGICAL LBIT,TOBE,REGO
      DIMENSION IREPA(*),WORK(LWORK)
C
      NCHUNK = 600
      IRNTYP = 0
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
C
C     Make array of boson irrep for each basis function
C
      CALL RPGEN(IREPA)
C
C     LL - integrals
C
 10   CONTINUE
      IF(.NOT.LBIT(IDFLAG,1).AND.LBIT(INTGEN,1)) THEN
C
C       Check for restart
C
cjkp+hjaaj: do not allow restart; we may have a new geometry
c       INQUIRE(FILE='DFLLSA',EXIST=TOBE)
c       REGO = TOBE
c       INQUIRE(FILE='DFLLSB',EXIST=TOBE)
c       REGO = TOBE.AND.REGO
c       IF(REGO) THEN
c         WRITE(LUPRI,'(A)')
c    &     '*** INFO *** : TWOGEN - LL-integrals from restart'
c         GOTO 20
c       ENDIF
cjkp+hjaaj end
C
C       Generate integrals
C
        I2TYP = 1
        CALL OPNFIL(LUINTA,'DFTWLL','UNKNOWN','TWOGE1')          
        IF(IPRTWO.GE.2)
     +     WRITE(LUPRI,'(1P,A,D9.2)') 'LL-int - threshhold:',THRS
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL TWOINT(WORK,LWORK,DUM,DUM,1,IDUM,IDUM,DUM,
     &              IDUM,IDUM,IRNTYP,MAXDIF,JATOM,
     &              .TRUE.,.TRUE.,.FALSE.,TKTIME,IPRTWO,
     &              IPRNTA,IPRNTB,IPRNTC,IPRNTD,RTNTWO,IDUM,I2TYP,
     &              IDUM,IDUM,DUM,DUM,DUM,IDUM,.TRUE.,
     &              IDUM,DUM)
        CALL TIMER('LL-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
C
C     SL - integrals
C
 20   CONTINUE
      IF(.NOT.LBIT(IDFLAG,2).AND.LBIT(INTGEN,2)) THEN
C
C       Check for restart
C
        INQUIRE(FILE='DFSLSA',EXIST=TOBE)
        REGO = TOBE
        INQUIRE(FILE='DFSLSB',EXIST=TOBE)
        REGO = TOBE.AND.REGO
C
C       IF one-center model 2,3,4 is used for geometry optimization,
C       do not generate integrals since these are one-center
C       and therefore unchanged during minimization.   /jkp
C
        IF(REGO) THEN
          IF ((ONECAP.AND.INTV1C.NE.1).AND.OPTIMI) THEN
            WRITE(LUPRI,'(/2A)')
     &      '** One-center LS-integrals unchanged in this geometry -->',
     &      '   Read from disk.'
            GOTO 30
cjkp+hjaaj: do not allow restart; we may have a new geometry
c         ELSE 
c            WRITE(LUPRI,'(A)')
c    &      '*** INFO *** : TWOGEN - SL-integrals from restart'
          ENDIF
        ENDIF
C
C       Generate integrals
C
        I2TYP  = 2
        THRBUF = THRS
        THRS   = THRS*THRFAC(1)
        IF(IPRTWO.GE.2)
     +     WRITE(LUPRI,'(1P,A,D9.2)') 'SL-int - threshhold:',THRS
        CALL OPNFIL(LUINTA,'DFTWSL','UNKNOWN','TWOGE1')          
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL TWOINT(WORK,LWORK,DUM,DUM,1,IDUM,IDUM,DUM,
     &              IDUM,IDUM,IRNTYP,MAXDIF,JATOM,
     &              .TRUE.,.TRUE.,.FALSE.,TKTIME,IPRTWO,
     &              IPRNTA,IPRNTB,IPRNTC,IPRNTD,RTNTWO,IDUM,I2TYP,
     &              IDUM,IDUM,DUM,DUM,DUM,IDUM,.TRUE.,
     &              IDUM,DUM)
        CALL TIMER('SL-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        THRS = THRBUF
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
C
C     SS - integrals
C
 30   CONTINUE
      IF(.NOT.LBIT(IDFLAG,3).AND.LBIT(INTGEN,3)) THEN
C
C       Check for restart
C
C        INQUIRE(FILE='DFSSTA',EXIST=TOBE)
C        REGO = TOBE
C        INQUIRE(FILE='DFSSTB',EXIST=TOBE)
C        REGO = TOBE.AND.REGO
        INQUIRE(FILE='DFSSSA',EXIST=TOBE)
        REGO = TOBE
        INQUIRE(FILE='DFSSSB',EXIST=TOBE)
        REGO = TOBE.AND.REGO
C
C       IF one-center model 2,3,4 is used for geometry optimization,
C       do not generate integrals since these are one-center
C       and therefore unchanged during minimization.   /jkp
C
        IF(REGO) THEN
          IF ((ONECAP.AND.INTV1C.NE.1).AND.OPTIMI) THEN
            WRITE(LUPRI,'(/2A)')
     &      '** One-center SS-integrals unchanged in this geometry -->',
     &      '   Read from disk.'
            GOTO 40
cjkp+hjaaj: do not allow restart; we may have a new geometry
c         ELSE 
c           WRITE(LUPRI,'(A)')
c    &      '*** INFO *** : TWOGEN - SS-integrals from restart'
          ENDIF
        ENDIF
C
C       Generate integrals
C
        I2TYP  = 3
        THRBUF = THRS
        THRS   = THRS*THRFAC(2)
        IF(IPRTWO.GE.2)
     +     WRITE(LUPRI,'(A,E9.2)') 'SS-int - threshhold:',THRS
        CALL OPNFIL(LUINTA,'DFTWSS','UNKNOWN','TWOGE1')          
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL TWOINT(WORK,LWORK,DUM,DUM,1,IDUM,IDUM,DUM,
     &              IDUM,IDUM,IRNTYP,MAXDIF,JATOM,
     &              .TRUE.,.TRUE.,.FALSE.,TKTIME,IPRTWO,
     &              IPRNTA,IPRNTB,IPRNTC,IPRNTD,RTNTWO,IDUM,I2TYP,
     &              IDUM,IDUM,DUM,DUM,DUM,IDUM,.TRUE.,
     &              IDUM,DUM)
        CALL TIMER('SS-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        THRS = THRBUF
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
 40   CONTINUE
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xx2mat */
      SUBROUTINE XX2MAT(LUINTA,LUINTB,NDIM,NNDIM,
     +  FMAT,DMAT,NDMAT,IFCTYP,IREPDM,BUF,IBUF,IND,IPRINT)
C*****************************************************************************
C
C     Insert LL/SS-integrals in Fock matrix
C
C     Written by T.Saue Mar 3 1995
C     Last revision : tsaue - Mar 3 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
      DIMENSION IBUF(*),BUF(*),FMAT(*),DMAT(*),IND(*),IFCTYP(*),
     &          IREPDM(*)
C
C     *********************************************
C     ***** Process integrals on file LUINTA ******
C     *********************************************
C
 10   CONTINUE
      REWIND LUINTA
C
C     Class D integrals:
C     ==================
C        I > K > J > L
C        I > K > L > L
C        I > J > K > L
C
      CALL DOIJKL(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class Ca integrals:
C     ===================
C        I = K > J > L
C        I > J = K > L
C        I > K > J = L
C
      CALL DOIJIL(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     *********************************************
C     ***** Process integrals on file LUINTB ******
C     *********************************************
C
 20   CONTINUE
      REWIND LUINTB
C
C     Class Cb integrals:
C     ===================
C        I > K = L > J
C        I > J > K = L
C        I = J > K > L
C
      CALL DOIJKK(LUINTB,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class Bb2 integrals:
C     ====================
C        I = K > J = L
C
      DO IREP = 1,NBSYM
        CALL DOIJIJ(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
      ENDDO
C
C     Class Ba1 integrals:
C     ====================
C       I > J = K = L
C
      CALL DOIJJJ(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
C
C     Class Ba2 integrals:
C     ====================
C        I = J = K > L
C
      CALL DOIJJJ(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
C
C     Class Bb1 integrals:
C     ====================
C        I = J > K = L
C
      CALL DOIIKK(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
C
C     Class A
C     =======
C
      CALL DOIIII(LUINTB,NDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IPRINT)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xy2mat */
      SUBROUTINE XY2MAT(LUINTA,LUINTB,
     +  FMAT,DMAT,NDMAT,IFCTYP,IREPDM,BUF,IBUF,IND,IPRINT)
C*****************************************************************************
C
C     Insert LL/SS-integrals in Fock matrix
C
C     Written by T.Saue Mar 3 1995
C     Last revision : tsaue - Mar 3 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
      DIMENSION IBUF(*),BUF(*),FMAT(*),DMAT(*),IND(*),
     &          IFCTYP(*),IREPDM(*)
C
C     *********************************************
C     ***** Process integrals on file LUINTA ******
C     *********************************************
C
      REWIND LUINTA
C
C     Class D integrals:
C     ==================
C       I > J, K > L
C
      CALL DOIJKL(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class B integrals:
C     ==================
C       I > J, K = L
C
      CALL DOIJKK(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     *********************************************
C     ***** Process integrals on file LUINTB ******
C     *********************************************
C
      REWIND LUINTB
C
C     Class C integrals:
C     ==================
C        I = J, K > L
C
      CALL DOIJKK(LUINTB,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class A integrals:
C     ==================
C        I = J, K = L
C
      NSLDIM = NTBAS(1)*NTBAS(2)
      CALL DOIIKK(LUINTB,NSLDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &              IREPDM,IBUF,BUF,IND,IPRINT)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijkl */
      SUBROUTINE DOIJKL(IUNIT,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &                  IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals with all indices unequal
C
C     Written by T.Saue Mar 5 1995
C     Last revision: Jan 7 1998 - jth (Linux cheat)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,*),DMAT(NBAST,NBAST,*) 
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),*),DMAT(NTBAS(0),NTBAS(0),*) 
#endif
      DIMENSION BUF(N2BBASX),IBUF(N2BBASX,NIBUF),IND(4,N2BBASX),
     &          IFCTYP(NDMAT),IREPDM(NDMAT)
      INTEGER A,B,C,D
#include "ibtfun.h"
C
   10   CONTINUE
        CALL REAINT(IUNIT,NBUF,IREPAB,N2BBASX,NIBUF,IBUF,BUF)
        IF(NBUF.EQ.-1) RETURN
C
C       Unpack indices
C       ==============
C
        IF (NIBUF .EQ. 1) THEN
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),3*NBITS),IBIT1)
            IND(2,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),2*NBITS),IBIT1)
            IND(3,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),  NBITS),IBIT1)
            IND(4,INTS) = IBTAND(       IBUF(INTS,1),         IBIT1)
          ENDDO
        ELSE
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),NBITS),IBIT1)
            IND(2,INTS) = IBTAND(       IBUF(INTS,1),       IBIT1)
            IND(3,INTS) = IBTAND(IBTSHR(IBUF(INTS,2),NBITS),IBIT1)
            IND(4,INTS) = IBTAND(       IBUF(INTS,2),       IBIT1)
          ENDDO
        END IF
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJKL: Integrals'
        WRITE(LUPRI,'(4I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),IND(3,I),IND(4,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,D,I)
              FMAT(C,D,I) = FMAT(C,D,I) + DINT*DMAT(A,B,I)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,D,I)
              FMAT(C,D,I) = FMAT(C,D,I) + DINT*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKL: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              GCD  = D2*BUF(INTS)*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GCD
              FMAT(B,A,I) = FMAT(B,A,I) + GCD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I) = FMAT(C,D,I) + GAB
              FMAT(D,C,I) = FMAT(D,C,I) + GAB
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              GCD  = D2*BUF(INTS)*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GCD
              FMAT(B,A,I) = FMAT(B,A,I) + GCD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I) = FMAT(C,D,I) + GAB
              FMAT(D,C,I) = FMAT(D,C,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKL: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJKL ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJKL: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      GOTO 10
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijil */
      SUBROUTINE DOIJIL(IUNIT,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &                  IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJIL
C     Indices are bit packed with the degenerate index as index 3.
C
C     Written by T.Saue Mar 4 1995
C     Last revision Mar 4 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D2 = 2.0D0,D4 = 4.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT) 
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT) 
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(N2BBASX),IBUF(N2BBASX,NIBUF),IND(3,N2BBASX)
      INTEGER A,B,D
#include "ibtfun.h"
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJIL: IFCTYP:'
        WRITE(LUPRI,'(5X,I5)') (IFCTYP(I),I=1,NDMAT)
      ENDIF        
   10   CONTINUE
        CALL REAINT(IUNIT,NBUF,IREPAB,N2BBASX,NIBUF,IBUF,BUF)
        IF(NBUF.EQ.-1) RETURN
C
C       Unpack indices
C       ==============
C
        IF (NIBUF .EQ. 1) THEN
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),3*NBITS),IBIT1)
            IND(2,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),2*NBITS),IBIT1)
            IND(3,INTS) = IBTAND(       IBUF(INTS,1),         IBIT1)
          ENDDO
        ELSE
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),NBITS),IBIT1)
            IND(2,INTS) = IBTAND(       IBUF(INTS,1),       IBIT1)
            IND(3,INTS) = IBUF(INTS,2)
          ENDDO
        END IF
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJIL: Integrals'
        WRITE(LUPRI,'(3I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),IND(3,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(A,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) + DINT*DMAT(A,B,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(A,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) + DINT*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJIL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJIL: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              GAD  = D2*BUF(INTS)*(DMAT(A,D,I) + DMAT(D,A,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GAD
              FMAT(B,A,I) = FMAT(B,A,I) + GAD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(A,D,I) = FMAT(A,D,I) + GAB
              FMAT(D,A,I) = FMAT(D,A,I) + GAB
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              GAD  = D2*BUF(INTS)*(DMAT(A,D,I) + DMAT(D,A,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GAD
              FMAT(B,A,I) = FMAT(B,A,I) + GAD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(A,D,I) = FMAT(A,D,I) + GAB
              FMAT(D,A,I) = FMAT(D,A,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJIL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJIL: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJIL ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJIL: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      GOTO 10
CC
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijkk */
      SUBROUTINE DOIJKK(IUNIT,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &                  IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJKK
C     Indices are bit packed with the degenerate index as index 3.
C
C     Written by T.Saue Mar 4 1995
C     Last revision Mar 4 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(N2BBASX),IBUF(N2BBASX,NIBUF),IND(3,N2BBASX)
      INTEGER A,B,C
#include "ibtfun.h"
C
 10     CONTINUE
        CALL REAINT(IUNIT,NBUF,IREPAB,N2BBASX,NIBUF,IBUF,BUF)
        IF(NBUF.EQ.-1) RETURN
C
C       Unpack indices
C       ==============
C
        IF (NIBUF .EQ. 1) THEN
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),3*NBITS),IBIT1)
            IND(2,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),2*NBITS),IBIT1)
            IND(3,INTS) = IBTAND(       IBUF(INTS,1),         IBIT1)
          ENDDO
        ELSE
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),NBITS),IBIT1)
            IND(2,INTS) = IBTAND(       IBUF(INTS,1),       IBIT1)
            IND(3,INTS) = IBUF(INTS,2)
          ENDDO
        END IF
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJKK: Integrals'
        WRITE(LUPRI,'(3I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),IND(3,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              DINT = D2*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,C,I)
              FMAT(C,C,I) = FMAT(C,C,I) + DINT*DMAT(A,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (1/4) (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              DINT = D2*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,C,I)
              FMAT(C,C,I) = FMAT(C,C,I) + DINT*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKK: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              GCC  = D2*BUF(INTS)*DMAT(C,C,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GCC
              FMAT(B,A,I) = FMAT(B,A,I) + GCC
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,C,I) = FMAT(C,C,I) + GAB
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              GCC  = D2*BUF(INTS)*DMAT(C,C,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GCC
              FMAT(B,A,I) = FMAT(B,A,I) + GCC
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,C,I) = FMAT(C,C,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKK: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJKK ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJKK: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      GOTO 10
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijjj */
      SUBROUTINE DOIJJJ(IUNIT,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJJJ
C
C     Written by T.Saue Mar 5 1995
C     Last revision: tsaue Mar 5 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#include "dgroup.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(NNDIM),IBUF(NNDIM),IND(2,NNDIM)
      INTEGER A,B
#include "ibtfun.h"
        CALL REAINT(IUNIT,NBUF,IREPAB,NNDIM,1,IBUF,BUF)
C
C       Unpack indices
C       ==============
C
        DO INTS = 1,NBUF
          IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS),NBITS),IBIT1)
          IND(2,INTS) = IBTAND(       IBUF(INTS),       IBIT1)
        ENDDO
C
C       Print integrals
C       ===============
C
        IF(IPRINT.GE.6) THEN
          WRITE(LUPRI,'(A)') '* DOIJJJ: Integrals'
          IF(NBUF.GT.0) WRITE(LUPRI,'(2I5,3X,F12.6)')
     &       (IND(1,I),IND(2,I),BUF(I),I=1,NBUF)
        ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + BUF(INTS)*DMAT(B,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,B,I) = FMAT(B,B,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(B,B,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
           ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              DINT = D2*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + DINT*DMAT(B,A,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJJJ ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJJJ: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              GBB  = BUF(INTS)*DMAT(B,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GBB
              FMAT(B,A,I) = FMAT(B,A,I) + GBB
              FMAT(B,B,I) = FMAT(B,B,I) 
     &                    + BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,B,I) = FMAT(B,B,I) 
     &                    - BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              GBB  = BUF(INTS)*DMAT(B,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) - GBB
              FMAT(B,A,I) = FMAT(B,A,I) - GBB
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,B,I) = FMAT(B,B,I) 
     &                    + D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              GBB  = D2*BUF(INTS)*DMAT(B,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GBB
              FMAT(B,A,I) = FMAT(B,A,I) + GBB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJJJ ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJJJ: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJJJ ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJJJ: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doiikk */
      SUBROUTINE DOIIKK(IUNIT,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IIKK
C
C     Written by T.Saue Mar 5 1995
C     Last revision: tsaue Mar 5 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(NNDIM),IBUF(NNDIM),IND(2,NNDIM)
      INTEGER A,B
#include "ibtfun.h"
C
      CALL REAINT(IUNIT,NBUF,IREPAB,NNDIM,1,IBUF,BUF)
C
C     Unpack indices
C     ==============
C
      DO INTS = 1,NBUF
        IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS),NBITS),IBIT1)
        IND(2,INTS) = IBTAND(       IBUF(INTS),       IBIT1)
      ENDDO
C
C     Print integrals
C     ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIIKK: Integrals'
        WRITE(LUPRI,'(2I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
       DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + BUF(INTS)*DMAT(A,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIKK: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + D2*BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + D2*BUF(INTS)*DMAT(A,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + D2*BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + D2*BUF(INTS)*DMAT(A,A,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIKK: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIIKK ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIIKK: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijij */
      SUBROUTINE DOIJIJ(IUNIT,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJIJ
C
C     Written by T.Saue Mar 5 1995
C     Sep 3 1996 - tsaue   Generalize to IFCTYP
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0,D3 = 3.0D0,DP5 = 0.5D0)
#include "dcbbas.h"
#include "dcbind.h"
#include "dgroup.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(NNDIM),IBUF(NNDIM),IND(2,NNDIM)
      INTEGER A,B
#include "ibtfun.h"
C
      CALL REAINT(IUNIT,NBUF,IREPAB,NNDIM,1,IBUF,BUF)
C
C     Unpack indices
C     ==============
C
      DO INTS = 1,NBUF
        IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS),NBITS),IBIT1)
        IND(2,INTS) = IBTAND(       IBUF(INTS),       IBIT1)
      ENDDO
C
C     Print integrals
C     ===============
C
      IF(IPRINT.GE.6.AND.NBUF.GT.0) THEN
        WRITE(LUPRI,'(A)') '* DOIJIJ: Integrals'
        WRITE(LUPRI,'(2I5,3X,F12.6)')
     &   (IND(1,I),IND(2,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C       Symmetric singlet Fock matrix
C       =============================
C       F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + D3*BUF(INTS)*DMAT(A,B,I)
              EINT = -DP5*BUF(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + EINT*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + EINT*DMAT(A,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(B,A,I)
              EINT = -DP5*BUF(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + EINT*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + EINT*DMAT(A,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + D4*BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJIJ ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
            CALL QUIT(
     &      'ERROR in DOIJIJ: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
#include"scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              GAB  = BUF(INTS)*DMAT(A,B,I)
              GBA  = BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) + D2*GAB + GBA
              FMAT(B,A,I) = FMAT(A,B,I) + D2*GBA + GAB
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
#include "scalar.h"
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I)+DMAT(B,A,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GAB
              FMAT(B,A,I) = FMAT(A,B,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJIJ ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &     'ERROR in DOIJIJ: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &    'DOIJIJ ERROR, specified IFCTYP not implemented yet',
     &    '              specified IFCTYP was',IFCTYP(I),
     &    '              for F,D matrix no.  ',I
          CALL QUIT(
     &    'ERROR in DOIJIJ: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doiiii */
      SUBROUTINE DOIIII(IUNIT,NDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IND,BUF,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IIII
C
C     Written by T.Saue Mar 5 1995
C     Sep 3 1996 - tsaue   Genrazlized to IFCTYP
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0,DP5=0.5D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION BUF(NDIM),IND(NDIM),IFCTYP(NDMAT),
     &          IREPDM(NDMAT)
      INTEGER A
      CALL REAINT(IUNIT,NBUF,IREPAB,NDIM,1,IND,BUF)
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIIII: Integrals'
        WRITE(LUPRI,'(I5,3X,F12.6)')
     &     (IND(I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + DP5*BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - DP5*BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(A,A,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIII ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIII: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + D2*BUF(INTS)*DMAT(A,A,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIII ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIII: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIIII ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIIII: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirscf */
      SUBROUTINE DIRSCF(FMAT,DSO,NFMAT,ISYMOP,IHRMOP,IFCKOP,
     &                  NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations
C
C     Written by T.Saue July 2 1995
C     tsaue -    Sep 2 1996 Generalized to NFMAT 
C                           and general operator symmetry IRPOP
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
#include "maxorb.h"
#include "aovec.h"
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),DSO(N2BBASX,NZ,NFMAT),
     &          ISYMOP(*),IHRMOP(*),IFCKOP(*),
     &          NPOS(*),WORK(LWORK)
#include "dgroup.h"
#include "dcbqal.h"
#include "dcbbas.h"
#include "blocks.h"
#include "cbihr2.h"
      EXTERNAL BDQALGEBRA
#include "ibtfun.h"
C
      CALL QENTER('DIRSCF')
#include "memint.h"
C
C     Memory allocation
      NDMAT   = NFMAT*4
      N2GAB   = NSYMBL*NSYMBL
      NDMRAO  = N2GAB*NDMAT
      CALL MEMGET('INTE',KIFC  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRD  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAB  ,N2GAB  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
C
C     Get G(ab) integrals from AOPROPER
C     =================================
C
      IF (SCRFCK .GT. D0) THEN
C
         IJOB   = 0
         ITYPE  = 0
         IGTYP  = 1
         MAXDIF = 0
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDIF,
     &        IPRINT,WORK(KGAB),WORK(KFREE),LFREE)
C
      END IF
C
      IF(NBSYM.GT.1) THEN
        CALL MEMGET('REAL',KFAO,N2BBASX*NDMAT,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KDAO,N2BBASX*NDMAT,WORK,KFREE,LFREE)
C
C       Transform density matrix from QO to AO
C       ======================================
C
        IOFF = 0
        DO I = 1,NFMAT
          IREP = ISYMOP(I)-1
          IPAR = JBTOF(IREP,1)
          DO IZ = 1,4
            IREPD = IRQMAT(IZ,IREP)
            IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)            
            IPQ   = IQTOPQ(IQ,IREP)
            IF(IPQ.GT.NZ) CALL QUIT('DIRSCF: IPQ.GT.NZ !')
            IF(IPQ.LE.0 ) CALL QUIT('DIRSCF: IPQ.LE.0 !')
            CALL DTSOAO(DSO(1,IPQ,I),WORK(KDAO+IOFF),
     &                 NTBAS(0),IREPD,IPRTWO)
            IOFF = IOFF + N2BBASX
          ENDDO
          IF(IPRINT.GE.5) THEN
            JOFF = KDAO + (I-1)*4*N2BBASX
            WRITE(LUPRI,'(A,I5)') '** Matrix no.',I
            CALL HEADER('Backtransformed hph. density matrix',-1)
            CALL PRQMAT(WORK(JOFF),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),4,IQMULT,LUPRI)
          ENDIF
        ENDDO
C
C       Construct Fock matrix
C       ======================
C
        CALL DZERO(WORK(KFAO),N2BBASX*NDMAT)
        CALL DIRSC1(WORK(KFAO),WORK(KDAO),NFMAT,
     &              WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &              WORK(KDMRAO),WORK(KGAB),NPOS,INTFLG,
     &              IPRINT,WORK(KFREE),LFREE)
C
C       Transform Fock matrix from AO to QO
C       ===================================
C
        IOFF = 0
        DO I = 1,NFMAT
          IREP = ISYMOP(I)-1
          IPAR = JBTOF(IREP,1)
          DO IZ = 1,4
            IREPD = IRQMAT(IZ,IREP)
            IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)            
            IPQ   = IQTOPQ(IQ,IREP)
            CALL MTAOSO(WORK(KFAO+IOFF),FMAT(1,IPQ,I),
     &                 NTBAS(0),IREPD,IPRTWO)
            IOFF = IOFF + N2BBASX
          ENDDO
        ENDDO
      ELSE
C
C       If there is no symmetry involved the density matrix DSO
C       is sent to DIRSC1 module:
C
        CALL DIRSC1(FMAT,DSO,NFMAT,
     &              WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &              WORK(KDMRAO),WORK(KGAB),NPOS,INTFLG,
     &              IPRINT,WORK(KFREE),LFREE)
      ENDIF
      CALL MEMREL('DIRSCF.DIRSC1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('DIRSCF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirsc1 */
      SUBROUTINE DIRSC1(FAO,DAO,NFMAT,
     &                  IFCTYP,IREPDM,ISYMOP,IHRMOP,IFCKOP,
     &                  DMRAO,GABRAO,NPOS,INTFLG,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2 = 2.0D0,D10 = 10.0D0)
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
C
C Used from COMMON blocks:
C   CCOM  : THRS,NHTYP,KHK
C   DHFIOU: LUINTA
C   SHELLS: KMAX
C   DCBBAS: NTBAS(0),N2BBASX
C   DCBGEN: PARCAL
C
#include "ccom.h"
#include "dcbgen.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "shells.h"
#include "blocks.h"
#include "dcbqal.h"
#if defined (VAR_MPI)
#include "infpar.h"
#endif
      LOGICAL PERTUR, EXPECT, UNDIFF, SPNORB, DIRFCK, DISTRI, SQ12EL,
     &        LONDON, SUSCEP, DDFOCK, ADISTR, LBIT, FIRST
      CHARACTER SECTID*12, TOTTID*12, WALLTID*12, CPUTID*12
      DIMENSION FAO(N2BBASX*4,NFMAT),DAO(*),
     &          IFCTYP(4,NFMAT),IREPDM(4,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),
     &          NPOS(*),DMRAO(*),GABRAO(*),WORK(LWORK),DNTSKP(8)
      EXTERNAL BDQALGEBRA
C
      CALL QENTER('DIRSC1')
#include "memint.h"
C
Ctec Is FIRST needed???
      FIRST = .TRUE.
      NDMAT = NFMAT*4
      CALL SETFCK(IFCTYP,IREPDM,NFMAT,4,ISYMOP,IHRMOP,IFCKOP,IPRINT)
      CALL DZERO(DNTSKP,8)
      IF (SCRFCK.GT.D0) THEN
C
C       Make reduced density matrix for screening
C       ==========================================
C       (Reduced GAB is already present in GABRAO)
C
        CALL MKDRAO(DAO,DMRAO,NDMAT,WORK,LWORK,IPRINT)
      ENDIF
C
C     ****************************
C     **** Process integrals *****
C     ****************************
C
      IRNTYP = 3
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      KPOS   = 1
      IF (PARCAL) THEN
         CALL SCRSTA('MPIHeader',DUM,DUM)
      ELSE
         CALL SCRSTA('Header',DUM,DUM)
      END IF
C
C initialize MTOTTK
C
C     LL - integrals
C
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      IF(LBIT(INTFLG,1).AND.LBIT(IDFLAG,1)) THEN
        I2TYP = 1
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( 1 )
C
           CALL GETTIM(CPU1,WALL1)
           CALL HER_PARDRV(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DUM,DNTSKP)
C
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA('AOfock:LL',DNTSKP,WALL)
           IF (IPRINT .GT. 0) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &             TOTTID
              IF (IPRINT .GT. 2) THEN
                 CPUTID = SECTID(CPU)
                 WALLTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &             CPUTID,'/',WALLTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DUM,DNTSKP,.TRUE.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA('AOfock:LL',DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
        IF(IPRINT.GE.5) THEN
           DO I = 1,NFMAT
              CALL HEADER('Skeleton LL-FOCK',-1)
              WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
              CALL PRQMAT(FAO(1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0
     &             ),4,IQMULT,LUPRI)
           ENDDO
        ENDIF
      ENDIF
      KPOS = KPOS + MTOTTK
C
C     SL - integrals
C
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      IF(LBIT(INTFLG,2).AND.LBIT(IDFLAG,2)) THEN
        I2TYP  = 2
        THRBUF = THRS
        THRS   = THRS*THRFAC(1)
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( 1 )
C
           CALL GETTIM(CPU1,WALL1)
           CALL HER_PARDRV(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DUM,DNTSKP)
C
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA('AOfock:SL',DNTSKP,WALL)
           IF (IPRINT .GT. 0) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &         '>>>> Total wall time used in HER_PARDRV so far   :',
     &         TOTTID
              IF (IPRINT .GT. 2) THEN
                 CPUTID = SECTID(CPU)
                 WALLTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &         CPUTID,'/',WALLTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DUM,DNTSKP,.TRUE.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA('AOfock:SL',DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
        IF(IPRINT.GE.5) THEN
          DO I = 1,NFMAT
            CALL HEADER('Skeleton LL+SL-FOCK',-1)
            WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
            CALL PRQMAT(FAO(1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  4,IQMULT,LUPRI)
          ENDDO
        ENDIF
        THRS = THRBUF
        CALL FLSHFO(LUPRI)
      ENDIF
      KPOS = KPOS + MTOTTK
C
C     SS - integrals
C
      MTOTTK = NSMLBL*(NSMLBL+1)/2
      IF(LBIT(INTFLG,3).AND.LBIT(IDFLAG,3)) THEN
        I2TYP  = 3
        THRBUF = THRS
        THRS   = THRS*THRFAC(2)
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( 1 )
C
           CALL GETTIM(CPU1,WALL1)
           CALL HER_PARDRV(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DUM,DNTSKP)
C
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA('AOfock:SS',DNTSKP,WALL)
           IF (IPRINT .GT. 0) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &             TOTTID
              IF (IPRINT .GT. 2) THEN
                 CPUTID = SECTID(CPU)
                 WALLTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &         CPUTID,'/',WALLTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DUM,DNTSKP,.TRUE.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA('AOfock:SS',DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
        IF(IPRINT.GE.5) THEN
          DO I = 1,NFMAT
            CALL HEADER('Skeleton LL+SL+SS-FOCK',-1)
            WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
            CALL PRQMAT(FAO(1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  4,IQMULT,LUPRI)
          ENDDO
        ENDIF
        CALL FLSHFO(LUPRI)
        THRS = THRBUF
      ENDIF
      KPOS = KPOS + MTOTTK
C
C     Symmetrize Fock matrix
C     ======================
C
      CALL SKLFCK(FAO,WORK,LWORK,IPRINT,.TRUE.,.FALSE.,.FALSE.,
     &            .FALSE.,.TRUE.,0,.FALSE.,NDMAT,IREPDM,IFCTYP,IDUM)
      IF(IPRINT.GE.5) THEN
        DO I = 1,NFMAT
          CALL HEADER('Symmetrized skeleton LL+SL+SS-FOCK',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(FAO(1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                4,IQMULT,LUPRI)
        ENDDO
      ENDIF
C
      CALL QEXIT('DIRSC1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck conscf */
      SUBROUTINE CONSCF(FSO,DSO,NFMAT,
     &                  IFCKOP,ISYMOP,IHRMOP,INTFLG,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for conventional SCF
C
C     Written by T.Saue July 2 1995
C     Sep 2 1996 tsaue - Allocate buffer IND for unpacking of indices
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C   DHFIOU: LUINT
C   INDPCK: NIBUF
C
#include "dcbind.h"
#include "dcbgen.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION FSO(*),DSO(*),
     &          IFCKOP(NFMAT),ISYMOP(NFMAT),IHRMOP(NFMAT),
     &          WORK(LWORK)
C
#include "memint.h"
C
C     Memory allocation
C
      NDMAT = NZ*NFMAT
      CALL MEMGET('INTE',KIFC  ,NDMAT        ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRD  ,NDMAT        ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF  ,N2BBASX      ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBUF ,N2BBASX*NIBUF,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIND  ,N2BBASX*4    ,WORK,KFREE,LFREE)
      CALL CONSC1(FSO,DSO,NFMAT,
     &            WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &            INTFLG,WORK(KBUF),WORK(KIBUF),
     &            WORK(KIND),IPRINT,WORK(KFREE),LFREE)
C
C     Memory deallocation
C
      CALL MEMREL('CONSCF',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CONSC1 */
      SUBROUTINE CONSC1(FSO,DSO,NFMAT,
     &                  IFCTYP,IREPDM,ISYMOP,IHRMOP,IFCKOP,
     &                  INTFLG,BUF,IBUF,IND,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Conventional SCF
C
C     Written by T.Saue Feb 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbind.h"
#include "dcbgen.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      LOGICAL LBIT
      DIMENSION FSO(N2BBASX,*),DSO(N2BBASX,*),
     &          IFCTYP(NZ,NFMAT),IREPDM(NZ,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),
     &          BUF(*),IBUF(*),IND(*),
     &          WORK(LWORK)
C
      LUINTA = 14
      LUINTB = 15
      NDMAT  = NZ*NFMAT
      CALL SETFCK(IFCTYP,IREPDM,NFMAT,NZ,ISYMOP,IHRMOP,IFCKOP,IPRINT)
C
C     LL-integrals
C     ============
C
      IF(LBIT(INTFLG,1).AND.(.NOT.LBIT(IDFLAG,1))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: LL-integrals:'
        CALL TIMER('START ',TIMSTR,TIMEND)
        NNDIM = NTBAS(1)*(NTBAS(1)-1)/2
        CALL OPNFIL(LUINTA,'DFLLSA','OLD','LL2SMT')
        CALL OPNFIL(LUINTB,'DFLLSB','OLD','LL2SMT')
        CALL XX2MAT(LUINTA,LUINTB,NTBAS(1),NNDIM,
     &       FSO,DSO,NDMAT,IFCTYP,IREPDM,
     &       BUF,IBUF,IND,IPRINT)
        CLOSE(LUINTA,STATUS='KEEP')
        CLOSE(LUINTB,STATUS='KEEP')
        CALL TIMER('CONSCF(LL)',TIMSTR,TIMEND)
      ENDIF
C
C     SL-integrals
C     ============
C
      IF(LBIT(INTFLG,2).AND.(.NOT.LBIT(IDFLAG,2))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: SL-integrals:'
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL OPNFIL(LUINTA,'DFSLSA','OLD','SL2SMT')
        CALL OPNFIL(LUINTB,'DFSLSB','OLD','SL2SMT')
        CALL XY2MAT(LUINTA,LUINTB,FSO,DSO,NDMAT,IFCTYP,IREPDM,
     &              BUF,IBUF,IND,IPRINT)
        CLOSE(LUINTA,STATUS='KEEP')
        CLOSE(LUINTB,STATUS='KEEP')
        CALL TIMER('CONSCF(SL)',TIMSTR,TIMEND)
      ENDIF
C
C     SS-integrals
C     ============
C
      IF(LBIT(INTFLG,3).AND.(.NOT.LBIT(IDFLAG,3))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: SS-integrals:'
        CALL TIMER('START ',TIMSTR,TIMEND)
        NNDIM = NTBAS(2)*(NTBAS(2)-1)/2
        CALL OPNFIL(LUINTA,'DFSSSA','OLD','SS2TMT')
        CALL OPNFIL(LUINTB,'DFSSSB','OLD','SS2SMT')
        CALL XX2MAT(LUINTA,LUINTB,NTBAS(2),NNDIM,
     &      FSO,DSO,NDMAT,IFCTYP,IREPDM,
     &      BUF,IBUF,IND,IPRINT)
        CLOSE(LUINTA,STATUS='KEEP')
        CLOSE(LUINTB,STATUS='KEEP')
        CALL TIMER('CONSCF(SS)',TIMSTR,TIMEND)
      ENDIF
C
C     Gaunt-integrals
C     ===============
C     to be programmed
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C*/Deck setfck */
      SUBROUTINE SETFCK(IFCTYP,IREPDM,NFMAT,MZ,ISYMOP,IHRMOP,IFCKOP,
     &                  IPRINT)
C*****************************************************************************
C
C     This routine will assign IFCTYP, indicating Fock matrices
C     to construct based on information about
C       hermiticiy of matrices     IHRMOP 0 - no symmetry
C                                         1 - Hermitian
C                                        -1 - anti-Hermitian
C       point group symmetry       ISYMOP 
C       Fock matrix type           IFCKOP 1, 2, or 3 - see IPFCK below
C
C     IFCTYP = XY
C       X - symmetry about the diagonal  - based on IPHRM
C         0 - G: (General) no symmetry
C         1 - S: Symmetric
C         2 - A: Anti-symmetric
C         3 - M: (Mixed) some blocks symmetric, other blocks anti-symmetric
C       Y - two-electron contributions   - based on IPFCK
C         1 - Coulomb
C         2 - Exchange
C         3 - Coulomb+Exchange
C
C     For the moment we will assume closed-shell Fock-matrices
C
C     Written by T.Saue Sep 4 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"    
C
#include "dcbqal.h"
#include "pgroup.h"
#include "dgroup.h"
#include "dcbham.h"
C     External variables
      EXTERNAL BDQALGEBRA
      DIMENSION IFCTYP(MZ,NFMAT),IREPDM(MZ,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT)
C     Local variables
      PARAMETER (D0 = 0.0D0)
      SAVE IPFCK,HTYP
      DIMENSION IPFCK(4,3),IRFCK(4),IQFCK(4)
      CHARACTER HTYP(0:3)*1
      DATA HTYP /'G','S','A','M'/
      DATA IPFCK/3,2,2,2, 1,0,0,0, 2,2,2,2/
#include "ibtfun.h"
C
      IF (IPRINT.GE.2) THEN
        CALL HEADER('Output from SETFCK',-1)
        WRITE(LUPRI,'(A,F7.4)') '* Exchange factor: ',HFXFAC
      ENDIF
      CALL IZERO(IFCTYP,MZ*NFMAT)
      DO I = 1,NFMAT
        IREP = ISYMOP(I)-1
C       IREP is always 0 (totally symmetric) for Hartre-Fock
C       (defined in SETDHF(dirac/dirset.F) and GETFC1(prp/pamset.F))
C       and it is the irrep of the perturbation operator in response
C       (defined in SETGMT(prp/pamrvc.F)). /hjaaj Apr 2002
        IHRM = IHRMOP(I)
        IPAR = JBTOF(IREP,1)
        IFCK = IFCKOP(I)
C
C       First determine contributions
C
        DO IZ = 1,4
          IRFCK(IZ) = IRQMAT(IZ,IREP)
          IY        = IPFCK(IZ,IFCK)
          IC        = MOD(IY,2)
          IE        = IY - IC
          IH        = IHQMAT(IZ,IHRM)
C         An anti-symmetric matrix gives no Coulomb contributions
          IF (IH.EQ.2) IC = 0
C         No exchange contributions if HFXFAC = 0
          IF (HFXFAC.EQ.D0) IE = 0
          IY   = IC + IE
          IF (MZ.LT.4) THEN
C           ... we are using .SOFOCK or non-direct (CONSCF)
            IQFCK(IZ)     = IQMULT(1,JQBAS(IRFCK(IZ),IPAR),IZ)
            IPQ           = IQTOPQ(IQFCK(IZ),IREP)
C
            IYOLD = MOD(IFCTYP(IPQ,I),10)
            IXOLD = IFCTYP(IPQ,I)/10
            IFCTYP(IPQ,I) = 10*IBTOR(IXOLD,IH) + IBTOR(IYOLD,IY)
C           ... hjaaj apr 2002: with IBTOR we are sure we include
C           all the requested types. If some blocks are symmetric
C           and some are antisymmetric we end with the new type
C           of IX = 3 (IFCTYP = 10*IX + IY), which is in a way a
C           special case of IX = 0 (unsymmetric). /hjaaj Apr 2002
C
            IREPDM(IPQ,I) = IREP
C           ... MZ.lt.4 means that SOFOCK has been requested, and
C           the density matrix will be of mixed symmetry. We
C           use IREPDM to tell FCKOUT(abacus/her2out.F) which
C           irrep gives symmetric matrix elements and thus
C           Coulomb contributions. /hjaaj Apr 2002
C
          ELSE
C           ... MZ.eq.4: we are using DIRSCF (not .SOFOCK)
            IQFCK(IZ)     = IZ
            IFCTYP(IZ,I)  = 10*IH+IY
            IREPDM(IZ,I) =  IRFCK(IZ)
          ENDIF
        ENDDO
C
C     Print section
C
        IF(IPRINT.GE.2) THEN
          WRITE(LUPRI,'(A,I5,3X,A,3X,A,A3)') 
     &      '* Fock matrix no.',I,' - ',
     &      'Operator irrep: ',REP(IREP)
          WRITE(LUPRI,'(A,4(2X,A3))')
     &      'Matrix symmetry:',(REP(IRFCK(IZ)),IZ=1,4)
          IF(MZ.LT.4) THEN
            WRITE(LUPRI,'(A,4(3X,A1,1X))')
     &        'Packing        :',(QUNIT(IQFCK(IZ)),IZ=1,4)
            WRITE(LUPRI,'(A,4(3X,A1,1X))')
     &        'Eff.   packing :',
     &        (QUNIT(IQTOPQ(IQFCK(IZ),IREP)),IZ=1,4)
          ENDIF
            WRITE(LUPRI,'(A)') 
     &      '    Mat.  IFC  SYM  Coulomb Exchange Irrep'
          DO IZ = 1,MZ
            IY = MOD(IFCTYP(IZ,I),10)
            IH = (IFCTYP(IZ,I)-IY)/10
            IC = MOD(IY,2)
            IE = IY - IC
            IQ = IZ
            IF(MZ.LT.4) IQ = IPQTOQ(IZ,IREP)
            WRITE(LUPRI,'(4X,I1,3A1,I5,3X,A1,6X,I1,7X,I1,6X,A3)')
     &        IZ,'(',QUNIT(IQ),')',
     &        IFCTYP(IZ,I),HTYP(IH),IC,IE/2,REP(IREPDM(IZ,I))
          ENDDO
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirsof */
      SUBROUTINE DIRSOF(FMAT,DSO,NFMAT,ISYMOP,IHRMOP,IFCKOP,
     &                  NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations in SO-basis
C
C     Written by T.Saue August 28 1996
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
#include "maxorb.h"
#include "aovec.h"
Ctec Is FIRST needed??
      LOGICAL FIRST
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),DSO(N2BBASX,NZ,NFMAT),
     &          IFCKOP(NFMAT),NPOS(*),
     &          IHRMOP(NFMAT),ISYMOP(NFMAT),WORK(LWORK)
#include "dgroup.h"
#include "cbihr2.h"
#include "dcbbas.h"
#include "blocks.h"
C
C      SAVE FIRST
      FIRST = .TRUE.
      CALL QENTER('DIRSOF')
#include "memint.h"
C
C     Memory allocation
C
      NDMAT  = NFMAT*NZ
      N2GAB  = NSYMBL*NSYMBL
      NDMRAO = N2GAB*NDMAT
      NDMRSO = MAXSHL*MAXSHL*NDMAT
      CALL MEMGET('INTE',KIFC  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRD  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGAB  ,N2GAB  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDMRSO,NDMRSO ,WORK,KFREE,LFREE)
C
C     Get G(ab) integrals from AOPROPER
C     =================================
C
      IF (SCRFCK .GT. D0) THEN
         IJOB   = 0
         ITYPE  = 0
         IGTYP  = 1
         MAXDIF = 0
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDIF,
     &        IPRINT,WORK(KGAB),WORK(KFREE),LFREE)
C
      END IF
C
      CALL DIRSO1(FMAT,DSO,NFMAT,
     &            WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &            WORK(KDMRAO),WORK(KDMRSO),WORK(KGAB),
     &            NPOS,INTFLG,IPRINT,WORK(KFREE),LFREE)
      CALL MEMREL('DIRSCF.DIRSC1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('DIRSOF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirso1 */
      SUBROUTINE DIRSO1(FSO,DSO,NFMAT,
     &                  IFCTYP,IREPDM,ISYMOP,IHRMOP,IFCKOP,
     &                  DMRAO,DMRSO,GABRAO,NPOS,INTFLG,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations in SO-basis
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2 = 2.0D0,D4 = 4.0D0,D8 = 8.0D0,D0=0.0D0)
#include "dummy.h"
#include "maxaqn.h"
C
#include "maxorb.h"
#include "aovec.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "ccom.h"
#include "blocks.h"
#if defined (VAR_MPI)
#include "infpar.h"
#endif
      LOGICAl LBIT
      CHARACTER SECTID*12, INTTID*12, TOTTID*12, CPUTID*12, WALLTID*12
      DIMENSION FSO(N2BBASX,*),DSO(N2BBASX,*),
     &          IFCTYP(*),IREPDM(*),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),
     &          DMRAO(*),DMRSO(*),GABRAO(*),
     &          DNTSKP(8),NPOS(*),WORK(LWORK)
C
      CALL QENTER('DIRSO1')
#include "memint.h"
C
C     Initialization
C
      NDMAT = NFMAT*NZ
      CALL DZERO(DNTSKP,8)
      CALL SETFCK(IFCTYP,IREPDM,NFMAT,NZ,ISYMOP,IHRMOP,IFCKOP,IPRINT)
C
C     Make reduced density matrix for screening
C     =========================================
      IF(SCRFCK.GT.D0) THEN
C       in SO-basis
        CALL MKDRSO(DSO,DMRSO,NDMAT,WORK,LWORK,IPRINT)
C       in AO-basis
        CALL MKDRQO(DSO,DMRAO,NFMAT,ISYMOP,WORK,LWORK,IPRINT)
      ENDIF
C
C     ****************************
C     **** Process integrals *****
C     ****************************
C
      IRNTYP = 9
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      KPOS   = 1
      IF (PARCAL) THEN
         CALL SCRSTA('MPIHeader',DUM,DUM)
      ELSE
         CALL SCRSTA('Header',DUM,DUM)
      END IF
C
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      IF(LBIT(INTFLG,1).AND.LBIT(IDFLAG,1)) THEN
        I2TYP = 1
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( 1 )
C
           CALL GETTIM(CPU1,WALL1)
           CALL HER_PARDRV(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DMRSO,DNTSKP)
C
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA('SOfock:LL',DNTSKP,WALL)
           IF (IPRINT .GT. 0) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &             TOTTID
              IF (IPRINT .GT. 5) THEN
                 CPUTID = SECTID(CPU)
                 WALLTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &         CPUTID,'/',WALLTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DMRSO,DNTSKP,.TRUE.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA('SOfock:LL',DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
      ENDIF
      KPOS = KPOS + MTOTTK
C
C     SL - integrals
C
      MTOTTK = NLRGBL*(NLRGBL+1)/2
      IF(LBIT(INTFLG,2).AND.LBIT(IDFLAG,2)) THEN
        I2TYP  = 2
        THRBUF = THRS
        THRS   = THRS*THRFAC(1)
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( 1 )
C
           CALL GETTIM(CPU1,WALL1)
           CALL HER_PARDRV(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DMRSO,DNTSKP)
C
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA('SOfock:SL',DNTSKP,WALL)
           IF (IPRINT .GT. 0) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &             TOTTID
              IF (IPRINT .GT. 2) THEN
                 CPUTID = SECTID(CPU)
                 WALLTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &         CPUTID,'/',WALLTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DMRSO,DNTSKP,.TRUE.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA('SOfock:SL',DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
        THRS = THRBUF
        CALL FLSHFO(LUPRI)
      ENDIF
      KPOS = KPOS + MTOTTK
C
C     SS - integrals
C
      MTOTTK = NSMLBL*(NSMLBL+1)/2
      IF(LBIT(INTFLG,3).AND.LBIT(IDFLAG,3)) THEN
        I2TYP  = 3
        THRBUF = THRS
        THRS   = THRS*THRFAC(2)
        CALL TIMER('START ',TIMSTR,TIMEND)
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( 1 )
C
           IF (IPRINT .GT. 0) CALL GETTIM(CPU1,WALL1)
C
           CALL HER_PARDRV(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DMRSO,DNTSKP)
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA('SOfock:SS',DNTSKP,WALL)
C
           IF (IPRINT .GT. 0) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &             TOTTID
C     
              IF (IPRINT .GT. 2) THEN
                 INTTID = SECTID(CPU)
                 WRITE(LUPRI,'(A,A12)')
     &            '>>>> CPU  time used in HER_PARDRV last iteration :',
     &            INTTID
                 INTTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12)')
     &            '>>>> Wall time used in HER_PARDRV last iteration :',
     &            INTTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DMRSO,DNTSKP,.TRUE.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA('SOfock:SS',DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
        CALL FLSHFO(LUPRI)
        THRS = THRBUF
      ENDIF
      KPOS = KPOS + MTOTTK
C     
C     Scale matrices
C
      DO I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IF(IX.EQ.0) THEN
          CALL DSCAL(N2BBASX,D2,FSO(1,I),1)
        ELSE
          CALL DSCAL(N2BBASX,D4,FSO(1,I),1)
        ENDIF
      ENDDO
      CALL QEXIT('DIRSO1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Mkdrqo */
      SUBROUTINE MKDRQO(DQO,DMRQO,NFMAT,ISYMOP,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Make reduced density matrices in AO-basis starting from
C     density matrices in QO-basis
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "dgroup.h"
#include "dcbqal.h"
#include "dcbbas.h"
#include "blocks.h"
      DIMENSION DQO(N2BBASX,NZ,NFMAT),
     &          DMRQO(NSYMBL,NSYMBL,NZ,NFMAT),
     &          ISYMOP(NFMAT),WORK(LWORK)
      EXTERNAL BDQALGEBRA      
C
#include "ibtfun.h"
      CALL QENTER('MKDRQO')
#include "memint.h"
C
C     Initialize
C
      NDMAT = NFMAT*NZ
      N2DRAO = NSYMBL*NSYMBL*NDMAT
      CALL DZERO(DMRQO,N2DRAO)
C
C     Make index array from AOs to blocks
C
      CALL MEMGET('INTE',KINRAO,NTBAS(0),WORK,KFREE,LFREE)
      CALL PTRDAO(WORK(KINRAO))
C
C     Make reduced AO matrix
C
      IF(NBSYM.GT.1) THEN
        CALL MEMGET('REAL',KDBF,N2BBASX,WORK,KFREE,LFREE)
        DO I = 1,NFMAT
          IREP = ISYMOP(I)-1
          IPAR = JBTOF(IREP,1)
          DO IZ = 1,4
            IREPDM = IRQMAT(IZ,IREP)
            IQ     = IQMULT(1,JQBAS(IREPDM,IPAR),IZ)
            IQP    = IQTOPQ(IQ,IREP)
            IF(IQP.GT.NZ) CALL QUIT('MKDRQO: IQP.GT.NZ !')
            IF(IQP.LE.0 ) CALL QUIT('MKDRQO: IQP.LE.0 !')
            CALL DTSOAO(DQO(1,IQP,I),WORK(KDBF),
     &                 NTBAS(0),IREPDM,IPRINT)
            CALL GATMAT(0,NTBAS(0),NTBAS(0),WORK(KDBF),
     &              DMRQO(1,1,IQP,I),WORK(KINRAO),NSYMBL)
          ENDDO
        ENDDO
      ELSE
        DO I = 1,NFMAT
          DO IZ = 1,NZ
            CALL GATMAT(0,NTBAS(0),NTBAS(0),DQO(1,IZ,I),
     &              DMRQO(1,1,IZ,I),WORK(KINRAO),NSYMBL)
          ENDDO
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) CALL PRDRAO(DMRQO,NDMAT,WORK(KINRAO))
C
      CALL MEMREL('MKDRQO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('MKDRQO')
      RETURN
      END
