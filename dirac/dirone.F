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
C VERSION : $Revision: 1.78.2.2 $
C DATE    : $Date: 2004/09/01 08:08:49 $
C FILE    : dirone.F
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onegen */
      SUBROUTINE ONEGEN(WORK,LWORK)
C*****************************************************************************
C
C     This routine goes through the list of one-electron integrals specified
C     in the COMMON block XCBCLS. The integrals will be put on the file
C     AOPROPER.
C
C     Written by T.Saue - May 1996
C     Last revision: May 19 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "mxcent.h"
      DIMENSION WORK(LWORK)
      CHARACTER COMB*4
C
#include "dcbgen.h"
#include "dcbprl.h"
#include "dcbham.h"
#include "dcbcls.h"
#include "cbihr1.h"
#include "cbiher.h"
C
C     Loop through list of class labels; calculate only the active ones(CLSCAL)
C     =========================================================================
C
      DO 10 I = 1,NPRPCLS
      IF(.NOT.CLSCAL(I)) GOTO 10
Chj     CLSCMB and PDOINT is now '+', '-', or '0',
Chj     corresponding to a factor of +1, -1, or 0 on that block.
        DO J = 1,4
          IF (CLSCMB(I)(J:J) .EQ. '0') THEN
             COMB(J:J) = 'F'
          ELSE
             COMB(J:J) = 'T'
          END IF
        END DO  
C
C     Nuclear attraction integrals (molecular field)
C
      IF    (CLSINT(I).EQ.'MOLFLD ') THEN
        IF (NOSMLV) THEN
           CLSCMB(I)(4:4) = '0'
           PDOINT(I)(4:4) = '0'
           COMB(4:4) = 'F'
        END IF
        CALL NUCATT(WORK,LWORK,COMB,IPRONE)
C
C     Contribution from beta-matrix
C
      ELSEIF(CLSINT(I).EQ.'BETAMAT') THEN
        CALL BETAMT(COMB,WORK,LWORK,IPRONE)
C
C     All other integrals
C
      ELSE
        CALL OP1GEN(CLSINT(I),COMB,IORDCL(I),WORK,LWORK,IPRONE)
      ENDIF
      CLSCAL(I) = .FALSE.
   10 CONTINUE
C
      CLOSE(LU1INT,STATUS='KEEP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck op1gen */
      SUBROUTINE OP1GEN(WORD,COMB,IORDER,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the generation of one-electron integrals
C     to be placed on file AOPROPER.
C
C     Written by Trond Saue - May 1996
C
C     Last revision May 19 1996
C
C  MI(+TEC): minor changes for LAO based RPA, dec.2002
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
      LOGICAL DOINT(4),TOFILE
      CHARACTER WORD*7,COMB*4,LABINT(5*MXCOOR)*8,MTFORM*6
      DIMENSION WORK(LWORK)
C
      CALL QENTER('OP1GEN')
#include "memint.h"
C     Write integrals to file
      TOFILE   = .TRUE.
      NCOMP    = 0
      MAXTYP = 5*MXCOOR
      IF ((WORD.EQ.'ELFGRDC').OR.(WORD.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE IF (WORD.EQ.'DSO    ') THEN
CMI ... to be fixed for the future LAO implementation...
         NINTAD = 9*NUCDEP*NUCDEP
      ELSE
         NINTAD = MAXTYP
      END IF
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
C
C     Memory allocation
C
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      NCOMP = 0
      IF (WORD.EQ.'S1MAGR ') THEN 
         TRIANG = .FALSE.
         MTFORM = 'SQUARE'
      ELSE
         TRIANG = .TRUE.
         MTFORM = 'TRIANG'
      END IF
      CALL PR1IN1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,DUMMY,NCOMP,TOFILE,MTFORM, 
     &            DOINT)
C
C     Memory deallocation
C
      CALL MEMREL('OP1GEN',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('OP1GEN')
C
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'OP1GEN: Class label:',WORD
      WRITE(LUPRI,'(A,A)') 'OP1GEN ERROR: Not able to read COMB =',COMB
      CALL QUIT('OP1GEN: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck betamt */
      SUBROUTINE BETAMT(COMB,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the generation of beta-matrix integrals
C     to be placed on file AOPROPER.
C
C     Written by Trond Saue - May 1996
C     Last revision May 19 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dummy.h"
      LOGICAL DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      DIMENSION WORK(LWORK)
#include "dcbgen.h"
#include "dcbbas.h"
C
      CALL QENTER('BETAMT')
#include "memint.h"
C
C     Initialization
C
      NCOMP    = 1
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
C
C     Memory allocation
C
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      CALL RGETINT('OVERLAP',WORK(KONEIN),DOINT,NCOMP,IDUMMY,IDUMMY,
     &            'TRIANG',.FALSE.,IPRINT,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Small component overlap integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
      ENDIF
C
C     Generate integral labels
C
      LABEL = 'BETAMAT '
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(WORK(KONEIN),NNBBASX,LABEL,RTNLBL,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('BETAMT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('BETAMT')
C
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'BETAMT ERROR: Not able to read COMB =',COMB
      CALL QUIT('BETAMT: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucatt */
      SUBROUTINE NUCATT(WORK,LWORK,COMB,IPRINT)
C*****************************************************************************
C
C     Generate nuclear attraction integrals
C
C     Written by T.Saue May 1995
C     Last revision May 19 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAl DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4,OMITVNUC(2)*4
      DIMENSION WORK(LWORK)
#include "dcbbas.h"
C
      CALL QENTER('NUCATT')
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
      OMITVNUC(1) = 'FFFF'
      OMITVNUC(2) = 'FFFF'
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,OMITVNUC,IPRINT)
      IF(IPRINT.GE.4) THEN
        CALL HEADER('Nuclear attraction integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),1,LUPRI)
      ENDIF
C
C     Generate integral labels
C
      LABEL = 'MOLFIELD'
      CALL GETDAT(RTNLBL(1),RTNLBL(2))
      RTNLBL(2)(1:2) = 'SY'
      WRITE(RTNLBL(2)(3:4),'(I2)') 1
      RTNLBL(2)(5:8) = COMB
C
C     Write integrals to file
C
      CALL WRTPRO(WORK(KONEIN),NNBBASX,LABEL,RTNLBL,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('NUCATT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('NUCATT')
      RETURN
 1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
      WRITE(LUPRI,'(A,A)') 'NUCATT: Not able to read COMB =',COMB
      CALL QUIT('NUCATT: Not able to read COMB')
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck nucat1 */
      SUBROUTINE NUCAT1(OP1INT,WORK,LWORK,DOINT,OMITVNUC,IPRINT)
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC(2)
      DIMENSION OP1INT(NNBBASX),WORK(LWORK)
#include "dcbbas.h"
#include "nuclei.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('NUCAT1')
#include "memint.h"
C
      CALL MEMGET('REAL',KSTHMA,3*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACIN,2*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
C
C     PROPTY = .FALSE.
C     MAXDIF = 0
C     DIFDIP = .FALSE.
C     SECDER = .FALSE.
C     NOPV   = .FALSE.
C     NOCV   = .FALSE.
C
      CALL ONEDR1(WORK(KSTHMA),DUMMY,DUMMY,WORK(KFACIN),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KNCENT),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),DUMMY,WORK(KFREE),LWRK,
     &            IPRINT,.FALSE.,0,.TRUE.,.TRUE.,.FALSE.,NTBAS(0),
     &            NNBBASX,NNBBAST,DOINT,OMITVNUC,.FALSE.)
C
C     Symmetry unpack nuclear attraction integrals
C
      KNUCS = KSTHMA + NNBBASX
      CALL SYMUPK(WORK(KNUCS),OP1INT,1,NNBBASX)
C
C     Memory deallocation
C
      CALL MEMREL('NUCAT1.ONEDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('NUCAT1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onefck */
      SUBROUTINE ONEFCK(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C     Generate one-electron Fock matrix, using the one-electron operators
C     specified in COMMON block CBIHAM.
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "cbihr1.h"
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('ONEFCK')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     ************************************************
C     ****** Construct one-electron Fock matrix ******
C     ************************************************
C
      CALL ONEFC1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK(KFREE),LFREE,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('ONEFCK.ONEFC1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('ONEFCK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck onefc1 */
      SUBROUTINE ONEFC1(FOCK,FTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C     Generate one-electron Fock matrix, using the one-electron operators
C     specified in COMMON block CBIHAM.
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
C
      LOGICAL FIRST(NZ),TOBE,BARNUC
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
#include "memint.h"
C
C     Initialize
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL LSET(NZ,.TRUE.,FIRST)
C
C     Nuclear attraction integrals
C     ============================
C
      IOP = IPMOLFLD
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK(KFREE),LFREE,IPRHAM)
      IF (SMLV1C) THEN
         CALL CATCORR(FTRI,DUMMY,DUMMY,WORK,LFREE,IPRHAM)
         CALL PR1CEX1 (FTRI,FIRST,FOCK,WORK,KFREE,LFREE,IPRHAM)
      ENDIF
C
C     Beta matrix
C     ===========
C
      IOP = IPBETAMT
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK(KFREE),LFREE,IPRHAM)
C
C     Kinetic energy
C     ==============
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK(KFREE),LFREE,IPRHAM)
C
C     Finite field: 
C     Loop through list of one-electron operators in CBIHAM
C     =====================================================
C
      DO I = 1,N1OPER
        IOP = IPR1OP(I)
        CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,
     &              WORK(KFREE),LFREE,IPRHAM)
      ENDDO
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,NZ,WORK(KFREE),LFREE)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER('One-electron Fock matrix',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(FOCK(I2BASX(I,I)+1,1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prpmao */
      SUBROUTINE PRPMAO(LU1INT,IOP,PMAT,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C     Generates lower triangular matrix of one-electron integrals
C     corresponding to property IOP of COMMON block XCBPRP
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,PTOL = 1.0D-12)
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "dcbxpr.h"
#include "dcbprl.h"
C
      CHARACTER MXFORM*6, FMT*6, RTNLBL(2)*8,JCHR*1
      LOGICAL FIRST(*), DOINT(4)
      DIMENSION PMAT(NNBBASX,*),OP1INT(NNBBASX),WORK(LWORK), IDOINT(4)
C
      EXTERNAL BDQALGEBRA
C
C
C     Extract information about property
C     ==================================
C
      IOPTYP = IPRPTYP(IOP)
      ITIMOP = IPRPTIM(IOP)
      IBRPOP = IPRPSYM(IOP)-1
      NCMP   = MCMP(IOPTYP)
      IF(IPRINT.GE.3) THEN
        CALL HEADER('Output from PRPMAO',-1)
        WRITE(LUPRI,'(A,A16)')
     &    'Generating property matrix for property ',PRPNAM(IOP)
      ENDIF
C
C     Loop over components
C     ====================
C
      DO 10 I = 1,NCMP
        PFAC  = FACPRP(I,IOP)
        IF(ABS(PFAC).LT.PTOL) GOTO 10
C
C       Extract information about component
C
C       SCALAR OPERATOR:
C       ----------------
C       KLBL - pointer to list of labels
C       ISYM - scalar operator is symmetric(1)/antisymmetric(2) about diagonal
C       IBREP - boson irrep of scalar operator
C
C       MATRIX OPERATOR:
C       ----------------
C
C       IM4  - pointer to matrix operator
C       IFAC - sign of coefficient in linear combination
C       IFAS - sign of matrix operator
C         Anti-symmetric scalar operators are provided with an imaginary phase 
C         in order to retain Hermiticity
C       IMAT - quaternion position of matrix operator
C
        KLBL  = IPRPLBL(I,IOP)
C
C       Determine what matrix the operator enters
C       =========================================
C
C       Is this a LL+SS operator (IC=1) or a LS+SL operator (IC=2) ?
        IF(PDOINT(KLBL)(2:3) .EQ. '00') THEN
          IC = 1
        ELSE IF((PDOINT(KLBL)(1:1) .EQ. '0') .AND.
     &          (PDOINT(KLBL)(4:4) .EQ. '0')) THEN 
          IC = 2
        ELSE
          CALL QUIT('PRPMAO: illegal PDOINT!')
        ENDIF
        IBRP = IPRLREP(KLBL)
        IH   = IPRLTYP(KLBL)
        IM4  = JM4(I,IOPTYP)
        IZ   = JM4POS(IM4)
        IQ   = IQMULT(1,JQBAS(IBRP,IC),IZ)
        IQP  = IQTOPQ(IQ,IBRPOP)
        IF(IQP.LE.0) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A)') 'IQP(',IQP,').LE.0'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          WRITE(LUPRI,'(A,7I5)')
     &        'IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP :',
     &         IM4,IBRP,IH,IZ,IQ,IBRPOP,IQP
          WRITE(LUPRI,'(A)') 'IQTOPQ(1:4,0:NBSYM-1):'
          DO IZ = 1,4
             WRITE(LUPRI,'(8I5)') (IQTOPQ(IZ,J),J=0,NBSYM-1)
          END DO
          WRITE(LUPRI,'(A)') 'IPQTOQ(1:NZ,0:NBSYM-1):'
          DO IZ = 1,NZ
             WRITE(LUPRI,'(8I5)') (IPQTOQ(IZ,J),J=0,NBSYM-1)
          END DO
          CALL QUIT('PRPMAO: IQP.LE.0 !!!')
        END IF
        IF(IQP.GT.NZ) THEN
          WRITE(LUPRI,'(/A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A)') 'IQP(',IQP,').GT.NZ(',
     &          NZ,')'
          CALL WRIXPR(0,IOP)
          WRITE(LUPRI,'(A,I5)') 'Q: Scalar :',JQBAS(IBRP,IC)
          WRITE(LUPRI,'(A,I5)') 'Q: Matrix :',IZ
          CALL QUIT('PRPMAO: IQP.GT.NZ !!!')
        ENDIF
        IF(IPRINT.GE.4) THEN
          WRITE(LUPRI,'(A,I5,3X,A8)') 'Component nr.',I,PRPLBL(KLBL)
          FMT = MXFORM(PFAC,20)
          WRITE(LUPRI,'(A,4X,'//FMT//',A,A4)')
     &      'Factor: ',PFAC,' and PDOINT:',PDOINT(KLBL)
          WRITE(LUPRI,'(A,I5)') 'Integrals enter matrix :',IQP
          WRITE(LUPRI,'(A,L1)') 'FIRST = ',FIRST(IQP)
        ENDIF
C
C       Get integral from file
C
        CALL PRPREA(LU1INT,PRPLBL(KLBL),RTNLBL,PDOINT(KLBL),
     &              OP1INT,NNBBASX,IPRINT)
        IF(IPRINT.GE.8) THEN
           WRITE(LUPRI,'(A)') '* Integrals from AOPROPER'
           CALL OUTPAK(OP1INT,NTBAS(0),1,LUPRI)
        ENDIF        
C
C       Determine factor to scale matrix with
C       =====================================
C
        IFAC  = JCOM(I,IOPTYP)*JM4FAS(IM4)
        PFAC  = FACPRP(I,IOP)*IFAC
        IF(PFAC.NE.D1) THEN
          CALL DSCAL(NNBBASX,PFAC,OP1INT,1)
          IF(IPRINT.GE.8) THEN
             WRITE(LUPRI,'(A)') '* Scaled integrals'
             CALL OUTPAK(OP1INT,NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF        
C
C       Phase transformation/PDOINT factor if necessary
C
Chj     PDOINT is now '+', '-', or '0',
Chj     corresponding to a factor of +1, -1, or 0 on that block.
        NM1 = 0
        DO J = 1,4
           JCHR = PDOINT(KLBL)(J:J)
           IF (JCHR .EQ. '0') THEN
              IDOINT(J) = 0
           ELSE IF (JCHR .EQ. '+') THEN
              IDOINT(J) = 1
           ELSE IF (JCHR .EQ. '-') THEN
              NM1 = 1
              IDOINT(J) = -1
           ELSE
              CALL QUIT('Illegal factor in PDOINT')
           END IF
        END DO
C
        IF (NZ.LT.4 .OR. NM1.GT.0) THEN
          CALL QNPHASE(IZ,IBRP,OP1INT,IDOINT)
          IF(IPRINT.GE.8) THEN
             WRITE(LUPRI,'(A/3A)')
     &          '* Phase and/or PDOINT transformed integrals:',
     &          '  ( PDOINT = ',PDOINT(KLBL),' )'
             CALL OUTPAK(OP1INT,NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF
C
C       Determine where integrals go and add them on
C
Chj     IH = IHQMAT(IQ,1)
Chj       IF(IH.EQ.1) THEN
Chj         CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
Chj       ELSEIF(IH.EQ.2) THEN
Chj         CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
Chj       ENDIF
Chj
        IF(FIRST(IQP)) THEN
Chj       CALL DSPTSI(NTBAS(0),OP1INT,PMAT(1,IQP))
          CALL DCOPY(NNBBASX,OP1INT,1,PMAT(1,IQP),1)
          FIRST(IQP) = .FALSE.
        ELSE
Chj       CALL DSPPSI(NTBAS(0),OP1INT,PMAT(1,IQP))
          CALL DAXPY(NNBBASX,D1,OP1INT,1,PMAT(1,IQP),1)
        ENDIF
 10   CONTINUE
      IF(IPRINT.GE.8) THEN
        DO IZ = 1,NZ
          WRITE(LUPRI,'(A,I2,A)') '* Property matrix (IZ=',IZ,')'
          IF (FIRST(IZ)) THEN
             WRITE(LUPRI,'(/A)') '  --- zero matrix ---'
          ELSE
             CALL OUTPAK(PMAT(1,IZ),NTBAS(0),1,LUPRI)
Chj          CALL OUTPUT(PMAT(1,IZ),1,NTBAS(0),1,NTBAS(0),
Chj  &                   NTBAS(0),NTBAS(0),1,LUPRI)
          END IF
        ENDDO
      ENDIF
CLuuk:  Here is the right place to generalize the spinfree formalism to
Cluuk:  properties. Add a call to a generalized version of SPFAO or
Cluuk:  TBOSON (and backtransformation) here. This will make the 
Cluuk:  property matrices spinfree once and for all ! (but)
CLuuk:  In addition you may need to explicity zero out the 2e-so
CLuuk:  contribution when doing second or higher order calculations
      RETURN
C
 1000 CONTINUE
C
C     Not able to read IDOINT information from PDOINT
C
      WRITE(LUPRI,'(A,A)') 'PRPMAO: Class label:',PRPLBL(KLBL)
      WRITE(LUPRI,'(A,A)') 'PRPMAO ERROR: Not able to read COMB =',
     &                     PDOINT(KLBL)
      CALL QUIT('PRPMAO: Not able to read COMB')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdin */
      SUBROUTINE LOWDIN(VMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construction of the transformation matrix for Lowdin
C        ortogonalization
C
C     METHOD:
C        In this subroutine V+ is computed.  As the basis set is
C        real it reduces to the transpose:
C
C           V+ = TRN(V) = t(s**(-0.5))
C
C     Written by T.Saue
C     Last revision: Nov 23 1994 - tsaue
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbbas.h"
      DIMENSION VMAT(*),WORK(LWORK)
C
      CALL QENTER('LOWDIN')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,MXBBAS ,WORK,KFREE,LFREE)
      IF(ISPHTR.NE.0) THEN
        NCSPH = 0
        DO LL = 1,NHTYP
          NCSPH = NCSPH + LL*(LL+1)*(2*LL-1)/2
        ENDDO
        NICSP = NHTYP
      ELSE
        NCSPH = 0
        NICSP = 0
      ENDIF
      CALL MEMGET('INTE',KICSP,NICSP ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCSPH,NCSPH ,WORK,KFREE,LFREE)
C
C     Make Lowdin matrix
C     ==================
C
      CALL LOWDI1(VMAT,WORK(KSMAT),WORK(KEIG),WORK(KCSPH),
     &            WORK(KICSP),WORK(KFREE),LFREE,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('LOWDIN.LOWDI1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
      CALL QEXIT('LOWDIN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdi1 */
      SUBROUTINE LOWDI1(VMAT,SMAT,EIG,CSPH,ICSP,WORK,LWORK,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construction of the transformation matrix for Lowdin
C        ortogonalization
C
C     METHOD:
C        In this subroutine V+ is computed.  As the basis set is
C        real it reduces to the transpose:
C
C           V+ = TRN(V) = t(s**(-0.5))
C
C        Written by T.Saue
C        Last revision: Nov 23 1994 - tsaue
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
      CHARACTER COMP(2)*1
      DIMENSION VMAT(*),SMAT(*),EIG(*),CSPH(*),ICSP(*),WORK(LWORK)
      LOGICAL LBIT
C
      CALL QENTER('LOWDI1')
      IF(IPRINT.GE.5) CALL TITLER('Output from LOWDI1','*',103)
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     Get overlap matrix packed on fermion irreps
C     ===========================================
C
      CALL GTOVLT(SMAT,D1,WORK,LWORK)
C
C     If spherical transformations, get general transformation matrix
C     ===============================================================
C
      IF(ISPHTR.NE.0) CALL GETSPH(CSPH,ICSP,NHTYP,WORK,LWORK,IPRINT)
C
C     Initialize Lowdin matrix
C     ========================
C
      CALL DZERO(VMAT,N2BBASX)
      NBRP = 4/NZ
      ISOFF = 1
      IVOFF = 1
      CALL HEADER('Generating Lowdin matrix:',-1)
      DO IFRP = 1,NFSYM
        JSOFF = ISOFF
        JVOFF = IVOFF
C
C       Large components
C       ----------------
C
        IC = 1
        NFORB(IFRP,IC) = 0
        IP    = MOD(IFRP+IC,NFSYM) + 1
        DO JSYM = 1,NBRP
          NBORB(JSYM,IFRP,IC) = 0
          IREP = JFSYM(JSYM,IP) - 1
          IF(NBBAS(IREP,IC).GT.0) THEN
            IF(LBIT(ISPHTR,IC)) THEN
              CALL SPHLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                    NFBAS(IFRP,0),NEFF,STOL(IC),EIG,
     &                    CSPH,ICSP,NHTYP,WORK,LWORK,IPRINT)
              NDEL  = NBBAS(IREP,IC) - NEFF
            ELSE
              CALL LOWGEN(SMAT(JSOFF),NBBAS(IREP,IC),
     &                    VMAT(JVOFF),NEFF,NFBAS(IFRP,0),
     &                    STOL(IC),IPRINT,EIG,WORK,LWORK)
              NDEL  = NBBAS(IREP,IC) - NEFF
              WRITE(LUPRI,'(3X,A1,3X,A3,3X,A,I10)') 
     &          COMP(IC),REP(IREP),'* Deleted: ',NDEL
            ENDIF
C
C           Zero columns of Lowdin matrix that are to be deleted
C
            IF(NDEL.GT.0) THEN
              J = JVOFF + NFBAS(IFRP,0)*NEFF 
              DO I = 1,NDEL
                CALL DZERO(VMAT(J),NBBAS(IREP,IC))
                J = J + NFBAS(IFRP,0)
              ENDDO
            ENDIF
            JSOFF = JSOFF + (NFBAS(IFRP,0)+1)*NBBAS(IREP,IC)
            JVOFF = JVOFF + NFBAS(IFRP,0)*NEFF + NBBAS(IREP,IC)
            NBORB(JSYM,IFRP,IC) = NEFF
            NFORB(IFRP,IC) = NFORB(IFRP,IC) + NEFF
          ENDIF
        ENDDO
C
C       Small components
C       ----------------
C
        IC = 2
        NFORB(IFRP,IC) = 0
        IP    = MOD(IFRP+IC,NFSYM) + 1
        DO JSYM = 1,NBRP
          NBORB(JSYM,IFRP,IC) = 0
          IREP = JFSYM(JSYM,IP) - 1
          IF(NBBAS(IREP,IC).GT.0) THEN
            IF(LBIT(ISPHTR,IC)) THEN
              CALL RKBLOW(IREP,IC,SMAT(JSOFF),VMAT(JVOFF),
     &                    NFBAS(IFRP,0),NEFF,STOL(IC),EIG,
     &                    CSPH,ICSP,NHTYP,WORK,LWORK,IPRINT)
              NDEL  = NBBAS(IREP,IC) - NEFF
            ELSE
              CALL LOWGEN(SMAT(JSOFF),NBBAS(IREP,IC),
     &                    VMAT(JVOFF),NEFF,NFBAS(IFRP,0),
     &                    STOL(IC),IPRINT,EIG,WORK,LWORK)
              NDEL  = NBBAS(IREP,IC) - NEFF
              WRITE(LUPRI,'(3X,A1,3X,A3,3X,A,I10)') 
     &           COMP(IC),REP(IREP),'* Deleted: ',NDEL
            ENDIF
C
C           Zero columns of Lowdin matrix that are to be deleted
C
            IF(NDEL.GT.0) THEN
              J = JVOFF + NFBAS(IFRP,0)*NEFF 
              DO I = 1,NDEL
                CALL DZERO(VMAT(J),NBBAS(IREP,IC))
                J = J + NFBAS(IFRP,0)
              ENDDO
            ENDIF
            NBORB(JSYM,IFRP,IC) = NEFF
            NFORB(IFRP,IC) = NFORB(IFRP,IC) + NEFF
            JSOFF = JSOFF + (NFBAS(IFRP,0)+1)*NBBAS(IREP,IC)
            JVOFF = JVOFF + NFBAS(IFRP,0)*NEFF+NBBAS(IREP,IC)
          ENDIF
          NBORB(JSYM,IFRP,0) = NBORB(JSYM,IFRP,1) + NBORB(JSYM,IFRP,2)
        ENDDO
        NFORB(IFRP,0) = NFORB(IFRP,1) + NFORB(IFRP,2)
        IF ( IPRINT .GE. 5) THEN
           WRITE(LUPRI,'(A,I2)') '* Lowdin matrix for ifrep = ',IFRP
           CALL OUTPUT(VMAT(IVOFF),1,NFBAS(IFRP,0),1, NFORB(IFRP,0),
     &                   NFBAS(IFRP,0),NFORB(IFRP,0),1,LUPRI)
        END IF
        ISOFF = ISOFF + NFBAS(IFRP,0)*NFBAS(IFRP,0)
        IVOFF = IVOFF + NFBAS(IFRP,0)*NFORB(IFRP,0)
      ENDDO
      CALL QEXIT('LOWDI1')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rgetint */
      SUBROUTINE RGETINT(WORD,OP1INT,DOINT,NCOMP,IORDER,NPQUAD,
     &                  MXFORM,TOFILE,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for getting one-electron property integrals
C     Written by T.Saue October 1994
C     Last revision : Nov 15 1994
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "shells.h"
#include "symmet.h"
#include "nuclei.h"
      DIMENSION OP1INT(*),WORK(LWORK)
      CHARACTER WORD*7,MXFORM*6,LABINT(3*MXCOOR)*8
      LOGICAL DOINT(2,2),PROPRI,TRIANG,TOFILE
C
      CALL QENTER('RGETINT')
#include "memint.h"
C
      MAXTYP = 3*MXCOOR
      IF ((WORD.EQ.'ELFGRDC').OR.(WORD.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
      ELSE
         NINTAD = MAXTYP
      END IF
      PROPRI   = .FALSE.
      TRIANG   = .TRUE.
C
C     Memory allocation
C
      CALL MEMGET('INTE',KINTRP,MAXTYP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINTAD,NINTAD,WORK,KFREE,LFREE)
C
C     Calculate integrals
C
      CALL PR1IN1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,OP1INT,NCOMP,TOFILE,MXFORM,DOINT)
C
C     Memory deallocation
C
      CALL MEMREL('RGETINT.PR1IN1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RGETINT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vfsort */
      SUBROUTINE VFSORT(NLEFF,NSEFF,VMAT,VBUF,WORK,LWORK)
C*****************************************************************************
C
C     Sort Lowdin matrix on fermion symmetries
C
C     Written by T.Saue Sep 12 1995
C     Last revision: Sep 12 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C Used from COMMON blocks:
C  DCBBAS: NTBAS(0)
C
#include "dcbbas.h"
      DIMENSION VMAT(*),VBUF(*),WORK(LWORK)
      CALL QENTER('VFSORT')
#include "memint.h"
C
C     Memory allocation
      CALL MEMGET('INTE',KSORT,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL VFSOR1(NLEFF,NSEFF,VMAT,VBUF,WORK(KSORT))
C     Memory release
      CALL MEMREL('VFSORT.VFSOR1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('VFSORT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck vfsor1 */
      SUBROUTINE VFSOR1(NLEFF,NSEFF,VMAT,VBUF,ISORT)
C*****************************************************************************
C
C     Sort Lowdin matrix on fermion symmetries
C
C     Written by T.Saue Sep 12 1995
C     Last revision: Sep 12 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcbbas.h"
#include "dcborb.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION VBUF(NBAST,*)
#else
      DIMENSION VBUF(NTBAS(0),*)
#endif
      DIMENSION VMAT(*),ISORT(*),NVEC(2)
C
C
C     Fermion symmetry 1
C     ==================
C
      CALL IZERO(NFORB,6)
      IOFF    = 1
      N2      = 0
      NVEC(1) = NLEFF
      NVEC(2) = NSEFF
      JJ = 0
      DO IC = 1,2
        DO J = 1,NVEC(IC)
          JJ = JJ + 1
          IND = IDAMAX(NTBAS(0),VBUF(1,JJ),1)
          IF(IND.LE.NFBAS(1,0)) THEN
            NFORB(1,IC) = NFORB(1,IC) + 1
            CALL DCOPY(NFBAS(1,0),VBUF(1,JJ),1,VMAT(IOFF),1)
            IOFF = IOFF + NFBAS(1,0)
          ELSE
            N2          = N2 + 1
            NFORB(2,IC) = NFORB(2,IC) + 1
            ISORT(N2)   = JJ
          ENDIF
        ENDDO
      ENDDO
C
C     Fermion  irrep 2
C     ================
C
      DO J = 1,N2
        CALL DCOPY(NFBAS(2,0),VBUF(NFBAS(1,0)+1,ISORT(J)),1,
     &             VMAT(IOFF),1)
        IOFF = IOFF + NFBAS(2,0)
      ENDDO
      NFORB(1,0) = NFORB(1,1)+NFORB(1,2)
      NFORB(2,0) = NFORB(2,1)+NFORB(2,2)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck freemt */
      SUBROUTINE FREEMT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for the construction of the one-electron free particle matrix
C     in SO-basis
C     Trond Saue, October 1995 - Oslo
C     Last revision: May 21 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION FOCK(*),WORK(LWORK)
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      CALL QENTER('FREEMT')
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     ****** Construct Fock matrix ******
C
      CALL FREEM1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &              WORK(KFREE),LFREE,IPRINT)
C     Memory deallocation
      CALL MEMREL('FREEMT.FREEM1',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('FREEMT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck freem1 */
      SUBROUTINE FREEM1(FOCK,FTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Driver for the construction of the one-electron free particle matrix
C     in SO-basis
C     Trond Saue, October 1995 - Oslo
C
C-----------------------------------------------------------------------------
C
C     The free particle matrix has the following form in the C1 case:
C
C             La      Sa      Lb      Sb
C     La  |   0    | -icDz |   0   | -icD- |
C     Sa  | -icDz  |  -C2  | -icD- |   0   |
C     Lb  |   0    | -icD+ |   0   | -icDz |
C     Sb  | -icD+  |   0   | -icDz |  -C2  |
C
C
C     The integrals:
C        S  - overlap integrals                  symmetric
C        C2 - (2c^2)S                            symmetric
C        D  - dipole velocity integrals          antisymmetric
C            D+ = Dx + iDy
C            D- = Dx - iDy
C     The ordering of components is chosen to reveal the following
C     matrix structure which is due to time reversal symmetry:
C
C     |  A  |  B  |
C     | -B* |  A* |
C
C     In this routine the full matrices A and B are constructed:
C
C            A(real)   |   A(imag)   |    B(real)   |  B(imag)
C     ------------------------------------------------------------
C           La      Sa     La    Sa  |  Lb     Sb     Lb     Sb
C     La |   0  |   0  |   0  | -cDz |    0  | -cDy |   0  | -cDx |
C     Sa |   0  |  -C2 | -cDz |   0  |  -cDy |   0  | -cDx |   0  |
C
C     It should be noted that Ar is symmetric, whereas Ai,Br and Bi
C     are antisymmetric.
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbgen.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),OP1INT(NNBBASX),
     &          WORK(LWORK)
C
C     Initialize
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL LSET(NZ,.TRUE.,FIRST)
C
C     Beta matrix
C     ===========
C
      IOP = IPBETAMT
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
C
C     Kinetic energy
C     ==============
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
      IF(IPRHAM.GE.5) THEN
        CALL HEADER('Unsorted one-electron Fock matrix',-1)
        CALL PRQMAT(FOCK,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,NZ,WORK,LWORK)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER('One-electron free particle matrix',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(FOCK(I2BASX(I,I)+1,1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gmotra */
      SUBROUTINE GMOTRA(WORK,LWORK)
C*****************************************************************************
C
C     This matrix provides initial transformation matrix from
C     SO- to MO-basis. With no restrictions on kinetic balance
C     the Lowdin matrix is used. With restricted balance, the
C     Lowdin matrix is used to diagonalize the free-particle
C     matrix. Unphysical solutions are then deleted from the set
C     of vectors and the remaining vectors used as transformation matrix.
C
C     Written by T.Saue, October 1995
C     Last revision May 21 1996 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C   DCBBAS: N2BBASX
C   PAMIOU: LUTMAT
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr1.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION WORK(LWORK)
      CALL QENTER('GMOTRA')
#include "memint.h"
C
C     By default we do not partition the space. This may be
C     modified if we have linear symmetry or approximate the
C     Hamiltonian : Levy-Leblond, Spinfree, ..
C
      SUB_BL = .FALSE.
C
C     Make one-center projection matrix
C     =================================
C
      IF (ONECAP.AND.(INTV1C.EQ.1.AND.ICTLV1C(2).EQ.1)) 
     &    CALL MK1CEX (WORK,KFREE,LFREE,IPRGEN)
C
C     Get Lowdin matrix
C     =================
C
      CALL MEMGET('REAL',KVMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL LOWDIN(WORK(KVMAT),WORK(KFREE),LFREE,IPRONE)
C
C     Unrestricted kinetic balance; use Lowdin matrix
C     ===============================================
C
C     LV : The unrestricted procedure should be revised so that
C     the large components are also sorted on boson irreps !!!
C
      CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)
      IF(.NOT.URKBAL.AND.NTBAS(2).EQ.0) THEN
        URKBAL = .TRUE.
        WRITE(LUPRI,'(A)') 
     &   '*** WARNING *** : No small bases. Setting URKBAL = T !'
      ENDIF
      IF(URKBAL) THEN
        NZT = 1
        DO IFRP = 1,NFSYM
          NESH(IFRP)  = NFORB(IFRP,1)
          NPSH(IFRP)  = NFORB(IFRP,2)
          NORB(IFRP)  = NESH(IFRP) + NPSH(IFRP)
          NTMO(IFRP)  = NORB(IFRP)
        ENDDO
        KTMAT = KVMAT
      ELSE
C
C     Charge conjugation symmetry
C     ===========================
C
        CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL DELPOS(WORK(KTMAT),WORK(KVMAT),WORK(KTBUF),WORK(KEIG),
     &            WORK(KFREE),LFREE)
      ENDIF
      N2TMT   = 0
      N2TMOTQ = 0
      DO IFRP = 1,NFSYM
        NISHMF(IFRP) = NISH(IFRP)
        I2TMT(IFRP) = N2TMT
        N2TMT = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
        I2TMOT(IFRP) = N2TMOTQ
        N2TMOTQ = N2TMOTQ + NTMO(IFRP)*NTMO(IFRP)*NZ
      ENDDO
C
C     Construct right index transformed overlap matrix Saomo = Sao C
C
      CALL MKSAOMO(WORK(KTBUF),WORK(KTMAT),IPRONE,WORK(KFREE),LFREE)
C
C     If requested, reduce transformation matrix by projection with
C     respect to a set of reference orbitals
C
      IF(PROJEC) THEN
        CALL PSIPRJ(WORK(KTMAT),WORK(KTBUF),IPRONE,WORK,KFREE,LFREE)
        N2TMOTQ = 0
        DO IFRP = 1,NFSYM
          I2TMOT(IFRP) = N2TMOTQ
          N2TMOTQ = N2TMOTQ + NTMO(IFRP)*NTMO(IFRP)*NZ
          NORB(IFRP) = NTMO(IFRP) + NFRO(IFRP)
        ENDDO
        CALL SETDC2(0)
      ENDIF
C
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
      IF (.NOT.SUB_BL) THEN
         DO IFRP = 1, NFSYM
            N_SUB_BL(IFRP) = 4 / NZ
            DO IB = 1, N_SUB_BL(IFRP) 
              NORB_SUB(IB,IFRP,0) = NBORB(IB,IFRP,0)
              NORB_SUB(IB,IFRP,1) = NBORB(IB,IFRP,1)
              NORB_SUB(IB,IFRP,2) = NBORB(IB,IFRP,2)
              IREP = JFSYM(IB,IFRP) - 1
              ID_SUB_BL(IB,IFRP) = IREP
            ENDDO
         ENDDO
      ENDIF
C
      IF (LEVYLE.OR.SPINFR) THEN
         SUB_BL = .TRUE.
      ENDIF
C
C     Output
C
      CALL TITLER('Orbital dimensions','*',125)
      IF ( NFSYM .EQ. 1 ) THEN
         WRITE(LUPRI,'(A,I5)')
     &        'No. of electronic orbitals (NESH): ',NESH(1)
         WRITE(LUPRI,'(A,I5)')
     &        'No. of positronic orbitals (NPSH): ',NPSH(1)
         WRITE(LUPRI,'(A,I5)')
     &        'Total no. of orbitals      (NORB): ',NORB(1)
      ELSE
         WRITE(LUPRI,'(35X,5(A5,1X,I1,1X))') 
     &        ('Irrep',I,I=1,NFSYM),'Sum'
         J = 0
         DO I = 1, NFSYM
            J = J + NESH(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X))')
     &        'No. of electronic orbitals (NESH): ',
     &        (NESH(I),I=1,NFSYM),J
         J = 0
         DO I = 1, NFSYM
            J = J + NPSH(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X))')
     &        'No. of positronic orbitals (NPSH): ',
     &        (NPSH(I),I=1,NFSYM),J
         J = 0
         DO I = 1, NFSYM
            J = J + NORB(I)
         END DO
         WRITE(LUPRI,'(A,5(I5,3X)/)')
     &        'Total no. of orbitals      (NORB): ',
     &        (NORB(I),I=1,NFSYM),J
      END IF
C
C
C     Write transformation matrices
C
      CALL OPNFIL(LUTMAT,'DFTMAT','UNKNOWN','GMOTRA')
      CALL WRITT(LUTMAT,N2TMT,WORK(KTMAT))
      CALL WRITT(LUTMAT,N2TMT,WORK(KTBUF))
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     Memory deallocation
C
      CALL MEMREL('GMOTRA',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('GMOTRA')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck delpos */
      SUBROUTINE DELPOS(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Diagonalize free particle matrix and remove unphysical solutions
C     Or transform to basis for Levy-Leblond Hamiltonian
C
C     Written by T.Saue, October 1995
C     Modified by L. Visscher, Jan 20, 1998 to include Levy-Leblond 
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM2 = -2.0D0,D0 = 0.0D0)
C
C Used from COMMON blocks
C  CBIPAM:
C  DHFCAL:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
#include "pgroup.h"
#include "cbihr1.h"
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
#include "memint.h" 
C
      FAC = DM2*CVAL*CVAL
      NSYM = 4/NZ
C
C     Generate free particle (or bare nucleus) matrix in MO-basis
C     ===========================================================
C
      IF(VEXTPJ) THEN
        CALL ONEFCK(TBUF,IPRONE,WORK,LWORK)
      ELSE
        CALL FREEMT(TBUF,IPRONE,WORK,LWORK)
      ENDIF
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NFORB(I,0).EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                         NFORB(I,0),NFORB(I,0),
     &            TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRGEN)
        IF (IPRHAM.GE.6) THEN
           IF(VEXTPJ) THEN
             CALL HEADER('Bare nucleus matrix in orthon. basis',-1)
           ELSE
             CALL HEADER('Free particle matrix in orthon. basis',-1)
           ENDIF
           WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
           CALL PRQMAT(TMAT(IOFT),NFORB(I,0),NFORB(I,0),NFORB(I,0),
     &                 NFORB(I,0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   10 CONTINUE
C
C     The Levy-Leblond equation:
C     ==========================
C
C       Set up basis for Levy-Leblond equation as:
C
C              [                (Psi(L), 0 ]
C        Psi = [                           ]
C              [(2mc)^-1 sigma.p Psi(L), 0 ]
C
      IF(LEVYLE) THEN
        NZT   = NZ
        FAC2  = -D1/FAC
        IBUF  = 1
        IMAT  = 1
        DO 20 I = 1,NFSYM
          NESH(I)  = NFORB(I,1)
          NPSH(I)  = 0
          NORB(I)  = NESH(I) + NPSH(I)
          NTMO(I)  = NORB(I)
          IF(NORB(I).EQ.0) GOTO 20
          CALL DZERO(TBUF(IBUF),NFORB(I,0)*NORB(I)*NZT)
          JBUF = IBUF
          DO J = 1, NFORB(I,1)
            TBUF(JBUF) = D1
            JBUF = JBUF + NFORB(I,0) + 1
          ENDDO
          DO IZ = 1, NZ
            JBUF = IBUF + NFORB(I,1)
            JMAT = IMAT + NFORB(I,1)
            DO J = 1, NFORB(I,1)
              CALL DCOPY (NFORB(I,2),TMAT(JMAT),1,TBUF(JBUF),1)
              CALL DSCAL (NFORB(I,2),FAC2,TBUF(JBUF),1)
              JMAT = JMAT + NFORB(I,0)
              JBUF = JBUF + NFORB(I,0)
            ENDDO
            IF (IPRHAM.GE.6) THEN
             CALL HEADER('Transformation matrix in orthon. basis',-1)
             CALL PRQMAT(TBUF(IBUF),NFORB(I,0),NORB(I),
     &                   NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
            IMAT = IMAT + NFORB(I,0) * NFORB(I,0)
            IBUF = IBUF + NFORB(I,0) * NORB(I)
          ENDDO
C         Finish by adjusting the number of positron solutions to zero
          DO ISYM = 1, NSYM
            NBORB(ISYM,I,2) = 0
            NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
          ENDDO
   20   CONTINUE
C
C     Generate transformation matrix
C     ==============================
C
        IBUF = 1
        IMAT = 1
        IOFV = 1
        DO 30 I = 1,NFSYM
          IF(NORB(I).EQ.0) GOTO 30
          CALL BCKTR1(TMAT(IMAT),NFBAS(I,0),NFORB(I,1),
     &                TBUF(IBUF),NFORB(I,0),NFORB(I,1),
     &                NFORB(I,0),NZ,
     &                NFORB(I,1),1,NFBAS(I,0),
     &                VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &                IPRGEN)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
          IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
          IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   30   CONTINUE
      ELSE      
C
C       The modified Dirac eqation:
C       ===========================
C
C       Basis as in Levy-Leblond but without fixed ratio
C
C
C               [                (Psi(L),                      0  ]
C         Psi = [                                                 ]
C               [                      0, (2mc)^-1 sigma.p Psi(L) ]
        FAC2 = D1/FAC/DM2
        NZT  = NZ
        IMAT = 1
        IBUF = 1
        DO 40 I = 1,NFSYM
          NESH(I)  = NFORB(I,1)
          NPSH(I)  = NFORB(I,1)
          NORB(I)  = NESH(I) + NPSH(I)
          NTMO(I)  = NORB(I)
          IF(NORB(I).EQ.0) GOTO 40
C
C         Prepare small component part
C         ----------------------------
C
          LWRK = MXBBAS*MXBBAS*NSYM
          CALL MEMGET('REAL',KSMT,LWRK,WORK,KFREE,LFREE)
          LWRK = NFORB(I,1) * NFORB(I,1)
          CALL MEMGET('REAL',KVMT,LWRK,WORK,KFREE,LFREE)
          IVMT  = KVMT
          IBOL  = 0
          NEFF  = 0
          DO ISYM = 1,NSYM
          IF(NBORB(ISYM,I,1).GT.0) THEN
C
C           Make (sigma.p)(sigma.p) = p2
C
            IREP = JFSYM(ISYM,I) - 1
            JSL = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
            JLS = IMAT + NFORB(I,1) * NFORB(I,0) + IBOL
            CALL QGEMM(NBORB(ISYM,I,1),NBORB(ISYM,I,1),NFORB(I,2),D1,
     &          'N','N',IPQTOQ(1,0),TMAT(JLS),NFORB(I,0),NFORB(I,0),NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &          D0,IPQTOQ(1,0),WORK(KSMT),
     &          NBORB(ISYM,I,1),NBORB(ISYM,I,1),1)
            IF (IPRHAM.GE.7) THEN
              WRITE(LUPRI,'(A,A3)') ' * Boson irrep: ',REP(IREP)
              CALL HEADER('Kinetic energy matrix in orthon. basis',-1)
              CALL PRQMAT(WORK(KSMT),NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                    NBORB(ISYM,I,1),NBORB(ISYM,I,1),
     &                    1,IPQTOQ(1,0),LUPRI)
            ENDIF
            CALL LOWGEN(WORK(KSMT),NBORB(ISYM,I,1),WORK(IVMT),NEFFI,
     &                  NBORB(ISYM,I,1),STOL(2),IPRONE,EIG,
     &                  WORK(KFREE),LFREE)
            IF (NEFFI.LT.NBORB(ISYM,I,1)) WRITE(LUPRI,'(A)') 
     &            '*** WARNING *** Linear dependence in RKB basis'
            NEFF = NEFF + NEFFI
            NBORB(ISYM,I,2) = NEFFI
             IF (IPRHAM.GE.7) THEN
              CALL HEADER('Ext.transf. matrix in orthon. basis',-1)
              CALL PRQMAT(WORK(IVMT),NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                    NBORB(ISYM,I,1),NBORB(ISYM,I,2),
     &                    1,IPQTOQ(1,0),LUPRI)
             ENDIF
            IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
            IBOL = IBOL + NBORB(ISYM,I,1)
          ELSE
            NBORB(ISYM,I,2) = 0
          ENDIF
          NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
          ENDDO
          NPSH(I)  = NEFF
          NORB(I)  = NESH(I) + NPSH(I)
          NTMO(I)  = NORB(I)
C
          CALL DZERO(TBUF(IBUF),NFORB(I,0)*NORB(I)*NZ)
          IBO  = 0
          IBOL = 0
          IVMT = KVMT
          DO ISYM = 1,NSYM
          IF(NBORB(ISYM,I,1).GT.0) THEN
            IREP = JFSYM(ISYM,I) - 1
C
C           Small component part
C
            JBUF = IBUF + NFORB(I,0)*IBO + NFORB(I,1)
            JSL  = IMAT + IBOL * NFORB(I,0) + NFORB(I,1)
            CALL QGEMM(NFORB(I,2),NBORB(ISYM,I,2),NBORB(ISYM,I,1),D1,
     &           'N','N',IPQTOQ(1,0),TMAT(JSL),NFORB(I,0),NFORB(I,0),NZ,
     &           'N','N',IPQTOQ(1,0),WORK(IVMT),
     &                   NBORB(ISYM,I,1),NBORB(ISYM,I,2),1,
     &           D0,IPQTOQ(1,0),TBUF(JBUF),NFORB(I,0),NORB(I),NZ)
            IBO = IBO + NBORB(ISYM,I,2)
C
C           Large component part
C
            JBUF = IBUF + NFORB(I,0)*IBO + IBOL
            DO J = 1, NBORB(ISYM,I,1)
              TBUF(JBUF) = D1
              JBUF = JBUF + NFORB(I,0) + 1
            ENDDO
            IBO = IBO + NBORB(ISYM,I,1)
C
C           Update offsets
C
            IBOL = IBOL + NBORB(ISYM,I,1)
            IVMT = IVMT + NBORB(ISYM,I,1)*NBORB(ISYM,I,2)
          ENDIF
          ENDDO
          IF (IPRHAM.GE.6) THEN
           CALL HEADER('Transformation matrix in orthon. basis',-1)
           CALL PRQMAT(TBUF(IBUF),NFORB(I,0),NORB(I),
     &                 NFORB(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
          IMAT = IMAT + NFORB(I,0)*NFORB(I,0)*NZ
          IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
   40   CONTINUE
        CALL MEMREL('DELPOS',WORK,KSMT,KSMT,KFREE,LFREE)
C
        IF(VEXTPJ.OR.FREEPJ) THEN
C
C       Deleting positronic solutions of:
C       =================================
C         FREEPJ - free particle matrix
C         VEXTPJ - bare nucleus matrix          
C
          LWRK  = 0
Ctsaue:   These variables are later set by SETDC2, but needed here
C         in modified form....
          N2TMT  = 0
          NORBT  = 0
          N2ORBT = 0
          DO I = 1,NFSYM
            I2TMT(I) = N2TMT
            N2TMT    = N2TMT + NFORB(I,0)*NORB(I)*NZT
            IORB(I)  = NORBT
            NORBT    = NORBT  + NORB(I)
            I2ORBT(I)= N2ORBT*NZ
            N2ORBT   = N2ORBT + NFBAS(I,0)*NORB(I)
          ENDDO
          CALL MEMGET('REAL',KBUF,N2ORBT*NZ,WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KBEIG,NTBAS(0),WORK,KFREE,LFREE)
C
C         Transform free particle/bare nucleus matrix
C         to modified Dirac form; 
          IMAT = 1
          IBUF = 1
          IOFT = 1
          IOFV = 1
          DO 50 I = 1,NFSYM
            IF(NORB(I).EQ.0) GOTO 50
            IBUF = KBUF + I2ORBT(I)
            CALL QTRANS('AOMO','S',D0,NFORB(I,0),NFORB(I,0),
     &           NORB(I),NORB(I),
     &           TMAT(IOFT),NFORB(I,0),NFORB(I,0),NZ,IPQTOQ(1,0),
     &           WORK(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &           TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &           TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),NZT,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)
          IF (IPRHAM.GE.6) THEN
            IF(VEXTPJ) THEN
              CALL HEADER('Bare nucleus matrix in mod.MO-basis',-1)
            ELSE
              CALL HEADER('Free particle matrix in mod.MO-basis',-1)
            ENDIF
            CALL PRQMAT(WORK(IBUF),NORB(I),NORB(I),
     &                 NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
C         get full transformation as well !
            CALL BCKTR1(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                  TBUF(1+I2TMT(I)),NFORB(I,0),NORB(I),
     &                  NFORB(I,0),NZ,
     &                  NORB(I),1,NFBAS(I,0),
     &                  VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &                  IPRHAM)
            IF (IPRHAM.GE.6) THEN
              CALL HEADER('Transformation matrix',-1)
              CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
            ENDIF
            IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
            IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
            IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
            IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
 50       CONTINUE
C
C         Diagonalize
C
          CALL DFDIAG (WORK(KBUF),EIG,WORK(KBEIG),TBUF,
     &                 WORK,KFREE,LFREE)
C
C         Backtransform only electronic solutions
C
          N2TMT = 0
          IBUF = KBUF
          IMAT = 1
          NBUF = 0
          DO 60 I = 1,NFSYM
          IF(NESH(I).EQ.0) GOTO 60
          I2TMT(I) = N2TMT
          N2TMT = NFBAS(I,0)*NORB(I)*NZ
          IF(IPRHAM.GE.5) THEN
            CALL HEADER('Eigenvalues',-1)
            WRITE(LUPRI,'(I5,E16.8)') (J,EIG(IORB(I)+J),J=1,NORB(I))
          ENDIF
          CALL BCKTR1(WORK(IBUF),NFBAS(I,0),NESH(I),
     &              TBUF(I2ORBT(I)+1),NORB(I),NORB(I),
     &              NORB(I),NZ,NESH(I),(NPSH(I)+1),NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,
     &              IPRHAM)
          NPSH(I) = 0
          NORB(I) = NESH(I)
          NTMO(I) = NORB(I)
          NBUF = NFBAS(I,0)*NESH(I)*NZ
          CALL DCOPY(NBUF,WORK(IBUF),1,TMAT(IMAT),1)
          IBUF = IBUF + NBUF
          IMAT = IMAT + NBUF
C         Finish by adjusting the number of positron solutions to zero
          DO ISYM = 1, NSYM
            NBORB(ISYM,I,2) = 0
            NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
          ENDDO
 60       CONTINUE
          CALL MEMREL('DELPOS',WORK,KWORK,KWORK,KFREE,LFREE)
        ELSE      
C
C       Generate transformation matrix
C       ==============================
C
          IMAT = 1
          IBUF = 1
          IOFV = 1
          DO 70 I = 1,NFSYM
            IF(NFORB(I,0).EQ.0) GOTO 70
            CALL BCKTR1(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                  TBUF(IBUF),NFORB(I,0),NORB(I),
     &                  NFORB(I,0),NZ,
     &                  NORB(I),1,NFBAS(I,0),
     &                  VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &                  IPRGEN)
            IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
            IBUF = IBUF + NFORB(I,0)*NORB(I)*NZ
            IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   70     CONTINUE
        ENDIF
      ENDIF      
C
C     Print section
C
      CALL HEADER('Output from DELPOS',-1)
      IF(LEVYLE.OR.VEXTPJ.OR.FREEPJ) THEN
        WRITE(LUPRI,'(A)') '* All positronic solutions deleted !'
      ELSE
        WRITE(LUPRI,'(A)') '* Applied strict kinetic balance !'
      ENDIF
      IF(IPRHAM.GE.5) THEN
        CALL HEADER('Final transformation matrix ',-1)
        IMAT = 1
        DO 80 I = 1,NFSYM
        IF(NORB(I).GT.0) THEN
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        ENDIF
 80     CONTINUE
      ENDIF
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck prprea */
      SUBROUTINE PRPREA(LU1INT,PRPLBL,RTNLBL,PCOMB,PRPINT,NDIM,IPRINT)
C***********************************************************************
C
C     Read property integrals from file with label PRPLBL
C     and class blocks defined by COMB 
C     Return integrals in PRPINT and return labels in RTNLBL
C
C     Written by T.Saue, January 1995
C     Last revision: May 19 1996 - tsaue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLB3
      CHARACTER PRPLBL*8, RTNLBL(2)*8, PCOMB*4, COMB*4
      DIMENSION PRPINT(NDIM)
C
Chj   PDOINT is now '+', '-', or '0',
Chj   corresponding to a factor of +1, -1, or 0 on that block.
      DO J = 1,4
        IF (PCOMB(J:J) .EQ. '0' .OR. PCOMB(J:J) .EQ. 'F') THEN
           COMB(J:J) = 'F'
        ELSE
           COMB(J:J) = 'T'
        END IF
      END DO
C
      REWIND LU1INT
      IF(FNDLB3(PRPLBL,COMB,RTNLBL,LU1INT)) THEN
        CALL READT(LU1INT,NDIM,PRPINT)
      ELSE
        WRITE(LUPRI,'(A,I4)') '>>> Property '//PRPLBL//' '//COMB//
     &                     ' not found on LU1INT = ',LU1INT
        CALL QUIT('PRPREA:Property '//PRPLBL//' not found')
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndlb3 */
      LOGICAL FUNCTION FNDLB3(SRCLBL,COMB,RTNLBL,LU)
C***********************************************************************
C
C  18-Jan-1995 tsaue
C  1. Find MOLECULE label and block COMB on unformatted file LU
C       a)Cannot find label: return FALSE
C       b)Finds label      : continues
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLB2
      CHARACTER SRCLBL*8, RTNLBL(2)*8, COMB*4
      FNDLB3 = .FALSE.
   10 CONTINUE
      IF(FNDLB2(SRCLBL,RTNLBL,LU)) THEN
C     Found label SRCLBL on LU; check COMB
        IF(RTNLBL(2)(5:8).NE.COMB) GOTO 10
        FNDLB3 = .TRUE.
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtovlx */
      SUBROUTINE GTOVLX(SMAT,SSMETR,WORK,LWORK)
C*****************************************************************************
C
C     Get overlap matrix from file AOPROPER
C     The SS block is scaled by the factor SSMETR
C     No symmetry packing
C
C     Written by T.Saue May 26 1996
C     Last revision: Jan 20 1998 - L. Visscher
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(N2BBASX),WORK(LWORK)
      CHARACTER*8 RTNLBL(2)
      CALL QENTER('GTOVLX')
#include "memint.h"
C      
C     Allocate integral buffer
C     ========================
C
      CALL MEMGET('REAL',KSBUF,NNBBASX,WORK,KFREE,LFREE)
C
C     Read integral from file
C     =======================
C
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL PRPREA(LU1INT,'OVERLAP ',RTNLBL,'TFFT',WORK(KSBUF),NNBBASX,
     &            IPRONE)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Make full matrix
C     ================
C
      CALL DSPTSI(NTBAS(0),WORK(KSBUF),SMAT)
C
C     Memory deallocation
C     ===================
C
      CALL MEMREL('GTOVLX',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(SMAT,1,WORK,LWORK)
C
C     Scale SS-block
C     ==============
C
      IF(SSMETR.NE.D1) THEN
        DO I = 1,NFSYM
           JJ =  I2BASX(I,I) + (NFBAS(I,1) * NTBAS(0)) + NFBAS(I,1) + 1
           DO J = 1, NFBAS(I,2)
              CALL DSCAL (NFBAS(I,2),SSMETR,SMAT(JJ),1)
              JJ = JJ + NTBAS(0)
           ENDDO
        ENDDO
      ENDIF
C
      CALL QEXIT('GTOVLX')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gtovlt */
      SUBROUTINE GTOVLT(SMAT,SSMETR,WORK,LWORK)
C*****************************************************************************
C
C     Get overlap matrix from file AOPROPER
C     The SS block is scaled by the factor SSMETR
C     Packed on fermion irreps
C
C     Written by T.Saue May 26 1996
C     Last revision: Jan 20 1998 - L. Visscher
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "dcbgen.h"
#include "cbihr1.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "symmet.h"
      DIMENSION SMAT(*),WORK(LWORK)
      CHARACTER*8 RTNLBL(2)
      CALL QENTER('GTOVLT')
#include "memint.h"
C      
C     Allocate integral buffer
C     ========================
C
      CALL MEMGET('REAL',KSTRI,NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSBUF,N2BBASX,WORK,KFREE,LFREE)
C
C     Read integral from file
C     =======================
C
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL PRPREA(LU1INT,'OVERLAP ',RTNLBL,'TFFT',WORK(KSTRI),NNBBASX,
     &            IPRONE)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Make full matrix
C     ================
C
      CALL DSPTSI(NTBAS(0),WORK(KSTRI),WORK(KSBUF))
C
      CALL DZERO(SMAT,N2BAST)
      NBRP     = 4/NZ
      IOFF = 0
      DO I = 1,NFSYM
        DO IC = 1,2
          IP    = MOD(I+IC,NFSYM) + 1
          DO JSYM = 1,NBRP
            ISYM           = JFSYM(JSYM,IP)
            IREP           = ISYM - 1
            IF(NBBAS(IREP,IC).GT.0) THEN
              ISOFF = IOFF + (IBBAS(IREP,IC)-IBAS(I))*(NFBAS(I,0)+1)
              IUOFF = ICOS(ISYM,IC)*(NTBAS(0)+1)
              DO J = 1,NBBAS(IREP,IC)
                CALL DCOPY(NBBAS(IREP,IC),WORK(KSBUF+IUOFF),1,
     &                                    SMAT(1+ISOFF)    ,1)
                IUOFF = IUOFF + NTBAS(0)
                ISOFF = ISOFF + NFBAS(I,0)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        IOFF = IOFF + N2BAS(I)
      ENDDO
C
C     Scale SS-block
C     ==============
C
      IF(SSMETR.NE.D1) THEN
        IOFF = 1
        DO I = 1,NFSYM
          JJ = IOFF + NFBAS(I,1)*(NFBAS(I,0)+1)
          DO J = 1, NFBAS(I,2)
            CALL DSCAL (NFBAS(I,2),SSMETR,SMAT(JJ),1)
            JJ = JJ + NFBAS(I,0)
          ENDDO
          IOFF = IOFF + N2BAS(I)
        ENDDO
      ENDIF
C
C
C     Memory deallocation
C     ===================
C
      CALL MEMREL('GTOVLT',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('GTOVLT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck heinz */
      SUBROUTINE HEINZ(WORK,LWORK)
C*****************************************************************************
C
C     Solve one-electron systems with finite nuclei and project
C     out all free positronic solutions
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION WORK(LWORK)
      CALL QENTER('HEINZ')
#include "memint.h"
      CALL MEMGET('REAL',KCMO ,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFAO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFMO,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVBUF,N2BBASX ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2TMT   ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
      CALL HEINZ2(WORK(KCMO),WORK(KEIG),WORK(KFAO),WORK(KFMO),
     &            WORK(KFBUF),WORK(KVBUF),
     &            WORK(KTMAT),WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK(KFREE),LFREE,IPRGEN)
C
      CALL MEMREL('HEINZ',WORK,KWORK,KWORK,KFREE,LFREE)      
      CALL QEXIT('HEINZ')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck heinz2 */
      SUBROUTINE HEINZ2(CMO,EIG,FAO,FMO,FBUF,VBUF,TMAT,
     &                  FTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Solve one-electron systems with finite nuclei and project
C     out all free positronic solutions
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "codata.h"
      PARAMETER(DP5 = 0.5D0,D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      LOGICAL FIRST(NZ)
      DIMENSION CMO(*),EIG(*),FAO(N2BBASX,NZ),FMO(*),FBUF(N2BBASX,NZ),
     &          VBUF(N2BBASX),TMAT(*),OP1INT(*),FTRI(NNBBASX,NZ),
     &          WORK(LWORK)
C
C     Initialize
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      DO IZ = 1,NZ
        FIRST(IZ) = .TRUE.
      ENDDO
C
C     MAKE ONE-ELECTRON FOCK MATRIX MINUS NUCLEAR ATTRACTION
C     ======================================================
C
      CALL DZERO(FBUF,N2BBASXQ)
C
C     Beta matrix
C     ===========
C
      IOP = IPBETAMT
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
C
C     Kinetic energy
C     ==============
C
      IOP = IPKINERG
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FBUF(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FBUF(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FBUF(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FBUF,NZ,WORK,LWORK)
C
C     KEEP NUCLEAR ATTRACTION IN SEPARATE BUFFER
C
      CALL DZERO(VBUF,N2BBASX)
C
C     Nuclear attraction integrals
C     ============================
C
      FIRST(1) = .TRUE.
      IOP = IPMOLFLD
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
C
C     Make full matrix
C     ================
C
      CALL DSPTSI(NTBAS(0),FTRI(1,1),VBUF)
      CALL BUTOBS(VBUF,1,WORK,LWORK)
C
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Get MO-transformation matrix
C     ============================
C
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','HEINZ')
      CALL READT(LUTMAT,N2TMT,TMAT)
      CLOSE(LUTMAT,STATUS='KEEP')
C
      A = ALPHAC
      N =  1
      DN = N
      K = -1
      IF(K.GT.0) THEN
        L = K
        J2 = 2*L - 1
      ELSE
        L = -(K+1)
        J2 = 2*L + 1
      ENDIF
      NREL = N - ABS(K)
C
      IOFF = NPSH(1)*NFBAS(1,0)+1
      DO IZ = 1,137
C
        Z = IZ
        Z2 = Z*Z
        ZM2 = 1.0D0/Z2
C
C     Analytical value of one-electron energy
C
        GAMMA = SQRT(K**2-(Z*A)**2)
        APPN = SQRT((NREL+GAMMA)**2+(Z*A)**2)
        EREL = -Z**2/((1+SQRT(1-(Z*A/APPN)**2))*APPN**2)
        ENREL = -DP5*(Z/N)*(Z/N)
C
C     MAKE CURRENT FOCK MATRIX
C     ========================
C
        CALL DCOPY(N2BBASXQ,FBUF,1,FAO,1)      
        CALL DAXPY(N2BBASX,Z,VBUF,1,FAO(1,1),1)
C
C     SOLVE WITHOUT PROJECTION
C
        CALL QTRANS('AOMO','S',D0,NFBAS(1,0),NFBAS(1,0),
     &              NORB(1),NORB(1),
     &              FAO,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &              FMO,NORB(1),NORB(1),NZ,IPQTOQ(1,0),
     &              TMAT,NFBAS(1,0),NORB(1),NZ,IPQTOQ(1,0),
     &              TMAT,NFBAS(1,0),NORB(1),NZ,IPQTOQ(1,0),
     &              WORK,LWORK,IPRGEN)
        IF(DOJACO) THEN
          CALL RSJACO(NORB(1),NORB(1),NORB(1),FMO,
     &                EIG,1,1,0,CMO,WORK,LWORK)
         ELSE
           CALL QDIAG(NZ,NORB(1),FMO,NORB(1),NORB(1),
     &             EIG,1,CMO,NORB(1),
     &             NORB(1),WORK,LWORK,IERR)
         ENDIF
         EFULL = EIG(NPSH(1)+1)
C
C     SOLVE WITH PROJECTION
C
        CALL QTRANS('AOMO','S',D0,NFBAS(1,0),NFBAS(1,0),
     &                         NESH(1),NESH(1),
     &              FAO,NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &              FMO,NESH(1),NESH(1),NZ,IPQTOQ(1,0),
     &              TMAT(IOFF),NFBAS(1,0),NFORB(1,0),NZ,IPQTOQ(1,0),
     &              TMAT(IOFF),NFBAS(1,0),NFORB(1,0),NZ,IPQTOQ(1,0),
     &              WORK,LWORK,IPRGEN)
        IF(DOJACO) THEN
          CALL RSJACO(NESH(1),NESH(1),NESH(1),FMO,
     &                EIG,1,1,0,CMO,WORK,LWORK)
        ELSE
          CALL QDIAG(NZ,NESH(1),FMO,NESH(1),NESH(1),
     &             EIG,1,CMO,NESH(1),NESH(1),WORK,LWORK,IERR)
        ENDIF
        EPROJ = EIG(1)
        EINT  = EIG(2)
        EFULL = EFULL
        EPROJ = EPROJ
        EREL  = EREL
        ENREL = ENREL
        WRITE(LUPRI,'(I3,4F20.10)') IZ,EFULL,EREL,ENREL,EINT
        IVAL = 0
        CALL PREIGV(EIG,NESH(1),IVAL,WORK,LWORK)
        CALL FLSHFO(LUPRI)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reno */
      SUBROUTINE RENO(LDM,NDIM,SMAT)
#include "implicit.h"
#include "priunit.h"
      DIMENSION SMAT(LDM,NDIM)
      DO J = 1,NDIM
        DO I = 1,NDIM
          SMAT(I,J) = SMAT(I,J)/SQRT(SMAT(I,I)*SMAT(J,J))
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphlow */
      SUBROUTINE SPHLOW(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,
     &                  EIG,CSPH,ICSP,LMAX,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin matrix with spherical transformation embedded.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(*),WORK(*)
C
      CALL QENTER('SPHLOW')
#include "memint.h"
      NXYZ = (LMAX*(LMAX+1))/2
      NFUN  = (LMAX*(LMAX+1)*(LMAX+2))/6
      NCBUF = ((LMAX+2)*(LMAX+1)*(3*LMAX**2+6*LMAX+1)*LMAX)/60
      NCAR = NBBAS(IREP,IC)
      NMAT = NBBAS(IREP,IC)*NBBAS(IREP,IC)
      CALL MEMGET('REAL',KTMAT,NMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCBUF,NCBUF ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KICBF,3*LMAX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFUN,3*NFUN,WORK,KFREE,LFREE)
      CALL SPHLO1(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,EIG,NCAR,
     &            CSPH,ICSP,LMAX,NXYZ,NFUN,WORK(KTMAT),WORK(KCBUF),
     &            WORK(KICBF),WORK(KIFUN),
     &            WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('SPHLOW',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPHLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphlo1 */
      SUBROUTINE SPHLO1(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,EIG,NCAR,
     &                  CSPH,ICSP,LMAX,MXYZ,NFUN,TMAT,CBUF,ICBF,
     &                  IFUN,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin matrix with spherical transformation embedded.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D6=6.0D0,D4=4.0D0)
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "ccom.h"
#include "symmet.h"
#include "pgroup.h"
#include "shells.h"
#include "nuclei.h"
#include "dgroup.h"
#include "dcbrkb.h"
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(LMAX),IFUN(NFUN,3),
     &          TMAT(NCAR,NCAR),CBUF(*),ICBF(3,LMAX),
     &          WORK(*)
      DATA COMP/'L','S'/
#include "sdpre.h"
#include "ibtfun.h"
C
C**************************************************
C***  Generate spherical transformation matrix  ***
C**************************************************
C
C     Initialize
C     ==========
C
      IF(IC.EQ.1) THEN
        JA1 = 1
        JA2 = NLRGSH
      ELSEIF(IC.EQ.2) THEN
        JA1 = NLRGSH+1
        JA2 = KMAX
      ENDIF
      N2MAT = NCAR*NCAR
      CALL DZERO(TMAT,N2MAT)
C
C     Initial print section
C     =====================
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A,A3,2X,A1)') 
     &    'SPHLOW: Transformation matrix for ',REP(IREP),COMP(IC)
      ENDIF
C
C     Loop over shells
C     ================
C
      ISOFF = ICOS(IREP+1,IC)
      NAORB = NSTRT(JA1)
      INUC = 0
      IROW = 0
      ICOL = 0
      DO 10 JA = JA1,JA2
        JNUC = NCENT(JA)
C
C       New center; check out allowed symmetries
C       ----------------------------------------
        IF(JNUC.NE.INUC) THEN
          MULA = ISTBNU(JNUC)
          INUC = JNUC
          ILL  = 0
          IOFF = 0
        ENDIF
        LL  = NHKT(JA)
C
C       New L-value; extract relevant transformation matrix
C       ---------------------------------------------------
        IF(LL.NE.ILL) THEN
          ILL  = LL
          L    = LL - 1
          NLM  = 2*L + 1
          NXYZ = (LL*(LL+1))/2
          LOFF = NHKOFF(LL)
C
C         Cartesian functions
C
          KXYZ = 0
          DO ICMP = 1,NXYZ
            IVARB = IBTXOR(IREP,ISYMAO(LL,ICMP))
            IF(IBTAND(MULA,IVARB).EQ.0) THEN
              KXYZ              = KXYZ + 1
              IFUN(LOFF+KXYZ,1) = ICMP
              IFUN(LOFF+ICMP,3) = KXYZ
            ENDIF
            ICBF(1,LL) = KXYZ
          ENDDO
C
C         Solid harmonic functions
C
          KLM = 0
          DO ICMP = 1,NLM
            IVARB = IBTXOR(IREP,IREPLM(L,ICMP-LL))
            IF(IBTAND(MULA,IVARB).EQ.0) THEN
              KLM               = KLM + 1
              IFUN(LOFF+KLM,2)  = ICMP
            ENDIF
            ICBF(2,LL) = KLM
          ENDDO
C
C         Extract spherical transformation matrix
C
          ICBF(3,LL) = IOFF
          IF(KLM.GT.0) THEN
            JOFF = ICSP(LL)
            CALL EXTELM('N',CSPH(JOFF+1),NXYZ,NLM,
     &                      CBUF(IOFF+1),KXYZ,KXYZ,KLM,
     &                      IFUN(LOFF+1,1),IFUN(LOFF+1,2))
            JOFF = IOFF
            DO J = 1,KLM
              FAC = D1/DNRM2(KXYZ,CBUF(JOFF+1),1)
              CALL DSCAL(KXYZ,FAC,CBUF(JOFF+1),1)
              JOFF = JOFF + KXYZ
            ENDDO
            IOFF = IOFF + KXYZ*KLM
          ENDIF
          IF(IPRINT.GE.4) THEN
            KOFF = ICBF(3,LL)
            WRITE(LUPRI,'(A,I3)')
     &       ' * Center : ',JNUC,
     &       ' * L-value: ',L
            IF(KLM.GT.0) THEN
              CALL OUTPUT(CBUF(KOFF+1),1,KXYZ,1,KLM,KXYZ,KLM,1,LUPRI)
            ELSE
              WRITE(LUPRI,'(A)') ' - no functions - '
            ENDIF
          ENDIF
        ENDIF
C
C       Expand transformation matrix
C       ============================
C
        IF(KLM.EQ.0) GOTO 10
        KOFF = ICBF(3,LL)+1
        DO J = 1,KLM
          CALL DCOPY(KXYZ,CBUF(KOFF),1,
     &                TMAT(IROW+1,ICOL+J),1)   
          KOFF = KOFF + KXYZ
        ENDDO
        IROW = IROW + KXYZ
        ICOL = ICOL + KLM
   10 CONTINUE
C
C     Final print section
C
      IF(IPRINT.GE.3) THEN
        CALL PRSYMB(LUPRI,'-',40,0)
        CALL OUTPUT(TMAT,1,IROW,1,ICOL,NCAR,ICOL,1,LUPRI)
      ENDIF
C
C*********************************************************
C***  Transform overlap matrix to spherical basis and  ***
C***  generate Lowdin matrix                           ***
C*********************************************************
C
      NSPH = ICOL
      CALL QTRANS('AOMO','S',D0,NCAR,NCAR,NSPH,NSPH,
     &            SMAT,LDM ,NCAR,1,IPQTOQ(1,0),
     &            VMAT,LDM ,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRINT)
      CALL LOWGEN(VMAT,NSPH,SMAT,NEFF,LDM,STOL,IPRINT,
     &            EIG,WORK,LWORK)
      CALL DGEMM('N','N',NCAR,NEFF,NSPH,D1,TMAT,NCAR,
     &           SMAT,LDM,D0,VMAT,LDM)
      NA = NCAR - NSPH
      NB = NSPH - NEFF 
      ND = NA + NB
      WRITE(LUPRI,'(3X,A1,3X,A3,3X,3(A,I10),A)') 
     &           COMP(IC),REP(IREP),'* Deleted: ',ND,
     &           '(Proj: ',NA,', Lindep: ',NB,')'
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkblow */
      SUBROUTINE RKBLOW(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,
     &                  EIG,CSPH,ICSP,LMAX,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin matrix with RKB transformation embedded.
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(*),WORK(*)
C
      CALL QENTER('RKBLOW')
#include "memint.h"
      NXYZ = (LMAX*(LMAX+1))/2
      NFUN  = (LMAX*(LMAX+1)*(LMAX+2))/6
      NCBUF = ((LMAX+2)*(LMAX+1)*(3*LMAX**2+6*LMAX+1)*LMAX)/60
      IF(LMAX.GE.3) THEN
        L2    = LMAX - 2
        NCBUF = NCBUF + ((L2+2)*(L2+1)*(3*L2**2+6*L2+1)*L2)/60
      ENDIF
      NCAR = NBBAS(IREP,IC)
      NMAT = NBBAS(IREP,IC)*NBBAS(IREP,IC)
      CALL MEMGET('REAL',KTMAT,NMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCBUF,NCBUF ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KICBF,3*LMAX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFUN,3*NFUN,WORK,KFREE,LFREE)
      CALL RKBLO1(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,EIG,NCAR,
     &            CSPH,ICSP,LMAX,NXYZ,NFUN,WORK(KTMAT),WORK(KCBUF),
     &            WORK(KICBF),WORK(KIFUN),
     &            WORK(KFREE),LFREE,IPRINT)
      CALL MEMREL('RKBLOW',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RKBLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkblo1 */
      SUBROUTINE RKBLO1(IREP,IC,SMAT,VMAT,LDM,NEFF,STOL,EIG,NCAR,
     &                  CSPH,ICSP,LMAX,MXYZ,NFUN,TMAT,CRED,IRED,
     &                  IFUN,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     Generate Lowdin matrix with RKB transformation embedded.
C
C     Spherical transformation matrix:
C       CSPH - for each L value it has the dimension CSPH(NLM,NXYZ,*)
C       ICSP - gives offsets to each L value in spherical 
C              transformation matrix CSPH
C     For each center and L value a transformation matrix is extracted
C     from CSPH according to allowed symmetries.
C       IFUN(*,1) - points to Cartesian component in full list
C       IFUN(*,2) - points to spherical component in full list
C       IFUN(*,3) - points to Cartesian comp. in symmetry-reduced list
C     The symmetry-reduced spherical transformation matrix is
C     stored in CRED
C     For each L value
C       IRED(1,*) - gives number of symmetry-reduced Cartesians KXYZ
C       IRED(2,*) - gives number of symmetry-reduced sphericals KLM
C       IRED(3,*) - gives offsets to each L value
C     For each L value CRED has the dimension CRED(KXYZ,KLM)
C     The RKB part is adjoined to CRED and has the dimension
C        CRED(KRKB,KLMM2) where KRKB = KXYZ+KXYZM2
C     KXYZM2 and KLMM2 are the number of symmetry-reduced Cartesian
C     and spherical components, respectively, of the L-2 partner.
C
C      

C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D6=6.0D0,D4=4.0D0)
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "ccom.h"
#include "symmet.h"
#include "pgroup.h"
#include "shells.h"
#include "nuclei.h"
#include "dgroup.h"
#include "dcbrkb.h"
      LOGICAL LRKB
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(LMAX),IFUN(NFUN,3),
     &          TMAT(NCAR,NCAR),CRED(*),IRED(3,LMAX),
     &          WORK(*)
      DATA COMP/'L','S'/
#include "sdpre.h"
#include "ibtfun.h"
C
C**************************************************
C***  Generate spherical transformation matrix  ***
C**************************************************
C
C     Initialize
C     ==========
C
      IF(IC.EQ.1) THEN
        JA1 = 1
        JA2 = NLRGSH
      ELSEIF(IC.EQ.2) THEN
        JA1 = NLRGSH+1
        JA2 = KMAX
      ENDIF
      N2MAT = NCAR*NCAR
      CALL DZERO(TMAT,N2MAT)
C
C     Initial print section
C     =====================
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A,A3,2X,A1)') 
     &    'RKBLOW: Transformation matrix for ',REP(IREP),COMP(IC)
      ENDIF
C
C     Loop over shells
C     ================
C     (Shell number will be different from primitive number
C      if large component is contracted, thus the double
C      loop structure)
C
C     Offset in SO-basis:
      ISOFF = ICOS(IREP+1,IC)
C     Offset in AO-basis:
      NAORB = NSTRT(JA1)
C
      INUC = 0
      ICOL = 0
      JP1_TEST = 0
      DO 10 JA = JA1,JA2
         JP1  = JSTRT(JA) + 1
C
C        Skip primitives allready processed in another shell.
C
         IF ( JP1_TEST .EQ. JP1 ) GO TO 10
         JP1_TEST = JP1
         JP2  = JSTRT(JA) + NUCO(JA)
C
         JNUC = NCENT(JA)
C
C        New center; check out allowed symmetries
C        ----------------------------------------
         IF(JNUC.NE.INUC) THEN
            MULA = ISTBNU(JNUC)
            INUC = JNUC
            ILL  = 0
            IOFF = 0
         ENDIF
         LL  = NHKT(JA)
C
C        New L-value; extract relevant transformation matrix
C        ---------------------------------------------------
         IF(LL.NE.ILL) THEN
            ILL  = LL
            JRKB = 0
            L    = LL - 1
            NLM  = 2*L + 1
            NXYZ = (LL*(LL+1))/2
            LOFF = NHKOFF(LL)
C
C           Cartesian functions
C
            KXYZ = 0
            DO ICMP = 1,NXYZ
               IVARB = IBTXOR(IREP,ISYMAO(LL,ICMP))
               IF(IBTAND(MULA,IVARB).EQ.0) THEN
                  KXYZ              = KXYZ + 1
                  IFUN(LOFF+KXYZ,1) = ICMP
                  IFUN(LOFF+ICMP,3) = KXYZ
               ENDIF
               IRED(1,LL) = KXYZ
            ENDDO
C
C           Solid harmonic functions
C
            KLM = 0
            DO ICMP = 1,NLM
               IVARB = IBTXOR(IREP,IREPLM(L,ICMP-LL))
               IF(IBTAND(MULA,IVARB).EQ.0) THEN
                  KLM               = KLM + 1
                  IFUN(LOFF+KLM,2)  = ICMP
               ENDIF
               IRED(2,LL) = KLM
            ENDDO
C
C           Extract spherical transformation matrix
C
            IRED(3,LL) = IOFF
            IF(KLM.GT.0) THEN
               JOFF = ICSP(LL)
               CALL EXTELM('N',CSPH(JOFF+1),NXYZ,NLM,
     &              CRED(IOFF+1),KXYZ,KXYZ,KLM,
     &              IFUN(LOFF+1,1),IFUN(LOFF+1,2))
               JOFF = IOFF
               DO J = 1,KLM
                  FAC = D1/DNRM2(KXYZ,CRED(JOFF+1),1)
                  CALL DSCAL(KXYZ,FAC,CRED(JOFF+1),1)
                  JOFF = JOFF + KXYZ
               ENDDO
               IOFF = IOFF + KXYZ*KLM
            ENDIF
            IF(IPRINT.GE.4) THEN
               KOFF = IRED(3,LL)
               WRITE(LUPRI,'(A,I3)')
     &              ' * Center : ',JNUC,
     &              ' * L-value: ',L
               IF(KLM.GT.0) THEN
                 CALL OUTPUT(CRED(KOFF+1),1,KXYZ,1,KLM,KXYZ,KLM,1,LUPRI)
               ELSE
                 WRITE(LUPRI,'(A)') ' - no functions - '
               ENDIF
            ENDIF
         ENDIF
C
C        Expand transformation matrix
C        ============================
C
         IF(KLM.EQ.0) GOTO 10
C
C        Loop over primitives
C      
         DO 15 JP = JP1, JP2
C
C           A function can contribute to kinetic balance in two ways :
C           a. It balances a (L-1) large component function
C           b. It balances a (L+1) large component function
C
C           In case b it needs to be combined with a (L+2) function,
C           and we take this combination when we enter this (L+2) shell.
C           This means we can disgard the functions for the moment,
C           provided that case a. does not apply as well (as may happen
C           with family-type basis sets).
C
C           Look for (L-2) small function. If it does not exist, 
C           IRKB returns the pointer to (L-1) large component function
C           (case b). If this is zero it means that we have the
C           first function for case b and we may skip it.
C
            JB = IRKB(JP)
            IF (JB.EQ.0) GOTO 15
C
            IROW = NSTRT(JA + JP - JP1) + IFUN(LOFF+1,1)
C                        ^- to get correct shell number
C
C           JA is the pointer to the start of this series of uncontracted
C           shells, increase it by the distance from the start
C
            IROW = IPTSYM(IROW,IREP) - ISOFF
C
C           Fill in standard spherical transformation matrix
C
            KOFF = IRED(3,LL)+1
            DO J = 1,KLM
               CALL DCOPY(KXYZ,CRED(KOFF),1,
     &              TMAT(IROW,ICOL+J),1)   
               KOFF = KOFF + KXYZ
            ENDDO
            ICOL = ICOL + KLM
C
C           Possibly add RKB part
C
            LRKB = L.GE.2
            IF(LRKB) THEN
               LLM2   = LL-2
               KLMM2  = IRED(2,LLM2)
               LRKB   = KLMM2.GT.0
            ENDIF
            IF(LRKB) THEN
C
C           Possibly generate RKB matrix
C
C
               IF(JRKB.EQ.0) THEN
                  JRKB = 1
                  NXYZM2 = (LLM2*(LLM2+1))/2
                  NLMM2  = 2*LLM2-1
                  KXYZM2 = IRED(1,LLM2)                            
                  KRKB   = KXYZ + KXYZM2
                  IM2    = IOFF + KXYZ
                  LM2F   = NHKOFF(LLM2)
                  CALL DZERO(CRED(IOFF+1),KRKB*KLMM2)
C
C                 First extract spherical transformation matrix 
C                 for L-2 partner
                  JM2    = ICSP(LLM2)
                  CALL EXTELM('N',CSPH(JM2+1) ,NXYZM2,NLMM2,
     &                 CRED(IM2+1),KRKB,KXYZM2,KLMM2,
     &                 IFUN(LM2F+1,1),IFUN(LM2F+1,2))
C
C                 Then fill in corresponding number for r**2
C
                  DO K = 1,KXYZM2
                     IJ = IFUN(LM2F+K,1)
                     IA = INT(SQRT(SDPRE(2*IJ) + 0.25D0) + 0.4999D0)
                     IB = IJ - IA*(IA-1)/2
C**               x**2: (ia  ,ib  )
                     IJX  = (IA*(IA-1))/2 + IB
                     KOFF = IOFF + IFUN(LOFF+IJX,3)
                     CALL DAXPY(KLMM2,D1,CRED(IM2+K),KRKB,
     &                    CRED(KOFF),KRKB)
C**               y**2: (ia+2,ib  )
                     IJY = ((IA+2)*(IA+1))/2 + IB
                     KOFF = IOFF + IFUN(LOFF+IJY,3)
                     CALL DAXPY(KLMM2,D1,CRED(IM2+K),KRKB,
     &                    CRED(KOFF),KRKB)
C**               z**2: (ia+2,ib+2)
                     IJZ = ((IA+2)*(IA+1))/2 + IB + 2
                     KOFF = IOFF + IFUN(LOFF+IJZ,3)
                     CALL DAXPY(KLMM2,D1,CRED(IM2+K),KRKB,
     &                    CRED(KOFF),KRKB)
                  ENDDO
C
C                 Scale matrix for L-2 partner
C
                  FAC =  D6-D4*SDPRE(LL)
                  KOFF = IM2+1
                  DO K = 1,KLMM2
                     CALL DSCAL(KXYZM2,FAC,CRED(KOFF),1)
                     KOFF = KOFF + KRKB
                  ENDDO
C
C                 Normalize...
C
                  KOFF = IOFF+1
                  DO J = 1,KLMM2
                     FAC = D1/DNRM2(KRKB,CRED(KOFF),1)
                     CALL DSCAL(KRKB,FAC,CRED(KOFF),1)
                     KOFF = KOFF + KRKB
                  ENDDO
                  IF(IPRINT.GE.4) THEN
                     WRITE(LUPRI,'(A)') '* RKB part:'
                     CALL OUTPUT(CRED(IOFF+1),1,KRKB,1,KLMM2,
     &                    KRKB,KLMM2,1,LUPRI)
                  ENDIF
                  IOFF = IOFF + KRKB*KLMM2
               ENDIF
C
C              Copy extra functions for L
C            
               IRM2 = NSTRT(JB) + IFUN(LM2F+1,1)
               IRM2 = IPTSYM(IRM2,IREP) - ISOFF
               KOFF = IRED(3,LL)+KXYZ*KLM+1
               DO J = 1,KLMM2
                  CALL DCOPY(KXYZ,CRED(KOFF),1,
     &                 TMAT(IROW,ICOL+J),1)   
                  KOFF = KOFF + KXYZ
                  CALL DCOPY(KXYZM2,CRED(KOFF),1,
     &                 TMAT(IRM2,ICOL+J),1)   
                  KOFF = KOFF + KXYZM2
               ENDDO
               IROW = IROW + KXYZ
               ICOL = ICOL + KLMM2
            ENDIF
 15      CONTINUE
 10   CONTINUE
C
C     Final print section
C
      IF(IPRINT.GE.3) THEN
         CALL PRSYMB(LUPRI,'-',40,0)
         CALL OUTPUT(TMAT,1,NCAR,1,ICOL,NCAR,ICOL,1,LUPRI)
      ENDIF
C
C*********************************************************
C***  Transform overlap matrix to spherical basis and  ***
C***  generate Lowdin matrix                           ***
C*********************************************************
C
      NSPH = ICOL
      CALL QTRANS('AOMO','S',D0,NCAR,NCAR,NSPH,NSPH,
     &            SMAT,LDM ,NCAR,1,IPQTOQ(1,0),
     &            VMAT,LDM ,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            TMAT,NCAR,NSPH,1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRINT)
      CALL LOWGEN(VMAT,NSPH,SMAT,NEFF,LDM,STOL,IPRINT,
     &     EIG,WORK,LWORK)
      CALL DGEMM('N','N',NCAR,NEFF,NSPH,D1,TMAT,NCAR,
     &     SMAT,LDM,D0,VMAT,LDM)
      NA = NCAR - NSPH
      NB = NSPH - NEFF 
      ND = NA + NB
      WRITE(LUPRI,'(3X,A1,3X,A3,3X,3(A,I10),A)') 
     &     COMP(IC),REP(IREP),'* Deleted: ',ND,
     &     '(Proj: ',NA,', Lindep: ',NB,')'
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowgen */
      SUBROUTINE LOWGEN(SMAT,NS,VMAT,NV,LDM,STOL,IPRINT,
     &                  EIG,WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix
C     (Lowdin's canonical orthonormalization)
C
C     V_ij = U_ij/SQRT(s_j) where U_ij is the matrix that diagonalizes
C     the overlap matrix S and s_j is eigenvalue j.
C     Linear dependence is removed by eliminating columns of V
C     corresponding to eigenvalues of the overlap matric beyond
C     the given threshold.
C
C     Written by T.Saue Sep 11 1995
C     Last revision Apr 1 1998 - tsaue
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),WORK(LWORK)
C
      CALL QENTER('LOWGEN')
C
C     Diagonalize overlap matrix
C
      CALL RSJACO(LDM,NS,NS,SMAT,EIG,1,-1,0,VMAT,WORK,LWORK)
      NV = 0
      DO I = 1,NS
        IF(EIG(I).LE.STOL) GOTO 10
        NV  = NV + 1
        FAC = D1/SQRT(EIG(I))
        CALL DSCAL(NS,FAC,VMAT(1,I),1)
      ENDDO
 10   CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.2) WRITE(LUPRI,'(7X,A,1P,D8.1)')
     &   'LOWGEN: Smallest ".LINDEP" test value of a kept orbital:',
     &   EIG(NV)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Lowdin input matrix:',-1)
        CALL OUTPUT(SMAT,1,NS,1,NS,LDM,NS,1,LUPRI)
      ENDIF
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('LOWGEN: Eigenvalues in Lowdin diagonalization:',-1)
        CALL OUTPUT(EIG,1,NS,1,1,NS,1,1,LUPRI)
      ENDIF
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Canonical Lowdin matrix:',-1)
        CALL OUTPUT(VMAT,1,NS,1,NV,LDM,NS,1,LUPRI)
      ENDIF
      CALL QEXIT('LOWGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck tscorr */
      SUBROUTINE TSCORR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Generate SS nuclear attraction integrals to mimic
C     two-electron SS-integrals
C
C     Written by T.Saue July 1997
C     Last revision July 15 1997 - tsaue
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER(D1 = 1.0D0)
C
      DIMENSION WORK(LWORK),CBUF(2),GBUF(2),FMAT(*)
      LOGICAL DOINT(2,2)
      CHARACTER OMITVNUC(2)*4
#include "dgroup.h"
#include "dcbbas.h"
#include "nuclei.h"
C
      CALL QENTER('TSCORR')
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEMT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Buffer info
C
CTROND
      DO I = 1,NUCIND
        CBUF(I) = CHARGE(I)
        GBUF(I) = GNUEXP(I)
        GNUEXP(I) = 0.0D0
      ENDDO
      call quit('TSCORR not properly implemented in this version')
C     HJAaJ: this is only for a specific datomic!!! (TeSr ??)
C
      CHARGE(1) = -0.1812D0
      CHARGE(2) = -0.0832D0
C
C     Generate integrals
C
      DOINT(1,1) = .FALSE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.

      OMITVNUC(1) = 'TTTT'
      OMITVNUC(2) = 'FFFF'
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,OMITVNUC,IPRINT)
C
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('SS nuclear attraction integrals',-1)
          CALL PRQMAT(WORK(KONEMT),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      ENDIF
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
C
C     Memory deallocation
      CALL MEMREL('TSCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Buffer retrieve
C
      DO I = 1,NUCIND
        CHARGE(I) = CBUF(I) 
        GNUEXP(I) = GBUF(I) 
      ENDDO
C
      CALL QEXIT('TSCORR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bncorr */
      SUBROUTINE BNCORR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     BNCORR: Generate screened bare nucleus matrix
C     Jun 2003 Jesper Kielberg Pedersen (based on OLDBNCR)
C
C     Description : 
C     Extension of BNCORR that in a simple way takes the shell-structure
C     into account. For each center in a molecule we estimate the screening
C     neglected in the bare nucleus approximation by a sum of the 
C     contributions from each shells (n) in the atoms :
C           sum_(n) <X_A|(-Z(n)*G(n,a,r_C))/r_C|X_B> 
C           X=L,S ; G(a,r)=exp(-a*r^2)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL DOINT(2,2), FINISH
      CHARACTER*4 OMITVNUC
      DIMENSION WORK(LWORK),FMAT(*)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('BNCORR')
c
c     Start timer
c
      TIMSTR = SECOND() 
c
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEMT, N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN, NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KESHEL,7*NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KCHRG),NUCIND)
      CALL DZERO(WORK(KONEMT),N2BBASX)
      CALL DZERO(WORK(KESHEL),NUCIND)
C
C     Initialization (calculate LL and SS)
C
      DOINT(1,1) = .TRUE.
      IF (NOSMLV) THEN
         DOINT(2,2) = .FALSE.
      ELSE
         DOINT(2,2) = .TRUE.
      END IF
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
c
c   ... find number of occupied shells
c
      Z = D0
      DO I = 1,NUCIND
         Z = MAX(Z,CHARGE(I))
      ENDDO
      CALL BNCOR1(GNU,ZEFF,0,Z,NSHELL,0,DUMMY,NUCIND,.TRUE.)
c
c   ... loop over shells
c       in the first shell (N=1) we only place one electron to avoid
c       self-interaction. This electron is delocalized (to some degree)
c       in the zeroth run (N=0) by using a diffuse Gaussian (see BNCOR1)
c
      DO N = 0,NSHELL
c     ... In Slater's rules we write the electron configuration as
c         [1s] [2s2p] [3s3p][3d] [4s4p][4d][4f]... so for each shell we have
c          1     1        2            3           loops
         IF (N.LE.1) THEN
	    NSLOOP = 1
         ELSE
	    NSLOOP = N - 1
         ENDIF
c        ... we will not go beyond f-orbitals
         NSLOOP = MIN(NSLOOP,3)
         DO J = 1,NSLOOP
         ZTEST = D0
         DO I = 1,NUCIND
c      ... Get charges and exponents for the electrons in this group
           CALL BNCOR1(GNU,ZEFF,N,CHARGE(I),J,I,WORK(KESHEL),
     &          NUCIND,.FALSE.)
           ZTEST = MAX(ZTEST,ZEFF)
           WORK(KXGEXP-1+I) =  GNU
           WORK(KCHRG -1+I) = -ZEFF
         ENDDO
         IF(ZTEST.EQ.D0) GOTO 20
c
         IF(IPRINT .GE. 3) THEN
            IF(N.EQ.0) THEN
            WRITE(LUPRI,5) 
     &' *** Entering BNCORR  for <X_A|(-Z*G(a,r_C))/r_C|X_B> terms ***',
     &'     X=L,S ; G(a,r)=exp(-a*r^2)'
            WRITE(LUPRI,'(14X,A)')
     &      ' Nuclei      Charge(-Z)       Exponent(a)',
     &      '-----------------------------------------'
            ENDIF
            WRITE(LUPRI,'(14X,A,I2,A)')'(Loop over ',N,' . shell)'
            DO I = 1,NUCIND
               IF(WORK(KCHRG-1+I).LT.D0)
     &         WRITE(LUPRI,'(15X,A5,F15.2,F15.2)')
     &         NAMN(I),WORK(KCHRG-1+I),WORK(KXGEXP-1+I)
            ENDDO
            IF (N.EQ.NSHELL)
     &      WRITE(LUPRI,'(14X,A/)')
     &      '-----------------------------------------'
5        FORMAT(/,2(/,4X,A))
         ENDIF
C
C     Omit selected terms depending on OMITVNUC (See CATNUC for details)
C     In this version one-center terms are calculated with the same
C     charges/exponents as multi-center terms.
C
         OMITVNUC = 'FFFF'
c
c     ... Make list of nuclei.
c
         CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals (both LL and SS)
c
         CALL CATDR1(WORK(KONEMT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
c  ... end loop within shells
20    CONTINUE
         ENDDO
c  ... end loop over shells
      ENDDO
10    CONTINUE
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KONEMT),WORK(KONEIN),1,NNBBASX)
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.5) THEN
         CALL HEADER('BNCORR correction FMAT',-1)
         CALL OUTPUT(WORK(KONEMT),1,NTBAS(0),1,NTBAS(0),
     &        NTBAS(0),NTBAS(0),1,LUPRI)
      END IF
C
C     Add correction to Fock-matrix
C
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
      IF(IPRINT.GE.5) THEN
         CALL HEADER('BNCORR: screening corrected bare nucleus FMAT',-1)
         CALL PRQMAT(FMAT,NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Memory deallocation
      CALL MEMREL('BNCORR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      IF (IPRINT.GE.2)
     &  CALL TIMTXT('>>> Time used in BNCORR is',TIME,LUPRI)
      CALL QEXIT('BNCORR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bncor1 */
      SUBROUTINE BNCOR1(GNU,ZEFF,N,Z,J,INUC,SHELE,NNUC,MAXSHELL)
C*****************************************************************************
C
C     BNCOR1 :  Find effektive charges and Gaussian exponent 
C               (from Slater rules) for a given shell (N) for
C               a given atom (Z)
C
C     Jun 2003 Jesper Kielberg Pedersen (used for BNCORR)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D4 = 4.0D0,
     &          D6 = 6.0D0, D8 = 8.0D0, D10=10.0D0, D14=14.0D0,
     &          D18=18.0D0)
      LOGICAL MAXSHELL
      DIMENSION STO1G(7), SHELE(NNUC,7)
      SAVE STO1G
C     STO1G(n) is the exponent of the 1s-Gaussian that
C     best fits the ns Slater-type atomic orbitals of exponent 1.
C     (Robert F. Stewart, JCP 52(1),431 (Jan 1970))
C     This article provides exponent up to the n=5 shell. The exponents
C     for the last two shells have been generated from a fit :
C     (y = 0.28196 x^-1.5532)
      DATA STO1G /0.2709,0.1012,0.05297,0.03265,0.02217,0.01744,
     &            0.01373/
C
C     Table of s,p,d,f orbital occupations of the elements
C     (from which the electron-configuration can be deduced)
      DIMENSION NOCC(118*4)
      SAVE NOCC
      DATA NOCC 
C     Periode 1
     &          / 1, 0, 0, 0,
     &            2, 0, 0, 0,
C     Periode 2
     &            3, 0, 0, 0,
     &            4, 0, 0, 0,
     &            4, 1, 0, 0,
     &            4, 2, 0, 0,
     &            4, 3, 0, 0,
     &            4, 4, 0, 0,
     &            4, 5, 0, 0,
     &            4, 6, 0, 0,
C     Periode 3
     &            5, 6, 0, 0,
     &            6, 6, 0, 0,
     &            6, 7, 0, 0,
     &            6, 8, 0, 0,
     &            6, 9, 0, 0,
     &            6,10, 0, 0,
     &            6,11, 0, 0,
     &            6,12, 0, 0,
C     Periode 4
     &            7,12, 0, 0,
     &            8,12, 0, 0,
     &            8,12, 1, 0,
     &            8,12, 2, 0,
     &            8,12, 3, 0,
     &            8,12, 4, 0,
     &            8,12, 5, 0,
     &            8,12, 6, 0,
     &            8,12, 7, 0,
     &            8,12, 8, 0,
     &            8,12, 9, 0,
     &            8,12,10, 0,
     &            8,13,10, 0,
     &            8,14,10, 0,
     &            8,15,10, 0,
     &            8,16,10, 0,
     &            8,17,10, 0,
     &            8,18,10, 0,
C     Periode 5
     &            9,18,10, 0,
     &           10,18,10, 0,
     &           10,18,11, 0,
     &           10,18,12, 0,
     &           10,18,13, 0,
     &           10,18,14, 0,
     &           10,18,15, 0,
     &           10,18,16, 0,
     &           10,18,17, 0,
     &           10,18,18, 0,
     &           10,18,19, 0,
     &           10,18,20, 0,
     &           10,19,20, 0,
     &           10,20,20, 0,
     &           10,21,20, 0,
     &           10,22,20, 0,
     &           10,23,20, 0,
     &           10,24,20, 0,
C     Periode 6
     &           11,24,20, 0,
     &           12,24,20, 0,
     &           12,24,21, 0,
C    --> Lanthanides & actinides
     &           12,24,21, 1,
     &           12,24,20, 3,
     &           12,24,20, 4,
     &           12,24,20, 5,
     &           12,24,20, 6,
     &           12,24,20, 7,
     &           12,24,20, 8,
     &           12,24,20, 9,
     &           12,24,20,10,
     &           12,24,20,11,
     &           12,24,20,12,
     &           12,24,20,13,
     &           12,24,20,14,
     &           12,24,21,14,
C    <-- Lanthanides & actinides
     &           12,24,22,14, 12,24,23,14, 12,24,24,14,
     &           12,24,25,14, 12,24,26,14, 12,24,27,14,
     &           12,24,28,14, 12,24,29,14, 12,24,30,14,
     &           12,25,30,14, 12,26,30,14, 12,27,30,14,
     &           12,28,30,14, 12,29,30,14, 12,30,30,14,
C     Periode 7
     &           13,30,30,14,
     &           14,30,30,14,
     &           14,30,31,14,
C    --> Lanthanides & actinides
     &           14,30,32,14, 14,30,31,16, 14,30,31,17,
     &           14,30,31,18, 14,30,30,20, 14,30,30,21,
     &           14,30,31,21, 14,30,30,23, 14,30,30,24,
     &           14,30,30,25, 14,30,30,26, 14,30,30,27,
     &           14,30,30,28, 14,30,31,28,
C    <-- Lanthanides & actinides
     &           14,30,32,28, 14,30,33,28, 14,30,34,28,
     &           14,30,35,28, 14,30,36,28, 14,30,37,28,
     &           13,30,39,28, 13,30,40,28, 14,30,40,28,
     &           14,31,40,28, 14,32,40,28, 14,33,40,28,
     &           14,34,40,28, 14,35,40,28, 14,36,40,28/
      CALL QENTER('BNCOR1')
c
      IF (NINT(Z) .GT. 118) THEN
         WRITE(LUPRI,'(1X,2(A/))') 'FATAL : The bare nucleus '//
     &   'correction is only implemented for elemtents 1-118',
     &   '         You can disable BNCORR with the .NOBNCR keyword '//
     &   'in *DHFCAL'
         CALL QUIT('BNCORR : Unsupported element!')
      ENDIF
      IF (MAXSHELL) THEN
c    ... just find maximum number of shells we need to run over 
c        and return number in J-variable (used in BNCORR)
         NST = NOCC((NINT(Z)-1)*4 + 1)
         J = NST/2
         IF(MOD(NST,2).NE.0) J = J + 1
         GOTO 10
      ENDIF
C
      ZEFF = D0
c     ... If this is a ghost-atom (f.ex for Counter-Poise 
c         correction) we return (charge is zero in BNCORR)
      IF (Z .EQ. D0) GOTO 10
      IF (N.EQ.0) THEN
c     ...Always make 1 electron delocalized
         ZEFF = D1
         GNU  = 0.1D0
         SHELE(INUC,1) = SHELE(INUC,1) + D1
         GOTO 10
      ENDIF
c     ... If Hydrogen we are done
      IF (Z.EQ.D1) GOTO 10
c     ... special case case for the (2-1) 1s electron
      IF (N.EQ.1) THEN
         ZEFF = D1
         SEXP = Z - 0.35
         IF (Z.EQ.2) SEXP = SEXP + 0.05
         SHELE(INUC,1) = SHELE(INUC,1) + D1
         GOTO 20
      ENDIF
c     ... Generel case. We are past H and He and have thus accounted 
c         for two electron.
c     ... Find number of s-electrons for element (Z)
      NST = NOCC((NINT(Z)-1)*4 + 1)
      NPT = NOCC((NINT(Z)-1)*4 + 2)
      NDT = NOCC((NINT(Z)-1)*4 + 3)
      NFT = NOCC((NINT(Z)-1)*4 + 4)
c
c     ... Now figure out the occupation used for the Slater rules (N >= 2) :
c      [1s][2s2p][3s3p][3d][4s4p][4d][4f]...
c
c     =================================
c     How many s-electron in this shell
c     =================================
c
      NS = NST - (N-1)*D2
c
      IF (NS .LE. D0) THEN
         NS = D0
      ELSEIF (NS.GE.D2) THEN
         NS = D2
      ENDIF
c
c     =================================
c     How many p-electron in this shell
c     =================================
c
      NP = NPT - (N-2)*D6
c
      IF (NP .LE. D0) THEN
         NP = D0
      ELSEIF (NP.GE.D6) THEN
         NP = D6
      ENDIF
c
c     =================================
c     How many d-electron in this shell
c     =================================
c
      ND = NDT
      DO I  = 1,N-3
         ND = ND - D10
      ENDDO
c
      IF (ND .LE. D0) THEN
         ND = D0
      ELSEIF (ND.GE.D10) THEN
         ND = D10
      ELSE
         ND = NDT - (N-3)*D10
      ENDIF
c
c     =================================
c     How many f-electron in this shell
c     =================================
c
      NF = NFT
      DO I  = 1,N-4
         NF = NF - D14
      ENDDO
c
      IF (NF .LE. D0) THEN
         NF = D0
      ELSEIF (NF.GE.D14) THEN
         NF = D14
      ELSE
         NF = NFT - (N-4)*D14
      ENDIF
      IF (J.EQ.1) ZEFF = NS + NP
      IF (J.EQ.2) ZEFF = ND
      IF (J.EQ.3) ZEFF = NF
      SHELE(INUC,N) = SHELE(INUC,N) + ZEFF
      IF (ZEFF.EQ.D0) THEN
         GNU  = D0
         GOTO 10
      ENDIF
C
C     ... Finally use Slater's rules to estimate screening 
C         (and orbital exponent)
C
      IF (J.EQ.1) THEN
c        [nsnp]
         EMN = D0
         DO I = 1,N-2
            EMN = EMN + SHELE(INUC,I)
         ENDDO
         SEXP = (Z - ((ZEFF-1.0)*0.35 + SHELE(INUC,N-1)*0.85 + EMN))/N
      ELSE
c        [nd],[nf]
         EMN = D0
         DO I = 1,N-1
            EMN = EMN + SHELE(INUC,I)
         ENDDO
         SEXP = (Z - ((ZEFF-1.0)*0.35 + EMN))/N
      ENDIF
c
20    CONTINUE
c     ... Convert STO to GTO using 1 Gaussian (STO1G)
      GNU = D2*SEXP*SEXP*STO1G(N)
C     ... Make sure we dont get more delocalized than the first 
C         1s electron.
C      GNU = MAX(GNU,0.50D0)
10    CONTINUE
      CALL QEXIT ('BNCOR1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck oldbncr */
      SUBROUTINE OLDBNCR(FMAT,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     OLDBNCR: Generate screened bare nucleus matrix
C     Aug 2001 Hans Joergen Aa. Jensen (based on TSCORR)
C     Last revision: Sep 17 2002 - jkp (added timings, adapted for CATDR1)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0)
C
      LOGICAL DOINT(2,2)
      CHARACTER*4 OMITVNUC
      DIMENSION WORK(LWORK),FMAT(*)
      DIMENSION GNUTAB(103)
      SAVE GNUTAB
      DATA GNUTAB /2*0.1D0, 8*0.7D0, 8*1.0D0, 18*1.1D0, 18*1.2D0,
     &             49*2.0D0/
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "nuclei.h"
#include "symmet.h"
#include "cbiher.h"
#include "cbihr1.h"
#include "cbisol.h"
C
      CALL QENTER('OLDBNCR')
c
c     Start timer
c
      TIMSTR = SECOND() 
c
#include "memint.h"
C     Memory allocation
      CALL MEMGET('REAL',KONEMT, N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN, NNBBASX,WORK,KFREE,LFREE)
      LFCBA = NUCDEP**3
      CALL MEMGET('REAL',KFCBA ,  LFCBA ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCOORC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSIGNC,3*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KNCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJSYMC,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KJCENT,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KGEXP ,  NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KXGEXP,  NUCIND,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KCHRG ,  NUCIND,WORK,KFREE,LFREE)
      LWRK   = LWORK - KFREE + 1
      CALL DZERO(WORK(KCHRG),NUCIND)
c
c     ... set charges and gaussian exponents.
c     
C       TODO: check for point charges ...
C       for now we correct them as well
c
      IF(IPRINT .GE. 2) THEN
         WRITE(LUPRI,'(/,2(4X,A/))') 
     &' *** Entering OLDBNCR  for <X_A|(-Z*G(a,r_C))/r_C|X_B> terms ***',
     &'     X=L,S ; G(a,r)=exp(-a*r^2)'
C
         WRITE(LUPRI,'(14X,A)')
     &   ' Nuclei      Charge(-Z)       Exponent(a)'
         WRITE(LUPRI,'(14X,A)')
     &   '-----------------------------------------'
         DO I = 1,NUCIND
           ICHRG = NINT(CHARGE(I))
           ICHRG = ABS(ICHRG)
           IF (ICHRG .GT. 0) THEN
              WORK(KXGEXP-1+I) =  GNUTAB(ICHRG)
              WORK(KCHRG-1+I)  = -CHARGE(I)
           END IF
           WRITE(LUPRI,'(15X,A5,F15.2,F15.2)')
     &     NAMN(I),WORK(KCHRG-1+I),WORK(KXGEXP-1+I)
         ENDDO
         WRITE(LUPRI,'(14X,A/)')
     &   '-----------------------------------------'
      ENDIF
C
C     Initialization (calculate LL and SS)
C
      DOINT(1,1) = .TRUE.
      IF (NOSMLV) THEN
         DOINT(2,2) = .FALSE.
      ELSE
         DOINT(2,2) = .TRUE.
      END IF
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
C
C     Omit selected terms depending on OMITVNUC (See CATNUC for details)
C     In this version one-center terms are calculated with the same
C     charges/exponents as multi-center terms.
C
      OMITVNUC = 'FFFF'
c
c     ... Make list of nuclei.
c
      CALL CATNUC(NUCDEP,NCENTC,WORK(KCHRG),WORK(KFCBA),WORK(KGEXP),
     &            WORK(KXGEXP),WORK(KCOORC), WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),OMITVNUC,IPRINT)
c
c    ... Calculate integrals (both LL and SS)
c
      CALL DZERO(WORK(KONEMT),N2BBASX)
      CALL CATDR1(WORK(KONEMT),DUMMY,NCENTC,WORK(KFCBA),
     &            WORK(KCOORC),WORK(KSIGNC),WORK(KJSYMC),
     &            WORK(KJCENT),WORK(KGEXP),WORK(KFREE),LWRK,IPRINT,
     &            .FALSE.,0,NNBBASX,DOINT)
C
C     Symmetry unpack integrals
C
      CALL SYMUPK(WORK(KONEMT),WORK(KONEIN),1,NNBBASX)
      CALL DSPTSI(NTBAS(0),WORK(KONEIN),WORK(KONEMT))
      CALL BUTOBS(WORK(KONEMT),1,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.5) THEN
         CALL HEADER('OLDBNCR correction FMAT',-1)
         CALL OUTPUT(WORK(KONEMT),1,NTBAS(0),1,NTBAS(0),
     &        NTBAS(0),NTBAS(0),1,LUPRI)
      END IF
C
C     Add correction to Fock-matrix
C
      CALL DAXPY(N2BBASX,D1,WORK(KONEMT),1,FMAT,1)
      IF(IPRINT.GE.5) THEN
         CALL HEADER('OLDBNCR: screened corrected bare nucleus FMAT',-1)
         CALL PRQMAT(FMAT,NTBAS(0),NTBAS(0),
     &               NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
      END IF
C
C     Memory deallocation
      CALL MEMREL('OLDBNCR',WORK,KWORK,KWORK,KFREE,LFREE)
C
C
c     Stop timer
c
      TIMEND = SECOND()
      TIME   = TIMEND - TIMSTR
      IF (IPRINT.GE.2)
     &  CALL TIMTXT('>>> Time used in OLDBNCR is',TIME,LUPRI)
      CALL QEXIT('OLDBNCR')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck eonecap */
      SUBROUTINE EONECAP(FMAT,DMAT,WORK,LWORK,IQ)
C****************************************************************************
C     PURPOSE:
C        Calculate correction to the nuclear attraction energy 
C        from the one-center models
C
C      E = D(IJ)(*)F(IJ)
C
C     Sep 2001 : Jesper Kielberg Pedersen
C     Based on ERGCAL
C     Revised Juli 2002 : jkp
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D2=2.0D0)
C
#include "dcbdhf.h"
#include "dcbham.h"
#include "dcbbas.h"
      CHARACTER FMT*6,MXFORM*6
      DIMENSION FMAT(*), DMAT(*),WORK(LWORK)
C
      CALL QENTER('EONECAP')
C
      ECORR  = D0
C
C        Factor D2 comes from Time-Reversal Symmetry.
C
      IF(INTV1C.EQ.3.OR.(INTV1C.EQ.2.AND.ICTLV1C(2).EQ.0)) THEN
         IF (NOPEN .GE. 1) THEN
            CALL QUIT('EONECAP not implemented for open shell '//
     &                'for this ONECAP model yet.')
         END IF
      END IF
      ECORR = D2*DDOT(N2BBASX,DMAT,1,FMAT,1)
C
      CALL HEADER('Output from EONECAP',-1)
      FMT = MXFORM(ABS(ECORR),20)
      IF (IQ .EQ. 1) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'ONECAP correction to LL-Nuc att.  :  ',ECORR
      ELSE IF (IQ .EQ. 2) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'ONECAP correction to SS-Nuc att.  :  ',ECORR
      ELSE IF (IQ .EQ. 3) THEN
         WRITE(LUPRI,'(A,'//FMT//')')
     &  'Estimation of error on the energy in this model :  ',ECORR
      END IF
C
      CALL QEXIT('EONECAP')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck psiprj */
      SUBROUTINE PSIPRJ(TMAT,TINV,IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Eliminate a set of reference orbitals from the variational space;
C     this is done by projecting them out of the MO transformation matrix
C     Written by T.Saue March 26 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcbdhf.h"
C
      LOGICAL TOBE
      DIMENSION TMAT(*),TINV(*),WORK(*)
      DIMENSION NSTR(2,0:2,MPRJREF),KVEC(2,MPRJREF),
     &          KQ(2),KI(2),NVECS(2),NQ(2)
C
      KFRSAV = KFREE
      CALL TITLER('Reduction of variational space by projection',
     &     '*',125)
C
C     Find total number of reference orbitals
C     AND normal coefficients
      NORBT = 0
      NCMOT  = 0
      NVECT = 0
      DO IFRP = 1,NFSYM
        NORBT       = NORBT  + NORB(IFRP)
        ICMO(IFRP)  = NCMOT
        ICMOQ(IFRP) = NCMOT*NZ
        NCMO(IFRP)  = NFBAS(IFRP,0)*NORB(IFRP)
        NCMOQ(IFRP) = NCMO(IFRP)*NZ
        NCMOT       = NCMOT + NCMO(IFRP)
        NVECS(IFRP) = 0
        DO I = 1,NPRJREF
          CALL MEMGET('INTE',KVEC(IFRP,I),NORB(IFRP),WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = 1
          CALL NUMLST(VCPROJ(IFRP,I),WORK(KVEC(IFRP,I)),
     &              NORB(IFRP),-NPSH(IFRP),NESH(IFRP),
     &              IFRP,NSTR(IFRP,0,I),WORK(KFREE),LFREE)
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),IFRP,
     &                NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NVECS(IFRP) = NVECS(IFRP) + NSTR(IFRP,0,I)
          CALL MEMREL('PSIPRJ',WORK,KVEC(IFRP,I),KVEC(IFRP,I),
     &                KFREE,LFREE)
          CALL MEMGET('INTE',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
        NVECT = NVECT + NVECS(IFRP)
      ENDDO
      NCMOTQ  = NCMOT*NZ
      IF(NVECT.EQ.0) GOTO 30
C
C     Calculate dimensions of reference coefficients array
C
      NDIM  = 0
      DO IFRP = 1,NFSYM
        NQ(IFRP)         = NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        NDIM             = NDIM + NQ(IFRP)
        CALL MEMGET('INTE',KI(IFRP),NVECS(IFRP),WORK,KFREE,LFREE)
        CALL IZERO(WORK(KI(IFRP)),NVECS(IFRP))
      ENDDO       
      CALL MEMGET('REAL',KQ(1),NDIM,WORK,KFREE,LFREE)
      KQ(2)   = KQ(1) + NQ(1)
C
C
C     Select reference orbitals
C
      IOPT = 3
      CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
      IF(SPINFR.OR.LEVYLE) THEN
        IOPT = 11
        CALL MEMGET('REAL',KIBE,NORBT,WORK,KFREE,LFREE)
      ELSE
        KIBE = KFREE
      ENDIF
      DO IFRP = 1,NFSYM
        NQ(IFRP) = KQ(IFRP)
      ENDDO
      DO 10 I = 1,NPRJREF
        NVECT = 0
        DO IFRP = 1,NFSYM
          NVECT = NVECT + NSTR(IFRP,0,I)
        ENDDO
        IF(NVECT.EQ.0) THEN
          WRITE(LUPRI,'(A,A6)') 
     &       '* No coefficients selected for reference ',PRJFIL(I)
          GOTO 10
        ENDIF
        INQUIRE(FILE=PRJFIL(I),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
           WRITE(LUPRI,'(A,A6,A)') 
     &     'Reference coefficient file ',PRJFIL(I),' not present'
           CALL QUIT('PSIPRJ: No reference coefficients !')
        ENDIF
        CALL OPNFIL(LUCOEF,PRJFIL(I),'OLD','PSIPRJ')
        WRITE(LUPRI,'(A,A6)') 
     &       '* Coefficients for reference ',PRJFIL(I)
        IF(PROOWN) THEN
          CALL SELOWN(IOPT,WORK(KCMO),DUM,WORK(KIBE),
     &                I,LUCOEF,PRJFIL(I),NQ,IDUM,KI,
     &                NVECS,KVEC(1,I),NSTR(1,0,I),WORK,KFREE,LFREE)
          DO IFRP = 1, NFSYM
            IF(NSTR(IFRP,0,I).GT.0) THEN
              NQ(IFRP) = NQ(IFRP) + NSTR(IFRP,0,I)*NFBAS(IFRP,0)
            ENDIF
         ENDDO
        ELSE
          CALL REACMO(LUCOEF,PRJFIL(I),WORK(KCMO),DUM,WORK(KIBE),TOTERG,
     &                IOPT,WORK,KFREE,LFREE)
          DO IFRP = 1, NFSYM
            IF(NSTR(IFRP,0,I).GT.0) THEN
              CALL SELCFS (WORK(KCMO+ICMOQ(IFRP)),IFRP,WORK(NQ(IFRP)),
     &                       NVECS(IFRP),WORK(KVEC(IFRP,I)),
     &                       NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                       NFBAS(IFRP,0),NORB(IFRP))
              NQ(IFRP) = NQ(IFRP) + NSTR(IFRP,0,I)*NFBAS(IFRP,0)
              IF(SPINFR.OR.LEVYLE) THEN
C               Spin-free case: Fill array with boson irrep info
                CALL SELIBEIG(WORK(KIBE),IORB(IFRP),IFRP,
     &                         WORK(KI(IFRP)),WORK(KVEC(IFRP,I)),
     &                         NSTR(IFRP,2,I),NSTR(IFRP,1,I))
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        CLOSE(LUCOEF,STATUS='KEEP')
 10   CONTINUE
C
C     Throw away the full set of coefficients and eigenvalues
C     =======================================================
C
      CALL MEMREL('PSIPRJ.cfs',WORK,KCMO,KCMO,KFREE,LFREE)
C
C     If frozen orbitals, save to file
C     ================================
C
      IF(NFROT.GE.0) THEN
        CALL OPNFIL(LUCOEF,'DFFROZ','UNKNOWN','PSIPRJ')
        DO IFRP = 1,NFSYM
        IF(NFRO(IFRP).GT.0) THEN
          NBUF = NFBAS(IFRP,0)*NZ
          CALL MEMGET('INTE',KSEL,NVECS(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET('INTE',KIND,NFRO(IFRP),WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KBUF,NBUF,WORK,KFREE,LFREE)
          CALL WRIFRZ(LUCOEF,IFRP,WORK(KQ(IFRP)),WORK(KSEL),
     &                WORK(KBUF),WORK(KIND),NVECS(IFRP))
          CALL MEMREL('PSIPRJ.froz',WORK,KSEL,KSEL,KFREE,LFREE)       
        ENDIF
        ENDDO
        CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF
C
C     Reduce variational space by projection using the reference orbitals
C     ===================================================================
C     If more than one set of reference orbitals, they must be orthogonalized !
C
C     Do projection
      TOBE = .FALSE.
      N2TMT = 0
      DO 20 IFRP = 1,NFSYM
        WRITE(LUPRI,'(A,A3)') '* Fermion ircop ',FREP(IFRP)        
C       Transform reference coefficients to orthonormal basis
        NDEL = 0
        IF(NFRO(IFRP).EQ.0) GOTO 21
        NREF = NORB(IFRP)*NVECS(IFRP)*NZ
        NPRJ = NORB(IFRP)*NORB(IFRP)*NZ
        NBUF = MAX(NREF,NPRJ)
        CALL MEMGET('REAL',KBUF,NBUF,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBOS,NORB(IFRP),WORK,KFREE,LFREE)
        CALL QGEMM(NORB(IFRP),NVECS(IFRP),NFBAS(IFRP,0),D1,
     &            'H','N',IPQTOQ(1,0),TINV(I2TMT(IFRP)+1),
     &                    NFBAS(IFRP,0),NORB(IFRP),NZ,
     &            'N','N',IPQTOQ(1,0),WORK(KQ(IFRP)),
     &                    NFBAS(IFRP,0),NVECS(IFRP),NZ,
     &            D0,IPQTOQ(1,0),WORK(KBUF),NORB(IFRP),NVECS(IFRP),NZ)
        CALL DCOPY(NREF,WORK(KBUF),1,WORK(KQ(IFRP)),1)
C       Initialize unit matrix (transformation matrix in MO-basis)
        CALL DUNIT2(WORK(KBUF),NORB(IFRP),NORB(IFRP),NORB(IFRP),NZ)
        IF(SPINFR.OR.LEVYLE) THEN
C         Spin-free case: Fill array with boson irrep info
          CALL INIBOS(WORK(KIBOS),IFRP)
        ELSE
C         Spin-orbit coupling: Fill array with zeros..
          CALL ICOPY (NORB(IFRP),0,0,WORK(KIBOS),1) 
        ENDIF
C       Do projection
        CALL PRJORT(WORK(KBUF),WORK(KIBOS),
     &              WORK(KQ(IFRP)),WORK(KI(IFRP)),
     &              PRJTHR,THRMAX,NORB(IFRP),NORB(IFRP),
     &              NVECS(IFRP),NDEL,IPRHAM)
 21     CONTINUE
        WRITE(LUPRI,'(2X,A,I6)') '- Number of orbitals deleted: ',NDEL
        NTMO(IFRP) = NORB(IFRP)-NDEL
        IEFF       = I2TMT(IFRP)+1
        I2TMT(IFRP) = N2TMT
        N2TMT = N2TMT + NFBAS(IFRP,0)*NTMO(IFRP)*NZT
        IF(NDEL.GT.0) THEN
          IF(NFRO(IFRP).GT.0) THEN
            IF(NDEL.LT.NFRO(IFRP)) THEN
              WRITE(LUPRI,'(2X,A,I6,A)') 
     &        'is smaller than number of frozen orbitals',
     &        NFRO(IFRP),'!'
              WRITE(LUPRI,'(2X,A)')
     &        'You may rerun with PRINT.GE.2 under **HAMILTONIAN',
     &        'to see the norms of the remaining vectors.',
     &        'Based on this, you may set the projection threshold'
              WRITE(LUPRI,'(2X,A,E12.4)')
     &        'PRJTHR lower. Current value is ',PRJTHR
              CALL QUIT
     &        ('PSIPRJ:Fewer deleted orbitals than frozen ones!')
            ENDIF
          ENDIF
C         Backtransform to SO-basis
          NDIM = NFBAS(IFRP,0)*NORB(IFRP)*NZT
          CALL DCOPY(NDIM,TMAT(IEFF),1,TINV(IEFF),1)
          CALL PRJCNT(WORK(KBUF),TMAT(I2TMT(IFRP)+1),
     &                IFRP,NORB(IFRP),NFBAS(IFRP,0),
     &                NTMO(IFRP),WORK(KIBOS),TINV(IEFF))
CTROND: It is assumed that the orbitals projected out belong to the electronic
CTROND: space. This must be worked on !!!!
          NESH(IFRP) = NESH(IFRP)-NDEL+NFRO(IFRP)
          TOBE = .TRUE.       
        ELSEIF((I2TMT(IFRP)+1).NE.IEFF) THEN
C         Move transformation matrix if columns was deleted in previous
C         irrep ...
          NDIM = NFBAS(IFRP,0)*NTMO(IFRP)*NZT
          CALL DCOPY(NDIM,TMAT(IEFF),1,TMAT(I2TMT(IFRP)+1),1)
        ENDIF
        IF(NFRO(IFRP).GT.0) THEN
          CALL MEMREL('PSIPRJ.prj',WORK,KBUF,KBUF,KFREE,LFREE)
        ENDIF
 20   CONTINUE
      IF(IPRHAM.GE.5) THEN
        CALL HEADER('PSIPRJ: Final transformation matrix ',-1)
        DO 80 I = 1,NFSYM
        IF(NTMO(I).GT.0) THEN
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TMAT(I2TMT(I)+1),NFBAS(I,0),NTMO(I),
     &                NFBAS(I,0),NTMO(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
 80     CONTINUE
      ENDIF
      IF(TOBE) CALL MKSAOMO(TINV,TMAT,IPRINT,WORK(KFREE),LFREE)
 30   CONTINUE
      CALL MEMREL('PSIPRJ',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRJOUT */
      SUBROUTINE PRJOUT(VECS,VREF,THRLDP,THRMAX,SMAT,VBUF,IVECS,
     &                  NVDIM,NVECS,NVREF,NDEL,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1=-1.0D0)
      PARAMETER (T1MIN = 1.0D-8, THRRND = 1.0D-4)
C
#include "dgroup.h"
      DIMENSION VECS(NVDIM,NVECS,NZ),VREF(NVDIM,NVREF,NZ),
     &          VBUF(NVDIM,NZ),SMAT(NVDIM,NVDIM),IVECS(NVECS)
      DIMENSION P(4)
C
      THRMAX = D0
      NDEL   = 0
      N2VECS = NVDIM*NVECS
      N2VREF = NVDIM*NVREF
      IROUND = 0
      ITURN  = 0
 10   CONTINUE
      ITURN  = ITURN + 1
C
C     First project out reference vectors
C     ===================================
C
      DO IREF = 1,NVREF
        CALL QGEMM(NVDIM,1,NVDIM,D1,
     &         'N','N',IPQTOQ(1,0),SMAT,NVDIM,NVDIM,1,
     &         'N','N',IPQTOQ(1,0),VREF(1,IREF,1),NVDIM,NVREF,NZ,
     &              D0,IPQTOQ(1,0),VBUF,NVDIM,1,NZ)
        IF(IPRINT.GE.2) THEN
          CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VBUF,NVDIM,NZ,1,
     &           'N',IPQTOQ(1,0),VREF(1,IREF,1),N2VREF,NZ,1)
          VNRM=SQRT(VNRM)
          WRITE(LUPRI,'(A,I3,A,1P,D12.5)')
     &      'Reference vector no.',IREF,' Norm: ',VNRM
        ENDIF
        DO 100 IVEC = 1,NVECS
          IF(IVECS(IVEC).LT.0) GOTO 100
C
C         Form -<VREF|IVEC>
C
          CALL QDOT(P,IPQTOQ(1,0),NZ,NVDIM,
     &         'H',IPQTOQ(1,0),VBUF,NVDIM,NZ,1,
     &         'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
          DO IZ = 1,NZ
            P(IZ) = -P(IZ)
          ENDDO
C
C         Form VECS = VECS -VREF<VREF|VECS>
C
          CALL  QXAPY(NVDIM,'N','N',P,IPQTOQ(1,0),NZ,
     &                VREF(1,IREF,1),N2VREF,1,IPQTOQ(1,0),NZ,
     &                VECS(1,IVEC,1),N2VECS,1,IPQTOQ(1,0),NZ)
 100    CONTINUE
      ENDDO
C
C     Next commence internal orthogonalization
C     =========================================
C
      DO 200 IVEC = 1,NVECS
        IF(IVECS(IVEC).LT.0) GOTO 200
        CALL QGEMM(NVDIM,1,NVDIM,D1,
     &         'N','N',IPQTOQ(1,0),SMAT,NVDIM,NVDIM,1,
     &         'N','N',IPQTOQ(1,0),VECS(1,IVEC,1),NVDIM,NVECS,NZ,
     &              D0,IPQTOQ(1,0),VBUF,NVDIM,1,NZ)
C
C       Check size of norm
C       ==================
C
        CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VBUF,NVDIM,NZ,1,
     &           'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
        VNRM=SQRT(VNRM)
        THRMAX=MAX(THRMAX,VNRM)
        IF(VNRM.LT.THRLDP) THEN
          WRITE(LUPRI,'(A,I3,A,1P,D12.5)')
     &      'Vector no.',IVEC,' removed. Linear dependence: ',VNRM
          IVECS(IVEC) = -1
          NDEL = NDEL + 1
          GOTO 200
        ELSE
          IF(IPRINT.GE.2) WRITE(LUPRI,'(A,I3,A,1P,D12.5)')
     &      'Vector no.',IVEC,' Norm: ',VNRM
          IF(VNRM.LT.THRRND.AND.ITURN.EQ.1) IROUND=IROUND+1
          IF(VNRM.LT.T1MIN) THEN
            VNRM = D1/VNRM
            DO IZ = 1,NZ
              CALL DSCAL(NVDIM,VNRM,VECS(1,IVEC,IZ),1)
            ENDDO
            CALL QGEMM(NVDIM,1,NVDIM,D1,
     &         'N','N',IPQTOQ(1,0),SMAT,NVDIM,NVDIM,1,
     &         'N','N',IPQTOQ(1,0),VECS(1,IVEC,1),NVDIM,NVECS,NZ,
     &              D0,IPQTOQ(1,0),VBUF,NVDIM,1,NZ)
            CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VBUF,NVDIM,NZ,1,
     &           'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
          END IF
          VNRM = D1 / VNRM
          DO IZ = 1,NZ
            CALL DSCAL(NVDIM,VNRM,VECS(1,IVEC,IZ),1)
          ENDDO
        ENDIF
        CALL QGEMM(NVDIM,1,NVDIM,D1,
     &         'N','N',IPQTOQ(1,0),SMAT,NVDIM,NVDIM,1,
     &         'N','N',IPQTOQ(1,0),VECS(1,IVEC,1),NVDIM,NVECS,NZ,
     &              D0,IPQTOQ(1,0),VBUF,NVDIM,1,NZ)
        DO 300 JVEC = (IVEC+1),NVECS
          IF(IVECS(JVEC).LT.0) GOTO 300
C
C         Form -<IVEC|JVEC>
C
          CALL QDOT(P,IPQTOQ(1,0),NZ,NVDIM,
     &         'H',IPQTOQ(1,0),VBUF,NVDIM,NZ,1,
     &         'N',IPQTOQ(1,0),VECS(1,JVEC,1),N2VECS,NZ,1)
          DO IZ = 1,NZ
            P(IZ) = -P(IZ)
          ENDDO
C
C         Form JVEC = JVEC -IVEC<IVEC|JVEC>
C
          CALL  QXAPY(NVDIM,'N','N',P,IPQTOQ(1,0),NZ,
     &                VECS(1,IVEC,1),N2VECS,1,IPQTOQ(1,0),NZ,
     &                VECS(1,JVEC,1),N2VECS,1,IPQTOQ(1,0),NZ)
 300    CONTINUE
 200  CONTINUE
C
C     Check whether a new round of orthonormalization is
C     necessary
C     ==================================================
C
CTROND      IF(IROUND.GT.0) THEN
CTROND         IROUND = 0
         IF (ITURN.EQ.1) THEN
CTROND           WRITE(LUPRI,'(A)') 
CTROND     &     '* PRJORT: Second round of orthonormalization !'
           GO TO 10
         ENDIF
CTROND         CALL QUIT(
CTROND     &     'PRJOUT error: IROUND .gt. 0 and ITURN .gt. 1')
CTROND      END IF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PRJORT */
      SUBROUTINE PRJORT(VECS,IBVEC,VREF,IBREF,THRLDP,THRMAX,
     &                  NVDIM,NVECS,NVREF,NDEL,IPRINT)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,DM1=-1.0D0)
      PARAMETER (T1MIN = 1.0D-8, THRRND = 1.0D-4)
C
#include "dgroup.h"
      DIMENSION VECS(NVDIM,NVECS,NZ),VREF(NVDIM,NVREF,NZ),
     &          IBVEC(NVECS),IBREF(NVREF)
      DIMENSION P(4)
C
      IF(IPRINT.GE.7) THEN
         CALL HEADER('PRJORT: MO transformation matrix on input',-1)
         CALL PRQMAT(VECS,NVDIM,NVECS,NVDIM,NVECS,1,IPQTOQ(1,0),LUPRI)
      ENDIF
      THRMAX = D0
      NDEL   = 0
      N2VECS = NVDIM*NVECS
      N2VREF = NVDIM*NVREF
      IROUND = 0
      ITURN  = 0
 10   CONTINUE
      ITURN  = ITURN + 1
C
C     First project out reference vectors
C     ===================================
C
      DO IREF = 1,NVREF
        IF(IPRINT.GE.2) THEN
          CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VREF(1,IREF,1),N2VREF,NZ,1,
     &           'N',IPQTOQ(1,0),VREF(1,IREF,1),N2VREF,NZ,1)
          VNRM=SQRT(VNRM)
          WRITE(LUPRI,'(A,I3,A,1P,D12.5)')
     &      'Reference vector no.',IREF,' Norm: ',VNRM
        ENDIF
        DO 100 IVEC = 1,NVECS
          IF(IBVEC(IVEC).NE.IBREF(IREF)) GOTO 100
C
C         Form -<VREF|IVEC>
C
          CALL QDOT(P,IPQTOQ(1,0),NZ,NVDIM,
     &         'H',IPQTOQ(1,0),VREF(1,IREF,1),N2VREF,NZ,1,
     &         'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
          DO IZ = 1,NZ
            P(IZ) = -P(IZ)
          ENDDO
C
C         Form VECS = VECS -VREF<VREF|VECS>
C
          CALL  QXAPY(NVDIM,'N','N',P,IPQTOQ(1,0),NZ,
     &                VREF(1,IREF,1),N2VREF,1,IPQTOQ(1,0),NZ,
     &                VECS(1,IVEC,1),N2VECS,1,IPQTOQ(1,0),NZ)
 100    CONTINUE
      ENDDO
C
C     Next commence internal orthogonalization
C     =========================================
C
      DO 200 IVEC = 1,NVECS
        IF(IBVEC(IVEC).LT.0) GOTO 200
C
C       Check size of norm
C       ==================
C
        CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1,
     &           'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
        VNRM=SQRT(VNRM)
        THRMAX=MAX(THRMAX,VNRM)
        IF(VNRM.LT.THRLDP) THEN
          WRITE(LUPRI,'(3X,A,I3,A,1P,D12.5)')
     &      'Vector no.',IVEC,' removed. Linear dependence: ',VNRM
          IBVEC(IVEC) = -1-IBVEC(IVEC)
          NDEL = NDEL + 1
          GOTO 200
        ELSE
          IF(IPRINT.GE.2) WRITE(LUPRI,'(A,I3,A,1P,D12.5)')
     &      'Vector no.',IVEC,' Norm: ',VNRM
          IF(VNRM.LT.THRRND.AND.ITURN.EQ.1) IROUND=IROUND+1
          IF(VNRM.LT.T1MIN) THEN
            VNRM = D1/VNRM
            DO IZ = 1,NZ
              CALL DSCAL(NVDIM,VNRM,VECS(1,IVEC,IZ),1)
            ENDDO
            CALL QDOT(VNRM,IPQTOQ(1,0),1,NVDIM,
     &           'H',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1,
     &           'N',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1)
          END IF
          VNRM = D1 / VNRM
          DO IZ = 1,NZ
            CALL DSCAL(NVDIM,VNRM,VECS(1,IVEC,IZ),1)
          ENDDO
        ENDIF
        DO 300 JVEC = (IVEC+1),NVECS
          IF(IBVEC(JVEC).NE.IBVEC(IVEC)) GOTO 300
C
C         Form -<IVEC|JVEC>
C
          CALL QDOT(P,IPQTOQ(1,0),NZ,NVDIM,
     &         'H',IPQTOQ(1,0),VECS(1,IVEC,1),N2VECS,NZ,1,
     &         'N',IPQTOQ(1,0),VECS(1,JVEC,1),N2VECS,NZ,1)
          DO IZ = 1,NZ
            P(IZ) = -P(IZ)
          ENDDO
C
C         Form JVEC = JVEC -IVEC<IVEC|JVEC>
C
          CALL  QXAPY(NVDIM,'N','N',P,IPQTOQ(1,0),NZ,
     &                VECS(1,IVEC,1),N2VECS,1,IPQTOQ(1,0),NZ,
     &                VECS(1,JVEC,1),N2VECS,1,IPQTOQ(1,0),NZ)
 300    CONTINUE
 200  CONTINUE
      IF (ITURN.LT.3) THEN
        GO TO 10
      ENDIF
      IF(IPRINT.GE.7) THEN
         CALL HEADER('PRJORT: MO transformation matrix on output',-1)
         CALL PRQMAT(VECS,NVDIM,NVECS,NVDIM,NVECS,1,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck prjcnt */
      SUBROUTINE PRJCNT(VMO,VAO,IFRP,NMO,NAO,NEFF,IVEC,TMAT)
C***********************************************************************
C
C     After a projection step, backtransform to SO-basis
C     VMO is Mo transformation matrix after transformation.
C     The columns for which IVEC=1 are backtransformed to
C     VAO, the transformation matrix in AO-basis, using
C     the transformation matrix TMAT
C
C     Written by T.Saue Jan 31 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION VMO(NMO,NMO,NZ),VAO(NAO,NEFF,NZ),IVEC(*),
     &          TMAT(NAO,NMO,NZ)
      DIMENSION IBOS(0:7)
      NBRP = 4/NZ
      DO ISYM = 1, NBRP
        IREP = JFSYM(ISYM,IFRP) - 1
        IBOS(IREP) = ISYM
      ENDDO
      II = 0
      DO I = 1,NMO
      IF(IVEC(I).GE.0) THEN
        II = II + 1
        CALL QGEMM(NAO,1,NMO,D1,
     &          'N','N',IPQTOQ(1,0),TMAT,NAO,NMO,NZ,
     &          'N','N',IPQTOQ(1,0),VMO(1,I,1), NMO,NMO,NZ,
     &          D0,IPQTOQ(1,0),VAO(1,II,1),NAO,NEFF,NZ)
      ELSE
        IF(SPINFR.OR.LEVYLE) THEN
C         tsaue: I assume that electronic orbitals are projected out...
          IREP = -1-IVEC(I)
          ISYM = IBOS(IREP)
          NBORB(ISYM,IFRP,1) = NBORB(ISYM,IFRP,1) - 1
        ENDIF
      ENDIF
      ENDDO
      DO ISYM = 1, NBRP
        NBORB(ISYM,IFRP,0) = NBORB(ISYM,IFRP,1) + NBORB(ISYM,IFRP,2)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck prjcnt2 */
      SUBROUTINE PRJCNT2(IVEC,IFRP,NBAS,TMAT,NMAT,TBUF,NBUF)
C***********************************************************************
C
C     After a projection step, update dimension arrays in 
C     COMMON blocks DCBBAS and DCBORB
C     Written by T.Saue Jan 31 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
      DIMENSION IVEC(*),TMAT(NBAS,NMAT,NZ),TBUF(NBAS,NBUF,NZ)
      NSYM = 4/NZ
      JVEC = 0
      JBUF = 0
      DO ISYM = 1,NSYM
        DO IC = 2,1,-1
          NVEC = NBORB(ISYM,IFRP,IC)
          IF(NVEC.GT.0) THEN
            NBORB(ISYM,IFRP,IC)=ISUM(NVEC,IVEC(JVEC+1),1)
            DO II = 1,NVEC
              JVEC = JVEC + 1
              IF(IVEC(JVEC).EQ.1) THEN
                JBUF = JBUF + 1
                DO IZ = 1,NZ
                  CALL DCOPY(NBAS,TMAT(1,JVEC,IZ),1,TBUF(1,JBUF,IZ),1)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        NBORB(ISYM,IFRP,0) = NBORB(ISYM,IFRP,1) + NBORB(ISYM,IFRP,2)
      ENDDO
      NESH(IFRP) = ISUM(NSYM,NBORB(1,IFRP,1),1)
      NPSH(IFRP) = ISUM(NSYM,NBORB(1,IFRP,2),1)
      NORB(IFRP) = NESH(IFRP) + NPSH(IFRP)
      RETURN
      END
C
C*****************************************************************************
C  /* Deck mk1cex */
      SUBROUTINE MK1CEX (WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construct one-center expansion coefficients
C        
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION WORK(*)
C
      CALL QENTER('MK1CEX')
      KSAVE = KFREE
C
      LDM = NTBAS(0)
      LDM2 = LDM * LDM
      CALL MEMGET('REAL',KSMAT1,LDM2,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT2,LDM2,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAO,LDM,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KINDAOS,LDM,WORK,KFREE,LFREE)
      NIND = 2 * NUCDEP + 2
      CALL MEMGET('INTE',KINDNAO,NIND,WORK,KFREE,LFREE)
C
C     Get overlap matrix
C     ==================
C
      CALL GTOVLX(WORK(KSMAT1),D1,WORK(KFREE),LFREE)
C
C     Backtransform overlap matrix to AO-basis
C     ========================================
C
      CALL MTBSBU(WORK(KSMAT1),WORK(KSMAT2))
      CALL MTSOAO(WORK(KSMAT2),WORK(KSMAT1),NTBAS(0),0,IPRINT)
C
C     Get index arrays to sort on expansion centers (nuclei)
C     ======================================================
C
      CALL INDAOC (WORK(KINDAO),WORK(KINDAOS),WORK(KINDNAO))
C
C     Sort the overlap matrix to this order
C     =====================================
C
      CALL REORMT (WORK(KSMAT1),WORK(KSMAT2),LDM,LDM,
     &             WORK(KINDAO),WORK(KINDAO))
C
C     Construct the projection operator and put it in SMAT1
C     =====================================================
C
      CALL MK1CEX2 (LDM,WORK(KSMAT2),WORK(KSMAT1),WORK(KINDNAO),
     &              WORK,KFREE,LFREE,IPRINT)
C
C     Sort the projection operator back to the standard AO order
C     ========================================================
C
      CALL REORMT (WORK(KSMAT1),WORK(KSMAT2),LDM,LDM,
     &             WORK(KINDAOS),WORK(KINDAOS))
C
C     Write the projection operator to file for later use
C     ===================================================
C
      CALL OPNFIL(LUPMAT,'DFPMAT','UNKNOWN','MK1CEX')
      CALL WRITT(LUPMAT,N2BBASX,WORK(KSMAT2))
      CLOSE(LUPMAT,STATUS='KEEP')
C
      CALL MEMREL('MK1CEX',WORK,KSAVE,KSAVE,KFREE,LFREE)
C
      CALL QEXIT('MK1CEX')
      RETURN
      END
C*****************************************************************************
C  /* Deck mk1cex2 */
      SUBROUTINE MK1CEX2 (LDM,STOT,PTOT,INDNAO,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     PURPOSE:
C        Construct one-center expansion coefficients
C        
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "pgroup.h"
#include "dcbbas.h"
#include "mxcent.h"
#include "nuclei.h"
      DIMENSION STOT(LDM,LDM),PTOT(LDM,LDM)
      DIMENSION WORK(*)
      DIMENSION INDNAO(2,*)
C
      CALL QENTER('MK1CEX2')
      KSAVE = KFREE
      CALL DZERO (PTOT,LDM*LDM)
C
C     Loop over large/small component
C     Only the small part is considered because the
C     one-center approximation is not used in the 
C     large component.
C 
      DO IC = 2, 2
         NB0 = INDNAO(IC,1)+1
C
C        Loop over expansion centers
C
         DO IEXPC = 1, NUCDEP
           NB1 = INDNAO(IC,IEXPC)+1
           NB2 = INDNAO(IC,IEXPC+1)
           NLOC = NB2-NB1+1
           NV = NLOC*LDM
           NW = NLOC*NLOC
           CALL MEMGET('REAL',KEIG,NLOC,WORK,KFREE,LFREE)
           CALL MEMGET('REAL',KV,NV,WORK,KFREE,LFREE)
           CALL MEMGET('REAL',KW,NW,WORK,KFREE,LFREE)
C
C          Make the local projection operator
C
           CALL MK1CEX3 (STOT(NB1,NB0),STOT(NB1,NB1),PTOT(NB1,NB0),LDM,
     &                   NTBAS(IC),NLOC,
     &                   WORK(KV),WORK(KW),WORK(KEIG),STOL(IC),
     &                   IPRINT,WORK,KFREE,LFREE)
           CALL MEMREL('MK1CEX2',WORK,KEIG,KEIG,KFREE,LFREE)
C
C          This operator contains also the projection on the functions of the same
C          center which should be (close to) the unit matrix. Make this part zero.
C
C          TODO : check that this part is indeed the unit matrix.
C
           DO I = NB1, NB2
              CALL DZERO(PTOT(NB1,I),NLOC)
           ENDDO
         ENDDO
      ENDDO
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Projection operator in AOC basis',-1)
        CALL OUTPUT(PTOT,1,LDM,1,LDM,LDM,LDM,1,LUPRI)
      ENDIF
C
      CALL MEMREL('MK1CEX2',WORK,KSAVE,KSAVE,KFREE,LFREE)
      CALL QEXIT('MK1CEX2')
      RETURN
      END
C*****************************************************************************
C  /* Deck mk1cex3 */
      SUBROUTINE MK1CEX3 (STOT,SLOC,PTOT,LDM,NTOT,NLOC,V,W,EIG,STOL,
     &                    IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C     PURPOSE:
C        Make projection operator for this center. 
C        We reduce the basis if the overlap matrix on this
C        center is nearly singular. This may happen as we work in
C        the untransformed non-orthogonal AO basis. STOL should be
C        consistent with the values used in LOWDIN !!
C        
C        On input :
C           STOT : The complete overlap matrix
C           SLOC : The overlap matrix for the current center
C           LDM : Leading dimension (should be NBAST)
C           NTOT : The total number of functions to be projected 
C                  (now all large or small functions in the molecule
C                   but generalization to domains is easy)
C           NLOC : The number of functions on the current center
C           STOL : Linear dependency treshold
C
C        On output :
C           PTOT : The updated projection operator matrix (added the
C                  rows for the current center)
C
C        NB : STOT and PTOT should be positioned such that the current
C             center is the first row !
C 
C        Scratch arrays : V, W, EIG
C 
C        Written by L.Visscher
C        Last revision: Aug 16 2001 - lv
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
      DIMENSION WORK(*)
      DIMENSION STOT(LDM,*),PTOT(LDM,*)
      DIMENSION SLOC(LDM,NLOC),V(LDM,NLOC),W(NLOC,NLOC),EIG(NLOC)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local overlap matrix',-1)
        CALL OUTPUT(SLOC,1,NLOC,1,NLOC,LDM,NLOC,1,LUPRI)
      ENDIF
C
C     Form (Sloc)**-1/2 * U
C
      CALL LOWGEN(SLOC,NLOC,V,NV,LDM,STOL,IPRINT,EIG,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local Lowdin matrix',-1)
        CALL OUTPUT(V,1,NLOC,1,NV,LDM,NLOC,1,LUPRI)
      ENDIF
C
C     Form (Sloc)**-1/2 * U * U+ * (Sloc)**-1/2 = (Sloc)**-1
C
      CALL DGEMM('N','T',NLOC,NLOC,NV,D1,V,LDM,V,LDM,D0,W,NLOC)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local inverse overlap matrix',-1)
        CALL OUTPUT(W,1,NLOC,1,NLOC,NLOC,NLOC,1,LUPRI)
      ENDIF
C
C     Make (Sloc)**-1  * Stot : the projection operator  
C
      CALL DGEMM('N','N',NLOC,NTOT,NLOC,D1,W,NLOC,STOT,LDM,D0,PTOT,LDM)
C
      IF(IPRINT.GE.7) THEN
        CALL HEADER('Local projection operator',-1)
        CALL OUTPUT(PTOT,1,NLOC,1,NTOT,LDM,NLOC,1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX1 */
      SUBROUTINE PR1CEX1 (FTRI,FIRST,FOCK,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     Correct for neglect of multicenter integrals in symmetry packed SO matrix
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 16, 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbqal.h"
#include "dcbgen.h"
#include "dcbham.h"
C
      LOGICAL FIRST(NZ)
      DIMENSION FOCK(N2BBASX,NZ),FTRI(NNBBASX,NZ),
     &          WORK(*)
C
C     Make full matrix
C     ================
C
      DO IZ = 1,NZ
        IF(FIRST(IZ)) THEN
          CALL DZERO(FOCK(1,IZ),N2BBASX)
        ELSE
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSPTSI(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DAPTGE(NTBAS(0),FTRI(1,IZ),FOCK(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
C     Allocate work space to store the projection operator and
C     a work array.
C
      CALL MEMGET('REAL',KPROJ,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF1,N2BBASX*4,WORK,KFREE,LFREE)
C
C     Get the projection operator
C
      CALL OPNFIL(LUPMAT,'DFPMAT','OLD','PR1CEX1')
      CALL READT(LUPMAT,N2BBASX,WORK(KPROJ))
      CLOSE(LUPMAT,STATUS='KEEP')
C
C     Add contributions from missing multicenter integrals
C     ====================================================
C
      IREP = 0
      CALL PR1CEX2 (NTBAS(0),IREP,FOCK,WORK(KPROJ),WORK(KBUF1),IPRINT)
C
C     Release the memory
C
      CALL MEMREL('PR1CEX1',WORK,KPROJ,KPROJ,KFREE,LFREE)
C
C     Pack matrix
C     ===========
C
      DO IZ = 1,NZ
        IF(.NOT.FIRST(IZ)) THEN
          IQ = IPQTOQ(IZ,0)
          IH = IHQMAT(IQ,1)
          IF(IH.EQ.1) THEN        
            CALL DSITSP(NTBAS(0),FOCK(1,IZ),FTRI(1,IZ))
          ELSEIF(IH.EQ.2) THEN
            CALL DGETAP(NTBAS(0),FOCK(1,IZ),FTRI(1,IZ))
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX2 */
      SUBROUTINE PR1CEX2 (NB,IREP,FOCK,PROJ,BUF1,IPRINT)
C***********************************************************************
C     Correct for neglect of multicenter integrals in SO matrix
C
C     A basis function |q> can be projected on the set of functions
C     |r> on center A using the projection operators P(r,q) 
C     made in MK1CEX1. 
C
C     We use this to correct the input matrix F(p,q) that contains
C     only the one-center contributions (F1) because the multicenter
C     contributions (F2) were neglected.
C 
C     F(p,q) = F1(p,q) + F2(p,q)
C
C     FR(p,q) = F1(p,r) * P(r,q)
C     FL(p,q) = P+(p,r) * F1(r,q) 
C
C     F(p,q) is approximated by F1(p,q) + 1/2 (FL(p,q) + FR(p,q)) 
C     We use the average of the left and right projected matrices to
C     preserve the permutation symmetry.
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 16, 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbqal.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, HALF=0.5D0)
      DIMENSION FOCK(NB,NB,NZ),PROJ(NB,NB),BUF1(NB,NB,4)
C
C     Transform Fock matrix to AO-basis
C     ==================================
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX2: Fock matrix in SO basis',-1)
         CALL PRQMAT(FOCK,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IPAR = JBTOF(IREP,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTSOAO(FOCK(1,1,IPQ),BUF1(1,1,IZ),NB,IREPD,IPRINT)
      ENDDO
C
      DO IZ = 1, 4
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX2: Fock matrix in AO basis',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,1,LUPRI)
         ENDIF
C
C        Construct the right projected fock matrix and store it in FOCK
C
         CALL DGEMM('N','N',NB,NB,NB,D1,BUF1(1,1,IZ),NB,PROJ,NB,
     &               D0,FOCK,NB)
C
C        Construct and add the left projected fock matrix
C
         CALL DGEMM('T','N',NB,NB,NB,D1,PROJ,NB,BUF1(1,1,IZ),NB,
     &               D1,FOCK,NB)
C
C        Add the correction to the Fock matrix
C 
         IF (IPRINT.GE.10) THEN
            CALL HEADER('PR1CEX2 : Projection operator',-1)
            CALL OUTPUT(PROJ,1,NB,1,NB,NB,NB,1,LUPRI)
            CALL HEADER('PR1CEX2: Correction to Fock matrix (*2)',-1)
            CALL OUTPUT(FOCK,1,NB,1,NB,NB,NB,1,LUPRI)
         ENDIF
C
         NB2 = NB * NB
         CALL DAXPY (NB2,HALF,FOCK,1,BUF1(1,1,IZ),1)
C
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX2: Corrected Fock matrix',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,1,LUPRI)
         ENDIF
      ENDDO
C     
C     Transform Fock matrix to SO-basis
C     =================================
C
      CALL DZERO(FOCK,NB2*NZ)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTAOSO(BUF1(1,1,IZ),FOCK(1,1,IPQ),NB,IREPD,IPRINT)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX2: Corrected Fock matrix in SO basis',-1)
         CALL PRQMAT(FOCK,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PR1CEX3 */
      SUBROUTINE PR1CEX3 (NB,IREP,DENS,PROJ,BUF1,IPRINT)
C***********************************************************************
C     Correct for neglect of multicenter integrals in SO matrix
C     Note that this routines resembles PR1CEX2 that is used for
C     Fock matrix (contributions). This one is good for densities
C     where the transpose of P needs to be used and the transformation
C     to AO basis requires a different routine.
C
C     A basis function |q> can be projected on the set of functions
C     |r> on center A using the projection operators P(r,q) 
C     made in MK1CEX1. 
C
C     We use this to correct the input matrix D(q,p) to map multicenter
C     contributions to the calculated one-center integrals
C 
C     D(q,p) = D1(q,p) + D2(q,p)
C
C     DR(q,p) = D1(q,r) * P+(r,p)
C     DL(q,p) = P(q,r) * D1(r,p) 
C
C     D(q,p) is approximated by F1(q,p) + 1/2 (DL(q,p) + DR(q,p)) 
C     We use the average of the left and right projected matrices to
C     preserve the permutation symmetry.
C
C     Written by L. Visscher, august 2001
C     Last revision: Aug 17, 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbqal.h"
C
      PARAMETER (D0=0.D0, D1=1.D0, HALF=0.5D0)
      DIMENSION DENS(NB,NB,NZ),PROJ(NB,NB),BUF1(NB,NB,4)
C
C     Transform density matrix to AO-basis
C     =====================================
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX3: Density matrix in SO basis',-1)
         CALL PRQMAT(DENS,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      IPAR = JBTOF(IREP,1)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL DTSOAO(DENS(1,1,IPQ),BUF1(1,1,IZ),NB,IREPD,IPRINT)
      ENDDO
C
      DO IZ = 1, 4
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX3: Density matrix in AO basis',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,1,LUPRI)
         ENDIF
C
C        Construct the right projected fock matrix and store it in DENS
C
         CALL DGEMM('N','T',NB,NB,NB,D1,BUF1(1,1,IZ),NB,PROJ,NB,
     &               D0,DENS,NB)
C
C        Construct and add the left projected fock matrix
C
         CALL DGEMM('N','N',NB,NB,NB,D1,PROJ,NB,BUF1(1,1,IZ),NB,
     &               D1,DENS,NB)
C
C        Add the correction to the density matrix
C 
         IF (IPRINT.GE.10) THEN
            CALL HEADER('PR1CEX3 : Projection operator',-1)
            CALL OUTPUT(PROJ,1,NB,1,NB,NB,NB,1,LUPRI)
            CALL HEADER('PR1CEX3: Correction to density matrix (*2)',-1)
            CALL OUTPUT(DENS,1,NB,1,NB,NB,NB,1,LUPRI)
         ENDIF
C
         NB2 = NB * NB
         CALL DAXPY (NB2,HALF,DENS,1,BUF1(1,1,IZ),1)
C
         IF (IPRINT.GE.8) THEN
            CALL HEADER('PR1CEX3: Corrected density matrix',-1)
            CALL OUTPUT(BUF1(1,1,IZ),1,NB,1,NB,NB,NB,1,LUPRI)
         ENDIF
      ENDDO
C     
C     Transform density matrix to SO-basis
C     ====================================
C
      CALL DZERO(DENS,NB2*NZ)
      DO IZ = 1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL DTAOSO(BUF1(1,1,IZ),DENS(1,1,IPQ),NB,IREPD,IPRINT)
      ENDDO
C
      IF (IPRINT.GE.5) THEN
         CALL HEADER('PR1CEX3: Corrected density matrix in SO basis',-1)
         CALL PRQMAT(DENS,NB,NB,NB,NB,NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C
      RETURN
      END
C*****************************************************************************
C  /* Deck glinsm */
      SUBROUTINE GLINSM(WORK,LWORK)
C*****************************************************************************
C
C     This routine makes an additional transformation from the
C     matrix supplied by GMOTRA to a basis that is diagonal in
C     one component of angular momentum (j in general, l for spinfree
C     and Levy-Leblond). This can be used for systems that posess
C     linear symmetry.
C
C     Written by L. Visscher october/november 2003
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <dcbgen.h>
#include <dcbham.h>
#include <cbihr1.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION WORK(LWORK)
      CALL QENTER('GLINSM')
#include <memint.h>
C
      IF (URKBAL) THEN
         WRITE (LUPRI,1000) "unrestricted kinetic balance"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      IF (LEVYLE.OR.SPINFR) THEN
         WRITE (LUPRI,1000) "spinfree or Levy-Leblond"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
      IF (ZORA) THEN
         WRITE (LUPRI,1000) "ZORA"
         LINEAR = .FALSE.
         CALL QEXIT('GLINSM')
         RETURN
      ENDIF
C
C     Get old transformation matrix
C     =============================
C
      CALL MEMGET('REAL',KVMAT,N2TMT,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL OPNFIL(LUTMAT,'DFTMAT','OLD','PREDHF')
      CALL READT(LUTMAT,N2TMT,WORK(KVMAT))
C
      CALL LINSYM(WORK(KTMAT),WORK(KVMAT),WORK(KTBUF),WORK(KEIG),
     &            WORK(KFREE),LFREE)
C
C     Construct right index transformed overlap matrix Saomo = Sao C
C
      CALL MKSAOMO(WORK(KTBUF),WORK(KTMAT),IPRONE,WORK(KFREE),LFREE)
C
C     Write the new transformation matrix.
C
      REWIND (LUTMAT)
      CALL WRITT(LUTMAT,N2TMT,WORK(KTMAT))
      CALL WRITT(LUTMAT,N2TMT,WORK(KTBUF))
      CLOSE(LUTMAT,STATUS='KEEP')
C
C     We need to diagonalize blockwise in the SCF step 
C
      SUB_BL = .TRUE.
C
C     Memory deallocation
      CALL MEMREL('GLINSM',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('GLINSM')
      RETURN
C
 1000 FORMAT (//" WARNING : linear symmetry is not implemented for ",
     &A," yet" /"           continue calculation in lower symmetry")
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck linsym */
      SUBROUTINE LINSYM(TMAT,VMAT,TBUF,EIG,WORK,LWORK)
C*****************************************************************************
C
C     Construct matrix J(x,y,z) and transform to eigenvectors of this operator
C     Routine based on DELPOS and parts of the ZORA diagonalization routine
C
C     Written by L.Visscher, october 2003
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0=0.D0, D1=1.D0, D2=2.D0)
C
C Used from COMMON blocks
C  CBIPAM:
C  DHFCAL:
C  DCBBAS:
C  DCBORB:
C  DGROUP:
#include <dcbgen.h>
#include <dcbham.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
#include <cbihr1.h>
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      DIMENSION IQBF(4),IQW(4)
#include <memint.h>
C
      CALL HEADER('Output from LINSYM',-1)
      NBRP = 4 / NZ
C
C     Construct mj matrix
C
      CALL AMJMAT(TBUF,IPRONE,WORK,LWORK)
C
C     Transform to MO-basis
C     =====================
C
C Note : It's tempting to use I2ORBT to find the start of a block, but this is not
C possible as this array is not yet set (is done in dirset). Took me half an hour to
C figure this one out (LV).
C
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 10
        IOFTB = I2BASX(I,I) + 1
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &            NORB(I),NORB(I),
     &            TBUF(IOFTB),NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),
     &            TMAT(IOFT),NORB(I),NORB(I),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),
     &            VMAT(IOFV),NFBAS(I,0),NORB(I),1,IPQTOQ(1,0),
     &            WORK,LWORK,IPRGEN)
        IOFT = IOFT + NORB(I)*NORB(I)
        IOFV = IOFV + NORB(I)*NFBAS(I,0)
   10 CONTINUE
C
C     Initialize TBUF to zero (we will fill in parts later and its easier to initialize
C     the rest here). Also initialize the array that counts the number of functions for 
C     each Mj value.
C
      CALL DZERO (TBUF,N2TMT)
      CALL IZERO (NORB_SUB,MAX_SUB_BL*3*2)
      WRITE (LUPRI,1000)
C
      IOFT = 1
      DO 20 I = 1,NFSYM
C
        IF(NORB(I).EQ.0) GOTO 20
C         
C       Set dimensions
C
        NL = NESH(I) 
        NS = NPSH(I)
        NT = NORB(I)
        NLL = NL * NL
        NSS = NS * NS
        NTT = NT * NT
        IEIGP = IORB(I) + 1
        IEIGE = IEIGP + NS
C
C       Allocate memory for the diagonalization. We want to
C       have the LL and SS blocks in separate blocks.
C
        CALL MEMGET('REAL',KFLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KFSS,NSS,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVLL,NLL,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KVSS,NSS,WORK,KFREE,LFREE)
C
C       Extract LL block of Jz matrix
C
        CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,0,0,1,1,1,NZ,
     &               WORK(KFLL),NL,NL,1,
     &               TMAT(IOFT),NT,NT,NZ)
C
C       Find eigenvalues of Jz
C
        IF (NL.GT.0) THEN
           CALL QDIAG(1,NL,WORK(KFLL),NL,NL,
     &                EIG(IEIGE),1,
     &                WORK(KVLL),NL,NL,
     &                WORK(KFREE),LFREE,IERR)
        ENDIF
C
C       Extract SS block of Jz matrix
C
        CALL EXTRSB (1,I,NBORB,4,NBRP,NBORB,4,NBRP,0,0,2,2,1,NZ,
     &               WORK(KFSS),NS,NS,1,
     &               TMAT(IOFT),NT,NT,NZ)
C
        IF (NS.GT.0) THEN
           CALL QDIAG(1,NS,WORK(KFSS),NS,NS,
     &                EIG(IEIGP),1,
     &                WORK(KVSS),NS,NS,
     &                WORK(KFREE),LFREE,IERR)
        ENDIF
C
C       Count the Mj-values.
C       We will either find a negative or a positive value for each MJ.
C
        MJMAX = 0
        DO IX = IEIGE, IEIGE+NL-1
           MJ = NINT(D2*EIG(IX))
           MJABS = (ABS(MJ)+1)/2
           ID_SUB_BL(MJABS,I)  = MJ
           NORB_SUB(MJABS,I,1) = NORB_SUB(MJABS,I,1) + 1
           NORB_SUB(MJABS,I,0) = NORB_SUB(MJABS,I,0) + 1
           MJMAX = MAX(MJABS,MJMAX)
        ENDDO
C
        DO IX = IEIGP, IEIGP+NS-1
           MJ = NINT(D2*EIG(IX))
           MJABS = (ABS(MJ)+1)/2
           IF (MJABS.GT.MAX_SUB_BL) CALL QUIT ("increase max_sub_bl")
           ID_SUB_BL(MJABS,I)  = MJ
           NORB_SUB(MJABS,I,2) = NORB_SUB(MJABS,I,2) + 1
           NORB_SUB(MJABS,I,0) = NORB_SUB(MJABS,I,0) + 1
           MJMAX = MAX(MJABS,MJMAX)
        ENDDO
C
        N_SUB_BL(I) = MJMAX
C
C       Write information to output
C
        DO MJABS = 1, MJMAX
           WRITE (LUPRI,1010) 3-2*I,ID_SUB_BL(MJABS,I),
     &                       (NORB_SUB(MJABS,I,IC),IC=0,2)
        ENDDO
C
C       We need to be careful now. In making NORB_SUB we chose to
C       use the absolute value of MJ because this gives the most logical
C       block structure. The eigenvectors and solutions are, however,
C       still ordered by the signed value. Reorder them first before
C       continuing.
C
        DO IX = IEIGP, IEIGP+NS-1
           EIG(IX) = ABS(EIG(IX))
        ENDDO
C
C       The SELBOS routine does just the reordering that we want,
C       use work(kfss) as scratch, and allocate a few other small
C       scratch arrays. The IBEIG part will not be used but we will
C       initialize it to prevent problems with uninitialized data while
C       sorting it.
C
        CALL MEMGET('INTE',KINDX,NS,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIDUM,NS,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEI,NS,WORK,KFREE,LFREE)
        MJ = 0
        CALL ICOPY (NS,MJ,0,WORK(KIBEI),1)
        CALL SELBOS (NS,WORK(KINDX),EIG(IEIGP),
     &               WORK(KIBEI),WORK(KVSS),
     &               WORK(KFSS),WORK(KIDUM))
        CALL MEMREL('LINSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       Same thing for the positive energy solutions
C
        DO IX = IEIGE, IEIGE+NL-1
           EIG(IX) = ABS(EIG(IX))
        ENDDO
        CALL MEMGET('INTE',KINDX,NL,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIDUM,NL,WORK,KFREE,LFREE)
        CALL MEMGET('INTE',KIBEI,NL,WORK,KFREE,LFREE)
        MJ = 0
        CALL ICOPY (NL,MJ,0,WORK(KIBEI),1)
        CALL SELBOS (NL,WORK(KINDX),EIG(IEIGE),
     &               WORK(KIBEI),WORK(KVLL),
     &               WORK(KFLL),WORK(KIDUM))
        CALL MEMREL('LINSYM',WORK,KINDX,KINDX,KFREE,LFREE)
C
C       Put the MO vectors in TBUF
C
C       We now have the solutions ordered on on the right index 
C       by absolute mj-value and at the left index on boson irreps.
C       Keep this ordering when inserting into the full matrix so
C       that the new transformation matrix (generated below) will
C       give an ordering on MJ values.
C
        CALL EXTRSB (2,I,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(I),0,0,1,1,1,NZ,
     &               WORK(KVLL),NL,NL,1,
     &               TBUF(IOFT),NT,NT,NZ)
        CALL EXTRSB (2,I,NBORB,4,NBRP,
     &               NORB_SUB,MAX_SUB_BL,N_SUB_BL(I),0,0,2,2,1,NZ,
     &               WORK(KVSS),NS,NS,1,
     &               TBUF(IOFT),NT,NT,NZ)
C
C       Done; release the extra memory used for diagonalization
C
        CALL MEMREL('LINSYM',WORK,KFLL,KFLL,KFREE,LFREE)
C
C       Update offset
C
        IOFT = IOFT + NORB(I)*NORB(I)
C
   20 CONTINUE
C
C     Generate transformation matrix in non-orthogonal basis
C     ======================================================
C
      IOFT = 1
      IOFU = 1
      IOFV = 1
      DO 30 I = 1,NFSYM
        IF(NORB(I).EQ.0) GOTO 30
        IF (IPRHAM.GE.5) THEN
           CALL HEADER('Transformation matrix in orthon. basis',-1)
           CALL PRQMAT(TBUF(IOFU),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),1,IPQTOQ(1,0),LUPRI)
        ENDIF
C
C       Backtransform all solutions
C       ===========================
C
        CALL BCKTRA(TMAT(IOFT),NFBAS(I,0),NORB(I),
     &              TBUF(IOFU),NORB(I),NORB(I),
     &              NORB(I),1,
     &              NORB(I),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NORB(I),1,
     &              IPRHAM)
        IOFT = IOFT + NFBAS(I,0)*NORB(I)
        IOFU = IOFU + NORB(I)*NORB(I)
        IOFV = IOFV + NORB(I)*NFBAS(I,0)
   30 CONTINUE
C
      CALL FLSHFO(LUPRI)
      RETURN
 1000 FORMAT 
     & (/' Parity   MJ  Functions(total) Functions(LC) Functions(SC)')
 1010 FORMAT (T5,I2,T9,I3,'/2',T20,I5,T33,I5,T46,I5)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amjmat */
      SUBROUTINE AMJMAT(FOCK,IPRINT,WORK,LWORK)
C*****************************************************************************
C     Generate one-electron matrix of Jz
C
C     Written by Luuk Visscher, oct 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <cbihr1.h>
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('AMJMAT')
#include <memint.h>
C
C     Memory allocation
C
      CALL MEMGET('REAL',KFTRI ,NZ*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KONEIN,   NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KFIRST,NZ        ,WORK,KFREE,LFREE)
C
C     ************************************************
C     ****** Construct one-electron Fock matrix ******
C     ************************************************
C
      CALL AMJMA1(FOCK,WORK(KFTRI),WORK(KONEIN),WORK(KFIRST),
     &            WORK,KFREE,LFREE,IPRINT)
C
C     Memory deallocation
      CALL MEMREL('AMJMAT',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('AMJMAT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck amjma1 */
      SUBROUTINE AMJMA1(FOCK,FTRI,OP1INT,FIRST,WORK,KFREE,LFREE,IPRINT)
C*****************************************************************************
C     Generate one-electron matrix of Jz
C
C     Written by Luuk Visscher, oct 1999
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbgen.h>
#include <dcbham.h>
C
      PARAMETER (HALF=0.5D0, D1=1.0D0)
      LOGICAL FIRST
      DIMENSION FOCK(N2BBASX),FTRI(NNBBASX),OP1INT(NNBBASX),
     &          WORK(*)
C
C     Sanity check : linear symmetry requires use of real groups (C2v or D2h)
C
      IF (NZ.NE.1) CALL QUIT ("Linear symmetry requires NZ=1")
C
#if defined (SYS_IBM)
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
         OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
C
C     Get integrals of component of angular momentum that we want to
C     quantize on.
C
      FIRST = .TRUE.
      IOP =  IPANGMOM
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK(KFREE),LFREE,IPRHAM)
      CALL DSPTSI(NTBAS(0),FTRI,FOCK)
C     
C     Get the spin-z operator
C     =======================
C     
      FIRST = .TRUE.
      IOP =  IPSPNMOM
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK(KFREE),LFREE,IPRHAM)
      CALL MEMGET('REAL',KBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL DSPTSI(NTBAS(0),FTRI,WORK(KBUF))
C     
C     Make the Jz operator (in case of spinfree we make Lz)
C     =====================================================
C     
      IF (.NOT.(SPINFR.OR.LEVYLE))
     &   CALL DAXPY (N2BBASX,D1,WORK(KBUF),1,FOCK,1)
C
C     Release memory & close AOPROPER
C
      CALL MEMREL('AMJMA1',WORK,KBUF,KBUF,KFREE,LFREE)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(FOCK,1,WORK(KFREE),LFREE)
C
C     Output section
C     ==============
C
      IF(IPRHAM.GE.4) THEN
        CALL HEADER('One-electron MJ matrix',-1)
        DO 10 I = 1,NFSYM
          IF(NFBAS(I,0).EQ.0) GOTO 10
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL PRQMAT(FOCK(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI)
   10   CONTINUE
      ENDIF
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Wrifrz */
      SUBROUTINE WRIFRZ(IUNIT,IFRP,VECS,IVEC,BUF,IBUF,NVEC)
C***********************************************************************
C
C     Dump frozen orbitals to file
C     Written by T. Saue May 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
      DIMENSION VECS(*),BUF(*),IVEC(*),IBUF(*)
C
C     Get index array for frozen orbitals; 
C
      NTOT=1
      CALL NUMLS3(VCFROZ(IFRP),IVEC,NFBAS(IFRP,0),
     &            -NPSH(IFRP),NESH(IFRP),NTOT,NZERO)
C
C     Consistency check
C
      IF(NFRO(IFRP).NE.(NTOT-NZERO)) THEN
        CALL QUIT('WRIFRZ:Inconsistency !')
      ENDIF
C
C     Selected orbitals that are not frozen have index zero;
C     they must be taken out and pointers made to the remaining ones
C     Write compressed index array to file
C
C     We recompute NASHMFT as this value may change if there are frozen open shells
C
      NASHMFT = 0
      II = 0
      DO I = 1,NTOT
      IP = IVEC(I)
      IF(IP.NE.0) THEN
        II = II + 1
        IBUF(II)=I        
C       Positronic position
        IF(IP.LT.0.AND.IP.GE.(-NPSH(IFRP))) THEN
          CALL QUIT('WRIFRZ: Positronic frozen orbitals not allowed !')
          IVEC(II)=NPSH(IFRP)+1+IP
C       Electronic position
        ELSEIF (IP.GT.0.AND.IP.LE.NESH(IFRP)) THEN
          IVEC(II)=NPSH(IFRP)+IP
C         Check if orbital belongs to closed or open shells.
C         If so, subtract from occupation number
          IF(IP.GT.0.AND.IP.LE.NISH(IFRP)) THEN
            NISHMF(IFRP) = NISHMF(IFRP) - 1
          ELSE
            I1 = NISH(IFRP)
            DO IOPEN = 1,NOPEN
              I2 = I1 + NACSH(IFRP,IOPEN)
              IF(IP.GT.I1.OR.IP.LE.I2) THEN
                NACSHMF(IFRP,IOPEN) = NACSHMF(IFRP,IOPEN) - 1
                GOTO 20
              ENDIF
              I1 = I2
            ENDDO
 20         CONTINUE
          ENDIF
        ENDIF
      ENDIF
      ENDDO
      DO IOPEN = 1,NOPEN
        NASHMFT = NASHMFT + NACSHMF(IFRP,IOPEN)
      ENDDO
      CALL WRITI(IUNIT,NFRO(IFRP),IVEC)
C
C     Write coefficients to file
C
      IDIM = NFBAS(IFRP,0)*NVEC
      JDIM = NFBAS(IFRP,0)*NZ
      DO II = 1,NFRO(IFRP)
        I = IBUF(II)
        IOFF = NFBAS(IFRP,0)*(I-1)+1
        JOFF = 1
        DO IZ = 1,NZ
          CALL DCOPY(NFBAS(IFRP,0),VECS(IOFF),1,BUF(JOFF),1)
          IOFF = IOFF + IDIM
          JOFF = JOFF + NFBAS(IFRP,0)
        ENDDO
        CALL WRITT(IUNIT,JDIM,BUF)
      ENDDO
C
      RETURN
      END
