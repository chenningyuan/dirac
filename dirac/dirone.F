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
C FILE    : dirone.F
C SOURCE  : $Source$
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
#include <implicit.h>
#include <priunit.h>
C
#include <mxcent.h>
      DIMENSION WORK(LWORK)
C
#include <dcbgen.h>
#include <dcbprl.h>
#include <dcbham.h>
#include <dcbcls.h>
#include <cbihr1.h>
#include <cbiher.h>
C
C     Loop through list of class labels; calculate only the active ones(CLSCAL)
C     =========================================================================
C
      DO 10 I = 1,NPRPCLS
      IF(.NOT.CLSCAL(I)) GOTO 10
C
C     Nuclear attraction integrals (molecular field)
C
      IF    (CLSINT(I).EQ.'MOLFLD ') THEN
        IF (NOSMLV) CLSCMB(I) = 'TFFF'
        IF (NOSMLV) PDOINT(I) = 'TFFF'
        CALL NUCATT(WORK,LWORK,CLSCMB(I),IPRONE)
C
C     Contribution from beta-matrix
C
      ELSEIF(CLSINT(I).EQ.'BETAMAT') THEN
        CALL BETAMT(CLSCMB(I),WORK,LWORK,IPRONE)
C
C     All other integrals
C
      ELSE
        CALL OP1GEN(CLSINT(I),CLSCMB(I),IORDCL(I),WORK,LWORK,IPRONE)
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
C     Last revision May 19 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dummy.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
      LOGICAL DOINT(4),TOFILE
      CHARACTER WORD*7,COMB*4,LABINT(3*MXCOOR)*8
      DIMENSION WORK(LWORK)
#include <dcbgen.h>
#include <nuclei.h>
#include <symmet.h>
#include <cbiher.h>
C
      CALL QENTER('OP1GEN')
#include <memint.h>
C     Write integrals to file
      TOFILE   = .TRUE.
      NCOMP    = 0
      MAXTYP = 3*MXCOOR
      IF ((WORD.EQ.'ELFGRDC').OR.(WORD.EQ.'ELFGRDS')) THEN
         NINTAD = 9*NUCIND*(MAXREP+1)
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
      CALL PR1IN1(WORK(KFREE),LFREE,WORK(KINTRP),WORK(KINTAD),
     &            LABINT,WORD,IORDER,NPQUAD,TRIANG,
     &            PROPRI,IPRINT,DUMMY,NCOMP,TOFILE,'TRIANG',
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
#include <implicit.h>
#include <priunit.h>
C
#include <dummy.h>
      LOGICAL DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      DIMENSION WORK(LWORK)
#include <dcbgen.h>
#include <dcbbas.h>
C
      CALL QENTER('BETAMT')
#include <memint.h>
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
      CALL GETINT('OVERLAP',WORK(KONEIN),DOINT,NCOMP,IDUMMY,IDUMMY,
     &            'TRIANG',.FALSE.,IPRINT,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Small component overlap integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),NTBAS(0),LUPRI)
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
#include <implicit.h>
#include <priunit.h>
C
      LOGICAl DOINT(4)
      CHARACTER LABEL*8,RTNLBL(2)*8,COMB*4
      DIMENSION WORK(LWORK)
#include <dcbbas.h>
C
      CALL QENTER('NUCATT')
#include <memint.h>
C     Memory allocation
      CALL MEMGET('REAL',KONEIN,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate integrals
C
      READ(COMB,'(4L1)',ERR=1000) (DOINT(I),I=1,4)
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,IPRINT)
      IF(IPRINT.GE.6) THEN
        CALL HEADER('Nuclear attraction integrals',-1)
        CALL OUTPAK(WORK(KONEIN),NTBAS(0),NTBAS(0),LUPRI)
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
      SUBROUTINE NUCAT1(OP1INT,WORK,LWORK,DOINT,IPRINT)
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <mxcent.h>
C
      LOGICAL DOINT(2,2)
      DIMENSION OP1INT(NNBBASX),WORK(LWORK)
#include <dcbbas.h>
#include <nuclei.h>
#include <cbiher.h>
#include <cbihr1.h>
#include <cbisol.h>
C
      CALL QENTER('NUCAT1')
#include <memint.h>
C
      CALL MEMGET('REAL',KSTHMA,3*NNBBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFACIN,  NUCDEP,WORK,KFREE,LFREE)
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
     &            NNBBASX,NNBBAST,DOINT,.FALSE.)
C
C     Memory deallocation
C
      CALL MEMREL('NUCAT1.ONEDR1',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Symmetry unpack nuclear attraction integrals
C
      KNUCS = KSTHMA + NNBBASX
      CALL SYMUPK(WORK(KNUCS),OP1INT,1,NNBBASX)
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <cbihr1.h>
      DIMENSION FOCK(*),WORK(LWORK)
C
      CALL QENTER('ONEFCK')
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbgen.h>
#include <dcbham.h>
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
C     Nuclear attraction integrals
C     ============================
C
      IOP = IPMOLFLD
      CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
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
C     Finite field: 
C     Loop through list of one-electron operators in CBIHAM
C     =====================================================
C
      DO I = 1,N1OPER
        IOP = IPR1OP(I)
        CALL PRPMAO(LU1INT,IOP,FTRI,OP1INT,FIRST,WORK,LWORK,IPRHAM)
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
      CALL BUTOBS(FOCK,NZ,WORK,LWORK)
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
      SUBROUTINE PRPMAO(LU1INT,IOP,PTRI,OP1INT,FIRST,WORK,LWORK,IPRINT)
C*****************************************************************************
C     Generates lower triangular matrix of one-electron integrals
C     corresponding to property IOP of COMMON block XCBPRP
C
C     Written by Trond Saue May 1996
C     Last revision: May 19 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,PTOL = 1.0D-12)
C
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbqal.h>
#include <dcbxpr.h>
#include <dcbprl.h>
C
      CHARACTER MXFORM*6,FMT*6
      LOGICAL FIRST(*),DOINT(4)
      DIMENSION PTRI(NNBBASX,*),OP1INT(NNBBASX),WORK(LWORK)
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
      IF(IOPTYP.EQ.1) THEN
        IC = 1
      ELSE
        IC = 2
      ENDIF
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
        IM4  = JM4(I,IOPTYP)
        IBRP = IPRLREP(KLBL)
        IH   = IPRLTYP(KLBL)
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
          WRITE(LUPRI,'(A,4X,'//FMT//')')
     &      'Factor: ',PFAC
          WRITE(LUPRI,'(A,I5)') 'Integrals enter matrix :',IQP
          WRITE(LUPRI,'(A,L1)') 'FIRST = ',FIRST(IQP)
        ENDIF
C
C       Get integral from file
C
        CALL PRPREA(LU1INT,PRPLBL(KLBL),PDOINT(KLBL),
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
C       Phase transformation if necessary
C
        IF(NZ.LT.4) THEN
          READ(PDOINT(KLBL),'(4L1)',ERR=1000) (DOINT(J),J=1,4)
          CALL QNPHASE(IZ,IBRP,OP1INT,DOINT)
          IF(IPRINT.GE.8) THEN
             WRITE(LUPRI,'(A)') '* Phase transformed integrals:'
             CALL OUTPAK(OP1INT,NTBAS(0),1,LUPRI)
          ENDIF
        ENDIF
C
C       Determine where integrals go and add them on
C
        IF(FIRST(IQP)) THEN
          CALL DCOPY(NNBBASX,OP1INT,1,PTRI(1,IQP),1)
          FIRST(IQP) = .FALSE.
        ELSE
          CALL DAXPY(NNBBASX,D1,OP1INT,1,PTRI(1,IQP),1)
        ENDIF
 10   CONTINUE
      IF(IPRINT.GE.8) THEN
        DO IZ = 1,NZ
          WRITE(LUPRI,'(A,I2,A)') '* Triangular prop.matrix(IZ=',IZ,')'
          IF (FIRST(IZ)) THEN
             WRITE(LUPRI,'(/A)') '  --- zero matrix ---'
          ELSE
             CALL OUTPAK(PTRI(1,IZ),NTBAS(0),1,LUPRI)
          END IF
        ENDDO
      ENDIF
      RETURN
C
 1000 CONTINUE
C
C     Not able to read DOINT information from PDOINT
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
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
C
#include <ccom.h>
#include <dcbgen.h>
#include <dcbbas.h>
      DIMENSION VMAT(*),WORK(LWORK)
C
      CALL QENTER('LOWDIN')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
#include <dummy.h>
#include <maxaqn.h>
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include <ccom.h>
#include <dcbgen.h>
#include <dcbham.h>
#include <dcborb.h>
#include <dgroup.h>
#include <pgroup.h>
#include <dcbbas.h>
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
     &                     STOL(IC),IPRINT,EIG,WORK,LWORK)
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
C  /* Deck getint */
      SUBROUTINE GETINT(WORD,OP1INT,DOINT,NCOMP,IORDER,NPQUAD,
     &                  MXFORM,TOFILE,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for getting one-electron property integrals
C     Written by T.Saue October 1994
C     Last revision : Nov 15 1994
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
#include <shells.h>
#include <symmet.h>
#include <nuclei.h>
      DIMENSION OP1INT(*),WORK(LWORK)
      CHARACTER WORD*7,MXFORM*6,LABINT(3*MXCOOR)*8
      LOGICAL DOINT(2,2),PROPRI,TRIANG,TOFILE
C
      CALL QENTER('GETINT')
#include <memint.h>
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
      CALL MEMREL('GETINT.PR1IN1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('GETINT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mkovlp */
      SUBROUTINE MKOVLP(WORK,LWORK)
C*****************************************************************************
C
C     Generate sorted overlap matrix
C
C     Written by T.Saue Sep 1995
C     Last revision: Sep 11 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C Used from COMMON blocks
C  DCBBAS: N2BBASX
C  PAMIOU: LUOVLP
C
#include <dcbbas.h>
#include <dcbgen.h>
      DIMENSION WORK(LWORK)
#include <memint.h>
C     Memory allocation
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSBUF,NNBBASX,WORK,KFREE,LFREE)
C
C     Generate sorted overlap matrix
C     ==============================
C
      CALL MKOVL1(WORK(KSMAT),WORK(KSBUF),WORK(KFREE),LFREE)
C
C     Write overlap matrix to file
C     ============================
C
      CALL OPNFIL(LUOVLP,'DFOVLP','UNKNOWN','MKOVLP')
      CALL WRITT(LUOVLP,N2BBASX,WORK(KSMAT))
      CLOSE(LUOVLP,STATUS='KEEP')
C
C     Memory deallocation
C     ===================
C
      CALL MEMREL('MKOVLP.MKOVL1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mkovl1 */
      SUBROUTINE MKOVL1(SMAT,SBUF,WORK,LWORK)
C*****************************************************************************
C
C     Generate sorted overlap matrix
C
C     Written by T.Saue Sep 1995
C     Last revision: Sep 11 1995 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
C Used from COMMON blocks:
C  DCBBAS: NTBAS(0)
C
#include <dcbgen.h>
#include <dcbham.h>
#include <dcbbas.h>
#include <cbihr1.h>
      DIMENSION SMAT(*),SBUF(*),WORK(LWORK)
C
C     Read triangular batch of overlap integrals from file
C     ====================================================
C
#if defined (SYS_IBM)
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='APROPER')
#else
      OPEN (LU1INT,STATUS='UNKNOWN',FORM='UNFORMATTED',
     &         FILE='AOPROPER')
#endif
      CALL PRPREA(LU1INT,'OVERLAP ','TFFT',SBUF,NNBBASX,IPRONE)
      CLOSE(LU1INT,STATUS = 'KEEP')
C     
C     Make full matrix
C     ================
C
      CALL DSPTSI(NTBAS(0),SBUF,SMAT)
C
      IF(IPRONE.GE.5) THEN
        CALL HEADER('Unsorted overlap matrix',-1)
        CALL OUTPUT(SMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),
     &              NTBAS(0),1,LUPRI)
      ENDIF
C
C     Reindex overlap matrix to sorted basis
C     ======================================
C
      CALL BUTOBS(SMAT,1,WORK,LWORK)
      IF(IPRONE.GE.4) THEN
        CALL HEADER('Overlap matrix',-1)
        CALL OUTPUT(SMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),
     &              NTBAS(0),1,LUPRI)
      ENDIF
C
      RETURN
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
#include <implicit.h>
#include <priunit.h>
C Used from COMMON blocks:
C  DCBBAS: NTBAS(0)
C
#include <dcbbas.h>
      DIMENSION VMAT(*),VBUF(*),WORK(LWORK)
      CALL QENTER('VFSORT')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
C Used from COMMON blocks:
C  DCBBAS: NTBAS,NFORB
C
#include <dcbbas.h>
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
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION FOCK(*),WORK(LWORK)
#include <dcbgen.h>
#include <cbihr1.h>
#include <dcbbas.h>
#include <dgroup.h>
C
      CALL QENTER('FREEMT')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
#include <dcbham.h>
#include <dcbgen.h>
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
#include <implicit.h>
#include <priunit.h>
C
C Used from COMMON blocks:
C   DCBBAS: N2BBASX
C   PAMIOU: LUTMAT
#include <dcbgen.h>
#include <dcbham.h>
#include <cbihr1.h>
#include <dcbdhf.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      DIMENSION WORK(LWORK)
      CALL QENTER('GMOTRA')
#include <memint.h>
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
C
C     Output
C
      CALL HEADER('Output from GMOTRA',-1)
      IF ( NFSYM .EQ. 1 ) THEN
         WRITE(LUPRI,'(40X,A)') 'Irrep 1'
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
C     Write the transformation matrix.
C
      N2TMT = 0
      DO IFRP = 1,NFSYM
        I2TMT(IFRP) = N2TMT
        N2TMT = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZT
      ENDDO
      CALL OPNFIL(LUTMAT,'DFTMAT','UNKNOWN','GMOTRA')
      CALL WRITT(LUTMAT,N2TMT,WORK(KTMAT))
C
C     Construct right index transformed overlap matrix Saomo = Sao C
C
      IF(SPINFR.OR.DIISAO.OR.NAELEC.GT.0) THEN
        ISIZE = 0
        DO I = 1,NFSYM
          ISIZE = ISIZE + NZ*NORB(I)*NFBAS(I,0)
        END DO
C       ISIZE .le. N2BBASXQ .le. size of DTMP
        CALL MKSAOMO(WORK(KTBUF),WORK(KTMAT),IPRONE,WORK(KFREE),LFREE)
        CALL WRITT(LUTMAT,ISIZE,WORK(KTBUF))
      ENDIF
C
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
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,DM2 = -2.0D0,D0 = 0.0D0)
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
#include <pgroup.h>
#include <cbihr1.h>
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
#include <memint.h> 
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
        FACSL = - D1 / FAC
        IBUF  = 1
        IMAT  = 1
        DO 20 I = 1,NFSYM
          NESH(I)  = NFORB(I,1)
          NPSH(I)  = 0
          NORB(I)  = NESH(I) + NPSH(I)
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
              CALL DSCAL (NFORB(I,2),FACSL,TBUF(JBUF),1)
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
          CALL BCKTRA(TMAT(IMAT),NFBAS(I,0),NFORB(I,1),
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
        NZT  = NZ
        IMAT = 1
        IBUF = 1
        DO 40 I = 1,NFSYM
          NESH(I)  = NFORB(I,1)
          NPSH(I)  = NFORB(I,1)
          NORB(I)  = NESH(I) + NPSH(I)
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
C         Transform to free particle/bare nucleus matrix
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
            CALL BCKTRA(TMAT(IMAT),NFBAS(I,0),NORB(I),
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
          CALL BCKTRA(WORK(IBUF),NFBAS(I,0),NESH(I),
     &              TBUF(I2ORBT(I)+1),NORB(I),NORB(I),
     &              NORB(I),NZ,NESH(I),(NPSH(I)+1),NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZT,
     &              IPRHAM)
          NPSH(I) = 0
          NORB(I) = NESH(I)
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
            CALL BCKTRA(TMAT(IMAT),NFBAS(I,0),NORB(I),
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
          WRITE(LUPRI,'(A,I5)') '* Fermion ircop no.',I
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
 80     CONTINUE
      ENDIF
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck prprea */
      SUBROUTINE PRPREA(LU1INT,PRPLBL,COMB,PRPINT,NDIM,IPRINT)
C*****************************************************************************
C
C     Read property integrals from file
C
C     Written by T.Saue, January 1995
C     Last revision: May 19 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      LOGICAL FNDLB3
      CHARACTER PRPLBL*8, RTNLBL(2)*8, COMB*4
      DIMENSION PRPINT(NDIM)
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fndlb3 */
      LOGICAL FUNCTION FNDLB3(SRCLBL,COMB,RTNLBL,LU)
C*****************************************************************************
C
C  18-Jan-1995 tsaue
C  1. Find MOLECULE label and block COMB on unformatted file LU
C       a)Cannot find label: return FALSE
C       b)Finds label      : continues
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
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
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
C
#include <dcbgen.h>
#include <cbihr1.h>
#include <dcbbas.h>
#include <dgroup.h>
      DIMENSION SMAT(N2BBASX),WORK(LWORK)
      CALL QENTER('GTOVLX')
#include <memint.h>
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
      CALL PRPREA(LU1INT,'OVERLAP ','TFFT',WORK(KSBUF),NNBBASX,IPRONE)
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
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0)
#include <mxcent.h>
#include <maxaqn.h>
#include <mxorb.h>
C
#include <dcbgen.h>
#include <cbihr1.h>
#include <dcbbas.h>
#include <dgroup.h>
#include <symmet.h>
      DIMENSION SMAT(*),WORK(LWORK)
      CALL QENTER('GTOVLT')
#include <memint.h>
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
      CALL PRPREA(LU1INT,'OVERLAP ','TFFT',WORK(KSTRI),NNBBASX,IPRONE)
      CLOSE(LU1INT,STATUS = 'KEEP')
C
C     Make full matrix
C     ================
C
      CALL DSPTSI(NTBAS(0),WORK(KSTRI),WORK(KSBUF))
C
      CALL DZERO(SMAT,N2BBASX)
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <dcbham.h>
#include <dcbbas.h>
#include <dgroup.h>
      DIMENSION WORK(LWORK)
      CALL QENTER('HEINZ')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
#include <alphac.h>
      PARAMETER(DP5 = 0.5D0,D0 = 0.0D0)
C
#include <dcbgen.h>
#include <dcbham.h>
#include <dcbbas.h>
#include <dcborb.h>
#include <dgroup.h>
      LOGICAL FIRST(NZ)
      DIMENSION CMO(*),EIG(*),FAO(N2BBASX,NZ),FMO(*),FBUF(N2BBASX,NZ),
     &          VBUF(N2BBASX),TMAT(*),OP1INT(*),FTRI(NNBBASX,NZ),
     &          WORK(LWORK),IQBF(4)
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
C  /* Deck reno */
      SUBROUTINE RENO(LDM,NDIM,SMAT)
#include <implicit.h>
#include <priunit.h>
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
#include <implicit.h>
#include <priunit.h>
C
C
#include <dcbbas.h>
#include <dgroup.h>
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(*),WORK(*)
C
      CALL QENTER('SPHLOW')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D6=6.0D0,D4=4.0D0)
C
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
C
#include <ccom.h>
#include <symmet.h>
#include <pgroup.h>
#include <shells.h>
#include <nuclei.h>
#include <dgroup.h>
#include <dcbrkb.h>
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(LMAX),IFUN(NFUN,3),
     &          TMAT(NCAR,NCAR),CBUF(*),ICBF(3,LMAX),
     &          WORK(*)
      DATA COMP/'L','S'/
#include <sdpre.h>
#include <ibtfun.h>
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
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <dgroup.h>
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),
     &          CSPH(*),ICSP(*),WORK(*)
C
      CALL QENTER('RKBLOW')
#include <memint.h>
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
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D6=6.0D0,D4=4.0D0)
C
#include <mxcent.h>
#include <mxorb.h>
#include <maxaqn.h>
C
#include <ccom.h>
#include <symmet.h>
#include <pgroup.h>
#include <shells.h>
#include <nuclei.h>
#include <dgroup.h>
#include <dcbrkb.h>
      LOGICAL LRKB
      CHARACTER COMP(2)*1
      DIMENSION SMAT(LDM,*),VMAT(LDM,*),EIG(*),IQBF(4),
     &          CSPH(*),ICSP(LMAX),IFUN(NFUN,3),
     &          TMAT(NCAR,NCAR),CRED(*),IRED(3,LMAX),
     &          WORK(*)
      DATA COMP/'L','S'/
#include <sdpre.h>
#include <ibtfun.h>
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
C           Look for (L-1) large function;
C           skip if no partner found
C
            JL = IRKB(JP)
            IF(JL.EQ.0) GOTO 15
            IROW = NSTRT(JA + JP - JP1) + IFUN(LOFF+1,1)
C                        ^- to get correct shell number
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
               JB = IRKB(JL)
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
C     Written by T.Saue Sep 11 1995
C     Last revision Apr 1 1998 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
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
#include <implicit.h>
#include <priunit.h>
#include <mxcent.h>
      PARAMETER(D1 = 1.0D0)
C
      LOGICAL DOINT(2,2)
      DIMENSION WORK(LWORK),CBUF(2),GBUF(2),FMAT(*)
#include <dgroup.h>
#include <dcbbas.h>
#include <nuclei.h>
C
      CALL QENTER('TSCORR')
#include <memint.h>
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
      CHARGE(1) = -0.1812D0
      CHARGE(2) = -0.0832D0
C
C     Generate integrals
C
      DOINT(1,1) = .FALSE.
      DOINT(2,2) = .TRUE.
      DOINT(1,2) = .FALSE.
      DOINT(2,1) = .FALSE.
      CALL NUCAT1(WORK(KONEIN),WORK(KFREE),LFREE,
     &            DOINT,IPRINT)
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
