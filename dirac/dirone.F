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
        IF(IQP.GT.NZ) THEN
          WRITE(LUPRI,'(A)') '***** PRPMAO ERROR *****'
          WRITE(LUPRI,'(A)') PRPLBL(KLBL)
          WRITE(LUPRI,'(A,I2,A,I2,A)') 'IQP(',IQP,').GT.MZ(',
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
          WRITE(LUPRI,'(A,I5)') 'Integrals enters matrix :',IQP
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
          CALL OUTPAK(PTRI(1,IZ),NTBAS(0),1,LUPRI)
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
C
#include <dcbbas.h>
      DIMENSION VMAT(*),WORK(LWORK)
C
      CALL QENTER('LOWDIN')
#include <memint.h>
C
C     Memory allocation
C
      CALL MEMGET('REAL',KVBUF,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSMAT,N2BBASX,WORK,KFREE,LFREE)
C
C     Make Lowdin matrix
C     ==================
C
      CALL LOWDI1(VMAT,WORK(KSMAT),WORK(KVBUF),
     &            WORK(KFREE),LFREE,IPRINT)
C
C     Memory deallocation
C
      CALL MEMREL('LOWDIN.LOWDI1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('LOWDIN')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowdi1 */
      SUBROUTINE LOWDI1(VMAT,SMAT,VBUF,WORK,LWORK,IPRINT)
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
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
#include <dcbgen.h>
#include <dgroup.h>
#include <pgroup.h>
#include <dcbbas.h>
      CHARACTER COMP(2)*1
      DIMENSION VMAT(*),SMAT(*),VBUF(*),WORK(LWORK)
      LOGICAL LBIT
C
      CALL QENTER('LOWDI1')
      IF(IPRINT.GE.5) CALL TITLER('Output from LOWDI1','*',103)
      COMP(1) = 'L'
      COMP(2) = 'S'
C
C     Get overlap matrix and transform to AO-basis
C     ==================================================
C
      CALL GTOVLP(SMAT,WORK,LWORK)
      IF(NBSYM.GT.1) THEN
        CALL MTBSBU(SMAT,VMAT)
        CALL MTSOAO(VMAT,SMAT,NTBAS(0),0,IPRINT)
      ENDIF
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Overlap matrix',-1)
        CALL OUTPUT(SMAT,1,NTBAS(0),1,NTBAS(0),NTBAS(0),NTBAS(0),
     &              1,LUPRI)
      ENDIF
C
C     Initialize Lowdin matrix
C     ========================
C
      CALL DZERO(VMAT,N2BBASX)
C
C     Large components
C     ================
C
      NLEFF = 0
      IF(NTBAS(1).EQ.0) GOTO 10
      IF(LBIT(ISPHTR,1)) THEN
        CALL SPHLOW(1,SMAT,VMAT,NTBAS(0),NLEFF,STOL(1),NDEL,
     &              IPRINT,WORK,LWORK)
      ELSE
        CALL LOWGEN(1,SMAT,VMAT,NTBAS(0),NTBAS(1),NLEFF,
     &              STOL(1),IPRINT,WORK,LWORK)
        NDEL = NTBAS(1) - NLEFF
      ENDIF
      IF(NDEL.GT.0) THEN
         WRITE(LUPRI,'(A,I5)')
     &  '*** WARNING *** Eigenvalues below threshold removed from'//
     &  ' large component overlap matrix:',NDEL
         IF (ISPHTR.EQ.3) WRITE(LUPRI,'(A)')
     &   '*** CARAMBA *** THIS MAY VERY WELL NOT GIVE CORRECT RKB!!'
      ENDIF
C     
C     Small components
C     ================
C
 10   CONTINUE
      NSEFF = 0
      IF(NTBAS(2).EQ.0) GOTO 20
      IOFF = NTBAS(0)*NTBAS(1) + NTBAS(1) + 1
      JOFF = NTBAS(0)*NLEFF    + NTBAS(1) + 1
      IF(ISPHTR.EQ.3) THEN
         CALL RKBLOW(2,SMAT(IOFF),VMAT(JOFF),NTBAS(0),NSEFF,
     &        STOL(2),NDEL,IPRINT,WORK,LWORK)
      ELSEIF(ISPHTR.EQ.2) THEN
         CALL SPHLOW(2,SMAT(IOFF),VMAT(JOFF),NTBAS(0),NSEFF,
     &        STOL(2),NDEL,IPRINT,WORK,LWORK)
      ELSE
         CALL LOWGEN(2,SMAT(IOFF),VMAT(JOFF),NTBAS(0),NTBAS(2),
     &        NSEFF,STOL(2),IPRINT,WORK,LWORK)
         NDEL = NTBAS(2) - NSEFF
      ENDIF
      IF(NDEL.GT.0) THEN
         WRITE(LUPRI,'(A,I5)')
     &  '*** WARNING *** Eigenvalues below threshold removed from'//
     &  ' small component overlap matrix:',NDEL
         IF (ISPHTR.EQ.3) WRITE(LUPRI,'(A)')
     &  '*** CARAMBA *** THIS MAY VERY WELL NOT GIVE CORRECT RKB!!'
      ENDIF
C
C     Backsort to SO-basis
C     ====================
C
 20   CONTINUE
      IF(NBSYM.EQ.1) THEN
        NFORB(1,1) = NLEFF
        NFORB(1,2) = NSEFF
        NFORB(1,0) = NFORB(1,1) + NFORB(1,2)
        N2VMT      = NFORB(1,0)*NFBAS(1,0)
      ELSE
        NTEFF = NLEFF + NSEFF
        CALL VTAOSO(VMAT,VBUF,NTBAS(0),NTEFF,IPRINT)
        CALL VTBUBS(VBUF,VMAT,NTEFF)
C
C       Shuffle on fermion symmetries
C       =============================
C
        IF(NFSYM.EQ.1) THEN
          NFORB(1,1) = NLEFF
          NFORB(1,2) = NSEFF
          NFORB(1,0) = NFORB(1,1) + NFORB(1,2)
          N2VMT      = NFORB(1,0)*NFBAS(1,0)
        ELSE
          CALL DCOPY(NTEFF*NTBAS(0),VMAT,1,VBUF,1)
          CALL VFSORT(NLEFF,NSEFF,VMAT,VBUF,WORK,LWORK)
          N2VMT = NFORB(1,0)*NFBAS(1,0) + NFORB(2,0)*NFBAS(2,0)
        ENDIF
      ENDIF
C
C     Output section
C     ==============
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('Lowdin matrix',-1)
        IOFF = 1
        DO 100 I = 1,NFSYM
          IF(NFORB(I,0).EQ.0) GOTO 100
          WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
          CALL OUTPUT(VMAT(IOFF),1,NFBAS(I,0),1, NFORB(I,0),
     &                NFBAS(I,0),NFORB(I,0),1,LUPRI)
          IOFF = IOFF + NFBAS(I,0)*NFORB(I,0)
  100   CONTINUE
      ENDIF
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
C  /* Deck ltest */
      SUBROUTINE LTEST(WORK,LWORK)
#include <implicit.h>
#include <priunit.h>
#include <pgroup.h>
#include <dgroup.h>
#include <dcbdhf.h>
#include <dcbbas.h>
      DIMENSION WORK(LWORK)
      CALL QENTER('LTEST')
#include <memint.h>
      CALL MEMGET('REAL',KSSO,N2BBASX,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSAO,N2BBASX,WORK,KFREE,LFREE)
C
C     Get overlap matrix and transform to unsorted basis
C     ==================================================
C
      CALL GTOVLP(WORK(KSSO),WORK(KFREE),LFREE)
      CALL BSTOBU(WORK(KSSO),1,WORK(KFREE),LFREE)
C
C     Transform overlap matrix to AO-basis
C     ====================================
C
      CALL MTSOAO(WORK(KSSO),WORK(KSAO),NTBAS(0),0,10)
      CALL DZERO(WORK(KSSO),N2BBASX)
      CALL MTAOSO(WORK(KSAO),WORK(KSSO),NTBAS(0),0,10)
      CALL MEMREL('LTEST',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('LTEST')
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
C     Last revision: Jan 7 1998 - jth LINUX
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C Used from COMMON blocks:
C  DCBBAS: NTBAS,NFORB
C
#include <dcbbas.h>
C
#ifdef VAR_G77
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
      IF(.NOT.URKBAL.AND.NTBAS(2).EQ.0) THEN
        URKBAL = .TRUE.
        WRITE(LUPRI,'(A)') 
     &   '*** WARNING *** : No small bases. Setting URKBAL = T !'
      ENDIF
      IF(URKBAL) THEN
        NZT = 1
        N2TMT = 0
        DO I = 1,NFSYM
          N2TMT = N2TMT + NFBAS(I,0)*NFORB(I,0)
        ENDDO
        CALL OPNFIL(LUTMAT,'DFTMAT','UNKNOWN','GMOTRA')
        CALL WRITT(LUTMAT,N2TMT,WORK(KVMAT))
        CLOSE(LUTMAT,STATUS='KEEP')
      ELSE
C
C     Charge conjugation symmetry
C     ===========================
C
        CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)
        CALL DELPOS(WORK(KTMAT),WORK(KVMAT),WORK(KTBUF),WORK(KEIG),
     &            WORK(KFREE),LFREE)
        N2TMT = 0
        NFAC = 2
        IF(FREEPJ.OR.VEXTPJ) THEN
          DO I = 1,NFSYM
            N2TMT = N2TMT + NFBAS(I,0)*NFORB(I,1)
          ENDDO
        ELSE
          DO I = 1,NFSYM
            NDIM  = NFORB(I,1) + MIN(NFORB(I,2),NFORB(I,1))
            N2TMT = N2TMT + NFBAS(I,0)*NDIM
          ENDDO
        ENDIF
        N2TMT = NZ*N2TMT
        CALL OPNFIL(LUTMAT,'DFTMAT','UNKNOWN','GMOTRA')
        CALL WRITT(LUTMAT,N2TMT,WORK(KTMAT))
        CLOSE(LUTMAT,STATUS='KEEP')
      ENDIF
C
C     Memory deallocation
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
C
C     Written by T.Saue, October 1995
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
#include <cbihr1.h>
C
      DIMENSION TMAT(*),TBUF(*),VMAT(*),EIG(*),WORK(LWORK)
      DIMENSION IQBF(4)
C
C     Construct free particle matrix
C
      IF(VEXTPJ) THEN
        CALL ONEFCK(TBUF,IPRONE,WORK,LWORK)
      ELSE
        CALL FREEMT(TBUF,IPRONE,WORK,LWORK)
      ENDIF
C
C       Transform to MO-basis
C       =====================
C
      IOFT = 1
      IOFV = 1
      DO 10 I = 1,NFSYM
        IF(NFORB(I,0).EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                         NFORB(I,0),NFORB(I,0),
     &              TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),
     &              TMAT(IOFT),NFORB(I,0),NFORB(I,0),NBF,
     &              IQBF,
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,IPQTOQ(1,0),
     &              WORK,LWORK,IPRGEN)
        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** DELPOS: NBF.NE.NZ !'
   10 CONTINUE
C
C       Diagonalize Fock-matrix in MO-basis
C       ===================================
C
      IOFT = 1
      IEIG = 1
      DO 20 I = 1,NFSYM
        IF(NFORB(I,0).EQ.0) GOTO 20
        IF(DOJACO) THEN
          CALL RSJACO(NFORB(I,0),NFORB(I,0),NFORB(I,0),TMAT(IOFT),
     &                EIG(IEIG),1,1,0,TBUF(IOFT),WORK,LWORK)
        ELSE
          CALL QDIAG(NZ,NFORB(I,0),TMAT(IOFT),NFORB(I,0),NFORB(I,0),
     &             EIG(IEIG),1,TBUF(IOFT),NFORB(I,0),
     &             NFORB(I,0),WORK,LWORK,IERR)
        ENDIF
        IOFT = IOFT + NFORB(I,0)*NFORB(I,0)*NZ
        IEIG = IEIG + NFORB(I,0)
   20 CONTINUE
C
C     Generate transformation matrix
C     ==============================
C
      IOFT = 1
      IOFU = 1
      IOFV = 1
      NFAC = 2
      IF(FREEPJ.OR.VEXTPJ) NFAC = 1
      DO 30 I = 1,NFSYM
        IF(NFORB(I,0).EQ.0) GOTO 30
        NTDIM = NFAC*NFORB(I,1)
C
C       Backtransform accepted positronic solutions
C       ===========================================
C
        IF(FREEPJ.OR.VEXTPJ) THEN
          IOFTE = IOFT
        ELSE
          CALL BCKTRA(TMAT(IOFT),NFBAS(I,0),NTDIM,
     &              TBUF(IOFU),NFORB(I,0),NFORB(I,0),
     &              NFORB(I,0),NZ,IPQTOQ(1,0),
     &              NFORB(I,1),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &              IPQTOQ(1,0),IPRGEN)
          NDIM  = MIN(NFORB(I,1),NFORB(I,2))
          IOFTE = IOFT + NFBAS(I,0)*NDIM
        ENDIF
        IOFUE = IOFU + NFORB(I,0)*NFORB(I,2)
C
C       Backtransform electronic solutions
C       ===========================================
C
        CALL BCKTRA(TMAT(IOFTE),NFBAS(I,0),NTDIM,
     &              TBUF(IOFUE),NFORB(I,0),NFORB(I,0),
     &              NFORB(I,0),NZ,IPQTOQ(1,0),
     &              NFORB(I,1),1,NFBAS(I,0),
     &              VMAT(IOFV),NFBAS(I,0),NFORB(I,0),1,
     &              IPQTOQ(1,0),IPRGEN)
        IOFT = IOFT + NFBAS(I,0)*NTDIM*NZ
        IOFU = IOFU + NFORB(I,0)*NFORB(I,0)*NZ
        IOFV = IOFV + NFORB(I,0)*NFBAS(I,0)
   30 CONTINUE
C
C     Print section
C
      FAC = DM2*CVAL*CVAL
      CALL HEADER('Output from DELPOS',-1)
      IF(FREEPJ.OR.VEXTPJ) THEN
        WRITE(LUPRI,'(A)') '* All positronic solutions deleted !'
      ELSE
        IEIG = 0
        DO I = 1,NFSYM
          IOFF = IEIG + NFORB(I,1)
          NDIM = MIN(NFORB(I,1),NFORB(I,2))
          NDEL = NFORB(I,2) - NDIM
          WRITE(LUPRI,'(3X,A,I2,A,I5)') '*** Fermion corep ',I,
     &         ', Number of spinors deleted: ',NDEL
          IF(IPRGEN.GE.0) THEN
             IF(((DABS(EIG(IOFF+1)-FAC)).GT.1D-07).
     &            OR.((DABS(EIG(IOFF+NDEL)-FAC)).GT.1D-07)) THEN
             WRITE(LUPRI,'(/3X,A/)') 
     &            '*** Eigenvalues of deleted spinors:'
             NLST = 5
             IF(NDEL.LT.5) NLST = NDEL
             WRITE(LUPRI,'(6X,I5,2F18.8)')
     &            (J,EIG(IOFF+J),EIG(IOFF+J)-FAC,J=1,NLST)
             IF(NDEL.GT.15) THEN
                WRITE(LUPRI,'(15X,A)') ('.', J=1,3)
                WRITE(LUPRI,'(6X,I5,2F18.8)')
     &               (J,EIG(IOFF+J),EIG(IOFF+J)-FAC,J=NDEL-5,NDEL)
             ELSE
                WRITE(LUPRI,'(6X,I5,2F18.8)')
     &               (J,EIG(IOFF+J),EIG(IOFF+J)-FAC,J=NLST+1,NDEL)
             ENDIF
             ENDIF
          ENDIF
          IF(IPRGEN.GE.2) THEN
             WRITE(LUPRI,'(/3X,A/)') 
     &            '*** Eigenvalues of deleted spinors:'
             WRITE(LUPRI,'(6X,2F18.8)')
     &            (EIG(IOFF+J),EIG(IOFF+J)-FAC,J=1,NDEL)
          ENDIF
          IF(IPRGEN.GE.3) THEN
            CALL HEADER('DELPOS: Charge conjugation: ',-1)
            WRITE(LUPRI,'(5X,2A20)') '* Positrons:','* Electrons :'
            JJ = IEIG+NFORB(I,0)+1
            WRITE(LUPRI,'(I5,3F20.10)')
     &        (J,EIG(IEIG+J),EIG(JJ-J),(EIG(IEIG+J)+EIG(JJ-J)-FAC),
     &         J=1,NFORB(I,0))
          ENDIF
          IEIG = IEIG + NFORB(I,0)
        ENDDO
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
C  /* Deck gtovlp */
      SUBROUTINE GTOVLP(SMAT,WORK,LWORK)
C*****************************************************************************
C
C     Get overlap matrix from file AOPROPER
C
C     Written by T.Saue May 26 1996
C     Last revision: May 26 1996 - tsaue
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbgen.h>
#include <cbihr1.h>
#include <dcbbas.h>
      DIMENSION SMAT(N2BBASX),WORK(LWORK)
      CALL QENTER('GTOVLP')
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
      CALL MEMREL('GTOVLP',WORK,KWORK,KWORK,KFREE,LFREE)
C
C     Reindex to sorted basis
C     =======================
C
      CALL BUTOBS(SMAT,1,WORK,LWORK)
C
      CALL QEXIT('GTOVLP')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck heinz */
      SUBROUTINE HEINZ(WORK,LWORK)
C*****************************************************************************
C
C     Solve one-electron systems with finite nuclei and projecting
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
C     Solve one-electron systems with finite nuclei and projecting
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
     &              FMO,NORB(1),NORB(1),NBF,IQBF,
     &              TMAT,NFBAS(1,0),NORB(1),NZ,IPQTOQ(1,0),
     &              TMAT,NFBAS(1,0),NORB(1),NZ,IPQTOQ(1,0),
     &              WORK,LWORK,IPRGEN)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** HEINZ2: NBF.NE.NZ !'
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
     &              FMO,NESH(1),NESH(1),NBF,IQBF,
     &              TMAT(IOFF),NFBAS(1,0),NFORB(1,0),NZ,IPQTOQ(1,0),
     &              TMAT(IOFF),NFBAS(1,0),NFORB(1,0),NZ,IPQTOQ(1,0),
     &              WORK,LWORK,IPRGEN)
        IF(NBF.NE.NZ) WRITE(LUPRI,'(A)')
     &     '*** WARNING *** HEINZ2: NBF.NE.NZ !'
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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphlow */
      SUBROUTINE SPHLOW(LCOMP,SMAT,VMAT,LDM,NEFF,STOL,NDEL,IPRINT,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix with spherical transformation embedded
C
C     Written by T.Saue Sep 22 1995
C     Last revision Sep 22 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
      DIMENSION SMAT(*),VMAT(*),WORK(LWORK)
C
      CALL QENTER('SPHLOW')
#include <memint.h>
C     Memory allocation
      NDIM = NTBAS(LCOMP)
      NSPH = NSPHCM(LCOMP)
      N2SPH = NSPH*NSPH
      CALL MEMGET('REAL',KSSPH,N2SPH,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVSPH,N2SPH,WORK,KFREE,LFREE)
      CALL SPHLO1(LCOMP,SMAT,VMAT,NDIM,NEFF,LDM,STOL,NDEL,
     &            WORK(KSSPH),WORK(KVSPH),
     &            NSPH,IPRINT,WORK(KFREE),LFREE)
      CALL MEMREL('SPHLOW.SPHLO1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('SPHLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck sphlo1 */
      SUBROUTINE SPHLO1(LCOMP,SMAT,VMAT,NDIM,NEFF,LDM,STOL,NDEL,
     &                  SSPH,VSPH,NSPH,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix with spherical transformation embedded
C
C     Written by T.Saue Sep 22 1995
C     Last revision Sep 22 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION SMAT(LDM,NDIM) ,VMAT(LDM,NDIM),
     &          SSPH(NSPH,NSPH),VSPH(NSPH,NSPH),
     &          WORK(LWORK)
C
C     Transform overlap matrix to spherical basis
C
      CALL MTAOCS(SMAT,LDM,SSPH,NSPH,LCOMP,WORK,LWORK,IPRINT)
C
C     Generate spherical Lowdin matrix
C
      CALL LOWGEN(LCOMP,SSPH,VSPH,NSPH,NSPH,NEFF,STOL,IPRINT,
     &            WORK,LWORK)
      NDEL = NSPH - NEFF
C
C     Backtransform rows of Lowdin matrix
C      --> totally symmetric contributions are projected out
C
      CALL VTCSAO(VSPH,NSPH,VMAT,LDM,NEFF,LCOMP,IPRINT)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowgen */
      SUBROUTINE LOWGEN(LCOMP,SMAT,VMAT,LDM,NDIM,NEFF,STOL,IPRINT,
     &     WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix
C
C     Written by T.Saue Sep 11 1995
C     Last revision Sep 19 1995 - tsaue
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
      DIMENSION SMAT(LDM,NDIM),VMAT(LDM,NDIM),WORK(LWORK)
C
      CALL QENTER('LOWGEN')
#include <memint.h>
C     Memory allocation
      CALL MEMGET('REAL',KEIG ,NDIM  ,WORK,KFREE,LFREE)
      CALL LOWGE1(LCOMP,SMAT,VMAT,NDIM,NEFF,LDM,STOL,WORK(KEIG),
     &            IPRINT,WORK(KFREE),LFREE)
      CALL MEMREL('LOWGEN.LOWGE1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('LOWGEN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lowge1 */
      SUBROUTINE LOWGE1(LCOMP,SMAT,VMAT,NDIM,NEFF,LDM,STOL,EIG,IPRINT,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix
C
C     Written by T.Saue Sep 11 1995
C     Last revision Sep 19 1995
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D1 = 1.0D0,DM1 = -1.0D0)
C
C Used from COMMON blocks
C  DCBBAS: NTBAS,NFBAS
C  DHFCAL: STOL
C  DGROUP: NFSYM,NBSYM
C
#include <dgroup.h>
      DIMENSION SMAT(LDM,NDIM),VMAT(LDM,NDIM),EIG(NDIM),
     &          WORK(LWORK)
      NTOT = (NDIM-1)*LDM + NDIM
C
C     Diagonalize overlap matrix
C
C     IJOB   =  1: Eigenvalues and eigenvectors
C     IORDER = -1: Sort eigenvalues and eigenvectors in descending order
C     IPACK  =  0: Full overlap matrix given
C
      CALL RSJACO(LDM,NDIM,NDIM,SMAT,EIG,1,-1,0,VMAT,WORK,LWORK)
C      Avoid Householder; it mixes symmetries !!!
C      CALL RS(LDM,NDIM,SMAT,EIG,1,VMAT,BUF(1,1),BUF(1,2),IERR)
C
C     DLIMIT Should be set to ca. 1.0D-05 in final version!
C
      DLIMIT = 1.0D-02
      IF((IPRINT.GE.0).AND.(EIG(NDIM).LT.DLIMIT)) THEN
         CALL HEADER('Generating AO-MO Transformation',-1)
         WRITE(LUPRI,'(/3X,A,X,E6.1)') 
     &        '*** Eigenvalues of overlap matrix below', DLIMIT
         IF (LCOMP.EQ.1) THEN
            WRITE(LUPRI,'(3X,A/)') 
     &        '    Large component basis:'
         ELSEIF(LCOMP.EQ.2) THEN
            WRITE(LUPRI,'(3X,A/)') 
     &        '    Small component basis:'
         ELSE
            CALL QUIT('LOWGEN: Undefined LCOMP')
         ENDIF
         DO I = NDIM,1,-1
            IF (EIG(I).LT.DLIMIT) WRITE(LUPRI,'(I5,F16.9)') I,EIG(I)
         ENDDO
      ENDIF
      IF(IPRINT.GE.4) THEN
        CALL HEADER('LOWGEN: All Eigenvalues of overlap matrix:',-1)
        WRITE(LUPRI,'(I5,F16.9)') (I,EIG(I),I=1,NDIM)
      ENDIF
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Eigenvectors of overlap matrix:',-1)
        CALL OUTPUT(VMAT,1,NDIM,1,NDIM,LDM,NDIM,1,LUPRI)
      ENDIF
      NEFF = 0
      DO J = 1,NDIM
        IF(EIG(J).LE.STOL) GOTO 10
        NEFF = NEFF + 1
        FAC   = D1/SQRT(EIG(J))
        CALL DSCAL(NDIM,FAC,VMAT(1,J),1)
      ENDDO
 10   CONTINUE
C
C     Zero the columns that are to be deleted
C
      NDEL = NDIM - NEFF
      IF(NDEL.GE.1) THEN
        CALL HEADER('WARNING from LOWGEN',-1)
        IF(IPRINT.GE.3) THEN
           WRITE(LUPRI,'(I5,A,1P,E9.3)')
     &     NDEL,' eigenvalues of overlap matrix below threshold ',
     &     STOL
          WRITE(LUPRI,'(3X,A)')
     &      'Eigenvectors of AO-overlap matrix deleted:'
          CALL OUTPUT(VMAT,(NEFF+1),NDIM,1,NDIM,LDM,NDIM,1,LUPRI)
        ENDIF
        DO J = 1,NDEL
          CALL DZERO(VMAT(1,NEFF+J),NDIM)
        ENDDO
      ENDIF
      IF(IPRINT.GE.6) THEN
        CALL HEADER('LOWGEN: Lowdin matrix:',-1)
        CALL OUTPUT(VMAT,1,NDIM,1,NEFF,LDM,NDIM,1,LUPRI)
      ENDIF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rkblow */
      SUBROUTINE RKBLOW(LCOMP,SMAT,VMAT,LDM,NEFF,STOL,NDEL,IPRINT,
     &                  WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix for S components with RKB transformation embedded
C
C                          J.K. Laerdahl May 08 1996
C                          Modification of T. Saue's SPHLOW
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
C
#include <dcbbas.h>
#include <maxaqn.h>
#include <mxorb.h>
#include <dcbrkb.h>
      DIMENSION SMAT(*),VMAT(*),WORK(LWORK)
C
      CALL QENTER('RKBLOW')
#include <memint.h>
C     Memory allocation
      NDIM = NTBAS(LCOMP)
      NRKB = NRKBSF
      N2RKB = NRKB*NRKB
      CALL MEMGET('REAL',K1RKB,N2RKB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',K2RKB,N2RKB,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',NTRRKB,NRKB*NDIM,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',NTRRKB2,NRKB*NDIM,WORK,KFREE,LFREE)
      CALL RKBLO1(LCOMP,SMAT,VMAT,NDIM,NEFF,LDM,STOL,NDEL,
     &            WORK(K1RKB),WORK(K2RKB),
     &            WORK(NTRRKB),WORK(NTRRKB2),
     &            NRKB,IPRINT,WORK(KFREE),LFREE)
      CALL MEMREL('RKBLOW.RKBLO1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('RKBLOW')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rkblo1 */
      SUBROUTINE RKBLO1(LCOMP,SMAT,VMAT,NDIM,NEFF,LDM,STOL,NDEL,
     &                  SRKB,VRKB,TRNMA1,TRNMA2,NRKB,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Generate Lowdin matrix with RKB transformation embedded
C     Written by J.K. Laerdahl 1996
C
C*****************************************************************************
#include <implicit.h>
#include <priunit.h>
      DIMENSION SMAT(LDM,NDIM) ,VMAT(LDM,NDIM),
     &     SRKB(NRKB,NRKB),VRKB(NRKB,NRKB),
     &     TRNMA1(NRKB,NDIM),TRNMA2(NDIM,NRKB),
     &     WORK(LWORK)
      PARAMETER(D1=1.0D0, D0=0.0D0)
C     
C     Construct Transformation matrices
C     
      CALL RKBTRN(TRNMA1,TRNMA2,NRKB,NDIM,LCOMP,IPRINT)
C
C     Transform overlap matrix to RKB basis
C     
#include <memint.h>
      NWRK = NRKB*NDIM
      CALL MEMGET('REAL',KW1,NWRK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KW2,NWRK,WORK,KFREE,LFREE)
      CALL DGEMM('N','N',NRKB,NDIM,NDIM,D1,TRNMA1,NRKB,
     &               SMAT,LDM,D0,WORK(KW1),NRKB)
C
C     Transpose and transform columns
C
      CALL MTRSP(NRKB,NDIM,WORK(KW1),NRKB,WORK(KW2),NDIM)
      CALL DGEMM('N','N',NRKB,NRKB,NDIM,D1,TRNMA1,NRKB,
     &               WORK(KW2),NDIM,D0,WORK(KW1),NRKB)
C
C     Backtranspose to obtain final matrix
C
      CALL MTRSP(NRKB,NRKB,WORK(KW1),NRKB,SRKB,NRKB)
C
      IF (IPRINT.GT.4) THEN
         CALL HEADER('SRKB: Transformed matrix:',-1)
         CALL OUTPUT(SRKB,1,NRKB,1,NRKB,NRKB,NRKB,1,LUPRI)
      ENDIF
C
C     Generate RKB Lowdin matrix
C      
      CALL LOWGEN(LCOMP,SRKB,VRKB,NRKB,NRKB,NEFF,STOL,IPRINT,
     &     WORK,LWORK)
      NDEL = NRKB - NEFF
C     
C     Backtransform rows of Lowdin matrix
C     
      CALL DGEMM('N','N',NDIM,NRKB,NRKB,D1,TRNMA2,NDIM,
     &     VRKB,NRKB,D0,VMAT,LDM)
C     
      CALL MEMREL('RKBLO1.RKBTRN',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rkbtrn */
      SUBROUTINE RKBTRN(TRNMA1,TRNMA2,ND1,NDIM,LCOMP,IPRINT)
#include <implicit.h>
#include <priunit.h>
      PARAMETER(D6=6.0D0,D4=4.0D0,D1=1.0D0)
#include <maxaqn.h>
#include <mxorb.h>
#include <mxcent.h>
#include <nuclei.h>
#include <shells.h>
#include <dcbbas.h>
#include <dcbrkb.h>
      LOGICAL FOUND
      DIMENSION TRNMA1(ND1,NDIM),TRNMA2(NDIM,ND1)
C
      CALL DZERO(TRNMA1,ND1*NDIM)
C     
      IF (IPRINT.GE.5) CALL HEADER('RKBTRN: Building TRNMA1',1)
      IOFFI = 0
      IOFFJ = 0
C
C     Loop over Shells of S comp functions
C
      DO 210 ISHELL = NLRGSH+1,KMAX
         L     = NHKT(ISHELL)
         NLM   = 2*(L-1) + 1
         NRKB  = NLM
         NXYZ  = L*(L+1)/2
         NDEG   = NUCDEG(NCENT(ISHELL))
         IF (L.GE.3) THEN
            NRKB = NRKB + 2*(L-2)-1
            NXYZM2 = (L-2)*(L-1)/2
            NLMM2  = 2*(L-2)-1
C     
C     Scaling factor for S comp L-2 partner 
C     Independent of Alpha because of
C     normalization constant
C
            DLM2SC = (-D4*DFLOAT(L)+D6)
         ENDIF
         ICNTL = ISTRT(ISHELL) + 1
         FOUND = .FALSE.
         IXX = 0
 50      IXX = IXX + 1
         IF (IXX.GT.IRKBLE) GOTO 55
         IF (ICNTL.EQ.IRKBLS(IXX)) THEN
            FOUND = .TRUE.
            GOTO 55
         ENDIF
         GOTO 50
 55      CONTINUE
C
C     Loop over degenerate centers
C
         DO IDEG = 1,NDEG
            IF (FOUND) THEN
               NRKB = 0
               GOTO 111
            ENDIF
            DO I = 1, NRKB
               DO J = 1, NXYZ
                  TRNMA1(IOFFI+I,IOFFJ+J) = CRKB(ISRKB(L)
     &                 +NRKB*(J-1)+I-1)
                  IF(IPRINT.GE.10)
     &                 WRITE(LUPRI,'(A,2I5,F15.8)') '#', I,J, 
     &                 TRNMA1(IOFFI+I,IOFFJ+J)
               ENDDO
C     
C     Contract with L-2 partner
C     
               IF (I.GT.NLM) THEN
                  IDEGP = (IDEG-1)*KHKT(IRKBPT(ISHELL,2))
                  LM2 = I - NLM
                  DO IM2 = 1, NXYZM2
                     ELEMNT = CSP(ISPADR(L-2)+(IM2-1)*NLMM2+LM2-1)
                     IF(IPRINT.GE.10) THEN
                        WRITE(LUPRI,'(A,2I5,F15.8,7I6)') '###', 
     &                       I,IM2,ELEMNT, ISPADR(L-2),NXYZM2,LM2, 
     &                       IOFFI+I,IOFFJ+J,IRKBPT(ISHELL,1),IDEGP
                     ENDIF
                     IPUT=IRKBPT(ISHELL,1)+IDEGP+IM2-NTBAS(1)-1 
                     TRNMA1(IOFFI+I,IPUT) = DLM2SC*ELEMNT
                  ENDDO
               ENDIF
            ENDDO
 111        CONTINUE
            IOFFI = IOFFI + NRKB
 200        IOFFJ = IOFFJ + NXYZ
         ENDDO
 210     CONTINUE
C     
C     Print section
C
      IF(IPRINT.GE.5) THEN
         CALL HEADER('TRNMA1:',-1)
         CALL OUTPUT(TRNMA1,1,ND1,1,NDIM,ND1,NDIM,1,LUPRI)
      ENDIF
C
C      It is NOT necessary to Normalize matrix
C
      CALL MTRSP(ND1,NDIM,TRNMA1,ND1,TRNMA2,NDIM) 
      IF(IPRINT.GE.6) THEN
         CALL HEADER('TRNMA2:',-1)
         CALL OUTPUT(TRNMA2,1,NDIM,1,ND1,NDIM,ND1,1,LUPRI)
      ENDIF
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
  1000 CONTINUE
C
C     Not able to read DOINT information from COMB
C
c     There is no goto 1000 is this routine./Feb 9 1998 - jth
c      WRITE(LUPRI,'(A,A)') 'TSCORR: Not able to read COMB =',COMB
c      CALL QUIT('TSCORR: Not able to read COMB')
      END
C
