*
* 2 electron contributions with mk=1
*
* The generic form is
*
* sigma(Ia,Ib) = sum(ijkl) <Ia!a+ia!Ja> <Ib!a+jb a kb a lb!Jb>
*                (ij!kl) C(Ja,Jb)
* The current code goes as
*
* Loop over batches of Ib
*
* set up C(Ka,I,Jb) = sum(Ja) <Ia!a+ ia!Ja> C(Ja,Jb>
*
* Set up excitations <Ib!a+jb !Kb>, <Kb! a kb !Mb> <Mb! a lb !Jb>
*
      SUBROUTINE DMK1TWO(NIA,IASM,IATP,IAGRP,IAOC,
     &                   NIB,IBSM,IBTP,IBGRP,IBOC,
     &                   NJA,JASM,JATP,JAGRP,JAOC,
     &                   NJB,JBSM,JBTP,JBGRP,JBOC,
     &          I_AC,JKL_AC,NGAS,SB,CB,
     &          NALPTS,NBEPTS,MAXK,
     &          SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &          XINT,NSMOB,NSMST,IUSEAB,
     &          CJRES,SIRES,SCLFAC,NTEST,
     &          IUSE_PH,IPHGAS,IBU,JBU,KBU,LBU,
     &          XIJILS)
*
* Note on symmetry handling : I assume that double group symmetry is
* in use, and use info in DBGSYM
*
* =====
* Input
* =====
*
* NIA(NIB,NJA,NJB)         : Number of IA strings
* IASM(IBSM,JASM,JBSM)     : Symmetry of IA strings
* IATP(IBTP,JATP,JBTP)     : Supergroup of strings IA (relative to start of type!)
* IAGRP(IBGRP,JAGRP,JBGRP) : Type (number of electrons) of IA strings
* IAOCC(IBGRP,JAOCC,JBOCC) : Occupation of  each GAS space for IA strings
*
* I_AC = 1, JKL_AC(1-3) = 2,2,1 => a ia a+jb a+kb alb
* I_AC = 2, JKL_AC(1-3) = 2,1,1 => a+ia a+jb a kb alb
*
* NALPTS : Number of alpha spinors per sym and type
* NALPTS : Number of beta  spinors per sym and type
* NGAS : Number of GAS orbital spaces
* CB   : Input C block
* ADSXA : sym of a+, a+a => sym of a
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX : Number of symmetries of orbitals,strings,
*       single excitations
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* ======
* Output
* ======
* SB : updated sigma block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* XINT  : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1997
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. Symmetry info
#include "symm.inc"
*. General input
      INTEGER NALPTS(MXPNGAS,*),NBEPTS(MXPNGAS,*)
*
*.Input
      DIMENSION CB(*)
      INTEGER I_AC(1),JKL_AC(3)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
*
*.Local arrays
      PARAMETER (MXNX = 256)
      DIMENSION ITP(256),JKLTP(256,3)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*
      INTEGER MBOC(20),KBOC(20),LBOC(20), KBSPGP(20),MBSPGP(20)

*. Arrays for reorganization
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C     CHARACTER*6 DBGSYM
C     COMMON/DBGSYM/IDBGMULT(8,8),IPGMULT(4,4),
C    &              NDBGIRR,IMOSF_SP(8,2),ISPSF_MO(8,2),
C    &              IADJSYM(8),INVELM(8)
* Contents of DBGSYM
*    IDBGMULT : Group multiplication table  in double group
*    IPGMULT  : Group multiplication table  in point  group
*    NDBGIRR  : Number of irreps in double group
*    IMOSF_SP : Translation table from MO symmetry to spinor symmetry
*               IMOSF_SP(IOBSM,IAB) Gives the spinor symmetry
*               of an orbital of symmetry IOBSM times an alpha (IAB=1)
*               or and beta spin-function (IAB = 2 )
*    ISPSF_MO : Translation table from spinor symmetry to MO symmetry
*               ISPSF(ISPSM,IAB) gives the symmetry of an MO that
*               combined with an alpha(IAB=1) or beta(IAB=2)
*               spin function give a spinor with symmetry ISPSM
*
*   IADJSYM   : IADJSYM(ISPSYM) is the adjoint symmetry of ISPSYM
*   INVELM    : Inverse elements
*
* The irreps are organized in the usual way with the boson irreps first
*
      CALL QENTER('DMK1TWO')
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ==================== '
        WRITE(6,*) '  SIGMA_MK1 speaking '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' Occupation of groups (IA,IB,JA,JB) '
        WRITE(6,*)
        CALL IWRTMA(IAOCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        CALL IWRTMA(IBOCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        CALL IWRTMA(JAOCC,1,NGAS,1,NGAS)
        WRITE(6,*)
        CALL IWRTMA(JBOCC,1,NGAS,1,NGAS)
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup ( Group number for each gas space )
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,JAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,JBGRP,JBSPGP)
*. Symmetry of allowed excitations
      IISM = IDBGMULT(IADJSYM(IASM),JASM)
      IF(I_AC(1).EQ.1) THEN
        ISM = IADJSYM(IISM)
      ELSE
        ISM = IISM
      END IF
      JKLSM = IDBGMULT(IADJSYM(IBSM),JBSM)
      IF(ISM.EQ.0.OR.JKLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM ISM ',IASM,JASM,ISM
        write(6,*) ' IBSM JBSM JKLSM ',IBSM,JBSM,JKLSM
      END IF
*.Types of Excitations  that connects the two strings
C     GXFTP(NOP,IAC,ILOCC,IROCC,NGAS,NX,IX,MXNX)
      CALL GXFTP(1,I_AC,IAOC,JAOC,NGAS,NIX,ITP,MXNX)
      CALL GXFTP(3,JKL_AC,IBOC,JBOC,NGAS,NJKLX,JKLTP,MXNX)
      IF(NIX.EQ.0.OR.NJKLX.EQ.0) GOTO 9999
      DO 2001 IX = 1, NIX
*
        ITYP = ITP(IX)
        KAFRST = 1
        NI = NALPTS(ITYP,ISM)
*
*
*. Generate creation- or annihilation- mappings <Ia!a+ ia! Ja>
*. Ordered according to Ia so we set up tables <Ja! a ia !Ia>
        IF(I_AC(1).EQ.1) THEN
          IAC_CONJ = 2
        ELSE
          IAC_CONJ = 1
        END IF
*
        CALL ADAST_GAS(ISM,ITYP,NGAS,JASPGP,JASM,
     &       I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFAC,IAC_CONJ)
*. No compression at the moment
        NKAEFF = NKASTR
*. Loop over batches of KA strings (KA strings = IA strings)
        NKABTC = NKAEFF/MAXK
        IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
        DO 1801 IKABTC = 1, NKABTC
           KABOT = (IKABTC-1)*MAXK + 1
           KATOP = MIN(KABOT+MAXK-1,NKAEFF)
           LKABTC = KATOP-KABOT+1
*. Obtain C(ka,i,JB) for Ka in batch
            DO II = 1, NI
              CALL GET_CKAJJB(CB,NI,NJA,CJRES,LKABTC,NJB,
     &             II,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC)
            FACS = 1.0D0
*
            DO 2000 JKL = 1, NJKLX
              JTYP = JKLTP(JKL,1)
              KTYP = JKLTP(JKL,2)
              LTYP = JKLTP(JKL,3)
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' JTYP, KTYP, LTYP', JTYP,KTYP, LTYP
              END IF
*
              DO 1940 JSM = 1, NSMOB
              DO 1930 KSM = 1, NSMOB
*
                IFIRST = 1
*
                IF(JKL_AC(1).EQ.1) THEN
                  JJSM = IADJSYM(JSM)
                ELSE
                  JJSM = JSM
                END IF
*
                IF(JKL_AC(2).EQ.1) THEN
                  KKSM = IADJSYM(KSM)
                ELSE
                  KKSM = KSM
                END IF
*
                JKSM = IDBGMULT(JJSM,KKSM)
                LLSM = IDBGMULT(INVELM(JKSM),JKLSM)
                IF(JKL_AC(3).EQ.1) THEN
                  LSM = IADJSYM(LLSM)
                ELSE
                  LSM = LLSM
                END IF
*
                IF(LSM.EQ.0) GOTO 1930
                NJ = NBEPTS(JTYP,JSM)
                NK = NBEPTS(KTYP,KSM)
                NL = NBEPTS(LTYP,LSM)
                IF(NJ.EQ.0.OR.NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
* we will set the following mappings up
*   <KB! a j (conj) !IB>
*   <MB! a k (conj) !KB>
*   <JB! a l (conj) !MB>, drivenm by the strings to the right
*. Find symmetry, group and type of  KB and MB
                KBSM = IDBGMULT(JJSM,IBSM)
                MBSM = IDBGMULT(KKSM,KBSM)
*
                CALL ICOPVE(IBOC,KBOC,NGAS)
                IF(JKL_AC(1).EQ.1) THEN
                  KBOC(JTYP) = KBOC(JTYP) + 1
                  KBTP = IBTP + 2
                  JCONJ = 2
                ELSE
                  KBOC(JTYP) = KBOC(JTYP) - 1
                  KBTP = IBTP - 2
                  JCONJ = 1
                END IF
*
                CALL ICOPVE(KBOC,MBOC,NGAS)
                IF(JKL_AC(2).EQ.1) THEN
                  MBOC(JTYP) = KBOC(JTYP) + 1
                  MBTP = KBTP + 2
                  KCONJ = 2
                ELSE
                  MBOC(JTYP) = KBOC(JTYP) - 1
                  MBTP = KBTP - 2
                  KCONJ = 1
                END IF
*
                IF(JKL_AC(3).EQ.1) THEN
                  LCONJ = 2
                ELSE
                  LCONJ = 1
                END IF
*   <KB! a j (conj) !IB>
                CALL ADAST_GAS(JSM,JTYP,NGAS,KBSPGP,KBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,JKL_AC(1))
                IF(NKBSTR.EQ.0) GOTO 1930
*   <MB! a k (conj) !KB>
                CALL ADAST_GAS(KSM,KTYP,NGAS,MBSPGP,MBSM,
     &               I3,XI3S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,JKL_AC(2))
                IF(NKBSTR.EQ.0) GOTO 1930
*   <JB! a l (conj) !MB>
                CALL ADAST_GAS(LSM,LTYP,NGAS,JBSPGP,JBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,JKL_AC(3))
* Fetch Integrals as (i j |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                CALL GETINT_DBG(XINT,ISM,ITYP,IBU,JSM,JTYP,JBU,
     &                           KSM,KTYP,KBU,LSM,LTYP,LBU,
     &                           IRI,IXCHNG,0,0,ICOUL)
*
* S(Ka,Ib) = sum(i,j,k,l,Jb)<Ib!a jb a kb a lb!Jb>C(Ka,i,Jb)*(ij!kl)
*
                CALL SKICKJ_MK1(SB,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,NKBSTR,NKBSTR,XINT,NI,NJ,NK,NL,
     &               NKBSTR,NKBSTR,I2,XI2S,I3,XI3S,I4,XI4S,IKORD,
     &               FACS,MXTSOB,XIJILS)
C               SKICKJ_MK1(SKII,CKJJ,NKA,NIB,NJB,NKB,NMB,XIJKL,
C    &                  NI,NJ,NK,NL,MAXI,MAXK,MAXM
C    &                  KBIB,XKBIB,MBKB,XMBKB,JBMB,XJBMB,IKORD,
C    &                  FACS,MXTSOB,XIJILS)

 1930         CONTINUE
*             ^ End of loop over KSM
 1940       CONTINUE
*          .^ End of loop over ISM
 2000       CONTINUE
*           ^ End of loop over KLTYP
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL QEXIT('DMK1TWO')
      RETURN
      END
