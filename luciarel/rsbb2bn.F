      SUBROUTINE RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                   JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                   IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                   SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                   NOBPTS,       MAXK,
     &                   SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                   XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                   CJRES,SIRES,SCLFAC,NTEST,
     &                   NSEL2E,ISEL2E,IUSE_PH,IPHGAS)
*
* Combined alpha-beta double excitation
* contribution from given C block to given S block
*. If IUSAB only half the terms are constructed
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha  strings in sigma
* IBSM,IBTP : Symmetry and type of beta   strings in sigma
* JASM,JATP : Symmetry and type of alpha  strings in C
* JBSM,JBTP : Symmetry and type of beta   strings in C
* NIA,NIB : Number of alpha-(beta-) strings in sigma
* NJA,NJB : Number of alpha-(beta-) strings in C
* IAGRP : String group of alpha strings
* IBGRP : String group of beta strings
* IAEL1(3) : Number of electrons in RAS1(3) for alpha strings in sigma
* IBEL1(3) : Number of electrons in RAS1(3) for beta  strings in sigma
* JAEL1(3) : Number of electrons in RAS1(3) for alpha strings in C
* JBEL1(3) : Number of electrons in RAS1(3) for beta  strings in C
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
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + CJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks
* October 96 : New routines for accessing annihilation information
*             Cleaned and shaved, only IROUTE = 3 option active
* October   97 : allowing for N-1/N+1 switch
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*)
*
      INTEGER ISEL2E(*)
*.Input
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION XINT(*)
      DIMENSION CJRES(*),SIRES(*)
*
      PARAMETER(MXTSOB = 40)
      DIMENSION H(MXTSOB*MXTSOB)
*.Local arrays
      DIMENSION ITP(20),JTP(20),KTP(20),LTP(20)
      DIMENSION IOP_TYP(2),IOP_AC(2),IOP_REO(2)
*
      DIMENSION IJ_TYP(2),IJ_DIM(2),IJ_REO(2),IJ_AC(2),IJ_SYM(2)
      DIMENSION KL_TYP(2),KL_DIM(2),KL_REO(2),KL_AC(2),KL_SYM(2)
*
      DIMENSION IASPGP(20),IBSPGP(20),JASPGP(20),JBSPGP(20)
*. Arrays for reorganization
      DIMENSION NADDEL(6),IADDEL(4,6),IADOP(4,6),ADSIGN(6)
C    &          SIGNREO,NADOP,NADDEL,IADDEL,ADSIGN)
*
      CALL QENTER('RS2B ')
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' RSBB2BN speaking '
        WRITE(6,*) ' =============== '
      END IF
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Groups defining each supergroup
      CALL GET_SPGP_INF(IATP,IAGRP,IASPGP)
      CALL GET_SPGP_INF(JATP,IAGRP,JASPGP)
      CALL GET_SPGP_INF(IBTP,IBGRP,IBSPGP)
      CALL GET_SPGP_INF(JBTP,IBGRP,JBSPGP)
*
*. Symmetry of allowed excitations
      IJSM = STSTSX(IASM,JASM)
      KLSM = STSTSX(IBSM,JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP2_GAS(NKLTYP,KTP,LTP,NGAS,IBOC,JBOC,IPHGAS)
      CALL SXTYP2_GAS(NIJTYP,ITP,JTP,NGAS,IAOC,JAOC,IPHGAS)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
*
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Should N-1 or N+1 projection be used for alpha strings
          IJ_TYP(1) = ITYP
          IJ_TYP(2) = JTYP
          IJ_AC(1)  = 2
          IJ_AC(2) =  1
          NOP = 2
          IF(IUSE_PH.EQ.1) THEN
            CALL ALG_ROUTERX(IAOC,JAOC,NOP,IJ_TYP,IJ_AC,IJ_REO,
     &           SIGNIJ)
          ELSE
*. Enforced a+ a
            IJ_REO(1) = 1
            IJ_REO(2) = 2
            SIGNIJ = 1.0D0
          END IF
*. Two choices here :
*  1 : <Ia!a+ ia!Ka><Ja!a+ ja!Ka> ( good old creation mapping)
*  2 :-<Ia!a  ja!Ka><Ja!a  ia!Ka>  + delta(i,j)
C?        WRITE(6,*) ' RSBB2BN : IOP_REO : ', (IOP_REO(II),II=1,2)
          IF(IJ_REO(1).EQ.1.AND.IJ_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
            IJAC = 2
            SIGNIJ2 = SCLFAC
*
            IJ_DIM(1) = NI
            IJ_DIM(2) = NJ
            IJ_SYM(1) = ISM
            IJ_SYM(2) = JSM
            IJ_TYP(1) = ITYP
            IJ_TYP(2) = JTYP
*
            NOP1   = NI
            IOP1SM = ISM
            IOP1TP = ITYP
            NOP2   = NJ
            IOP2SM = JSM
            IOP2TP = JTYP
          ELSE
*. Terra Nova, annihilation map
            IJAC = 1
            SIGNIJ2 = -SCLFAC
*
            IJ_DIM(1) = NJ
            IJ_DIM(2) = NI
            IJ_SYM(1) = JSM
            IJ_SYM(2) = ISM
            IJ_TYP(1) = JTYP
            IJ_TYP(2) = ITYP
*
            NOP1   = NJ
            IOP1SM = JSM
            IOP1TP = JTYP
            NOP2   = NI
            IOP2SM = ISM
            IOP2TP = ITYP
          END IF
*
*. Generate creation- or annihilation- mappings for all Ka strings
*
*. For operator connecting to |Ka> and |Ja> i.e. operator 2
          CALL ADAST_GAS(IJ_SYM(2),IJ_TYP(2),NGAS,JASPGP,JASM,
     &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SIGNIJ2,IJAC)
C         CALL ADAST_GAS(JSM,JTYP,JATP,JASM,IAGRP,
C    &         I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFACS,IJ_AC)
*. For operator connecting |Ka> and |Ia>, i.e. operator 1
          CALL ADAST_GAS(IJ_SYM(1),IJ_TYP(1),NGAS,IASPGP,IASM,
     &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJAC)
C         CALL ADAST_GAS(ISM,ITYP,NGAS,IASPGP,IASM,
C    &         I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE,IJ_AC)
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
              CALL COMPRS2LST(I1,XI1S,IJ_DIM(2),I3,XI3S,IJ_DIM(1),
     &                        NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
*
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, IJ_DIM(2)
              CALL GET_CKAJJB(CB,IJ_DIM(2),NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            MXACJ=MAX(MXACJ,NIB*LKABTC*IJ_DIM(1),NJB*LKABTC*IJ_DIM(2))
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*IJ_DIM(1))
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
              IF(NTEST.GE.100) THEN
                WRITE(6,*) ' KTYP, LTYP', KTYP, LTYP
              END IF
*. Should this group of excitations be included
              IF(NSEL2E.NE.0) THEN
               IAMOKAY=0
               IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP)THEN
                 DO JSEL2E = 1, NSEL2E
                   IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY = 1
                 END DO
               END IF
               IF(IAMOKAY.EQ.0) GOTO 2000
              END IF
*
              KL_TYP(1) = KTYP
              KL_TYP(2) = LTYP
              KL_AC(1)  = 2
              KL_AC(2) =  1
              NOP = 2
              IF(IUSE_PH.EQ.1) THEN
                CALL ALG_ROUTERX(IBOC,JBOC,NOP,KL_TYP,KL_AC,KL_REO,
     &               SIGNKL)
              ELSE
*. Enforced a+ a
                KL_REO(1) = 1
                KL_REO(2) = 2
                SIGNKL = 1.0D0
              END IF
*
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*
                IF(KL_REO(1).EQ.1.AND.KL_REO(2).EQ.2) THEN
*. Business as usual i.e. creation map
                  KLAC = 2
                  KL_DIM(1) = NK
                  KL_DIM(2) = NL
                  KL_SYM(1) = KSM
                  KL_SYM(2) = LSM
                  KL_TYP(1) = KTYP
                  KL_TYP(2) = LTYP
                ELSE
*. Terra Nova, annihilation map
                  KLAC = 1
                  KL_DIM(1) = NL
                  KL_DIM(2) = NK
                  KL_SYM(1) = LSM
                  KL_SYM(2) = KSM
                  KL_TYP(1) = LTYP
                  KL_TYP(2) = KTYP
                END IF
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP)
     &          IKORD = 1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
*. currently we are using creation mappings for kl
*. (Modify to use ADAST later )
                CALL ADAST_GAS(KL_SYM(2),KL_TYP(2),NGAS,JBSPGP,JBSM,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,SIGNKL,KLAC)
C               CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
C    &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADAST_GAS(KL_SYM(1),KL_TYP(1),NGAS,IBSPGP,IBSM,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE,KLAC)
C               CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
C    &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (iop2 iop1 |  k l )
*
                IXCHNG = 0
                ICOUL = 1
                CALL LGETINT(XINT,IJ_TYP(2),IJ_SYM(2),
     &               IJ_TYP(1),IJ_SYM(1),
     &               KL_TYP(1),KL_SYM(1),KL_TYP(2),KL_SYM(2),IXCHNG,
     &               0,0,ICOUL)
C               CALL LGETINT(XINT,JTYP,JSM,ITYP,ISM,KTYP,KSM,
C    &                     LTYP,LSM,IXCHNG,0,0,ICOUL)
*
* S(Ka,j,Ib) = sum(k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                IROUTE = 3
                CALL SKICKJ(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,IJ_DIM(1),IJ_DIM(2),
     &               KL_DIM(1),KL_DIM(2),
     &               NKBSTR,I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE,MXTSOB)
 1930         CONTINUE
*             ^ End of loop over KSM
 2000       CONTINUE
*           ^ End of loop over KLTYP
*
*. Scatter out from s(Ka,Ib,i)
*
            DO II = 1, IJ_DIM(1)
              CALL ADD_SKAIIB(SB,IJ_DIM(1),NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*.        ^End of loop over partitioning of alpha strings
 1940   CONTINUE
*       ^ End of loop over ISM
 2001 CONTINUE
*     ^ End of loop over IJTYP
*
 9999 CONTINUE
*
*
      CALL QEXIT('RS2B ')
      RETURN
      END
