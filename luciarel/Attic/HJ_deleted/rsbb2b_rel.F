      SUBROUTINE RSBB2B_REL(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                  JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,
     &                  ADSXA,STSTSX,
     &                  NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &                  ITSOB,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                  XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,IUSEAB,
     &                  ICJKAIB,CJRES,SIRES,S2,SCLFAC,NTESTG,
     &                  NSEL2E,ISEL2E,IRI,SCR)
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
* ICJKAIB =1 =>  construct C(Ka,Jb,j) and S(Ka,Ib,i) as intermediate
*                 matrices in order to reduce overhead
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
* H : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
* Feb 92 : Loops restructured ; Generation of I2,XI2S moved outside
* October 1993 : IUSEAB added
* January 1994 : Loop restructured + ICJKAIB introduced
* February 1994 : Fetching and adding to transposed blocks
*
*. October 96 : New routines for accessing annihilation information
*.              Cleaned and shaved, only IROUTE = 3 option active
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),
     &        NOBPTS2(MXPNGAS,*),IOBPTS2(MXPNGAS,*),
     &        ITSOB(*)
*
      INTEGER ISEL2E(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
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
      DIMENSION S2(*)
      dimension SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS),KTP(MXPNGAS),LTP(MXPNGAS)
      INTEGER IREO(4),IILEN(4)
*
CTF
* Added these nice things ...
#include "symm.inc"
*
      CALL QENTER('RS2B ')
      NTESTL = 0
      NTEST = MAX(NTESTL,NTESTG)
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' RSBB2B speaking '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        WRITE(6,*) ' Initial sigma block '
        CALL WRTMAT(SB,NIB,NIA,NIB,NIA)

      END IF

      IROUTE = 3
* IROUTE = 1 : Normal (i.e. old) route,
* IROUTE = 2 : New route with j first
* IROUTE = 3 : C(Ka,j,Jb)
*
*
*. Symmetry of allowed excitations
      IJSM = IDBGMULT(IADJSYM(IASM),JASM)
      KLSM = IDBGMULT(IADJSYM(IBSM),JBSM)
      IF(IJSM.EQ.0.OR.KLSM.EQ.0) GOTO 9999
      IF(NTEST.GE.600) THEN
        write(6,*) ' IASM JASM IJSM ',IASM,JASM,IJSM
        write(6,*) ' IBSM JBSM KLSM ',IBSM,JBSM,KLSM
      END IF
*.Types of SX that connects the two strings
      CALL SXTYP_GAS_REL(NKLTYP,KTP,LTP,NGAS,MXPNGAS,IBOC,JBOC)
      CALL SXTYP_GAS_REL(NIJTYP,ITP,JTP,NGAS,MXPNGAS,IAOC,JAOC)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
*
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = IDBGMULT(ISM,IADJSYM(IJSM))
          IF(NTESTL.GE.2000)
     &    WRITE(6,*) ' Loop 1940, ISM, KSM ', ISM,JSM
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all K strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS_REL(JSM,JTYP,JATP,JASM,IAGRP,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE = 1.0D0
          CALL ADSTN_GAS_REL(ISM,ITYP,IATP,IASM,IAGRP,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,ONE)
*. Compress list to common nonvanishing elements
          IDOCOMP = 1
          IF(IDOCOMP.EQ.1) THEN
C             COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
              CALL COMPRS2LST(I1,XI1S,NJ,I3,XI3S,NI,NKASTR,NKAEFF)
          ELSE
              NKAEFF = NKASTR
          END IF

*. Loop over batches of KA strings
          NKABTC = NKAEFF/MAXK
          IF(NKABTC*MAXK.LT.NKAEFF) NKABTC = NKABTC + 1
C         IF(NKABTC.GT.1) THEN
C           WRITE(6,*) ' RSBB2B : NKABTC .GT. 1 '
C           WRITE(6,*) ' I am not prepared for this '
C           STOP 'RSBB2B :  NKABTC .GT. 1'
C         END IF
          DO 1801 IKABTC = 1, NKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*
            ZERO = 0.0D0
            CALL SETVEC(SIRES,ZERO,NIB*LKABTC*NI)
            FACS = 1.0D0
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
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
              DO 1930 KSM = 1, NSMOB
                IFIRST = 1
                LSM = IDBGMULT(KSM,IADJSYM(KLSM))
                IF(NTEST.GE.2000)
     &          WRITE(6,*) ' Loop 1930, KSM, LSM = ',KSM,LSM
                IF(LSM.EQ.0) GOTO 1930
                KOFF = IOBPTS2(KTYP,KSM)
                LOFF = IOBPTS2(LTYP,LSM)
                NK = NOBPTS2(KTYP,KSM)
                NL = NOBPTS2(LTYP,LSM)
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
                CALL ADSTN_GAS_REL(LSM,LTYP,JBTP,JBSM,IBGRP,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS_REL(KSM,KTYP,IBTP,IBSM,IBGRP,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*
* Fetch Integrals as (j i k l )
*
*. We want the integrals out in the order (J I ! K L).
*. However, for complex integrals we cannot just switch the
* indeces around, so obtain the correct integrals and
* reorder subsequently
*
* 1 : ( I J ! K_ L_)
                IXCHNG = 0
                ICOUL = 1
                CALL GETINT_DBG(SCR,ITYP,ISM,2,JTYP,JSM,2,KTYP,KSM,1,
     &           LTYP,LSM,1,IRI,IXCHNG,0,0,ICOUL)
*. Reorder to ( J I ! K L)
                 IREO(1) = 2
                 IREO(2) = 1
                 IREO(3) = 3
                 IREO(4) = 4
                 IILEN(1) = NI
                 IILEN(2) = NJ
                 IILEN(3) = NK
                 IILEN(4) = NL
C      REO4MT(XIN,XOUT,LIN,IREO,FACIN,FACOUT)
                 ONE = 1.0D0
                 ZERO = 0.0D0
                 CALL REO4MT(SCR,XINT,IILEN,IREO,ONE,ZERO)
*. 2 : (K_ J !L I_)
C                CALL GETINT_DBG(SCR,KTYP,KSM,1,JTYP,JSM,2,LTYP,LSM,1,
C    &                           ITYP,ISM,2,IRI,IXCHNG,0,0,ICOUL)
*. Reorder to (J I ! K L)
C                IREO(1) = 3
C                IREO(2) = 1
C                IREO(3) = 4
C                IREO(4) = 2
C                IILEN(1) = NK
C                IILEN(2) = NJ
C                IILEN(3) = NL
C                IILEN(4) = NI
C                CALL REO4MT(SCR,XINT,IILEN,IREO,ONE,ONE)
*
                 if (NTEST.ge.1000) then
                   WRITE(6,'(A,8I3)') ' Sym and type for ijkl ',
     &             ISM,ITYP,JSM,JTYP,KSM,KTYP,LSM,LTYP
                   WRITE(6,*) ' XINT as (J I ! K L ) '
                   CALL WRTMAT(XINT,NI*NJ,NK*NL,NI*NJ,NK*NL)
                 end if
*
* S(Ka,j,Ib) = sum(k,l,Jb)<Ib!a+kba lb!Jb>C(Ka,j,Jb)*(ji!kl)
*
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' C(Ka,j,JB) as C(Kaj,JB) '
                  CALL WRTMAT(CJRES,LKABTC*NJ,NJB,LKABTC*NJ,NJB)
                END IF
*
                CALL SKICKJ_REL(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,XINT,NI,NJ,NK,NL,NKBSTR,
     &               I4,XI4S,I2,XI2S,IKORD,
     &               FACS,IROUTE,MXTSOB)
                FACS = 1.0D0
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Updated S(Ka,i,IB) as S(Kai,IB) '
                  CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
*
*
 1930         CONTINUE
 2000       CONTINUE
*
*. Scatter out from s(Ka,Ib,i)
*
                IF(NTEST.GE.1000) THEN
                  WRITE(6,*) ' Final S(Ka,i,IB) as S(Kai,IB) '
                  CALL WRTMAT(SIRES,LKABTC*NI,NIB,LKABTC*NI,NIB)
                END IF
            DO II = 1, NI
              CALL ADD_SKAIIB(SB,NI,NIA,SIRES,LKABTC,NIB,II,
     &             I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
*
      CALL QEXIT('RS2B ')
      RETURN
      END
