      SUBROUTINE GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
     &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
     &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &                  SB,CB,ADSXA,STSTSX,MXPNGAS,
     &                  NOBPTS,IOBPTS,MAXK,
     &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
     &                  NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,IUSEAB,
     &                  CJRES,SIRES,NORB,NTESTG,SCLFAC)
*
* alpha-beta contribution to two-particle density matrix
* from given c-block and s-block.
*
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
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I2, XI2S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* X : Space for block of two-electron integrals
*
* Jeppe Olsen, Fall of 1996
*
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. General input
      INTEGER ADSXA(MXPOBS,MXPOBS),STSTSX(NSMST,NSMST)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*)
*.Input
      DIMENSION CB(*),SB(*)
*. Output
      DIMENSION RHO2(*)
*.Scratch
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION X(*)
      DIMENSION CJRES(*),SIRES(*)
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS),KTP(MXPNGAS),LTP(MXPNGAS)
*
      CALL QENTER('GSD2B')
      NTESTL = 000
      NTEST = MAX(NTESTL,NTESTG)
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================== '
        WRITE(6,*) ' GSBBD2B speaking '
        WRITE(6,*) ' ================== '
      END IF
C?    WRITE(6,*) ' NJAS NJB = ',NJA,NJB
C?    WRITE(6,*) ' IAGRP IBGRP = ', IAGRP,IBGRP
C?    WRITE(6,*) ' MXPNGAS = ', MXPNGAS
C?    WRITE(6,*) ' NSMOB = ', NSMOB
      IROUTE = 3
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
      CALL SXTYP_GAS(NKLTYP,KTP,LTP,NGAS,MXPNGAS,IBOC,JBOC)
      CALL SXTYP_GAS(NIJTYP,ITP,JTP,NGAS,MXPNGAS,IAOC,JAOC)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 9999
      DO 2001 IJTYP = 1, NIJTYP
        ITYP = ITP(IJTYP)
        JTYP = JTP(IJTYP)
        DO 1940 ISM = 1, NSMOB
          JSM = ADSXA(ISM,IJSM)
          IF(JSM.EQ.0) GOTO 1940
          KAFRST = 1
          if(ntest.ge.1500) write(6,*) ' ISM JSM ', ISM,JSM
          IOFF = IOBPTS(ITYP,ISM)
          JOFF = IOBPTS(JTYP,JSM)
          NI = NOBPTS(ITYP,ISM)
          NJ = NOBPTS(JTYP,JSM)
          IF(NI.EQ.0.OR.NJ.EQ.0) GOTO 1940
*. Generate annihilation mappings for all Ka strings
*. a+j!ka> = +/-/0 * !Ja>
          CALL ADSTN_GAS(JSM,JTYP,JATP,JASM,IAGRP,
     &                   I1,XI1S,NKASTR,IEND,IFRST,KFRST,KACT,
     &                   SCLFAC)
*. a+i!ka> = +/-/0 * !Ia>
          ONE    = 1.0D0
          CALL ADSTN_GAS(ISM,ITYP,IATP,IASM,IAGRP,
     &                   I3,XI3S,NKASTR,IEND,IFRST,KFRST,KACT,
     &                   ONE   )
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
          DO 1801 IKABTC = 1, NKABTC
C?          write(6,*) ' Batch over kstrings ', IKABTC
            KABOT = (IKABTC-1)*MAXK + 1
            KATOP = MIN(KABOT+MAXK-1,NKAEFF)
            LKABTC = KATOP-KABOT+1
*. Obtain C(ka,J,JB) for Ka in batch
            DO JJ = 1, NJ
              CALL GET_CKAJJB(CB,NJ,NJA,CJRES,LKABTC,NJB,
     &             JJ,I1(KABOT+(JJ-1)*NKASTR),
     &             XI1S(KABOT+(JJ-1)*NKASTR))
            END DO
*. Obtain S(ka,i,Ib) for Ka in batch
            DO II = 1, NI
              CALL GET_CKAJJB(SB,NI,NIA,SIRES,LKABTC,NIB,
     &             II,I3(KABOT+(II-1)*NKASTR),
     &             XI3S(KABOT+(II-1)*NKASTR))
            END DO
*
            DO 2000 KLTYP = 1, NKLTYP
              KTYP = KTP(KLTYP)
              LTYP = LTP(KLTYP)
*
              DO 1930 KSM = 1, NSMOB
                LSM = ADSXA(KSM,KLSM)
                IF(LSM.EQ.0) GOTO 1930
C?              WRITE(6,*) ' Loop 1930, KSM LSM ',KSM,LSM
                KOFF = IOBPTS(KTYP,KSM)
                LOFF = IOBPTS(LTYP,LSM)
                NK = NOBPTS(KTYP,KSM)
                NL = NOBPTS(LTYP,LSM)
*. If IUSEAB is used, only terms with i.ge.k will be generated so
                IKORD = 0
                IF(IUSEAB.EQ.1.AND.ISM.GT.KSM) GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.LT.KTYP)
     &          GOTO 1930
                IF(IUSEAB.EQ.1.AND.ISM.EQ.KSM.AND.ITYP.EQ.KTYP) IKORD=1
*
                IF(NK.EQ.0.OR.NL.EQ.0) GOTO 1930
*. Obtain all connections a+l!Kb> = +/-/0!Jb>
                ONE = 1.0D0
                CALL ADSTN_GAS(LSM,LTYP,JBTP,JBSM,IBGRP,
     &               I2,XI2S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE   )
                IF(NKBSTR.EQ.0) GOTO 1930
*. Obtain all connections a+k!Kb> = +/-/0!Ib>
                CALL ADSTN_GAS(KSM,KTYP,IBTP,IBSM,IBGRP,
     &               I4,XI4S,NKBSTR,IEND,IFRST,KFRST,KACT,ONE)
                IF(NKBSTR.EQ.0) GOTO 1930
*
*. Update two-electron density matrix
*  Rho2b(ij,kl) =  Sum(ka)S(Ka,i,Ib)<Ib!Eb(kl)!Jb>C(Ka,j,Jb)
*
                ZERO = 0.0D0
                CALL SETVEC(X,ZERO,NI*NJ*NK*NL)
*
C               WRITE(6,*) ' Before call to ABTOR2'
                CALL ABTOR2(SIRES,CJRES,LKABTC,NIB,NJB,
     &               NKBSTR,X,NI,NJ,NK,NL,NKBSTR,
     &               I4,XI4S,I2,XI2S,IKORD)
*. contributions to Rho2(ij,kl) has been obtained, scatter out
C?              WRITE(6,*) ' Before call to ADTOR2'
C?              WRITE(6,*) ' RHO2B (X) matrix '
C?              call wrtmat(x,ni*nj,nk*nl,ni*nj,nk*nl)
                CALL ADTOR2(RHO2,X,2,
     &                NI,IOFF,NJ,JOFF,NK,KOFF,NL,LOFF,NORB)
C?              write(6,*) ' updated density matrix '
C?              call prsym(rho2,NORB*NORB)

 1930         CONTINUE
 2000       CONTINUE
 1801     CONTINUE
*. End of loop over partitioning of alpha strings
 1940   CONTINUE
 2001 CONTINUE
*
 9999 CONTINUE
*
*
      CALL QEXIT('GSD2B')
      RETURN
      END
