      SUBROUTINE RSSBCB2(IASM,IATP,IOCPTA,IBSM,IBTP,IOCTPB,
     &                  JASM,JATP,JBSM,JBTP,NGAS,
     &                  IAOC,IBOC,JAOC,JBOC,NAEL,NBEL,
     &                  IJAGRP,IJBGRP,SB,CB,JDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  NIA,NIB,NJA,NJB,MXPOBS,IDC,PS,
     &                  ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                  MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,IPRNT,IHAPR,
     &                  IPTSPC,JPTSPC,IJOP,NNSEL2E,ISEL2E,SCLFAC,
     &                  IUSE_PH,IPHGAS,I_RES_AB,
     &                  MXTSOB,SOMESCR,XIJILS)
**************************************************************
* WARNING!!
*  The last line of parameters has not yet been accounted for !!
*  (&     MXTSOB,SOMESCR,XIJILS)
*  Pass to relevant subroutines, allocate and modify !!
* WARNING!!
**************************************************************
*
* Contributions to sigma block (iasm iatp, ibsm ibtp ) from
* C block (jasm jatp , jbsm, jbtp)
*
* =====
* Input
* =====
*
* IASM,IATP : Symmetry and type of alpha strings in sigma
* IBSM,IBTP : Symmetry and type of beta  strings in sigma
* JASM,JATP : Symmetry and type of alpha strings in C
* JBSM,JBTP : Symmetry and type of beta  strings in C
* NGAS      : Number of active spaces in calculation
* IAOC,IBOC : Number of electrons in each AS for sigma supergroups
* JAOC,JBOC : Number of electrons in each AS for C     supergroups
* NAEL : Number of alpha electrons
* NBEL : Number of  beta electrons
* IJAGRP    : IA and JA belongs to this group of strings
* IJBGRP    : IB and JB belongs to this group of strings
* CB : Input c block
* IDOH2 : = 0 => no two electron operator
* IDOH2 : = 1 =>    two electron operator
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
*
* IHAPR : .ne. 0 implies thatt the exact Hamiltonian shoulf not be uses
* In the case IPTSPC and JPTSPC defined the perturbation spaces
* a nonvanishing perturbation is allowed inside each subspace.
* The actual type of approximate Hamiltonian in each subspace is defined by
* IHFORM
* NNSEL2E : Only selected 2e terms will be included
* ISEL2E : orbital spaces in which 2e terms are included
*          (Currently : all indeces identical )
*
* ======
* Output
* ======
* SB : fresh sigma block
*
* =======
* Scratch
* =======
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* I1, XI1S   : at least MXSTSO : Largest number of strings of given
*              type and symmetry
* C2 : Must hold largest STT block of sigma or C
*
* XINT : Scratch space for integrals.
*
* Jeppe Olsen , Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
*. Output
      DIMENSION CB(*),SB(*)
*. Scratch
      DIMENSION SSCR(*),CSCR(*),I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*)
      DIMENSION C2(*)
      DIMENSION CJRES(*),SIRES(*)
      DIMENSION IBLOCK(8)
      DIMENSION IPHGAS(*)
*. For H(apr)
      DIMENSION ISEL2E(*)
*
      NTESTO= NTEST
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
*
C?    IF(.NOT. (IATP.EQ.3.AND.IBTP.EQ.1) )THEN
C?      WRITE(6,*) ' Enforced return from RSSBCB2 '
C?      RETURN
C?    END IF
C?    IF(.NOT.(JATP.EQ.4.AND.JBTP.EQ.1)) THEN
C?      WRITE(6,*) ' Enforced return from RSSBCB2 '
C?      RETURN
C?    END IF
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' ==============================='
        WRITE(6,*) ' RSSBCB2 :  C block (transposed)'
        WRITE(6,*) ' ================================'
        CALL WRTMAT(CB,NJB,NJA,NJB,NJA)
        WRITE(6,*) ' ======================================='
        WRITE(6,*) ' RSSBCB2 : Initial  S block(transposed) '
        WRITE(6,*) ' ======================================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
        WRITE(6,*) ' Overall scalefactor ',SCLFAC
        WRITE(6,*) ' IHAPR,JDOH2 = ', IHAPR,JDOH2
        WRITE(6,*) ' IUSE_PH,I_RES_AB = ', IUSE_PH,I_RES_AB
      END IF
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' IAOC and IBOC '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' JAOC and JBOC  : '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IASM IATP JASM JATP ', IASM,IATP,JASM,JATP
        WRITE(6,*) ' IBSM IBTP JBSM JBTP ', IBSM,IBTP,JBSM,JBTP
        WRITE(6,*) ' NAEL NBEL ', NAEL, NBEL
      END IF
* Should the corresponding Hamiltonian matrix block be
* calculated exactly or approximately
      IF(IHAPR.NE.0) THEN
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
        IF(NTEST.GE. 20)
     &  WRITE(6,*) ' RSSBCBN : ', NNSEL2E, ISEL2E(1)
        NSEL2E = NNSEL2E
        IF(IMZERO.NE.0) GOTO 9999
      ELSE
*. Operator specified by input
        IAPRLEV =-1
        IDOH2 = JDOH2
        IDIAG = 0
        NSEL2E = 0
      END IF
      IF(NTEST.GE. 20)
     &WRITE(6,*) ' IHAPR, IDIAG IDOH2 ' , IHAPR,IDIAG, IDOH2
*
*
      IF(IDC.EQ.2.AND.IATP.EQ.IBTP.AND.IASM.EQ.IBSM .AND.
     &   I_RES_AB.EQ.0.AND.JASM.EQ.JBSM.AND.JATP.EQ.JBTP) THEN
*. Diagonal sigma block, use alpha-beta symmetry to reduce computations.
        IUSEAB = 1
      ELSE
        IUSEAB = 0
      END IF
*
      IF(IDIAG.EQ.0) THEN
*
* Calculate block exactly
*
      IF(I_RES_AB.NE.1.AND.IUSEAB.EQ.0.
     &   AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* =============================
* Sigma beta beta contribution
* =============================
*
* Sigma aa(IA,IB) = sum(i.gt.k,j.gt.l)<IB!Eb(ij)Eb(kl)!JB>
*                 * ((ij!kl)-(il!kj)) C(IA,JB)
*                 + sum(ij) <IB!Eb(ij)!JB> H(ij) C(IA,JB)
*.One electron part
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)

        IF(NBEL.GE.0) THEN
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB before RSBB1E'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB1E'
C         CALL RSBB1E(IBSM,IBTP,IOCTPB,JBSM,JBTP,IOCTPB,
C    &         IJBGRP,NIA,NGAS,IBOC,JBOC,
C    &         SB,CB,ADSXA,SXSTST,STSTSX,
C    &         NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
C    &         SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
C    &         NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
C    &         NIB,MXSXST,CJRES(1),CJRES(1+MXSXBL),MOCAA,SCLFAC,
C    &         IUSE_PH,IPHGAS,NTEST)
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB after RSBB1E'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
        END IF
        IF(IDOH2.NE.0.AND.NBEL.GE.0) THEN
*. Two electron part
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB2A'
          CALL RSBB2A(IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,NIB,
     &                NGAS,IBOC,JBOC,SB,CB,
     &                ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &                NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &                CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC,NTEST,0,0,
     &                IUSE_PH,IPHGAS)
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB after RSBB2a'
             call wrtmat(sb,nib,nia,nib,nia)
           END IF
        END IF
        CALL TRPMT3(SB,NIA,NIB,C2)
        CALL COPVEC(C2,SB,NIA*NIB)
        CALL TRPMT3(CB,NJA,NJB,C2)
        CALL COPVEC(C2,CB,NJA*NJB)
      END IF
*
* =============================
* Sigma alpha beta contribution
* =============================
*
      IF(IDOH2.NE.0.AND.NAEL.GE.0.AND.NBEL.GE.0) THEN
        IF(NTEST.GE.101)
     &  WRITE(6,*) ' I am going to call RSBB2B'
        IIITRNS = 1
        IF(IIITRNS.EQ.1) THEN
*. Call advice routine
C     ADVICE_SIGMA(IAOCC,IBOCC,JAOCC,JBOCC,ITERM,LADVICE)
           CALL ADVICE_SIGMA(IAOC,IBOC,JAOC,JBOC,1,LADVICE)
*. LADVICE = 2 => implies transpose
           IF(LADVICE.EQ.2) THEN
             JJJTRNS = 1
           ELSE
             JJJTRNS = 0
           END IF
        ELSE
           JJJTRNS = 0
        END IF
*
        IF (JJJTRNS.EQ.0) THEN
          CALL RSBB2BN(IASM,IATP,IBSM,IBTP,NIA,NIB,
     &         JASM,JATP,JBSM,JBTP,NJA,NJB,
     &         IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
     &         NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS)
*
         ELSE IF ( JJJTRNS.EQ.1) THEN
*. well lets give the transpose routine some more practice : Transpose back
          CALL TRPMT3(SB,NIB,NIA,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
*
          CALL TRPMT3(CB,NJB,NJA,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
C         WRITE(6,*) ' RSSBCB2 : Transpose path choosen'
*
          CALL RSBB2BN(IBSM,IBTP,IASM,IATP,NIB,NIA,
     &         JBSM,JBTP,JASM,JATP,NJB,NJA,
     &         IJBGRP,IJAGRP,NGAS,IBOC,IAOC,JBOC,JAOC,
     &         SB,CB,ADSXA,STSTSX,MXPNGAS,
     &         NOBPTS,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &         XINT,NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         IUSEAB,CJRES,SIRES,SCLFAC,NTEST,0,0,IUSE_PH,IPHGAS)

*. Transpose ( To compensate later transposition )
          CALL TRPMT3(SB,NIA,NIB,C2)
          CALL COPVEC(C2,SB,NIA*NIB)
          CALL TRPMT3(CB,NJA,NJB,C2)
          CALL COPVEC(C2,CB,NJA*NJB)
        END IF
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB after RSBB2B, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB after RSBB2b'
             call wrtmat(sb,nia,nib,nia,nib)
           END IF
      END IF
*
* =============================
* Sigma alpha alpha contribution
* =============================
*
      IF(I_RES_AB.NE.-1.AND.
     &   NAEL.GE.0.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* alpha single excitation
*
        IF(NTEST.GE.101)
     &  WRITE(6,*) ' I am going to call RSBB1E (last time )'
C       CALL RSBB1E(IASM,IATP,IOCTPA,JASM,JATP,IOCTPA,
C    &                   IJAGRP,NIB,
C    &                   NGAS,IAOC,JAOC,
C    &                   SB,CB,
C    &                   ADSXA,SXSTST,STSTSX,
C    &                   NOBPTS,IOBPTS,
C    &                   ITSOB,MAXI,MAXK,
C    &                   SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
C    &                   NSMOB,NSMST,NSMSX,MXPOBS,MOCAA,
C    &                   NIA,MXSXST,CJRES(1),CJRES(1+MXSXBL),
C    &                   MOCAA,SCLFAC,IUSE_PH,IPHGAS,NTEST)
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB transposed after RSBB1, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
           IF(NTEST.GE.500) THEN
             WRITE(6,*) ' SB transposed  after RSBB1E'
             call wrtmat(SB,nib,nia,nib,nia)
           END IF
*
* alpha double excitation
*
        IF(IDOH2.NE.0.AND.NAEL.GE.0) THEN
          IF(NTEST.GE.101)
     &    WRITE(6,*) ' I am going to call RSBB2A (last time )'
          CALL RSBB2A(IASM,IATP,JASM,JATP,IJAGRP,NIB,NIA,
     &         NGAS,IAOC,JAOC,
     &         SB,CB,
     &         ADSXA,DXSTST,STSTDX,SXDXSX,MXPNGAS,
     &         NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &         SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &         NSMOB,NSMST,NSMSX,NSMDX,MXPOBS,
     &         CJRES,SIRES,MXSXST,MXSXBL,MOCAA,SCLFAC,NTEST,0,0,
     &         IUSE_PH,IPHGAS)
        END IF
*
           IF(NTEST.GE.101) THEN
             WRITE(6,*) ' SB transposed after RSBB2A, first element '
             call wrtmat(sb,1,1    ,nia,nib)
           END IF
        IF(NTEST.GE.500) THEN
          WRITE(6,*) ' SB after RSBB2A'
          call wrtmat(sb,nia,nib,nia,nib)
        END IF
      END IF
*
      ELSE IF (IDIAG.EQ.1) THEN
*
*. Diagonal approxiation
*
       IBLOCK(1) = IATP
       IBLOCK(2) = IBTP
       IBLOCK(3) = IASM
       IBLOCK(4) = IBSM
       IBLOCK(5) = 1
       IBLOCK(6) = 1
       IF(IDOH2.EQ.0) THEN
         I12 = 1
       ELSE
         I12 = 2
       END IF
C?     WRITE(6,*) ' IDOH2, I12 ', IDOH2,I12
       ITASK = 2
       FACTOR = 0.0D0
*. Input is in det basis
       IIDC = 1
*. Well, we are working with transposed matrices so
       CALL TRPMT3(CB,NIB,NIA,C2)
C?     WRITE(6,*) ' DIATERM2_GAS will be called '
       CALL DIATERM2_GAS(FACTOR,ITASK,C2,1,IBLOCK,1,0,I12,IIDC)
       CALL TRPMT3(C2,NIA,NIB,CB)
       ONE = 1.0D0
       IF(IUSEAB.EQ.0) THEN
         FACTOR = 1.0D0*SCLFAC
       ELSE
         FACTOR = 0.5D0*SCLFAC
       END IF
       CALL VECSUM(SB,SB,CB,ONE,FACTOR,NIA*NIB)
      END IF
*
 9999 CONTINUE
      IF(IHAPR.NE.0) THEN
*. Clean up
        CALL HMATAPR(IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &               IPTSPC,JPTSPC,IJOP,IJOP,IIF,JDOH2,IDOH2,
     &               IMZERO,IDIAG)
      END IF
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' ==================================='
        WRITE(6,*) ' RSSBCB : Final S block (transposed)'
        WRITE(6,*) ' ==================================='
        CALL WRTMAT(SB,NIB,NIA,NIB,NIA)
      END IF
      NTESTO = NTEST
      RETURN
      END
