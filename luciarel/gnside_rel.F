      SUBROUTINE GNSIDE_REL(ISD,IRIS,IRIC,ICA,ICB,IAA,IAB,
     &           IAOC,IBOC,JAOC,JBOC,
     &           NIA,NIB,NJA,NJB,
     &           TB,IXTP,
     &           SB,CB,ISSM,ICSM,
     &           I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &           TBSUB,CJRES,SIRES,MAXLB,
     &           ICA_STR,ICB_STR,IAA_STR,IAB_STR,
     &           KJA,XKJA,KJB,XKJB,KIA,XKIA,KIB,XKIB,
     &           KJAD,XKJAD,KJBD,XKJBD,KIAD,XKIAD,KIBD,XKIBD,
     &           ITDIAG,ITABTRNSP,IC_RESTRICT,IS_RESTRICT,
     &           SIGN_EXTERNAL,IPRS,IPRD)
*
* Sigma routine and density routine for for general operator 
*
* sum(Ica,Icb,Iaa,Iab) T(Ica,Icb,Iaa,Iab) Ica+ Icb+ Iaa Iab
*
* contribution from C vector with occupation JAOC,IBOC
* to Sigma vector with occupation IAOC, IBOC
*
* All symmetryblocks of C and Sigma are treated 
*
* ITDIAG = 1 => Diagonal T-block 
*               Only symmetryblocks with 
*                  (ICASM,IAASM).GE.(ICBSM,IABSM)
*               are included, and blocks with 
*                  (ICASM,IAASM) = (ICBSM,IABSM)
*               are given in diagonal form 
*               (Iaa.ge.Iab, if Iaa=Iab, Ica .ge. Icb)
* ITABTRSNP = 1 => alpha-beta transposed block, input block is 
*               T(ICB,ICA,IAB,IAA)
* IX_RESTRICT = 1: BLocks in X are restricted to those with IASM.GE.IBSM
* IX_RESTRICT =-1: BLocks in X are restricted to those with IASM.LT.IBSM
*
* SIGN_EXTERNAL : A sign multiplying all contributions
*                 ( And he said : Mighty Lord, Give me a sign)
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
*
* Jeppe Olsen, August 1999    
* Timo Fleig, June 2001, modified for relativistic case
*             cont´d, 11th of July, Mallorca, around noon (I´m crazy)
*                     15th of July
*                     August 2001
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "orbinp.inc"
#include "csm.inc"
#include "symm.inc"
#include "cgas.inc"
#include "ctcc.inc"
*     
*.Input
      DIMENSION CB(*)
      INTEGER NIA(*),NIB(*),NJA(*),NJB(*)
      INTEGER IAOC(*),IBOC(*),JAOC(*),JBOC(*)
      INTEGER ICA(*),ICB(*),IAA(*),IAB(*)
*.Input or Output
      DIMENSION SB(*)
      dimension TB(*)
*.Local Scratch
      INTEGER KAOC(MXPNGAS),KBOC(MXPNGAS)
      INTEGER NIKAINTM(MXPOBS,MXPLCCOP),NIKBINTM(MXPOBS,MXPLCCOP)
      INTEGER NJKAINTM(MXPOBS,MXPLCCOP),NJKBINTM(MXPOBS,MXPLCCOP)
      INTEGER NKA(MXPOBS),NKB(MXPOBS)
*. 
      INTEGER NICA(MXPOBS),NICB(MXPOBS),NIAA(MXPOBS),NIAB(MXPOBS)
      INTEGER ICA_EXP(MXNOP),ICB_EXP(MXNOP),
     &        IAA_EXP(MXNOP),IAB_EXP(MXNOP)
      INTEGER ICAGP(MXPNGAS),ICBGP(MXPNGAS)
      INTEGER IAAGP(MXPNGAS),IABGP(MXPNGAS)
      integer KAGP(MXPNGAS),KBGP(MXPNGAS)
      integer JAGP(MXPNGAS),JBGP(MXPNGAS)
      INTEGER IAC_AR(MXNOP)
      INTEGER IB_CA(MXNDGIRR,MXNDGIRR,MXPLCCOP),
     &        IB_CB(MXNDGIRR,MXNDGIRR,MXPLCCOP)
      INTEGER IB_AA(MXNDGIRR,MXNDGIRR,MXPLCCOP),
     &        IB_AB(MXNDGIRR,MXNDGIRR,MXPLCCOP)
      INTEGER IB_T(MXNDGIRR,MXNDGIRR,MXNDGIRR)
      INTEGER IB_C(MXNDGIRR),IB_S(MXNDGIRR)
*. Scratch through input 
      DIMENSION I1(*),XI1S(*),I2(*),XI2S(*)
      DIMENSION I3(*),XI3S(*),I4(*),XI4S(*)
      DIMENSION TBSUB(*)
      INTEGER KJA(MAXLB),KJB(MAXLB),KIA(MAXLB),KIB(MAXLB)
      DIMENSION XKJA(MAXLB),XKJB(MAXLB),XKIA(MAXLB),XKIB(MAXLB)
      INTEGER KJAD(MAXLB),KJBD(MAXLB),KIAD(MAXLB),KIBD(MAXLB)
      DIMENSION XKJAD(MAXLB),XKJBD(MAXLB),XKIAD(MAXLB),XKIBD(MAXLB)
*. ... Must hold excitations for all intermediate strings of given sym 
*. and all orbitals of given type
      INTEGER ICA_STR(*),ICB_STR(*),IAA_STR(*),IAB_STR(*)
*     ... Must hold strings of given sym for T ops.
      DIMENSION CJRES(*),SIRES(*)
*
      COMMON/KKKDUMMY/LEN_C2,LEN_S2
      COMMON/CMXCJ/MXCJ
*
      CALL QENTER('GENSIG')
*
      NTESTL = 0000 
      NTEST = MAX(NTESTL,max(IPRS,IPRD))
*
      ISKIP = 0
      LEN_C = LEN_TT_BLOCK(ICSM,NJA,NJB,NSMST)
      XCNORM = INPROD(CB,CB,LEN_C)
      IF(XCNORM.EQ.0.0D0) ISKIP = 1
      if (ISKIP.EQ.1) then
        if (NTEST.ge.1000) 
     &      write(6,*) 'SKIPPING gnside_rel run. LEN_C = ',LEN_C
        goto 9999
      end if
*
      if (NTEST.ge.100) write(6,*) 'Coupling with type ',IXTP
*
      IF(NTEST.GE.500) THEN
        WRITE(6,*) ' ================ '
        WRITE(6,*) ' GNSIDE speaking '
        WRITE(6,*) ' ================ '
        write(6,*) ' Real / imaginary contribution :'
        write(6,*) ' Left , Right ',IRIS,IRIC
        WRITE(6,*)
        WRITE(6,*) ' IAOC and IBOC in GENSIG'
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' JAOC and JBOC in GENSIG'
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Type of operator '
        CALL WRT_TP_GENOP(ICA,ICB,IAA,IAB,NGAS)
        IF(ISD.EQ.1) THEN
          WRITE(6,*) ' Sigma generation '
        ELSE 
          WRITE(6,*) ' Density generation '
      END IF
*
      WRITE(6,*) ' Input C vector '
      CALL WRTVH1(CB,ICSM,NJA,NJB,NSMST,0)
      if (ISD.eq.2) then
         WRITE(6,*) ' Input S vector '
         CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
       end if
      END IF
*
*. A few constants
      IONE = 1
      ZERO = 0.0D0
      ONE = 1.0D0
*. Symmetry of T
      IOPSM = IDBGMULT(INVELM(ICSM),ISSM)
*. Offsets for symmetryblocks of C and S
*. Total length of C and sigma blocks
      LEN_S = LEN_TT_BLOCK(ISSM,NIA,NIB,NSMST)
*
      ISKIP = 0
      IF(ISD.EQ.2) THEN
       XSNORM = INPROD(SB,SB,LEN_S)
       IF(XSNORM.EQ.0.0D0) ISKIP = 1
      END IF
      IF(ISKIP.EQ.1) GOTO 9999
*. Offset of C and sigma routines
      CALL IOFF_SYMBLK_MAT(NSMST,NJA,NJB,ICSM,IB_C,IC_RESTRICT)
      CALL IOFF_SYMBLK_MAT(NSMST,NIA,NIB,ISSM,IB_S,IS_RESTRICT)
*. types of strings in T in groupnotation
      CALL OCC_TO_GRP(ICA,ICAGP,1)
      CALL OCC_TO_GRP(ICB,ICBGP,1)
      CALL OCC_TO_GRP(IAA,IAAGP,1)
      CALL OCC_TO_GRP(IAB,IABGP,1)
C          OCC_TO_GRP(IOCC,IGRP,IWAY)
*
* Why is the norm of the C vector commented out here? JEPPE
C     XCNORM = INPROD(CB,CB,LEN_C)
C     IF(XCNORM.EQ.0.0D0) GOTO 9999
*. Operators in T in expanded form
      CALL REF_OP(ICA,ICA_EXP,NCA_OP,NGAS,1,MXNOP)
      CALL REF_OP(ICB,ICB_EXP,NCB_OP,NGAS,1,MXNOP)
      CALL REF_OP(IAA,IAA_EXP,NAA_OP,NGAS,1,MXNOP)
      CALL REF_OP(IAB,IAB_EXP,NAB_OP,NGAS,1,MXNOP)
C          REF_OP(IOPGAS,IOP,NOP,NGAS,IWAY)
*. Sign for bringing operators and amplitudes into new order
      NPERM =  NAA_OP*NCB_OP 
     &       + NAA_OP*(NAA_OP-1)/2
     &       + NAB_OP*(NAB_OP-1)/2
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNXXX  =  -1.0D0*SIGN_EXTERNAL
      ELSE
        SIGNXXX = SIGN_EXTERNAL
      END IF
*. Type of Ka and Kb
      CALL CCEX_OCC_OCC(JAOC,KAOC,NGAS,1,IAA,IKA_ZERO)
      CALL CCEX_OCC_OCC(JBOC,KBOC,NGAS,1,IAB,IKB_ZERO)
* Check that annihilation has not acted upon empty GAS:
      IF(IKA_ZERO.EQ.0.AND.IKB_ZERO.EQ.0) THEN
* Get group of K mapped occupations
       call occ_to_grp(KAOC,KAGP,1)
       call occ_to_grp(KBOC,KBGP,1)
* Number of strings in this supergroup (new DBG routine)
       do KSM = 1,NSMST,1
         call nst_spgrp_dbg(NGAS,KAGP,KSM,NSMST,NKA(KSM),NKADST,1)
         call nst_spgrp_dbg(NGAS,KBGP,KSM,NSMST,NKB(KSM),NKBDST,2)
       end do
*. Number of strings in excitation ops
       do ISYM = 1,NSMST,1
         call nst_spgrp_dbg(NGAS,ICAGP,ISYM,NSMST,NICA(ISYM),NCADST,1)
         call nst_spgrp_dbg(NGAS,ICBGP,ISYM,NSMST,NICB(ISYM),NCBDST,2)
*. All arrays refer to strings of creation operators of given occ.
*. So do not bother with difference between anni and crea here.
         call nst_spgrp_dbg(NGAS,IAAGP,ISYM,NSMST,NIAA(ISYM),NAADST,1)
         call nst_spgrp_dbg(NGAS,IABGP,ISYM,NSMST,NIAB(ISYM),NABDST,2)
       end do
*. Offset for symmetryblocks of T
       IF(ITABTRNSP.EQ.0) THEN
         CALL Z_TCC_OFF(IB_T,NICA,NICB,NIAA,NIAB,IOPSM,NSMST,ITDIAG)
       ELSE
         CALL Z_TCC_OFF(IB_T,NICB,NICA,NIAB,NIAA,IOPSM,NSMST,ITDIAG)
       END IF
C             Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
*. Creation path from C/S strings to K strings
*. Mappings from !I> to !K> are annihilations mappings 
       NMXOP_T = MAX(NCA_OP,NCB_OP,NAA_OP,NAB_OP)
       DO JOP = 1, NMXOP_T
         IAC_AR(JOP) = 1
       END DO
       ONE = 1.0D0
* CA : KA => IA mappings in 1
C           MAP_EXSTR(IOP,IAC,NOP,IREFOC,IX,SX,NK,IB,SCLFAC,IAB)
       IAB_MAP = 1
       CALL MAP_EXSTR(ICA_EXP,IAC_AR,NCA_OP,IAOC,I1,XI1S,
     &                NIKAINTM,IB_CA,ONE,IAB_MAP)
* CB : KB => IB     
       IAB_MAP = 2
       CALL MAP_EXSTR(ICB_EXP,IAC_AR,NCB_OP,IBOC,I2,XI2S,
     &                NIKBINTM,IB_CB,ONE,IAB_MAP)
* AA : KA => JA
       IAB_MAP = 1
       CALL MAP_EXSTR(IAA_EXP,IAC_AR,NAA_OP,JAOC,I3,XI3S,
     &                NJKAINTM,IB_AA,ONE,IAB_MAP)
* AB : KB => JB
       IAB_MAP = 2
       CALL MAP_EXSTR(IAB_EXP,IAC_AR,NAB_OP,JBOC,I4,XI4S,
     &                NJKBINTM,IB_AB,ONE,IAB_MAP)
*. Loop over symmetry-blocks of C and Sigma
       DO IASM = 1, NSMST 
        IBSM = IDBGMULT(ISSM,INVELM(IASM))
*
        IS_INCLUDED = 1
        IF(IS_RESTRICT.EQ.1.AND.IASM.LT.IBSM) IS_INCLUDED = 0
        IF(IS_RESTRICT.EQ.-1.AND.IASM.GE.IBSM) IS_INCLUDED = 0
*
        IF(IS_INCLUDED.EQ.1) THEN
        DO JASM = 1, NSMST
         JBSM = IDBGMULT(ICSM,INVELM(JASM))
*
         IC_INCLUDED = 1
         IF(IC_RESTRICT.EQ.1.AND.JASM.LT.JBSM) IC_INCLUDED = 0
         IF(IC_RESTRICT.EQ.-1.AND.JASM.GE.JBSM) IC_INCLUDED = 0
*
         IF(IC_INCLUDED.EQ.1) THEN
         IF(NTEST.GE.1000) WRITE(6,'(A,4I4)') 
     &   ' IASM, IBSM, JASM, JBSM ', IASM,IBSM,JASM,JBSM
         DO KASM = 1, NSMST
*. IAA_SM is symmetry of ANNIHILATOR string
          IAA_SM = IDBGMULT(KASM,INVELM(JASM))
          ICA_SM = IDBGMULT(IASM,INVELM(KASM))
* To obtain the number of annnihilator strings of a given sym
*. we find the number of creator strings of conjugated symmetry
          LIAA_SM = NIAA(INVELM(IAA_SM))
          LICA_SM = NICA(ICA_SM)
          LKA_SM = NKA(KASM)
*
          DO KBSM = 1, NSMST
*
           IAB_SM = IDBGMULT(KBSM,INVELM(JBSM))
           ICB_SM = IDBGMULT(IBSM,INVELM(KBSM))
*
           LIAB_SM = NIAB(INVELM(IAB_SM))
           LICB_SM = NICB(ICB_SM)
*
           LKB_SM = NKB(KBSM)
*
           ICBOFF = IB_C(JASM)
           ISBOFF = IB_S(IASM)
*
           IF(LIAA_SM*LICA_SM*LIAB_SM*LICB_SM*LKA_SM*LKB_SM.NE.0)
     &     THEN
*
             IF(NTEST.GE.700) THEN
               WRITE(6,*) ' ICA_SM, ICB_SM, IAA_SM,IAB_SM =',
     &                      ICA_SM, ICB_SM, IAA_SM,IAB_SM
             END IF 
             CALL FLUSH(6)
*
             IF(NTEST.GE.1000) 
     &       WRITE(6,'(A,2I4)') ' KASM,KBSM ', KASM, KBSM
             IF(NTEST.GE.1000) 
     &       WRITE(6,*) ' LKA_SM, LKB_SM', LKA_SM, LKB_SM
             IF(NTEST.GE.1000) 
     &       WRITE(6,'(A,4I3)')      
     &       'LIAA_SM, LICA_SM,LIAB_SM,LICB_SM',
     &        LIAA_SM, LICA_SM,LIAB_SM,LICB_SM
*
*. The operator strings of ICA, ICB, IAA, IAB 
            IUB = 1
            CALL GETSTR2_TOTSM_SPGP(IUB,ICAGP,NGAS,ICA_SM,NCA_OP,
     &           LICA_SM2,ICA_STR,NORBT,0,0,0)
            IUB = 2
            CALL GETSTR2_TOTSM_SPGP(IUB,ICBGP,NGAS,ICB_SM,NCB_OP,
     &           LICB_SM2,ICB_STR,NORBT,0,0,0)
            IUB = 1
* We are generating strings of creation operators so 
            IAA_SM_ADJ = INVELM(IAA_SM)
            CALL GETSTR2_TOTSM_SPGP(IUB,IAAGP,NGAS,IAA_SM_ADJ,NAA_OP,
     &           LIAA_SM2,IAA_STR,NORBT,0,0,0)
            IUB = 2
            IAB_SM_ADJ = INVELM(IAB_SM)
            CALL GETSTR2_TOTSM_SPGP(IUB,IABGP,NGAS,IAB_SM_ADJ,NAB_OP,
     &           LIAB_SM2,IAB_STR,NORBT,0,0,0)
* Batching                   
           LEN_KA_BAT = MAXLB
           LEN_KB_BAT = MAXLB
*
           IF(LKA_SM.LT.MAXLB) THEN
             LEN_KA_BAT = LKA_SM
             LEN_KB_BAT = MIN(LKB_SM,MAXLB*MAXLB/LKA_SM)
           END IF
*
           IF(LKB_SM.LT.MAXLB) THEN
             LEN_KB_BAT = LKB_SM
             LEN_KA_BAT = MIN(LKA_SM,MAXLB*MAXLB/LKB_SM)
           END IF
*
           NKA_BT = LKA_SM/LEN_KA_BAT
           IF(NKA_BT*LEN_KA_BAT. LT. LKA_SM) NKA_BT = NKA_BT + 1
           NKB_BT = LKB_SM/LEN_KB_BAT
           IF(NKB_BT*LEN_KB_BAT. LT. LKB_SM) NKB_BT = NKB_BT + 1
*
           NCA_BT = LICA_SM/MAXLB
           IF(NCA_BT*MAXLB.LT.LICA_SM) NCA_BT = NCA_BT + 1
           NCB_BT = LICB_SM/MAXLB
           IF(NCB_BT*MAXLB.LT.LICB_SM) NCB_BT = NCB_BT + 1
           NAA_BT = LIAA_SM/MAXLB
           IF(NAA_BT*MAXLB.LT.LIAA_SM) NAA_BT = NAA_BT + 1
           NAB_BT = LIAB_SM/MAXLB
           IF(NAB_BT*MAXLB.LT.LIAB_SM) NAB_BT = NAB_BT + 1
*. Loop over batches
           DO IKA_BT = 1, NKA_BT
           DO IKB_BT = 1, NKB_BT
            IKA_B = (IKA_BT-1)*LEN_KA_BAT + 1
            IKA_E = MIN(LKA_SM,IKA_B + LEN_KA_BAT - 1)
            NKA_B = IKA_E - IKA_B + 1
            NKA_B0 = NKA_B
            IKB_B = (IKB_BT-1)*LEN_KB_BAT + 1
            IKB_E = MIN(LKB_SM,IKB_B +  LEN_KB_BAT - 1)
            NKB_B = IKB_E - IKB_B + 1
            NKB_B0 = NKB_B
*. Loop over annihilation batches
            DO IAA_BT = 1, NAA_BT
             IAA_B = (IAA_BT-1)*MAXLB + 1
             IAA_E = MIN(LIAA_SM,IAA_B + MAXLB -1)
             NAA_B = IAA_E-IAA_B+1
*. KA => JA mapping
             IUB_KJ = 1
             CALL K_TO_J_TOT(KJAD,XKJAD,KASM,IKA_B,IKA_E,
     &                       IAA_STR,IAA_SM,NAA_B,IAA_B,IAA_E,
     &                       I3,XI3S,IB_AA,NJKAINTM,NAA_OP,IZEROKAJA,
     &                       IUB_KJ)
            IF(IZEROKAJA.EQ.0) THEN
            DO IAB_BT = 1, NAB_BT
*B(eginning) and E(nd) of each Batch
             IAB_B = (IAB_BT-1)*MAXLB + 1
             IAB_E = MIN(LIAB_SM,IAB_B + MAXLB -1)
             NAB_B = IAB_E-IAB_B+1
*. KB => JB mapping
             IUB_KJ = 2
             CALL K_TO_J_TOT(KJBD,XKJBD,KBSM,IKB_B,IKB_E,
     &                       IAB_STR,IAB_SM,NAB_B,IAB_B,IAB_E,
     &                       I4,XI4S,IB_AB,NJKBINTM,NAB_OP,IZEROKBJB,
     &                       IUB_KJ)
            IF(IZEROKBJB.EQ.0) THEN
*. Loop over creation batches 
            DO ICA_BT = 1, NCA_BT
             ICA_B = (ICA_BT-1)*MAXLB + 1
             ICA_E = MIN(LICA_SM,ICA_B + MAXLB -1)
             NCA_B = ICA_E-ICA_B+1
*. KA => IA mapping
             IUB_KJ = 1
             CALL K_TO_J_TOT(KIAD,XKIAD,KASM,IKA_B,IKA_E,
     &                       ICA_STR,ICA_SM,NCA_B,ICA_B,ICA_E,
     &                       I1,XI1S,IB_CA,NIKAINTM,NCA_OP,IZEROKAIA,
     &                       IUB_KJ)
            DO ICB_BT = 1, NCB_BT
             ICB_B = (ICB_BT-1)*MAXLB + 1
             ICB_E = MIN(LICB_SM,ICB_B + MAXLB -1)
             NCB_B = ICB_E-ICB_B+1
*. KB => IB mapping
             IUB_KJ = 2
             CALL K_TO_J_TOT(KIBD,XKIBD,KBSM,IKB_B,IKB_E,
     &                       ICB_STR,ICB_SM,NCB_B,ICB_B,ICB_E,
     &                       I2,XI2S,IB_CB,NIKBINTM,NCB_OP,IZEROKBIB,
     &                       IUB_KJ)
             IF(IZEROKBJB.EQ.0.AND.IZEROKAJA.EQ.0.AND.
     &          IZEROKBIB.EQ.0.AND.IZEROKAIA.EQ.0     ) THEN
*. Compress Ka and Kb strings to active subset 
C       COMPRS2LST(I1,XI1,N1,I2,XI2,N2,NKIN,NKOUT)
*. Ka
             CALL COMPRS2LST_B(KIAD,XKIAD,NCA_B,KJAD,XKJAD,NAA_B,
     &                         NKA_B0,NKA_B,KIA,XKIA,KJA,XKJA)
*. Kb
             CALL COMPRS2LST_B(KIBD,XKIBD,NCB_B,KJBD,XKJBD,NAB_B,
     &                         NKB_B0,NKB_B,KIB,XKIB,KJB,XKJB)
             IM_ACTIVE = 1
* Offset and form of T-coefficients
             IF(ITABTRNSP.EQ.0) THEN
               ITOFF = IB_T(ICA_SM,ICB_SM,INVELM(IAA_SM))
               IF(ITOFF.LT.0) THEN
                 ITOFF = -ITOFF
                 ITRNSP = 1
               ELSE
                 ITRNSP = 0
               END IF
             ELSE
               ITOFF = IB_T(ICB_SM,ICA_SM,INVELM(IAB_SM))
               IF(ITOFF.LT.0) THEN
                 ITOFF = -ITOFF
CE               ITRNSP = 1
                 ITRNSP = 0
               ELSE
CE               ITRNSP = 0
                 ITRNSP = 1
               END IF
             END IF
             IF(ITDIAG.EQ.1.AND.ICA_SM.EQ.ICB_SM.AND.
     &          IAA_SM.EQ.IAB_SM) THEN
                IDIAG = 1
             ELSE
                IDIAG = 0
             END IF
*
             IF(ISD.EQ.1) THEN
               if ((IRIS.eq.1.and.IRIC.eq.2).or.
     &             (IRIS.eq.2.and.IRIC.eq.1)) then
                 ITOFF = ITOFF + LEN_T_VEC
               end if
               CALL TCC_SUBBLK(TB(ITOFF),TBSUB,1,
     &                         LICA_SM,NCA_B,ICA_B,LICB_SM,NCB_B,ICB_B,
     &                         LIAA_SM,NAA_B,IAA_B,LIAB_SM,NAB_B,IAB_B,
     &                         ITRNSP,IDIAG,NTEST)
             ELSE
*. Zero density subblock
              ZERO = 0.0D0
              CALL SETVEC(TBSUB,ZERO,NCA_B*NCB_B*NAA_B*NAB_B)
             END IF
*  
* ==========================
*. Obtain C(Ka,Kb,Iaa,Iab)  
* ==========================
*
*. C(Ka,Kb,Iaa,Iab) = sum(Ja,Jb) <Ka!Iaa!Ja><Kb!Iab!Jb> C(Ja,Jb)
             CALL C_TO_CKK(CB(ICBOFF),CJRES,1,
     &                     NJA(JASM),NJB(JBSM),
     &                     NKA_B,NKB_B,NAA_B,NAB_B,
     &                     KJA,XKJA,KJB,XKJB)
*
             IF(ISD.EQ.2) THEN
               CALL C_TO_CKK(SB(ISBOFF),SIRES,1,
     &                       NIA(IASM),NIB(IBSM),
     &                       NKA_B,NKB_B,NCA_B,NCB_B,
     &                       KIA,XKIA,KIB,XKIB)
             END IF
* =================================================================
* For sigma :
* SIRES(Ka,Kb,Ica,Icb) = CJRES(Ka,Kb,Iaa,Iab)*T(Ica,Icb,Iaa,Iab)
* For densi :
* T(Ica,Icb,Iaa,Iab)   = SIRES(Ka,Kb,Ica,Icb)*CJRES(Ka,Kb,Ica,Icb)
* =================================================================
*
             LKAB = NKA_B*NKB_B
             LCAB = NCA_B*NCB_B 
             LAAB = NAA_B*NAB_B
*
             FACTORC = 0.0D0
             FACTORAB = 1.0D0*SIGNXXX
*
             IF(ISD.EQ.1) THEN
               CALL MATML7(SIRES,CJRES,TBSUB,
     &         LKAB,LCAB,LKAB,LAAB,LCAB,LAAB,FACTORC,FACTORAB,2)
*
               IF(NTEST.GE.1000) THEN
                 write(6,*)
                 write(6,*) ' Excitation type ',IXTP
                 write(6,*) ' IRIS,IRIC,ITOFF',IRIS,IRIC,ITOFF
                 WRITE(6,*) ' Updated S(Ka,Kb,*)  from MATML7 '
                 CALL WRTMAT(SIRES,LKAB,LCAB,LKAB,LCAB)
                 write(6,*)
               END IF
             ELSE
               FACTORC = 0.0D0
               CALL MATML7(TBSUB,SIRES,CJRES,
     &         LCAB,LAAB,LKAB,LCAB,LKAB,LAAB,FACTORC,FACTORAB,1) 
               IF(NTEST.GE.100) THEN
                 write(6,*)
                 write(6,*) '========================'
                 write(6,*) ' UPDATED DENSITY        '
                 write(6,*) '========================'
                 write(6,*) ' Excitation type ',IXTP
                 write(6,*) ' IRIS,IRIC,ITOFF',IRIS,IRIC,ITOFF
                 WRITE(6,*) ' Updated TBSUB '
                 CALL WRTMAT(TBSUB,LCAB,LAAB,LCAB,LAAB)
                 write(6,*)
               END IF
             END IF
* =================================================================
* For sigma :
* S(Ia,IB) = SI(Ia,Ib) + <Ka!O+_ca!Ia><Kb!O+cb!Ib>SIRES(Ka,Kb,Ica,Icb)
* For densi :
* Just scatter out to complete density block
* =================================================================
*
             IF(ISD.EQ.1) THEN
               CALL C_TO_CKK(SB(ISBOFF),SIRES,
     &                       2,NIA(IASM),NIB(IBSM),
     &                       NKA_B,NKB_B,NCA_B,NCB_B,
     &                       KIA,XKIA,KIB,XKIB)
               if (NTEST.ge.100) then
                 write(6,*)
                 write(6,*) '========================'
                 write(6,*) ' UPDATED SIGMA VECTOR   '
                 write(6,*) '========================'
                 write(6,*) 'ISBOFF : ',ISBOFF
                 call wrtmat(SB(ISBOFF),NIA(IASM),NIB(IBSM),
     &                                  NIA(IASM),NIB(IBSM))
               end if
             ELSE IF(ITRNSP.EQ.0) THEN
C?             WRITE(6,*) ' Densi-call to TCC, ITOFF  = ', ITOFF
               if ((IRIS.eq.1.and.IRIC.eq.2).or.
     &             (IRIS.eq.2.and.IRIC.eq.1)) then
                 ITOFF = ITOFF + LEN_T_VEC
               end if
               CALL TCC_SUBBLK(TB(ITOFF),TBSUB,2,
     &                  LICA_SM,NCA_B,ICA_B,LICB_SM,NCB_B,ICB_B,
     &                  LIAA_SM,NAA_B,IAA_B,LIAB_SM,NAB_B,IAB_B,
     &                  ITRNSP,IDIAG,NTEST)
             END IF
            END IF
*           ^ End if all maps are nonvanishings
            END DO
            END DO
*           ^ End of loop over batches of creation   operators
            END IF
*           ^ End if batch of beta annihilation was nontrivial
            END DO
            END IF
*           ^ End if batch of alpha annihilation was nontrivial
            END DO
*           ^ End of loop over batches of annihilation operators
           END DO
           END DO
*          ^ End of loop over batches of Ka,Kb
           END IF
*          ^ End if symmetry combinations have nonvanishing dimensions
         END DO
         END DO
*        ^ End of loop over KASM, KBSM
        END IF
*       ^ End if symmetryblock of C should be included
        END DO
*       ^ End of loop over JASM
       END IF
*      ^ End if symmetryblock of S should be included
       END DO
*      ^ End of loop over IASM
      END IF
*     ^ End if Ka and Kb are nontrivial strings
      IF(NTEST.GE.300) THEN 
        IF(ISD.EQ.1) THEN
          WRITE(6,*) ' Updated sigma block '
C              WRTVH1(H,IHSM,NRPSM,NCPSM,NSMOB,ISYM)
          CALL WRTVH1(SB,ISSM,NIA,NIB,NSMST,0)
        ELSE IF (ISD.EQ.2) THEN
          WRITE(6,*) ' Updated density block '
          CALL WRT_TCC_BLK(TB,1,NICA,NICB,NIAA,NIAB,NSMST)
C              WRT_TCC_BLK(TCC,ITCC_SM,NCA,NCB,NAA,NAB,NSMST)
        END IF
      END IF
*
 9999 CONTINUE
*
      CALL QEXIT('GENSIG')
      RETURN
      END
