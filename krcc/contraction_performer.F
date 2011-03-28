************************************************************************
*                                                                      *
* These are primarily the files from Jeppe's NEWCCV code modified to   *
* only perform the contraction of zero or one commutator instead of up *
* to four. Hence simplified.                                           *
*                                                                      *
* Lasse 2011                                                           *
*                                                                      *
************************************************************************
      SUBROUTINE EMTHET_COM_KRCC(T,ITSM,CC_VEC_FNC,NCOMMU,IHM,
     &           NSPOBEX,ISPOBEX_TP,
C next two lines needs to be checked
     &           IBSPOBEX,LSPOBEX,NOBEX,
     &           IBSOX_FOR_OX,NSOX_FOR_OX,ISOX_FOR_OX,
C end the check
     &           NHTP,IHTP,IHINDEX,IHSM,
     &           ILTTWO,I2DIATP,
! Info on contractions
     &           NOP,IOP,
     &           NCONT,IOPRES,
     &           IOPT,
     &           IPOSCON,
     &           IPERM,
     &           IPREFAC,
     &           IPOSCONPERM,
     &           ISTARTOP,IENDOP,
     &           ITOP,NTOP,
     &           IRESOP,NRESOP)!,
C    &           IHOP,IHM)
*
* New version by Lasse 2009. Replaces old version by Jeppe
*
* With precomputed contractions to eliminate overhead
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "crun.inc"
*. Specific input
*. Info on included spinorbital excitations
      INTEGER ISPOBEX_TP(4*NGAS,NSPOBEX+1),LSPOBEX(NSPOBEX+1)
      INTEGER IBSPOBEX(NSPOBEX+1)
*. Orbital excitations to spin-orbital excitations
      INTEGER IBSOX_FOR_OX(*),NSOX_FOR_OX(*),ISOX_FOR_OX(*)
*. Hamilton operator types
      INTEGER IHTP(4*NGAS,NHTP), IHINDEX(4,NHTP), IHSIGN(NHTP)
*. Input CC vector
      DIMENSION T(*)
*. New input
      INTEGER IOP(NOP),IPOSCON(NOP),IPREFAC(4,NOP),IPERM(NOP)
      INTEGER IPOSCONPERM(4,NOP)
      INTEGER IOPRES(NCONT),IOPT(NCONT)
      INTEGER ISTARTOP(NOP),IENDOP(NOP)
      INTEGER ITOP(4*NGAS,NTOP),IRESOP(4*NGAS,NRESOP)
C     INTEGER IHOP(4*NGAS,NOP)
*. Output CC vector function
      DIMENSION CC_VEC_FNC(*)
*
      NTEST = 100
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Information from EMTHET_COM '
        WRITE(6,*) ' ============================='
        WRITE(6,*) ' Input T-amplitudes '
        CALL WRTMAT(T,1,N_CC_AMP,1,N_CC_AMP)
      END IF
*
      IUNIOP = NSPOBEX+1
      IBUNIOP = IBSPOBEX(IUNIOP)

*
      IDUM = -1
      XDUM = -1.0D0
      ONE = 1.0D0
      IONE = 1
*
*. Start by zeroing the CC-vector function
      ZERO = 0.0D0
      CALL SETVEC(CC_VEC_FNC,ZERO,N_CC_AMP+1)
*. New initialization
*. Loop over blocks of Hamilton operator
      DO IH = 1, NOP !NHTP
* First is this a Hamiltonian or an Intermediate operator?
        IF(IHM.EQ.1) THEN
* For the Hamiltonian operators
          IHOP = IOP(IH)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Info for Hamiltonian op number ',IHOP
            WRITE(6,*) ' Absolute address for op ', IHOP
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,IHOP),1)
          END IF
        ELSE IF(IHM.EQ.2) THEN
* For the intermediate operators
          IHOP = IH
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Info for Intermediate op number ',IHOP
            WRITE(6,*) ' Absolute address for op ', IHOP
            CALL WRT_SPOX_TP_CC_KRCC(IHTP(1,IHOP),1)
          END IF
        ELSE 
          WRITE(6,*) ' Lasse forgot to set IHM '
          WRITE(6,*) ' Please do it for him '
          STOP ' Check error message in output '
        END IF
*
C Notice ISPIN is passed down via common block to tell that we now shall
C look for a spin-orbit integral
C
C Not sure to keep or if it will only apply for the Hamiltonian
C operators!
C
        ISPINORBIT = 0
        IF(IHOP.GE.ILTTWO+1) THEN
           ISPINORBIT = 1
        ELSE IF(IHOP.GE.I2DIATP) THEN
           ISPINORBIT = 2
        ELSE
           ISPINORBIT = 0
        END IF
*
*. About signs : signs are all generated in HTYPE_TO_ED_REL, so p.t
*
* But this is causing me problems. So i will try to fix it!
*
        FACTORH = 1.0D0
*
* Zero commutators Still needs fixing!!!!
*
        IF(NCOMMU.EQ.0) THEN
C         IHCTP = NHT_CONN0(IHOP)
          CALL HCT1234_REL(IHTP(1,IHOP),IHINDEX(1,IHOP),0,
     &         ISPOBEX_TP(1,IUNIOP),ISPOBEX_TP(1,IUNIOP),
     &         ISPOBEX_TP(1,IUNIOP),ISPOBEX_TP(1,IUNIOP),
     &         IONE,IONE,IONE,IONE,IHSM,
     &         ONE,ONE,ONE,ONE,
     &         CC_VEC_FNC(IBSPOBEX(IHCTP)),
     &         FACTORH,LCCFBLK)
        END IF
*
* One commutator
*
        IF(NCOMMU.EQ.1) THEN
* Loop over the number of possible contractions
          IF(IPOSCON(IH).GE.1) THEN
            DO I =1,IPERM(IH)
              IF(NTEST.GE.100) THEN
                WRITE(6,*) 'Number of possible contractions ',
     &          IPOSCONPERM(I,IH),' for contraction number ',I
                WRITE(6,*) 'Permutation factor for contraction ',
     &          IPREFAC(I,IH)
              END IF
              IF(IPOSCONPERM(I,IH).GE.1) THEN
                DO J=1,IPOSCONPERM(I,IH)
                  ICOUNT = ICOUNT + 1
                  IF(NTEST.GE.100) THEN
                   WRITE(6,*) ' Operator contracted with ',IOPT(ICOUNT)
                   CALL WRT_SPOX_TP_CC_KRCC(ITOP(1,IOPT(ICOUNT)),1)
                   WRITE(6,*) ' Operators contracted to ',IOPRES(ICOUNT)
                   CALL WRT_SPOX_TP_CC_KRCC(IRESOP(1,IOPRES(ICOUNT)),1)
                  END IF
                  IT1 = IOPT(ICOUNT) !operator contracted with
                  IHCTP = IOPRES(ICOUNT) !resultant operator
                  IT1_B = IBSPOBEX(IT1)
                  LCCFBLK = LSPOBEX(IHCTP)
             print*,'IT1,IHCTP,IT1_B,LCCFBLK',IT1,IHCTP,IT1_B,LCCFBLK
              CALL HCT1234_REL(IHTP(1,IHOP),IHINDEX(1,IHOP),1,
     &             ISPOBEX_TP(1,IT1),
     &             ITSM,IHSM,
     &             T(IT1_B),
     &             CC_VEC_FNC(IBSPOBEX(IHCTP)),
     &             LCCFBLK)
                END DO
              END IF
            END DO
          ELSE
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' No possible contractions '
            END IF
          END IF
        END IF
*
      END DO
*     ^ End of loop over H-operators
*
      RETURN
      END
*
      FUNCTION FAC_IN_BCH(NT,IT1,IT2,IT3,IT4)
*
* Determine factor for 1-4 operators in 
* BCH expansion of commutating operators. Permutation factors 
* are included, so only one term for combination is assumed 
* calculated  
*
* Jeppe Olsen, April 2003
*              Modified (bug for NT = 4) found Nov. 2004
*
#include "implicit.inc"
*
      IF(NT.EQ.0) THEN
        FAC = 1.0D0
      ELSE IF (NT.EQ.1) THEN
        FAC = 1.0D0
      ELSE IF (NT. EQ. 2) THEN
        IF(IT1.NE.IT2) THEN
          FAC = 1.0D0
        ELSE 
          FAC = 0.5D0
        END IF
      ELSE IF (NT. EQ. 3) THEN
        IF(IT1.NE.IT2.AND.IT2.NE.IT3.AND.IT1.NE.IT3) THEN
          FAC = 1.0D0
        ELSE IF (IT1.EQ.IT2.AND.IT2.EQ.IT3) THEN 
          FAC = 1.0D0/6.0D0
        ELSE 
          FAC = 1.0D0/2.0D0
        END IF
      ELSE IF (NT.EQ.4) THEN
        IF(IT1.NE.IT2.AND.IT1.NE.IT3.AND.IT1.NE.IT4.AND.
     &     IT2.NE.IT3.AND.IT2.NE.IT4.AND.IT3.NE.IT4     )THEN
*. Four different
          FAC = 1.0D0
        ELSE IF(IT1.EQ.IT2.AND.IT2.EQ.IT3.AND.IT3.EQ.IT4) THEN
*. Four identical 
          FAC = 1.0D0/24.0D0
        ELSE IF((IT1.EQ.IT2.AND.IT2.EQ.IT3) .OR. 
     &          (IT1.EQ.IT2.AND.IT2.EQ.IT4) .OR.
     &          (IT1.EQ.IT3.AND.IT3.EQ.IT4) .OR.
     &          (IT2.EQ.IT3.AND.IT3.EQ.IT4)     ) THEN
*. Three identical 
         FAC =  1.0D0/6.0D0
        ELSE IF ((IT1.EQ.IT2.AND.IT3.EQ.IT4) .OR. 
     &           (IT1.EQ.IT3.AND.IT2.EQ.IT4) .OR.
     &           (IT1.EQ.IT4.AND.IT2.EQ.IT3)    )  THEN
*. Two pairs of identical numbers pairs 
         FAC = 1.0D0/4.0D0
        ELSE 
*. Two identical numbers, the last two differs 
C        FAC = 1.0D0/6.0D0
         FAC = 1.0D0/2.0D0
        END IF
      END IF
*     ^ End of switch between different number of operators
      FAC_IN_BCH = FAC
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        IF(NT.EQ.0) THEN 
          WRITE(6,*) ' No T-ops, FAC == ', FAC
        ELSE IF (NT.EQ.1) THEN
          WRITE(6,*) ' One T-op, FAC = ', FAC
        ELSE IF (NT.EQ.2) THEN
          WRITE(6,*) ' T1, T2 , FAC = ', IT1,IT2, FAC
        ELSE IF (NT.EQ.3) THEN
          WRITE(6,*) 'T1, T2, T3, FAC = ', IT1,IT2,IT3, FAC
        ELSE IF (NT.EQ.4) THEN
          WRITE(6,*) ' T1, T2, T3, T4, FAC = ',  
     &                IT1,IT2,IT3,IT4,FAC
        END IF
      END IF
*
      RETURN
      END 
*
      SUBROUTINE HCT1234_REL(IHCAAB,IHINDEX,NTOP,IT1,
     &                   IT1SM,IHSM,T1,
     &                   HT1234,LHT1234)
* 
* All possible contractions of a given type of the Hamiltonian 
* with upto four T-operators of given type. 
*
*
* Jeppe Olsen, April 2003
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
*
*=========
*. Input 
*=========
*
*.. H-operator in CAAB form
      INTEGER IHCAAB(4*NGAS),IHINDEX(4)
*.. T-operators : type and coefficients
      INTEGER IT1(4*NGAS)!,IT2(4*NGAS),IT3(4*NGAS),IT4(4*NGAS)
      DIMENSION T1(*)!, T2(*), T3(*), T4(*), HT1234(*)
*
*. Local scratch
**. Assuming that H contains atmost 2 particle operator 
      PARAMETER(MXPCONN = 4*4*4*4)
      INTEGER ICONN(4,MXPCONN)
      INTEGER IHEXP(3,4)
C     INTEGER ID4(4*MXPNGAS),ID2(4*MXPNGAS),ID3(4*MXPNGAS)
      INTEGER ID1(4*MXPNGAS), IEX(4*MXPNGAS)
      INTEGER IORD(4),IEXD1234(4)
C     INTEGER IT4_REO(4*MXPNGAS),IT2_REO(4*MXPNGAS),IT3_REO(4*MXPNGAS)
C     INTEGER IT1_REO(4*MXPNGAS)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' ======================='
        WRITE(6,*) ' Welcome to HCT1234_REL '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' H-operator in action (CAAB)'
        CALL WRT_SPOX_TP_CC_KRCC(IHCAAB,1)
        WRITE(6,*) ' IT1SM ',IT1SM
        IF(NTOP.GE.1) THEN
          WRITE(6,*) ' T operator '
          CALL WRT_SPOX_TP_CC_KRCC(IT1,1)
        END IF
        IF (NTEST.GE.1000) THEN
          WRITE(6,*) ' HCT1234 : Input block of CC vector function '
          CALL WRTMAT(HT1234,1,LHT1234,1,LHT1234)
        END IF
      END IF
*. Reform H operator to expanded form - assumed in CONTR_POS
C     REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,IDOREO,IREO)
C     CALL REFORM_HTYP(IHCAAB,IHEXP,NHOP,1,NGAS,1,IHINDEX)
*. Obtain in CAAB form 
C May still be used Lasse
      CALL REFORM_HTYP_REL(IHCAAB,IHEXP,NHOP,1,NGAS,0,IHINDEX)
C      stop 'after REFORM_HTYP_REL'
*. Find all possible contractions of H operator with the NTOP T operators
C This is already found on the outside so we know which contraction is
C supposed be performed so no need to find it! 
C     CALL CONTR_POS_REL(IHEXP,NHOP,NTOP,IT1,IT2,IT3,IT4,
C    &               IONLY_CONN,NCONN,ICONN)
C Will here first for convenience replace ICONN array
      CALL ISETVC(ICONN(1,1),NTOP,4)
C No loop needed since it is one well defined contraction comming down.
C     DO JCONN = 1, NCONN
*. Obtain contraction order, and the individual deexcitation operators
* Will keep the sign here! At least for now.
      CALL HTYPE_TO_ED_REL(IHCAAB,IHINDEX,NHOP,ICONN(1,1),NTOP,
     &                   ID1,ID2,ID3,ID4,IEX,IEXD1234,IORD,ISIGN_DE)
C     stop 'after HTYPE_TO_ED_REL'
*. Reorder t-operators so they correspond to new order
*. I will use 4 operators, as the last (inactive have just been 
*. defined to be unit operator. This is done as OPCT1234 
*. currently requires four ops and four symmetries
* No need to reorder operators since there is only one anyway!
* And they come in the way they should be contracted!
*
C        ISIGN_DE = ISIGN_DE * FACX
C        print*,'ISIGN_DE * FACX',ISIGN_DE,FACX
C     CALL REO_CAABS_REL(IT1,IT2,IT3,IT4,4,
C    &                 IT1_REO,IT2_REO,IT3_REO,IT4_REO,
C    &                 IT1SM,IT2SM,IT3SM,IT4SM,
C    &                 IT1SM_REO,IT2SM_REO,IT3SM_REO,IT4SM_REO,
C    &                 IORD, NGAS)
*.   And then do the contraction 
*
* There are 24 different ways of making the contractions so, I must
* have 24 cases, and insert the appropriate T-operators, 
* If I do this without problems, I will give myself a beer ..
* 
* No need for beer if you only have one way. Lasse
*
*.1234 Always keep the same order now!!! 
* Commented out no longer needed things (good vocabulary today)
      CALL OPCT1234M_REL(IEX,ID1,!ID2,ID3,ID4,
     &                   IT1,!IT1_REO!,IT2_REO,IT3_REO,IT4_REO,
     &                   T1,!T2,T3,T4,
     &                   HT1234,
     &                   IT1SM_REO,!IT2SM_REO,IT3SM_REO,IT4SM_REO,
     &                   IHSM,IEXD1234,ISIGN_DE)
C     OPCT1234M(IOEX,IO1DX,IO2DX,IO3DX,IO4DX,
C    &                    IT1,IT2,IT3,IT4,T1,T2,T3,T4,OT1234,
C    &                    IT1SM,IT2SM,IT3SM,IT4SM,IOPSM,IOPINDX,
C    &                    ISIGN,FACX)
C      END DO
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' HCT1234 : Output block of CC vector function '
         CALL WRTMAT(HT1234,1,LHT1234,1,LHT1234)
      END IF
      NTEST = 0
*
      RETURN
      END 
*
      SUBROUTINE REFORM_HTYP_REL(IHCAAB,IHOP,NOP,IWAY,NGAS,IDOREO,IREO)
*
* Reform different forms of Hamiltonian operator 
*
* IWAY = 1 : IH_CAAB,IHINDX => IHOP
*      = 2 : IH_CAAB <= IHOP
* 
* If IDOREO = 1, reordering of indeces is done according to IREO
*
*. Note : Operators are ordered in CA CB AA AB order
*
* Jeppe Olsen, July 16, 2001
*
#include "implicit.inc"
*. Input
      INTEGER IREO(*)
*. Input and otput
      INTEGER IHCAAB(NGAS,4), IHOP(4,3)
*. Local scratch
      INTEGER JHOP(4,3)
*
      IF(IWAY.EQ.1) THEN
* HCAAB => HOP
       LOP = 0
       ICAAB = 0
       DO ICA = 1, 2
        DO IAB = 1, 2
          ICAAB = (ICA-1)*2 + IAB
          DO IGAS = 1, NGAS
            L = IHCAAB(IGAS,ICAAB)
            DO IOP = 1, L
              LOP = LOP + 1
              JHOP(LOP,1) = IGAS
              JHOP(LOP,2) = ICA
              JHOP(LOP,3) = IAB
            END DO
          END DO
        END DO
       END DO
       NOP = LOP
*. Reordering or transfer
       DO IOP = 1, NOP
         IF(IDOREO.EQ.0) THEN
           IOPEFF = IOP
         ELSE
           IOPEFF = IREO(IOP)
         END IF
         DO IND = 1, 3
           IHOP(IOP,IND) = JHOP(IOPEFF,IND)
         END DO
       END DO
      ELSE 
*. HOP => HCAAB
*. Nop must here be input
       IZERO = 0
       CALL ISETVC(IHCAAB,IZERO,4*NGAS)
       DO IOP = 1, NOP
         IGAS = IHOP(IOP,1)
         ICA  = IHOP(IOP,2)
         IAB  = IHOP(IOP,3)
         ICAAB = (ICA-1)*2+IAB
         IHCAAB(IGAS,ICAAB) = IHCAAB(IGAS,ICAAB) + 1
       END DO
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reform of Htype '
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' CAAB => OP reform '
        ELSE
         WRITE(6,*) ' OP => CAAB reform'
        END IF
        IF(IDOREO.NE.0) THEN
          WRITE(6,*) ' The reorder array '
          CALL IWRTMA(IREO,1,NOP,1,NOP)
        END IF 
        WRITE(6,*) ' CAAB form '
        CALL WRT_SPOX_TP_CC_KRCC(IHCAAB,1)
        WRITE(6,*) ' Operator form '
C            WRT_CNTR(ICONT,NCONT,LDIM)
        CALL WRT_CNTR3(IHOP,NOP,4)
      END IF
*
      RETURN
      END
*
      SUBROUTINE CONTR_POS_REL(IHOP,NHOP,NTOP,IT1,IT2,IT3,IT4,
     &                     IONLY_CONN,NCONN,ICONN)
*
* A general operator IHOP containing NHOP operators are given.
*.Find possible contractions with the NTOP excitation 
*.operators IT1, IT2, ..
* Only connections where all  deexcitation operators in IHOP are 
* contracted with excitation operators are included
*
*. IF IONLY_CONN = 1, then connected terms are included
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      INTEGER IHOP(4,3)
*. GAS : 1, ICA : 2, IAB : 3
      INTEGER IT1(4*NGAS),IT2(4*NGAS),IT3(4*NGAS),IT4(4*NGAS)
*. Output
      INTEGER ICONN(4,*)
*. ICONN(JOP,JCONN) : Which T-operator is operator JOP in IHOP
*.                    connected with in connection JCONN
*. ICONN(JOP,JCONN) = 0 => No connection of operator JOP in IHOP
*
*. Local scratch 
      INTEGER ISCR(4), IMAXVAL(4)
      INTEGER IT(4*MXPNGAS,4)
      INTEGER JCONN(4)
      INTEGER IDE(4)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input to CONTR_POS_REL' 
        WRITE(6,*) ' Form of general operator : '
        CALL WRT_CNTR3(IHOP,NHOP,4)
        IF(NTOP.GE.1) CALL WRT_SPOX_TP_CC_KRCC(IT1,1)
        IF(NTOP.GE.2) CALL WRT_SPOX_TP_CC_KRCC(IT2,1)
        IF(NTOP.GE.3) CALL WRT_SPOX_TP_CC_KRCC(IT3,1)               
        IF(NTOP.GE.4) CALL WRT_SPOX_TP_CC_KRCC(IT4,1)
      END IF
*. Check if some of the operators are equivalent, i.e. 
*. have identical CA, AB, IGAS. 
*. Currently I am only checking whether operators 1 and 2 are 
*. identical and whether operators 3 and 4 are identical 
      N12EQV = 1
      DO ICOMP = 1, 3
        IF(IHOP(1,ICOMP).NE.IHOP(2,ICOMP)) N12EQV = 0
      END DO
*
      N34EQV = 0
      IF(NHOP.EQ.4) THEN
        N34EQV = 1
        DO ICOMP = 1, 3
          IF(IHOP(3,ICOMP).NE.IHOP(4,ICOMP)) N34EQV = 0
        END DO
      END IF
*. Use IPHGAS1 to classify operators as excitation or deexcitation operators
*. Should be modified for open shell reference 
*
*  Could this not just be changed to IPHGAS_AB ? Lasse. + more changes
*
      DO JOP = 1, NHOP
        IF((IPHGAS1(IHOP(JOP,1)).EQ.1.AND.IHOP(JOP,2).EQ.2).OR.
     &     (IPHGAS1(IHOP(JOP,1)).EQ.2.AND.IHOP(JOP,2).EQ.1) ) THEN
           IDE(JOP) = 1
        ELSE
           IDE(JOP) = 2
        END IF
      END DO
C     WRITE(6,*) ' The IDE array ', (IDE(I),I=1,NHOP)
*. A given connection is defined by NHOP integers between 1 and 
*. NTOP. Loop over these integers
      CALL ISETVC(IMAXVAL,NTOP,NHOP)
      IFIRST = 1
      IZERO = 0
      NCONN = 0
 1000 CONTINUE
        IF(IFIRST.EQ.1) THEN
          CALL ISETVC(ISCR,IZERO,NHOP)
          NONEW = 0
          IFIRST = 0
        ELSE
C              NXTNUM2(INUM,NELMNT,MINVAL,MAXVAL,NONEW)
          CALL NXTNUM2(ISCR,NHOP,IZERO,IMAXVAL,NONEW)
        END IF
C       WRITE(6,*) ' Next suggested  connection '
C       CALL IWRTMA(ISCR,1,NHOP,1,NHOP)
        IF(NONEW.EQ.0) THEN
*. Check if this connection is nonvanishing
          CALL ICOPVE(IT1,IT(1,1),4*NGAS)
          IF(NTOP.GE.2) CALL ICOPVE(IT2,IT(1,2),4*NGAS)
          IF(NTOP.GE.3) CALL ICOPVE(IT3,IT(1,3),4*NGAS)
          IF(NTOP.GE.4) CALL ICOPVE(IT4,IT(1,4),4*NGAS)
*
          I_AM_OKAY = 1
          DO IOP = 1, NHOP
            JTOP = ISCR(IOP)
            JGAS = IHOP(IOP,1)
            JCA =  IHOP(IOP,2)
            JAB =  IHOP(IOP,3)
*. Well, a creation operator in IHOP should connect with 
*. an annihilation operator in IT and vice versa so
            IF(JCA.EQ.1) THEN
             JTCA = 2
            ELSE
             JTCA = 1
            END IF
*
            JCAAB = (JTCA-1)*2+JAB
            J_AD = JGAS + (JCAAB-1)*NGAS
            IF(JTOP.GT.0) THEN
              IT(J_AD,JTOP) = IT(J_AD,JTOP)-1
              IF(IT(J_AD,JTOP).LT.0)  I_AM_OKAY = 0
            END IF
          END DO
C         WRITE(6,*) ' I_AM_OKAY for this connection ', I_AM_OKAY
*. If required, check whether term is connected
          IF(IONLY_CONN.EQ.1.AND.I_AM_OKAY.EQ.1) THEN 
            CALL ISETVC(JCONN,IZERO,NTOP)
            DO JHOP = 1, NHOP
              IF(ISCR(JHOP).NE.0) JCONN(ISCR(JHOP)) = 1
            END DO
            ICONNECTED = 1
            DO JTOP = 1, NTOP
              IF(JCONN(JTOP).EQ.0) ICONNECTED = 0
            END DO
          ELSE
            ICONNECTED = 1
          END IF
C         WRITE(6,*) ' ICONNECTED for this connection ', ICONNECTED
*         ^ End of check for connectedness is required
*. Check to see ensure that all deexcitation operators are 
*. connected
          DO JHOP = 1, NHOP
            IF(IDE(JHOP).EQ.1.AND.ISCR(JHOP).EQ.0) 
     &      I_AM_OKAY = 0    
C           PRINT*,'IDE(JHOP),ISCR(JHOP)',IDE(JHOP),ISCR(JHOP)   
          END DO
C         WRITE(6,*) ' Check , All deexcitations connected ', I_AM_OKAY
*. Check that for equivalent operators the contracted operators 
*. occur in ascending order. This is to ensure that 
*. contractions splitting equivalent operators only 
*. are included once
          IF(N12EQV.EQ.1.AND.ISCR(1).LT.ISCR(2)) I_AM_OKAY = 0
C         PRINT*,'I_AM_OKAY AFTER N12EQV',I_AM_OKAY
          IF(NHOP.EQ.4) THEN
            IF(N34EQV.EQ.1.AND.ISCR(3).LT.ISCR(4)) I_AM_OKAY = 0
          END IF
*
C         print*,'I_AM_OKAY,ICONNECTED',I_AM_OKAY,ICONNECTED
          IF(I_AM_OKAY.EQ.1.AND.ICONNECTED.EQ.1) THEN
*. Well a new connection has been born 
             NCONN = NCONN + 1
             CALL ICOPVE(ISCR,ICONN(1,NCONN),NHOP)
          END IF
      GOTO 1000
        END IF
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Output from CONTR_POS_REL '
*
        WRITE(6,*) ' Number of obtained connections = ', NCONN
        WRITE(6,*) ' Obtained connections : '
        CALL IWRTMA(ICONN,NHOP,NCONN,4,NCONN)
      END IF
*
      RETURN
      END
*
      SUBROUTINE HTYPE_TO_ED_REL(IHTYPE,IHINDEX,NHOP,ICONT,
     &           NTOP,ID1,ID2,ID3,ID4,IEX,IEXD1234,IORD,ISIGN_DE)
*    
*. A Hamilton operator (IHTYPE) is given as a set of C/A operators  
*. (in CA CB AA AB order), and a ordering array IHINDEX.
*. A contraction is given (ICONT), telling which T-operator 
*. each operator in H should be connected to.
*. Obtain 1 : The order in which the contraction should be done 
*.        2 : The individual excitations and deexcitations in CAAB form
*.        3 : Sign for bringing CA CB AA AB operator into EX D1 D2 D3 D4
*.            form
*
* Sign for bringing CA CB AA AB operator into EX D1 D2 D3 D4 form
*
* Jeppe Olsen, July 16, 2001, Finished Easter 2003 ( April)
*              Debugging initiated Nov. 2004 
*. Some phase modifications, Nov. 2004 ... 
*  ( I understand why somebody asks higher powers for signs,
*    they are difficult to get yourself ...)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
*. Input
      INTEGER IHTYPE(NGAS,4), IHINDEX(4),ICONT(4)
*. Output 
      INTEGER ID1(NGAS,4),ID2(NGAS,4),ID3(NGAS,4),ID4(NGAS,4)
      INTEGER IEX(NGAS,4),IEXD1234(4)
      INTEGER IORD(4)
* IORD(IOP) : New operator IOP is old operator IORD(IOP) ( new => old)
*. Local scratch 
      INTEGER IHOP(4,3)
      INTEGER IDOP(4,3,4),IEXOP(4,3), LD(4) 
C     INTEGER IT_AR(4)
      INTEGER IEX_INDEX(4),ID_INDEX(4,4)
*
      IZERO = 0
*
*. Reform IHTYPE from CAAB form to a string of NHOP operators
C     REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,IDOREO,IREO)
C     CALL REFORM_HTYP(IHTYPE,IHOP,NHOP,1,NGAS,1,IHINDEX)
C     CALL REFORM_HTYP(IHTYPE,IHOP,NHOP,1,NGAS,1,IHINDEX)
*. Keep CAAB ordering ..
      CALL REFORM_HTYP_REL(IHTYPE,IHOP,NHOP,1,NGAS,0,IHINDEX)
*. Obtain the deexcitations/excitations D1,D2,D3,D4,EX (they will now
*. be in CAAB order)
      DO I = 1, 4
        LD(I) = 0
      END DO
      NEX = 0
      DO IOP = 1, NHOP
        INDEX = ICONT(IOP)
*.      ^ Index tells to which T-operator an operator is connected, 0=>disc.
        IF(INDEX.EQ.0) THEN
*. Excitation part
          NEX = NEX + 1
          IEXOP(NEX,1) = IHOP(IOP,1)
          IEXOP(NEX,2) = IHOP(IOP,2)
          IEXOP(NEX,3) = IHOP(IOP,3)
C*. Original index for this operator
C         IEX_INDEX(NEX) = IHINDEX(IOP)
*. index for this operator in IHOP
          IEX_INDEX(NEX) = IOP
        ELSE 
*.Deexcitation part
          LD(INDEX) = LD(INDEX) + 1
          I = LD(INDEX)
          IDOP(I,1,INDEX) = IHOP(IOP,1)
          IDOP(I,2,INDEX) = IHOP(IOP,2)
          IDOP(I,3,INDEX) = IHOP(IOP,3)
C         ID_INDEX(I,INDEX) = IHINDEX(IOP)
          ID_INDEX(I,INDEX) = IOP
        END IF
      END DO
*. Obtain the execution order
C     CONTR_ORD4(NDOP,ID1,ID2,ID3,ID4,
C    &           LD1,LD2,LD3,LD4,IORD) 
       CALL CONTR_ORD4_REL(NTOP,IDOP(1,1,1),IDOP(1,1,2),
     &                 IDOP(1,1,3),IDOP(1,1,4),
     &                 LD(1),LD(2),LD(3),LD(4),IORD)
*. I want IORD to contain 4 well-defined terms, also 
*. when there are less than 4 operators so
       IF(NTOP.LT.4) IORD(4) = 4
       IF(NTOP.LT.3) IORD(3) = 3
       IF(NTOP.LT.2) IORD(2) = 2
       IF(NTOP.LT.1) IORD(1) = 1
*...   (^I guess I do not need LT.0)
*. On return IORD : New => Old order
*. Set up IEXD1234 : Index in operator EXD1D2D3D4 to original order
      DO JEX = 1, NEX
        IEXD1234(JEX) = IEX_INDEX(JEX)
C        print*,'IEX_INDEX(JEX)',IEX_INDEX(JEX)
      END DO
      LEN = NEX
*. Loop over deexcitations in NEW order
      DO JDOP = 1, NTOP
        JDOP_OLD = IORD(JDOP)
C?      WRITE(6,*) 'JDOP, JDOP_OLD', JDOP, JDOP_OLD
C?      WRITE(6,*) ' LD(JDOP_OLD) ',  LD(JDOP_OLD)
        DO JOP = 1, LD(JDOP_OLD)
          LEN = LEN + 1
          IEXD1234(LEN) = ID_INDEX(JOP,JDOP_OLD)
C          print*,'ID_INDEX(JOP,JDOP_OLD)',ID_INDEX(JOP,JDOP_OLD)
        END DO
      END DO
*. IEXD1234 refers to indeces in IHOP, modify so they 
*. correspond to indeces in the original operator
      DO JOP = 1, NHOP
*. Set up the deexcitations numbered according to execution order 
        IEXD1234(JOP) = IHINDEX(IEXD1234(JOP))
C        print*,'IHINDEX(IEXD1234(JOP))',IHINDEX(IEXD1234(JOP))
C        print*,'IHINDEX',IHINDEX(JOP)
      END DO
*. Excitation in CAAB form
      IDUM = 0
      CALL REFORM_HTYP_REL(IEX,IEXOP,NEX,2,NGAS,0,IDUM)
C          REFORM_HTYP(IHCAAB,IHOP,NOP,IWAY,NGAS,0,IDUM)
*. First operator to be contracted
      CALL REFORM_HTYP_REL(ID1,IDOP(1,1,IORD(1)),LD(IORD(1)),2,
     &                     NGAS,0,IDUM)
*. Second op
      IF(NTOP.GE.2) THEN 
       CALL REFORM_HTYP_REL(ID2,IDOP(1,1,IORD(2)),LD(IORD(2)),2,
     &                      NGAS,0,IDUM)
      ELSE
       CALL ISETVC(ID2,IZERO,4*NGAS)
      END IF
*. third  op
      IF(NTOP.GE.3) THEN 
       CALL REFORM_HTYP_REL(ID3,IDOP(1,1,IORD(3)),LD(IORD(3)),2,
     &                      NGAS,0,IDUM)
      ELSE
       CALL ISETVC(ID3,IZERO,4*NGAS)
      END IF
*. fourth  op
      IF(NTOP.GE.4) THEN
       CALL REFORM_HTYP_REL(ID4,IDOP(1,1,IORD(4)),LD(IORD(4)),2,
     &                      NGAS,0,IDUM)
      ELSE
       CALL ISETVC(ID4,IZERO,4*NGAS)
      END IF
*. Sign for bringing H operator from CAAB order to EXD1234 order
C IPERM_PARITY(IPERM,NELMNT)
C      CALL IWRTMA(IEXD1234,1,NHOP,1,NHOP)
      ISIGN_DE = IPERM_PARITY(IEXD1234,NHOP)
C      CALL IWRTMA(IEXD1234,1,NHOP,1,NHOP)
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Output from HTYPE_TO_ED_REL ' 
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Contraction of ops. (ICONT) '
        CALL IWRTMA(ICONT,1,NHOP,1,4)
        WRITE(6,*) ' Suggested order of contraction : '
        CALL IWRTMA(IORD,1,NTOP,1,NTOP)
        WRITE(6,*) ' IEX,D1,D2,D3,D4 : '
        CALL WRT_SPOX_TP_CC_KRCC(IEX,1)
        CALL WRT_SPOX_TP_CC_KRCC(ID1,1)
        IF(NTOP.GE.2) 
     &  CALL WRT_SPOX_TP_CC_KRCC(ID2,1)
        IF(NTOP.GE.3) 
     &  CALL WRT_SPOX_TP_CC_KRCC(ID3,1)
        IF(NTOP.GE.4) 
     &  CALL WRT_SPOX_TP_CC_KRCC(ID4,1) 
        WRITE(6,*) ' IEXD1234 : new to old index array '
        CALL IWRTMA(IEXD1234,1,NHOP,1,NHOP)
        WRITE(6,*) ' ISIGN_DE = ',ISIGN_DE
      END IF
*
      RETURN
      END
*
      SUBROUTINE REO_CAABS_REL(IN1,IN2,IN3,IN4,NIN,
     &                     IOUT1,IOUT2,IOUT3,IOUT4,
     &                     INSM1,INSM2,INSM3,INSM4,
     &                     IOUTSM1,IOUTSM2,IOUTSM3,IOUTSM4,
     &                     IREO, NGAS)
*
* Reorder 0-4 CAAB operators
*
* Jeppe Olsen, April 2003
*
#include "implicit.inc"
#include "mxpdim.inc"
*. Input
      INTEGER IN1(*),IN2(*),IN3(*),IN4(*)
*. New => old ( reordered => original)
      INTEGER IREO(*)
*. Output 
      INTEGER IOUT1(*),IOUT2(*),IOUT3(*),IOUT4(*)
*. Local scratch
      INTEGER ICAAB(4*MXPNGAS)
*
      ISMX = -1
      DO IOP = 1, NIN
        IF(IOP.EQ.1) THEN 
           CALL ICOPVE(IN1,ICAAB,4*NGAS)
           ISMX = INSM1
        END IF
        IF(IOP.EQ.2) THEN 
           CALL ICOPVE(IN2,ICAAB,4*NGAS)
           ISMX = INSM2
        END IF
        IF(IOP.EQ.3) THEN 
           CALL ICOPVE(IN3,ICAAB,4*NGAS)
           ISMX = INSM3
        END IF
        IF(IOP.EQ.4) THEN 
           CALL ICOPVE(IN4,ICAAB,4*NGAS)
           ISMX = INSM4
        END IF
*. What operator is IOP going into ?
*. Compiler warnings ...
        IOP_REO = -2206
        DO JOP = 1, NIN
          IF(IREO(JOP).EQ.IOP) IOP_REO = JOP
        END DO
        IF(IOP_REO.EQ.1) CALL ICOPVE(ICAAB,IOUT1,4*NGAS)
        IF(IOP_REO.EQ.2) CALL ICOPVE(ICAAB,IOUT2,4*NGAS)
        IF(IOP_REO.EQ.3) CALL ICOPVE(ICAAB,IOUT3,4*NGAS)
        IF(IOP_REO.EQ.4) CALL ICOPVE(ICAAB,IOUT4,4*NGAS)
*
        IF(IOP_REO.EQ.1) IOUTSM1 = ISMX
        IF(IOP_REO.EQ.2) IOUTSM2 = ISMX
        IF(IOP_REO.EQ.3) IOUTSM3 = ISMX
        IF(IOP_REO.EQ.4) IOUTSM4 = ISMX
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Reordering of T-operators '
        WRITE(6,*) ' Reorder array, new => old '
        CALL IWRTMA(IREO,1,NIN,1,NIN)
        WRITE(6,*)
        WRITE(6,*) ' Input T-operators : '
        WRITE(6,*)
        IF(NIN .GE.1) THEN
          WRITE(6,*) ' T1_IN : '
          CALL WRT_SPOX_TP_CC_KRCC(IN1,1)
        END IF
        IF(NIN .GE.2) THEN
          WRITE(6,*) ' T2_IN : '
          CALL WRT_SPOX_TP_CC_KRCC(IN2,1)
        END IF
        IF(NIN .GE.3) THEN
          WRITE(6,*) ' T3_IN : '
          CALL WRT_SPOX_TP_CC_KRCC(IN3,1)
        END IF
        IF(NIN .GE.4) THEN
          WRITE(6,*) ' T4_IN : '
          CALL WRT_SPOX_TP_CC_KRCC(IN4,1)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' Output T-operators : '
        WRITE(6,*)
        IF(NIN .GE.1) THEN
          WRITE(6,*) ' T1_OUT : '
          CALL WRT_SPOX_TP_CC_KRCC(IOUT1,1)
        END IF
        IF(NIN .GE.2) THEN
          WRITE(6,*) ' T2_OUT : '
          CALL WRT_SPOX_TP_CC_KRCC(IOUT2,1)
        END IF
        IF(NIN .GE.3) THEN
          WRITE(6,*) ' T3_OUT : '
          CALL WRT_SPOX_TP_CC_KRCC(IOUT3,1)
        END IF
        IF(NIN .GE.4) THEN
          WRITE(6,*) ' T4_OUT : '
          CALL WRT_SPOX_TP_CC_KRCC(IOUT4,1)
        END IF
        WRITE(6,*) ' Output symmetries : '
        IF(NIN.EQ.1) WRITE(6,*) IOUTSM1
        IF(NIN.EQ.2) WRITE(6,*) IOUTSM1, IOUTSM2
        IF(NIN.EQ.3) WRITE(6,*) IOUTSM1, IOUTSM2, IOUTSM3
        IF(NIN.EQ.4) WRITE(6,*) IOUTSM1, IOUTSM2, IOUTSM3, IOUTSM4
      END IF
*
      RETURN
      END
*
      SUBROUTINE OPCT1234M_REL(IOEX,IO1DX,!IO2DX,IO3DX,IO4DX,
     &                    IT1,!IT2,IT3,IT4,
     &                    T1,!T2,T3,T4,
     &                    OT1234,
     &                    IT1SM,!IT2SM,IT3SM,IT4SM,
     &                    IOPSM,IOPINDX,ISIGN)
*
* Master routine for contracting operator with upto 4 operators
*
* Jeppe Olsen, April 2003 ( I hope it will be working before May ..)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
C and the local scratch 
C     COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
C    &               KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTOCC7,
C    &               KLSTOCC8,
C    &               KLSTREO,
C    &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
C    &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
C    &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
C    &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
C    &               KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
C    &               KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
C    &               KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
C    &               KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
C    &               KIX5_CA,KSX5_CA,KIX5_CB,KSX5_CB,
C    &               KIX5_AA,KSX5_AA,KIX5_AB,KSX5_AB,
C    &               KIX6_CA,KSX6_CA,KIX6_CB,KSX6_CB,
C    &               KIX6_AA,KSX6_AA,KIX6_AB,KSX6_AB,
C    &               KIX7_CA,KSX7_CA,KIX7_CB,KSX7_CB,
C    &               KIX7_AA,KSX7_AA,KIX7_AB,KSX7_AB,
C    &               KIX8_CA,KSX8_CA,KIX8_CB,KSX8_CB,
C    &               KIX8_AA,KSX8_AA,KIX8_AB,KSX8_AB,
C    &               KLTSCR1,KLTSCR2,KLTSCR3,
C    &               KLTSCR4,KLTSCR5,KLTSCR6,
C    &               KLOPSCR,
C    &               KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
C    &               KLSMD1,KLSMD2,KLSMD3,KLSMD4,
C    &               KLSMEX,KLSMK1,KLSMK2,KLSML1,
C    &               KLNMD1,KLNMD2,KLNMD3,KLNMD4,
C    &               KLNMEX,KLNMK1,KLNMK2,KLNML1,
C    &               KLOCK1, KLOCK2, KLOCK3,KLOCK4, 
C    &               KLOCL1, KLOCL2, KLOCL3,KLOCL4, KL_IBF, 
C    &               KLEXEORD
C from MEM
      COMMON/CC_SCR2/KLZ,KLZSCR,KLSTOCC1,KLSTOCC2,
     &               KLSTOCC3,KLSTOCC4,KLSTOCC5,KLSTOCC6,KLSTOCC7,
     &               KLSTOCC8,KLSTREO,
     &               KIX1_CA,KSX1_CA,KIX1_CB,KSX1_CB,
     &               KIX1_AA,KSX1_AA,KIX1_AB,KSX1_AB,
     &               KIX2_CA,KSX2_CA,KIX2_CB,KSX2_CB,
     &               KIX2_AA,KSX2_AA,KIX2_AB,KSX2_AB,
     &               KIX3_CA,KSX3_CA,KIX3_CB,KSX3_CB,
     &               KIX3_AA,KSX3_AA,KIX3_AB,KSX3_AB,
     &               KIX4_CA,KSX4_CA,KIX4_CB,KSX4_CB,
     &               KIX4_AA,KSX4_AA,KIX4_AB,KSX4_AB,
     &               KIX5_CA,KSX5_CA,KIX5_CB,KSX5_CB,
     &               KIX5_AA,KSX5_AA,KIX5_AB,KSX5_AB,
     &               KIX6_CA,KSX6_CA,KIX6_CB,KSX6_CB,
     &               KIX6_AA,KSX6_AA,KIX6_AB,KSX6_AB,
     &               KIX7_CA,KSX7_CA,KIX7_CB,KSX7_CB,
     &               KIX7_AA,KSX7_AA,KIX7_AB,KSX7_AB,
     &               KIX8_CA,KSX8_CA,KIX8_CB,KSX8_CB,
     &               KIX8_AA,KSX8_AA,KIX8_AB,KSX8_AB,
     &               KLTSCR1,KLTSCR2,KLTSCR3,
Classe     &               KLTSCR4,KLTSCR5,KLTSCR6,
     &               KLOPSCR,
     &               KLIOD1_ST,KLIOD2_ST,KLIOEX_ST,
     &               KLSMD1,KLSMD2,KLSMD3,KLSMD4,
     &               KLSMEX,KLSMK1,KLSMK2,KLSML1,
     &               KLNMD1,KLNMD2,KLNMD3,KLNMD4,
     &               KLNMEX,KLNMK1,KLNMK2,KLNML1,
     &               KLOCK1, KLOCK2, KLOCK3,KLOCK4, 
     &               KLOCL1, KLOCL2, KLOCL3,KLOCL4, KL_IBF, 
     &               KLEXEORD
      CALL OPCT1234_REL(IOEX,IO1DX,!IO2DX,IO3DX,IO4DX,
     &           IT1,!IT2,IT3,IT4,
     &           T1,!T2,T3,T4,
     &           OT1234,
     &           IT1SM,!IT2SM,IT3SM,IT4SM,
     &           IOPSM,
     &           LCCBD12,LCCB,     
     &    WORK(KIX1_CA),WORK(KSX1_CA),WORK(KIX1_CB),WORK(KSX1_CB),
     &    WORK(KIX1_AA),WORK(KSX1_AA),WORK(KIX1_AB),WORK(KSX1_AB),
     &    WORK(KIX2_CA),WORK(KSX2_CA),WORK(KIX2_CB),WORK(KSX2_CB),
     &    WORK(KIX2_AA),WORK(KSX2_AA),WORK(KIX2_AB),WORK(KSX2_AB),
     &    WORK(KIX3_CA),WORK(KSX3_CA),WORK(KIX3_CB),WORK(KSX3_CB),
     &    WORK(KIX3_AA),WORK(KSX3_AA),WORK(KIX3_AB),WORK(KSX3_AB),
     &    WORK(KIX4_CA),WORK(KSX4_CA),WORK(KIX4_CB),WORK(KSX4_CB),
     &    WORK(KIX4_AA),WORK(KSX4_AA),WORK(KIX4_AB),WORK(KSX4_AB),
     &    WORK(KIX5_CA),WORK(KSX5_CA),WORK(KIX5_CB),WORK(KSX5_CB),
     &    WORK(KIX5_AA),WORK(KSX5_AA),WORK(KIX5_AB),WORK(KSX5_AB),
     &    WORK(KIX6_CA),WORK(KSX6_CA),WORK(KIX6_CB),WORK(KSX6_CB),
     &    WORK(KIX6_AA),WORK(KSX6_AA),WORK(KIX6_AB),WORK(KSX6_AB),
     &    WORK(KIX7_CA),WORK(KSX7_CA),WORK(KIX7_CB),WORK(KSX7_CB),
     &    WORK(KIX7_AA),WORK(KSX7_AA),WORK(KIX7_AB),WORK(KSX7_AB),
     &    WORK(KIX8_CA),WORK(KSX8_CA),WORK(KIX8_CB),WORK(KSX8_CB),
     &    WORK(KIX8_AA),WORK(KSX8_AA),WORK(KIX8_AB),WORK(KSX8_AB),
     &    WORK(KLOCK1),WORK(KLOCK2),WORK(KLOCK3),
     &    WORK(KLOCK4),WORK(KLOCL1),WORK(KLOCL2),WORK(KLOCL3),
     &    WORK(KLOCL4),
     &    WORK(KLSTOCC1),WORK(KLSTOCC2),WORK(KLSTOCC3),
     &    WORK(KLSTOCC4),WORK(KLSTOCC5),WORK(KLSTOCC6),
     &    WORK(KLSTOCC7),
     &    WORK(KLTSCR1),WORK(KLTSCR2),WORK(KLTSCR3),
Classe     &    WORK(KLTSCR4),WORK(KLTSCR5),WORK(KLTSCR6),
     &    WORK(KLOPSCR),
     &    WORK(KLSMD1),WORK(KLSMD2),WORK(KLSMD3),WORK(KLSMD4),
     &    WORK(KLSMEX),WORK(KLSMK1),WORK(KLSML1),
     &    WORK(KLNMD1),WORK(KLNMD2),WORK(KLNMD3),WORK(KLNMD4),
     &    WORK(KLNMEX),WORK(KLNMK1),WORK(KLNML1),
     &    IOPINDX,
     &    WORK(KLZ),WORK(KLZSCR),WORK(KLSTREO),
     &    ISIGN,N_TDL_MAX)
*
      RETURN
      END
*
      SUBROUTINE OPCT1234_REL(IOEX,IO1DX,!IO2DX,IO3DX,IO4DX,
     &           IT1,!IT2,IT3,IT4,
     &           T1,!T2,T3,T4,
     &           OT1234,
     &           IT1SM,!IT2SM,IT3SM,IT4SM,
     &           IOPSM,
     &           LDB,LB,
     &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,
     &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,
     &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &           IX3_CA,SX3_CA,IX3_CB,SX3_CB,
     &           IX3_AA,SX3_AA,IX3_AB,SX3_AB,
     &           IX4_CA,SX4_CA,IX4_CB,SX4_CB,
     &           IX4_AA,SX4_AA,IX4_AB,SX4_AB,
     &           IX5_CA,SX5_CA,IX5_CB,SX5_CB,
     &           IX5_AA,SX5_AA,IX5_AB,SX5_AB,
     &           IX6_CA,SX6_CA,IX6_CB,SX6_CB,
     &           IX6_AA,SX6_AA,IX6_AB,SX6_AB,
     &           IX7_CA,SX7_CA,IX7_CB,SX7_CB,
     &           IX7_AA,SX7_AA,IX7_AB,SX7_AB,
     &           IX8_CA,SX8_CA,IX8_CB,SX8_CB,
     &           IX8_AA,SX8_AA,IX8_AB,SX8_AB,
     &           IOC_K1,IOC_K2,IOC_K3,IOC_K4,
     &           IOC_L1,IOC_L2,IOC_L3,IOC_L4,
     &           ISTR_D1,ISTR_D2,ISTR_D3,ISTR_D4,ISTR_EX,ISTR_K,ISTR_L,
     &           TSCR1,TSCR2,TSCR3,
Classe     &           TSCR4,TSCR5,TSCR6,
     &           OPSCR,
     &           ISM_CAAB_D1,ISM_CAAB_D2,ISM_CAAB_D3,ISM_CAAB_D4,
     &           ISM_CAAB_EX,ISM_CAAB_K,ISM_CAAB_L,
     &           INM_CAAB_D1,INM_CAAB_D2,INM_CAAB_D3,INM_CAAB_D4,
     &           INM_CAAB_EX,INM_CAAB_K,INM_CAAB_L,
     &           IEXD1234_INDX,
     &           IZ,IZSCR,ISTREO,
     &           ISIGNG,N_TDL_MAX)
*
* Contract indeces of Operator O with indeces of 
* excitation operator T1,T2, T3, T4
*
* Operator O is defined by an operator part IOEX, and 
* deexcitation parts IO1DX,IO2DX,IO3DX,IO4DX.
*
* LDB : Batch size for D1, D2, D3, D4
* LB2   : Batch size for for remaining expansions
*
* Jeppe Olsen, May 2000, Finished March 2003 ( Well it became July 2003)
*                        and even November 2004
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "wrkspc.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "ctcc.inc"
#include "symm.inc"
*. Input   
      INTEGER IOEX(NGAS,4),IO1DX(NGAS,4)
C     INTEGER IO2DX(NGAS,4),IO3DX(NGAS,4),IO4DX(NGAS,4)
      INTEGER IT1(NGAS,4)!,IT2(NGAS,4),IT3(NGAS,4),IT4(NGAS,4)
      DIMENSION T1(*)!,T2(*),T3(*),T4(*)
*. Map of index EXD1234 operator to original index of Hamiltonian
      INTEGER IEX1234_INDX(4)
*. Output
      DIMENSION OT1234(*)
*
*. Local scratch
*
*. Occupation of conjugated operators
      INTEGER IO1DX_DAG(4*MXPNGAS)!,IO2DX_DAG(4*MXPNGAS)
C     INTEGER IO3DX_DAG(4*MXPNGAS),IO4DX_DAG(4*MXPNGAS)
*. Occupation of OP T1 T2
C     INTEGER IOC_OT1T2(NGAS,4)
*
*. Occupation of gasspaces for various strings 
*. Number of strings per sym of the various string supergroups
      INTEGER NOEX(8,4),NO1DX(8,4),NO2DX(8,4),NO3DX(8,4),NO4DX(8,4)
      INTEGER NT1(8,4), NT2(8,4), NT3(8,4), NT4(8,4)
      INTEGER NK1(8,4),NK2(8,4),NK3(8,4),NK4(8,4)
      INTEGER NL1(8,4),NL2(8,4),NL3(8,4),NL4(8,4)
C     INTEGER NOT1T2(8,4)
*. Offsets of strings with given sym for the strings of various CC ops
      INTEGER IBOEX(8,4),IBO1DX(8,4)!,IBO2DX(8,4),IBO3DX(8,4)
C     INTEGER IBO4DX(8,4)
*. Offset in operators to strings with given sym
      INTEGER IBT1_TCC(8,8,8)!,IBT2_TCC(8,8,8)
C     INTEGER IBT3_TCC(8,8,8),IBT4_TCC(8,8,8)
      INTEGER IBL1_TCC(8,8,8)!,IBL2_TCC(8,8,8)
C     INTEGER IBL3_TCC(8,8,8),IBL4_TCC(8,8,8)
C     INTEGER IBOT1T2_TCC(8,8,8)
*
C     INTEGER IBOT1T2_TCC(8,8,8)
*
      INTEGER IB_D1K1(8,8,4)!,IB_D2K2(8,8,4) 
C     INTEGER IB_D3K3(8,8,4),IB_D4K4(8,8,4)
      INTEGER IB_EXK1(8,8,4)!,IB_L1K2(8,8,4)
C     INTEGER IB_L2K3(8,8,4),IB_L3K4(8,8,4)
*
*
*. Scratch through parameter list. 
*
*. IX1_* : Number of operators in operator * Largest C or A string
      DIMENSION IX1_CA(*),SX1_CA(*),IX1_CB(*),SX1_CB(*)
      DIMENSION IX1_AA(*),SX1_AA(*),IX1_AB(*),SX1_AB(*)
*. IX2, SX2
      DIMENSION IX2_CA(*),SX2_CA(*),IX2_CB(*),SX2_CB(*)
      DIMENSION IX2_AA(*),SX2_AA(*),IX2_AB(*),SX2_AB(*)
*. IX3, SX3
      DIMENSION IX3_CA(*),SX3_CA(*),IX3_CB(*),SX3_CB(*)
      DIMENSION IX3_AA(*),SX3_AA(*),IX3_AB(*),SX3_AB(*)
*. IX4, SX4
      DIMENSION IX4_CA(*),SX4_CA(*),IX4_CB(*),SX4_CB(*)
      DIMENSION IX4_AA(*),SX4_AA(*),IX4_AB(*),SX4_AB(*)
*. IX5, SX5
      DIMENSION IX5_CA(*),SX5_CA(*),IX5_CB(*),SX5_CB(*)
      DIMENSION IX5_AA(*),SX5_AA(*),IX5_AB(*),SX5_AB(*)
*. IX6, SX6
      DIMENSION IX6_CA(*),SX6_CA(*),IX6_CB(*),SX6_CB(*)
      DIMENSION IX6_AA(*),SX6_AA(*),IX6_AB(*),SX6_AB(*)
*. IX7, SX7
      DIMENSION IX7_CA(*),SX7_CA(*),IX7_CB(*),SX7_CB(*)
      DIMENSION IX7_AA(*),SX7_AA(*),IX7_AB(*),SX7_AB(*)
*. IX8, SX8
      DIMENSION IX8_CA(*),SX8_CA(*),IX8_CB(*),SX8_CB(*)
      DIMENSION IX8_AA(*),SX8_AA(*),IX8_AB(*),SX8_AB(*)
*
*. ISTR_D1, ISTR_K, ISTR_L : for occupations of strings  
      INTEGER ISTR_D1(*)!, ISTR_D2(*), ISTR_D3(*), ISTR_D4(*)
      INTEGER ISTR_EX(*)
      INTEGER ISTR_K(*), ISTR_L(*) 
*. For occupation of intermediate strings
      INTEGER IOC_K1(NGAS,4)!,IOC_K2(NGAS,4)
C     INTEGER IOC_K3(NGAS,4),IOC_K4(NGAS,4)
      INTEGER IOC_L1(NGAS,4)!,IOC_L2(NGAS,4)
C     INTEGER IOC_L3(NGAS,4),IOC_L4(NGAS,4)
*. For intermediates with both strings batched
*. Changed from squared to cubed by Lasse
C      DIMENSION TSCR1(LDB*LDB*LDB),TSCR2(LDB*LDB*LDB)
C      DIMENSION TSCR3(LDB*LDB*LDB)
      DIMENSION TSCR1(*),TSCR2(*)
      DIMENSION TSCR3(*)
*. For an intermediate with only one string batched. Now allocated
*. inside since problems with dimensioning. Lasse
C      DIMENSION TSCR4(N_TDL_MAX), TSCR5(N_TDL_MAX),TSCR6(N_TDL_MAX)
*. For a batch of coefficients for Operator
*. Operator accessed as OP(EX,D2,D1) so
C      DIMENSION OPSCR(LDB*LDB*LB)
      DIMENSION OPSCR(*)
*. For part of Hamiltonian
      INTEGER ISM_CAAB_D1(4,*), INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_D2(4,*), INM_CAAB_D2(4,*)
      INTEGER ISM_CAAB_D3(4,*), INM_CAAB_D3(4,*)
      INTEGER ISM_CAAB_D4(4,*), INM_CAAB_D4(4,*)
      INTEGER ISM_CAAB_EX(4,*), INM_CAAB_EX(4,*)
*. For CC operators
      INTEGER ISM_CAAB_K(4,*), INM_CAAB_K(4,*)
      INTEGER ISM_CAAB_L(4,*), INM_CAAB_L(4,*)
*
      INTEGER IZ(*), IZSCR(*), ISTREO(*)
      INTEGER IGRP_TEST(MXPNGAS)
*. Collecting info on largest T(D,L) Block
      COMMON/CNTDL_MAX/NTDL_MAX_ACT,ID2_MAX_ACT(MXPNGAS*4),
     &      IL1_MAX_ACT(MXPNGAS*4),IEX_MAX_ACT(MXPNGAS*4),
     &      IT1_MAX_ACT(MXPNGAS*4),IT2_MAX_ACT(MXPNGAS*4)
*
* The story goes as
* Loop over symmetry of D4
* Loop over batches  of D4
*  Loop over symmetry of D3
*  Loop over batches  of D3
*   Loop over symmetry of D2
*   Loop over batches  of D2
*
* Part 1 : O, T1 => OT1(d2,d3,d4,l1) in TSCR4
*
*    Loop over symmetry of D1
*    Loop over batches  of D1
*
*    Loop over batches of Ex
*      Fetch O as O(d1,d2,d3,d4,ex)
*      Loop over batches of K1
*        T1(I) => T1(d1,k1)
*        OT1(d2,d3,d4,ex,k1) = sum(d1) O(d1,d2,d3,d4,ex)T1(d1,k1)
*        OT1(d2,d3,d4,ex,k1) => OT1(d2,d3,d4,l1)
*      End of loop over batches of K1
*    End of loop over bathces of Ex
*    End of loop over batches of d1
*    End of loop over symmetry of D1
* ( We now have OT1(l1,d4,d3,d2), all l1 !, d2,d3,d4 in batch)    
*
* Part 2 : OT1 * T2 => OT12(d3,d4,l2) in TSCR5
*
*    Loop over batches of k2
*      T2(i) => T2(d2,k2)
*      Loop over batches of l1
*        OT12(d3,d4,l1,k2) = OT1(d2,d3,d4,l1)*T2(d2,k2) 
*        OT12(d3,d4,l1,k2) => OT12(d3,d4,l2)
*      End of loop over batches of l1
*    End of loop over batches of K2
*   End of loop over batches  of D2
*   End of loop over symmetry of D2
*
* Part 3 : OT12 * T3 => OT123(d4,l3) in TSCR6
*
*   Loop over batches of K3
*     T3(I) => T2(D3,K3)
*     Loop over batches of L2
*       OT123(d4,l2,k3) = sum(d3) OT12(d3,d4,l2) T2(d3,k3) 
*       OT123(d4,l2,k3) => OT123(d4,l3)
*     End of loop over batches of L2
*   End of loop over batches of K3
*  End of loop over batches  of D3
*  End of loop over symmetry of D3
*
* Part 4 : OT123 * T4 => OT1234(I)
*
*  Loop over batches of K4
*    T4(I) => T4(d4,k4)
*    Loop over batches of L3
*     OT1234(l3,k4) = OT123(d4,l3)T4(d4,k4)
*     OT1234(l3,k4) => OT1234(I)
*    End of loop over batches of L3
*  End of loop over batches of K4
* End of loop over batches of D4
* End of loop over symmetry of D4
*
* - That's all she wrote
*
C     CALL MEMMAR(IDUM,IDUM,'MARK  ',IDUM,'OPCTT ')
*
      ZERO = 0.0D0
      ONE = 1.0D0
*
      NTEST = 100 
      IF(NTEST.GE.100 )THEN
        WRITE(6,*) ' OPCT1234_REL entered '
        WRITE(6,*) ' ==================== '
        WRITE(6,*) ' Ex Deex1 Deex2 Deex3 Deex4 ' 
        WRITE(6,*)
        CALL WRT_SPOX_TP_CC_KRCC(IOEX,1)
        CALL WRT_SPOX_TP_CC_KRCC(IO1DX,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IO2DX,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IO3DX,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IO4DX,1)
        WRITE(6,*) 
        WRITE(6,*) ' Form of T1 T2 T3 T4'
        CALL WRT_SPOX_TP_CC_KRCC(IT1,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IT2,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IT3,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IT4,1)
*
      END IF
      IF(NTEST.GE.100) 
     &WRITE(6,*) ' N_TDL_MAX in OPCT1T2 ', N_TDL_MAX
      SIGNG = DFLOAT(ISIGNG)
*. Symmetry of final type  
C     IF(NTEST.GE.100) WRITE(6,*) ' IT1SM, IT2SM, IT3SM, IT4SM : ',
C    &             IT1SM, IT2SM, IT3SM, IT4SM
      IF(NTEST.GE.100) WRITE(6,*) ' IT1SM : ',IT1SM
C     IT12SM = IDBGMULT(IT1SM,IT2SM)
C     IT123SM= IDBGMULT(IT12SM,IT3SM)
C     IT1234SM= IDBGMULT(IT123SM,IT4SM)
C     IL4SM= IDBGMULT(IOPSM,IT1234SM)
      IL4SM= IDBGMULT(IOPSM,IT1SM)
C     IF(NTEST.GE.100) WRITE(6,*) 'IT12SM,IT123SM,IT1234SM,IL4SM ',
C    &            IT12SM,IT123SM,IT1234SM,IL4SM
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Symmetry of OT1234 = ', IL4SM
      END IF
*
* =================================================
* Occupation of the various gas spaces for strings 
* =================================================
*
*. K1 = O1DX  I1 
      CALL NEW_CAAB_OC(IOC_K1,IT1,IO1DX,1,1,NGAS)
*. L1 = OEX K1 
      CALL NEW_CAAB_OC(IOC_L1,IOC_K1,IOEX,2,1,NGAS)
*. K2 = O2DX I2
C     CALL NEW_CAAB_OC(IOC_K2,IT2,IO2DX,1,1,NGAS)
*. L2 = K2 L1
C     CALL NEW_CAAB_OC(IOC_L2,IOC_K2,IOC_L1,2,1,NGAS)
*. K3 = O3DX I3
C     CALL NEW_CAAB_OC(IOC_K3,IT3,IO3DX,1,1,NGAS)
*. L3 = K3 L2
C     CALL NEW_CAAB_OC(IOC_L3,IOC_K3,IOC_L2,2,1,NGAS)
*. K4 = O2DX I2
C     CALL NEW_CAAB_OC(IOC_K4,IT4,IO4DX,1,1,NGAS)
*. L4 = K4 L3
C     CALL NEW_CAAB_OC(IOC_L4,IOC_K4,IOC_L3,2,1,NGAS)
*. (L4 is final type of OT1234)
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Occupation of K1, K2, K3, K4 '
        CALL WRT_SPOX_TP_CC_KRCC(IOC_K1,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IOC_K2,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IOC_K3,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IOC_K4,1)
        WRITE(6,*) ' Occupation of L1, L2, L3, L4 '
        CALL WRT_SPOX_TP_CC_KRCC(IOC_L1,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IOC_L2,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IOC_L3,1)
C       CALL WRT_SPOX_TP_CC_KRCC(IOC_L4,1)
      END IF
*
* =============================================================
*. Obtain symmetry-dimensions and -offsets for various strings
* =============================================================
*
      DO I_CAAB = 1, 4
       call evenodd(IEO,I_CAAB)
       IF(IEO.EQ.1) THEN
         IUB = 1
       ELSE
         IUB = 2
       END IF
       IF(I_CAAB.GE.3) THEN
         ICA = 2
       ELSE
         ICA = 1
       END IF
C NST_SPGP_CC_REL should be checked!
C check with nst_spgrp_dbg from sigden.
       CALL NST_SPGP_CC_REL(IOEX(1,I_CAAB),NOEX(1,I_CAAB),IUB,ICA)
       CALL ZBASE(NOEX(1,I_CAAB),IBOEX(1,I_CAAB),NSMST)
*. D1
       CALL NST_SPGP_CC_REL(IO1DX(1,I_CAAB),NO1DX(1,I_CAAB),IUB,ICA)
       CALL ZBASE(NO1DX(1,I_CAAB),IBO1DX(1,I_CAAB),NSMST)
*. D2
C      CALL NST_SPGP_CC_REL(IO2DX(1,I_CAAB),NO2DX(1,I_CAAB),IUB,ICA)
C      CALL ZBASE(NO2DX(1,I_CAAB),IBO2DX(1,I_CAAB),NSMST)
*. D3
C      CALL NST_SPGP_CC_REL(IO3DX(1,I_CAAB),NO3DX(1,I_CAAB),IUB,ICA)
C      CALL ZBASE(NO3DX(1,I_CAAB),IBO3DX(1,I_CAAB),NSMST)
*. D4
C      CALL NST_SPGP_CC_REL(IO4DX(1,I_CAAB),NO4DX(1,I_CAAB),IUB,ICA)
C      CALL ZBASE(NO4DX(1,I_CAAB),IBO4DX(1,I_CAAB),NSMST)
*. T1
       CALL NST_SPGP_CC_REL(IT1(1,I_CAAB),NT1(1,I_CAAB),IUB,ICA)
*. T2
C      CALL NST_SPGP_CC_REL(IT2(1,I_CAAB),NT2(1,I_CAAB),IUB,ICA)
*. T3
C      CALL NST_SPGP_CC_REL(IT3(1,I_CAAB),NT3(1,I_CAAB),IUB,ICA)
*. T4
C      CALL NST_SPGP_CC_REL(IT4(1,I_CAAB),NT4(1,I_CAAB),IUB,ICA)
*. K1
       CALL NST_SPGP_CC_REL(IOC_K1(1,I_CAAB),NK1(1,I_CAAB),IUB,ICA)
*. K2
C      CALL NST_SPGP_CC_REL(IOC_K2(1,I_CAAB),NK2(1,I_CAAB),IUB,ICA)
*. K3
C      CALL NST_SPGP_CC_REL(IOC_K3(1,I_CAAB),NK3(1,I_CAAB),IUB,ICA)
*. K4
C      CALL NST_SPGP_CC_REL(IOC_K4(1,I_CAAB),NK4(1,I_CAAB),IUB,ICA)
*. L1   
       CALL NST_SPGP_CC_REL(IOC_L1(1,I_CAAB),NL1(1,I_CAAB),IUB,ICA)
*. L2   
C      CALL NST_SPGP_CC_REL(IOC_L2(1,I_CAAB),NL2(1,I_CAAB),IUB,ICA)
*. L3   
C      CALL NST_SPGP_CC_REL(IOC_L3(1,I_CAAB),NL3(1,I_CAAB),IUB,ICA)
*. L4   
C      CALL NST_SPGP_CC_REL(IOC_L4(1,I_CAAB),NL4(1,I_CAAB),IUB,ICA)
*
      END DO
*. We now have the various dimensions, so we can write T1,T2,T3,T4 if 
*. required
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Input T1 block '
        CALL WRT_TCC_BLK_CC_REL(T1,IT1SM,NT1(1,1),NT1(1,2),NT1(1,3),
     &                   NT1(1,4),NSMST)
C       WRITE(6,*) ' Input T2 block '
C       CALL WRT_TCC_BLK_CC_REL(T2,IT2SM,NT2(1,1),NT2(1,2),NT2(1,3),
C    &                   NT2(1,4),NSMST)
C       WRITE(6,*) ' Input T3 block '
C       CALL WRT_TCC_BLK_CC_REL(T3,IT3SM,NT3(1,1),NT3(1,2),NT3(1,3),
C    &                   NT3(1,4),NSMST)
C       WRITE(6,*) ' Input T4 block '
C       CALL WRT_TCC_BLK_CC_REL(T4,IT4SM,NT4(1,1),NT4(1,2),NT4(1,3),
C    &                   NT4(1,4),NSMST)
C       WRITE(6,*) ' Input OT1234 block '
C       CALL WRT_TCC_BLK_CC_REL(OT1234,IL4SM,NL4(1,1),NL4(1,2),NL4(1,3),
C    &        NL4(1,4),NSMST)
        WRITE(6,*) ' Input OT1234 block '
        CALL WRT_TCC_BLK_CC_REL(OT1234,IL4SM,NL1(1,1),NL1(1,2),NL1(1,3),
     &        NL1(1,4),NSMST)
      END IF
      stop 'for now again'
*
*
* ============================
*. Offsets to T1, T2, T3, T4
* ============================
*
*. T1
      CALL Z_TCC_OFF2_REL(IBT1_TCC,LEN_T1,NT1(1,1),NT1(1,2),NT1(1,3),
     &                NT1(1,4),IT1SM,NSMST)
*. T2
C     CALL Z_TCC_OFF2_REL(IBT2_TCC,LEN_T2,NT2(1,1),NT2(1,2),NT2(1,3),
C    &                NT2(1,4),IT2SM,NSMST)
*. T3
C     CALL Z_TCC_OFF2_REL(IBT3_TCC,LEN_T3,NT3(1,1),NT3(1,2),NT3(1,3),
C    &                NT3(1,4),IT3SM,NSMST)
*. T4
C     CALL Z_TCC_OFF2_REL(IBT4_TCC,LEN_T4,NT4(1,1),NT4(1,2),NT4(1,3),
C    &                NT4(1,4),IT4SM,NSMST)
*. The following has been moved as symmetry of L1 - L4 is not known here
* This is because of the loops over symmetry later. And L1-L4 is a product of these. Lasse
*. L1
C     CALL Z_TCC_OFF2(IBL1_TCC,LEN_L1,NL1(1,1),
C    &     NL1(1,2),NL1(1,3),NL1(1,4),IL1,NSMST)
*. L2
C     CALL Z_TCC_OFF2(IBL2_TCC,LEN_L2,NL2(1,1),
C    &     NL2(1,2),NL2(1,3),NL2(1,4),IL2,NSMST)
*. L4
C     CALL Z_TCC_OFF2(IBL3_TCC,LEN_L3,NL3(1,1),
C    &     NL3(1,2),NL3(1,3),NL3(1,4),IL3,NSMST)
*. L4
C     CALL Z_TCC_OFF2(IBL4_TCC,LEN_L4,NL4(1,1),
C    &     NL4(1,2),NL4(1,3),NL4(1,4),IL4,NSMST)
*
*
* ========================
*.  D1{\dagger} K1  => T1 mapping
* =========================
*
*. Obtain D1{\dagger}
      CALL CONJ_CAAB(IO1DX,IO1DX_DAG,NGAS,SP_D1K1)
      CALL T1T2_TO_T12_MAP_REL(IO1DX_DAG,IOC_K1,IT1,
     &     IX1_CA,SX1_CA,IX1_CB,SX1_CB,IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &     IB_D1K1,ISTR_D1,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D1K1)
*
* ======================
* Ex K1 => L1 mapping
* ======================
*
      CALL T1T2_TO_T12_MAP_REL(IOEX,IOC_K1,IOC_L1,
     &     IX2_CA,SX2_CA,IX2_CB,SX2_CB,IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &     IB_EXK1,ISTR_K,ISTR_L,ISTREO,IZ,IZSCR,SIGN_EXK1)
*
* ========================
*.  D2{\dagger} K2  => T2 mapping
* =========================
*
*. Obtain D2{\dagger}
C     CALL CONJ_CAAB(IO2DX,IO2DX_DAG,NGAS,SP_D2K2)
C     CALL T1T2_TO_T12_MAP_REL(IO2DX_DAG,IOC_K2,IT2,
C    &     IX3_CA,SX3_CA,IX3_CB,SX3_CB,IX3_AA,SX3_AA,IX3_AB,SX3_AB,
C    &     IB_D2K2,ISTR_D2,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D2K2)
*
* ==========
* L1 K2 => L2
* ==========
*
C     CALL T1T2_TO_T12_MAP_REL(IOC_L1,IOC_K2,IOC_L2,
C    &     IX4_CA,SX4_CA,IX4_CB,SX4_CB,IX4_AA,SX4_AA,IX4_AB,SX4_AB,
C    &     IB_L1K2,ISTR_L,ISTR_K,ISTREO,IZ,IZSCR,SIGN_L1K2)
*
* ========================
*.  D3{\dagger} K3  => T3 mapping
* =========================
*
*. Obtain D3{\dagger}
C     CALL CONJ_CAAB(IO3DX,IO3DX_DAG,NGAS,SP_D3K3)
C     CALL T1T2_TO_T12_MAP_REL(IO3DX_DAG,IOC_K3,IT3,
C    &     IX5_CA,SX5_CA,IX5_CB,SX5_CB,IX5_AA,SX5_AA,IX5_AB,SX5_AB,
C    &     IB_D3K3,ISTR_D3,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D3K3)
*
* ==========
* L2 K3 => L3
* ==========
*
C     CALL T1T2_TO_T12_MAP_REL(IOC_L2,IOC_K3,IOC_L3,
C    &     IX6_CA,SX6_CA,IX6_CB,SX6_CB,IX6_AA,SX6_AA,IX6_AB,SX6_AB,
C    &     IB_L2K3,ISTR_L,ISTR_K,ISTREO,IZ,IZSCR,SIGN_L2K3)
* ========================
*.  D4{\dagger} K4  => T4 mapping
* =========================
*
*. Obtain D4{\dagger}
C     CALL CONJ_CAAB(IO4DX,IO4DX_DAG,NGAS,SP_D4K4)
C     CALL T1T2_TO_T12_MAP_REL(IO4DX_DAG,IOC_K4,IT4,
C    &     IX7_CA,SX7_CA,IX7_CB,SX7_CB,IX7_AA,SX7_AA,IX7_AB,SX7_AB,
C    &     IB_D4K4,ISTR_D4,ISTR_K,ISTREO,IZ,IZSCR,SIGN_D4K4)
*
* ==========
* L3 K4 => L4
* ==========
*
C     CALL T1T2_TO_T12_MAP_REL(IOC_L3,IOC_K4,IOC_L4,
C    &     IX8_CA,SX8_CA,IX8_CB,SX8_CB,IX8_AA,SX8_AA,IX8_AB,SX8_AB,
C    &     IB_L3K4,ISTR_L,ISTR_K,ISTREO,IZ,IZSCR,SIGN_L3K4)
*

*. And the individual strings : all symmetries constructed
*. D1 strings 
C Check STR_CAAB_REL for bugs
      CALL STR_CAAB_REL(IO1DX,ISTR_D1)
*. D2 strings 
C     CALL STR_CAAB_REL(IO2DX,ISTR_D2)
*. D3 strings 
C     CALL STR_CAAB_REL(IO3DX,ISTR_D3)
*. D4 strings 
C     CALL STR_CAAB_REL(IO4DX,ISTR_D4)
*. Ex strings 
      CALL STR_CAAB_REL(IOEX,ISTR_EX)
*. Batch length : Should allow to vary 
      LD1_BAT = LDB
C     LD2_BAT = LDB
C     LD3_BAT = LDB
C     LD4_BAT = LDB
      LEX_BAT = LDB
      LKL_BAT = LDB 
*. Calculate length of TSCR4, TSCR5 and TSCR6. Memory saving by Lasse
*. Do not allocate more than what is needed
*. Assumes at the moment static batching
      NL1_TOT_MAX = 0
C     NL2_TOT_MAX = 0
C     NL3_TOT_MAX = 0
C     ND2_TOT_MAX = 0
C     ND3_TOT_MAX = 0
C     ND4_TOT_MAX = 0
      DO IISYM = 1,NSMST
        NL1_TOT_INT = LEN_TCCBLK_REL(NL1(1,1),NL1(1,2),NL1(1,3),
     &                               NL1(1,4),IISYM,NSMST)
C       NL2_TOT_INT = LEN_TCCBLK_REL(NL2(1,1),NL2(1,2),NL2(1,3),
C    &                               NL2(1,4),IISYM,NSMST)
C       NL3_TOT_INT = LEN_TCCBLK_REL(NL3(1,1),NL3(1,2),NL3(1,3),
C    &                               NL3(1,4),IISYM,NSMST)
        NL1_TOT_MAX = MAX(NL1_TOT_MAX,NL1_TOT_INT)
C       NL2_TOT_MAX = MAX(NL2_TOT_MAX,NL2_TOT_INT)
C       NL3_TOT_MAX = MAX(NL3_TOT_MAX,NL3_TOT_INT)
C        print*,'NL1_TOT_INT,NL2_TOT_INT,NL3_TOT_INT',
C     &          NL1_TOT_INT,NL2_TOT_INT,NL3_TOT_INT
C       ND2_TOT_INT = LEN_TCCBLK_REL(NO2DX(1,1),NO2DX(1,2),NO2DX(1,3),
C    &                               NO2DX(1,4),IISYM,NSMST)
C       ND3_TOT_INT = LEN_TCCBLK_REL(NO3DX(1,1),NO3DX(1,2),NO3DX(1,3),
C    &                               NO3DX(1,4),IISYM,NSMST)
C       ND4_TOT_INT = LEN_TCCBLK_REL(NO4DX(1,1),NO4DX(1,2),NO4DX(1,3),
C    &                               NO4DX(1,4),IISYM,NSMST)
C       ND2_TOT_MAX = MAX(ND2_TOT_MAX,ND2_TOT_INT)
C       ND3_TOT_MAX = MAX(ND3_TOT_MAX,ND3_TOT_INT)
C       ND4_TOT_MAX = MAX(ND4_TOT_MAX,ND4_TOT_INT)
C        print*,'ND2_TOT_INT,ND3_TOT_INT,ND4_TOT_INT',
C     &          ND2_TOT_INT,ND3_TOT_INT,ND4_TOT_INT
      END DO
CC     print*,'NL1_TOT_MAX,NL2_TOT_MAX,NL3_TOT_MAX',
CC    &        NL1_TOT_MAX,NL2_TOT_MAX,NL3_TOT_MAX
CC     print*,'ND2_TOT_MAX,ND3_TOT_MAX,ND4_TOT_MAX',
CC    &        ND2_TOT_MAX,ND3_TOT_MAX,ND4_TOT_MAX
CC     ND2_BAT_MAX = ND2_TOT_MAX/LD2_BAT
CC     IF(ND2_BAT_MAX*LD2_BAT.LT.ND2_TOT_MAX) ND2_BAT_MAX=ND2_BAT_MAX+1
CC     ND3_BAT_MAX = ND3_TOT_MAX/LD3_BAT
CC     IF(ND3_BAT_MAX*LD3_BAT.LT.ND3_TOT_MAX) ND3_BAT_MAX=ND3_BAT_MAX+1
CC     ND4_BAT_MAX = ND4_TOT_MAX/LD4_BAT
CC     IF(ND4_BAT_MAX*LD4_BAT.LT.ND4_TOT_MAX) ND4_BAT_MAX=ND4_BAT_MAX+1
CC     ID2_START_MIN = 1
CC     ID3_START_MIN = 1
CC     ID4_START_MIN = 1
CC     ID2_START_MAX = (ND2_BAT_MAX - 1) * LD2_BAT + 1
CC     ID3_START_MAX = (ND3_BAT_MAX - 1) * LD3_BAT + 1
CC     ID4_START_MAX = (ND4_BAT_MAX - 1) * LD4_BAT + 1
C     IF(ND2_TOT_MAX.GE.LD2_BAT) THEN
C       ID2_BATLEN_MAX = LD2_BAT
C     ELSE
C       ID2_BATLEN_MAX = MAX(1,ND2_TOT_MAX)
C     END IF
C     IF(ND3_TOT_MAX.GE.LD3_BAT) THEN
C       ID3_BATLEN_MAX = LD3_BAT
C     ELSE
C       ID3_BATLEN_MAX = MAX(1,ND3_TOT_MAX)
C     END IF
C     IF(ND4_TOT_MAX.GE.LD4_BAT) THEN
C       ID4_BATLEN_MAX = LD4_BAT
C     ELSE
C       ID4_BATLEN_MAX = MAX(1,ND4_TOT_MAX)
C     END IF
C     LEN4 = NL1_TOT_MAX*ID2_BATLEN_MAX*ID3_BATLEN_MAX*ID4_BATLEN_MAX
C     LEN5 = NL2_TOT_MAX*ID3_BATLEN_MAX*ID4_BATLEN_MAX
C     LEN6 = NL3_TOT_MAX*ID4_BATLEN_MAX
C      print*,'LEN4,LEN5,LEN6,LDB',LEN4,LEN5,LEN6,LDB
C     CALL MEMMAR(KLTSCR4,LEN4,'ADDL  ',2,'TSCR4 ')
C     CALL MEMMAR(KLTSCR5,LEN5,'ADDL  ',2,'TSCR5 ')
C     CALL MEMMAR(KLTSCR6,LEN6,'ADDL  ',2,'TSCR6 ')
*. Symmetry of Op T1, Op T1 T2, Op T1 T2 T3, Op T1 T2 T3 T4
*
* Check Op for annihilation symmetry
*
* IOPT1SM now final symmetry
      IOPT1SM    = IDBGMULT(IOPSM,IT1SM)
C     IOPT12SM   = IDBGMULT(IOPT1SM,IT2SM)
C     IOPT123SM  = IDBGMULT(IOPT12SM,IT3SM)
C     IOPT1234SM = IDBGMULT(IOPT123SM,IT4SM)
*
* Needs to be changed Lasse 2011
*
      NL4_TOT = LEN_TCCBLK_REL(NL4(1,1),NL4(1,2),NL4(1,3),
     &                 NL4(1,4),IOPT1234SM,NSMST)
*
C     DO ID4SM = 1, NSMST
C correct symmetry ? yes
C      K4SM = IDBGMULT(ID4SM,IT4SM)
*. By def L4 = final operator so 
C      L4SM = IOPT1234SM
C      NL4_TOT = LEN_TCCBLK_REL(NL4(1,1),NL4(1,2),NL4(1,3),
C    &                NL4(1,4),L4SM,NSMST)
C      ID4SM_NEW = 1
*. Offsets for L4 with sym L4SM
C      CALL Z_TCC_OFF2_REL(IBL4_TCC,LEN_L4,NL4(1,1),
C    &      NL4(1,2),NL4(1,3),NL4(1,4),L4SM,NSMST)
*. Number of D4 strings with actual symmetry
C      ND4_TOT = LEN_TCCBLK_REL(NO4DX(1,1),NO4DX(1,2),NO4DX(1,3),
C    &               NO4DX(1,4),ID4SM,NSMST)
*. Number of K4 strings with actual symmetry
C      NK4_TOT = LEN_TCCBLK_REL(NK4(1,1),NK4(1,2),NK4(1,3),
C    &                 NK4(1,4),K4SM,NSMST)
*. Number of D4 batches 
C      ND4_BAT = ND4_TOT/LD4_BAT
C      IF(ND4_BAT*LD4_BAT.LT.ND4_TOT) ND4_BAT=ND4_BAT+1
*. And loop over batches of D4
C      DO ID4_BAT = 1, ND4_BAT
C       ID4_START = (ID4_BAT-1)*LD4_BAT + 1
C       ID4_STOP  = MIN(ND4_TOT,ID4_START+LD4_BAT-1)
C       ID4_BATLEN = ID4_STOP-ID4_START+1
C       IF(NTEST.GE.100) THEN
C         WRITE(6,*) ' ID4SM, ID4_BAT, ID4_BATLEN = ',
C    &                 ID4SM, ID4_BAT, ID4_BATLEN
C       END IF
*. Generate D4 strings for given sym and batch
C       CALL ISMNM_FOR_TCC_BAT_REL(NO4DX,ISM_CAAB_D4,INM_CAAB_D4,ID4SM,
C    &       ID4_BATLEN,ID4SM_NEW,
C    &       ISM_C1_D4,ISM_CA1_D4,ISM_AA1_D4,
C    &       INM_AB1_D4,INM_AA1_D4,INM_CA1_D4,INM_CB1_D4,
C    &       ISM_CINI_D4,ISM_CAINI_D4,ISM_AAINI_D4,      
C    &       INM_ABINI_D4,INM_AAINI_D4,INM_CAINI_D4,INM_CBINI_D4,0)
C       ID4SM_NEW = 0
* The results of this loop is the generation of OT123(d4,l3) in TSCR6
C
C Not sure about symmetry
C
C       L3SM = IDBGMULT(IOPT123SM,IADJSYM(ID4SM))
*. Number of L3 strings with actual symmetry
C       NL3_TOT = LEN_TCCBLK_REL(NL3(1,1),NL3(1,2),NL3(1,3),
C    &                 NL3(1,4),L3SM,NSMST)
*. Offsets for L3 with sym L3SM
C       CALL Z_TCC_OFF2_REL(IBL3_TCC,LEN_L3,NL3(1,1),
C    &       NL3(1,2),NL3(1,3),NL3(1,4),L3SM,NSMST)
C       LENGTH = ID4_BATLEN*NL3_TOT
C       IF(LENGTH.GT.LEN6) THEN
C          WRITE(6,*) ' N_TDL_MAX too small '
C          WRITE(6,*) ' Current, required = ', LEN6, LENGTH
C          WRITE(6,*) 'ID4_BATLEN,NL3_TOT = ', ID4_BATLEN, NL3_TOT
C          STOP ' N_TDL_MAX too small '
C       END IF
C       CALL SETVEC(WORK(KLTSCR6),ZERO,LENGTH)
*. Loop over D3
C       DO ID3SM = 1, NSMST
C sym ok
C        K3SM = IDBGMULT(ID3SM,IT3SM)
C        IK3SM_NEW = 1
C        ID3SM_NEW = 1
C        IF(NTEST.GE.100)  WRITE(6,*) ' ID3SM,K3SM = ',ID3SM, K3SM
*. Number of D3 strings with actual symmetry
         
C        ND3_TOT = LEN_TCCBLK_REL(NO3DX(1,1),NO3DX(1,2),NO3DX(1,3),
C    &                 NO3DX(1,4),ID3SM,NSMST)
*. Number of K3 strings with actual symmetry
C        NK3_TOT = LEN_TCCBLK_REL(NK3(1,1),NK3(1,2),NK3(1,3),
C    &                 NK3(1,4),K3SM,NSMST)
*. Number of D3 batches 
C        IF(NK3_TOT.GT.0 .AND. NL3_TOT.GT.0) THEN
C        ND3_BAT = ND3_TOT/LD3_BAT
C        IF(ND3_BAT*LD3_BAT.LT.ND3_TOT) ND3_BAT=ND3_BAT+1
*. And loop over batches of D3
C      IF(NTEST.GE.100) WRITE(6,*) 'ND3_TOT, ND3_BAT ',ND3_TOT,ND3_BAT
C        DO ID3_BAT = 1, ND3_BAT
C         IF(NTEST.GE.100) WRITE(6,*) ' ID3_BAT = ', ID3_BAT
C         ID3_START = (ID3_BAT-1)*LD3_BAT + 1
C         ID3_STOP  = MIN(ND3_TOT,ID3_START+LD3_BAT-1)
C         ID3_BATLEN = ID3_STOP-ID3_START+1
C         IF(NTEST.GE.100) THEN
C           WRITE(6,*) ' ID3SM, ID3_BAT, ID3_BATLEN = ',
C    &                   ID3SM, ID3_BAT, ID3_BATLEN
C         END IF
*. Generate D3 strings for given sym and batch
C         CALL ISMNM_FOR_TCC_BAT_REL(NO3DX,ISM_CAAB_D3,INM_CAAB_D3,
C    &         ID3SM,
C    &         ID3_BATLEN,ID3SM_NEW,
C    &         ISM_C1_D3,ISM_CA1_D3,ISM_AA1_D3,
C    &         INM_AB1_D3,INM_AA1_D3,INM_CA1_D3,INM_CB1_D3,
C    &         ISM_CINI_D3,ISM_CAINI_D3,ISM_AAINI_D3,      
C    &         INM_ABINI_D3,INM_AAINI_D3,INM_CAINI_D3,INM_CBINI_D3,0)
C         ID3SM_NEW = 0
*. In the loops inside, we will construct OT1T2(d3,d4,l2) in TSCR5 
*  by looping over symmetries and batches of D2 and D1. 
C
C Again not sure about Lsym
C
C         ID34SM = IDBGMULT(ID3SM,ID4SM)
C         L2SM   = IDBGMULT(IOPT12SM,IADJSYM(ID34SM))
*. Number of L2 strings with sym L2SM
C         NL2_TOT = LEN_TCCBLK_REL(NL2(1,1),NL2(1,2),NL2(1,3),
C    &                  NL2(1,4),L2SM,NSMST)
*. Offsets for L2 with sym L2SM
C         CALL Z_TCC_OFF2_REL(IBL2_TCC,LEN_L2,NL2(1,1),
C    &         NL2(1,2),NL2(1,3),NL2(1,4),L2SM,NSMST)
C         LENNY = ID3_BATLEN*ID4_BATLEN*NL2_TOT
C         IF(NTEST.GE.100) WRITE(6,*) ' L2SM, NL2_TOT ',
C    &                                  L2SM, NL2_TOT
C         IF(LENNY .GT. LEN5 ) THEN
*. Problem, too small length for TSCR5 
C            WRITE(6,*) ' Dimension of TSCR5 too small '
C            WRITE(6,*) ' Required and allocated : ',LENNY, LEN5 
C            WRITE(6,*) 'NL2_TOT, ID3_BATLEN, ID4_BATLEN = ',
C    &                   NL2_TOT, ID3_BATLEN, ID4_BATLEN
C            STOP ' Dimension of TSCR5 too small '
C         END IF
C         CALL SETVEC(WORK(KLTSCR5),ZERO,ID3_BATLEN*ID4_BATLEN*NL2_TOT)
C         DO ID2SM = 1, NSMST
C          K2SM = IDBGMULT(ID2SM,IT2SM)
*. The symmetry of L1 obtained from OT1(l1,d4,d3,d2) has symmetry OPT1SM
C
C Again not sure about Lsym
C
           ID234SM = IDBGMULT(ID2SM,ID34SM)
           L1SM = IDBGMULT(IOPT1SM,IADJSYM(ID234SM))
*. Number of L1 strings with sym L1SM
           NL1_TOT = LEN_TCCBLK_REL(NL1(1,1),NL1(1,2),NL1(1,3),
     &                   NL1(1,4),L1SM,NSMST)
*. Offsets for L1 with sym L1SM
           CALL Z_TCC_OFF2_REL(IBL1_TCC,LEN_L1,NL1(1,1),
     &          NL1(1,2),NL1(1,3),NL1(1,4),L1SM,NSMST)
*
           IK2SM_NEW = 1
           ID2SM_NEW = 1
           IF(NTEST.GE.100)  WRITE(6,*) ' ID2SM,K2SM = ',ID2SM, K2SM
*. Number of D2 strings with actual symmetry
           ND2_TOT = LEN_TCCBLK_REL(NO2DX(1,1),NO2DX(1,2),NO2DX(1,3),
     &                   NO2DX(1,4),ID2SM,NSMST)
           IF(NTEST.GE.100) WRITE(6,*) ' ND2_TOT = ', ND2_TOT 
*. Number of K2 strings with actual symmetry
           NK2_TOT = LEN_TCCBLK_REL(NK2(1,1),NK2(1,2),NK2(1,3),
     &                   NK2(1,4),K2SM,NSMST)
           IF(NK2_TOT.GT.0.AND.ND2_TOT.NE.0) THEN
*. Number of D2 batches 
           ND2_BAT = ND2_TOT/LD2_BAT
           IF(ND2_BAT*LD2_BAT.LT.ND2_TOT) ND2_BAT=ND2_BAT+1
*. And loop over batches of D2
       IF(NTEST.GE.100) WRITE(6,*) 'ND2_TOT, ND2_BAT ', ND2_TOT, ND2_BAT
           DO ID2_BAT = 1, ND2_BAT
*. Generate D3 strings for given sym and batch
            ID2_START = (ID2_BAT-1)*LD2_BAT + 1
            ID2_STOP  = MIN(ND2_TOT,ID2_START+LD2_BAT-1)
            ID2_BATLEN = ID2_STOP-ID2_START+1
            IF(NTEST.GE.100) THEN
              WRITE(6,*) ' ID2SM, ID2_BAT, ID2_BATLEN = ',
     &                     ID2SM, ID2_BAT, ID2_BATLEN
            END IF
*. Generate D2 strings for given sym and batch
            CALL ISMNM_FOR_TCC_BAT_REL(NO2DX,ISM_CAAB_D2,INM_CAAB_D2,
     &           ID2SM,
     &           ID2_BATLEN,ID2SM_NEW,
     &           ISM_C1_D2,ISM_CA1_D2,ISM_AA1_D2,
     &           INM_AB1_D2,INM_AA1_D2,INM_CA1_D2,INM_CB1_D2,
     &           ISM_CINI_D2,ISM_CAINI_D2,ISM_AAINI_D2,      
     &           INM_ABINI_D2,INM_AAINI_D2,INM_CAINI_D2,INM_CBINI_D2,0)
            ID2SM_NEW = 0
*
* Thinks this should be new start Lasse
*
* ==============================================================
* Part 1 O1(d1,d2,d3,d4,ex) T1(i) => OT1(d2,d3,d4,l1) in TSCR4
* ==============================================================
*
* Part 1 : Obtain in TSCR4 OT1(d4,d3,d2,l1) for all l1 of correct sym
*
* Length of OT1(d4,d3,d2,l1)
            LEN_OT1 = NL1_TOT*ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
            IF(LEN_OT1 .GT. LEN4 ) THEN
*. Problem, too small length for TSCR4 
             WRITE(6,*) ' Dimension of TSCR4 too small '
             WRITE(6,*) ' Required and allocated : ',LEN_OT1, LEN4 
             WRITE(6,*) 'NL1_TOT,ID2_BATLEN, ID3_BATLEN, ID4_BATLEN = ',
     &                   NL1_TOT,ID2_BATLEN, ID3_BATLEN, ID4_BATLEN
             WRITE(6,*) ' Ex, D2,D3,D4 and L1 : '
             CALL WRT_SPOX_TP_CC_KRCC(IOEX,1)
             CALL WRT_SPOX_TP_CC_KRCC(IO2DX,1)
             CALL WRT_SPOX_TP_CC_KRCC(IO3DX,1)
             CALL WRT_SPOX_TP_CC_KRCC(IO4DX,1)
             CALL WRT_SPOX_TP_CC_KRCC(IOC_L1,1)
             WRITE(6,*) ' T1, T2, T3, T4 '
             CALL WRT_SPOX_TP_CC_KRCC(IT1,1)
             CALL WRT_SPOX_TP_CC_KRCC(IT2,1)
             CALL WRT_SPOX_TP_CC_KRCC(IT3,1)
             CALL WRT_SPOX_TP_CC_KRCC(IT4,1)
             STOP ' Dimension of TSCR4 too small '
            END IF
*
            ZERO = 0.0D0
            CALL SETVEC(WORK(KLTSCR4),ZERO,LEN_OT1)
*. Loop 1 can be realized in a number of different ways. 
*. The present approach minimizes the number of times 
*  integrals must be fetched 
*
            DO ID1SM = 1, NSMST
*. Obtain exsm from OP(ex,d4,d3,d2,d1) have symmetry IOPSM
             IF(NTEST.GE.100) WRITE(6,*) ' ID1SM = ', ID1SM
C
C Not sure about EXsym
C
             ID1234SM   = IDBGMULT(ID1SM,ID234SM)
             IEXSM = IDBGMULT(IOPSM,IADJSYM(ID1234SM))
             K1SM = IDBGMULT(IT1SM,ID1SM)
*. Number of K1 strings with sym K1SM
             NK1_TOT = LEN_TCCBLK_REL(NK1(1,1),NK1(1,2),NK1(1,3),
     &                 NK1(1,4),K1SM,NSMST)
             IF(NTEST.GE.100) WRITE(6,*) 'K1SM, NK1_TOT ,IEXSM = ', 
     &                                    K1SM, NK1_TOT ,IEXSM
*. Number of excitations strings in O with this symmetry
             NEX_TOT = LEN_TCCBLK_REL(NOEX(1,1),NOEX(1,2),NOEX(1,3),
     &                 NOEX(1,4),IEXSM,NSMST)
*. Number of D1 strings with given symmetry
             ND1_TOT = LEN_TCCBLK_REL(NO1DX(1,1),NO1DX(1,2),NO1DX(1,3),
     &              NO1DX(1,4),ID1SM,NSMST)
             IF(ND1_TOT.GT.0.AND.NEX_TOT.GT.0.AND.NK1_TOT.GT.0) THEN
*. Number of batches of excitation part of O
             NEX_BAT = NEX_TOT/LEX_BAT
             IF(NEX_BAT*LEX_BAT.LT.NEX_TOT) NEX_BAT = NEX_BAT + 1
*. Number of D1 batches
             ND1_BAT = ND1_TOT/LD1_BAT
             IF(ND1_BAT*LD1_BAT.LT.ND1_TOT) ND1_BAT=ND1_BAT+1
*. Loop over D1 batches
             ID1SM_NEW = 1
       IF(NTEST.GE.100) WRITE(6,*) 'ND1_TOT, ND1_BAT ', ND1_TOT, ND1_BAT
             DO ID1_BAT = 1, ND1_BAT
              IF(NTEST.GE.100) WRITE(6,*) ' ID1_BAT = ', ID1_BAT
              ID1_START = (ID1_BAT-1)*LD1_BAT + 1
              ID1_STOP  = MIN(ND1_TOT,ID1_START+LD1_BAT-1)
              ID1_BATLEN = ID1_STOP-ID1_START+1
              IF(NTEST.GE.100) 
     &        WRITE(6,*) ' ID1_STOP, ID1_START, ID1_BATLEN ',
     &                     ID1_STOP, ID1_START, ID1_BATLEN
*. Generate ID1 strings for given sym and batch
              CALL ISMNM_FOR_TCC_BAT_REL(NO1DX,ISM_CAAB_D1,INM_CAAB_D1,
     &             ID1SM,ID1_BATLEN,ID1SM_NEW,
     &             ISM_C1_D1,ISM_CA1_D1,ISM_AA1_D1,
     &             INM_AB1_D1,INM_AA1_D1,INM_CA1_D1,INM_CB1_D1,
     &             ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,      
     &             INM_ABINI_D1,INM_AAINI_D1,INM_CAINI_D1,INM_CBINI_D1,
     &             0)
              ID1SM_NEW = 0
*
* all excitation operators in one shot. 
* As we only have limited ISM, INM, the integrals (d1,d2,d3,d4,ex)
* are obtained in batches 
              IEXSM_NEW = 1
              DO IEX_BAT = 1, NEX_BAT
               IF(NTEST.GE.100) WRITE(6,*) ' IEX_BAT = ', IEX_BAT
               IEX_START = (IEX_BAT-1)*LB + 1
               IEX_STOP  = MIN(NEX_TOT,IEX_START+LB-1)
               IEX_BATLEN = IEX_STOP - IEX_START + 1
*. Generate IEX strings for given sym and batch
               CALL ISMNM_FOR_TCC_BAT_REL(NOEX,ISM_CAAB_EX,INM_CAAB_EX,
     &              IEXSM,IEX_BATLEN,IEXSM_NEW,
     &              ISM_C1_EX,ISM_CA1_EX,ISM_AA1_EX,
     &              INM_AB1_EX,INM_AA1_EX,INM_CA1_EX,INM_CB1_EX,
     &              ISM_CINI_EX,ISM_CAINI_EX,ISM_AAINI_EX,      
     &              INM_ABINI_EX,INM_AAINI_EX,INM_CAINI_EX,INM_CBINI_EX,
     &              0)
               IEXSM_NEW = 0
*. Length of EXD2D3D4
               LEN_EXD234 = IEX_BATLEN*ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
*. Obtain integrals OP(D1,D2,D3,D4,EX)
                 IOFF = 1 +
     &           (IEX_START-1)*ID1_BATLEN*ID2_BATLEN*ID3_BATLEN*
     &                         ID4_BATLEN
*
                 CALL GET_OPINT4_REL(OPSCR(IOFF),
     &           IO4DX,ID4_BATLEN,INM_CAAB_D4(1,1),ISM_CAAB_D4(1,1),
     &           ISTR_D4,IBO4DX,
     &           IO3DX,ID3_BATLEN,INM_CAAB_D3(1,1),ISM_CAAB_D3(1,1),
     &           ISTR_D3,IBO3DX,
     &           IO2DX,ID2_BATLEN,INM_CAAB_D2(1,1),ISM_CAAB_D2(1,1),
     &           ISTR_D2,IBO2DX,
     &           IO1DX,ID1_BATLEN,INM_CAAB_D1(1,1),ISM_CAAB_D1(1,1),
     &           ISTR_D1,IBO1DX,
     &           IOEX,IEX_BATLEN,INM_CAAB_EX(1,1),ISM_CAAB_EX(1,1),
     &           ISTR_EX,IBOEX,IEXD1234_INDX)
*
              END DO
*            ^ End of loop over batches of excitation operators OEX
*. We now have OP(d1,d2,d3,d4,ex) for all ex 
*. OP(d1,d2,d3,d4,ex)*T1(I) => OPT(d2,d3,d4,l1) in TSCR4
*
              SIGNX = SP_D1K1*SIGN_D1K1
              SIGN_LK = SIGN_EXK1
              SIGNG = FLOAT(ISIGNG)
              IF(NTEST.GE.100) WRITE(6,*) ' SIGNX, SIGNG, SIGN_LK = ', 
     /                     SIGNX, SIGNG, SIGN_LK
              NDUM = ID2_BATLEN*ID3_BATLEN*ID4_BATLEN
              IF(NTEST.GE.100) 
     &        WRITE(6,*) ' OP(d1,d2,d3,d4,ex)*T1(I) => OPT(d2,d3,d4,l1)'
C              WRITE(6,*) ' SIGNX, SIGNG, SIGN_LK = ',
C     /                     SIGNX, SIGNG, SIGN_LK
              CALL OT_T_REL(OPSCR,T1,WORK(KLTSCR4),
     &           NDUM,K1SM,IEXSM,ID1SM,
     &           NK1,NOEX,NT1,NO1DX,NL1,LB,NK1_TOT,NEX_TOT,LEN_T1,
     &           ID1_BATLEN,
     &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,I_NM_CAAB_L,
     &           ISM_CINI_D1,ISM_CAINI_D1,ISM_AAINI_D1,      
     &           INM_ABINI_D1,INM_AAINI_D1,INM_CBINI_D1,INM_CAINI_D1,
     &           SIGNX,SIGNG,SIGN_LK, 
     &           IB_D1K1,IB_EXK1,
     &           IX1_CA,SX1_CA,IX1_CB,SX1_CB,
     &           IX1_AA,SX1_AA,IX1_AB,SX1_AB,
     &           IX2_CA,SX2_CA,IX2_CB,SX2_CB,
     &           IX2_AA,SX2_AA,IX2_AB,SX2_AB,
     &           IBT1_TCC, IBL1_TCC,
     &           TSCR1,TSCR2,TSCR3)
*
             END DO
*             ^ End of loop over D1 batch
              IF(NTEST.GE.100) THEN 
                WRITE(6,*) ' OT1(D4D3D2,L1) matrix after ID1SM=  ',ID1SM
                ND432 = ID2_BATLEN*ID3_BATLEN*ID4_BATLEN 
                print*,'ND432',ND432
                CALL WRTMAT(WORK(KLTSCR4),ND432,NL1_TOT,ND432,NL1_TOT)
              END IF
*
            END IF
*           ^ End if there were nonvanishing dimensions
            END DO
*           ^ End of loop over symmetries of D1
*
            IF(NTEST.GE.100) THEN 
              WRITE(6,*) ' OT1(D4D3D2,L1) matrix '
              ND432 = ID2_BATLEN*ID3_BATLEN*ID4_BATLEN 
              print*,'ND432',ND432
              CALL WRTMAT(WORK(KLTSCR4),ND432,NL1_TOT,ND432,NL1_TOT)
            END IF
*
* Think this should be the end!!!
*
*
* =========================================================
* Part 2 OT1(d2,d3,d4,l1) T2(i) => OT12(d3,d4,l2) in TSCR5
* =========================================================
*
*  We have OT1(d4,d3,d2,L1) for d4,d3,d2 in batch and all L1 of given symmetry
*  Construct OT12(d4,d3,l2) for d4,d3 in batch and all l2
*
C           LEND34 = ID3_BATLEN*ID4_BATLEN
C           SIGNX = SIGN_D2K2*SP_D2K2
C           SIGNY = 1.0D0
C           SIGNLK = SIGN_L1K2
CM          CALL SETVEC(TSCR5,ZERO,LEND34*NL2_TOT)
C           NTEST = 00
C           IF(NTEST.GE.100) 
C    &      WRITE(6,*) '  OT1(d2,d3,d4,l1) T2(i) => OT12(d3,d4,l2) '
C           NTEST = 00
C           CALL OT_T_REL(WORK(KLTSCR4),T2,WORK(KLTSCR5),LEND34,
C    &       K2SM,L1SM,ID2SM,
C    &       NK2,NL1,NT2,NO2DX,NL2,LB,NK2_TOT,NL1_TOT,LEN_T2,
C    &       ID2_BATLEN,
C    &       ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &       ISM_CINI_D2,ISM_CAINI_D2,ISM_AAINI_D2,
C    &       INM_ABINI_D2,INM_AAINI_D2,INM_CBINI_D2,INM_CAINI_D2,
C    &       SIGNX,SIGNY,SIGNLK,
C    &       IB_D2K2,IB_L1K2,
C    &       IX3_CA,SX3_CA,IX3_CB,SX3_CB,
C    &       IX3_AA,SX3_AA,IX3_AB,SX3_AB,
C    &       IX4_CA,SX4_CA,IX4_CB,SX4_CB,
C    &       IX4_AA,SX4_AA,IX4_AB,SX4_AB,
C    &       IBT2_TCC,IBL2_TCC,TSCR1,TSCR2,TSCR3)  
           END DO
*          ^ End of loop over D2 batches
          END IF 
*         ^ End of NK2_TOT*ND2_TOT .ne. 0
C         END DO
*         ^ End of loop over symmetries of D2 batches
C         IF(NTEST.GE.100) THEN 
C           WRITE(6,*) ' OT12(D4D3,L2) matrix '
C           ND43 = ID3_BATLEN*ID4_BATLEN 
C           CALL WRTMAT(WORK(KLTSCR5),ND43,NL2_TOT,ND43,NL2_TOT)
C         END IF
* We now have OPT12(d3,d4,l2) in TSCR5
*
* Part 3 : OPT123(d4,l3) in TSCR6
*
C         SIGNX = SIGN_D3K3*SP_D3K3
C         SIGNY = 1.0D0
C         SIGNLK = SIGN_L2K3
C         CALL SETVEC(TSCR6,ZERO,ID4_BATLEN*NL3_TOT)
C           NTEST = 00
C         IF(NTEST.GE.100) 
C    &    WRITE(6,*) '  OT12(d3,d4,l2) T3(i) => OT123(d4,l3) '
C           NTEST = 00
C         CALL OT_T_REL(WORK(KLTSCR5),T3,WORK(KLTSCR6),ID4_BATLEN,
C    &       K3SM,L2SM,ID3SM,
C    &       NK3,NL2,NT3,NO3DX,NL3,LB,NK3_TOT,NL2_TOT,LEN_T3,
C    &       ID3_BATLEN,
C    &       ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &       ISM_CINI_D3,ISM_CAINI_D3,ISM_AAINI_D3,
C    &       INM_ABINI_D3,INM_AAINI_D3,INM_CBINI_D3,INM_CAINI_D3,
C    &       SIGNX,SIGNY,SIGNLK,
C    &       IB_D3K3,IB_L2K3,
C    &       IX5_CA,SX5_CA,IX5_CB,SX5_CB,
C    &       IX5_AA,SX5_AA,IX5_AB,SX5_AB,
C    &       IX6_CA,SX6_CA,IX6_CB,SX6_CB,
C    &       IX6_AA,SX6_AA,IX6_AB,SX6_AB,
C    &       IBT3_TCC,IBL3_TCC,TSCR1,TSCR2,TSCR3)  
C        END DO
*        ^ End of loop over D3 batches 
C       END IF
*.      ^ End if there are K3 and L3 strings
C       END DO
*       ^ End of loop over D3 symmetry
C       IF(NTEST.GE.100) THEN 
C          WRITE(6,*) ' OT123(D4,L3) matrix '
C          ND4 = ID4_BATLEN 
C          CALL WRTMAT(WORK(KLTSCR6),ND4,NL3_TOT,ND4,NL3_TOT)
C        END IF
*
* Part 4 : OT123(d4,l3)*T4 => OT1234(I)
C       SIGNX = SIGN_D4K4*SP_D4K4
C       SIGNY = 1.0D0
C       SIGNLK = SIGN_L3K4
C       IF(NTEST.GE.100) 
C    &  WRITE(6,*) '  OT123(d4,l1) T4(i) => OT1234(l4) '
C       CALL OT_T_REL(WORK(KLTSCR6),T4,OT1234,1,
C    &       K4SM,L3SM,ID4SM,
C    &       NK4,NL3,NT4,NO4DX,NL4,LB,NK4_TOT,NL3_TOT,LEN_T4,
C    &       ID4_BATLEN,
C    &       ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
C    &       ISM_CINI_D4,ISM_CAINI_D4,ISM_AAINI_D4,
C    &       INM_ABINI_D4,INM_AAINI_D4,INM_CBINI_D4,INM_CAINI_D4,
C    &       SIGNX,SIGNY,SIGNLK,
C    &       IB_D4K4,IB_L3K4,
C    &       IX7_CA,SX7_CA,IX7_CB,SX7_CB,
C    &       IX7_AA,SX7_AA,IX7_AB,SX7_AB,
C    &       IX8_CA,SX8_CA,IX8_CB,SX8_CB,
C    &       IX8_AA,SX8_AA,IX8_AB,SX8_AB,
C    &       IBT4_TCC,IBL4_TCC,TSCR1,TSCR2,TSCR3)  
C      END DO
*      ^ End of D4 batches 
C     END DO
*     ^ End of D4 symmetries
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Updated OT1234 block '
        CALL WRT_TCC_BLK_CC_REL(OT1234,IL4SM,NL4(1,1),NL4(1,2),NL4(1,3),
     &       NL4(1,4),NSMST)
      END IF
*
C     CALL MEMMAR(IDUM,IDUM,'FLUSM ',IDUM,'OPCTT ')
      RETURN
      END
*
      SUBROUTINE CONTR_ORD4_REL(NDOP,ID1,ID2,ID3,ID4,
     &           LD1,LD2,LD3,LD4,IORD) 
*
* NDOP contraction operators ID1, ID2, ... are given. 
* Find which order these contraction operators should be 
* applied 
*
* Jeppe Olsen, July15, 2001, debugging initialized nov. 2004
*
* Initial version : Contraction containing most elements
* are contracted first
* 
* Output : IORD : New operator I is original operator IORD(I)

*
#include "implicit.inc"
*. Input : Contraction operators 
      DIMENSION ID1(4,3),ID2(4,3),ID3(4,3),ID4(4,3)
*. Output 
      INTEGER IORD(4)
*. Local scratch
      INTEGER LEN_AR(4),INO(4), ISCR(4)
*
*. Length of different contractions
                    CALL DIM_CNTR_REL(ID1,LD1,4,LEN_AR(1))
      IF(NDOP.GE.2) CALL DIM_CNTR_REL(ID2,LD2,4,LEN_AR(2))
      IF(NDOP.GE.3) CALL DIM_CNTR_REL(ID3,LD3,4,LEN_AR(3))
      IF(NDOP.GE.4) CALL DIM_CNTR_REL(ID4,LD4,4,LEN_AR(4))
*. Multiply with -1 so numerically largest numbers 
*. will come first in sort
      DO I = 1, NDOP
       LEN_AR(I) = - LEN_AR(I)
      END DO
C           ORDINT(IINST,IOUTST,NELMNT,INO,IPRNT)
       CALL ORDINT(LEN_AR,ISCR,NDOP,INO,0)   
*. INO now gives new order => old order which is what we want so
      CALL ICOPVE(INO,IORD,NDOP)
*. Fusk, enforce order 1 2 3 4
C     WRITE(6,*) ' Order of contraction is set to 1234'
C     IORD(1) = 1
C     IORD(2) = 2
C     IORD(3) = 3
C     IORD(4) = 4
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Contraction operators  : '
C  WRT_CNTR(ICONT,NCONT,LDIM)
       CALL WRT_CNTR3(ID1,LD1,4)
       IF(NDOP.GE.2) CALL WRT_CNTR3(ID2,LD2,4)
       IF(NDOP.GE.3) CALL WRT_CNTR3(ID3,LD3,4)
       IF(NDOP.GE.4) CALL WRT_CNTR3(ID4,LD4,4)
       WRITE(6,*) ' Order of contraction for operators '
       CALL IWRTMA(IORD,1,NDOP,1,NDOP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE OT_T_REL(OT1,T2,OT1T2,NDUM,
     &           ISM_K,ISM_L,ISM_D,
     &           NK,NL,NI,ND,NLK,LEN_KLBAT,NK_TOT,NL_TOT,NI_TOT,ND_TOT,
     &           ISM_CAAB_K,INM_CAAB_K,ISM_CAAB_L,INM_CAAB_L,
     &           ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,
     &           INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
     &           SIGNX,SIGNG,SIGNLK,
     &           IB_DK,IB_LK,
     &           IX_DK_CA,SX_DK_CA,IX_DK_CB,SX_DK_CB,
     &           IX_DK_AA,SX_DK_AA,IX_DK_AB,SX_DK_AB,
     &           IX_LK_CA,SX_LK_CA,IX_LK_CB,SX_LK_CB,
     &           IX_LK_AA,SX_LK_AA,IX_LK_AB,SX_LK_AB,
     &           IBT2_TCC,IBI_TCC,TSCR1,TSCR2,TSCR3)
* Input / Argument list 
* =====================
*
* OT1 : OT1 coefficients
*  T2 : T2  coefficients
* OT1T2 (OUTPUT) : OT1T2 coefficients
* NDUM : Number of dummy indeces 
* ISM_K : Symmetry of K in T2(D,K)
* ISM_L : Symmetry of L-operator in OT1(IDUM,D,L)
* ISM_D : Symmetry of  D in T2(D,K)
* NK : Number of CA,CB,AA,AB operators per sym for K-strings in T2(D,K)
* NL : Number of CA,CB,AA,AB operators per sym for L-strings in  OT1(IDUM,D,L)
* NI : Number of CA,CB,AA,AB operators per sym for I-strings in  T2(I)
* ND : Number of CA,CB,AA,AB operators per sym for D-strings in OT1(IDUM,D,L)
* NLK : Number of CA,CB,AA,AB operators per sym for LK-strings in OT1T2(IDUM,LK)
* LEN_KLBAT : Length of batches over L and K
* NK_TOT : Total number of K-strings in T2(D,K)
* NI_TOT : Total number of I-strings in T2(I)
* ND_TOT : Total number of D-strings in T2(D,K) ( number of D's in batch)
* ISM_CAAB_K, INM_CAAB_K, ISM_CAAB_L, INM_CAAB_L (SCRATCH) : scratch 
* arrays for CAAB of K and L
* ISM_CINI_L,ISM_CAINI_L,ISM_AAINI_L,INM_ABINI_L,INM_AAINI_L,INM_CBINI_L,
* INM_CAINI_L : Initial value
* SIGNX, SIGNG, SIGNLK : General sign factors 
* IB_DK(IDSM,IKSM,ICAAB) Offset to mappings for D and K ops with given sym 
*                        for ICAAB part of operator
* IB_LK(ILSM,IKSM,ICAAB) Offset to mappings for L and K ops with given sym 
*                        for ICAAB part of operator
* IX_DK_CA, SX_DK_CA, : (D,K) = > I mapping for CA
* IX_DK_CB, SX_DK_CB, : (D,K) = > I mapping for CB
* IX_DK_AA, SX_DK_AA, : (D,K) = > I mapping for AA
* IX_DK_AB, SX_DK_AB, : (D,K) = > I mapping for AB
* IX_LK_CA, SX_LK_CA  : (L,K) => LK mapping for CA
* IX_LK_CB, SX_LK_CB  : (L,K) => LK mapping for CB
* IX_LK_AA, SX_LK_AA  : (L,K) => LK mapping for AA
* IX_LK_AB, SX_LK_AB  : (L,K) => LK mapping for AB
* IBT2_TCC : Offsets to symmetry blocks in T2
* IBI_TCC : Offset to symmetry blocks in T1
* TSCR1, TSCR2,TSCR3 (SCRATCH)
* 

*
* Include next operator in collected operator 
*
* OT1(D,IDUM,L)T2(I) => OT12(IDUM,LK) for all LK,I,L, D in batch
*
* Loop over batches of K
*  T2(I) => T2(D,K)
*  Loop over batches of L
*   OT1T2(IDUM,L,K) = sum(D) OT1(D,IDUM,L)*T2(D,K)  
*   OT1T2(IDUM,L,K) => OT1T2(IDUM,LK)
*  End of loop over L batches 
* End of loop over K batches
*
* Jeppe Olsen, July 14, 2001, completed and debugged July 2003
*                             (Well it became Nov. 2004)
*
#include "implicit.inc"
*. Number of strings per symmetry
      INTEGER NK(8,4),NL(8,4),NI(8,4),ND(8,4),NLK(8,4)
*. Maps 
      INTEGER IB_DK(8,8,4)
      DIMENSION IX_DK_CA(*),SX_DK_CA(*),IX_DK_CB(*),SX_DK_CB(*)
      DIMENSION IX_DK_AA(*),SX_DK_AA(*),IX_DK_AB(*),SX_DK_AB(*)
      INTEGER IB_LK(8,8,4)
      DIMENSION IX_LK_CA(*),SX_LK_CA(*),IX_LK_CB(*),SX_LK_CB(*)
      DIMENSION IX_LK_AA(*),SX_LK_AA(*),IX_LK_AB(*),SX_LK_AB(*)
*. Offset 
      INTEGER IBT2_TCC(8,8,8),IBI_TCC(8,8,8)
*. Input
      DIMENSION OT1(*), T2(*)
*. Output 
      DIMENSION OT1T2(*)
*. Scratch through input 
      INTEGER ISM_CAAB_K(*),INM_CAAB_K(*)
      INTEGER ISM_CAAB_L(*),INM_CAAB_L(*)
      DIMENSION TSCR1(*),TSCR2(*),TSCR3(*)
*
#include "csm.inc"
*
*
*
      NTEST = 00
      ZERO = 0.0D0
      ONE  = 1.0D0
*
C      print*,'NDUM',NDUM
C      print*,'SIGNG',SIGNG
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' OT_T entered '
        WRITE(6,*) ' NK_TOT, NL_TOT, ND_TOT, NI_TOT = ',
     &               NK_TOT, NL_TOT, ND_TOT, NI_TOT
      END IF
*. Number of K-batches
C?    WRITE(6,*) 'NK_TOT, LEN_KLBAT = ', NK_TOT, LEN_KLBAT
      NK_BAT = NK_TOT/LEN_KLBAT
      IF(NK_BAT*LEN_KLBAT.LT.NK_TOT) NK_BAT = NK_BAT + 1
*. Number of batches of L
      NL_BAT = NL_TOT/LEN_KLBAT
      IF(NL_BAT*LEN_KLBAT.LT.NL_TOT) NL_BAT = NL_BAT + 1
*. And loop over batches of K 
      IKSM_NEW = 1
      DO IK_BAT = 1, NK_BAT
        IK_START = (IK_BAT-1)*LEN_KLBAT + 1
        IK_STOP  = MIN(NK_TOT,IK_START+LEN_KLBAT-1)
        IK_BATLEN = IK_STOP-IK_START+1 
        IF(NTEST.GE.100) THEN
          WRITE(6,*) ' IK_BAT, IK_START, IK_STOP = ',
     &                 IK_BAT, IK_START, IK_STOP 
        END IF
*. Generate K strings for given sym and batch
        CALL ISMNM_FOR_TCC_BAT_REL(NK,ISM_CAAB_K,INM_CAAB_K,ISM_K,
     &       IK_BATLEN,IKSM_NEW,
     &       ISM_C1_K,ISM_CA1_K,ISM_AA1_K,
     &       INM_AB1_K,INM_AA1_K,INM_CA1_K,INM_CB1_K,     
     &       ISM_CINI_K,ISM_CAINI_K,ISM_AAINI_K,      
     &       INM_ABINI_K,INM_AAINI_K,INM_CAINI_K,INM_CBINI_K,1)
        IKSM_NEW = 0
*. Obtain T2(D,K) in TSCR1
        IF(NTEST.GE.100) WRITE(6,*) '  T2(I) => T2(D,K) reordering'
C       SIGNX = SIGN_D2K2*SP_D2K2
        CALL TI_TO_TOKBN_REL(NSMST,
     &       ISM_CINI_D,ISM_CAINI_D,ISM_AAINI_D,ISM_D,
     &       INM_ABINI_D,INM_AAINI_D,INM_CBINI_D,INM_CAINI_D,
     &       ND_TOT,ND,
     &       ISM_CINI_K,ISM_CAINI_K,ISM_AAINI_K,ISM_K,
     &       INM_ABINI_K,INM_AAINI_K,INM_CBINI_K,INM_CAINI_K,
     &       IK_BATLEN,NK,
     &       IX_DK_CA,SX_DK_CA,IB_DK(1,1,1),
     &       IX_DK_CB,SX_DK_CB,IB_DK(1,1,2),
     &       IX_DK_AA,SX_DK_AA,IB_DK(1,1,3),
     &       IX_DK_AB,SX_DK_AB,IB_DK(1,1,4),
     &       NI(1,1),NI(1,2),NI(1,3),NI(1,4),
     &       TSCR1,T2,IBT2_TCC,1,1,SIGNX,1)
C        print*,'T2',T2(1)
*. And loop over batches of L
        ILSM_NEW = 1
        DO IL_BAT = 1, NL_BAT
          IL_START = (IL_BAT-1)*LEN_KLBAT + 1
          IL_STOP  = MIN(NL_TOT,IL_START+LEN_KLBAT-1)
          IF(NTEST.GE.100) WRITE(6,*) ' IL_BAT, IL_START, IL_STOP ',
     &                 IL_BAT, IL_START, IL_STOP
          IL_BATLEN = IL_STOP-IL_START+1 
*. Generate L strings for given sym and batch
          CALL ISMNM_FOR_TCC_BAT_REL(NL,ISM_CAAB_L,INM_CAAB_L,ISM_L,
     &         IL_BATLEN,ILSM_NEW,
     &         ISM_C1_L,ISM_CA1_L,ISM_AA1_L,
     &         INM_AB1_L,INM_AA1_L,INM_CA1_L,INM_CB1_L,     
     &         ISM_CINI_L,ISM_CAINI_L,ISM_AAINI_L,      
     &         INM_ABINI_L,INM_AAINI_L,INM_CAINI_L,INM_CBINI_L,1)
          ILSM_NEW = 0
*. OT12(IDUM,L,K) = sum(D) OT1(D,IDUM,L)*T2(D,K) in TSCR3
*. Signg is introduced here !!
          NR_OT12 = NDUM*IL_BATLEN
          NC_OT12 = IK_BATLEN
          NR_OT1 = ND_TOT 
          NC_OT1 = NR_OT12
          NR_T2  = ND_TOT
          NC_T2 = IK_BATLEN
CE        IOFF = (IL_START - 1)*ND_TOT*NL_TOT + 1
          IOFF = (IL_START - 1)*ND_TOT*NDUM + 1
C?        WRITE(6,*) ' IOFF in OT_T ', IOFF
          IF(NTEST.GE.100) THEN
            WRITE(6,*) 'OT1(IOFF) block '
            CALL WRTMAT(OT1(IOFF),NR_OT1,NC_OT1,NR_OT1,NC_OT1)
            WRITE(6,*) ' TSCR1 '
            CALL WRTMAT(TSCR1,NR_T2,NC_T2,NR_T2,NC_T2)
          END IF
          ZERO = 0.0D0
          CALL MATML7(TSCR3,OT1(IOFF),TSCR1, NR_OT12, NC_OT12 ,
     &         NR_OT1,NC_OT1,NR_T2,NC_T2,ZERO,SIGNG,1)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' TSCR fresh from MATML7'
            CALL WRTMAT(TSCR3,NR_OT12,NC_OT12,NR_OT12,NC_OT12)
          END IF
*. Expand OT1T2(IDUM,L,K) to OT1T2(IDUM,I) 
          IF(NTEST.GE.100) 
     &    WRITE(6,*) ' OT1T2(L,K) =>  OT1T2(I) reordering '
C?        WRITE(6,*) ' SIGNLK before call to TI_TO_TOKBN(2) ', SIGNLK
          CALL TI_TO_TOKBN_REL(NSMST,
     &         ISM_CINI_L,ISM_CAINI_L,ISM_AAINI_L,ISM_L,
     &         INM_ABINI_L,INM_AAINI_L,INM_CBINI_L,INM_CAINI_L,
     &         IL_BATLEN,NL,
     &         ISM_CINI_K,ISM_CAINI_K,ISM_AAINI_K,ISM_K,
     &         INM_ABINI_K,INM_AAINI_K,INM_CBINI_K,INM_CAINI_K,
     &         IK_BATLEN,NK,
     &         IX_LK_CA,SX_LK_CA,IB_LK(1,1,1),
     &         IX_LK_CB,SX_LK_CB,IB_LK(1,1,2),
     &         IX_LK_AA,SX_LK_AA,IB_LK(1,1,3),
     &         IX_LK_AB,SX_LK_AB,IB_LK(1,1,4),
     &         NLK(1,1),NLK(1,2),NLK(1,3),NLK(1,4),
     &         TSCR3,OT1T2,IBI_TCC,2,NDUM,SIGNLK,0)
        END DO
*       ^ End of loop over batches of L
      END DO
*.    ^ End of loop over batches of K
*
      RETURN
      END
*
      FUNCTION IPERM_PARITY(IPERM,NELMNT)
*
* Find sign required to bring permutation IPERM 
* of the first NELMNT integers into order.
*
* KISS version
*
* Jeppe Olsen, Oct. 2000
*
#include "implicit.inc"
*. Input
      INTEGER IPERM(NELMNT)
*. Local scratch 
      PARAMETER(MAXELMNT=100)
      INTEGER ISCR(MAXELMNT)
*
      IF(NELMNT.GT.MAXELMNT) THEN
        WRITE(6,*) ' IPERM_PARITY in trouble, NELMNT > MAXELMNT'
        WRITE(6,*) ' NELMNT, MAXELMNT = ',  NELMNT, MAXELMNT
        STOP       ' IPERM_PARITY in trouble, NELMNT > MAXELMNT'
      END IF
*. Ensure that all elements from 1 to NELMNT are included
      I_AM_OKAY = 1
      DO I = 1, NELMNT
        IFOUND = 0
        DO IELMNT = 1, NELMNT
          IF(IPERM(IELMNT).EQ.I) IFOUND = 1
        END DO
        IF(IFOUND.EQ.0) I_AM_OKAY = 0
      END DO
*
      IF(I_AM_OKAY.EQ.0) THEN
        WRITE(6,*) ' Illegal input to IPERM_PARITY'
        CALL IWRTMA(IPERM,1,NPERM,1,NPERM)
        STOP       ' Illegal input to IPERM_PARITY'
      END IF
*
      CALL ICOPVE(IPERM,ISCR,NELMNT)
      LPERM = 0
      DO IELMNT = 1, NELMNT
*. Find IELMNT in ISCR (Cannot occur before IELMNT)
        KELMNT = 0
        DO JELMNT = IELMNT, NELMNT
          IF(ISCR(JELMNT).EQ.IELMNT) KELMNT = JELMNT
        END DO
        DO JELMNT = KELMNT, IELMNT+1, -1
          ISCR(JELMNT) = ISCR(JELMNT-1)
        END DO
        ISCR(IELMNT) = IELMNT
        LPERM = LPERM + KELMNT - IELMNT
C?      WRITE(6,*) ' Updated list for IELMNT =', IELMNT
C?      CALL IWRTMA(ISCR,1,NELMNT,1,NELMNT)
      END DO
      IF(MOD(LPERM,2).EQ.0) THEN
        ISIGN = 1
      ELSE
        ISIGN = -1
      END IF
*
      IPERM_PARITY= ISIGN
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
       WRITE(6,*) ' Permutation : '
       CALL IWRTMA(IPERM,1,NELMNT,1,NELMNT)
       WRITE(6,*) ' Parity of permutation = ', ISIGN
      END IF
*
      RETURN
      END
*
      SUBROUTINE NEW_CAAB_OC(IOCC_L,IOCC_R,IOC_OP,ICE,ILR,NGAS)
*
* ICE = 1 contraction map 
*        IOCC_L = IOC_OP contracted with IOCC_R
* ICE = 2 Excitation map 
*        IOCC_L = IOC_OP multiplied with IOCC_R 
*
* ILR = 1 : Find IOCC_L
* ILR = 2 : Find IOCC_R
*
#include "implicit.inc"
*. Input or output
      INTEGER IOCC_L(NGAS,4),IOC_OP(NGAS,4)
      INTEGER IOCC_R(NGAS,4)
*
      IONE = 1
      MONE = -1
*
* Contraction map
*
      IF(ICE.EQ.1) THEN
*. CA_L = CA_R - AA_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,1),IOCC_R(1,1),IOC_OP(1,3),IONE,MONE,
     &                NGAS)
        ELSE IF (ILR.EQ.2) THEN
          CALL IVCSUM(IOCC_R(1,1),IOCC_L(1,1),IOC_OP(1,3),IONE,IONE,
     &               NGAS)
        END IF
*. CB_OUT = CB_IN -AB_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,2),IOCC_R(1,2),IOC_OP(1,4),IONE,MONE,
     &                NGAS)
        ELSE 
          CALL IVCSUM(IOCC_R(1,2),IOCC_L(1,2),IOC_OP(1,4),IONE,IONE,
     &                NGAS)
        END IF
*. AA_OUT = AA_IN - CA_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,3),IOCC_R(1,3),IOC_OP(1,1),IONE,MONE,
     &                NGAS)
         ELSE 
          CALL IVCSUM(IOCC_R(1,3),IOCC_L(1,3),IOC_OP(1,1),IONE,IONE,
     &                NGAS)
        END IF
*. AB_OUT = AB_IN - CB_OP
        IF(ILR.EQ.1) THEN
          CALL IVCSUM(IOCC_L(1,4),IOCC_R(1,4),IOC_OP(1,2),IONE,MONE,
     &                NGAS)
        ELSE 
          CALL IVCSUM(IOCC_R(1,4),IOCC_L(1,4),IOC_OP(1,2),IONE,IONE,
     &                NGAS)
        END IF
      ELSE IF (ICE.EQ.2) THEN
*
* Excitation map
*
* IL_CA = IR _CA + IOP_CA
        
        IF(ILR.EQ.1) THEN
          DO ICAAB = 1, 4
           CALL IVCSUM(IOCC_L(1,ICAAB),IOCC_R(1,ICAAB),IOC_OP(1,ICAAB),
     &          IONE,IONE,NGAS)
          END DO
        ELSE IF( ILR.EQ.2) THEN
          DO ICAAB = 1, 4
           CALL IVCSUM(IOCC_R(1,ICAAB),IOCC_L(1,ICAAB),IOC_OP(1,ICAAB),
     &          IONE,MONE,NGAS)
          END DO
        END IF
      END IF
*     ^ End of ICE switch 
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 'NEW_CAAB_OC : ICE, ILR = ', ICE,ILR
        WRITE(6,*) ' IOCC_L,IOCC_OP,IOCC_R : '
        CALL IWRTMA(IOCC_L,NGAS,4,NGAS,4)
        CALL IWRTMA(IOC_OP,NGAS,4,NGAS,4)
        CALL IWRTMA(IOCC_R,NGAS,4,NGAS,4)
      END IF
*
      RETURN
      END
*
      SUBROUTINE ZBASE(NVEC,IVEC,NCLASS)
*
*  Some class division exists with NVEC(ICLASS) members in
*  class ICLASS.
*
*  Construct array IVEC(ICLASS) giving first element of
*  class ICLASS in full adressing
*
#include "implicit.inc"
      DIMENSION NVEC(*),IVEC(*)
*
      IVEC(1) = 1
      DO 100 ICLASS = 2,NCLASS
        IVEC(ICLASS) = IVEC(ICLASS-1)+NVEC(ICLASS-1)
  100 CONTINUE
*
      NTEST = 0
      IF( NTEST .NE. 0 ) THEN
        WRITE(6,'(A)') '  ZBASE : NVEC and IVEC '
        WRITE(6,'(A)') '  ===================== '
        CALL IWRTMA(NVEC,1,NCLASS,1,NCLASS)
        CALL IWRTMA(IVEC,1,NCLASS,1,NCLASS)
      END IF
*
      RETURN
      END
*
      FUNCTION LEN_TCCBLK_REL(NCA,NCB,NAA,NAB,ITSYM,NSMST)
*
* Find length of Coupled cluster excitation block,
* given types, total symmetry is ITSYM
*
* Jeppe Olsen, May of 2000
* Changed from NSMST**3 to NSMST**2 procedure, July 2002
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "multd2h.inc"
#include "symm.inc"
*. Input
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Local scratch
      INTEGER NC(8), NA(8)
*
*. Number of creation and annihilation strings per symmetry
      IZERO = 0
      CALL ISETVC(NC,IZERO,NSMST)
      CALL ISETVC(NA,IZERO,NSMST)
      DO IA_SM = 1, NSMST
        DO IB_SM = 1, NSMST
          IAB_SM = IDBGMULT(IA_SM,IB_SM)
C          IAB_SM = IDBGMULT(IA_SM,IADJSYM(IB_SM))
          NC(IAB_SM) = NC(IAB_SM) + NCA(IA_SM)*NCB(IB_SM)
C
          IAB_SM = IDBGMULT(IADJSYM(IA_SM),IADJSYM(IB_SM))
C          IAB_SM = IDBGMULT(IA_SM,IADJSYM(IB_SM))
C          IAB_SM = IDBGMULT(IA_SM,IB_SM)
C          NA(IAB_SM) = NA(IAB_SM) + NAA(IA_SM)*NAB(IB_SM)
      NA(IAB_SM) = NA(IAB_SM) + NAA(IADJSYM(IA_SM))*NAB(IADJSYM(IB_SM))
C      NA(IAB_SM) = NA(IAB_SM) + NAA(IA_SM)*NAB(IB_SM)
        END DO
      END DO
*. And number of strings with sym ITSYM
      LEN = 0
      DO ISM_C = 1, NSMST
C        ISM_A = IDBGMULT(ISM_C,ITSYM)
        ISM_A = IDBGMULT(ITSYM,IADJSYM(ISM_C))
*. Number of C and A strings with this sym
        LEN = LEN + NC(ISM_C)*NA(ISM_A)
      END DO
C
C Rewrite
C
C      IZERO = 0
C      CALL ISETVC(NC,IZERO,NSMST)
C      CALL ISETVC(NA,IZERO,NSMST)
C      DO ISM_C = 1, NSMST
C        ISM_A = IDBGMULT(ITSYM,INVELM(ISM_C))
*
      LEN_TCCBLK_REL = LEN
*
      RETURN
      END
*
      SUBROUTINE Z_TCC_OFF2_REL(IBT,LENGTH,NCA,NCB,NAA,NAB,ITSYM,NSMST)
*
* Offsets for symmetryblocks of TCC elements, sym of CA,CB,AA used
*
* Jeppe Olsen, Summer of 99
*
* Compared to Z_TCC_OFF : Total length added to argument list
#include "implicit.inc"
#include "mxpdim.inc"
#include "multd2h.inc"
#include "symm.inc"
*. Input
      INTEGER NCA(*),NCB(*),NAA(*),NAB(*)
*. Output
      INTEGER IBT(8,8,8)
*
      IOFF = 1
      DO ISM_C = 1, NSMST
C        ISM_A = IDBGMULT(ISM_C,ITSYM)
        ISM_A = IDBGMULT(ITSYM,INVELM(ISM_C))
        DO ISM_CA = 1, NSMST
          ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
          DO ISM_AA = 1, NSMST
            ISM_AB =  IDBGMULT(ISM_A,INVELM(ISM_AA))
            IBT(ISM_CA,ISM_CB,ISM_AA) = IOFF
            IOFF = IOFF +
     &      NCA(ISM_CA)*NCB(ISM_CB)*NAA(ISM_AA)*NAB(ISM_AB)
          END DO
        END DO
      END DO
      LENGTH = IOFF - 1
*
      RETURN
      END
*
      SUBROUTINE T1T2_TO_T12_MAP_REL(I1SPOBEX,I2SPOBEX,I12SPOBEX,
     &                       ICA_MAP,XCA_MAP,ICB_MAP,XCB_MAP,
     &                       IAA_MAP,XAA_MAP,IAB_MAP,XAB_MAP,
     &                       IB,I1OCC,I2OCC,I1REO,IZ,IZSCR,SIGNP) 
*
* Obtain mappings for T1T2 => T12 as 4 mappings for each CAAB component
*
* All creation and annihilation operators are assumed to commute.
*
* Jeppe Olsen, Oct 2000     
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "cgas.inc"
*
*. Input
*
      INTEGER I1SPOBEX(NGAS,4)
      INTEGER I2SPOBEX(NGAS,4)
      INTEGER I12SPOBEX(NGAS,4)
*
*. Output
*
*. T1*T2 => T12 mappings
      INTEGER ICA_MAP(*),ICB_MAP(*),IAA_MAP(*),IAB_MAP(*)               
      DIMENSION XCA_MAP(*),XCB_MAP(*),XAA_MAP(*),XAB_MAP(*)         
      INTEGER IB(8,8,4)
*
*. Local scratch 
*
      INTEGER NI1CA(MXPNGAS),NI1CB(MXPNGAS),
     &        NI1AA(MXPNGAS),NI1AB(MXPNGAS)
      INTEGER NI2CA(MXPNGAS),NI2CB(MXPNGAS),
     &        NI2AA(MXPNGAS),NI2AB(MXPNGAS)
      INTEGER NI12CA(MXPNGAS),NI12CB(MXPNGAS),
     &        NI12AA(MXPNGAS),NI12AB(MXPNGAS)
      INTEGER I1CA_EXP(100),I1CB_EXP(100),I1AA_EXP(100),I1AB_EXP(100)
      INTEGER I2CA_EXP(100),I2CB_EXP(100),I2AA_EXP(100),I2AB_EXP(100)
      INTEGER I12CA_EXP(100),I12CB_EXP(100),
     &        I12AA_EXP(100),I12AB_EXP(100)
C     INTEGER ICA(100)
      INTEGER I1OFF(8,8,8), I2OFF(8,8,8), I12OFF(8,8,8)
C     INTEGER IB_CA_SS(8,8),IB_CB_SS(8,8),IB_AA_SS(8,8),IB_AB_SS(8,8)
*
*. Scratch through parameter list
*
*. For occ of two sets of excitation strings and an reorder array 
      INTEGER I1OCC(*),I2OCC(*),I1REO(*)
*. For Z-matrix and its construction
      INTEGER IZ(*),IZSCR(*)
*
      IDUM = 0
      CALL MEMMAR(IDUM,IDUM,'MARK  ',IDUM,'TT_TM ')
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' Welcome to T1T2_TO_T12_MAP '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Excitation operator 1 :'
        WRITE(6,*)
C   WRT_SPOX_TP_CC_REL(IEX_TP,NEX_TP)
        CALL WRT_SPOX_TP_CC_KRCC(I1SPOBEX,1)
        WRITE(6,*) ' Excitation operator 2 :'
        WRITE(6,*)
        CALL WRT_SPOX_TP_CC_KRCC(I2SPOBEX,1)
        WRITE(6,*) ' Compound operator '
        CALL WRT_SPOX_TP_CC_KRCC(I12SPOBEX,1)
      END IF
*. Number of operators in each string 
      N1CA_OP = IELSUM(I1SPOBEX(1,1),NGAS)
      N1CB_OP = IELSUM(I1SPOBEX(1,2),NGAS)
      N1AA_OP = IELSUM(I1SPOBEX(1,3),NGAS)
      N1AB_OP = IELSUM(I1SPOBEX(1,4),NGAS)
*
      N2CA_OP = IELSUM(I2SPOBEX(1,1),NGAS)
      N2CB_OP = IELSUM(I2SPOBEX(1,2),NGAS)
      N2AA_OP = IELSUM(I2SPOBEX(1,3),NGAS)
      N2AB_OP = IELSUM(I2SPOBEX(1,4),NGAS)
*
      N12CA_OP = IELSUM(I12SPOBEX(1,1),NGAS)
      N12CB_OP = IELSUM(I12SPOBEX(1,2),NGAS)
      N12AA_OP = IELSUM(I12SPOBEX(1,3),NGAS)
      N12AB_OP = IELSUM(I12SPOBEX(1,4),NGAS)
*. Sign to bring O(ca1)O(cb1)O(aa1)O(ab1)
*.               O(ca2)O(cb2)O(aa2)O(ab2)
*. into          O(ca1)O(ca2)O(cb1)O(cb2)O(aa1)O(aa2)O(ab1)O(ab2)
*
      NPERM =  N1AB_OP*( N2CA_OP+N2CB_OP+N2AA_OP)
     &      +  N1AA_OP*( N2CA_OP+N2CB_OP )
     &      +  N1CB_OP*  N2CA_OP
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGNP = -1.0D0
      ELSE
        SIGNP = 1.0D0
      END IF

*
      MX_NOP = MAX(N12CA_OP,N12CB_OP,N12AA_OP,N12AB_OP)
*
C     CALL REF_OP_CC(ICA,ICA_EXP,NCA_OP,NGAS,1)
      CALL REF_OP_CC(I1SPOBEX(1,1),I1CA_EXP,N1CA_OP,NGAS,1)
      CALL REF_OP_CC(I1SPOBEX(1,2),I1CB_EXP,N1CB_OP,NGAS,1)
      CALL REF_OP_CC(I1SPOBEX(1,3),I1AA_EXP,N1AA_OP,NGAS,1)
      CALL REF_OP_CC(I1SPOBEX(1,4),I1AB_EXP,N1AB_OP,NGAS,1)
*
      CALL REF_OP_CC(I2SPOBEX(1,1),I2CA_EXP,N2CA_OP,NGAS,1)
      CALL REF_OP_CC(I2SPOBEX(1,2),I2CB_EXP,N2CB_OP,NGAS,1)
      CALL REF_OP_CC(I2SPOBEX(1,3),I2AA_EXP,N2AA_OP,NGAS,1)
      CALL REF_OP_CC(I2SPOBEX(1,4),I2AB_EXP,N2AB_OP,NGAS,1)
*
      CALL REF_OP_CC(I12SPOBEX(1,1),I12CA_EXP,N12CA_OP,NGAS,1)
      CALL REF_OP_CC(I12SPOBEX(1,2),I12CB_EXP,N12CB_OP,NGAS,1)
      CALL REF_OP_CC(I12SPOBEX(1,3),I12AA_EXP,N12AA_OP,NGAS,1)
      CALL REF_OP_CC(I12SPOBEX(1,4),I12AB_EXP,N12AB_OP,NGAS,1)
C      do i=1,4
C        print*,'I1SPOBEX,I2SPOBEX,I12SPOBEX',
C     &          I1SPOBEX(1,i),I2SPOBEX(1,i),I12SPOBEX(1,i)
C      end do
C      print*,N1CA_OP,N1CB_OP,N1AA_OP,N1AB_OP
C      print*,N2CA_OP,N2CB_OP,N2AA_OP,N2AB_OP
C      print*,N12CA_OP,N12CB_OP,N12AA_OP,N12AB_OP
*. Number of strings per sym
      IUB = 1
      ICA = 1
      CALL NST_SPGP_CC_REL(I1SPOBEX(1,1),NI1CA,IUB,ICA)
      IUB = 2
      ICA = 1
      CALL NST_SPGP_CC_REL(I1SPOBEX(1,2),NI1CB,IUB,ICA)
      IUB = 1
      ICA = 2
      CALL NST_SPGP_CC_REL(I1SPOBEX(1,3),NI1AA,IUB,ICA)
      IUB = 2
      ICA = 2
      CALL NST_SPGP_CC_REL(I1SPOBEX(1,4),NI1AB,IUB,ICA)
*
      IUB = 1
      ICA = 1
      CALL NST_SPGP_CC_REL(I2SPOBEX(1,1),NI2CA,IUB,ICA)
      IUB = 2
      ICA = 1
      CALL NST_SPGP_CC_REL(I2SPOBEX(1,2),NI2CB,IUB,ICA)
      IUB = 1
      ICA = 2
      CALL NST_SPGP_CC_REL(I2SPOBEX(1,3),NI2AA,IUB,ICA)
      IUB = 2
      ICA = 2
      CALL NST_SPGP_CC_REL(I2SPOBEX(1,4),NI2AB,IUB,ICA)
*
      IUB = 1
      ICA = 1
      CALL NST_SPGP_CC_REL(I12SPOBEX(1,1),NI12CA,IUB,ICA)
      IUB = 2
      ICA = 1
      CALL NST_SPGP_CC_REL(I12SPOBEX(1,2),NI12CB,IUB,ICA)
      IUB = 1
      ICA = 2
      CALL NST_SPGP_CC_REL(I12SPOBEX(1,3),NI12AA,IUB,ICA)
      IUB = 2
      ICA = 2
      CALL NST_SPGP_CC_REL(I12SPOBEX(1,4),NI12AB,IUB,ICA)
C      print*,'MXPNGAS',MXPNGAS
C      DO i=1,MXPNGAS
C      print*,'NI1CA,NI1CB,NI1AA,NI1AB',
C     &        NI1CA(i),NI1CB(i),NI1AA(i),NI1AB(i)
C      print*,'NI2CA,NI2CB,NI2AA,NI2AB',
C     &        NI2CA(i),NI2CB(i),NI2AA(i),NI2AB(i)
C      print*,'NI12CA,NI12CB,NI12AA,NI12AB',
C     &        NI12CA(i),NI12CB(i),NI12AA(i),NI12AB(i)
C      END DO
C      print*,I1SPOBEX(1,1),I1SPOBEX(1,2),I1SPOBEX(1,3),I1SPOBEX(1,4)
C      print*,I2SPOBEX(1,1),I2SPOBEX(1,2),I2SPOBEX(1,3),I2SPOBEX(1,4)
C      print*,I12SPOBEX(1,1),I12SPOBEX(1,2)
C      print*,I12SPOBEX(1,3),I12SPOBEX(1,4)
*. And offsets
C          Z_TCC_OFF(IBT,NCA,NCB,NAA,NAB,ITSYM,NSMST)
C     CALL Z_TCC_OFF2(I1OFF,LT1,I1SPOBEX(1,1),I1SPOBEX(1,2),
C    &     I1SPOBEX(1,3),I1SPOBEX(1,4),IT1SM,NSMST)   
C     CALL Z_TCC_OFF2(I2OFF,LT2,I2SPOBEX(1,1),I2SPOBEX(1,2),
C    &     I2SPOBEX(1,3),I2SPOBEX(1,4),IT2SM,NSMST)   
C     CALL Z_TCC_OFF2(I12OFF,LT11,I12SPOBEX(1,1),I12SPOBEX(1,2),
C    &     I12SPOBEX(1,3),I12SPOBEX(1,4),IT12SM,NSMST)   
*
*  ================
*. T1 * T2 mappings 
*  ================
C      WRITE(6,*) ' Memory check before STST_TO_ST_MAP'
C      CALL LMEMCHK
C      WRITE(6,*) ' Memcheck passed '
*
C          STST_TO_ST_MAP(IS1OC,IS2OC,IS12OC,
C    &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
C    &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO)
C      print*,'before STST_TO_ST_MAP_REL'
C
C     note must be adapted to rel formalism
C
C     ADDED IUB
C      print*,'STST_TO_ST_MAP_REL says please adapt me'
C      do i=1,4
C       print*,I1SPOBEX(1,i),I2SPOBEX(1,i),I12SPOBEX(1,i)
C      end do
C
C     IUB looks a little funny but works
C
*
*. CA
      IUB = 1 
      CALL STST_TO_ST_MAP_REL(I1SPOBEX(1,1),I2SPOBEX(1,1),
     &     I12SPOBEX(1,1), IB(1,1,1),ICA_MAP,XCA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO,IUB)
*. CB
      IUB = 2
      CALL STST_TO_ST_MAP_REL(I1SPOBEX(1,2),I2SPOBEX(1,2),
     &     I12SPOBEX(1,2), IB(1,1,2),ICB_MAP,XCB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO,IUB)
*. AA ( Jeppe Phase factor due to def of Annistrings must be added)
      IUB = 2
      CALL STST_TO_ST_MAP_REL(I1SPOBEX(1,3),I2SPOBEX(1,3),
     &     I12SPOBEX(1,3), IB(1,1,3),IAA_MAP,XAA_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO,IUB)
*. AB ( Jeppe Phase factor due to def of Annistrings must be added)
      IUB = 1
      CALL STST_TO_ST_MAP_REL(I1SPOBEX(1,4),I2SPOBEX(1,4),
     &     I12SPOBEX(1,4), IB(1,1,4),IAB_MAP,XAB_MAP,
     &     IZ,IZSCR,I1OCC,I2OCC,I1REO,IUB)
*
      CALL MEMMAR(IDUM,IDUM,'FLUSM ',IDUM,'TT_TM ')
      RETURN 
      END
*
      SUBROUTINE STR_CAAB_REL(ICAAB,ISTR_CAAB)
*
* Obtain strings of CAAB types given by ICAAB  for all symmetries 
*
* Jeppe Olsen, May 2000
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "ctcc.inc"
#include "lucinp.inc"
#include "csm.inc"

*. Input
      INTEGER ICAAB(NGAS,4)
*. Output
      INTEGER ISTR_CAAB(MX_ST_TSOSO_BLK_MX*NSMST,4)
C     INTEGER IBSTR_CAAB(NSMST,4)
*. Local scratch
      INTEGER IGRP_AR(MXPNGAS)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' STR_CAAB in action : '
        WRITE(6,*) ' ====================='
        WRITE(6,*)
        WRITE(6,*) 
     &  ' Strings will be generated for the four CAAB substrings of '
        CALL WRT_SPOX_TP_CC_KRCC(ICAAB,1)
      END IF
*
      IDUM = 0
      IBSTR_CA = 1
      IBSTR_CB = 1
      IBSTR_AA = 1
      IBSTR_AB = 1
*CA strings
C      print*,'before entering'
C      DO I=1,5
C       DO J=1,4
C         print*,'ISTR_CAAB(I,J)',ISTR_CAAB(I,J),I,J
C       END DO
C      END DO
      DO ISTSM = 1, NSMST
        CALL OCC_TO_GRP_CC_REL(ICAAB(1,1),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,1),NGAS) 
C        WRITE(6,*) ' NEL(CA) = ', NEL
        IUB = 1 
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_AR,NGAS,ISTSM,NEL,NSTR_CA,
     &      ISTR_CAAB((IBSTR_CA-1)*NEL+1,1), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,1) = IBSTR_CA
        IBSTR_CA = IBSTR_CA + NSTR_CA
*. CB strings
        CALL OCC_TO_GRP_CC_REL(ICAAB(1,2),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,2),NGAS) 
C        WRITE(6,*) ' NEL(CB) = ', NEL
        IUB = 2 
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_AR,NGAS,ISTSM,NEL,NSTR_CB,
     &       ISTR_CAAB((IBSTR_CB-1)*NEL+1,2), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,2) = IBSTR_CB
        IBSTR_CB = IBSTR_CB + NSTR_CB
*. AA strings
        CALL OCC_TO_GRP_CC_REL(ICAAB(1,3),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,3),NGAS) 
C        WRITE(6,*) ' NEL(AA) = ', NEL
        IUB = 1
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_AR,NGAS,ISTSM,NEL,NSTR_AA,
     &       ISTR_CAAB((IBSTR_AA-1)*NEL+1,3), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,3) = IBSTR_AA
        IBSTR_AA = IBSTR_AA + NSTR_AA
*. AB strings
        CALL OCC_TO_GRP_CC_REL(ICAAB(1,4),IGRP_AR,1)
        NEL = IELSUM(ICAAB(1,4),NGAS) 
C        WRITE(6,*) ' NEL(AB) = ', NEL
        IUB = 2
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,IGRP_AR,NGAS,ISTSM,NEL,NSTR_AB,
     &       ISTR_CAAB((IBSTR_AB-1)*NEL+1,4), NOCOB,0,IDUM,IDUM)
COLD    IBSTR_CAAB(ISTSM,4) = IBSTR_AB
        IBSTR_AB = IBSTR_AB + NSTR_AB
*
      END DO
*
C      print*,'after exit'
C      DO I=1,5
C       DO J=1,4
C         print*,'ISTR_CAAB(I,J)',ISTR_CAAB(I,J),I,J
C       END DO
C      END DO
      RETURN
      END
*
      SUBROUTINE NST_SPGP_CC_REL(IOCC,NSTFSM,IUB,ICA)
* Number of strings for given supergroup.
*.Input supergroup is defined by  occupation in each orb space
*
* Jeppe Olsen , March 99
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION IOCC(NGAS)
*. Output
      DIMENSION NSTFSM(NSMST)
*. General input
#include "mxpdim.inc"
#include "gasstr.inc"
#include "csm.inc"
#include "cgas.inc"
#include "symm.inc"
*. Scratch
      INTEGER ISM(MXPNGAS),MNSM(MXPNGAS),MXSM(MXPNGAS)
      INTEGER IGRP(MXPNGAS)
!
!
!
! Implement ICA change ism by iadjsym for anni. Done but not sure if there is bugs!!!
!
! Check if NSTFSMGP and NSTFSMGP2 are used correctly!!!
! 
*
      NTEST = 00
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' NST_SPGP_CC_REL is speaking '
        WRITE(6,*) ' ============================'
*
        WRITE(6,*) ' Occupation of supergroup '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
      END IF
C      PRINT*,'IUB,ICA',IUB,ICA
*. Occupation => groups
C  OCC_TO_GRP_CC(NSTR,IOCC,IGRP,IWAY)
      CALL OCC_TO_GRP_CC_REL(IOCC,IGRP,1)
      IZERO = 0
      CALL ISETVC(NSTFSM,IZERO,NSMST)
*. Max and Min for allowed symmetries
      DO IGAS = 1, NGAS
        MXSM(IGAS) = 1
        IF(IUB.EQ.2) THEN
          DO ISYM = 1, NSMST
            IF(ICA.EQ.1) THEN
              IF(NSTFSMGP2(ISYM,IGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
            ELSE
C              IF(NSTFSMGP2(IADJSYM(ISYM),IGRP(IGAS)) .NE. 0 ) THEN
              IF(NSTFSMGP2(ISYM,IGRP(IGAS)) .NE. 0 ) THEN
                MXSM(IGAS) = IADJSYM(ISYM)
C                MXSM(IGAS) = ISYM
              END IF
            END IF
C          print*,'IGRP(IGAS)',IGRP(IGAS)
C          print*,'NSTFSMGP2(ISYM,IGRP(IGAS)',NSTFSMGP2(ISYM,IGRP(IGAS))
          END DO
          MNSM(IGAS) = NSMST
          DO ISYM = NSMST,1, -1
            IF(ICA.EQ.1) THEN
              IF(NSTFSMGP2(ISYM,IGRP(IGAS)) .NE. 0 )  MNSM(IGAS) = ISYM
            ELSE
C          ISYM2 = IADJSYM(ISYM)
C              IF(NSTFSMGP2(IADJSYM(ISYM),IGRP(IGAS)) .NE. 0 ) THEN
              IF(NSTFSMGP2(ISYM,IGRP(IGAS)) .NE. 0 ) THEN
                MNSM(IGAS) = IADJSYM(ISYM)
C                MNSM(IGAS) = ISYM
              END IF
C          print*,'ISYM,IADJSYM(ISYM)',ISYM,IADJSYM(ISYM)
            END IF
C          print*,'IGRP(IGAS)',IGRP(IGAS)
C          print*,'NSTFSMGP2(ISYM,IGRP(IGAS)',NSTFSMGP2(ISYM,IGRP(IGAS))
          END DO
        ELSE
          DO ISYM = 1, NSMST
            IF(ICA.EQ.1) THEN
              IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) MXSM(IGAS) = ISYM
            ELSE
C              IF(NSTFSMGP(IADJSYM(ISYM),IGRP(IGAS)) .NE. 0 ) THEN
              IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) THEN
                MXSM(IGAS) = IADJSYM(ISYM)
C                MXSM(IGAS) = ISYM
              END IF
            END IF
C          print*,'IGRP(IGAS)',IGRP(IGAS)
C          print*,'NSTFSMGP(ISYM,IGRP(IGAS)',NSTFSMGP(ISYM,IGRP(IGAS))
          END DO
          MNSM(IGAS) = NSMST
          DO ISYM = NSMST,1, -1
            IF(ICA.EQ.1) THEN
              IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) MNSM(IGAS) = ISYM
            ELSE
C              IF(NSTFSMGP(IADJSYM(ISYM),IGRP(IGAS)) .NE. 0 ) THEN
              IF(NSTFSMGP(ISYM,IGRP(IGAS)) .NE. 0 ) THEN
                MNSM(IGAS) = IADJSYM(ISYM)
C                MNSM(IGAS) = ISYM
              END IF
            END IF
C          print*,'IGRP(IGAS)',IGRP(IGAS)
C          print*,'NSTFSMGP(ISYM,IGRP(IGAS)',NSTFSMGP(ISYM,IGRP(IGAS))
          END DO
        END IF
      END DO
*. Last space with more than one symmetry
      NGASL = 1
      DO IGAS = 1, NGAS
C        print*,'MXSM,MNSM,IGAS',MXSM(IGAS),MNSM(IGAS),IGAS
        IF(MXSM(IGAS).NE.MNSM(IGAS)) NGASL = IGAS
      END DO
*. First symmetry combination
      DO IGAS = 1, NGAS
         ISM(IGAS) = MNSM(IGAS)
      END DO
*. Loop over symmetries in each gas space
      IFIRST = 1
 1000 CONTINUE
      IF(IFIRST.EQ.1) THEN
        CALL ISETVC(ISM,1,NGAS)
        IFIRST = 0
        NONEW = 0
      ELSE
        CALL NXTNUM3(ISM,NGAS,MNSM,MXSM,NONEW)
      END IF
      IF(NONEW.EQ.0) THEN
*. Symmetry of current combination and number of strings in this supergroup
        ISMSPGP = ISM(1)
C        print*,'ISMSPGP',ISMSPGP
        IF(IUB.EQ.2) THEN
          IF(ICA.EQ.1) THEN
            NST = NSTFSMGP2(ISM(1),IGRP(1))
          ELSE
            NST = NSTFSMGP(ISM(1),IGRP(1))
          END IF
C        print*,'NST first',NST
          DO JGRP = 2, NGAS
C          print*,'ISM(JGRP)',ISM(JGRP)
            CALL SYMCOM_REL(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
            ISMSPGP = ISMSPGPO
            IF(ICA.EQ.1) THEN
              NST = NST * NSTFSMGP2(ISM(JGRP),IGRP(JGRP))
            ELSE
              NST = NST * NSTFSMGP(ISM(JGRP),IGRP(JGRP))
            END IF
          END DO
        ELSE
          NST = 1
          IF(ICA.EQ.1) THEN
            NST = NSTFSMGP(ISM(1),IGRP(1))
          ELSE
            NST = NSTFSMGP2(ISM(1),IGRP(1))
          END IF
C        print*,'NST first',NST
          DO JGRP = 2, NGAS
C          print*,'ISM(JGRP)',ISM(JGRP)
C          IF(ICA.EQ.2) ISM(JGRP) = IADJSYM(ISM(JGRP))
C          print*,'ISM(JGRP)',ISM(JGRP)
            CALL SYMCOM_REL(3,7,ISMSPGP,ISM(JGRP),ISMSPGPO)
            ISMSPGP = ISMSPGPO
            IF(ICA.EQ.1) THEN
              NST = NST * NSTFSMGP(ISM(JGRP),IGRP(JGRP))
            ELSE
              NST = NST * NSTFSMGP2(ISM(JGRP),IGRP(JGRP))
            END IF
          END DO
        END IF
        NSTFSM(ISMSPGP) =   NSTFSM(ISMSPGP) + NST
        GOTO 1000
      END IF
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*)
     &  ' Number of strings per symmetry for supergroup'
        CALL IWRTMA10(NSTFSM(1),1,NSMST,1,NSMST)
      END IF
*
      RETURN
      END
*
      SUBROUTINE ISMNM_FOR_TCC_BAT_REL(NS_CAAB,ISM_CAAB,INM_CAAB,ISM,
     &           LBAT,ISM_INI,
     &           ISM_C1,ISM_CA1,ISM_AA1,
     &           INM_AB1,INM_AA1,INM_CA1,INM_CB1,
     &           ISM_CINI,ISM_CAINI,ISM_AAINI,
     &           INM_ABINI,INM_AAINI,INM_CAINI,INM_CBINI,IONLY_LIM)     
*
* A Coupled cluster block is given with dimensions N_CAAB(ISM,I_CAAB) for 
* the four substrings.
* Obtain symmetries and string numbers of the four substrings for  
* a batch of strings with total sym ISM. 
* The last string in previous batch is defined 
* by  ISM_C1,ISM_CA1,ISM_AA1,
*     INM_AB1,INM_AA1,INM_CA1,INM_CB1
*
* IF IONLY_LIM = 1, only the limits are returned
*
* If ISM_INI = 1, this is first batch with this symmetry
* 
*
* Jeppe Olsen, May of 2001  
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "orbinp.inc"
#include "cc_exc.inc"
#include "symm.inc"
*. Specific input
      INTEGER NS_CAAB(8,4)             
*. Output
      INTEGER ISM_CAAB(4,*), INM_CAAB(4,*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) '  ISMNM_FOR_TCC_BAT speaking : '
        WRITE(6,*)
      END IF
*
      IF(ISM_INI.EQ.1) THEN
*. Initialize start indeces
       ISM_C1 = 1
       ISM_CA1 = 1
       ISM_AA1 = 1
*
       INM_AB1 = 1
       INM_AA1 = 1
       INM_CB1 = 1
       INM_CA1 = 0
      END IF
*. Save initial indeces
      ISM_CINI  =  ISM_C1
      ISM_CAINI =  ISM_CA1
      ISM_AAINI =  ISM_AA1 
*
      INM_ABINI =  INM_AB1
      INM_AAINI =  INM_AA1 
      INM_CBINI =  INM_CB1 
      INM_CAINI =  INM_CA1
*
      IT =  0
      INI_LOOP = 1
      DO ISM_C = ISM_C1, NSMST
        ISM_A = IDBGMULT(ISM,INVELM(ISM_C)) 
        IF(INI_LOOP .EQ. 1 ) THEN
          ISM_CA2 = ISM_CA1
        ELSE
          ISM_CA2 = 1
        END IF
        DO ISM_CA = ISM_CA2, NSMST
          ISM_CB = IDBGMULT(ISM_C,INVELM(ISM_CA))
          IF(INI_LOOP .EQ. 1 ) THEN
            ISM_AA2 = ISM_AA1
          ELSE
            ISM_AA2 = 1
          END IF
          DO ISM_AA = ISM_AA2, NSMST
            ISM_AB =  IDBGMULT(ISM_A,INVELM(ISM_AA))
*. Loop over operators as  matrix (I_CA,I_CB,I_AA,I_AB)
            NSTR_CA = NS_CAAB(ISM_CA,1)
            NSTR_CB = NS_CAAB(ISM_CB,2)
            NSTR_AA = NS_CAAB(ISM_AA,3)
            NSTR_AB = NS_CAAB(ISM_AB,4)
*
            IF(INI_LOOP.EQ.1) THEN
             INM_AB2 = INM_AB1
            ELSE
             INM_AB2 = 1
            END IF
            DO I_AB = INM_AB2, NSTR_AB
             IF(INI_LOOP.EQ.1) THEN
              INM_AA2 = INM_AA1
             ELSE
              INM_AA2 = 1
             END IF
             DO I_AA = INM_AA2, NSTR_AA
              IF(INI_LOOP.EQ.1) THEN
                INM_CB2 = INM_CB1
              ELSE 
                INM_CB2 = 1
              END IF
              DO I_CB = INM_CB2, NSTR_CB
               IF(INI_LOOP.EQ.1) THEN
                 INM_CA2 = INM_CA1
               ELSE 
                 INM_CA2 = 0
               END IF
               DO I_CA = INM_CA2 + 1, NSTR_CA
                IT = IT + 1
*
                IF(IONLY_LIM.EQ.0) THEN
                  INM_CAAB(1,IT) = I_CA
                  INM_CAAB(2,IT) = I_CB
                  INM_CAAB(3,IT) = I_AA
                  INM_CAAB(4,IT) = I_AB
*
                  ISM_CAAB(1,IT) = ISM_CA
                  ISM_CAAB(2,IT) = ISM_CB
                  ISM_CAAB(3,IT) = ISM_AA
                  ISM_CAAB(4,IT) = ISM_AB
                END IF
                  IF(IT.EQ.LBAT) THEN
*. Save last indeces for me
                    ISM_C1  = ISM_C
                    ISM_CA1 = ISM_CA
                    ISM_AA1 = ISM_AA 
*
                    INM_AB1 = I_AB 
                    INM_AA1 = I_AA
                    INM_CB1 = I_CB
                    INM_CA1 = I_CA
                    GOTO 1001
                   END IF
*                  ^ End if last element
               END DO
               INI_LOOP = 0
              END DO
              INI_LOOP = 0
             END DO
             INI_LOOP = 0
            END DO
            INI_LOOP = 0
*           ^ End of loop over elements of block
          END DO
          INI_LOOP = 0
*         ^ End of loop over ISM_AA
        END DO
        INI_LOOP = 0
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
 1001 CONTINUE
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*) 
        WRITE(6,*) ' String numbers for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(INM_CAAB,4,LBAT,4,LBAT)
        WRITE(6,*)
        WRITE(6,*) ' Symmetry for ca, cb, aa, ab    '
        WRITE(6,*) ' ==================================== '
        CALL IWRTMA(ISM_CAAB,4,LBAT,4,LBAT)
        WRITE(6,*)
      END IF
*
      RETURN
      END
*
      SUBROUTINE DIM_CNTR_REL(ICONT,NCONT,LDIM,NCNTR)
*
* Find dimension, NCNTR,  
* of contraction operator ICONT having NCONT operators 
*
*. General input
#include "implicit.inc"
#include "mxpdim.inc"
#include "orbinp.inc"
*. Specific input
      INTEGER ICONT(LDIM,4)
*. Local scratch, assuming contraction is atmost 2-e oprator 
      INTEGER IUSED(4)
*
      IZERO = 0
      CALL ISETVC(IUSED,IZERO,NCONT)
      N = 1
      DO IOP = 1, NCONT
      IF(IUSED(IOP).EQ.0) THEN
*. other operators of the same GAS, SPIN AND CA ?
        NOP = 1
        DO JOP = IOP+1,NCONT
          IF(ICONT(IOP,1).EQ.ICONT(JOP,1).AND.
     &       ICONT(IOP,2).EQ.ICONT(JOP,2).AND.
     &       ICONT(IOP,3).EQ.ICONT(JOP,3)     )THEN
               NOP = NOP + 1
               IUSED(JOP) = 1
          END IF
        END DO
        LOB = NOBPT(ICONT(IOP,1))
        IDIV = 1
        DO I = 1, NOP
          N = N*(LOB+1-I)
          IDIV = IDIV * I
        END DO
        N = N/IDIV
      END IF
      END DO
*
      NCNTR = N
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Gasspace  Spin    Cr/An '
        WRITE(6,*) ' ========================'
        DO JCONT = 1, NCONT
          WRITE(6,'(I4,6X,I2,7X,I2)') 
     &    ICONT(JCONT,1),ICONT(JCONT,2),ICONT(JCONT,3)
        END DO
        WRITE(6,*) ' Number of contractions = ', NCNTR
      END IF
*
      RETURN
      END 
*
      SUBROUTINE TI_TO_TOKBN_REL(NSMST,
     &           IOP_SM_C1,IOP_SM_CA1,IOP_SM_AA1,IOP_SM,
     &           IOP_NM_AB1,IOP_NM_AA1,IOP_NM_CB1,IOP_NM_CA1,
     &           NOP_BAT,NOP_CAAB,
     &           KK_SM_C1,KK_SM_CA1,KK_SM_AA1,KK_SM,
     &           KK_NM_AB1,KK_NM_AA1,KK_NM_CB1,KK_NM_CA1,
     &           NK_BAT,NKK_CAAB,
     &           IKJ_CA_MAP,SKJ_CA_MAP,IBKJ_CA,
     &           IKJ_CB_MAP,SKJ_CB_MAP,IBKJ_CB,
     &           IKJ_AA_MAP,SKJ_AA_MAP,IBKJ_AA,
     &           IKJ_AB_MAP,SKJ_AB_MAP,IBKJ_AB,
     &           NISFSM_CA,NISFSM_CB,NISFSM_AA,NISFSM_AB,
     &           TOPK,TI,IB_TI,ISG,LDUM,SIGNI,IOPDAG)
* A set of operators OP_CA, OP_CB, OP_AA, OP_AB is given
* OP_NM_CAAB, OP_SM_CAAB contains a number of such operator quadruplets
*
*
* ISG = 1 :
*   T(IDUM,OP,K) :=  Sign(OP,K)* T(IDUM,I)
* ISG = 2 :
*   T(IDUM,I) := T(IDUM,OP,K) + Sign(OP,K)*T(IDUM,I)
*
* If IOPDAG = 1 IOP_NM, IOP_SM contains info for daggered operator
* instead of operator itself
*
* Version with complete precomputed maps OP * K => I
*
* Jeppe Olsen, October 2000 (TOKB version)
*              June 2001    TOKBN (Turbo version)
*
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "multd2h.inc"
#include "symm.inc"
*. Input
*. Operator
      INTEGER NOP_CAAB(8,4)
*. Number of K-strings per symmetry
      INTEGER NKK_CAAB(8,4)
*. Offset to strings with given symmetry
C     INTEGER IBOP_CA(8),IBOP_CB(8),IBOP_AA(8)
*. Offset to block with given sym in TI
      INTEGER IB_TI(8,8,8)
*. Number of strings in I per symmetry
      INTEGER NISFSM_CA(8),NISFSM_CB(8),NISFSM_AA(8),NISFSM_AB(8)
*. Kstrings
COLD  INTEGER K_SM_CAAB(4,NK_BAT),K_NM_CAAB(4,NK_BAT)
*. Info on OP * K to I mappings
      INTEGER IKJ_CA_MAP(*),IBKJ_CA(8,8)
      INTEGER IKJ_CB_MAP(*),IBKJ_CB(8,8)
      INTEGER IKJ_AA_MAP(*),IBKJ_AA(8,8)
      INTEGER IKJ_AB_MAP(*),IBKJ_AB(8,8)
*
      DIMENSION SKJ_CA_MAP(*),SKJ_CB_MAP(*),SKJ_AA_MAP(*),SKJ_AB_MAP(*)
*. Input or output
C     DIMENSION TOPK(LDUM,NOP_BAT,NK_BAT)
C     DIMENSION TI(LDUM,*)
      DIMENSION TOPK(*),TI(*)
      COMMON/ROLLO/NINNER,NINNER1,NINNER2,NINNER3,NINNER4,NINNER5,
     &             NINNER6,NINNER7,NINNER8,NINNER9,NINNER10
*

      NTEST = 00
      MX_TI_COL = 0
      I_CAAB_TI_MX = 0
*.
C      do i=1,8
C       do j=1,4
C        print*,'NOP_CAAB,i,j',NOP_CAAB(i,j),i,j
C       end do
C      end do
C      print*,'ISG,IOPDAG',ISG,IOPDAG
      IF(ISG.EQ.1) THEN
        ZERO = 0.0D0
        CALL SETVEC(TOPK,ZERO,LDUM*NOP_BAT*NK_BAT)
      END IF
*. Loop over I-operators
      INI_LOOP_OP = 1
      IOP = 0
      DO IOP_SM_C = IOP_SM_C1, NSMST
        IOP_SM_A = IDBGMULT(IOP_SM,INVELM(IOP_SM_C))
        IF(INI_LOOP_OP .EQ. 1 ) THEN
          IOP_SM_CA2 = IOP_SM_CA1
        ELSE
          IOP_SM_CA2 = 1
        END IF
        DO IOP_SM_CA = IOP_SM_CA2, NSMST
          IOP_SM_CB = IDBGMULT(IOP_SM_C,INVELM(IOP_SM_CA))
          IF(INI_LOOP_OP .EQ. 1 ) THEN
            IOP_SM_AA2 = IOP_SM_AA1
          ELSE
            IOP_SM_AA2 = 1
          END IF
          DO IOP_SM_AA = IOP_SM_AA2, NSMST
            IOP_SM_AB =  IDBGMULT(IOP_SM_A,INVELM(IOP_SM_AA))
*. Loop over operators as  matrix (I_CA,I_CB,I_AA,I_AB)
            NOP_STR_CA = NOP_CAAB(IOP_SM_CA,1)
            NOP_STR_CB = NOP_CAAB(IOP_SM_CB,2)
            NOP_STR_AA = NOP_CAAB(IOP_SM_AA,3)
            NOP_STR_AB = NOP_CAAB(IOP_SM_AB,4)
            IF(NOP_STR_CA*NOP_STR_CB*NOP_STR_AA*NOP_STR_AB.NE.0)THEN
            NINNER10 = NINNER10 + 1
*
            IF(INI_LOOP_OP.EQ.1) THEN
             IOP_NM_AB2 = IOP_NM_AB1
            ELSE
             IOP_NM_AB2 = 1
            END IF
            DO IOP_NM_AB = IOP_NM_AB2, NOP_STR_AB
                NINNER9 = NINNER9 + 1
             IF(INI_LOOP_OP.EQ.1) THEN
              IOP_NM_AA2 = IOP_NM_AA1
             ELSE
              IOP_NM_AA2 = 1
             END IF
             DO IOP_NM_AA = IOP_NM_AA2, NOP_STR_AA
                NINNER8 = NINNER8 + 1
              IF(INI_LOOP_OP.EQ.1) THEN
                IOP_NM_CB2 = IOP_NM_CB1
              ELSE
                IOP_NM_CB2 = 1
              END IF
              DO IOP_NM_CB = IOP_NM_CB2, NOP_STR_CB
                NINNER7 = NINNER7 + 1
               IF(INI_LOOP_OP.EQ.1) THEN
                 IOP_NM_CA2 = IOP_NM_CA1
               ELSE
                 IOP_NM_CA2 = 0
               END IF
               DO IOP_NM_CA = IOP_NM_CA2 + 1, NOP_STR_CA
                NINNER6 = NINNER6 + 1
                IOP = IOP + 1
*. A major part of the following could be moved outside
                IF(IOPDAG.EQ.0) THEN
                  IIOP_SM_CA = IOP_SM_CA
                  IIOP_NM_CA = IOP_NM_CA
*
                  IIOP_SM_CB = IOP_SM_CB
                  IIOP_NM_CB = IOP_NM_CB
*
                  IIOP_SM_AA = IOP_SM_AA
                  IIOP_NM_AA = IOP_NM_AA
*
                  IIOP_SM_AB = IOP_SM_AB
                  IIOP_NM_AB = IOP_NM_AB
C      print*,'IIOP_SM_CA,IIOP_NM_CA,IIOP_SM_CB,IIOP_NM_CB',
C     &        IIOP_SM_CA,IIOP_NM_CA,IIOP_SM_CB,IIOP_NM_CB
C      print*,'IIOP_SM_AA,IIOP_NM_AA,IIOP_SM_AB,IIOP_NM_AB',
C     &        IIOP_SM_AA,IIOP_NM_AA,IIOP_SM_AB,IIOP_NM_AB
                ELSE
Cori                  IIOP_SM_CA = IOP_SM_AA
Cori                  IIOP_NM_CA = IOP_NM_AA
*
Cori                  IIOP_SM_CB = IOP_SM_AB
Cori                  IIOP_NM_CB = IOP_NM_AB
*
Cori                  IIOP_SM_AA = IOP_SM_CA
Cori                  IIOP_NM_AA = IOP_NM_CA
*
Cori                  IIOP_SM_AB = IOP_SM_CB
Cori                  IIOP_NM_AB = IOP_NM_CB
C
C comment: Why does IOP_SM_CA come out totally symmetric?
C Is this operator not present? Or is a bug in assigning its symmetry?
C Non-relativistic there is no difference between adjoining and flipping c and a
C with regards to symmetry. 
C The flipping of strings is a total inversion from c to a 
C Here the symmetry should also be adjoint otherwise the symmetry will not be correct
C
                  IIOP_SM_CA = INVELM(IOP_SM_AA)
                  IIOP_NM_CA = INVELM(IOP_NM_AA)
                  IIOP_NM_CA = IOP_NM_AA
*
                  IIOP_SM_CB = INVELM(IOP_SM_AB)
                  IIOP_NM_CB = INVELM(IOP_NM_AB)
                  IIOP_NM_CB = IOP_NM_AB
*
                  IIOP_SM_AA = INVELM(IOP_SM_CA)
                  IIOP_NM_AA = INVELM(IOP_NM_CA)
                  IIOP_NM_AA = IOP_NM_CA
*
                  IIOP_SM_AB = INVELM(IOP_SM_CB)
                  IIOP_NM_AB = INVELM(IOP_NM_CB)
                  IIOP_NM_AB = IOP_NM_CB
C      print*,'IIOP_SM_CA,IIOP_NM_CA,IIOP_SM_CB,IIOP_NM_CB',
C     &        IIOP_SM_CA,IIOP_NM_CA,IIOP_SM_CB,IIOP_NM_CB
C      print*,'IIOP_SM_AA,IIOP_NM_AA,IIOP_SM_AB,IIOP_NM_AB',
C     &        IIOP_SM_AA,IIOP_NM_AA,IIOP_SM_AB,IIOP_NM_AB
                END IF
*
                IF(IOPDAG.EQ.0) THEN
C?                  LOP_AA = NOP_CAAB(IIOP_SM_CA,1)
C?                  LOP_AB = NOP_CAAB(IIOP_SM_CB,2)
C?                  LOP_CA = NOP_CAAB(IIOP_SM_AA,3)
C?                  LOP_CB = NOP_CAAB(IIOP_SM_AB,4)
                  LOP_CA = NOP_CAAB(IIOP_SM_CA,1)
                  LOP_CB = NOP_CAAB(IIOP_SM_CB,2)
                  LOP_AA = NOP_CAAB(IIOP_SM_AA,3)
                  LOP_AB = NOP_CAAB(IIOP_SM_AB,4)
C      print*,'LOP_CA,LOP_CB,LOP_AA,LOP_AB',LOP_CA,LOP_CB,LOP_AA,LOP_AB
                ELSE
C Not sure if this is correct Lasse
                  LOP_AA = NOP_CAAB(IOP_SM_CA,1)
                  LOP_AB = NOP_CAAB(IOP_SM_CB,2)
                  LOP_CA = NOP_CAAB(IOP_SM_AA,3)
                  LOP_CB = NOP_CAAB(IOP_SM_AB,4)
C      print*,'LOP_CA,LOP_CB,LOP_AA,LOP_AB',LOP_CA,LOP_CB,LOP_AA,LOP_AB
C                  LOP_CA = NOP_CAAB(IIOP_SM_CA,3)
C                  LOP_CB = NOP_CAAB(IIOP_SM_CB,4)
C                  LOP_AA = NOP_CAAB(IIOP_SM_AA,1)
C                  LOP_AB = NOP_CAAB(IIOP_SM_AB,2)
C      print*,'LOP_CA,LOP_CB,LOP_AA,LOP_AB',LOP_CA,LOP_CB,LOP_AA,LOP_AB
                END IF
*. Loop over K-strings
      INI_LOOP_KK = 1
      KSTR = 0
      DO KK_SM_C = KK_SM_C1, NSMST
        KK_SM_A = IDBGMULT(KK_SM,INVELM(KK_SM_C))
        IF(INI_LOOP_KK .EQ. 1 ) THEN
          KK_SM_CA2 = KK_SM_CA1
        ELSE
          KK_SM_CA2 = 1
        END IF
        DO KK_SM_CA = KK_SM_CA2, NSMST
          KK_SM_CB = IDBGMULT(KK_SM_C,INVELM(KK_SM_CA))
*
          NKK_STR_CA = NKK_CAAB(KK_SM_CA,1)
          NKK_STR_CB = NKK_CAAB(KK_SM_CB,2)
*
          IBOPK_CA = IBKJ_CA(IIOP_SM_CA, KK_SM_CA)
C
C Need to look at symmetry for I_CA_SM and I_CB_SM
C
          I_CA_SM = IDBGMULT(IIOP_SM_CA,KK_SM_CA)
*
          IBOPK_CB = IBKJ_CB(IIOP_SM_CB, KK_SM_CB)
          I_CB_SM = IDBGMULT(IIOP_SM_CB,KK_SM_CB)
*
          IF(INI_LOOP_KK .EQ. 1 ) THEN
            KK_SM_AA2 = KK_SM_AA1
          ELSE
            KK_SM_AA2 = 1
          END IF
          IF( NKK_STR_CA* NKK_STR_CB.EQ.0) KK_SM_AA2 = NSMST + 1
          DO KK_SM_AA = KK_SM_AA2, NSMST
            KK_SM_AB =  IDBGMULT(KK_SM_A,INVELM(KK_SM_AA))
*
            IBOPK_AA = IBKJ_AA(IIOP_SM_AA, KK_SM_AA)
C            print*,'IIOP_SM_AA,KK_SM_AA,IBOPK_AA',
C     &              IIOP_SM_AA,KK_SM_AA,IBOPK_AA
C
C Need to look at symmetry for I_AA_SM
C
            I_AA_SM = IDBGMULT(IIOP_SM_AA,KK_SM_AA)
*
            IBOPK_AB = IBKJ_AB(IIOP_SM_AB, KK_SM_AB)
*. Offset to block in TI with this symmetry combination
            IBTI_SSSS = IB_TI(I_CA_SM,I_CB_SM,I_AA_SM)
*
*. Loop over operators as  matrix (CA,CB,AA,AB)
            NKK_STR_AA = NKK_CAAB(KK_SM_AA,3)
            NKK_STR_AB = NKK_CAAB(KK_SM_AB,4)
C           IF(NKK_STR_CA*NKK_STR_CB*NKK_STR_AA*NKK_STR_AB.NE.0)THEN

*
            IF(INI_LOOP_KK.EQ.1) THEN
             KK_NM_AB2 = KK_NM_AB1
            ELSE
             KK_NM_AB2 = 1
            END IF
*. Tired of IF/ END IF today
            IF(NKK_STR_AA*NKK_STR_AB.EQ.0) THEN
              KK_NM_AB2 = NKK_STR_AB + 1
            ELSE
              NINNER5 = NINNER5 + 1
            END IF
            DO KK_NM_AB = KK_NM_AB2, NKK_STR_AB
          NINNER4 = NINNER4 + 1
          I_AB = IKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
          S_AB = SKJ_AB_MAP(IBOPK_AB-1+(KK_NM_AB-1)*LOP_AB+IIOP_NM_AB)
C      print*,'I_AB,S_AB,IBOPK_AB,KK_NM_AB,LOP_AB,IIOP_NM_AB',
C     &        I_AB,S_AB,IBOPK_AB,KK_NM_AB,LOP_AB,IIOP_NM_AB
             IF(INI_LOOP_KK.EQ.1) THEN
              KK_NM_AA2 = KK_NM_AA1
             ELSE
              KK_NM_AA2 = 1
             END IF
             DO KK_NM_AA = KK_NM_AA2, NKK_STR_AA
          NINNER3 = NINNER3 + 1
          I_AA = IKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
          S_AA = SKJ_AA_MAP(IBOPK_AA-1+(KK_NM_AA-1)*LOP_AA+IIOP_NM_AA)
C      print*,'I_AA,S_AA,IBOPK_AA,KK_NM_AA,LOP_AA,IIOP_NM_AA',
C     &        I_AA,S_AA,IBOPK_AA,KK_NM_AA,LOP_AA,IIOP_NM_AA
              IF(INI_LOOP_KK.EQ.1) THEN
                KK_NM_CB2 = KK_NM_CB1
              ELSE
                KK_NM_CB2 = 1
              END IF
*
              I_CAAB_TI00 = IBTI_SSSS-1+
     &        (I_AB-1)*
     &        NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)*NISFSM_AA(I_AA_SM)
     &       +(I_AA-1)*NISFSM_CA(I_CA_SM)*NISFSM_CB(I_CB_SM)

*
              IKJ_CB_ADR = IBOPK_CB-1+(KK_NM_CB2-2)*LOP_CB+IIOP_NM_CB
              DO KK_NM_CB = KK_NM_CB2, NKK_STR_CB
               NINNER2 = NINNER2 + 1
               IKJ_CB_ADR =  IKJ_CB_ADR + LOP_CB
*
               I_CB = IKJ_CB_MAP(IKJ_CB_ADR)
               S_CB = SKJ_CB_MAP(IKJ_CB_ADR)
               SIGN123 = S_CB*S_AA*S_AB*SIGNI
*
               IF(INI_LOOP_KK.EQ.1) THEN
                 KK_NM_CA2 = KK_NM_CA1
               ELSE
                 KK_NM_CA2 = 0
               END IF
*.
C               print*,'I_AA,I_AB,I_CB',I_AA,I_AB,I_CB
               IF(I_AA*I_AB*I_CB.EQ.0) THEN
                 KSTR = KSTR + NKK_STR_CA - KK_NM_CA2
                 IF(KSTR.GE.NK_BAT) GOTO 2001
               ELSE
                 I_CAAB_TI0 =  I_CAAB_TI00
     &          +(I_CB-1)*NISFSM_CA(I_CA_SM) - 1
*
                 KJ_ADR = IBOPK_CA-1 + (KK_NM_CA2-1)*LOP_CA+IIOP_NM_CA
                 IF(LDUM.EQ.1) THEN
*
                 IADR_TOPK0 = (KSTR-1)*NOP_BAT + IOP-1
*
                 N_CA_EFF = NKK_STR_CA - KK_NM_CA2
                 IF(KSTR+N_CA_EFF.GT.NK_BAT) N_CA_EFF = NK_BAT - KSTR
                 KSTR = KSTR + N_CA_EFF
*
                 IF(ISG.EQ.1) THEN
C                DO KK_NM_CA = KK_NM_CA2 + 1, NKK_STR_CA
                 DO KK_NM_CA = 1,  N_CA_EFF
C                 KSTR = KSTR + 1
                  IADR_TOPK0 = IADR_TOPK0 + NOP_BAT
                  KJ_ADR = KJ_ADR + LOP_CA
                  I_CA = IKJ_CA_MAP(KJ_ADR)
*
                  IF(I_CA .NE. 0) THEN
*. Adress of I_CB, I_CA, I_AA, I_AB
                    SIGN = SKJ_CA_MAP(KJ_ADR)*SIGN123
                    TOPK(1+IADR_TOPK0) = SIGN*TI(1+I_CAAB_TI0+I_CA)
                  END IF
*                 ^ End if Istrings was nonvanishing
C               IF(KSTR.EQ.NK_BAT) GOTO 2001
               END DO
                IF(KSTR.EQ.NK_BAT) GOTO 2001
*              ^ End of loop over KK_NM_CA
               ELSE
*              ^ ISG switch

C                DO KK_NM_CA = KK_NM_CA2 + 1, NKK_STR_CA
                 DO KK_NM_CA = 1,  N_CA_EFF
C                 KSTR = KSTR + 1
                  IADR_TOPK0 = IADR_TOPK0 + NOP_BAT
C                  print*,'KJ_ADR,LOP_CA',KJ_ADR,LOP_CA
                  KJ_ADR = KJ_ADR + LOP_CA
                  I_CA = IKJ_CA_MAP(KJ_ADR)
C                  print*,'I_CA,KJ_ADR',I_CA,KJ_ADR
*
                  IF(I_CA .NE. 0) THEN
*. Adress of I_CB, I_CA, I_AA, I_AB
                    SIGN = SKJ_CA_MAP(KJ_ADR)*SIGN123
                    TI(1 + I_CAAB_TI0 + I_CA) =
     &              TI(1 + I_CAAB_TI0 + I_CA)
     &             +SIGN*TOPK(1+IADR_TOPK0)
C      print*,'TI,TOPK,SIGN1',
C     &    TI(1 + I_CAAB_TI0 + I_CA),TOPK(1+IADR_TOPK0),SIGN
                  END IF
*                 ^ End if Istrings was nonvanishing
*
C               IF(KSTR.EQ.NK_BAT) GOTO 2001
               END DO
                IF(KSTR.EQ.NK_BAT) GOTO 2001
*              ^ End of loop over KK_NM_CA

               END IF
*              ^ ISF switch
                 ELSE
*                ^ If LDUM .eq .1
                 IADR_TOPK0 = (KSTR-1)*NOP_BAT*LDUM + (IOP-1)*LDUM
                 DO KK_NM_CA = KK_NM_CA2 + 1, NKK_STR_CA
                  NINNER1 = NINNER1 + 1
                  KSTR = KSTR + 1
                  IADR_TOPK0 = IADR_TOPK0 + NOP_BAT*LDUM
                  KJ_ADR = KJ_ADR + LOP_CA
                  I_CA = IKJ_CA_MAP(KJ_ADR)
*
                  IF(I_CA .NE. 0) THEN
                    NINNER = NINNER + LDUM
*. Adress of I_CB, I_CA, I_AA, I_AB
                    SIGN = SKJ_CA_MAP(KJ_ADR)*SIGN123
*
C                   I_CAAB_TI = I_CAAB_TI0 + I_CA
C                   IADR_TI0 = (I_CAAB_TI0 + I_CA)*LDUM
                    IF(ISG.EQ.1) THEN
                      IF(LDUM.EQ.1) THEN
                       IADR_TI0 = I_CAAB_TI0 + I_CA
                       TOPK(1+IADR_TOPK0) = SIGN*TI(1+IADR_TI0)
                      ELSE
                       DO I = 1, LDUM
                         IADR_TI0 = (I_CAAB_TI0 + I_CA)*LDUM
                         TOPK(I+IADR_TOPK0) = SIGN*TI(I+IADR_TI0)
                       END DO
                      END IF
*                     ^ End if LDUM = 1
                    ELSE IF (ISG.EQ.2) THEN
C                     print*,'LDUM',LDUM
                     IF(LDUM.EQ.1) THEN
                      IADR_TI0 = I_CAAB_TI0 + I_CA
                      TI(1+IADR_TI0) =
     &                TI(1+IADR_TI0)+SIGN*TOPK(1+IADR_TOPK0)
C       print*,'TI,TOPK,SIGN2',TI(1+IADR_TI0),TOPK(1+IADR_TOPK0),SIGN
                     ELSE
                      DO I = 1, LDUM
                       IADR_TI0 = (I_CAAB_TI0 + I_CA)*LDUM
                       TI(I+IADR_TI0) =
     &                 TI(I+IADR_TI0)+SIGN*TOPK(I+IADR_TOPK0)
C       print*,'TI,TOPK,SIGN3',TI(1+IADR_TI0),TOPK(1+IADR_TOPK0),SIGN
                      END DO
                    END IF
*                   ^ End if LDUM = 1
                  END IF
*                 ^ End of scatter/gather switch
               END IF
*              ^ End if Istrings was nonvanishing
*
                IF(KSTR.EQ.NK_BAT) GOTO 2001
               END DO
*              ^ End of loop over KK_NM_CA
               END IF
*              ^ End if LDUM = 1 switch

               END IF
*              ^ End if I_AA*I_AB*I_CB .NE. 0
               INI_LOOP_KK = 0
              END DO
*              ^ End of loop over KK_NM_CB
              INI_LOOP_KK = 0
             END DO
             INI_LOOP_KK = 0
            END DO
*           ^ End of loop over elements of block
C           END IF
*           ^ End if number of K-strings was nonvanishing
            INI_LOOP_KK = 0
          END DO
          INI_LOOP_KK = 0
*         ^ End of loop over ISM_AA
        END DO
        INI_LOOP_KK = 0
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
 2001 CONTINUE
* ^ End of loop over K-strings
*
                IF(IOP.EQ.NOP_BAT) GOTO 1001
               END DO
               INI_LOOP_OP = 0
              END DO
              INI_LOOP_OP = 0
             END DO
             INI_LOOP_OP = 0
            END DO
            INI_LOOP_OP = 0
*           ^ End of loop over elements of block
            END IF
*           ^ End if number of I strings was nonvanishing
          INI_LOOP_OP = 0
          END DO
*         ^ End of loop over ISM_AA
        END DO
        INI_LOOP_OP = 0
*        ^ End of loop over ISM_CA
      END DO
*     ^ End of loop over ISM_C
 1001 CONTINUE
*
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' TI_TO_TOKB speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' NOP_BAT, NK_BAT, LDUM = ', NOP_BAT, NK_BAT, LDUM
        WRITE(6,*) ' Largest column number used = ',  I_CAAB_TI_MX
        WRITE(6,*) ' TI as TI(Idum,I)  '
        CALL WRTMAT(TI,LDUM,I_CAAB_TI_MX,LDUM,I_CAAB_TI_MX)
        WRITE(6,*) ' TOPK as TOPK(IdumIop,Kstr) '
        CALL WRTMAT(TOPK,LDUM*NOP_BAT,NK_BAT,LDUM*NOP_BAT,NK_BAT)
      END IF
C      print*,'TI,TOPK',TI(1),TOPK(1)
*
      RETURN
      END
*
      SUBROUTINE STST_TO_ST_MAP_REL(IS1OC,IS2OC,IS12OC,
     &           IBS1S2,IS1S2_TO_S12,XS1S2_TO_S12,
     &           IZ,IZSCR,IS1_STR,IS2_STR,IS12_REO,IUB)
*
* Consider string multiplication S1*S2 => S12 
* Find the mapping between strings             
* Maps are organized as matrices (S1,S2) for each symmetry of S1 and S2
*
* Jeppe Olsen, May 1, 2000
*
* New version reducing number of calls to GET_STR2_TOTSM..
* a rainy day at HNIE, July 2002
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "multd2h.inc"
#include "csm.inc"
#include "lucinp.inc"
#include "orbinp.inc"
#include "symm.inc"
*. Input :  number of operators in each gasspace
C     INTEGER IS1OP(NGAS),IS2OP(NGAS),IS12OP(NGAS)
      INTEGER IS1OC(*),IS2OC(*),IS12OC(*)
*
*. Local scratch
*
      INTEGER IS1GRP(NGAS),IS2GRP(NGAS),IS12GRP(NGAS)
      INTEGER ISTR_OUT(100)
*
*. Scratch space through parameter list
*
      INTEGER IZ(*),IZSCR(*)
      INTEGER IS1_STR(*),IS2_STR(*),IS12_REO(*)
*     ^ Should hold largest list of strings with given CAAB
*. Local scratch 
      INTEGER LS1_STR(8) , IBS1_STR(8)
      INTEGER LS2_STR(8) , IBS2_STR(8)
C     INTEGER LS12_STR(8), IBS12_STR(8)
*
*. Output :
*
*. Offset to mappings for given sym of S1 and S2
      INTEGER IBS1S2(8,8)
*. And the mappings 
      INTEGER IS1S2_TO_S12(*)
      DIMENSION XS1S2_TO_S12(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' STST_TO_ST_MAP entered '
        WRITE(6,*) ' ======================='
        WRITE(6,*)
        WRITE(6,*) ' Occupation of S1, S2 and S12 '
        CALL IWRTMA(IS1OC,1,NGAS,1,NGAS)
        CALL IWRTMA(IS2OC,1,NGAS,1,NGAS)
        CALL IWRTMA(IS12OC,1,NGAS,1,NGAS)
        WRITE(6,*) ' IUB = ', IUB
      END IF
*
*. Occupation to group notation
      CALL OCC_TO_GRP_CC_REL(IS1OC,IS1GRP,1)
      CALL OCC_TO_GRP_CC_REL(IS2OC,IS2GRP,1)
      CALL OCC_TO_GRP_CC_REL(IS12OC,IS12GRP,1)
*. Number of operators 
      NS1OP = IELSUM(IS1OC,NGAS)
      NS2OP = IELSUM(IS2OC,NGAS)
      NS12OP = IELSUM(IS12OC,NGAS)
*. Z array for 12 strings 
      CALL WEIGHT_SPGP(IZ,NGAS,IS12OC,NOBPT,IZSCR,0)
*. Set up I12 Reorder array for all symmetries
      IOFF = 1
C      print*,'NSMST',NSMST
      DO IS12SM = 1, NSMST
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,
     &       IS12GRP,NGAS,IS12SM,NS12OP,NS12STR,
     &       IS1_STR,NOCOB,1,IZ,IS12_REO)
      END DO
*. Set up all I1 strings
      IOFF = 1
      DO IS1SM = 1, NSMST
        IBS1_STR(IS1SM) = IOFF
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,
     &       IS1GRP,NGAS,IS1SM,NS1OP,NS1STR,
     &       IS1_STR(IOFF),NOCOB,0,0,0)              
        LS1_STR(IS1SM) = NS1STR
C        print*,'LS1_STR(IS1SM)',LS1_STR(IS1SM)
        IOFF = IOFF + NS1STR*NS1OP
      END DO
*. And I2 strings
      IOFF = 1
      DO IS2SM = 1, NSMST
        IBS2_STR(IS2SM) = IOFF
        CALL GETSTR2_TOTSM_SPGP_REL(IUB,
     &       IS2GRP,NGAS,IS2SM,NS2OP,NS2STR,
     &       IS2_STR(IOFF),NOCOB,0,0,0)              
        LS2_STR(IS2SM) = NS2STR
C        print*,'LS2_STR(IS2SM)',LS2_STR(IS2SM) 
        IOFF = IOFF + NS2STR*NS2OP
      END DO
*
      IBOFF = 1
      DO IS12SM = 1, NSMST
        DO IS1SM = 1, NSMST
          IS2SM = IDBGMULT(IADJSYM(IS1SM),IS12SM)
          IBS1S2(IS1SM,IS2SM) = IBOFF
          NS1STR = LS1_STR(IS1SM)
          NS2STR = LS2_STR(IS2SM)
          IBS1STR = IBS1_STR(IS1SM)
          IBS2STR = IBS2_STR(IS2SM)
          IF(NTEST.GE.100) 
     &    WRITE(6,*) ' NS1STR, NS2STR = ', NS1STR, NS2STR
          DO IS2 = 1, NS2STR
          DO IS1 = 1, NS1STR
* S1 * S2 => S12 ( all strings are considered to be creation strings
            IS1_OFF =  IBS1STR + (IS1-1)*NS1OP
            IS2_OFF =  IBS2STR + (IS2-1)*NS2OP
            CALL CRAN_STR_REL(IS1_STR(IS1_OFF),IDUM,NS1OP,0,
     &                    IS2_STR(IS2_OFF),NS2OP,ISTR_OUT,
     &                    ISIGN,IZERO_STR)
            IJ = IBOFF -1 +(IS2-1)*NS1STR + IS1
            IF(IZERO_STR.NE.1) THEN
             INUM = ISTRNM(ISTR_OUT,NOCOB,NS12OP,IZ,IS12_REO,1)
             IS1S2_TO_S12(IJ) = INUM
             XS1S2_TO_S12(IJ) = DFLOAT(ISIGN)
            ELSE
             IS1S2_TO_S12(IJ) = 0
             XS1S2_TO_S12(IJ) = 0.0              
            END IF
          END DO
          END DO
*         ^ End of loop over IS1, IS2
          IBOFF = IBOFF + NS1STR*NS2STR
        END DO
      END DO
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' S1S2 => S12 mapping, adress and signs' 
        WRITE(6,*)  ' ====================================='
        DO IS1SM = 1, NSMST
        DO IS2SM = 1, NSMST
          NS1 = LS1_STR(IS1SM)
          NS2 = LS2_STR(IS2SM)
          IB = IBS1S2(IS1SM,IS2SM)
          WRITE(6,*) ' Sym of S1 and S2 : ', IS1SM, IS2SM
          WRITE(6,*)
          CALL IWRTMA(IS1S2_TO_S12(IB),NS1,NS2,NS1,NS2)
          CALL WRTMAT(XS1S2_TO_S12(IB),NS1,NS2,NS1,NS2)
        END DO
        END DO
      END IF
*
      RETURN
      END     
*
      SUBROUTINE CRAN_STR_REL(ICR,IAN,NCR,NAN,ISTR_IN,NEL_IN,
     &                    ISTR_OUT,ISIGN,IZERO_STR)
*
* ISTR_OUT = ISIGN* ICR IAN ISTR_IN
*
* Where ICR is a string of creation operators and IAN is a string 
* of annihilation operators.
*
* Input string is assumed to be given in ascending order, 
* and output string will be delivered with orbitals in 
* ascending order
*
*. Initial version, I hope it is not for mission critical routines 
* (could be speeded up)
*
* Jeppe Olsen, March 2000
*
* Change of phase of annihilations strings, Oct2000
*
#include "implicit.inc"
#include "newccp.inc"
*. Input
      INTEGER ICR(NCR),IAN(NAN)
      INTEGER ISTR_IN(NEL_IN)
*. Output
      INTEGER ISTR_OUT(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CRAN_STR speaking '
        WRITE(6,*) ' =================='
        WRITE(6,*) ' Input string '
        CALL IWRTMA(ISTR_IN,1,NEL_IN,1,NEL_IN) 
        WRITE(6,*) ' Annihilation string '
        CALL IWRTMA(IAN,1,NAN,1,NAN)
        WRITE(6,*) ' Creation string '
        CALL IWRTMA(ICR,1,NCR,1,NCR)
      END IF
*. Make sure that annihilation strings are properly increasing
C?    DO JAN = 2, NAN
C?      IF(IAN(JAN).LE.IAN(JAN-1)) THEN
C?        WRITE(6,*) ' CRAN confused, strange annihilation string :' 
C?        CALL IWRTMA(IAN,1,NAN,1,NAN)
C?        STOP ' CRAN confused, strange annihilation string'
C?      END IF
C?    END DO
*
      NEL_OUT = NEL_IN - NAN + NCR
*
      IZERO_STR = 0
      ISIGN = 1.0D0
      CALL ICOPVE(ISTR_IN,ISTR_OUT,NEL_IN)
*. Annihilate  : IAN(1) IAN(2) .... !STR_IN>
      DO IANNI = 1, NAN
        IFOUND = 0
        DO IEL = 1, NEL_IN-IANNI+1
C?        WRITE(6,*) ' CRAN : IANNI IEL ISTR IAN ',
C?   &    IANNI,IEL,ISTR_OUT(IEL),IAN(NAN-IANNI+1)
          IF(ISTR_OUT(IEL).EQ.IAN(NAN-IANNI+1)) THEN 
            ISIGN = ISIGN*(-1)**(IEL-1)
            IFOUND = 1
            DO JEL = IEL, NEL_IN-IANNI
              ISTR_OUT(JEL) = ISTR_OUT(JEL+1)
             END DO
          END IF
        END DO
        IF(IFOUND.EQ.0) THEN
*. orbital to be annihilated not found, output string is zero
          IZERO_STR = 1
          GOTO 1001 
        END  IF
      END DO
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Annihilated string '
        CALL IWRTMA(ISTR_OUT,NEL_IN-NAN,1,NEL_IN-NAN,1)
      END IF
*. Creation part 
      DO ICREA = 1, NCR
*. Place to insert orbital 
        ICR_ORB = ICR(NCR-ICREA+1)
        IPLACE = 1
        DO IEL = 1, NEL_IN-NAN + ICREA-1
          IF(ISTR_OUT(IEL).EQ.ICR_ORB) THEN
*. Electron is already around, zero
            IZERO_STR = 1
            GOTO 1001 
          END IF
*
          IF(IEL.LT.NEL_IN-NAN + ICREA-1) THEN
            IF(ISTR_OUT(IEL).LT.ICR_ORB.AND.
     &         ISTR_OUT(IEL+1).GT.ICR_ORB   ) THEN
               IPLACE = IEL+1
            END IF
          ELSE IF(IEL.EQ.NEL_IN-NAN+ICREA-1) THEN
            IF(ISTR_OUT(IEL).LT. ICR_ORB   ) THEN
              IPLACE = IEL + 1
            END IF
          END IF
        END DO
        ISIGN  = ISIGN*(-1)**(IPLACE-1)
        DO IEL = NEL_IN-NAN+ICREA,IPLACE + 1, -1
          ISTR_OUT(IEL) = ISTR_OUT(IEL-1)
        END DO
        ISTR_OUT(IPLACE) = ICR_ORB
      END DO
*
 1001 CONTINUE
*. A bit on the sign : In LUCIA the order of the annihilation is actually 
* descending !, find permutation sign
      IF(I_USE_NEWCCP.EQ.0.AND.NAN.GT.1) THEN
        NPERM = NAN*(NAN-1)/2
        ISIGN = ISIGN*(-1)**NPERM
      END IF
*
      IF(NTEST.GE.100) THEN 
        IF(IZERO_STR.EQ.0) THEN
          WRITE(6,*) ' Output string '
          CALL IWRTMA(ISTR_OUT,1,NEL_OUT,1,NEL_OUT)
          WRITE(6,*) ' ISIGN = ', ISIGN
        ELSE 
          WRITE(6,*) ' Vanishing string '
        END IF
      END IF
*
      RETURN
      END
*
      SUBROUTINE OCC_TO_GRP_CC_REL(IOCC,IGRP,IWAY)
*. Translate between occupation and group labels of supergroup
*
* IWAY = 1 => OCC to group
*      = 2 => Group to occ
*
* Number of strings in supergroup with given symmetry and number of 
* electrons in each space
*
* Jeppe Olsen, March 1999
*
#include "implicit.inc"
*. General input
#include "mxpdim.inc"
#include "cgas.inc" 
#include "gasstr.inc"
*. Local scratch
C     INTEGER JGRP(MXPNGAS)
*. Specific input/output
      INTEGER IOCC(NGAS),IGRP(NGAS)
*
      NTEST = 00
*
      IF(IWAY.EQ.1) THEN 
*. Occupation => Group number 
        DO IOBTP = 1, NGAS
          JJGRP = 0
          DO KGRP = IBGPSTR(IOBTP), IBGPSTR(IOBTP) + NGPSTR(IOBTP)-1
            IF(NELFGP(KGRP).EQ.IOCC(IOBTP)) JJGRP = KGRP
          END DO
          IGRP(IOBTP) = JJGRP
*
          IF(JJGRP.EQ.0) THEN
           WRITE(6,*) ' Group not included in list '
           WRITE(6,*) ' GAS space with problem : ', IOBTP
           WRITE(6,*) ' Input occupations : '
           CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
           print*,'KGRP range ',IBGPSTR(IOBTP), IBGPSTR(IOBTP) 
     &                          + NGPSTR(IOBTP)-1
           STOP       ' Group not included in list '
          END IF
        END DO
*
      ELSE
*. Group => Occupation 
        DO IOBTP = 1, NGAS
          IOCC(IOBTP) = NELFGP(IGRP(IOBTP))
        END DO
      END IF
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' Occupation and corresponding group array '
        CALL IWRTMA(IOCC,1,NGAS,1,NGAS)
        CALL IWRTMA(IGRP,1,NGAS,1,NGAS) 
      END IF
*
      RETURN
      END 
*
      SUBROUTINE GET_OPINT4_REL(OPSCR,
     &IOD4X,LD4,INM_CAAB_D4,ISM_CAAB_D4,ISTR_D4,IBO4DX,
     &IOD3X,LD3,INM_CAAB_D3,ISM_CAAB_D3,ISTR_D3,IBO3DX,
     &IOD2X,LD2,INM_CAAB_D2,ISM_CAAB_D2,ISTR_D2,IBO2DX,
     &IOD1X,LD1,INM_CAAB_D1,ISM_CAAB_D1,ISTR_D1,IBO1DX,
     &IOEX ,LEX,INM_CAAB_EX,ISM_CAAB_EX,ISTR_EX,IBOEX,
     &IEXD1234_INDX)
*
* Fetch batch of operator integrals ordered as O(D1,D2,D3,D4,EX) 
*
* Connected with standard integral input
*
* Jeppe Olsen, May of 2000
*
#include "implicit.inc"
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "ctcc.inc"
#include "cgas.inc"
#include "csm.inc"
#include "orbinp.inc"
C
#include "wrkspc.inc"
#include "glbbas.inc"
*. Input 
*.. Occupation of gas spaces         
      INTEGER IOD4X(NGAS,4),IOD3X(NGAS,4)
      INTEGER IOD2X(NGAS,4),IOD1X(NGAS,4),IOEX(NGAS,4)
*. Index in EXD1D2D3D4 to original index in H
      INTEGER IEXD1234_INDX(4)
*.. The sym and number of the included strings 
      INTEGER ISM_CAAB_EX(4,*),INM_CAAB_EX(4,*)
      INTEGER ISM_CAAB_D1(4,*),INM_CAAB_D1(4,*)
      INTEGER ISM_CAAB_D2(4,*),INM_CAAB_D2(4,*)
      INTEGER ISM_CAAB_D3(4,*),INM_CAAB_D3(4,*)
      INTEGER ISM_CAAB_D4(4,*),INM_CAAB_D4(4,*)
*.. The actual occupation of the various strings
      INTEGER ISTR_D1(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D2(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D3(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_D4(MX_ST_TSOSO_BLK_MX*NSMST,4)
      INTEGER ISTR_EX(MX_ST_TSOSO_BLK_MX*NSMST,4)
*.. Start of strings with given symmetry
      INTEGER IBO1DX(8,4),IBO2DX(8,4),IBO3DX(8,4), IBO4DX(8,4)
      INTEGER IBOEX(8,4)
* 
*. Local scratch
*
C     INTEGER IOCC(4),ICREA(4),IANNI(4)
      character*4 TYPE
      INTEGER JD1STR(4),JD2STR(4),JD3STR(4),JD4STR(4),JEXSTR(4)
      INTEGER IJKL_EDD(4), IJKL_ORIG(4)
      INTEGER NOP_D1_CAAB(4),NOP_D2_CAAB(4),NOP_EX_CAAB(4)
      INTEGER NOP_D3_CAAB(4),NOP_D4_CAAB(4)
C
      INTEGER IKPERM(4),IKFAC
*. Output
      DIMENSION OPSCR(*)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) '  GET_OPINT speaking : '
        WRITE(6,*) ' ===================== '
        WRITE(6,*)
        WRITE(6,*) ' IOEX, IOD4X, IOD3X, IOD2X, IOD1X  : '
        CALL WRT_SPOX_TP_CC_KRCC(IOEX,1)
        CALL WRT_SPOX_TP_CC_KRCC(IOD4X,1)
        CALL WRT_SPOX_TP_CC_KRCC(IOD3X,1)
        CALL WRT_SPOX_TP_CC_KRCC(IOD2X,1)
        CALL WRT_SPOX_TP_CC_KRCC(IOD1X,1)
        WRITE(6,*)
      END IF
*
*. Number of C/A A/B indeces 
*. D1
      NOP_D1_CA = IELSUM(IOD1X(1,1),NGAS)
      NOP_D1_CB = IELSUM(IOD1X(1,2),NGAS)
      NOP_D1_AA = IELSUM(IOD1X(1,3),NGAS)
      NOP_D1_AB = IELSUM(IOD1X(1,4),NGAS)
      NOP_D1 = NOP_D1_CA+NOP_D1_CB+NOP_D1_AA+NOP_D1_AB
      NOP_D1_CAAB(1) = NOP_D1_CA
      NOP_D1_CAAB(2) = NOP_D1_CB
      NOP_D1_CAAB(3) = NOP_D1_AA
      NOP_D1_CAAB(4) = NOP_D1_AB
*. D2
      NOP_D2_CA = IELSUM(IOD2X(1,1),NGAS)
      NOP_D2_CB = IELSUM(IOD2X(1,2),NGAS)
      NOP_D2_AA = IELSUM(IOD2X(1,3),NGAS)
      NOP_D2_AB = IELSUM(IOD2X(1,4),NGAS)
      NOP_D2 = NOP_D2_CA+NOP_D2_CB+NOP_D2_AA+NOP_D2_AB
      NOP_D2_CAAB(1) = NOP_D2_CA
      NOP_D2_CAAB(2) = NOP_D2_CB
      NOP_D2_CAAB(3) = NOP_D2_AA
      NOP_D2_CAAB(4) = NOP_D2_AB
*. D3
      NOP_D3_CA = IELSUM(IOD3X(1,1),NGAS)
      NOP_D3_CB = IELSUM(IOD3X(1,2),NGAS)
      NOP_D3_AA = IELSUM(IOD3X(1,3),NGAS)
      NOP_D3_AB = IELSUM(IOD3X(1,4),NGAS)
      NOP_D3 = NOP_D3_CA+NOP_D3_CB+NOP_D3_AA+NOP_D3_AB
      NOP_D3_CAAB(1) = NOP_D3_CA
      NOP_D3_CAAB(2) = NOP_D3_CB
      NOP_D3_CAAB(3) = NOP_D3_AA
      NOP_D3_CAAB(4) = NOP_D3_AB
*. D4
      NOP_D4_CA = IELSUM(IOD4X(1,1),NGAS)
      NOP_D4_CB = IELSUM(IOD4X(1,2),NGAS)
      NOP_D4_AA = IELSUM(IOD4X(1,3),NGAS)
      NOP_D4_AB = IELSUM(IOD4X(1,4),NGAS)
      NOP_D4 = NOP_D4_CA+NOP_D4_CB+NOP_D4_AA+NOP_D4_AB
*
      NOP_D4_CAAB(1) = NOP_D4_CA
      NOP_D4_CAAB(2) = NOP_D4_CB
      NOP_D4_CAAB(3) = NOP_D4_AA
      NOP_D4_CAAB(4) = NOP_D4_AB
*
      NOP_EX_CA = IELSUM(IOEX(1,1),NGAS)
      NOP_EX_CB = IELSUM(IOEX(1,2),NGAS)
      NOP_EX_AA = IELSUM(IOEX(1,3),NGAS)
      NOP_EX_AB = IELSUM(IOEX(1,4),NGAS)
      NOP_EX = NOP_EX_CA+NOP_EX_CB+NOP_EX_AA+NOP_EX_AB
      NOP_EX_CAAB(1) = NOP_EX_CA
      NOP_EX_CAAB(2) = NOP_EX_CB
      NOP_EX_CAAB(3) = NOP_EX_AA
      NOP_EX_CAAB(4) = NOP_EX_AB
*
      NCREA_ALPHA  = NOP_D1_CA + NOP_D2_CA + NOP_EX_CA
     &             + NOP_D3_CA + NOP_D4_CA
      NCREA_BETA   = NOP_D1_CB + NOP_D2_CB + NOP_EX_CB
     &             + NOP_D3_CB + NOP_D4_CB
      NANNI_ALPHA  = NOP_D1_AA + NOP_D2_AA + NOP_EX_AA
     &             + NOP_D3_AA + NOP_D4_AA
      NANNI_BETA   = NOP_D1_AB + NOP_D2_AB + NOP_EX_AB
     &             + NOP_D3_AB + NOP_D4_AB
*
C?    WRITE(6,*) ' NANNI_ALPHA, NANNI_BETA = ', 
C?   &             NANNI_ALPHA, NANNI_BETA
*
      NCREA = NCREA_ALPHA + NCREA_BETA
      NANNI = NANNI_ALPHA + NANNI_BETA
*
      NALPHA = NANNI_ALPHA + NCREA_ALPHA
      NBETA  = NANNI_BETA  + NCREA_BETA
*
      NOP = NCREA + NANNI
C      print*,'NALPHA',NALPHA
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Original index '
        CALL IWRTMA(IEXD1234_INDX,1,NOP,1,NOP)
      END IF
************************************************************
      INT = 0
      DO JEX = 1, LEX
C?    WRITE(6,*) ' EX string : '
      CALL GET_CCSTR_FROM_LIST 
     &     (JEXSTR,INM_CAAB_EX(1,JEX),ISM_CAAB_EX(1,JEX),
     &      NOP_EX_CAAB,ISTR_EX(1,1),ISTR_EX(1,2),ISTR_EX(1,3),
     &      ISTR_EX(1,4),IBOEX(1,1),IBOEX(1,2),IBOEX(1,3),
     &      IBOEX(1,4) )
      DO JD4 = 1, LD4
      CALL GET_CCSTR_FROM_LIST 
     &     (JD4STR,INM_CAAB_D4(1,JD4),ISM_CAAB_D4(1,JD4),
     &      NOP_D4_CAAB,ISTR_D4(1,1),ISTR_D4(1,2),ISTR_D4(1,3),
     &      ISTR_D4(1,4),IBO4DX(1,1),IBO4DX(1,2),IBO4DX(1,3),
     &      IBO4DX(1,4) )
      DO JD3 = 1, LD3
      CALL GET_CCSTR_FROM_LIST 
     &     (JD3STR,INM_CAAB_D3(1,JD3),ISM_CAAB_D3(1,JD3),
     &      NOP_D3_CAAB,ISTR_D3(1,1),ISTR_D3(1,2),ISTR_D3(1,3),
     &      ISTR_D3(1,4),IBO3DX(1,1),IBO3DX(1,2),IBO3DX(1,3),
     &      IBO3DX(1,4) )
      DO JD2 = 1, LD2
C?    WRITE(6,*) ' D2 string : '
      CALL GET_CCSTR_FROM_LIST 
     &     (JD2STR,INM_CAAB_D2(1,JD2),ISM_CAAB_D2(1,JD2),
     &      NOP_D2_CAAB,ISTR_D2(1,1),ISTR_D2(1,2),ISTR_D2(1,3),
     &      ISTR_D2(1,4),IBO2DX(1,1),IBO2DX(1,2),IBO2DX(1,3),
     &      IBO2DX(1,4) )
      DO JD1 = 1, LD1
*. Occupation of J1 string in JD1STR
C?    WRITE(6,*) ' D1 string : '
      CALL GET_CCSTR_FROM_LIST 
     &     (JD1STR,INM_CAAB_D1(1,JD1),ISM_CAAB_D1(1,JD1),
     &      NOP_D1_CAAB,ISTR_D1(1,1),ISTR_D1(1,2),ISTR_D1(1,3),
     &      ISTR_D1(1,4),IBO1DX(1,1),IBO1DX(1,2),IBO1DX(1,3),
     &      IBO1DX(1,4) )
*. Indeces of operator EXD1D2 
      DO IOP_EX = 1, NOP_EX
        IJKL_EDD(IOP_EX) = JEXSTR(IOP_EX)
      END DO
      DO IOP_D1 = 1, NOP_D1
        IJKL_EDD(NOP_EX+IOP_D1) = JD1STR(IOP_D1)
      END DO
      DO IOP_D2 = 1, NOP_D2
        IJKL_EDD(NOP_EX+NOP_D1+IOP_D2) = JD2STR(IOP_D2)
      END DO
      DO IOP_D3 = 1, NOP_D3
        IOP_EFF = NOP_EX+NOP_D1+NOP_D2 + IOP_D3
        IJKL_EDD(IOP_EFF) = JD3STR(IOP_D3)
      END DO
      DO IOP_D4 = 1, NOP_D4
        IOP_EFF = NOP_EX+NOP_D1+NOP_D2 + NOP_D3 + IOP_D4
        IJKL_EDD(IOP_EFF) = JD4STR(IOP_D4)
      END DO
*. Original order
      DO I = 1,4
        IJKL_ORIG(I) = 0
      END DO
      DO IOP = 1, NOP   
        IJKL_ORIG(IEXD1234_INDX(IOP)) = IJKL_EDD(IOP)
      END DO
*
      IF(NTEST.GE.100) THEN
      WRITE(6,*) ' IJKL_EDD and IJKL_ORIG '
        CALL IWRTMA(IJKL_EDD,1,NOP,1,NOP)
        CALL IWRTMA(IJKL_ORIG,1,NOP,1,NOP)    
      END IF
*
      INT = INT + 1
C
C fast version only
C        print*,'in fast version only'
        INOINT = 0
        IF(NOP.EQ.2.AND.(NALPHA.EQ.2.OR.NBETA.EQ.2)) THEN
          IOFF = NTOOB*(IJKL_ORIG(1) - 1) + IJKL_ORIG(2) - 1
          OPSCR(INT) = WORK(KFI + IOFF)
        END IF
        IF(NOP.EQ.4) THEN
        IF(NALPHA.EQ.4.OR.NALPHA.EQ.0.OR.ISPINORBIT.EQ.1) THEN
          IF(IJKL_ORIG(1).EQ.IJKL_ORIG(2).OR.
     &       IJKL_ORIG(3).EQ.IJKL_ORIG(4)) THEN
            INOINT = 1
          END IF
        END IF
        IF(NOP.EQ.4.AND.INOINT.EQ.0) THEN
        IJKL_ORIG(1) = ABS(IJKL_ORIG(1))
        IJKL_ORIG(2) = ABS(IJKL_ORIG(2))
        IJKL_ORIG(3) = ABS(IJKL_ORIG(3))
        IJKL_ORIG(4) = ABS(IJKL_ORIG(4))
        CALL GET_INT2(XOUT,WORK(KT_CC),IJKL_ORIG(1),IJKL_ORIG(4),
     &                                 IJKL_ORIG(2),IJKL_ORIG(3),
     &                                 NALPHA)
        ELSE
          XOUT = 0.0D0
        END IF
        OPSCR(INT) = XOUT
        END IF
      NTEST = 00
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' INT, OPSCR(INT) = ', INT, OPSCR(INT)
        print*,'NALPHA,NBETA',NALPHA,NBETA
      END IF
      NTEST = 0
*.    ^ End of switch between different integral types
      END DO
      END DO
      END DO
      END DO
      END DO
*
      NTEST = 00
      IF(NTEST.GE.100) THEN 
         WRITE(6,*)
         WRITE(6,*) ' Output matrix from GET_OPINT as X(D4D3D2EX,D1)' 
         WRITE(6,*) ' ============================================= '
         WRITE(6,*)
         CALL WRTMAT(OPSCR,LD4*LD3*LD2*LEX,LD1,LD4*LD3*LD2*LEX,LD1)
      END IF
*
      RETURN
      END 
*
      SUBROUTINE GET_CCSTR_FROM_LIST(ISTRING,INUM_CAAB,ISYM_CAAB,
     &           NOP_CAAB,LIST_CA,LIST_CB,LIST_AA,LIST_AB,
     &           IB_CA,IB_CB,IB_AA,IB_AB)
* Obtain CC string defined by INUM_CAAB, ISYM_CAAB
*
* Jeppe Olsen, May 2000 is running out
*
#include "implicit.inc"
*. Input
      INTEGER INUM_CAAB(4),ISYM_CAAB(4)
      INTEGER NOP_CAAB(4)
      INTEGER LIST_CA(*),LIST_CB(*),LIST_AA(*),LIST_AB(*)
      INTEGER IB_CA(*),IB_CB(*),IB_AA(*),IB_AB(*)
*. Output
      INTEGER ISTRING(*)
*. CA part
C      print*,'printout'
C      do I=1,4
C        print*,'INUM_CAAB,ISYM_CAAB,NOP_CAAB,I',
C     &          INUM_CAAB(I),ISYM_CAAB(I),NOP_CAAB(I),I
C      end do
C      print*,'entering GET_CCSTR_FROM_LIST'
      IOFF = 1 
      IF(NOP_CAAB(1).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(1),ISYM_CAAB(1),
     &                       NOP_CAAB(1),LIST_CA,IB_CA)
      IOFF = IOFF + NOP_CAAB(1)
      END IF
*. CB 
      IF(NOP_CAAB(2).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(2),ISYM_CAAB(2),
     &                       NOP_CAAB(2),LIST_CB,IB_CB)
      IOFF = IOFF + NOP_CAAB(2)
      END IF
*. AA
      IF(NOP_CAAB(3).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(3),ISYM_CAAB(3),
     &                       NOP_CAAB(3),LIST_AA,IB_AA)
      IOFF = IOFF + NOP_CAAB(3)
      END IF
*. AA
      IF(NOP_CAAB(4).NE.0) THEN
      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(4),ISYM_CAAB(4),
     &                       NOP_CAAB(4),LIST_AB,IB_AB)
      IOFF = IOFF + NOP_CAAB(4)
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        NOP_TOT = IOFF - 1
        WRITE(6,*) ' CAAB string from GET_CCSTR ..'
        CALL IWRTMA(ISTRING,1,NOP_TOT,1,NOP_TOT)
      END IF
*
      RETURN
      END
*
      SUBROUTINE CONJ_CAAB(ICAAB_IN,ICAAB_OUT,NGAS,SIGN)
*
* Conjugate ICAAB_IN to obtain ICAAB_OUT
* 
*. Jeppe Olsen, Oct 2000
#include "implicit.inc"
*.Input 
      INTEGER ICAAB_IN(NGAS,4)
*. Output
      INTEGER ICAAB_OUT(NGAS,4)
*. AA_out is obtain by conjugating CA_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,1), ICAAB_OUT(1,3),NGAS)
      CALL ICOPVE(ICAAB_IN(1,1), ICAAB_OUT(1,3),NGAS)
*. AB_out is obtain by conjugating CB_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,2), ICAAB_OUT(1,4),NGAS)
      CALL ICOPVE(ICAAB_IN(1,2), ICAAB_OUT(1,4),NGAS)
*. CA_out is obtain by conjugating AA_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,3), ICAAB_OUT(1,1),NGAS)
      CALL ICOPVE(ICAAB_IN(1,3), ICAAB_OUT(1,1),NGAS)
*. CB_out is obtain by conjugating AB_IN
CE    CALL REV_ORD_IARR(ICAAB_IN(1,4), ICAAB_OUT(1,2),NGAS)
      CALL ICOPVE(ICAAB_IN(1,4), ICAAB_OUT(1,2),NGAS)
*. Directly the conjugated operator is 
*  O(ab){\dag}O(aa){\dag}O(cb){\dag}O(ca){\dag}
*. Sign required to change  
*  O(aa){\dag}O(ab){\dag}O(ca){\dag}O(cb){\dag}
      NCA = IELSUM(ICAAB_IN(1,1),NGAS)
      NCB = IELSUM(ICAAB_IN(1,2),NGAS)
      NAA = IELSUM(ICAAB_IN(1,3),NGAS)
      NAB = IELSUM(ICAAB_IN(1,4),NGAS)
*
      NPERM = NCA*NCB + NAA*NAB
*. Sign required to bring the individual strings into ascending order
      NPERM = NPERM + 
     &        NCA*(NCA-1)/2+NCB*(NCB-1)/2+NAA*(NAA-1)/2+NAB*(NAB-1)/2
*
      IF(MOD(NPERM,2).EQ.1) THEN
        SIGN = -1.0D0
      ELSE
        SIGN = 1.0D0
      END IF
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' CAAB operator and conjugated CAAB operator : '
        CALL WRT_SPOX_TP_CC_KRCC(ICAAB_IN,1)
        CALL WRT_SPOX_TP_CC_KRCC(ICAAB_OUT,1)
        WRITE(6,*) ' Sign = ', sign
      END IF
*
      RETURN
      END
*
      SUBROUTINE REF_OP_CC(IOPGAS,IOP,NOP,NGAS,IWAY)
* An operatorstring may be specifed as
*
* IOPGAS : Number of operators per GASspace
* IOP    : GASpace of each operator
*
* Transform between these two form
*
* Iway = 1 : IOPGAS => IOP
* Iway = 2 : IOP    => IOPGAS
*
* Jeppe Olsen, Summer of 99
*
#include "implicit.inc"
*. Input/Output
      INTEGER IOPGAS(NGAS),IOP(NOP)
*
      IF(IWAY.EQ.1) THEN
        JOP = 0
        DO JGAS = 1, NGAS
          LJGAS = IOPGAS(JGAS)
          DO JJOP = 1, LJGAS
            JOP = JOP + 1
            IOP(JOP) = JGAS
          END DO
        END DO
        NOP = JOP
      ELSE
        DO JGAS = 0, NGAS
          JOP = 0
          DO JJOP = 1, NOP
            IF(IOP(JJOP).EQ.JGAS) JOP = JOP +1 
          END DO
          IOPGAS(JGAS) = JOP
        END DO
      END IF
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        IF(IWAY.EQ.1) THEN
         WRITE(6,*) ' IOPGAS => IOP '
        ELSE
         WRITE(6,*) ' IOP => IOPGAS '
        END IF
        WRITE(6,*) ' IOPGAS and IOP '
        CALL IWRTMA(IOPGAS,1,NGAS,1,NGAS)
        CALL IWRTMA(IOP,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE GET_STR_FROM_LIST(ISTRING,INUM,ISYM,NOP,LIST,IB)
C
C      CALL GET_STR_FROM_LIST(ISTRING(IOFF),INUM_CAAB(1),ISYM_CAAB(1),
C     &                       NOP_CAAB(1),LIST_CA,IB_CA)
*
* Obtain string INUM of sym ISYM from complete list of strings
*
* Jeppe Olsen, Still May 2000
*
#include "implicit.inc"
*.Input
      INTEGER LIST(NOP,*),IB(*)
*. Output 
      INTEGER ISTRING(NOP)
*
      INUM_ABS = IB(ISYM)-1+INUM
C      print*,'INUM_ABS',INUM_ABS
C      print*,'IB,ISYM,,INUM',IB(ISYM),ISYM,INUM
C      print*,'LIST(1,INUM_ABS),NOP',LIST(1,INUM_ABS),NOP
      CALL ICOPVE(LIST(1,INUM_ABS),ISTRING,NOP)
*
C?    IF(NOP.GT.0) WRITE(6,*) ' GET_STR... LIST(1,1) = ', LIST(1,1)
C?    IF(NOP.GT.0) WRITE(6,*) ' INUM_ABS = ', INUM_ABS
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' String number and sym = ', INUM,ISYM
        WRITE(6,*) ' Corresponding string : '
        CALL IWRTMA(ISTRING,1,NOP,1,NOP)
      END IF
*
      RETURN
      END
*
      SUBROUTINE WRT_CNTR3(ICONT,NCONT,LDIM)
*
* Write contraction operator ICONT in the form used in VNEWCCV
*
#include "implicit.inc"
      INTEGER ICONT(LDIM,3)
*
      WRITE(6,*)
     &' Information about operators to be contracted'
      WRITE(6,*) ' Gasspace  Cr/An   Spin  '
      WRITE(6,*) ' ========================'
      DO JCONT = 1, NCONT
        WRITE(6,'(I4,6X,I2,7X,I2,7X,I2)')
     &  ICONT(JCONT,1),ICONT(JCONT,2),ICONT(JCONT,3)
      END DO
*
      RETURN
      END
      SUBROUTINE WRT_3SM_ARRAY(ISYM, NSYM, MXPSYM)
*
* Write array containing three symmetry labels 
*
#include "implicit.inc"
*
      INTEGER ISYM(MXPSYM,MXPSYM,MXPSYM)
*
      DO I3 = 1, NSYM
        WRITE(6,*) ' Value of third index = ', I3
        CALL IWRTMA(ISYM(1,1,I3),NSYM,NSYM,MXPSYM,MXPSYM)
      END DO
*
      RETURN
      END
*
*
*
*
*
*
