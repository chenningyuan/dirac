*
      SUBROUTINE RSBB2A(ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,NSCOL,
     &                  NGAS,ISOC,ICOC,
     &                  SB,CB,
     &                  ADSXA,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,ITSOB,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  RIKSX,RKJSX,MXSXST,MXSXBL,IMOC,SCLFAC,NTESTG,
     &                  NSEL2E,ISEL2E,IRI,IBU,SCR)
*
* two electron excitations on column strings
*
* Updated for double group calculation (IRI_INT, IBU added)
*
* =====
* Input
* =====
*
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NSCOL : Number of columns in S block
* ISEL1(3) : Number of electrons in RAS1(3) for S block
* ICEL1(3) : Number of electrons in RAS1(3) for C block
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* DXSTST : Sym of sx,!st> => sym of sx !st>
* STSTDX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* SXDXSX : Symmetry of SX1,SX1*SX2 => symmetry of SX2
* NTSOB  : Number of orbitals per type and symmetry
* IBTSOB : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* NSMOB,NSMST,NSMSX,NSMDX : Number of symmetries of orbitals, strings,
*                           single excitations, double excitations
* MAXI   : Largest number of 'spectator strings' treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
*
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
* XINT : Space for two electron integrals
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS),DXSTST(NSMDX,NSMST)
      INTEGER STSTDX(NSMST,NSMST)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
*
      INTEGER ISEL2E(*)
C     INTEGER NTSOB(3,*),IBTSOB(3,*),ITSOB(*)
*.Input
      DIMENSION CB(*)
      INTEGER ISOC(NGAS),ICOC(NGAS)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),XINT(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
      DIMENSION RIKSX(MXSXBL,4),RKJSX(MXSXBL,4)
      dimension SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*
*.Local arrays
      DIMENSION ITP(256),JTP(256),KTP(256),LTP(256)
*
      INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
CTF   PARAMETER(MXTSOB = 40)
*
#include "symm.inc"
*
*   All important symmetry variables not on pointers
*
C     character*6 DBGSYM
C     integer   NIRR_PN,DDIM,NIRR_DG,DOUGRP,
C    &          SYDI,SYMDIM,ITOTDM,
C    &          IDIMGRP,NSDGIRR,NSPGIRR,
C    &          N_SPINR_EXP,NGASPI,NGSOBPI,IMAXDIM,
C    &          ITRIA,ITRIADB,ITRIASUM,ITRIADBSUM
C     dimension IDBGMULT(8,8),IPGMULT(4,4),IMOSF_SP(8,2),
C    &          ISPSF_MO(8,2),IADJSYM(8),INVELM(8)
C     dimension SYDI(8),SYMDIM(4),
C    &          NGSOBPI(8,16),
C    &          ITRIA(8),ITRIADB(4),ITOTDM(8)
C     common /DGSYMM/NIRR_PN,DDIM,NIRR_DG,DOUGRP,
C    &          SYDI,SYMDIM,ITOTDM,
C    &          IDIMGRP,NSDGIRR,NSPGIRR,
C    &          N_SPINR_EXP,NGASPI,NGSOBPI,IMAXDIM,
C    &          ITRIA,ITRIADB,ITRIASUM,ITRIADBSUM,
C    &          IDBGMULT,IPGMULT,IMOSF_SP,ISPSF_MO,
C    &          IADJSYM,INVELM

*
* DDIM	 	: double dimension (of SYMDIM(x))
* NIRR_DG	: number of double group irreps
* DOUGRP	: double group
* NIRR_PN	: number of spatial irreps
* SYDI		: symmetry dimensions of spatial irreps
* SYMDIM	: two SYDI's added (depending on double group)
* IDIMGRP	: dimension group for int trafo (1 SYMDIM)
* NSDGIRR	: number of sub double group irreps
* NSPGIRR	: number of sub point group irreps
* N_SPINR_EXP   : number of spinor expansions
* NGASPI	: number of GAS spaces (in picasso)
* NGSOBPI	: matrix with number orbitals in GAS and symmetries
* IMAXDIM	: largest dimension of spatial irreps (number from ordering
*                    in Koster, Dimmock, Wheeler, Statz)
* ITRIA         : triangle lengths of point group blocks
* ITRIASUM      : sum of these
* ITRIADB       : triangle lengths of two connected blocks in double group
* ITOTDM        : total orbital space dimensions per spatial irrep
* IDBGMULT      : Double Group multiplication table
* IPGMULT       : Point group mult. table
* IMOSF_SP      : Translation table from MO symmetry to spinor symmetry
*                 IMOSF_SP(IOBSM,IAB) Gives the spinor symmetry
*                 of an orbital of symmetry IOBSM times an alpha (IAB=1)
*                 or and beta spin-function (IAB = 2 )
* ISPSF_MO      : Translation table from spinor symmetry to MO symmetry
*                 ISPSF(ISPSM,IAB) gives the symmetry of an MO that
*                 combined with an alpha(IAB=1) or beta(IAB=2)
*                 spin function give a spinor with symmetry ISPSM
*
* IADJSYM       : IADJSYM(ISPSYM) is the adjoint symmetry of ISPSYM
* INVELM        : Inverse elements
*
* The irreps are organized in the usual way with the boson irreps first
*
      DIMENSION IACAR(2),ITPAR(2)
      CALL QENTER('RS2A')
*
      IMOC = 0
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' RSBB2A speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' ISCSM and ICCSM = ', ISCSM,ICCSM
        WRITE(6,*) ' IBU = ', IBU
C?      WRITE(6,*) ' MXSXST MXSXBL = ', MXSXST,MXSXBL
        WRITE(6,*) ' RSBB2A : IMOC ', IMOC
        WRITE(6,*) ' ISOC and ICOC : '
        CALL IWRTMA(ISOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ICOC,1,NGAS,1,NGAS)
        WRITE(6,*) ' NTESTG and NSEL2E', NTESTG, NSEL2E
      END IF
*
      IFRST = 1
      JFRST = 1
*
*.Types of DX that connects the two strings
*
C     IDXSM = STSTDX(ISCSM,ICCSM)
      IDXSM = IDBGMULT(ISCSM,INVELM(ICCSM))
      IF(IDXSM.EQ.0)  GOTO 2001
*. Connecting double excitations
      CALL DXTYP_GAS(NDXTYP,ITP,JTP,KTP,LTP,NGAS,ISOC,ICOC)
      DO 2000 IDXTYP = 1, NDXTYP
        ITYP = ITP(IDXTYP)
        JTYP = JTP(IDXTYP)
        KTYP = KTP(IDXTYP)
        LTYP = LTP(IDXTYP)
        NIJKL1 = 0
        IF(ITYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(JTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(KTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(LTYP.EQ.1) NIJKL1 = NIJKL1+1
        IF(NIJKL1.EQ.0) CALL QENTER('BB2A0')
        IF(NIJKL1.EQ.1) CALL QENTER('BB2A1')
        IF(NIJKL1.EQ.2) CALL QENTER('BB2A2')
        IF(NIJKL1.EQ.3) CALL QENTER('BB2A3')
        IF(NIJKL1.EQ.4) CALL QENTER('BB2A4')
*. Integrals included ??
        IF(NSEL2E.NE.0) THEN
          IAMOKAY=0
          IF(ITYP.EQ.JTYP.AND.ITYP.EQ.KTYP.AND.ITYP.EQ.LTYP) THEN
            DO JSEL2E = 1, NSEL2E
              IF(ISEL2E(JSEL2E).EQ.ITYP)IAMOKAY=1
            END DO
          END IF
          IF(IAMOKAY.EQ.0) GOTO 2000
        END IF
*. Largest possible number of orbital pairs
        MI = 0
        MJ = 0
        MK = 0
        ML = 0
        DO IOBSM = 1, NSMST
          MI = MAX(MI,NOBPTS(ITYP,IOBSM))
          MJ = MAX(MJ,NOBPTS(JTYP,IOBSM))
          MK = MAX(MK,NOBPTS(KTYP,IOBSM))
          ML = MAX(ML,NOBPTS(LTYP,IOBSM))
        END DO
        MXPAIR = MAX(MI*MK,MJ*ML)
C       WRITE(6,*) ' MXPAIR  = ', MXPAIR
*. Largest posssible
*. Symmetry of allowed Double excitation, loop over excitations
        DO 1950 IKOBSM = 1, NSMOB
C         JLOBSM = SXDXSX(IKOBSM,IDXSM)
          JLOBSM = IDBGMULT(IDXSM,INVELM(IKOBSM))
C         WRITE(6,*) ' IDXSM IKOBSM JLOBSM ',IDXSM,IKOBSM,JLOBSM
          IF(JLOBSM.EQ.0) GOTO 1950
*. types + symmetries defined => K strings are defined
          KFRST = 1
          K2FRST = 1
*
*. Number of batchs of symmetry pairs IK
*
          LENGTH = 0
          NIKBT = 0
          NBLK = 0
          NBLKT = 0
          DO ISM = 1, NSMOB
C           KSM = ADSXA(ISM,IKOBSM)
            KSM = IDBGMULT(IKOBSM,INVELM(ISM))
            NI = NOBPTS(ITYP,ISM)
            NK = NOBPTS(KTYP,KSM)
            IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
              NIK = NI*(NI+1)/2
C           ELSE IF(ITYP.GT.KTYP.OR.ITYP.EQ.KTYP.AND.ISM.GT.KSM) THEN
            ELSE IF( IOBPTS(ITYP,ISM).GT.IOBPTS(KTYP,KSM)) THEN
              NIK = NI*NK
            ELSE
              NIK = 0
            END IF
            IF(NIK.NE.0) THEN
              NBLKT = NBLKT + 1
              IF(LENGTH+NIK .GT. MXPAIR) THEN
*. The present batch is complete
                NIKBT = NIKBT+1
                IKBT(1,NIKBT) = NBLKT - NBLK
                IKBT(2,NIKBT) = NBLK
                IKBT(3,NIKBT) = LENGTH
                LENGTH = 0
                NBLK = 0
              END IF
              NBLK = NBLK + 1
              LENGTH = LENGTH + NIK
              IKSMBT(1,NBLKT) = ISM
              IKSMBT(2,NBLKT) = KSM
            END IF
          END DO
*. The last batch
          IF(NBLK.NE.0) THEN
            NIKBT = NIKBT+1
            IKBT(1,NIKBT) = NBLKT - NBLK + 1
            IKBT(2,NIKBT) = NBLK
            IKBT(3,NIKBT) = LENGTH
          END IF
*.
          IF(NTEST.GE.2000) THEN
            WRITE(6,*) ' ITYP, KTYP, IKOBSM,  NIKBT = ',
     &                   ITYP, KTYP, IKOBSM,  NIKBT
            WRITE(6,*) ' IKBT : Offset, number, length '
            DO JIKBT = 1, NIKBT
              WRITE(6,'(3i3)') (IKBT(II,JIKBT), II = 1, 3)
            END DO
            WRITE(6,*) ' IKSMBT '
            CALL IWRTMA(IKSMBT,2,NBLKT,2,8)
          END IF
*
*. Number of batchs of symmetry pairs JL
*
          LENGTH = 0
          NJLBT = 0
          NBLK = 0
          NBLKT = 0
          DO JANISM = 1, NSMOB
C           LSM = ADSXA(JSM,JLOBSM)
            LANISM = IDBGMULT(JLOBSM,INVELM(JANISM))
*. JANI and LANI are annihilation symmetries, creation symmetries
            JSM = IADJSYM(JANISM)
            LSM = IADJSYM(LANISM)
C           WRITE(6,*) 'JLOBSM,JANISM,LANISM, JSM,LSM',
C    &                 JLOBSM,JANISM,LANISM, JSM,LSM

*
            NJ = NOBPTS(JTYP,JSM)
            NL = NOBPTS(LTYP,LSM)
            IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
              NJL = NJ*(NJ+1)/2
C           ELSE IF(JTYP.GT.LTYP.OR.JTYP.EQ.LTYP.AND.JSM.GT.LSM) THEN
            ELSE IF (IOBPTS(JTYP,JSM).GT.IOBPTS(LTYP,LSM) ) THEN
              NJL = NJ*NL
            ELSE
              NJL = 0
            END IF
            IF(NJL.NE.0) THEN
              NBLKT = NBLKT + 1
              IF(LENGTH+NJL .GT. MXPAIR) THEN
*. The present batch is complete
                NJLBT = NJLBT+1
                JLBT(1,NJLBT) = NBLKT - NBLK
                JLBT(2,NJLBT) = NBLK
                JLBT(3,NJLBT) = LENGTH
                LENGTH = 0
                NBLK = 0
              END IF
              NBLK = NBLK + 1
              LENGTH = LENGTH + NJL
              JLSMBT(1,NBLKT) = JSM
              JLSMBT(2,NBLKT) = LSM
            END IF
          END DO
*. The last batch
          IF(NBLK.NE.0) THEN
            NJLBT = NJLBT+1
            JLBT(1,NJLBT) = NBLKT - NBLK + 1
            JLBT(2,NJLBT) = NBLK
            JLBT(3,NJLBT) = LENGTH
          END IF

*.
          IF(NTEST.GE.2000) THEN
            WRITE(6,*) ' JTYP, LTYP, JLOBSM,  NJLBT = ',
     &                   JTYP, LTYP, JLOBSM,  NJLBT
            WRITE(6,*) ' JLBT : Offset, number, length '
            DO JJLBT = 1, NJLBT
              WRITE(6,'(3i3)') (JLBT(II,JJLBT), II = 1, 3)
            END DO
            WRITE(6,*) ' JLSMBT '
            CALL IWRTMA(JLSMBT,2,NBLKT,2,8)
          END IF
*
*. Loop over batches of IK strings
          DO 1940 IKBTC = 1, NIKBT
          IF(NTEST.GE.1000) WRITE(6,*) ' IKBTC = ', IKBTC
*. Loop over batches of JL strings
          DO 1930 JLBTC = 1, NJLBT
*
CMOVED        KFRST = 1
CMOVED        K2FRST = 1
*
              IF(IMOC.EQ.0) THEN
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
              IFIRST = 1
*. Loop over batches of I strings
              NPART = NROW/MAXI
              IF(NPART*MAXI.NE.NROW) NPART = NPART + 1
              IF(NTEST.GE.2000)
     &        write(6,*) ' NROW, MAXI NPART ',NROW,MAXI,NPART
              DO 1801 IPART = 1, NPART
                IBOT = 1+(IPART-1)*MAXI
                ITOP = MIN(IBOT+MAXI-1,NROW)
                NIBTC = ITOP-IBOT+1
*.Loop over batches of intermediate strings
                KBOT = 1- MAXK
                KTOP = 0
 1800           CONTINUE
                  KBOT = KBOT + MAXK
                  KTOP = KTOP + MAXK
*
                  IONE = 1
                  JLBOFF = 1
                  NJLT = JLBT(3,JLBTC)
                  DO JLPAIR = 1, JLBT(2,JLBTC)
                    JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                    LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                    NJ = NOBPTS(JTYP,JSM)
                    NL = NOBPTS(LTYP,LSM)
                    IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                      NJL = NJ*(NJ+1)/2
                      JLSM = 1
                    ELSE
                      NJL = NJ * NL
                      JLSM = 0
                    END IF
C                   WRITE(6,*) ' JSM LSM JLSM ' , JSM, LSM,JLSM

*
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain all double excitations from this group of K strings
CT                CALL QENTER('ADADS')
                  I12C = 1
                  K12C = 1
                  IONE = 1
                  CALL ADADST_GAS(IONE,JSM,JTYP,NJ,
     &                            IONE,LSM,LTYP,NL,
     &                        ICCTP,ICCSM,IGRP,
     &                        KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                        JFRST,KFRST,I12C,K12C,SCLFAC)
                  JFRST = 0
                  KFRST = 0
*
CT                CALL QEXIT('ADADS')
C                 WRITE(6,*) ' NKBTC after ADADST_GAS',NKBTC
                  IF(NKBTC.EQ.0) GOTO 1930
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
CT                CALL QENTER('MATCG')
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
*.CB(IA,KB,jl) = +/-C(IA,a+la+jIA)
                    JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                            I1(1,I1JL),XI1S(1,I1JL))
                    END IF
                  END DO
CT                CALL QEXIT ('MATCG')
*
                  JLBOFF = JLBOFF + NJL
                  END DO
*. ( End of loop over jl in batch )
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals xint(ik,jl) = (ij!kl)-(il!kj)
                  IF(IFIRST.EQ.1) THEN
                    IXCHNG = 1
* Obtain integrals in ik batch
                    NIKT = IKBT(3,IKBTC)
                    NJLT = JLBT(3,JLBTC)
                    JLOFF = 1
                    DO JLPAIR = 1, JLBT(2,JLBTC)
                    IKOFF = 1
                    DO IKPAIR = 1, IKBT(2,IKBTC)
*
                      ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                      KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                      JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                      LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
*
                      IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                        IKSM = 1
                        NIK = NOBPTS(ITYP,ISM)*(NOBPTS(ITYP,ISM)+1)/2
                      ELSE
                        IKSM = 0
                        NIK = NOBPTS(ITYP,ISM)*NOBPTS(KTYP,KSM)
                      END IF
*
                      IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                        JLSM = 1
                        NJL = NOBPTS(JTYP,JSM)*(NOBPTS(JTYP,JSM)+1)/2
                      ELSE
                        JLSM = 0
                        NJL = NOBPTS(JTYP,JSM)*NOBPTS(LTYP,LSM)
                      END IF

                      ICOUL = 0
                      CALL GETINT_DBG(SCR,ITYP,ISM,IBU,JTYP,JSM,IBU,
     &                                KTYP,KSM,IBU,LTYP,LSM,IBU,
     &                                IRI,IXCHNG,IKSM,JLSM,ICOUL)

C                     CALL GETINT(SCR,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,
C    &                          LSM,IXCHNG,IKSM,JLSM,ICOUL)
                      DO JL = 1, NJL
                        CALL COPVEC(SCR((JL-1)*NIK+1),
     &                             XINT((JLOFF-1+JL-1)*NIKT+IKOFF),NIK)
                      END DO
                      IKOFF = IKOFF + NIK
                    END DO
                    JLOFF = JLOFF + NJL
                    END DO
                  END IF
                  IFIRST = 0
*.and now ,to the work
                  LIKB = NIBTC*NKBTC
                  IF(NTEST.GE.3000) THEN
                   WRITE(6,*) ' Integral block '
                   CALL WRTMAT(XINT,NIKT,NJLT,NIKT,NJLT)
                  END IF
                  IF(NTEST.GE.3000) THEN
                    WRITE(6,*) ' CSCR matrix '
                    CALL WRTMAT(CSCR,LIKB,NJLT,LIKB,NJLT)
                  END IF
*
                  FACTORC = 0.0D0
                  FACTORAB = 1.0D0
                  CALL MATML7(SSCR,CSCR,XINT,
     &                        LIKB,NIKT,LIKB,NJLT,NIKT,NJLT,
     &                        FACTORC,FACTORAB,2)
                  IF(NTEST.GE.3000) THEN
                    WRITE(6,*) ' SSCR matrix '
                    CALL WRTMAT(SSCR,LIKB,NIKT,LIKB,NIKT)
                    WRITE(6,*) ' Mem check '
                    CALL LMEMCHK
                  END IF
* ============================
* Loop over ik and scatter out
* ============================
*. Generate double excitations from K strings
*. I strings connected with K strings in batch <I!a+i a+k!K)
                  IF(NTEST.GE.3000) THEN
                    WRITE(6,*) ' S block second call to ADADST_GAS'
                    CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                  END IF
                  I12C = 2
*
                  IONE = 1
                  IKBOFF = 1
                  DO IKPAIR = 1, IKBT(2,IKBTC)
                    ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                    KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                    NI = NOBPTS(ITYP,ISM)
                    NK = NOBPTS(KTYP,KSM)
                    IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                      NIK = NI*(NI+1)/2
                      IKSM = 1
                    ELSE
                      NIK = NI * NK
                      IKSM = 0
                    END IF
CT                  CALL QENTER('ADADS')
                    IF(IFRST.EQ.1) KFRST = 1
                    ONE = 1.0D0
                    CALL ADADST_GAS(IONE,ISM,ITYP,NI,
     &                              IONE,KSM,KTYP,NK,
     &                              ISCTP,ISCSM,IGRP,
     &                              KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                              IFRST,KFRST,I12C,K12C,ONE)
*
                    IFRST = 0
                    K2FRST = 0
CT                  CALL QEXIT ('ADADS')
*
CT                  CALL QENTER('MATCS')
                    I = 0
                    K = 1
                    IF(NTEST.GE.3000) THEN
                      WRITE(6,*) ' S block before MATCA '
                      CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                    END IF
                    DO IK = 1, NIK
                      CALL NXTIJ(I,K,NI,NK,IKSM,NONEW)
                      IKOFF = (K-1)*NI + I
                      ISBOFF = 1+(IKBOFF-1+IK-1)*NIBTC*NKBTC
                      IF(IKSM.EQ.1.AND.I.EQ.k) THEN
* a+ i a+i gives trivially zero
                      ELSE
                        CALL MATCAS(SSCR(ISBOFF),SB,NIBTC,NROW,IBOT,
     &                              NKBTC,I1(1,IKOFF),XI1S(1,IKOFF))
                      END IF
                    END DO
CT                  CALL QEXIT ('MATCS')
                    IKBOFF = IKBOFF + NIK
                    IF(NTEST.GE.3000) THEN
                     WRITE(6,*) ' Updated S block'
                     CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                    END IF
                  END DO
*
                IF(KEND.EQ.0) GOTO 1800
*. End of loop over partitionings of resolution strings
 1801         CONTINUE
                ELSE IF (IMOC.EQ.1) THEN
*
* =========================================================================
*                    Use MOC method
* =========================================================================
*
* Fetch two electron integrals as (i,k,j,l)
                  IXCHNG = 1
                    ICOUL = 0
                    CALL LGETINT(XINT,ITYP,ISM,JTYP,JSM,KTYP,KSM,LTYP,
     &                          LSM,IXCHNG,0,0,ICOUL     )


*. Generate off diagonal I => K single excitations
              CALL SXSTR(ISCSM,ISCTP,IGRP,LSM,LTYP,
     &                   KSM,KTYP,MXSXST,RIKSX(1,1),RIKSX(1,2),
     &                   RIKSX(1,3),RIKSX(1,4),I1,1,NTEST)
*. Generate off-diagonal K => J single excitations
      IACAR(1) = 2
      IACAR(2) = 1
*. We are actually using conjugated operators so
      ITPAR(1) = LTYP
      ITPAR(2) = KTYP
*     ITPAR(1) = KTYP
*     ITPAR(2) = LTYP
      CALL NEWTYP(ISCTP,IACAR,ITPAR,2,KSTTYP)
C     NEWTYP(INSPGP,IACOP,ITPOP,NOP,OUTSPGP)
C      SYMCOM(ITASK,IOBJ,I1,I2,I12)
      CALL SYMCOM(3,0,KLSM,ISCSM,KSTRSM)
C?    WRITE(6,*) ' KSTRSM KSTTYP = ', KSTRSM,KSTTYP
*. Generate off diagonal K=> J single excitations
              CALL SXSTR(KSTRSM,KSTTYP,IGRP,JSM,JTYP,
     &                   ISM,ITYP,MXSXST,RKJSX(1,1),RKJSX(1,2),
     &                   RKJSX(1,3),RKJSX(1,4),XI1S,1,NTEST)
C?     write(6,*) ' memory check in after second SXSTR '
C?    write(6,*) ' MXSXST 3 ', MXSXST
C?     call memchk
*. Combine single excitations to double excitations
*. Loop over columns trings of sigma
       DO ISTR = 1, NSCOL
*. Double excitations out from DXFSTR
       CALL DXFSTR(RIKSX(1,1),RIKSX(1,2),RIKSX(1,3),RIKSX(1,4),
     &             I1,
     &             RKJSX(1,1),RKJSX(1,2),RKJSX(1,3),RKJSX(1,4),
     &             XI1S,
     &             NI,NJ,NK,NL,
     &             IOFF,JOFF,KOFF,LOFF,
     &             ISTR,MXSXST,XINT,
     &             I2,XI2S,NEX,NTEST)
C?     write(6,*) ' memory check in after DXSTR        '
C?     call memchk
C?    write(6,*) ' MXSXST 4 ', MXSXST
*. Msaxpy update
        CALL MSAXPYN6(SB((ISTR-1)*NROW+1),XI2S,CB,NROW,
     &                NEX,I2)
C      write(6,*) ' memory check in RSBB2A '
C      call memchk
        END DO
              END IF
 1930       CONTINUE
 1940     CONTINUE
 1950   CONTINUE
        IF(NIJKL1.EQ.0) CALL QEXIT ('BB2A0')
        IF(NIJKL1.EQ.1) CALL QEXIT ('BB2A1')
        IF(NIJKL1.EQ.2) CALL QEXIT ('BB2A2')
        IF(NIJKL1.EQ.3) CALL QEXIT ('BB2A3')
        IF(NIJKL1.EQ.4) CALL QEXIT ('BB2A4')
 2000 CONTINUE
*
 2001 CONTINUE
      CALL QEXIT('RS2A ')
      RETURN
      END
