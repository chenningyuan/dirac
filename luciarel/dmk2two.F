      SUBROUTINE dmk2two(
     &           ISCSM,ISCTP,LSC,ISCSPGP,
     &           ISRSM,ISRTP,LSR,ISRSPGP,
     &           ICCSM,ICCTP,LCC,ICCSPGP,
     &           ICRSM,ICRTP,LCR,ICRSPGP,
     &           NGAS,ISCOC,ISROC,ICCOC,ICROC,
     &           SB,CB,
     &           ADSXA,
     &           NOBPTS_U,NOBPTS_B,MAXI,MAXK,
     &           SSCR,CSCR,I1,XI1S,I2,XI2S,
     &           XINT,
     &           NSMOB,NSMST,
     &           SCLFAC,NTESTG,IRI,IIDMK2,
     &           SCR)
*
* Delta MK = 2 ( or -2 ) contribution to sigma vector.
*
*
* The Delta MK=2 term is (bar is indicated by a following underscore )
*
* Sigma +2 (Ia,Ib) =
* 1/2 Sum(ijkl) <Ia Ib ! a+i a+k a l_ aj_! Ja Jb> C(ja,Jb) (i j_!k l_)
*
*
* Sigma -2 (Ia,Ib) =
* 1/2 Sum(ijkl) <Ia Ib ! a+i_ a+k_ a l aj! Ja Jb> C(ja,Jb) (i_ j!k_ l)
*
*
* IRI = 1 => Real part of operator
* IRI = 2 => Imaginary part of operator
*
* =====
* Input
* =====
*
* ISCSM,ISCTP,LSC,ISCOC,ISCSPGP
*         : Symmetry,type,length,occ, supergroup of sigma columns
* ISRSM,ISRTP,LSR,ISROC,ISRSPGP
*         : Symmetry,type,length,occ, supergroup of sigma rows
*
* ICCSM,ICCTP,LCC,ICCOC,ICCSPGP
*         : Symmetry,type,length,occ, supergroup of C     columns
* ICRSM,ICRTP,LRC,ICROC,ICRSPGP
*         : Symmetry,type,length,occ, supergroup of C     rows
*
* CB   : Input C block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* DXSTST : Sym of sx,!st> => sym of sx !st>
* STSTDX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
* SXDXSX : Symmetry of SX1,SX1*SX2 => symmetry of SX2
* NOBPTS_B : Number of barred spinors per type and symmetry
* NOBPTS_U : Number of unbarred spinors per type and symmetry
* IADJSYM     :Array symmetry of object => symmetry of adjoint object
* INVELM  : Inverse elements
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
* Jeppe Olsen, Summer of 97
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
*. General input
      INTEGER ADSXA(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS_U(MXPNGAS,*),NOBPTS_B(MXPNGAS,*)
*.Input
      INTEGER ISCOC(NGAS),ICCOC(NGAS),ISROC(NGAS),ICROC(NGAS)
#include "symm.inc"
      DIMENSION CB(*)
*.Output
      DIMENSION SB(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),XINT(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*),I2(MAXK,*),XI2S(MAXK,*)
*
      INTEGER IKBT(3,8),IKSMBT(2,8),JLBT(3,8),JLSMBT(2,8)
      dimension SCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
*
      DIMENSION IACAR(2),ITPAR(2)
      CALL QENTER('DMK2TWO')
      NTESTL = 000
      NTEST = MAX(NTESTG,NTESTL)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' SIGMA_MK2 speaking '
        WRITE(6,*) ' ==================='
*
        WRITE(6,*) ' alpha- and beta- supergroup for sigma: '
        CALL IWRTMA(ISCOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ISROC,1,NGAS,1,NGAS)
*
        WRITE(6,*) ' alpha- and beta- supergroup for sigma: '
        CALL IWRTMA(ISCOC,1,NGAS,1,NGAS)
        CALL IWRTMA(ISROC,1,NGAS,1,NGAS)
*
      END IF
*. Program runs formally as a DELTA MK = + 2 case.
* ( Switch A and be def around in call to get delta MK=-2
* and flag  in IIDMK2 )
      IPM = 1
*
      IBFRST = 1
      JBFRST = 1
      IAFRST = 1
      JAFRST = 1
*. Symmetry and types of connections
*. Symmetry of creation pairs
        IKOBSM = IDBGMULT(INVELM(IADJSYM(ISRSM)),ICRSM)
        JLOBSM = IDBGMULT(INVELM(IADJSYM(ICCSM)),ISCSM)
      IF(IKSM.EQ.0.OR.JLSM.EQ.0) GOTO 2001
*. Connecting creations/Annihilations
C     AATYP(NAATP,ITYP,JTYP,LOCC,ROCC, NGAS)
        CALL AATYP(NKLTYP,KTYP,LTYP,ICROC,ISROC,NGAS)
        CALL AATYP(NIJTYP,ITYP,JTYP,ISCOC,ICCOC,NGAS)
      IF(NIJTYP.EQ.0.OR.NKLTYP.EQ.0) GOTO 2001
*. Largest possible number of orbital pairs
      MI = 0
      MJ = 0
      MK = 0
      ML = 0
      DO IOBSM = 1, NSMST
        MI = MAX(MI,NOBPTS_B(ITYP,IOBSM),NOBPTS_U(ITYP,IOBSM))
        MJ = MAX(MJ,NOBPTS_B(JTYP,IOBSM),NOBPTS_U(JTYP,IOBSM))
        MK = MAX(MK,NOBPTS_B(KTYP,IOBSM),NOBPTS_U(KTYP,IOBSM))
        ML = MAX(ML,NOBPTS_B(LTYP,IOBSM),NOBPTS_U(LTYP,IOBSM))
      END DO
      MXPAIR = MAX(MI*MK,MJ*ML)
C     WRITE(6,*) ' MXPAIR  = ', MXPAIR
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
        KSM = ADSXA(ISM,IKOBSM)
          NI = NOBPTS_U(ITYP,ISM)
          NK = NOBPTS_U(KTYP,KSM)
*
        IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
          NIK = NI*(NI+1)/2
        ELSE IF(ITYP.GT.KTYP.OR.ITYP.EQ.KTYP.AND.ISM.GT.KSM) THEN
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
     &               ITYP, KTYP, IKOBSM,  NIKBT
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
      DO JSM = 1, NSMOB
        LSM = ADSXA(JSM,JLOBSM)
          NJ = NOBPTS_B(JTYP,JSM)
          NL = NOBPTS_B(LTYP,LSM)
*
        IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
          NJL = NJ*(NJ+1)/2
        ELSE IF(JTYP.GT.LTYP.OR.JTYP.EQ.LTYP.AND.JSM.GT.LSM) THEN
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
     &               JTYP, LTYP, JLOBSM,  NJLBT
        WRITE(6,*) ' JLBT : Offset, number, length '
        DO JJLBT = 1, NJLBT
          WRITE(6,'(3i3)') (JLBT(II,JJLBT), II = 1, 3)
        END DO
        WRITE(6,*) ' JLSMBT '
        CALL IWRTMA(JLSMBT,2,NBLKT,2,8)
      END IF
*
*. Loop over batches of IK symmetries
      DO 1940 IKBTC = 1, NIKBT
        IF(NTEST.GE.1000) WRITE(6,*) ' IKBTC = ', IKBTC
*. Loop over batches of JL symmetries
        DO 1930 JLBTC = 1, NJLBT
*
*
* =========================================================================
*                    Use N-2 projection method
* =========================================================================
*
* The evaluation proceeds as
*
* 1 : C(Ja,Ib,lj) = sum(Jb) <Ib!a l_ a j_!Jb> C(Ja,Jb)
* 2 : S(Ja,Ib,ik) = sum(jl) C(Ja,Ib,lj)[(ij_!kl_)-(il_!kj_)]
* 3 : S(Ia,Ib)    = S(Ia,Ib) + Sum(Ja)<Ia!a+i a+k!Ja> S(ja,Ib,Ik)
          IFIRST = 1
*. Loop over batches of I strings, that is beta strings of sigma (I_beta)
          NPART = LSC/MAXI
          IF(NPART*MAXI.NE.LSR ) NPART = NPART + 1
          IF(NTEST.GE.2000)
     &      write(6,*) ' LSR , MAXI NPART ',LSR ,MAXI,NPART
          DO 1801 IPART = 1, NPART
            IBOT = 1+(IPART-1)*MAXI
            ITOP = MIN(IBOT+MAXI-1,NROW)
            NIBTC = ITOP-IBOT+1
*.Loop over batches of Ja strings
            NPART_JA = LRC/MAXK
            IF(NPART_JA.LE.LRC) NPART_JA = NPART_JA + 1
            DO 1799 IPART_JA = 1, NPART_JA
              KBOT = (IPART_JA-1)*MAXK + 1
              KTOP = MIN(KBOT-1+MAXK,NRC)
              NKBTC = KTOP-KBOT+1
*
              IONE = 1
              JLBOFF = 1
              NJLT = JLBT(3,JLBTC)
              DO JLPAIR = 1, JLBT(2,JLBTC)
                JSM = JLSMBT(1,JLBT(1,JLBTC)-1+JLPAIR)
                LSM = JLSMBT(2,JLBT(1,JLBTC)-1+JLPAIR)
                  NJ = NOBPTS_B(JTYP,JSM)
                  NL = NOBPTS_B(LTYP,LSM)
                IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                  NJL = NJ*(NJ+1)/2
                  JLSM = 1
                ELSE
                  NJL = NJ * NL
                  JLSM = 0
                END IF

* 1 : C(Ja,Ib,lj) = sum(Jb) <Ib!a l_ a j_!Jb> C(Ja,Jb)
*. obtain cb(KB,IA,jl) = sum(JB)<KB!a lb a jb !IB>C(IA,JB)
*
*. Obtain double creations from this partitioning of Ib strings
                  I12 = 1
                  K12 = 1
                  IONE = 1
                  CALL ADADST_GAS(IONE,JSM,JTYP,NJ,
     &                     IONE,LSM,LTYP,NL,
     &                     ISCTP,ISCSM,ISCSPGP,
     &                     IBOT,ITOP,I1,XI1S,MAXI,NIBTC,IEND,
     &                     JBFRST,IBFRST,I12,K12,SCLFAC)
                  IBFRST = 0
                  JBFRST = 0
*. Loop over jl in TS classes
                  J = 0
                  L = 1
*
                  DO  IJL = 1, NJL
                    CALL NXTIJ(J,L,NJ,NL,JLSM,NONEW)
                    I1JL = (L-1)*NJ+J
*.CB(JA,IB,jl) = +/-C(JA,a+la+j IB)
                    JLOFF = (JLBOFF-1+IJL-1)*NKBTC*NIBTC+1
                    IF(JLSM.EQ.1.AND.J.EQ.L) THEN
*. a+j a+j gives trivially zero
                      ZERO = 0.0D0
                      CALL SETVEC(CSCR(JLOFF),ZERO,NKBTC*NIBTC)
                    ELSE
                      CALL MATCG(CB,CSCR(JLOFF),NROW,NIBTC,IBOT,NKBTC,
     &                        I1(1,I1JL),XI1S(1,I1JL))
                    END IF
                  END DO
                  JLBOFF = JLBOFF + NJL
              END DO
*.            ^ ( End of loop over jl in batch )
*==============================================
*. SSCR(I,K,ik) = CSR(I,K,jl)*((ij!kl)-(il!jk))
*===============================================
*.Obtain two electron integrals xint(ik,jl) = (ij_!kl_)-(il_!kj_)
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
                      NIK = NOBPTS_U(ITYP,ISM)*(NOBPTS_U(ITYP,ISM)+1)/2
                    ELSE
                      IKSM = 0
                      NIK = NOBPTS_U(ITYP,ISM)*NOBPTS_U(KTYP,KSM)
                    END IF
*
                    IF(JSM.EQ.LSM.AND.JTYP.EQ.LTYP) THEN
                      JLSM = 1
                      NJL = NOBPTS_B(JTYP,JSM)*(NOBPTS_B(JTYP,JSM)+1)/2
                    ELSE
                      JLSM = 0
                       NJL = NOBPTS_B(JTYP,JSM)*NOBPTS_B(LTYP,LSM)
                    END IF

                    ICOUL = 0
                    IXCHNG = 1
*. getint generalized so each index has three parameters
* IBU,ITYP,ISM, IBU = 1 => Barred spinor, IBU = 2 => unbarred spinor
*
                    IF(IIDMK2.EQ.2) THEN
                     IBU = 2
                     JBU = 1
                     KBU = 2
                     LBU = 1
                    ELSE
                     IBU = 1
                     JBU = 2
                     KBU = 1
                     LBU = 2
                    END IF
* (I J_ ! K L_) ( or (I_ J ! K_ L) as (IK!JL)
                    CALL GETINT_DBG(SCR,
     &                     ITYP,ISM,IBU,JTYP,JSM,JBU,
     &                     KTYP,KSM,KBU,LTYP,LSM,LBU,
     &                     IRI,IXCHNG,IKSM,JLSM,ICOUL)
                    DO JL = 1, NJL
                      CALL COPVEC(SCR((JL-1)*NIK+1),
     &                           XINT((JLOFF-1+JL-1)*NIKT+IKOFF),NIK)
                    END DO
                    IKOFF = IKOFF + NIK
                  END DO
*                 ^ End IK pair
                  JLOFF = JLOFF + NJL
                END DO
*               ^ End JL pair
              END IF
*             ^ End if IFIRST
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
     &                    LIKB,NIKT,LIKB,NJLT,NIKT,NJLT,
     &                    FACTORC,FACTORAB,2)
              IF(NTEST.GE.3000) THEN
                WRITE(6,*) ' SSCR matrix '
                CALL WRTMAT(SSCR,LIKB,NIKT,LIKB,NIKT)
                WRITE(6,*) ' Mem check '
                CALL LMEMCHK
              END IF
* ============================
* Loop over ik and scatter out
* ============================
* 3 : S(Ia,Ib)    = S(Ia,Ib) + Sum(Ja)<Ia!a+i a+k!Ja> S(ja,Ib,Ik)
              IF(NTEST.GE.3000) THEN
                WRITE(6,*) ' S block second call to ADADST_GAS'
                CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
              END IF
*
              IONE = 1
              IKBOFF = 1
              DO IKPAIR = 1, IKBT(2,IKBTC)
                ISM = IKSMBT(1,IKBT(1,IKBTC)-1+IKPAIR)
                KSM = IKSMBT(2,IKBT(1,IKBTC)-1+IKPAIR)
                  NI = NOBPTS_U(ITYP,ISM)
                  NK = NOBPTS_U(KTYP,KSM)
*
                IF(ISM.EQ.KSM.AND.ITYP.EQ.KTYP) THEN
                  NIK = NI*(NI+1)/2
                  IKSM = 1
                ELSE
                  NIK = NI * NK
                  IKSM = 0
                END IF
                I12 = 2
                K12 = 2
                ONE = 1.0D0
*. Generate double creations from Batch of Ja strings
                CALL ADADST_GAS(IONE,ISM,ITYP,NI,
     &                          IONE,KSM,KTYP,NK,
     &                          ICRTP,ICRSM,ICRSPGP,
     &                          KBOT,KTOP,I1,XI1S,MAXK,NKBTC,KEND,
     &                          IAFRST,JAFRST,I12,K12,ONE)
*
                IAFRST = 0
                JAFRST = 0
*. Well now, the scattering proceeds as
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
                    continue
                  ELSE
*. at the moment we use straight forward row scattering
C     SCTADRMT(ISCT,SSCT,AIN,AOUT,NROWIN,NROWOUT,NCOL,FACTOR)
                    ONE = 1.0D0
                    CALL SCTADRMT(I1(1,IKOFF),XI1S(1,IKOFF),
     &                   SSCR(ISBOFF),SB,NKBTC,NSR,NIBTC,ONE)
                  END IF
                END DO
*               ^ End of Loop over IK pairs
                IKBOFF = IKBOFF + NIK
                IF(NTEST.GE.3000) THEN
                  WRITE(6,*) ' Updated S block'
                  CALL WRTMAT(SB,NROW,NSCOL,NROW,NSCOL)
                END IF
              END DO
*             ^ End loop over IK pairs
*
 1799       END DO
*.          ^ End of loop over partitionings of Ja strings
 1801     CONTINUE
*.        ^ End of loop over partitionings of Ib strings
 1930   CONTINUE
*       ^ End loop over batches of JL symmetries
 1940 CONTINUE
*     ^ End loop over batches of IK symmetries
 2001 continue
*
      CALL QEXIT('DMK2TWO')
      RETURN
      END
