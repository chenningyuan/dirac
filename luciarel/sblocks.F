      SUBROUTINE SBLOCKS(ONEINTS,
     &                NSBLOCK,ISBLOCK,CB,SB,C2,
     &                ICOCOC,ICSMOS,ICBLTP,NACOB,NSSOA_S,NSSOB_S,
     &                NSSOA_C,NSSOB_C,
     &                NAEL_S,IAGRP_S,NBEL_S,IBGRP_S,
     &                NAEL_C,IAGRP_C,NBEL_C,IBGRP_C,
     &                IOCTPA_S,IOCTPB_S,NOCTPA_S,NOCTPB_S,
     &                IOCTPA_C,IOCTPB_C,NOCTPA_C,NOCTPB_C,
     &                NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,
     &                NOBPTS2,IOBPTS2,MXPNGAS,
     &                ITSOB,MAXIJ,MAXK,MAXI,LI,LC,LS,
     &                XINT,CSCR,SSCR,SXSTSM,STSTSX,STSTDX,
     &                SXDXSX,ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &                MXTSOB,SOMESCR,XIJILS,
     &                I1,XI1S,I2,XI2S,IDOH2,MXPOBS,
     &                PS,IPRNT,LUC,ICJKAIB,CJRES,SIRES,I3,XI3S,
     &                I4,XI4S,MXSXST,MXSXBL,MOCAA,MOCAB,IAPR,
     &                LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,IRESTRICT,
     &                ICONSPA,ICONSPB,SCLFAC,
     &                LUCBLK,IPERTOP,IH0INSPC,ISPINFREE,
CTF  &                IH0SPC,
     &                ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                IRELACI,IRIS,IRIC,I_DO_TRNSP,IELCI)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen , Winter of 1991
*
* Modified for NO-PAIR DBG CI , T.F + J.O, March 21 98
*
* =====
* Input
* =====
*
* NSBLOCK : Number of BLOCKS included
* ISBLOCK : Blocks included
*   ISBLOCK(1,*) : alpha type of block
*   ISBLOCK(2,*) : beta type of block
*   ISBLOCK(3,*) : sym of alpha in block
*   ISBLOCK(4,*) : Offset of block
*
* ICOCOC : Allowed type combinations for C
* ICSMOS : Symmetry array for C
* ICBLTP : Block types for C
* NACOB : Number of active orbitals
* NSSOA : Number of strings per type and symmetry for alpha strings
* NAEL  : Number of active alpha electrons
* NSSOB : Number of strings per type and symmetry for beta strings
* NBEL  : Number of active beta electrons
* NTSOB : Number of orbitals per type and symmetry
* NOBPTS: Orbitals of given type and symmetry
* IOBPTS: Offset for orbitals of given sym and type
*
* MAXIJ : Largest allowed number of orbital pairs treated simultaneously
* MAXK  : Largest number of N-2,N-1 strings treated simultaneously
* MAXI  : Max number of N strings treated simultaneously
*
* LI : Length of scratch array for integrals
* LC : Length of scratch array for C
* LS : Length of scratch array for S
* XINT : Scratch array for integrals
* CSCR : Scratch array for C vector
* SSCR : Scratch array for S vector
*
* CJRES,SIRES : Space for above matrices
* The C and S vectors are accessed through routines that
* either fetches/disposes symmetry blocks or
* Symmetry-occupation-occupation blocks
*
*
* Note : Prescreen + Restrict features eliminated
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Specific input
      INTEGER ISBLOCK(8,*)
#include "picasso.inc"
      dimension ONEINTS(INTSPACE)
*.General input
      INTEGER ICOCOC(NOCTPA_C,NOCTPB_C)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA_S(NSMST,*),NSSOB_S(NSMST,*)
      INTEGER NSSOA_C(NSMST,*),NSSOB_C(NSMST,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST),ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NOBPTS2(MXPNGAS,*),IOBPTS2(MXPNGAS,*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA_S,NOCTPA_S),ICONSPB(NOCTPB_S,NOCTPB_S)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION XINT(*),CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*),
     &          XI1S(*),XI2S(*),XI3S(*),XI4S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      dimension ISTRFL(12)
      DIMENSION SCLFAC(*)
*
      REAL * 8 INPROD
*
      CALL QENTER('SBLOCK_S')
C     WRITE(6,*) ' LMEMCHK at start of SBLOCKS'
C     CALL LMEMCHK
      NTEST = 000
      NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' SBLOCKS speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' Number of sigma blocks to be calculated ',
     &  NSBLOCK
        WRITE(6,*) ' TTSS for each ACTIVE sigma block'
          DO IBLOCK = 1, NSBLOCK
            IF(ISBLOCK(1,IBLOCK).GT.0)
     &      WRITE(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          END DO
          WRITE(6,*) ' IDC PS ', IDC,PS
          WRITE(6,*) ' IDOH2 = ',IDOH2
*
          WRITE(6,*) ' NAEL_C, NBEL_C, NAEL_S,NBEL_S',
     &                 NAEL_C, NBEL_C, NAEL_S,NBEL_S
      END IF
*
      NTEST = 100
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Initial C vector '
        CALL WRTVCD(CB,LUC,1,-1)
      END IF
      NTEST = 00
* ===========================
* 1 : Arrays for accessing C
* ============================
*. Find batches of C - strings
      CALL PART_CIV2(IDC,ICBLTP,NSSOA_C,NSSOB_C,NOCTPA_C,NOCTPB_C,NSMST,
     &               LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &               I1CBLOCK,ICBLOCK,0,1)
*. Find the active blocks on LUC, store info in SCLFAC
C     FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
*     \/ Change to accomodate for ICBAT_RES
      CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
C     write(6,*) ' memchk after FIND_ACTIVE_BLOCKS'
C     CALL LMEMCHK
*
* Loop over batches over C blocks
      IF(IDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
C?    WRITE(6,*) ' IDOH2, MXEXC',IDOH2,MXEXC
      REWIND LUC
      IF(ICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
        if (NTEST.ge.1000) then
          write(6,*) 'Batch info from SBLOCKS:'
          write(6,*) 'JCBATCH = ',JCBATCH
          write(6,*) 'JCBAT_END = ',JCBAT_END
        end if
*
*. Read C blocks of this batch into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
        DO JJCBLOCK = 1, NJBLOCK
          if (NTEST.ge.1000) then
            write(6,*) 'Block info from SBLOCKS:'
            write(6,*) 'JJCBLOCK = ',JJCBLOCK
            write(6,*) 'NJBLOCK = ',NJBLOCK
          end if
*
          JBLOCK = I1CBLOCK(JCBATCH)-1+JJCBLOCK
          JATP = ICBLOCK(1,JBLOCK)
          JBTP = ICBLOCK(2,JBLOCK)
          JASM = ICBLOCK(3,JBLOCK)
          JBSM = ICBLOCK(4,JBLOCK)
          JOFF = ICBLOCK(5,JBLOCK)
          ISCALE = 1
*!!!
CTF
*  Caution: C and ICOOSC are not initialized at the moment.
*           But we read CB in from disc, so they are currently
*           NOT used in GSTTBL. If LUC is zero, they will be
*           used, i.e. C is input array for CB !
CTF
          CALL GSTTBL(C,CB(JOFF),JATP,JASM,JBTP,JBSM,ICOCOC,
     &                NOCTPA_C,NOCTPB_C,NSSOA_C,NSSOB_C,PS,
     &                ICOOSC,IDC,
     &                PL,LUC,C2,NSMST,ISCALE,SCLFAC(JBLOCK))
*
          if (NTEST.ge.1000) then
            NJATEST = NSSOA_C(JASM,JATP)
            NJBTEST = NSSOB_C(JBSM,JBTP)
            write(6,*) 'after GSTTBL in sblocks'
            write(6,*) 'JOFF is  ',JOFF
            call wrtmat(CB(JOFF),NJBTEST,NJATEST,
     &                  NJBTEST,NJATEST)
          end if
*
          SCLFAC(JBLOCK) = 1.0d0
        END DO
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
        DO 9000 ICBLK = I1CBLOCK(JCBATCH),
     &                I1CBLOCK(JCBATCH)-1+NJBLOCK
          if (NTEST.ge.1000) then
            write(6,*) 'Block info from SBLOCKS:'
            write(6,*) 'ICBLK = ',ICBLK
            write(6,*) 'I1CBLOCK(JCBATCH)-1+NJBLOCK = ',
     &                  I1CBLOCK(JCBATCH)-1+NJBLOCK
          end if
*
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          JASM = ICBLOCK(3,ICBLK)
          JBSM = ICBLOCK(4,ICBLK)
          ICOFF = ICBLOCK(5,ICBLK)
          NJA = NSSOA_C(JASM,JATP)
          NJB = NSSOB_C(JBSM,JBTP)
*
          IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
            CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                  LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
            DO 8765 IPERM = NPERM,1,-1
              if (NTEST.ge.1000) then
                write(6,*) 'Permutation info from SBLOCKS:'
                write(6,*) 'IPERM = ',IPERM
                write(6,*) 'NPERM = ',NPERM
              end if
*
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              NLLA = NSSOA_C(LLASM,LLATP)
              NLLB = NSSOB_C(LLBSM,LLBTP)
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 always corresponds to no transpose, so transpose!
              IF(IPERM.EQ.1) THEN
                IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM)
     &                   THEN
*. Diagonal blocks, Transposing corresponds to scaling
                  IF(PS.EQ.-1.0D0) THEN
                    CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                  END IF
                ELSE
*. offdiagonal blocks, explicit transposing
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                  CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
                  if (LSGN(IPERM).EQ.-1) then
                    write(6,*) 'LCOL and LROW undefined in sblocks.'
                    Call Abend2( 'Quitting.' )
CTF                 CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                  end if
                END IF
              END IF
*
              DO 10000 ISBLK = 1, NSBLOCK
                if (NTEST.ge.1000) then
                  write(6,*) 'Block info from SBLOCKS:'
                  write(6,*) 'ISBLK = ',ISBLK
                  write(6,*) 'NSBLOCK = ',NSBLOCK
                end if
*
                IF(ISBLOCK(1,ISBLK) .GT. 0 ) THEN
                  IATP = ISBLOCK(1,ISBLK)
                  IBTP = ISBLOCK(2,ISBLK)
                  IASM = ISBLOCK(3,ISBLK)
                  IBSM = ISBLOCK(4,ISBLK)
                  ISOFF = ISBLOCK(5,ISBLK)
                  NIA = NSSOA_S(IASM,IATP)
                  NIB = NSSOB_S(IBSM,IBTP)
*
                  IF(NIA*NIB.EQ.0) GOTO 10000
                  IF(NJA*NJB.EQ.0) GOTO 10000
*
                  CNORM = INPROD(CB(ICOFF),CB(ICOFF),NJA*NJB)
                  IF(CNORM.LT.1.0D-15) GOTO 10000
*
                  IF(NTEST.GE.100) THEN
                    WRITE(6,*) ' Next s block in batch : '
                    write(6,*) ' ISBLK IASM IBSM IATP IBTP'
                    write(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                  END IF
*
                  IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP.AND.
     &              ((LLBSM.GT.LLASM).OR.
     &              (LLASM.EQ.LLBSM).AND.(LLBTP.GT.LLATP)))
     &              GOTO 8764
*
                  IF(NTEST.GE.60) THEN
                    WRITE(6,*) ' RSSBCB will be called for '
                    WRITE(6,*) ' Sigma block : '
                    WRITE(6,*) ' ISOFF ', ISOFF
                    WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                    WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                    WRITE(6,*) ' C     block : '
                    WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                    WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                    WRITE(6,*) ' ICOFF ', ICOFF
                    WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                  END IF
*
                  XFAC = SCLFAC(ICBLK)
*
                  if (IRELACI.eq.0) then
*. Nonrelativistic Hamiltonian
                    CALL RSSBCB2(IASM,IATP,IOCTPA_S,
     &                 IBSM,IBTP,IOCTPB_S,
     &                 LLASM,LLATP,LLBSM,LLBTP,NGAS,
     &                 NELFSPGP(1,IATP+IOCTPA_S-1),
     &                 NELFSPGP(1,IBTP+IOCTPB_S-1),
     &                 NELFSPGP(1,LLATP+IOCTPA_S-1),
     &                 NELFSPGP(1,LLBTP+IOCTPB_S-1),
     &                 NAEL_S,NBEL_S,
     &                 IAGRP_S,IBGRP_S,
     &                 SB(ISOFF),CB(ICOFF),IDOH2,
     &                 ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                 NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
     &                 SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                 NSMOB,NSMST,NSMSX,NSMDX,
     &                 NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
     &                 ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                 MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
     &                 IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC,
     &                 MXTSOB,SOMESCR,XIJILS)
                  else if (IRELACI.eq.1) then
*. No-pair relativistic Hamiltonian
                    IATPABS = IATP +IOCTPA_S-1
                    IBTPABS = IBTP +IOCTPB_S-1
                    JATPABS = JATP +IOCTPA_C-1
                    JBTPABS = JBTP +IOCTPB_C-1
*
                    NTEST = 200
                    if (NTEST.ge.200) then
                      write(6,*) 'String type offsets in sblocks:'
                      write(6,*) 'IATPABS = ',IATPABS
                      write(6,*) 'IBTPABS = ',IBTPABS
                      write(6,*) 'JATPABS = ',JATPABS
                      write(6,*) 'JBTPABS = ',JBTPABS
                      write(6,*) 'IATP,IOCTPA_S',IATP,IOCTPA_S
                      write(6,*) 'IBTP,IOCTPB_S',IBTP,IOCTPB_S
                      write(6,*) 'JATP,IOCTPA_C',JATP,IOCTPA_C
                      write(6,*) 'JBTP,IOCTPB_C',JBTP,IOCTPB_C
                    end if
*
                    call rssbcbrel(ONEINTS,INTSPACE,
     &         IASM,IATP,IOCTPA_S,NELFSPGP(1,IATPABS),IAGRP_S,NIA,
     &         IBSM,IBTP,IOCTPB_S,NELFSPGP(1,IBTPABS),IBGRP_S,NIB,
     &         JASM,JATP,IOCTPA_C,NELFSPGP(1,JATPABS),IAGRP_C,NJA,
     &         JBSM,JBTP,IOCTPB_C,NELFSPGP(1,JBTPABS),IBGRP_C,NJB,
     &                  NGAS,SB(ISOFF),CB(ICOFF),IDOH2,
     &                  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                  NOBPTS,IOBPTS,NOBPTS2,IOBPTS2,
     &                  MXPNGAS,MAXI,MAXK,
     &                  SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
     &                  NSMOB,NSMST,NSMSX,NSMDX,
     &                  MXPOBS,IDC,
     &                  CJRES,SIRES,I3,XI3S,I4,XI4S,
     &                  MXSXBL,MXSXST,IPRNT,ISPINFREE,
     &                  IPTSPC,JPTSPC,IJOP,SCLFAC,
     &                  IRIS,IRIC,IELCI,
     &                  MXTSOB,SOMESCR,XIJILS)
C                   WRITE(6,*) 'before CB(19) = ', CB(19)
                  else
                    write(6,*) 'Something is wrong in SBLOCKS'
                    write(6,*) 'RELACI corrupted. IRELACI = ',IRELACI
                    Call Abend()
                  end if
*                 ^ End of switch relativistic/non-relativistic
 8764             CONTINUE        ! if great if block is true
                END IF
*               ^ end if ISBLOCK > 0
10000         CONTINUE
*.            ^ End of loop over sigma blocks
 8765       CONTINUE
*           ^ end loop start with transpose block
          END IF
*         ^ end  if SCLFAC > 0
 9000   CONTINUE
*.      ^ End of loop over C blocks in Batch
20000 CONTINUE
*.    ^End of loop over batches of C blocks
C     write(6,*) ' LMEMCHK after loop 20000'
C     CALL LMEMCHK
*
* Order
* The condition I_DO_TRNSP.EQ.1 is fulfilled if
*  IRIC.EQ.2.AND.IMK2_C.EQ.NMS2VAL
*  in rassg3. This means that the imaginary part (IRIC) of the
*  last MK2 value of the ket vectors (IMK2_C) for a given IRIS
*  and MK2_S is now being treated. After the sigma vector for
*  this partition is created, the final S block is reordered to
*  untransposed form (in the following).
*
      IF(I_DO_TRNSP.EQ.1) THEN
        DO  ISBLK = 1 , NSBLOCK
        IF(ISBLOCK(1,ISBLK).GT.0) THEN
          IATP = ISBLOCK(1,ISBLK)
          IBTP = ISBLOCK(2,ISBLK)
          IASM = ISBLOCK(3,ISBLK)
          IBSM = ISBLOCK(4,ISBLK)
          ISOFF  = ISBLOCK(5,ISBLK)
          ISOFFP = ISBLOCK(6,ISBLK)
          NIA = NSSOA_S(IASM,IATP)
          NIB = NSSOB_S(IBSM,IBTP)
*. Tranpose sigma block was obtained, transpose to obtain correct block
          CALL TRPMT3(SB(ISOFF),NIB,NIA,C2)
C?        WRITE(6,*) 'Transpose loop, NIA,NIB,ISOFF',NIA,NIB,ISOFF
          CALL COPVEC(C2,SB(ISOFF),NIA*NIB)
          IF(IDC.EQ.2.AND.IASM.EQ.IBSM.AND.IATP.EQ.IBTP)
     &    CALL TRPAD3(SB(ISOFF),PS,NIA)
        END IF
        END DO
      END IF
*
      NTEST = 50
      IF(NTEST.GE.50) THEN
        WRITE(6,*) ' output blocks from SBLOCKS '
        CALL WRTTTS(SB,ISBLOCK,NSBLOCK,
     &              NSMST,NOCTPA_S,NOCTPB_S,NSSOA_S,NSSOB_S,1)
      END IF
      NTEST = 0
*
C?    WRITE(6,*) ' LMEMCHK at END of SBLOCKS'
C?    CALL LMEMCHK
      CALL QEXIT('SBLOCK_S')
      RETURN
      END
