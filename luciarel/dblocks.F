      subroutine dblocks(NACOB,NACOB2,IRIL,IRIR,
     &                RHO1UUR,RHO1UUI,
     &                RHO1UBR,RHO1UBI,
     &                  RHO1BBR,RHO1BBI,
     &                  RHO1BUR,RHO1BUI,
     &                NSBLOCK,ISBLOCK,CB,SB,C2,
     &                RHO1S,
     &                ICOCOC,ICSMOS,ICBLTP,
     &                NSSOA_S,NSSOB_S,NSSOA_C,NSSOB_C,
     &                NAEL_S,IAGRP_S,NBEL_S,IBGRP_S,
     &                NAEL_C,IAGRP_C,NBEL_C,IBGRP_C,
     &                IOCTPA_S,IOCTPB_S,NOCTPA_S,NOCTPB_S,
     &                IOCTPA_C,IOCTPB_C,NOCTPA_C,NOCTPB_C,
     &                NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,
     &                NOBPTS2,IOBPTS2,MXPNGAS,
     &                ITSOB,MAXIJ,MAXK,MAXI,LC,
     &                LS,CSCR,SSCR,
     &                SXSTSM,STSTSX,STSTDX,SXDXSX,
     &                ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &                I1,XI1S,
     &                I2,XI2S,IDOH2,MXPOBS,
     &                PS,IPRNT,LUC,ICJKAIB,CJRES,
     &                SIRES,I3,XI3S,
     &                I4,XI4S,MXSXST,MXSXBL,
     &                MOCAA,MOCAB,IAPR,
     &                LCBLOCK,LECBLOCK,I1CBLOCK,ICBLOCK,
     &                IRESTRICT,ICONSPA,ICONSPB,
     &                SCLFAC,
     &                IPERTOP,IH0INSPC,
     &                ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                IRELACI,IDEN12,IFINMO)
*
* Direct GAS routine for density matrices
*
* Jeppe Olsen , Winter of 1991
*
* Modified for double group calculations , T.F + J.O, May 14, 1998
* cont'd      Sept. 1999  Timo Fleig
* cont'd      Feb.  2000  Timo Fleig
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
* LC : Length of scratch array for C
* LS : Length of scratch array for S
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
*.General input
      INTEGER ICOCOC(NOCTPA_C,NOCTPB_C)
      INTEGER ICSMOS(NSMST)
      INTEGER ICBLTP(*)
      INTEGER NSSOA_S(NSMST,*), NSSOB_S(NSMST,*)
      INTEGER NSSOA_C(NSMST,*), NSSOB_C(NSMST,*)
      INTEGER SXSTSM(NSMSX,NSMST),STSTSX(NSMST,NSMST)
      INTEGER STSTDX(NSMST,NSMST), ADSXA(MXPOBS,2*MXPOBS)
      INTEGER SXDXSX(2*MXPOBS,4*MXPOBS), ASXAD(MXPOBS,2*MXPOBS)
      INTEGER NOBPTS(MXPNGAS,*),IOBPTS(MXPNGAS,*),ITSOB(*)
      INTEGER NOBPTS2(MXPNGAS,*),IOBPTS2(MXPNGAS,*)
      INTEGER NELFSPGP(MXPNGAS,*)
      INTEGER ICONSPA(NOCTPA_S,NOCTPA_S), ICONSPB(NOCTPB_S,NOCTPB_S)
*.Scratch
      DIMENSION SB(*),CB(*),C2(*)
      DIMENSION CSCR(*),SSCR(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*),
     &          XI1S(*),XI2S(*),XI3S(*),XI4S(*)
      INTEGER   LCBLOCK(*),I1CBLOCK(*),ICBLOCK(8,*),LECBLOCK(*)
*
      DIMENSION CJRES(*),SIRES(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      dimension ISTRFL(12)
      DIMENSION SCLFAC(*)
*.
      REAL * 8 INPROD
*
C?    call qenter('DBLOC')
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.20) then
        write(6,*) ' ================='
        write(6,*) ' DBLOCKS speaking :'
        write(6,*) ' ================='
        write(6,*)
        write(6,*)
     &  ' Contributions from the following L blocks( TTSS) :'
          do IBLOCK = 1,NSBLOCK,1
            if (ISBLOCK(1,IBLOCK).gt.0)
     &      write(6,'(10X,4I3,2I8)') (ISBLOCK(II,IBLOCK),II=1,4)
          end do
      end if
*
      if (NTEST.ge.100) then
        write(6,*) ' C vector '
        call wrtvcd(CB,LUC,1,-1)
      end if
* ===========================
* 1 : Arrays for accessing C
* ============================
*. Find batches of C - strings
      CALL PART_CIV2(IDC,ICBLTP,NSSOA_C,NSSOB_C,NOCTPA_C,NOCTPB_C,NSMST,
     &              LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &              I1CBLOCK,ICBLOCK,0,1)
*
*. Find the active blocks on LUC, store info in SCLFAC
*  LUC is scratch unit; defined and provided with current
*  C batch in dassg3.
C     FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
*     \/ Change to accomodate for ICBAT_RES
      CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
      if (NTEST.ge.100) then
        write(6,*) ' memchk after FIND_ACTIVE_BLOCKS'
        CALL LMEMCHK
      end if
*
* Loop over batches over C blocks
CTF     MXEXC = 1
CTF   END IF
*
      REWIND LUC
      IF(ICBAT_RES.EQ.1) THEN
        if (NTEST.ge.200) then
          WRITE(6,*) ' Restricted set of C batches '
          WRITE(6,*) ' ICBAT_INI ICBAT_END', ICBAT_INI,ICBAT_END
        end if
        JCBAT_INI = ICBAT_INI
        JCBAT_END = ICBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = NCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END
*
*. Read C blocks of this batch into core
*
        ICOFF = 1
        NJBLOCK = LCBLOCK(JCBATCH)
        DO JJCBLOCK = 1, NJBLOCK
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
          SCLFAC(JBLOCK) = 1.0d0
        END DO
*
*. Loop over blocks of sigma and C in core and obtain  contribution from
*  given C block to given S block
*. Loop over C blocks
        DO 9000 ICBLK = I1CBLOCK(JCBATCH),
     &                I1CBLOCK(JCBATCH)-1+NJBLOCK
          JATP = ICBLOCK(1,ICBLK)
          JBTP = ICBLOCK(2,ICBLK)
          JASM = ICBLOCK(3,ICBLK)
          JBSM = ICBLOCK(4,ICBLK)
          if (NTEST.ge.500) then
            write(6,*) 'JATP = ',JATP
            write(6,*) 'JBTP = ',JBTP
            write(6,*) 'JASM = ',JASM
            write(6,*) 'JBSM = ',JBSM
          end if
          ICOFF = ICBLOCK(5,ICBLK)
          NJA = NSSOA_C(JASM,JATP)
          NJB = NSSOB_C(JBSM,JBTP)
*
          IF(SCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
            CALL PRMBLK(IDC,ISTRFL,JASM,JBSM,JATP,JBTP,PS,PL,
     &                  LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
            if (NTEST.ge.100) then
              write(6,*)
              write(6,'(A10,1X,I2)') 'NPERM =   ',NPERM
            end if
*. Start with transposed block
            DO 8765 IPERM = NPERM,1,-1
              if (NTEST.ge.100) then
                if (IPERM.ne.1) then
                  write(6,*) 'Doing permuted block!'
                end if
              end if
              LLASM = LASM(IPERM)
              LLBSM = LBSM(IPERM)
              LLATP = LATP(IPERM)
              LLBTP = LBTP(IPERM)
              if (NTEST.ge.500) then
                write(6,*) 'LLATP = ',LLATP
                write(6,*) 'LLBTP = ',LLBTP
                write(6,*) 'LLASM = ',LLASM
                write(6,*) 'LLBSM = ',LLBSM
              end if
              NLLA = NSSOA_C(LLASM,LLATP)
              NLLB = NSSOB_C(LLBSM,LLBTP)
*
*. The routines assumes on input that the blocks are transposed, so,
*. Initial tour, IPERM = 1 corresponds always to no transpose, so transpose!
              IF(IPERM.EQ.1) THEN
*
*  IDC is always 1, meaning that
*  PSSIGN.EQ.0.0D0.AND.PLSIGN.EQ.0.0D0
*  This is a relativistic code, no reference to S,L,Ms,Ml
*  quantum numbers whatsoever! No fuss!
                IF(IDC.EQ.2.AND.JATP.EQ.JBTP.AND.JASM.EQ.JBSM)
     &                   THEN
*
*. Diagonal blocks, Transposing corresponds to scaling
                  IF(PS.EQ.-1.0D0) THEN
                    CALL SCALVE(CB(ICOFF),PS,NJA*NJB)
                  END IF
                ELSE
*
*. offdiagonal blocks, explicit transposing
                  CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
                  CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
                  IF(LSGN(IPERM).EQ.-1)
     &            CALL SCALVE(CB(ICOFF),-1.0D0,LROW*LCOL)
                END IF
              END IF
*
              DO 10000 ISBLK = 1, NSBLOCK
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
                    WRITE(6,*) ' gsbbdndg will be called for '
                    WRITE(6,*) ' Left block : '
                    WRITE(6,*) ' ISOFF ', ISOFF
                    WRITE(6,*) ' ISBLK IASM IBSM IATP IBTP'
                    WRITE(6,'(5I5)')  ISBLK,IASM,IBSM,IATP,IBTP
                    WRITE(6,*) ' Right block : '
                    WRITE(6,*) ' ICBLK LLASM LLBSM LLATP LLBTP'
                    WRITE(6,'(5I5)')  ICBLK,LLASM,LLBSM,LLATP,LLBTP
                    WRITE(6,*) ' ICOFF ', ICOFF
                    WRITE(6,*) ' Overall scale',SCLFAC(ICBLK)
                  END IF
*
                  XFAC = SCLFAC(ICBLK)
*
                  if (IRELACI.eq.0) then
*. Nonrelativistic density
*  Not up to date. Do not use !!!
*
                    write(6,*) 'IRELACI = 0 in dblocks!'
                    write(6,*) 'Option not implemented yet!'
                    Call Abend2( 'Quitting.' )
*
CTF                 CALL GSDNBB2(IASM,IATP,IOCTPA_S,
CTF  &                 IBSM,IBTP,IOCTPB_S,
CTF  &                 LLASM,LLATP,LLBSM,LLBTP,NGAS,
CTF  &                 NELFSPGP(1,IATP+IOCTPA_S-1),
CTF  &                 NELFSPGP(1,IBTP+IOCTPB_S-1),
CTF  &                 NELFSPGP(1,LLATP+IOCTPA_S-1),
CTF  &                 NELFSPGP(1,LLBTP+IOCTPB_S-1),
CTF  &                 NAEL_S,NBEL_S,
CTF  &                 IAGRP_S,IBGRP_S,
CTF  &                 SB(ISOFF),CB(ICOFF),IDOH2,
CTF  &                 ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
CTF  &                 NOBPTS,IOBPTS,MXPNGAS,ITSOB,MAXI,MAXK,
CTF  &                 SSCR,CSCR,I1,XI1S,I2,XI2S,XINT,C2,
CTF  &                 NSMOB,NSMST,NSMSX,NSMDX,
CTF  &                 NIA,NIB,NLLA,NLLB,MXPOBS,IDC,PS,
CTF  &                 ICJKAIB,CJRES,SIRES,I3,XI3S,I4,XI4S,
CTF  &                 MXSXBL,MXSXST,MOCAA,MOCAB,IAPR,
CTF  &                 IPRNT,IPERTOP,IPTSPC,JPTSPC,IJOP,0,IDUM,XFAC)
                  else if (IRELACI.eq.1) then
*. No-pair relativistic Hamiltonian
                    IATPABS = IATP +IOCTPA_S-1
                    IBTPABS = IBTP +IOCTPB_S-1
                    JATPABS = JATP +IOCTPA_C-1
                    JBTPABS = JBTP +IOCTPB_C-1
*
                    if (NTEST.ge.200) then
                      write(6,*) 'String type offsets in dblocks:'
                      write(6,*) 'IATPABS = ',IATPABS
                      write(6,*) 'IBTPABS = ',IBTPABS
                      write(6,*) 'JATPABS = ',JATPABS
                      write(6,*) 'JBTPABS = ',JBTPABS
                    end if
*
Call relativistic one-body density with given pair of blocks
                    call gsdnbb2dg(IDEN12,
     &                 RHO1UUR,RHO1UUI,
     &                 RHO1UBR,RHO1UBI,
     &                  RHO1BBR,RHO1BBI,
     &                  RHO1BUR,RHO1BUI,
     &                 IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &                 NGAS,NELFSPGP(1,IATPABS),NELFSPGP(1,IBTPABS),
     &                 NELFSPGP(1,JATPABS),NELFSPGP(1,JBTPABS),
     &                 IAGRP_S,IBGRP_S,IAGRP_C,IBGRP_C,
     &                 SB(ISOFF),CB(ICOFF),C2,
     &                 ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &                 MXPNGAS,NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &                 MAXI,MAXK,
     &                 SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &                 NSMOB,NSMST,NSMSX,NSMDX,
     &                 NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,NACOB2,
     &                 RHO1S,SCLFAC,IRIL,IRIR)
                  else
                    write(6,*) 'Something is wrong in DBLOCKS'
                    write(6,*) 'IRELACI corrupted. IRELACI = ',IRELACI
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
*
C?    CALL QEXIT('DBLOC')
      RETURN
      END
