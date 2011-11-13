!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl

***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
*  Arrange density matrices for diagonalization
*  Input matrices are GAS-symmetry ordered.
*
*   1) Loop over output spinor symmetries (done outside!)
*    2) Loop over real/imaginary part (done outside!)
*     3) Loop over unbarred/barred type of left index
*      4) Loop over unbarred/barred type of right index
*         Pick relevant symmetry blocks from GAS blocks
*         and arrange them according to the spinor symmetry.
*
*   Timo Fleig
*    fixed bug in UB handling (05.11.2001)
*
      subroutine arr_dens(RHO1UU,RHO1UB,RHO1BU,RHO1BB,
     &                    NGSOB,NGSOB2,
     &                    NGSOB_HE,NGSOB2_HE,BLKOUT,
     &                    READEN,IDBG,ISMDIM,
     &                    ITDIM,ISPSM,NIRR,IGRP,
     &                    ISPSF_MO,ITOTDM,NGAS,IPRDEN)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*
      dimension RHOSYM(ISMDIM,ISMDIM,NIRR)
      dimension BLKOUT(*),READEN(ISMDIM**2)
      dimension NGSOB_HE(2,NGAS),NGSOB2_HE(2,NGAS),
     &          NGSOB(MXPOBS,MXPNGAS),NGSOB2(MXPOBS,MXPNGAS)
      dimension ISPSF_MO(MXNDGIRR,2),ITOTDM(MXNSPIR)
*
      NTESTL = 20
      NTEST = max(IPRDEN,NTESTL)
*
*  Copy matrix with number of unbarred spinors and number of
*  unbarred spinors (2) per (dgsym,GAS) to help array.
      do ISYM=1,2,1
        do IGAS=1,NGAS,1
          NGSOB_HE(ISYM,IGAS) = NGSOB(ISYM+NIRR+2*(IGRP-1),IGAS)
          NGSOB2_HE(ISYM,IGAS) = NGSOB2(ISYM+NIRR+2*(IGRP-1),IGAS)
        end do
      end do
*
      if (NTEST.ge.20) then
        write(6,*) 'Checking help matrices NGSOB(2)_HE: '
        call iwrtma(NGSOB_HE,NIRR,NGAS,NIRR,NGAS)
        call iwrtma(NGSOB2_HE,NIRR,NGAS,NIRR,NGAS)
      end if
*
      if (IDBG.ge.4.and.IDBG.le.6) then
        do IUB=1,2,1
          do JUB=1,2,1
* UB = 1 means unbarred index
* UB = 2 means barred index
            if (NTEST.ge.50) write(6,*) ' Processing UB block ',IUB,JUB
            IDIMR = ITOTDM(ISPSF_MO(ISPSM,IUB))
            IDIMC = ITOTDM(ISPSF_MO(ISPSM,JUB))
            if (IUB.eq.1.and.JUB.eq.1) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB_HE,
     &                   NGSOB_HE,IUB,JUB,TRDUM,TRDUM,RHO1UU,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = 1
                ICOLOFF = 1
              else if (IEOSS.eq.2) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            else if (IUB.eq.1.and.JUB.eq.2) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB_HE,
     &                   NGSOB2_HE,IUB,JUB,TRDUM,TRDUM,RHO1UB,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              else if (IEOSS.eq.2) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            else if (IUB.eq.2.and.JUB.eq.1) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB2_HE,
     &                   NGSOB_HE,IUB,JUB,TRDUM,TRDUM,RHO1BU,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = 1
              else if (IEOSS.eq.2) then
                IROWOFF = 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            else if (IUB.eq.2.and.JUB.eq.2) then
              call pick_active2(NGAS,NIRR,NGAS,NIRR,NGSOB2_HE,
     &                   NGSOB2_HE,IUB,JUB,TRDUM,TRDUM,RHO1BB,
     &                   BLKOUT,0,0,1,1,ITDIM,IDIMR,IDIMC,2,2,
     &                   NTEST)
              call evenodd(IEOSS,ISPSM)
              if (IEOSS.eq.1) then
                IROWOFF = ITOTDM(ISPSF_MO(ISPSM,IUB)) + 1
                ICOLOFF = ITOTDM(ISPSF_MO(ISPSM,JUB)) + 1
              else if (IEOSS.eq.2) then
                IROWOFF = 1
                ICOLOFF = 1
              end if
              call vectomat(BLKOUT,READEN,IDIMR*IDIMC,ISMDIM,ISMDIM,
     &                      IDIMR,IDIMC,IROWOFF,ICOLOFF)
            end if
          end do
        end do
      else if (IDBG.ge.7.and.IDBG.le.8) then
        write(6,*) ' Double group ',IDBG
        write(6,*) ' not inplemented yet in density calculation. '
        stop 'Quitting.'
      end if
*
      if (NTEST.ge.10) then
        write(6,*)
        write(6,*) 'Checking output block READEN in arr_dens.'
        call wrtmat(READEN,ISMDIM,ISMDIM,ISMDIM,ISMDIM)
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DASSG3(RHO1UUR,RHO1UUI,
     &                  RHO1UBR,RHO1UBI,
     &                  RHO1BBR,RHO1BBI,
     &                  RHO1BUR,RHO1BUI,
     &                  CB,SB,SCR,NBATS,LBATS,LEBATS,I1BATS,IBATS,
     &                  INONREL,IDEN12,IFINMO,IPRNT)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to DBLOCK
*
*               March 1998: Relativistic version
*                           real/imaginary loops,
*                           MK2 coupling
*
*               Density version in May 1998, TF and JO
*=========================================================
* OVERVIEW OVER TRANSPOSITIONS OF S AND C BLOCKS
*=========================================================
*
* The blocks are assumed transposed at the beginning, meaning
* a global transposition is needed in sblocks/dblocks.
* General rule: For processing a beta block, a transposition is needed.
*
* dblocks:      Offdiagonal blocks, explicit transposing
*               (Always done with PSSIGN and PLSIGN = 0 => IDC = 1
*               CALL TRPMT3(CB(ICOFF),NJA,NJB,C2)
*               CALL COPVEC(C2,CB(ICOFF),NJA*NJB)
*                - - - - - - - - - - - - - - - - - - - - - - -
* rssbcbrel:      <beta || beta>      TRP(S)  TRP(C)
*                                     do matrix operations
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <alpha || alpha>    do matrix operations
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <alpha || beta>     TRP(S)  TRP(C)
*                                     TRP(S)
*                                     do matrix operations
*                                     TRP(S)
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
*                 <beta || alpha>             TRP(C)
*                                     do matrix operations
*                                             TRP(C)
*                                     TRP(S)  TRP(C)
*                - - - - - - - - - - - - - - - - - - - - - - -
* sblocks:      Transpose final sigma block back.
*
*
* =====
* Input
* =====
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
*.Scratch
      DIMENSION SB(*),CB(*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cstate.inc"
#include "clunit.inc"
#include "cands.inc"
#include "gasstr.inc"
#include "strbas.inc"
*
      CALL QENTER('DASSG')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' DASSG3 speaking :'
        WRITE(6,*) ' ================='
        write(6,*)
        WRITE(6,*) ' DASSG3 : NBATS = ',NBATS
      END IF
*
*  LUSC62: logical unit right (C vector)
*  LUSC61: logical unit left  (S vector, identical)
*
      LBLK = -1
      CALL REWINE(LUSC61,-1)
*. Loop over vectors of sigma and C : MK2 space + Real/Imag
      DO IRIL = 1,2,1
        if (NTEST.ge.50) write(6,*) 'left r/i is  ',IRIL
        if (IRIL.eq.2) then
*. SKIP EOV mark between real and imaginary part of S
          CALL IFRMDS(IONEM,1,-1,LUSC61)
        end if
*
        JBATABS = 0
        do IMK2_L = 1,NMS2VAL,1
          MK2_L = MS2VAL(IMK2_L)
          IASTP = IST_FOR_DT(1,IMK2_L)
          IBSTP = IST_FOR_DT(2,IMK2_L)
*. Loop over batches of sigma for given IRIL and IMK2_L
          DO JBATS = 1,NBAT_MS2(IMK2_L),1
            JBATABS = JBATABS + 1
            if (NTEST.ge.20) then
              WRITE(6,*)
              WRITE(6,*)
     &        ' =============================================='
              WRITE(6,*)
     &        ' Start of construction from  Left batch ', JBATABS
              WRITE(6,*)
     &        ' =============================================='
              WRITE(6,*)
            end if
*. Read left batch in
            do ISBLK = 1,LBATS(JBATABS),1
              ISATP = IBATS(1,ISBLK)
              ISBTP = IBATS(2,ISBLK)
              ISASM = IBATS(3,ISBLK)
              ISBSM = IBATS(4,ISBLK)
              ISOFF = IBATS(5,ISBLK)
              LEN2  = IBATS(7,ISBLK)
*
              print*,'ISBLK ',ISBLK
              print*,'ISATP,ISBTP,ISASM,ISBSM,ISOFF',
     &                ISATP,ISBTP,ISASM,ISBSM,ISOFF
              XDUM = 0.0D0
              ISCALE = 1
*. Signs taken from CSTATE
              CALL GSTTBLD(SB(ISOFF),ISATP,ISASM,ISBTP,ISBSM,
     &              WORK(KNSTSO(IASTP)),WORK(KNSTSO2(IBSTP)),
     &              PSSIGN,IDC,PLSIGN,LUSC61,SCR,
     &              NSMST,ISCALE,XSCALE)
              print*,'TESTING L VECTOR FRACTION:'
              call wrtmat(SB(ISOFF),1,LEN2,1,LEN2)
            end do
*           ^ End of loop over S-blocks in batch
*
            call REWINE(LUSC62,-1)
            do IRIR = 1,2,1
              if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
              if (IRIR.EQ.2) then
*. Skip EOV mark between real and imaginary part
                call IFRMDS(IONEM,1,-1,LUSC62)
              end if
              do IMK2_R = 1,NMS2VAL,1
                MK2_R = MS2VAL(IMK2_R)
                IACTP = IST_FOR_DT(1,IMK2_R)
                IBCTP = IST_FOR_DT(2,IMK2_R)
*. Copy blocks with given MK2 to LUSC37
*  Explicit selection and access done in dblocks (analogy to sblocks)
                CALL REWINE(LUSC37,-1)
                NBL_C = NBLK_MS2(IMK2_R)
                CALL COPNBLKD(LUSC62,LUSC37,CB,NBL_C,0,LBLK)
C                    COPNBLKD(LUIN,LUOUT,SEGMNT,NBLOCK,IREW,LBLK)
*. End of vector mark
                CALL ITODS(-1,1,-1,LUSC37)
                if (NTEST.ge.20) then
                  WRITE(6,*)
     &            ' ********************************************* '
                  WRITE(6,*)
     &            ' *                                           * '
                  WRITE(6,*)
     &            ' * Start of contribution from given CISPACE  * '
                  WRITE(6,*)
     &            ' *                                           * '
                  WRITE(6,*)
     &            ' ********************************************* '
                  WRITE(6,*)
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of S :  MK2_L, IRIL =', MK2_L,IRIL
                  WRITE(6,'(4X,A,2I3)')
     &            ' MK2/IRI vector of C :  MK2_R, IRIR =', MK2_R,IRIR
                  WRITE(6,*)
                end if
*
C               CALL WRTVCD(CB,LUSC37,1,-1)
                IF(ABS(MK2_L-MK2_R).LE.2) THEN
*.    Connection
CTF  What about 2-particle densities?
CTF  Shouldn't ABS be .le. 4 instead?
*
*. Obtain contribution to batch sigma from all blocks with MK2_R
                  CALL DBLOCK(INONREL,IDEN12,IRIL,IRIR,
     &                        RHO1UUR,RHO1UUI,
     &                        RHO1UBR,RHO1UBI,
     &                        RHO1BBR,RHO1BBI,
     &                        RHO1BUR,RHO1BUI,
     &                        LBATS(JBATABS),IBATS(1,I1BATS(JBATABS)),
     &                        1,CB,SB,LUSC37,0,0,0,0)
*
                end if
*               ^End if interacting
              end do
*             ^End of loop over IMK2_R
            end do
*           ^End of loop over real/imaginary part of C
          end do
*         ^End of batches for given MK2_L
        end do
*       ^End of loop over IMK2_L
*. End of real/imag part of sigma => end of vector
      end do
*     ^End of loop over real/imag part of sigma
*
      call QEXIT('DASSG')
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine dbg_rho(C,HC,INONREL,IDEN12,JROOT)
*
* Outer routine for generation of density matrix
* for double group CI.
*
* The two vectors are supposed to reside on LUSC62 (C) and LUSC61 (S)
* C and HC are work space partitions called WORK(KVEC1) and
* WORK(KVEC2), respectively, above (in gasci).
*
* Generalized to allow for relativistic densities, T.F. + J.O., may 98
* cont'd ... fall 1998, T.F.
* cont'd ... august 1999, Timo Fleig
* cont'd ... december 1999, Timo Fleig
* cont'd ... june 2000, Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma
#include "cands.inc"
*
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
#include "symm.inc"
#include "clunit.inc"
#include "files_r.inc"
*
      common/glbref/ISPC,ISM,ITERGLB
*
      character*3 FILESTAT
*
      dimension EIGVDEN(MXNSPIPDMGRP)
      dimension NGSOB_HE(MXPIRR*MXPNGAS),
     &          NGSOB2_HE(MXPIRR*MXPNGAS)
*
      NTESTL = 100
      NTEST = max(NTESTL,IPRDEN)
*
      CALL QENTER('DBG_RH')
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'DBG_RH')
*
*. Arrays for partitioning of sigma
* (MS2 arrays are in CSTATE )
      CALL Z_BLKFO_REL(IDC,NMS2VAL,
     &                 ISPC,ISM,KLSLBT,KLSLEBT,KLSI1BT,KLSIBT,
     &                 KLSBLTP,NBATCH,NBLOCK,NBLK_MS2,IBLK_MS2,
     &                 NBAT_MS2,IBAT_MS2,IDUMMY,.FALSE.,0,0,IDUMMY)
*
*  Erase density matrices
      do IC=1,2,1
        do IB=1,2,1
          do IA=1,2,1
            do IER=1,NTOOB**2,1
              WORK(KRHO1_DBG(IA,IB,IC,JROOT)+IER-1) = 0.D0
            end do
          end do
        end do
      end do
*
      CALL DASSG3(WORK(KRHO1_DBG(2,2,1,JROOT)),
     &            WORK(KRHO1_DBG(2,2,2,JROOT)),
     &            WORK(KRHO1_DBG(2,1,1,JROOT)),
     &            WORK(KRHO1_DBG(2,1,2,JROOT)),
     &            WORK(KRHO1_DBG(1,1,1,JROOT)),
     &            WORK(KRHO1_DBG(1,1,2,JROOT)),
     &            WORK(KRHO1_DBG(1,2,1,JROOT)),
     &            WORK(KRHO1_DBG(1,2,2,JROOT)),
     &            C,HC,WORK(KVEC3),NBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &            WORK(KLSI1BT),WORK(KLSIBT),
     &            INONREL,IDEN12,IFINMO,IPRDEN)
*
*    WORK(KRHO1_DBG(A,B,C,D)), A=1,2  first  index barred,unbarred
*                              B=1,2  second index barred, unbarred
*                              C=1,2  real, imaginary part of density
*                              D=JROOT density for root JROOT
*
*  Check determined density matrices
      if (NTEST.ge.00) then
       write(6,*) 'Checking calculated density matrices for this root:'
       do I=1,2,1
c       I = 2
        do J=1,2,1
         do K=1,2,1
          write(6,*)
          write(6,*) 'Checking address on WORK array:'
          write(6,'(A10,I1,A1,I1,A1,I1,A1,I2,A2,I9)')
     &            'KRHO1_DBG(',I,',',J,',',K,',',JROOT,')=',
     &             KRHO1_DBG(I,J,K,JROOT)
          write(6,*)
          if (I.eq.1.and.J.eq.1.and.K.eq.1) write(6,*) 'BB den, real'
          if (I.eq.1.and.J.eq.1.and.K.eq.2) write(6,*) 'BB den, imag'
          if (I.eq.1.and.J.eq.2.and.K.eq.1) write(6,*) 'BU den, real'
          if (I.eq.1.and.J.eq.2.and.K.eq.2) write(6,*) 'BU den, imag'
          if (I.eq.2.and.J.eq.1.and.K.eq.1) write(6,*) 'UB den, real'
          if (I.eq.2.and.J.eq.1.and.K.eq.2) write(6,*) 'UB den, imag'
          if (I.eq.2.and.J.eq.2.and.K.eq.1) write(6,*) 'UU den, real'
          if (I.eq.2.and.J.eq.2.and.K.eq.2) write(6,*) 'UU den, imag'
          call wrtmat(WORK(KRHO1_DBG(I,J,K,JROOT)),
     &                NTOOB,NTOOB,NTOOB,NTOOB)
          write(6,*)
         end do
        end do
       end do
      end if
*
* Arrange density matrices for diagonalization
* and do the job. Write natural spinors to files.
*
      if (JROOT.eq.1) then
*  Currently we determine natural spinors only for the
*  lowest root. On the long run, the spinors should be
*  averaged over a selection of states (components).
        if (INATIT.eq.1) then
          FILESTAT = 'NEW'
        else if (INATIT.ge.2) then
          FILESTAT = 'OLD'
        end if
        open(unit=LUMOUTR,file=NLUMOUTR_MPI(1:LUMOUTRI_MPIL),
     &       status='UNKNOWN',form='FORMATTED')
        open(unit=LUMOUTI,file=NLUMOUTI_MPI(1:LUMOUTRI_MPIL),
     &       status='UNKNOWN',form='FORMATTED')
        call rewine(LUMOUTR,-1)
        call rewine(LUMOUTI,-1)
        ISUM = 0
        ISTORE = 0
        do ISYM=1,NIRR_PN,1
          do IGS=1,NGAS,1
            ISUM = ISUM + NGSSH(ISYM,IGS)
          end do
          if (ISUM.gt.ISTORE) ISTORE = ISUM
          ISUM = 0
        end do
        NFHELEN = ISTORE**2
        ISTODIM = (2*ISTORE)**2
c       ISTOTRI = ((2*ISTORE)**2 + 2*ISTORE) / 2
        ISTOTRI = (ISTORE**2 + ISTORE) / 2
      call memmar(KNBLK,NFHELEN,'ADDL',2,'TRHEMA')
      call memmar(KREADEN,ISTODIM*2,'ADDL',2,'READEN')
      call memmar(KREATRI,ISTOTRI*2,'ADDL',2,'REATRI')
        if (NTEST.ge.20) then
          write(6,'(A24,I8,A4)')
     &            'Alloc. full help matrix ',NFHELEN,' R*8'
          write(6,'(A24,I8,A4)')
     &            'Alloc. storing dimension',ISTODIM,' R*8'
          write(6,'(A24,I8,A4)')
     &            'Alloc. storing triangle ',ISTOTRI,' R*8'
        end if
*
        do ISYM=1,NIRR_PN,1
          ISPSM = ISYM + NIRR_PN
          if (NTEST.ge.10) then
            write(6,*) 
            write(6,*) ' ========================================'
            write(6,*) '  Determining natural spinors, symm.',ISPSM
            write(6,*) ' ========================================'
          end if
          if (ISYM.le.2) then
            IGRP = 1
            if (DOUGRP.eq.7.and.ISYM.eq.2) IGRP = 2
          else if (ISYM.ge.3) then
            IGRP = 2
          end if
          ITRIALEN = (SYMDIM(IGRP)**2 + SYMDIM(IGRP)) / 2
          ISCRLEN = (2*SYMDIM(IGRP)+1)*SYMDIM(IGRP)
      call memmar(KNATSPI,(SYMDIM(IGRP)*2)**2,
     &                'ADDL',2,'NATSPI')
      call memmar(KNSPIR,SYMDIM(IGRP)**2,'ADDL',2,'NSPINR')
      call memmar(KNSPII,SYMDIM(IGRP)**2,'ADDL',2,'NSPINI')
      call memmar(KSCRATCH,ISCRLEN,'ADDL',2,'SCRALE')
          if (NTEST.ge.20) then
            write(6,'(A24,I8,A4)')
     &           'Alloc. nat. spinor array',SYMDIM(IGRP)*2,' R*8'
            write(6,'(A24,I8,A4)')
     &           'Alloc. real NS array    ',SYMDIM(IGRP)**2,' R*8'
            write(6,'(A24,I8,A4)')
     &           'Alloc. imag. NS array   ',SYMDIM(IGRP)**2,' R*8'
          end if
*
          do IRI=1,2,1
            call arr_dens(WORK(KRHO1_DBG(2,2,IRI,JROOT)),
     &                    WORK(KRHO1_DBG(2,1,IRI,JROOT)),
     &                    WORK(KRHO1_DBG(1,2,IRI,JROOT)),
     &                    WORK(KRHO1_DBG(1,1,IRI,JROOT)),
     &                    NGSOB,NGSOB2,
     &                    NGSOB_HE,NGSOB2_HE,
     &                    WORK(KNBLK),
     &                    WORK(KREADEN+(IRI-1)*ISTODIM),
     &                    DOUGRP,SYMDIM(IGRP),
     &                    NTOOB,ISPSM,NIRR_PN,IGRP,
     &                    ISPSF_MO,ITOTDM,NGAS,IPRDEN)
            call fold(WORK(KREATRI+(IRI-1)*ISTOTRI),
     &                WORK(KREADEN+(IRI-1)*ISTODIM),
     &                SYMDIM(IGRP),1,ITRIALEN)
          end do
*
*  Diagonalize one-particle density matrix
*
          if (NTEST.ge.100) then
            write(6,*)
            write(6,*) 'Folded real density matrix:'
            call wrtmat(WORK(KREATRI),ISTOTRI,1,ISTOTRI,1)
            write(6,*)
            write(6,*) 'Folded imaginary density matrix:'
            call wrtmat(WORK(KREATRI+ISTOTRI),ISTOTRI,1,
     &                  ISTOTRI,1)
          end if
          call ceigen(WORK(KREATRI),WORK(KREATRI+ISTOTRI),
     &                SYMDIM(IGRP),EIGVDEN,WORK(KNATSPI),
     &                WORK(KSCRATCH))
*
          if (NTEST.ge.100) then
            write(6,*) '--------------------------------'
            write(6,*) 'Writing output natural spinors.'
            write(6,*) '--------------------------------'
            call wrtmat(WORK(KNATSPI),(SYMDIM(IGRP)*2)**2,1,
     &                   (SYMDIM(IGRP)*2)**2,1)
            write(6,*) '--------------------------------'
            write(6,*) '--------------------------------'
          end if
*
          if (NTEST.ge.50) then
            write(6,*) 'Natural spinor occupation numbers'
            write(6,*) 'Spinor symmetry: ',ISPSM
            call wrtmat(EIGVDEN,SYMDIM(IGRP),1,SYMDIM(IGRP),1)
          end if
*
*  Write this symmetry to file LUMOUT
*
          call wrt_natspi(WORK(KNATSPI),SYMDIM(IGRP),
     &                    LUMOUTR,LUMOUTI,ISPSM)
*
          write(6,*) 'Checking memory in SPISYM loop ...'
          call rmemchk
*
        end do
*       ^ End loop over spinor symmetries
*
        close (unit=LUMOUTR)
        close (unit=LUMOUTI)
*
      end if
*     ^ End if which root is being treated (to be modified)
*
*. Eliminate local memory
      call memmar(KDUM ,IDUM,'FLUSM ',2,'DBG_RH')
*
      CALL QEXIT('DBG_RH')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine dblock(INONREL,IDEN12,IRIL,IRIR,
     &                  RHO1UUR,RHO1UUI,
     &                  RHO1UBR,RHO1UBI,
     &                     RHO1BBR,RHO1BBI,
     &                     RHO1BUR,RHO1BUI,
     &                  NBLOCK,IBLOCK,IBOFF,CB,HCB,LUSC37,
     &                  IRESTRICT,
     &                  ICBAT_RES,ICBAT_INI,ICBAT_END)
*
* Calculate contribution to density matrices from a
* batch of left hand side blocks, all right hand side blocks.
* The left hand side blocks in action
* are the NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES.eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
* Generalized to NOPAIR Hamiltonian, T.F + J.O., March21 1998, 7.22
*                                    (not bad for a saturday)
*  cont'd Sept. 1999   Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
* IA/BS/CTP : Type (number of elecs in A/B,S/C)
#include "cands.inc"
*. Sigma blocks require
      integer IBLOCK(8,*)
*
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Used : NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
#include "cintfo.inc"
*
      integer ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      common/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      common/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'DBLOCK')
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRDEN)
*
      if (NTEST.ge.100) then
        write(6,*) ' DBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
        write(6,*) ' LUSC37 in DBLOCK ', LUSC37
      end if
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = IASTP
      IBTP = IBSTP
*
      JATP = IACTP
      JBTP = IBCTP
*. Number of supergroups
      NOCTPA_S = NOCTYP(IATP)
      NOCTPB_S = NOCTYP(IBTP)
*
      NOCTPA_C = NOCTYP(JATP)
      NOCTPB_C = NOCTYP(JBTP)
*. Offset for supergroups
      IOCTPA_S = IBSPGPFTP(IATP)
      IOCTPB_S = IBSPGPFTP(IBTP)
*
      IOCTPA_C = IBSPGPFTP(JATP)
      IOCTPB_C = IBSPGPFTP(JBTP)
*
      NAEL_S = NELEC(IATP)
      NBEL_S = NELEC(IBTP)
*
      NAEL_C = NELEC(JATP)
      NBEL_C = NELEC(JBTP)
*
*. connection matrices for supergroups, just st for for s-space
*  (only used when number of elecs are identical )
*
      call memmar(KCONSPA_S,NOCTPA_S**2,'ADDL  ',1,'CNSPAC')
      call memmar(KCONSPB_S,NOCTPB_S**2,'ADDL  ',1,'CNSPBC')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
*
      call SPGRPCON(IOCTPA_S,NOCTPA_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA_S),IPRCIX)
      call SPGRPCON(IOCTPB_S,NOCTPB_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB_S),IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      call memmar(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      call memmar(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      call STSTSM_REL(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MXSTBL = MXNSTR
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
      MAXIK = MAX(MAXI,MAXK)
* Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      do IOBTP = 1,NGAS,1
        do IOBSM = 1,NSMOB,1
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        end do
      end do
      MAXIJ = MXTSOB_L ** 2
*.Local scratch arrays for blocks of C and sigma
      if (ICISTR.le.2) then
        L0BLOCK = MXSB
      else if (ICISTR.eq.3) then
        L0BLOCK = MXSOOB
      end if
      LSCR1 = MAX(L0BLOCK,L2BLOCK)
Cold      LSCR1 = max(LSCR1,LCSBLK)
      if (IPRCIX.ge.2)
     &write(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*. Arrays giving allowed type combinations
      call memmar(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',2,'CIOIO ')
      call memmar(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
*. sigma needed for MXRESC
      call iaibcm_rel(ISSPC,IATP,IBTP,WORK(KSIOIO))
      call iaibcm_rel(ICSPC,JATP,JBTP,WORK(KCIOIO))
*. Arrays giving block type
      call memmar(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*.scratch space for projected matrices and a CI block
      LSCR2 = LCJRES
      if (IPRCIX.ge.2)
     & WRITE(6,*) ' Space for resolution matrices ',LSCR2
C  memory for blocks already allocated , so
      KC2 = KVEC3
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*  Scratch space for density multiply
      LRHO1 = NTOOB ** 2
      call memmar(KRHO1S,LRHO1,'ADDL  ',2,'RHO1S ')
CTF
*. Allocation of memory for annihilation/creation from supergroups
CTF setting this undefined variable to two:
      MXSXST = 2
      MXADKBLK = 2
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0,
     &            MXSTBL*MXTSOB_L)
      call memmar(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      call memmar(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      call memmar(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      call memmar(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      call memmar(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      call memmar(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      call memmar(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      call memmar(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
*. Type of blocks (all types set to 1 for now )
      call ZBLTP_REL(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
*.Some TTS arrays
      NTTS = MXNTTS
*. for partitioning of vector
      call memmar(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      call memmar(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      call memmar(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      call memmar(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      call memmar(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
*. Space for four blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (max(NAEL_S,NAEL_C,NBEL_S,NBEL_C)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (max(NAEL_S,NBEL_S,NAEL_C,NBEL_C)) * NOCOB
      do I1234=1,4,1
      call memmar(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      call memmar(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
      call memmar(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      end do
      call memmar(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
      call dblocks(NACOB,NACOB2,IRIL,IRIR,
     &            RHO1UUR,RHO1UUI,
     &            RHO1UBR,RHO1UBI,
     &                  RHO1BBR,RHO1BBI,
     &                  RHO1BUR,RHO1BUI,
     &            NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
     &            WORK(KRHO1S),
     &            WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            WORK(KNSTSO(JATP)),WORK(KNSTSO2(JBTP)),
     &            NAEL_S,IATP,NBEL_S,IBTP,
     &            NAEL_C,JATP,NBEL_C,JBTP,
     &            IOCTPA_S,IOCTPB_S,NOCTPA_S,NOCTPB_S,
     &            IOCTPA_C,IOCTPB_C,NOCTPA_C,NOCTPB_C,
     &            NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,
     &            NOBPTS2,IOBPTS2,MXPNGAS,
     &            ITSOB,MAXIJ,MAXK,MAXI,LSCR1,
     &            LSCR1,WORK(KCSCR),WORK(KSSCR),
     &            SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &            ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &            WORK(KI1),WORK(KXI1S),
     &            WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,
     &            PSSIGN,IPRDEN,LUSC37,ICJKAIB,WORK(KCJRES),
     &            WORK(KSIRES),WORK(KI3),WORK(KXI3S),
     &            WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
     &            MOCAA,MOCAB,IAPR,
     &            WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &            IRESTRICT,WORK(KCONSPA_S),WORK(KCONSPB_S),
     &            WORK(KLSCLFAC),
     &            IPERTOP,IH0INSPC,
     &            ICBAT_RES,ICBAT_INI,ICBAT_END,
     &            IRELACI,IDEN12,IFINMO)
*
*. Eliminate local memory
      IDUM = 0
      call memmar(KDUM,IDUM,'FLUSM ',2,'DBLOCK')
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
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
#include "ipoist8.inc"
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
      ITTSS_ORD = 1
      CALL PART_CIV2(IDC,ICBLTP,NSSOA_C,NSSOB_C,NOCTPA_C,NOCTPB_C,NSMST,
     &              LC,ICOCOC,ICSMOS,NCBATCH,LCBLOCK,LECBLOCK,
     &              I1CBLOCK,ICBLOCK,0,ITTSS_ORD)
*
*. Find the active blocks on LUC, store info in SCLFAC
*  LUC is scratch unit; defined and provided with current
*  C batch in dassg3.
C     FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
*     \/ Change to accomodate for ICBAT_RES
      CALL FIND_ACTIVE_BLOCKS(LUC,-1,SCLFAC,CB)
      if (NTEST.ge.100) then
        write(6,*) ' memchk after FIND_ACTIVE_BLOCKS'
        CALL RMEMCHK
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
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE DENSI2_REL(I12,RHO1,
CTF  &                      RHO2,
     &                      L,R,LUL,LUR)
*
* Density matrices between L and R
*
* I12 = 1 => only one-body-density
* I12 = 2 => one- and two-body-density matrices
*
* Jeppe Olsen,      Oct 94
* GAS modifications Aug 95
* Two body density added, '96
*
* Table-Block driven, June 97
*
* Modified for relativistic double group CI
*   Timo Fleig, April 1998
*
* Two-body density is stored as rho2(ijkl)=<l!e(ij)e(kl)-delta(jk)e(il)!r>
* ijkl = ij*(ij-1)/2+kl, ij.ge.kl
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
* =====
*.Input
* =====
*
*.Definition of L and R is picked up from CANDS
* with L being S and  R being C
*
*. Definition of c and sigma spaces
*  IA/BS/CTP : Type (number of elecs in A/B,S/C)
#include "cands.inc"
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "spinfo.inc"
#include "glbbas.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*./LUCINP/
*. Used : NSMOB
#include "lucinp.inc"
#include "clunit.inc"
*. Scratch for string information
      COMMON/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. Specific input
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*)
CTF   dimension RHO2(*)
*. Before I forget it :
      CALL QENTER('DENSI')
      IDUM = 0
      call memmar(KDUM,IDUM,'MARK ',IDUM,'DENSI ')
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO,NACOB ** 2)
      IF(I12.EQ.2)
     &CALL SETVEC(RHO2,ZERO,NACOB ** 2 *(NACOB**2+1)/2)
*
      WRITE(6,*) ' ISSPC ICSPC in DENSI2 ',ISSPC,ICSPC
      WRITE(6,*) ' IASTP IBSTP in DENSI2 ',IASTP,IBSTP
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      call memmar(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      call memmar(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM_REL(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*
* Info for this internal space
*. type of alpha and beta strings
* In principle same loops as in RASSG3 ...
c     do IRIL = 1,-1,-2
c        do IMK2_L = 1,NMS2VAL,1
c           MK2_L = MS2VAL(IMK2_L)
c           IALTP = IST_FOR_DT(1,IMK2_L)
c           IBLTP = IST_FOR_DT(2,IMK2_L)
c              do IRIR = 1,-1,-2
c                 do IMK2_R = 1,NMS2VAL,1
c                    etc.
c
*
*

            IATP = 1
            IBTP = 2
*. alpha and beta strings with an electron removed
            IATPM1 = 3
            IBTPM1 = 4
*. alpha and beta strings with two electrons removed
            IATPM2 = 5
            IBTPM2 = 6
*
            JATP = 1
            JBTP = 2
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offsets for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*
      ILSM = ISSM
      IRSM = ICSM
*. connection matrices for supergroups
      call memmar(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      call memmar(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX_REL(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX_REL(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX_REL(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX_REL(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX_REL(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX_REL(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRDEN.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
*. replace with MXINKA !!!
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
      MAXIJ = MXTSOB_L ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRDEN.GE.2) write(6,*) ' DENSI2 : MXSB MXTSOB_L MXSOOB ',
     &       MXSB,MXTSOB_L,MXSOOB
      IF(ICISTR.LE.2) THEN
        L0BLOCK = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        L0BLOCK = MXSOOB
      END IF
      LSCR1 = MAX(L0BLOCK,L2BLOCK)
Cold      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      IF(ICISTR.EQ.1) THEN
      call memmar(KCB,LSCR1,'ADDL  ',2,'KCB   ')
      call memmar(KSB,LSCR1,'ADDL  ',2,'KSB   ')
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      call memmar(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Arrays giving allowed type combinations '
      call memmar(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      call memmar(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &               IGSOCCX,
     &               NOCTPA,NOCTPB,
     &               ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &               NELFGP,MXPNGAS,NGAS,WORK(KSIOIO),IPRDIA)
      CALL IAIBCM_GAS(LCMBSPC(ISSPC),ICMBSPC(1,ISSPC),
     &               IGSOCCX,
     &               NOCTPA,NOCTPB,
     &               ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &               NELFGP,MXPNGAS,NGAS,WORK(KCIOIO),IPRDIA)
*. Scratch space for CJKAIB resolution matrices
      CALL MXRESC_REL(WORK(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                NSMST,NSTFSMSPGP,NSTFSMSPGP2,MXPNSMST,
     &                NSMOB,MXPNGAS,NGAS,NOBPTS,NOBPTS2,IPRCIX,MAXK,
     &                NELFSPGP,
     &                MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
      IF(IPRDEN.GE.2) THEN
        WRITE(6,*) ' DENSI12 :  : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                     MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
      LSCR12 = MAX(LSCR1,2*LSCR2)
*. It is assumed that the third block already has been allocated, so
      KC2 = KVEC3
COLD  CALL MEMMAR(KC2,LSCR12,'ADDL  ',2,'KC2   ')
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      call memmar(KI1,  LSCR3       ,'ADDL  ',1,'I1    ')
      call memmar(KI2,  LSCR3       ,'ADDL  ',1,'I2    ')
      call memmar(KI3,  LSCR3       ,'ADDL  ',1,'I3    ')
      call memmar(KI4,  LSCR3       ,'ADDL  ',1,'I4    ')
      call memmar(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
      call memmar(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
      call memmar(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ')
      call memmar(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ')
*. Arrays giving block type
      call memmar(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      call memmar(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        call memmar(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST_REL(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
      CALL ZBLTP_REL(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP))
      CALL ZBLTP_REL(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      call memmar(KRHO1S,NACOB ** 2,'ADDL  ',2,'RHO1S ')
*. For natural orbitals
      call memmar(KRHO1P,NACOB*(NACOB+1)/2,'ADDL  ',2,'RHO1P ')
      call memmar(KXNATO,NACOB **2,'ADDL  ',2,'XNATO ')
*. Natural orbitals in symmetry blocks
      call memmar(KRHO1SM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      call memmar(KXNATSM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      call memmar(KOCCSM,NACOB ,'ADDL  ',2,'RHO1S ')
*
*. Space for two  blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      call memmar(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*. Space for two blocks of string occupations and four arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      call memmar(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      DO K12 = 1, 2
      call memmar(KLOCSTR(K12),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO
      DO I1234 = 1, 4
      call memmar(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
      call memmar(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
*. Arrays for partitioning of Left vector = sigma
      NTTS = MXNTTS
      call memmar(KLLBTL ,NTTS  ,'ADDL  ',1,'LBT_L  ')
      call memmar(KLLEBTL,NTTS  ,'ADDL  ',1,'LEBT_L ')
      call memmar(KLI1BTL,NTTS  ,'ADDL  ',1,'I1BT_L ')
      call memmar(KLIBTL ,8*NTTS,'ADDL  ',1,'IBT_L  ')
      call memmar(KLSCLFCL,NTTS, 'ADDL  ',2,'SCLF_L')
      ITTSS_ORD = 1
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),
     &     WORK(KLI1BTL),WORK(KLIBTL),0,1)
*. Number of BLOCKS
        NBLOCKL = IFRMR_REL(WORK(KLI1BTL),1,NBATCHL)
     &          + IFRMR_REL(WORK(KLLBTL),1,NBATCHL) - 1
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      call memmar(KLLBTR ,NTTS  ,'ADDL  ',1,'LBT_R  ')
      call memmar(KLLEBTR,NTTS  ,'ADDL  ',1,'LEBT_R ')
      call memmar(KLI1BTR,NTTS  ,'ADDL  ',1,'I1BT_R ')
      call memmar(KLIBTR ,8*NTTS,'ADDL  ',1,'IBT_R  ')
      call memmar(KLSCLFCR,NTTS, 'ADDL  ',2,'SCLF_R')
      ITTSS_ORD = 1
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),
     &     NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),
     &     WORK(KLI1BTR),WORK(KLIBTR),0,1)
*. Number of BLOCKS
        NBLOCKR = IFRMR_REL(WORK(KLI1BTR),1,NBATCHR)
     &          + IFRMR_REL(WORK(KLLBTR),1,NBATCHR) - 1
        WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR
      WRITE(6,*) ' Memcheck in densi2 before GASDN2'
      CALL RMEMCHK

      IF(ICISTR.EQ.1) THEN
         WRITE(6,*) ' Sorry, ICISTR = 1 is out of fashion'
         WRITE(6,*) ' Switch to ICISTR = 2 - or reprogram '
         Call Abend2( ' DENSI2T : ICISTR = 1 in use ' )
      ELSE IF(ICISTR.GE.2) THEN
C       CALL GASDN2(I12,RHO1,
CTF  &       RHO2,
C    &       L,R,L,R,WORK(KC2),
C    &       WORK(KCIOIO),WORK(KSIOIO),ISMOST(1,ICSM),
C    &       ISMOST(1,ISSM),WORK(KCBLTP),WORK(KSBLTP),NACOB,
C    &       WORK(KNSTSO(IATP)),WORK(KISTSO(IATP)),
C    &       WORK(KNSTSO(IBTP)),WORK(KISTSO(IBTP)),
C    &       NAEL,IATP,NBEL,IBTP,IOCTPA,IOCTPB,NOCTPA,NOCTPB,
C    &       NSMST,NSMOB,NSMSX,NSMDX,MXPNGAS,NOBPTS,IOBPTS,
C    &       MAXK,MAXI,LSCR1,LSCR1,WORK(KCSCR),WORK(KSSCR),
C    &       SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
C    &       ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
C    &       WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
C    &       WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),WORK(KINSCR),
C    &       MXPOBS,IPRDEN,WORK(KRHO1S),LUL,LUR,
C    &       PSSIGN,PSSIGN,WORK(KRHO1P),WORK(KXNATO),
C    &       NBATCHL,WORK(KLLBTL),WORK(KLEBTL),WORK(KLI1BTL),
C    &       WORK(KLIBTL),
C    &       NBATCHR,WORK(KLLBTR),WORK(KLEBTR),WORK(KLI1BTR),
C    &       WORK(KLIBTR),WORK(KCONSPA),WORK(KCONSPB),
C    &       WORK(KLSCLFCL),WORK(KLSCLFCR))
C                        KLLBTR  KLLEBTR KLI1BTR KLIBTR
      END IF
      WRITE(6,*) ' Memcheck in densi2 after GASDN2'
      CALL RMEMCHK
*
* Natural Orbitals etc.
*
C     CALL NATORB_LCR(RHO1,NSMOB,NTOOBS,NACOBS,NINOBS,
C    &                IREOST,WORK(KXNATO),
C    &                WORK(KRHO1SM),WORK(KOCCSM),NACOB,
C    &                WORK(KRHO1P),IPRDEN)
      IF(IPRDEN.GE.5.AND.I12.EQ.2) THEN
        WRITE(6,*) ' Two-electron density '
        CALL PRSYM(RHO2,NACOB**2)
      END IF
*. Eliminate local memory
      call memmar(KDUM,IDUM,'FLUSM',IDUM,'DENSI ')
      CALL QEXIT('DENSI')
C     WRITE(6,*) ' Leaving DENSI '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GSBBD1DG(RHO1,NACOB,ISCSM,ISCTP,
     &     ICCSM,ICCTP,IGRP,NROW,
     &     NGAS,ISEL,ICEL,
     &     SB,CB,
     &     NOBPTS,IOBPTS,MAXI,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,
     &     NSMOB,RHO1S,
     &     SCALEF,SCLFAC,IPRNT)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
*
*  modified for relativistic double group calculation
*  Real/imaginary logic taken care of outside
*  cleaned and shaved ...
*  summer of '97 and april '98, Timo Fleig
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be calculated/updated
* NACOB : Number of active orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IGRP : String group of columns
* NROW : Number of rows in S and C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* CB   : Input C block
* SB   : Input S block
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of spinors per type and symmetry (this set)
* IOBPTS : base for spinors of given type and symmetry (this set)
* NSMOB : Number of symmetries of spinors
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1 : Calculated/updated density block
*        This is real or imaginary part of either the
*        barred-barred or the unbarred-unbarred partition
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXIJ*MAXI*MAXK, where MAXIJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Scratch space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*. General input
      INTEGER NOBPTS(MXPNGAS,*), IOBPTS(MXPNGAS,*)
*.Input
      INTEGER ISEL(NGAS),ICEL(NGAS)
*
* Double group symmetry info
#include "symm.inc"
*
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
*
*.Local arrays
      DIMENSION ITP(MXPNGAS),JTP(MXPNGAS)
*
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' ==================='
         WRITE(6,*) ' GSBBD1DG in action '
         WRITE(6,*) ' ==================='
         WRITE(6,*)
         WRITE(6,*) ' Occupation of active left strings '
         CALL IWRTMA(ISEL,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active Right strings '
         CALL IWRTMA(ICEL,1,NGAS,1,NGAS)
*
         WRITE(6,*) ' GSBBD1DG, sclfac ',SCLFAC
         WRITE(6,*) ' GSBBD1DG, scalef for r/i: ',SCALEF
      END IF
*
      IFRST = 1
      JFRST = 1
*
* Type of single excitations that connects the two column strings
      CALL SXTYP_GAS_REL(NSXTP,ITP,JTP,NGAS,MXPNGAS,ISEL,ICEL)
*.Symmetry of single excitation that connects IBSM and JBSM
      IJSM = IDBGMULT(ISCSM,IADJSYM(ICCSM))
      IF(NTEST.GE.1000)
     &     WRITE(6,*) ' ISCSM,ICCSM,IJSM ', ISCSM,ICCSM,IJSM
      IF(IJSM.EQ.0) GOTO 1001
      DO 900 IJTP = 1,NSXTP,1
         ITYP = ITP(IJTP)
         JTYP = JTP(IJTP)
         IF(NTEST.GE.1000) write(6,*) ' ITYP JTYP ', ITYP,JTYP
         DO 800 ISM = 1,NSMOB,1
*. new i and j so new intermediate strings
            KFRST = 1
*
            JSM = IDBGMULT(ISM,IADJSYM(IJSM))
            IF(JSM.EQ.0) GOTO 800
            IF(NTEST.GE.1000) write(6,*) ' ISM JSM ', ISM,JSM
            NIORB = NOBPTS(ITYP,ISM)
            NJORB = NOBPTS(JTYP,JSM)
            IBIORB = IOBPTS(ITYP,ISM)
            IBJORB = IOBPTS(JTYP,JSM)
            if (NTEST.ge.1000) then
               write(6,*) ' NIORB NJORB ', NIORB,NJORB
               write(6,*) ' offset left orb ', IBIORB
               write(6,*) ' offset right orb ', IBJORB
            end if
            IF(NIORB.EQ.0.OR.NJORB.EQ.0) GOTO 800
*
COLD. Loop over partitionings of the row strings
            NIPART = NROW/MAXI
            IF(NIPART*MAXI.NE.NROW) NIPART = NIPART + 1
*. Loop over partitionings of N-1 strings (K strings)
            KBOT = 1-MAXK
            KTOP = 0
 700        CONTINUE
            KBOT = KBOT + MAXK
            KTOP = KTOP + MAXK
*. Single excitation information independent of I strings
*
* Determine mapping <Jb|j_+|Kb>
            J12 = 1
            K12 = 1
* What are these good for? (relevant for info on common block HIDSCR)
            IF(JFRST.EQ.1) KFRST=1
            CALL ADST_GAS(IBJORB,NJORB,ICCTP,ICCSM,IGRP,KBOT,KTOP,
     &           I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,J12,K12,
     &           KKJACT,SCLFAC)
            JFRST = 0
            KFRST = 0
*. Determine mapping <Ib|i_+|Kb>
            I12 = 2
            IF(IFRST.EQ.1) KFRST = 1
            CALL ADST_GAS(IBIORB,NIORB,ISCTP,ISCSM,IGRP,KBOT,KTOP,
     &           I2,XI2S,MAXK,NKBTC,KEND,IFRST,KFRST,I12,K12,
     &           KKIACT,SCLFAC)
            IFRST = 0
            KFRST = 0
*. Appropriate place to start partitioning over I strings
*. Loop over partitionings of the row strings
            DO 701 IPART = 1, NIPART
               IBOT = (IPART-1)*MAXI+1
               ITOP = MIN(IBOT+MAXI-1,NROW)
               NIBTC = ITOP - IBOT + 1
*
* Obtain CSCR(I,K,JORB) = SUM(J) <K!j_+!J> . C(I,J)
*.Gather C Block
               DO JJORB = 1,NJORB
                  ICGOFF = 1 + (JJORB-1)*NKBTC*NIBTC
                  CALL MATCG(CB,CSCR(ICGOFF),NROW,NIBTC,IBOT,
     &                 NKBTC,I1(1,JJORB),XI1S(1,JJORB))
               END DO
*
* Obtain SSCR(I,K,IORB) = SUM(I) <I|i_+|K> . S(I,J)
*.Gather S Block
               DO IIORB = 1,NIORB
                  ISGOFF = 1 + (IIORB-1)*NKBTC*NIBTC
                  CALL MATCG(SB,SSCR(ISGOFF),NROW,NIBTC,IBOT,
     &                 NKBTC,I2(1,IIORB),XI2S(1,IIORB))
               END DO
               NKI = NKBTC*NIBTC
*
               FACTORC = 0.0D0
               FACTORAB = SCALEF
               CALL MATML7(RHO1S,SSCR,CSCR,NIORB,NJORB,NKI,NIORB,
     &              NKI,NJORB,FACTORC,FACTORAB,1)
* "1" means transpose matrix SSCR, so summation is over NKI,NKI as required!
*. Scatter out to complete matrix
               DO 610 JJORB = 1, NJORB
                  JORB = IBJORB-1+JJORB
                  DO 605 IIORB = 1, NIORB
                    IORB = IBIORB-1+IIORB
                    if (NTEST.ge.1000) then
                    write(6,'(A13,1X,F10.8,2X,A16,I5,A3,F10.8)')
     &                    'Adding RHO1S =',
     &                    RHO1S((JJORB-1)*NIORB+IIORB),
     &                    'to element RHO1(',
     &                    (JORB-1)*NACOB+IORB,')= ',
     &                    RHO1((JORB-1)*NACOB+IORB)
                    end if
                    RHO1((JORB-1)*NACOB+IORB) =
     &                   RHO1((JORB-1)*NACOB+IORB) +
     &                   RHO1S((JJORB-1)*NIORB+IIORB)
                    if (NTEST.ge.1000) then
                      write(6,'(3X,A3,F16.14)') ' = ',
     &                   RHO1((JORB-1)*NACOB+IORB)
                    end if
 605              CONTINUE
 610           CONTINUE

 701        CONTINUE
*. /\ end of this I partitioning
*.end of this K partitioning
            IF(KEND.EQ.0) GOTO 700
 800     CONTINUE
*.(end of loop over symmetries)
 900  CONTINUE
 1001 CONTINUE
*
C!    stop ' enforced stop in RSBBD1 '
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GSDNBB2DG(I12,
     &     RHO1UUR,RHO1UUI,
     &     RHO1UBR,RHO1UBI,
     &     RHO1BBR,RHO1BBI,
     &     RHO1BUR,RHO1BUI,
     &     IASM,IATP,IBSM,IBTP,JASM,JATP,JBSM,JBTP,
     &     NGAS,IAOC,IBOC,JAOC,JBOC,
     &     IAGRP,IBGRP,JAGRP,JBGRP,
     &     SB,CB,C2,
     &     ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX,
     &     MXPNGAS,NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &     MAXI,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,
     &     NSMOB,NSMST,NSMSX,NSMDX,
     &     NIA,NIB,NJA,NJB,MXPOBS,IPRNT,NACOB,NACOB2,
     &     RHO1S,SCLFAC,IRIL,IRIR)
*
* Contributions to density matrix from sigma block (iasm iatp, ibsm ibtp ) and
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
* NGAS : Number of As'es
* IAOC : Occpation of each AS for alpha strings in L
* IBOC : Occpation of each AS for beta  strings in L
* JAOC : Occpation of each AS for alpha strings in R
* JBOC : Occpation of each AS for beta  strings in R
* IAGRP    : IA belongs to this group of strings
* IBGRP    : IB belongs to this group of strings
* IBGRP    : IB belongs to this group of strings
* IBGRP    : IB belongs to this group of strings
* CB : Input c block
* SB : Input S block
* ADASX : sym of a+, a => sym of a+a
* ADSXA : sym of a+, a+a => sym of a
* SXSTST : Sym of sx,!st> => sym of sx !st>
* STSTSX : Sym of !st>,sx!st'> => sym of sx so <st!sx!st'>
*          is nonvanishing by symmetry
* DXSTST : Sym of dx,!st> => sym of dx !st>
* STSTDX : Sym of !st>,dx!st'> => sym of dx so <st!dx!st'>
*          is nonvanishing by symmetry
* MXPNGAS : Largest number of As'es allowed by program
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* IBORB  : Orbitals of given type and symmetry
* MAXI   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
*
* ======
* Output
* ======
* RHO1BBR/I : Real/Imaginary part of Barred-Barred     1-particle density
* RHO1UUR/I : ...................... Unbarred-Unbarred ..................
* RHO1UBR/I : ...................... Unbarred-Barred ....................
* RHO1BUR/I : ...................... Barred-Unbarred ....................
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
* RHO1S : Space for intermediate multiplication
*
*------------------------------------------------------------
* Jeppe Olsen , Winter of 1991
*------------------------------------------------------------
*
*------------------------------------------------------------
* modified for relativistic double group CI
*  Timo Fleig, April/May 1998
*              December 1998 / January 1999
*------------------------------------------------------------
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      INTEGER  ADSXA,SXSTST,STSTSX,DXSTST,STSTDX,SXDXSX
      integer DELTA_MK2
*. Input
      DIMENSION CB(*),SB(*)
      dimension IOBPTS(MXPNGAS,*),IOBPTS2(MXPNGAS,*)
      dimension NOBPTS(MXPNGAS,*),NOBPTS2(MXPNGAS,*)
*     . Output
      dimension RHO1UUR(*),RHO1UUI(*),RHO1UBR(*),RHO1UBI(*)
      dimension RHO1BBR(*),RHO1BBI(*),RHO1BUR(*),RHO1BUI(*)
*
*. Scratch
      DIMENSION SSCR(*),CSCR(*)
      DIMENSION I1(*),I2(*),I3(*),I4(*)
      dimension XI1S(*),XI2S(*),XI3S(*),XI4S(*)
      DIMENSION C2(*),RHO1S(*)
*
      CALL QENTER('GSDNBB2DG')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRNT)
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' GSDNBB2DG :  R block '
        WRITE(6,*) ' ====================='
        CALL WRTMAT(CB,NJA,NJB,NJA,NJB)
        WRITE(6,*) ' ====================='
        WRITE(6,*) ' GSDNBB2DG :  L block '
        WRITE(6,*) ' ====================='
        CALL WRTMAT(SB,NIA,NIB,NIA,NIB)
*
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in L '
        CALL IWRTMA(IAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in L '
        CALL IWRTMA(IBOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of alpha strings in R '
        CALL IWRTMA(JAOC,1,NGAS,1,NGAS)
        WRITE(6,*)
        WRITE(6,*) ' Occupation of beta  strings in R '
        CALL IWRTMA(JBOC,1,NGAS,1,NGAS)
*
        WRITE(6,'(A16,3I5)') ' MAXI,MAXK,NSMOB',MAXI,MAXK,NSMOB
*
        WRITE(6,*) 'SCLFAC =',SCLFAC
      END IF
      IACTIVE = 0
*
* Determine which part of this density block is being calculated
* real part or imaginary part
      if (IRIL.eq.IRIR) then
         if (NTEST.ge.200)
     &        write(6,*) 'Calculating real part of density matrix ...'
         IRIDEN = 1
         SCALEF = 1
      else if (IRIL.ne.IRIR) then
         if (NTEST.ge.200)
     &        write(6,*)
     $        'Calculating imaginary part of density matrix ...'
         IRIDEN = -1
         SCALEF = 1
      end if
*
*. Determine number of electrons in left and right strings and their type (UB)
      NLAEL = IELSUM(IAOC,NGAS)
      NLBEL = IELSUM(IBOC,NGAS)
      NRAEL = IELSUM(JAOC,NGAS)
      NRBEL = IELSUM(JBOC,NGAS)
*
      NLEL = NLAEL + NLBEL
      NREL = NRAEL + NRBEL
*
      if (NTEST.ge.200) then
         write(6,*) 'NLAEL =  ',NLAEL
         write(6,*) 'NLBEL =  ',NLBEL
         write(6,*) 'NRAEL =  ',NRAEL
         write(6,*) 'NRBEL =  ',NRBEL
         write(6,*) 'total left   : ',NLEL
         write(6,*) 'total right  : ',NREL
      end if
*
      if (NLEL.ne.NREL) then
         write(6,*) 'Total number of electrons ',
     &              'in left and right strings differs.'
         write(6,'(A25,4I4)') 'NLAEL,NLBEL,NRAEL,NRBEL  ',
     &                         NLAEL,NLBEL,NRAEL,NRBEL
         Call Abend2( 'gsdnbb2dg: problem with #els. in strings' )
      end if
*
      DELTA_MK2 = 2*(NLAEL-NRAEL)
      if (NTEST.ge.200)
     &     write(6,*) 'DELTA_MK2 =  ',DELTA_MK2
      IPM1 = DELTA_MK2 / 2
      if (NTEST.ge.200)
     &     write(6,*) 'IPM1 =  ',IPM1
*
      if (DELTA_MK2.eq.0) then
         if (NTEST.ge.100) then
            write(6,*) '............................'
            write(6,*) ' Entering DELTA_MK2 = 0 part'
            write(6,*) '............................'
         end if
         IF(NLBEL.GE.1.AND.IATP.EQ.JATP.AND.JASM.EQ.IASM) THEN
*
* ===================================
*  barred-barred contribution to RHO1
* ===================================
*
            if (IRIDEN.eq.1) then
               if (NTEST.ge.1000) then
                  write(6,*) 'gsbbd1dg will be called'
                  write(6,*) 'barred-barred contribution'
                  write(6,*) 'IRIDEN is  ',IRIDEN
                  write(6,*) 'REDUNDANT! Is being skipped ...'
               end if
*
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1BBR'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1BBR,NACOB2,
     &              IBSM,IBTP,JBSM,JBTP,IBGRP,NIA,
     &              NGAS,IBOC,JBOC,
     &              SB,CB,
     &              NOBPTS2,IOBPTS2,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            else if (IRIDEN.eq.-1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1BBI'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1BBI,NACOB2,
     &              IBSM,IBTP,JBSM,JBTP,IBGRP,NIA,
     &              NGAS,IBOC,JBOC,
     &              SB,CB,
     &              NOBPTS2,IOBPTS2,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            end if
         end if
*
         IF(NLAEL.GE.1.AND.IBTP.EQ.JBTP.AND.IBSM.EQ.JBSM) THEN
*
* =======================================
*  unbarred-unbarred contribution to RHO1
* =======================================
*
            CALL TRPMT3(CB,NJA,NJB,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMT3(SB,NIA,NIB,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
            if (IRIDEN.eq.1) then
               if (NTEST.ge.1000) then
                  write(6,*) 'gsbbd1dg will be called'
                  write(6,*) 'unbarred-unbarred contribution'
                  write(6,*) 'IRIDEN is  ',IRIDEN
               end if
*
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UUR'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1UUR,NACOB,
     &              IASM,IATP,JASM,JATP,IAGRP,NIB,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              NOBPTS,IOBPTS,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            else if (IRIDEN.eq.-1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UUI'
                 write(6,*) '////////////////////'
               end if
               CALL GSBBD1DG(RHO1UUI,NACOB,
     &              IASM,IATP,JASM,JATP,IAGRP,NIB,
     &              NGAS,IAOC,JAOC,
     &              SB,CB,
     &              NOBPTS,IOBPTS,MAXI,MAXK,
     &              SSCR,CSCR,I1,XI1S,I2,XI2S,
     &              NSMOB,RHO1S,
     &              SCALEF,SCLFAC,IPRNT)
            end if
         end if
*
C?      WRITE(6,*) ' GSBBD1DG was called '
C?      WRITE(6,*) ' Memory check '
C?      CALL RMEMCHK
*
* ================================
* beta-beta contribution to RHO2
* ================================
*
         IF(I12.EQ.2.AND.NBEL.GE.2) THEN
C?        WRITE(6,*) ' GSBBD2A will be called (beta)'
CTF       CALL GSBBD2A(RHO2,NACOB,IBSM,IBTP,JBSM,JBTP,IJBGRP,NIA,
CTF  &         NGAS,IBOC,JBOC,SB,CB,
CTF  &         ADSXA,SXSTST,STSTSX,SXDXSX,
CTF  &         NOBPTS,IOBPTS,MAXI,MAXK,
CTF  &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
CTF  &         NSMOB,NSMST,NSMSX,SCLFAC)
C?        WRITE(6,*) ' GSBBD2A was called '
*
C              GSBBD2A(RHO2,NACOB,ISCSM,ISCTP,ICCSM,ICCTP,IGRP,NROW,
C    &         NGAS,ISEL,ICEL,SB,CB,
C    &         ADSXA,SXSTST,STSTSX,SXDXSX,
C    &         NOBPTS,IOBPTS,MAXI,MAXK,
C    &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
C    &         NSMOB,NSMST,NSMSX)
         END IF
         IF(I12.EQ.2.AND.NAEL.GE.2) THEN
*
*
* ===================================
*  alpha-alpha contribution to RHO2
* ===================================
*
C?        WRITE(6,*) ' GSBBD2A will be called (alpha)'
CTF       CALL GSBBD2A(RHO2,NACOB,IASM,IATP,JASM,JATP,IJAGRP,NIB,
CTF  &         NGAS,IAOC,JAOC,SB,CB,
CTF  &         ADSXA,SXSTST,STSTSX,SXDXSX,
CTF  &         NOBPTS,IOBPTS,MAXI,MAXK,
CTF  &         SSCR,CSCR,I1,XI1S,I2,XI2S,X,
CTF  &         NSMOB,NSMST,NSMSX,SCLFAC)
C?      WRITE(6,*) ' GSBBD2A was called '
            CALL TRPMT3(CB,NJB,NJA,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMAT(SB,NIB,NIA,C2)
            CALL COPVEC(C2,SB,NIB*NIA)
         END IF
*
      else if (DELTA_MK2.eq.-2) then
         if (NTEST.ge.100) then
            write(6,*) '.............................'
            write(6,*) ' Entering DELTA_MK2 = -2 part'
            write(6,*) '.............................'
         end if

         if (NLAEL.ge.1.and.NLBEL.ge.1) then
*
* =====================================
*  barred-unbarred contribution to RHO1
* =====================================
*
            if (IRIDEN.eq.1) then
              if (NTEST.ge.1000) then
                 write(6,*) 'gsubd1dg will be called'
                 write(6,*) 'barred-unbarred contribution'
                 write(6,*) 'IRIDEN is  ',IRIDEN
                 write(6,*) 'IPM1 is  ',IPM1
                 write(6,*) 'REDUNDANT! Is being skipped ...'
              end if
*
              if (NTEST.ge.100) then
                write(6,*) '////////////////////'
                write(6,*) 'Writing to RHO1BUR'
                write(6,*) '////////////////////'
              end if
              call gsubd1dg(RHO1BUR,NACOB2,
     &             IASM,IBSM,IATP,IBTP,
     &             JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &             JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &             IAOC,IBOC,JAOC,JBOC,
     &             SB,CB,
     &             NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &             MAXI,MAXK,SSCR,CSCR,
     &             I1,XI1S,I2,XI2S,NSMOB,
     &             RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
           else if (IRIDEN.eq.-1) then
              if (NTEST.ge.100) then
                write(6,*) '////////////////////'
                write(6,*) 'Writing to RHO1BUI'
                write(6,*) '////////////////////'
              end if
              call gsubd1dg(RHO1BUI,NACOB2,
     &             IASM,IBSM,IATP,IBTP,
     &             JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &             JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &             IAOC,IBOC,JAOC,JBOC,
     &             SB,CB,
     &             NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &             MAXI,MAXK,SSCR,CSCR,
     &             I1,XI1S,I2,XI2S,NSMOB,
     &             RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
           end if
         end if
*
      else if (DELTA_MK2.eq.2) then
         if (NTEST.ge.100) then
            write(6,*) '............................'
            write(6,*) ' Entering DELTA_MK2 = 2 part'
            write(6,*) '............................'
         end if
*
* =====================================
*  unbarred-barred contribution to RHO1
* =====================================
*
         if (NTEST.ge.1000) then
            write(6,*) 'gsubd1dg will be called'
            write(6,*) 'unbarred-barred contribution'
            write(6,*) 'IRIDEN is  ',IRIDEN
            write(6,*) 'IPM1 is  ',IPM1
         end if
*
         if (IRIDEN.eq.1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UBR'
                 write(6,*) '////////////////////'
               end if
            call gsubd1dg(RHO1UBR,NACOB,
     &           IASM,IBSM,IATP,IBTP,
     &           JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &           JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &           IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,
     &           NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &           MAXI,MAXK,SSCR,CSCR,
     &           I1,XI1S,I2,XI2S,NSMOB,
     &           RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
         else if (IRIDEN.eq.-1) then
               if (NTEST.ge.100) then
                 write(6,*) '////////////////////'
                 write(6,*) 'Writing to RHO1UBI'
                 write(6,*) '////////////////////'
               end if
            call gsubd1dg(RHO1UBI,NACOB,
     &           IASM,IBSM,IATP,IBTP,
     &           JASM,JBSM,JATP,JBTP,IAGRP,IBGRP,
     &           JAGRP,JBGRP,NIA,NIB,NJA,NJB,NGAS,
     &           IAOC,IBOC,JAOC,JBOC,
     &           SB,CB,
     &           NOBPTS,NOBPTS2,IOBPTS,IOBPTS2,
     &           MAXI,MAXK,SSCR,CSCR,
     &           I1,XI1S,I2,XI2S,NSMOB,
     &           RHO1S,IPM1,SCALEF,SCLFAC,IPRNT)
         end if
*
      else
         write(6,*) 'DELTA_MK2 is  ',DELTA_MK2
         write(6,*) 'These cases are not debugged yet'
         write(6,*) 'I am skipping the relevant density routines'
      end if
*
      idelabs = abs(DELTA_MK2)
      if (IDELABS.ge.4) then
         write(6,*) 'CAUTION !!!'
         write(6,*) 'CAUTION !!!'
         write(6,*) 'Skipping RHO2 routines ...'
         write(6,*) 'CAUTION !!!'
         write(6,*) 'CAUTION !!!'
         goto 666
      end if
*
      if (DELTA_MK2.eq.-4) then
         if (NTEST.ge.100) then
            write(6,*) '.............................'
            write(6,*) ' Entering DELTA_MK2 = -4 part'
            write(6,*) '.............................'
         end if
*
* ===================================
*  alpha-beta contribution to RHO2
* ===================================
*
         IF(I12.EQ.2.AND.NAEL.GE.1.AND.NBEL.GE.1) THEN
*. Routine uses transposed blocks
            CALL TRPMT3(CB,NJA,NJB,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMT3(SB,NIA,NIB,C2)
            CALL COPVEC(C2,SB,NIA*NIB)
C?      WRITE(6,*) ' GSBBD2B will be called '
CTF     CALL GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
CTF  &                    JASM,JATP,JBSM,JBTP,NJA,NJB,
CTF  &                    IJAGRP,IJBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
CTF  &                    SB,CB,ADSXA,STSTSX,
CTF  &                    NOBPTS,IOBPTS,MAXK,
CTF  &                    I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
CTF  &                    NSMOB,NSMST,NSMSX,NSMDX,IUSEAB,
CTF  &                    SSCR,CSCR,NACOB,NTEST,SCLFAC)
C?      WRITE(6,*) ' GSBBD2B was called '
CTF  &
C     GSBBD2B(RHO2,IASM,IATP,IBSM,IBTP,NIA,NIB,
C    &                        JASM,JATP,JBSM,JBTP,NJA,NJB,
C    &                  IAGRP,IBGRP,NGAS,IAOC,IBOC,JAOC,JBOC,
C    &                  SB,CB,ADSXA,STSTSX,
C    &                  NOBPTS,IOBPTS,MAXK,
C    &                  I1,XI1S,I2,XI2S,I3,XI3S,I4,XI4S,X,
C    &                  NSMOB,NSMST,NSMSX,NSMDX,IUSEAB,
C    &                  CJRES,SIRES,NORB,NTEST)
            CALL TRPMT3(CB,NJB,NJA,C2)
            CALL COPVEC(C2,CB,NJA*NJB)
            CALL TRPMAT(SB,NIB,NIA,C2)
            CALL COPVEC(C2,SB,NIB*NIA)
         end if
      else if (DELTA_MK2.eq.4) then
         if (NTEST.ge.100) then
            write(6,*) '............................'
            write(6,*) ' Entering DELTA_MK2 = 4 part'
            write(6,*) '............................'
         end if
*
         write(6,*) 'DELTA_MK2 is  ',DELTA_MK2
         write(6,*) 'This part is still missing ...'
      END IF
*
 666  continue
      CALL QEXIT('GSDNBB2DG')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GSUBD1DG(RHO1,NACOB,
     &     ISCSMA,ISCSMB,ISCTPA,ISCTPB,
     &     ICCSMA,ICCSMB,ICCTPA,ICCTPB,
     &     IAGRP,IBGRP,JAGRP,JBGRP,
     &     NIROWA,NIROWB,NJROWA,NJROWB,
     &     NGAS,ISELA,ISELB,ICELA,ICELB,
     &     SB,CB,
     &     NOBPTS_U,NOBPTS_B,IOBPTS_U,IOBPTS_B,
     &     MAXL,MAXK,
     &     SSCR,CSCR,I1,XI1S,I2,XI2S,
     &     NSMOB,RHO1S,
     &     IPM1,SCALEF,SCLFAC,IPRNT)
*
* Contributions to one electron density matrix from column excitations
*
* GAS version, August 95 , Jeppe Olsen
*
*  modified for relativistic double group calculation
*  summer/fall of '97
*  april 1998, Timo Fleig
*
* =====
* Input
* =====
* RHO1  : One body density matrix to be calculated/updated
* NACOB : Number of active unbarred orbitals
* NACOB2 : Number of active barred orbitals
* ISCSM,ISCTP : Symmetry and type of sigma columns
* ICCSM,ICCTP : Symmetry and type of C     columns
* IAGRP,IBGRP : String group of S columns
* JAGRP,JBGRP : String group of C columns
* NIROWA,NIROWB : Number of rows in S block
* NJROWA,NJROWB : Number of rows in C block
* NGAS : Number of active spaces
* ISEL : Number of electrons per AS for S block
* ICEL : Number of electrons per AS for C block
* SB   : Input S block
* CB   : Input C block
* MXPNGAS : Max number of AS spaces ( program parameter )
* NOBPTS  : Number of orbitals per type and symmetry
* IOBPTS : base for orbitals of given type and symmetry
* NSMOB : Number of symmetries of spinors
* MAXL   : Largest Number of ' spectator strings 'treated simultaneously
* MAXK   : Largest number of inner resolution strings treated at simult.
* IPM1 : distinguish barred-unbarred from unbarred-barred part
*         +1 : operator i+ j_ (unbarred-barred)
*         -1 : operator i_+ j (barred-unbarred)
*
* ======
* Output
* ======
* RHO1 : Calculated/updated density block
*
* =======
* Scratch
* =======
*
* SSCR, CSCR : at least MAXLJ*MAXL*MAXK, where MAXLJ is the
*              largest number of orbital pairs of given symmetries and
*              types.
* I1, XI1S   : MAXK*Max number of orbitals of given type and symmetry
* I2, XI2S   : MAXK*Max number of orbitals of given type and symmetry
*              type and symmetry
* RHO1S : Space for one electron density
*
* Jeppe Olsen, Winter of 1991
* Updated for GAS , August '95
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
*. General input
      INTEGER NOBPTS_U(MXPNGAS,*),NOBPTS_B(MXPNGAS,*)
      integer IOBPTS_U(MXPNGAS,*),IOBPTS_B(MXPNGAS,*)
*.Input
      INTEGER ISELA(NGAS),ICELA(NGAS)
      integer ISELB(NGAS),ICELB(NGAS)
*
#include "symm.inc"
*
      DIMENSION CB(*),SB(*)
*.Output
      DIMENSION RHO1(*)
*.Scatch
      DIMENSION SSCR(*),CSCR(*),RHO1S(*)
      DIMENSION I1(MAXK,*),XI1S(MAXK,*)
      DIMENSION I2(MAXK,*),XI2S(MAXK,*)
*
*.Local arrays
      DIMENSION ITP(16),JTP(16)
*
      NTESTL = 0000
      NTEST = max(NTESTL,IPRNT)
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*)
         WRITE(6,*) ' ===================='
         WRITE(6,*) ' GSUBD1DG in action '
         WRITE(6,*) ' ===================='
         WRITE(6,*)
         WRITE(6,*) ' Occupation of active left unbarred strings '
         CALL IWRTMA(ISELA,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active Right unbarred strings '
         CALL IWRTMA(ICELA,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active left barred strings '
         CALL IWRTMA(ISELB,1,NGAS,1,NGAS)
         WRITE(6,*) ' Occupation of active Right barred strings '
         CALL IWRTMA(ICELB,1,NGAS,1,NGAS)
*
         WRITE(6,*) ' GSBBD1DG, sclfac ',SCLFAC
      END IF
*
      IFRST = 1
      JFRST = 1
*
*  determine C(Ja,Ib,j) = sum(Jb) <Jb| j_+ |Ib> C(Ja,Jb)
*  (if IPM1 = -1 determine C(Ia,Jb,j) = sum(Ja) <Ja| j+ |Ia> C(Ja,Jb) )
      if (IPM1.eq.1) then
         call adtyp(NSXTPB,JTYP,NGAS,ICELB,ISELB)
         call adtyp(NSXTPA,ITYP,NGAS,ISELA,ICELA)
      else if (IPM1.eq.-1) then
         call adtyp(NSXTPA,JTYP,NGAS,ICELA,ISELA)
         call adtyp(NSXTPB,ITYP,NGAS,ISELB,ICELB)
      else
         write(6,*) ' IPM1 out of range '
         write(6,*) ' IPM1 = ',IPM1
         goto 2000
      end if
      if (NSXTPA.eq.0 .or. NSXTPB.eq.0) goto 2000
*
* determine sym(J) for non-vanishing <Jb| j_+ |Ib>
*                                   (<Ja| j+ |Ia>)
* determine sym(I) for non-vanishing <Ia| i+ |Ja>
*                                   (<Ib| i_+ |Jb>)
      if (IPM1.eq.1) then
* transpose sigma block
         call trpmt3(SB,NIROWA,NIROWB,CSCR)
         call copvec(CSCR,SB,NIROWA*NIROWB)
         JSM = IDBGMULT(ICCSMB,IADJSYM(ISCSMB))
         ISM = IDBGMULT(ISCSMA,IADJSYM(ICCSMA))
         NJORB = NOBPTS_B(JTYP,JSM)
         NIORB = NOBPTS_U(ITYP,ISM)
      else
* transpose C block
         call trpmt3(CB,NJROWA,NJROWB,CSCR)
         call copvec(CSCR,CB,NJROWA*NJROWB)
         JSM = IDBGMULT(ICCSMA,IADJSYM(ISCSMA))
         ISM = IDBGMULT(ISCSMB,IADJSYM(ICCSMB))
         NJORB = NOBPTS_U(JTYP,JSM)
         NIORB = NOBPTS_B(ITYP,ISM)
      end if
      if (NTEST.ge.1000)
     &     write(6,*) ' JSM ISM ', JSM,ISM
*
* start with batches of K strings according to sum(Kb) <Jb| j_+ |Kb><Kb|1|Ib>
*                                             (sum(Ka) <Ja| j+ |Ka><Ka|1|Ia> )
      if (NJORB.eq.0 .or. NIORB.eq.0) goto 2000
      KBOT = 1 - MAXK
      KTOP = 0
 700  continue
      KBOT = KBOT + MAXK
      KTOP = KTOP + MAXK
      J12 = 2
      K12 = 2
      if (JFRST.eq.1) KFRST = 1
*
*  determine mapping <Jb| j_+ |Kb>
*                   ( <Ja| j+ |Ka> )
      if (IPM1.eq.1) then
         JBORB = IOBPTS_B(JTYP,JSM)
         if (NTEST.ge.1000) write(6,*) 'offset right orb ',JBORB
         call adst_gas(JBORB,NJORB,ICCTPB,ICCSMB,JBGRP,KBOT,KTOP,
     &        I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &        J12,K12,KKJACT,SCLFAC)
      else if (IPM1.eq.-1) then
         JBORB = IOBPTS_U(JTYP,JSM)
         if (NTEST.ge.1000) write(6,*) 'offset right orb ',JBORB
         call adst_gas(JBORB,NJORB,ICCTPA,ICCSMA,JAGRP,KBOT,KTOP,
     &        I1,XI1S,MAXK,NKBTC,KEND,JFRST,KFRST,
     &        J12,K12,KKJACT,SCLFAC)
      end if
*
      JFRST = 0
      KFRST = 0
* loop over batches of La, according to sum(La) <Ia| i+ |La><La|Ja>
      LBOT = 1 - MAXL
      LTOP = 0
 600  continue
      LBOT = LBOT + MAXL
      LTOP = LTOP + MAXL
      I12 = 1
      L12 = 1
      if (IFRST.eq.1) LFRST = 1
*
*  determine mapping <Ia| i+ |La>
*                   (<Ib| i_+ |Lb>)
      if (IPM1.eq.1) then
         IBORB = IOBPTS_U(ITYP,ISM)
         if (NTEST.ge.1000) write(6,*) 'offset left orb ',IBORB
         call adst_gas(IBORB,NIORB,ISCTPA,ISCSMA,IAGRP,LBOT,LTOP,
     &        I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &        I12,L12,LLIACT,SCLFAC)
      else
         IBORB = IOBPTS_B(ITYP,ISM)
         if (NTEST.ge.1000) write(6,*) 'offset left orb ',IBORB
         call adst_gas(IBORB,NIORB,ISCTPB,ISCSMB,IBGRP,LBOT,LTOP,
     &        I2,XI2S,MAXL,NLBTC,LEND,IFRST,LFRST,
     &        I12,L12,LLIACT,SCLFAC)
      end if
*
*  obtain CSCR(Ja,Kb,j) = sum(Jb) <Jb| j_+ |Kb><Kb|Ib> C(Ja,Jb)
*        (CSCR(Jb,Ka,j) = sum(Ja) <Ja| j+ |Ka><Ka|Ia> C(Jb,Ja) )
*
*  obtain SSCR(La,Ib,i) = sum(Ia) <Ia| i+ |La><La|Ja> C(Ib,Ia)
*        (SSCR(Lb,Ia,i) = sum(Ib) <Ib| i_+ |Lb><Lb|Jb> C(Ia,Ib) )
      if (NTEST.ge.1000) write(6,*) 'IPM1 =  ',IPM1
      if (IPM1.eq.1) then
         do JJORB = 1,NJORB,1
            ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
            call matcg(CB,CSCR(ICGOFF),NJROWA,NLBTC,LBOT,
     &           NKBTC,I1(1,JJORB),XI1S(1,JJORB))
         end do
         do IIORB = 1,NIORB,1
            ISGOFF = 1 + (IIORB-1)*NKBTC*NLBTC
            call matcg(SB,SSCR(ISGOFF),NIROWB,NKBTC,KBOT,
     &           NLBTC,I2(1,IIORB),XI2S(1,IIORB))
         end do
      else
         do JJORB = 1,NJORB,1
            ICGOFF = 1 + (JJORB-1)*NKBTC*NLBTC
            call matcg(CB,CSCR(ICGOFF),NJROWB,NLBTC,LBOT,
     &           NKBTC,I1(1,JJORB),XI1S(1,JJORB))
         end do
         do IIORB = 1,NIORB,1
            ISGOFF = 1 + (IIORB-1)*NKBTC*NLBTC
            call matcg(SB,SSCR(ISGOFF),NIROWA,NKBTC,KBOT,
     &           NLBTC,I2(1,IIORB),XI2S(1,IIORB))
         end do
      end if
*
      NKI = NKBTC*NLBTC
      FACTORC = 0.0D0
      FACTORAB = SCALEF
*
* Check transposing scheme before matrix multiply:
*  if IPM1.eq.1   SB has been transposed
*  if IPM1.eq.-1  CB has been transposed
* transpose back where required
*
      if (IPM1.eq.1) then
         call matml7(RHO1S,SSCR,CSCR,NIORB,NJORB,NIORB,NKI,
     &        NKI,NJORB,FACTORC,FACTORAB,0)
      else if (IPM1.eq.-1) then
         call trpmt3(SSCR,NIROWB,NIROWA,SB)
         call copvec(SB,SSCR,NIROWA*NIROWB)
         call trpmt3(CSCR,NJROWB,NJROWA,CB)
         call copvec(CB,CSCR,NJROWA*NJROWB)
         call matml7(RHO1S,SSCR,CSCR,NIORB,NJORB,NIORB,NKI,
     &        NKI,NJORB,FACTORC,FACTORAB,0)
      end if
*
*  Scatter out to complete matrix
*  We have 2 sets of spinors, and we are calculating density matrix
*  blocks for each of the four resulting blocks separately. This means
*  that the offsets and numbers of spinors are absolute values in each
*  of the blocks in question. The number of active orbitals corresponds
*  to the number of spinors in the row index of the block we are
*  calculating. This is taken care of in the form of separate variables
*  passed from above.
*
      do 900 JJORB = 1,NJORB,1
         JORB = JBORB - 1 + JJORB
         do 800 IIORB = 1,NIORB,1
            IORB = IBORB - 1 + IIORB
            if (NTEST.ge.1000) then
              write(6,'(A13,1X,F10.8,2X,A16,I5,A3,F10.8)')
     &              'Adding RHO1S =',
     &              RHO1S((JJORB-1)*NIORB+IIORB),
     &              'to element RHO1(',
     &              (JORB-1)*NACOB+IORB,')= ',
     &              RHO1((JORB-1)*NACOB+IORB)
            end if
            RHO1((JORB-1)*NACOB+IORB) = RHO1((JORB-1)*NACOB+IORB)
     &           + RHO1S((JJORB-1)*NIORB+IIORB)
            if (NTEST.ge.1000) then
                write(6,'(3X,A3,F16.14)') ' = ',
     &             RHO1((JORB-1)*NACOB+IORB)
            end if
 800     continue
 900  continue
*
      if (KEND.eq.0) goto 700
      if (LEND.eq.0) goto 600
*
 2000 continue
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Write natural (optimized) spinors to files
*   LUMOUTR (real part)
*   LUMOUTI (imaginary part)
*
      subroutine wrt_natspi(ONATSPI,IDIM,LUMOUTR,
     &                  LUMOUTI,ISPSM)
*
      implicit real*8 (A-H,O-Z)
#include "ipoist8.inc"
*
      dimension ONATSPI(2*IDIM,IDIM)
*
 100  format(4F18.14)
*
      write(6,*)
      write(6,*) 'Writing real part of natural spinors'
      write(6,*) 'to file NATSPIR.'
      write(LUMOUTR,'(A16,I2)') 'SPINOR SYMMETRY ',ISPSM
      write(LUMOUTI,'(A16,I2)') 'SPINOR SYMMETRY ',ISPSM
      do ICOL=1,IDIM,1
        write(LUMOUTR,'(A14,I4)') 'Spinor number ',ICOL
        write(LUMOUTR,100) (ONATSPI(I,ICOL),I=1,IDIM,1)
        write(LUMOUTI,'(A14,I4)') 'Spinor number ',ICOL
        write(LUMOUTI,100) (ONATSPI(I,ICOL),I=IDIM+1,2*IDIM,1)
      end do
*
      return
      end
