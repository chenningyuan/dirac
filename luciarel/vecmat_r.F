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
      SUBROUTINE EXPCIV_REL(ISM,ISPCIN,LUIN,
     &                  ISPCUT,LUUT,LBLK,
     &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
*
* Expand CI vector in CI space ISPCIN to CI vector in ISPCUT
* Input vector is supposed to be on LUIN
* Output vector will be placed on unit LUUT
*. If ICOPY .ne. 0 the output vectors will be copied to LUIN
*
* Storage form is defined by ICISTR
*
* Jeppe Olsen, February 1994
* GAS version August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "csm.inc"
#include "cgas.inc"
#include "gasstr.inc"

*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'EXPCIV')
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*
*. Allowed combinations of strings types for input and output
*. spaces
*
      call memmar(KLIABI,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABI')
      call memmar(KLIABU,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABU')
      CALL IAIBCM_GAS(LCMBSPC(ISPCIN),ICMBSPC(1,ISPCIN),IGSOCCX,
     &            NOCTPA,NOCTPB,
     &            ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,
     &            WORK(KLIABI),NTEST)
*
      CALL IAIBCM_GAS(LCMBSPC(ISPCUT),ICMBSPC(1,ISPCUT),IGSOCCX,
     &            NOCTPA,NOCTPB,
     &            ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,
     &            WORK(KLIABU),NTEST)
*
* type of each symmetry block ( full, lower diagonal, absent )
*
      call memmar(KLBLIN,NSMST,'ADDL  ',1,'KLBLIN')
      call memmar(KLBLUT,NSMST,'ADDL  ',1,'KLBLUT')
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLIN),IDUMMY)
        CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLUT),IDUMMY)

*
*. Number of dets etc per TTS block
*
      call memmar(KLTTSII,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
      call memmar(KLTTSNI,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLIN),NSMST,WORK(KLIABI),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSII),WORK(KLTTSNI),
     &            NCOMBI,0)
*
      call memmar(KLTTSIU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
      call memmar(KLTTSNU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLUT),NSMST,WORK(KLIABU),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSIU),WORK(KLTTSNU),
     &            NCOMBU,0)
*
*. Allocate memory for blocks of input and output space in
*  accordance with ICISTR
*
      IF(ICISTR.EQ.1) THEN
        LENGTHI = NCOMBI
        LENGTHU = NCOMBU
      ELSE IF (ICISTR .EQ. 2 ) THEN
        LENGTHI = MXSB
        LENGTHU = MXSB
      ELSE IF (ICISTR.EQ. 3 ) THEN
        LENGTHI = MXSOOB
        LENGTHU = MXSOOB
      END IF
*
      call memmar(KLBLI,LENGTHI,'ADDL  ',2,'KLBLI ')
      call memmar(KLBLU,LENGTHU,'ADDL  ',2,'KLBLU ')
*
*. and now : Let another subroutine complete the taks
*
      CALL REWINE(LUIN,-1)
      CALL REWINE(LUUT,-1)
*
*. Print for testing initial vectors out
*
C?    WRITE(6,*) ' Initial vectors in EXPCIV '
C?    WRITE(6,*) ' Initial vectors in EXPCIV '
C?    WRITE(6,*) ' Initial vectors in EXPCIV '
*
C?    DO IROOT = 1, NROOT
C?    WRITE(6,*) ' Root number ', IROOT
*
C?      IF(ICISTR.EQ.1) THEN
C?          CALL WRTMAT(WORK(KLBLI),1,NCOMBUT,1,NCOMBUT)
C?      ELSE
C?        CALL WRTVCD(WORK(KLBLI),LUIN,0,-1)
C?      END IF
C?    END DO
C?    CALL REWINE(LUIN,-1)
*
      write(6,*) ' ICISTR = ', ICISTR
      DO IROOT = 1, NROOT
*. Input vector should be first vector on file so
        IF(IROOT.EQ.1) THEN
          LLUIN = LUIN
        ELSE
          IF(ICISTR.EQ.1) THEN
            CALL REWINE(LUSCR,-1)
            CALL FRMDSC(WORK(KLBLI),NCOMBI,-1,LUIN,IMZERO,IAMPACK)
            CALL  TODSC(WORK(KLBLI),NCOMBI,-1,LUSCR)
            CALL REWINE(LUSCR,-1)
            LLUIN = LUSCR
          ELSE
*. With the elegance of an elephant
            CALL REWINE(LUSCR,-1)
            CALL REWINE(LUIN,-1)
            DO JROOT = 1, IROOT
              CALL REWINE(LUSCR,-1)
              CALL COPVCD_REL(LUIN,LUSCR,WORK(KLBLI),0,-1)
            END DO
            CALL REWINE(LUSCR,-1)
C?          write(6,*) ' Vector copied to LUSCR '
C?          CALL WRTVCD(WORK(KLBLI),LUSCR,1,-1)
C?          CALL REWINE(LUSCR,-1)
            LLUIN = LUSCR
          END IF
        END IF
*
        ITTSS_ORD = 2
        CALL EXPCIVS_REL(LLUIN,WORK(KLBLI),NCOMBI,WORK(KLTTSNI),
     &       WORK(KLTTSII),
     &       WORK(KLIABI),NOCTPA,NOCTPB,WORK(KLBLIN),
     &       LUUT,WORK(KLBLU),NCOMBU,WORK(KLTTSNU),
     &       WORK(KLTTSIU),WORK(KLIABU),
     &       WORK(KLBLUT),
     &       ICISTR,IDC,NSMST,
     &       LBLK,ITTSS_ORD)
*
      END DO
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINE(LUIN,-1)
        CALL REWINE(LUUT,-1)
        DO IROOT = 1, NROOT
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC(WORK(KLBLU),NCOMBU,-1,LUUT,IMZERO,IAMPACK)
            CALL  TODSC(WORK(KLBLU),NCOMBU,-1,LUIN)
          ELSE
            CALL COPVCD_REL(LUUT,LUIN,WORK(KLBLU),0,-1)
          END IF
        END DO
      END IF
*
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
*
C?    CALL REWINE(LUUT,-1)
C?    DO IROOT = 1, NROOT
C?    WRITE(6,*) ' Root number ', IROOT
*
C?      IF(ICISTR.EQ.1) THEN
C?          CALL WRTMAT(WORK(KLBLU),1,NCOMBUT,1,NCOMBUT)
C?      ELSE
C?        CALL WRTVCD(WORK(KLBLU),LUUT,0,-1)
C?      END IF
C?    END DO
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'EXPCV2')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE EXPCIVS_REL(LUIN,VECIN,NCOMBIN,NTTSIN,ITTSIN,IABIN,
     &                       NOCTPA,NOCTPB,IBLTPIN,
     &                       LUUT,VECUT,NCOMBUT,NTTSUT,ITTSUT,IABUT,
     &                       IBLTPUT,
     &                       ICISTR,IDC,NSMST,LBLK,ITTSS_ORD)
*
* Obtain those part of vector in cispace UT ,
* that can be obtained from terms in cispace IN
*
* Input vector on LUIN, Output vector in LUUT
* Output vector is supposed on start of vector
*
* LUIN is assumed to be single vector file,
* so rewinding will place vector on start of vector
*
* Both files are assumed on start of vector
*
* Jeppe Olsen, February 1994
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Input
      INTEGER IABIN(NOCTPA,NOCTPB),IABUT(NOCTPA,NOCTPB)
      INTEGER IBLTPIN(NSMST),IBLTPUT(NSMST)
      INTEGER NTTSIN(NOCTPA,NOCTPB,NSMST),ITTSIN(NOCTPA,NOCTPB,NSMST)
      INTEGER NTTSUT(NOCTPA,NOCTPB,NSMST),ITTSUT(NOCTPA,NOCTPB,NSMST)
      DIMENSION VECIN(*)
C     DIMENSION IATPUI(*),IBTPUI(*)
*. Output
      DIMENSION VECUT(*)
*
      WRITE(6,*) ' EXPCIVS in action'
      WRITE(6,*) '        Number of input  parameters',NCOMBIN
      WRITE(6,*) '        Number of output parameters',NCOMBUT
*
      IF(ICISTR.EQ.1) THEN
        CALL FRMDSC(VECIN,NCOMBIN,-1,LUIN,IMZERO,IAMPACK)
C       WRITE(6,*) ' Input vector '
C       CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
      END IF
*
*. Loop over TTS blocks of output vector
*
      IATPIN = 1
      IBTPIN = 0
      IASMIN = 1
*
      IATPUT = 1
      IBTPUT = 0
      IASMUT = 1
*
 1000 CONTINUE
*. Next output block
        CALL NXTBLK(IATPUT,IBTPUT,IASMUT,NOCTPA,NOCTPB,NSMST,
     &              IBLTPUT,IDC,NONEWUT,IABUT,ITTSS_ORD)
*. Corresponding input TTS block
        JATPIN = IATPUT
        JBTPIN = IBTPUT
        IF(IABIN(JATPIN,JBTPIN).EQ.0) THEN
          IZERO = 1
        ELSE
          IZERO = 0
        END IF
*
        NELMNT = NTTSUT(IATPUT,IBTPUT,IASMUT)
        IF(ICISTR.EQ.1) THEN
          IF(IZERO.EQ.0)
     &    IOFFIN = ITTSIN(JATPIN,JBTPIN,IASMUT)
          IOFFUT = ITTSUT(IATPUT,IBTPUT,IASMUT)
          IF(IZERO.EQ.0) THEN
            CALL COPVEC(VECIN(IOFFIN),VECUT(IOFFUT),NELMNT)
          ELSE
            ZERO = 0.0D0
            CALL SETVEC(VECUT(IOFFUT),ZERO,NELMNT)
          END IF
        END IF
*
        IF(ICISTR.NE.1.AND.NONEWUT.EQ.0) THEN
          IF(IZERO.EQ.1) THEN
C           ZERO = 0.0D0
C           CALL SETVEC(VECUT,ZERO,NELMNT)
            CALL ITODS(NELMNT,1,-1,LUUT)
            CALL ZERORC_REL(-1,LUUT)
C           CALL TODSC(VECUT,NELMNT,-1,LUUT)
          ELSE
*. Obtain input block
  999      CONTINUE
           CALL NXTBLK(IATPIN,IBTPIN,IASMIN,NOCTPA,NOCTPB,NSMST,
     &                 IBLTPIN,IDC,NONEWIN,IABIN,ITTSS_ORD)
           IF(NONEWIN.NE.0) THEN
             CALL REWINE(LUIN,-1)
             IATPIN = 1
             IBTPIN = 0
             IASMIN = 1
             GOTO 999
           ELSE IF (NONEWIN.EQ. 0 ) THEN
             CALL IFRMDS(LENGTH,1,-1,LUIN)
             CALL FRMDSC(VECIN,LENGTH,-1,LUIN,IMZERO,IAMPACK)
             IF(IATPIN.EQ.JATPIN.AND.IBTPIN.EQ.JBTPIN.AND.
     &          IASMIN.EQ.IASMUT) THEN
*. Correct block, save it
                  CALL ITODS(LENGTH,1,-1,LUUT)
                  IF(IMZERO.EQ.1) THEN
                    CALL ZERORC_REL(-1,LUUT)
                  ELSE
                    IF(IAMPACK.EQ.0) THEN
                      CALL TODSC(VECIN,LENGTH,-1,LUUT)
                    ELSE
                      CALL TODSCP(VECIN,LENGTH,-1,LUUT)
                    END IF
                  END IF
             ELSE
               GOTO 999
             END IF
           END IF
          END IF
        END IF
      IF(NONEWUT.EQ.0) GOTO 1000
*. End of loop over output blocks
      IF(ICISTR.EQ.1) THEN
        CALL TODSC(VECUT,NCOMBUT,-1,LUUT)
      ELSE
        CALL ITODS(-1,1,-1,LUUT)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' EXPCIVS Speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' ============ '
        WRITE(6,*) ' Input Vector '
        WRITE(6,*) ' ============ '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
        ELSE
          CALL WRTVCD(VECIN,LUIN,1,LBLK)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Output Vector '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECUT,1,NCOMBUT,1,NCOMBUT)
        ELSE
          CALL WRTVCD(VECUT,LUUT,1,LBLK)
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_BLOCKS_FROM_DISC_REL
     &           (LU,NBLOCK,IOFF,IBLOCK,NBLOCKT,C,IREW)
*
* Obtain blocks IOFF - IOFF + NBLOCK from file LU
*
* Jeppe Olsen, January 1996
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION C(*)
      DIMENSION IBLOCK(8,*)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "csm.inc"
#include "stinf.inc"
*
*. Local Scratch
*. Should be uncommented if NEW_OR_OLD is set to 2
*
C     PARAMETER (MXPNBLK = 33000)
C     DIMENSION ISCR(3*MXPNBLK)
*
C?    WRITE(6,*) ' ENTERING GET_BLOCKS_FROM_DISC'
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
CM    IREW = 1
      ISCAL = 0
      CALL GET_TTS_BATCHN(C,NBLOCK,IBLOCK(1,IOFF),NBLOCKT,IBLOCK,
     &                     NOCTPA,NOCTPB,NSMST,
     &                     WORK(KNSTSO(IATP)), WORK(KNSTSO2(IBTP)),
     &                     IDC,LU,IREW,ISCAL)

*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GSTTBLD(CTT,IATP,IASM,IBTP,IBSM,
     &                  NSASO,NSBSO,PSSIGN,IDC,
     &                  PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
*
* obtain  determinant block (iatp iasm, ibtp ibsm )
* from vector packed in combination format according to IDC
*
*. If ISCALE = 1, the routine scales and returns the block
*  in determinant normalization, and SCLFAC = 1.0D0
*
* If ISCALE = 0, the routine does not perform any overall
* scaling, and a scale factor is returned in SCLFAC
*
* IF ISCALE = 0, zero blocks are not set explicitly to zero,
* instead  zero is returned in SCLFAC
*
* ISCALE, SCLFAC added May 97
*
* Simplified version working only for vectors on disc
*
      IMPLICIT REAL*8  (A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION CTT(*),NSASO(NSMST, *),NSBSO(NSMST, *)
      DIMENSION SCR(*)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN 
        write(6,*) ' GSTTBLD, IATP,IASM,IBTP,IBSM,ISCALE'
        write(6,*)            IATP,IASM,IBTP,IBSM,ISCALE     
        WRITE(6,*) ' LUC = ', LUC
      END IF
* =================
* Read in from disc
* =================
      CALL IFRMDS(LBL,1,-1,LUC)
      IF(NTEST.GE.100) write(6,*) ' LBL = ', LBL
      IF(ISCALE.EQ.1) THEN
        CALL FRMDSC(SCR,LBL,-1,LUC,IMZERO,IAMPACK)
      ELSE
        NO_ZEROING = 1
        CALL FRMDSC2(SCR,LBL,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
      END IF
*
      IF(IMZERO.EQ.1.AND.ISCALE.EQ.0) THEN
        SCLFAC = 0.0D0
      ELSE
        NAST = NSASO(IASM,IATP)
        NBST = NSBSO(IBSM,IBTP)
        IF(LBL.NE.0) THEN
          CALL SDCMRF(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &         IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
        ELSE
          SCLFAC = 0.0D0
        END IF
      END IF
*
C?    WRITE(6,*) ' ISCALE and SCLFAC on return in GSTTBL',
C?   &ISCALE,SCLFAC
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE GSTTBLD2(CTT,LUC)
*
*     obtain  determinant block from disc
*
      IMPLICIT REAL*8  (A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION CTT(*)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' LUC = ', LUC
      END IF
*     =================
*     Read in from disc
*     =================
      CALL IFRMDS(LBL,1,-1,LUC)
      CALL FRMDSC(CTT,LBL,-1,LUC,IMZERO,IAMPACK)
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE INTBLK_IN_TO_OUT(IPSM,IQSM,IRSM,ISSM,I1234I,
     &           I12I,I34I,NLISTO,I1234O,I12O,I34O,
     &           ICCSMO,IMOSF_SP,ISPSF_MO,NIRRIN,NIRROUT,
     &           ISPITRA,ISPITRA_RC,
     &           NBLOCKO,IBLOCKO,LISTO,INDEXO,ISIGNO,IRIO,NCOMP)
*
* A given symmetry block of two-electron integrals
* is given by IPSM, IQSM, IRSM, ISSM.
* The permutational symmetry of the input integrals
* are defined by I1234I,I12I,I34I,ICCSMO
*
* NLISTO lists of output integrals are to be generated.
* These are determined by I1234O,I12O,I34O.
*
* Obtain the output blocks that this given input block
* will contribute to
*
* NBLOCKO : Number of output blocks
* IBLOCKO : Symmetry labels of output blocks
*           (in double group symmetry)
* ISIGNO  : Sign from permuting input symmetry block
*           (unused!)
* LISTO :  Output list of integral to which each block belongs
* INDEXO : Mapping of index in output block to
*          index in original input symmetry block
* IRIO   : Real or imaginary block
*
*
* The largest allowed number of output blocks is set by
* MXBLKO
*
* The blocks are - if required seperated into real and imaginary
* blocks
*
*
* Timo Fleig + Jeppe Olsen, April 1997
*
      IMPLICIT REAL(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "input.inc"
*
*. Input
      INTEGER I1234O(NLISTO),I12O(NLISTO),I34O(NLISTO)
      INTEGER ICCSMO(NLISTO)
      INTEGER IMOSF_SP(8,2),ISPSF_MO(8,2)
      INTEGER ISPITRA(4,NLISTO),ISPITRA_RC(NLISTO)
*
*. Note ISPITRA_RC(I) = 1 => This transformation is real
*                     = 2 => This transformation is imaginary
*                     = 3 => Contains both an imaginary and real part
*. Output
      INTEGER IBLOCKO(4,MXBLKO),ISIGNO(MXBLKO)
      INTEGER INDEXO (4,MXBLKO),LISTO(MXBLKO)
      INTEGER IRIO(MXBLKO)
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*) ' Welcome to INTBLK_IN_TO_OUT '
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,'(A,4I3)') ' IPSM IQSM IRSM ISSM for input block :',
     &               IPSM,IQSM,IRSM,ISSM
        WRITE(6,'(A,3I3)') ' I1234I,I12I,I34I :',I1234I,I12I,I34I
*
*
        WRITE(6,*) ' Number of output lists ', NLISTO
        WRITE(6,*) ' Output integral lists '
        WRITE(6,*)
        WRITE(6,*)
     &  ' List  I12  I34  I1234 ICC  SP1  SP2  SP3  SP4'
        WRITE(6,*)
     &  '=============================================='
        DO JLISTO = 1, NLISTO
        WRITE(6,'(1H , 9I5)')
     &   JLISTO,I12O(JLISTO),I34O(JLISTO),I1234O(JLISTO),
     &   ICCSMO(JLISTO),
     &   ISPITRA(1,JLISTO),ISPITRA(2,JLISTO),ISPITRA(3,JLISTO),
     &   ISPITRA(4,JLISTO)
        END DO
*
      END IF
*
      NBLOCKO = 0
*. Each input symmetry can in general correspond to two
*. output symmetries ( pointgroupXspinfunction=>double group).
      DO IPSF = 1,2
      DO IQSF = 1,2
*. Well here we are going to insert some resrtrictions.
* We will - for the moment assume that the basic integrals
* must have same spin components for particle 1 and 2.
*. This can be removed or modified later
Ctimo      IF(IPSF.EQ.IQSF) THEN
*
      DO IRSF = 1,2
      DO ISSF = 1,2
Ctimo      IF(IRSF.EQ.ISSF) THEN
        IPSMO = IMOSF_SP(IPSM,IPSF)
        IQSMO = IMOSF_SP(IQSM,IQSF)
        IRSMO = IMOSF_SP(IRSM,IRSF)
        ISSMO = IMOSF_SP(ISSM,ISSF)
        WRITE(6,'(1H ,A,4I3)')
     &  ' IPSMO IQSMO IRSMO ISSMO ',
     &    IPSMO,IQSMO,IRSMO,ISSMO
        IF(IPSMO.NE.0.AND.IQSMO.NE.0.AND.
     &  IRSMO.NE.0.AND.ISSMO.NE.0) THEN
*. Loop over all permutations of input block
          IF(I12I.EQ.0.OR.IPSM.EQ.IQSM) THEN
            J12MX = 1
          ELSE
            J12MX = 2
          END IF
*
          if(I34I.EQ.0.OR.IRSM.EQ.ISSM) then
            J34MX = 1
          ELSE
            J34MX = 2
          END IF
*
          IF(I1234I.EQ.1.OR.(IPSM.EQ.IRSM.AND.IQSM.EQ.ISSM)) THEN
            J1234MX = 1
          ELSE
            J1234MX = 2
          END IF
*
          DO J12 = 1, J12MX
            IF(J12.EQ.1) THEN
              JPSMO = IPSMO
              JQSMO = IQSMO
              J1 = 1
              J2 = 2
            ELSE
              JPSMO = IQSMO
              JQSMO = IPSMO
              J1 = 2
              J2 = 1
            END IF
            DO J34 = 1, J34MX
              IF(J34.EQ.1) THEN
               JRSMO = IRSMO
               JSSMO = ISSMO
               J3 = 3
               J4 = 4
              ELSE
               JRSMO = ISSMO
               JSSMO = IRSMO
               J3 = 4
               J4 = 3
              END IF
              DO J1234 = 1, J1234MX
                IF(J1234.EQ.1) THEN
                  KPSMO = JPSMO
                  KQSMO = JQSMO
                  KRSMO = JRSMO
                  KSSMO = JSSMO
*
                  K1 = J1
                  K2 = J2
                  K3 = J3
                  K4 = J4
*
                ELSE
                  KPSMO = JRSMO
                  KRSMO = JPSMO
                  KQSMO = JSSMO
                  KSSMO = JQSMO
*
                  K1 = J3
                  K2 = J4
                  K3 = J1
                  K4 = J2
                END IF
*
                IF(IPRINT.GE.2) THEN
                  WRITE(6,*) ' Permutation of input block'
                  WRITE(6,*) ' K1 K2 K3 K4',K1,K2,K3,K4
                  WRITE(6,*) ' KPSMO KQSMO KRSMO KSSMO ',
     &            KPSMO,KQSMO,KRSMO,KSSMO
                END IF
*
                DO JLISTO = 1, NLISTO
*. Does symmetries satisfy restrictions imposed on
* output list JLISTO
                  I_AM_OKAY = 1
                  IF(IPSMO.LT.IQSMO.AND.I12O(JLISTO).NE.0)
     &            I_AM_OKAY = 0
*. complex conjugation symmetry can be used to reduce
*. output lists to IPSMO .ge. IQSMO
                  IF(KPSMO.LT.KQSMO.AND.ICCSMO(JLISTO).NE.0)
     &            I_AM_OKAY = 0
                  IF(KRSMO.LT.KSSMO.AND.I34O(JLISTO).NE.0)
     &            I_AM_OKAY = 0
                  KPQSMO = (KPSMO-1)*NIRRIN+KQSMO
                  KRSSMO = (KRSMO-1)*NIRRIN+KSSMO
                  IF(KPQSMO.LT.KRSSMO.AND.I1234O(JLISTO).NE.0)
     &            I_AM_OKAY = 0
                  IF(I_AM_OKAY.EQ.1) THEN
*. Congratulations, you have passed all tests,
*  and you have now the right to be computed
*.
*
* One or two blocks due to real/imaginary separation
*
Ctimo                    IF(ISPITRA_RC(ISPITRA(1,JLISTO)).EQ.3.OR.
Ctimo     &                 ISPITRA_RC(ISPITRA(2,JLISTO)).EQ.3.OR.
Ctimo     &                 ISPITRA_RC(ISPITRA(3,JLISTO)).EQ.3.OR.
Ctimo     &                 ISPITRA_RC(ISPITRA(4,JLISTO)).EQ.3) THEN
*
*  Somethings wrong here. The following seems to be appropriate:
                    if (ISPITRA_RC(JLISTO).eq.3) then
                        IMIX = 1
                    ELSE
                        IMIX = 0
                    END IF
                    IF(IMIX.EQ.0) THEN
Ctimo                      NIM = 0
Ctimo                      IF(ISPITRA_RC(ISPITRA(1,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
Ctimo                      IF(ISPITRA_RC(ISPITRA(2,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
Ctimo                      IF(ISPITRA_RC(ISPITRA(3,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
Ctimo                      IF(ISPITRA_RC(ISPITRA(4,JLISTO)).EQ.1)
Ctimo     &                NIM = NIM + 1
                       NIM = 0
                       if (ISPITRA(1,JLISTO).eq.2) NIM = NIM + 1
                       if (ISPITRA(2,JLISTO).eq.2) NIM = NIM + 1
                       if (ISPITRA(3,JLISTO).eq.2) NIM = NIM + 1
                       if (ISPITRA(4,JLISTO).eq.2) NIM = NIM + 1
                    END IF
                    IF(IMIX.EQ.1) THEN
*  This transformation is complex
                      NCOMP = 2
                    ELSE
*  This transformation is real or imaginary
                      NCOMP = 1
                      IF(NIM.EQ.1.OR.NIM.EQ.3) THEN
*  This transformation is imaginary
                        IRI = 2
                      ELSE
*  This transformation is real
                        IRI = 1
                      END IF
                    END IF
                    DO JCOMP = 1, NCOMP
                      NBLOCKO = NBLOCKO + 1
                      IF(NBLOCKO.GT.MXBLKO) THEN
                        WRITE(6,*)
     &                  '  TROUBLE in INTBLK_IN_TO_OUT'
                        WRITE(6,*)
     &                  ' Too many output blocks '
                        WRITE(6,*)
     &                ' Raise MXBLKO from ', MXBLKO
                        Call Abend2( '  TROUBLE in INTBLK_IN_TO_OUT' )
                      END IF
                      IF(NCOMP.EQ.1) THEN
                        IRIO(NBLOCKO) = IRI
                      ELSE
                        IRIO(NBLOCKO) = JCOMP
                      END IF
                      LISTO(NBLOCKO) = JLISTO
*. symmetries of output block
                      IBLOCKO(1,NBLOCKO)=IPSMO
                      IBLOCKO(2,NBLOCKO)=IQSMO
                      IBLOCKO(3,NBLOCKO)=IRSMO
                      IBLOCKO(4,NBLOCKO)=ISSMO
*. Index of input block to index of output block
                      INDEXO(1,NBLOCKO) = K1
                      INDEXO(2,NBLOCKO) = K2
                      INDEXO(3,NBLOCKO) = K3
                      INDEXO(4,NBLOCKO) = K4
                    END DO
*                   /\ end  of loop over real/imag components
                  END IF
*                      (end i_am_okay check)
                END DO
*               /\ End of loop over output lists of integrals
*              given spatial symmetry
              END DO
*             /\ End of J1234 loop
            END DO
*           /\ End  of J34 Loop
          END DO
*         /\ End of J12 Loop
*
        END IF
*       /\ End if for nonvanishing symmetries
Ctimo      END IF
*     /\ End if for identity of spincomponents for particle 2
      END DO
      END DO
Ctimo      END IF
*     /\ End if for identity of spincomponents for particle 1
      END DO
      END DO
*     /\ End of loop over the spinor symmetries associated with
*     given spatial symmetry
*
*
      IF(IPRINT.GE.2) THEN
        WRITE(6,*) ' Output from INTBLK_IN_TO_OUT'
        WRITE(6,*) ' ============================'
        WRITE(6,*)
        WRITE(6,*) ' Number of output blocks', NBLOCKO
        WRITE(6,*)
        WRITE(6,*)
     &  ' Block  list   symmetries     indeces      real/imag  '
        WRITE(6,*)
     &  ' ======================================================='
        DO JBLOCKO = 1, NBLOCKO
          WRITE(6,'(1H ,I4,5x,I2,5X,4I2,5x,4I2,9X,I2)')
     &    JBLOCKO,LISTO(JBLOCKO),(IBLOCKO(I,JBLOCKO),I=1,4),
     &    (INDEXO(I,JBLOCKO),I=1,4),IRIO(JBLOCKO)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      FUNCTION IOFF_BLOCK_MAT(LROW,LCOL,IBLK)
*
* Offset of block IBLK in blocked matrix
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION LROW(*), LCOL(*)
*
      IOFF = 1
      DO JBLK = 1, IBLK-1
        IOFF = IOFF + LROW(JBLK)*LCOL(JBLK)
      END DO
*
      IOFF_BLOCK_MAT = IOFF
C?    WRITE(6,*) ' Offset for block ', IBLK, ' is ', IOFF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE IOFF_SYMBLK_MAT(NSMST,NA,NB,ITOTSM,IOFF,IRESTRICT)
*
* Offset for symmetry blocked matrix. 
* Row-symmetry used as primary index
*
* Jeppe Olsen, Summer of 99
*               July 2000 : IRESTRICT added
* Modified for relativistic case, Timo Fleig, July 2001
*
* IRESTRICT = 0 : No restrictions of symmetry blocks
*           = 1 : Only blocks with ASM .GE. BSM are included
*           =-1 : Only blocks with ASM .LT. BSM are included
#include "implicit.inc"
#include "ipoist8.inc"
*. General input  
#include "mxpdim.inc"
#include "symm.inc"
*. Specific input 
      INTEGER NA(NSMST),NB(NSMST)
*. Output
      INTEGER IOFF(NSMST)
*
      IB = 1
      DO IASM = 1, NSMST
        IBSM = IDBGMULT(INVELM(IASM),ITOTSM)
        INCLUDED = 0
        IF(IRESTRICT.EQ.0) INCLUDED = 1
        IF(IRESTRICT.EQ.1.AND.IASM.GE.IBSM) INCLUDED = 1
        IF(IRESTRICT.EQ.-1.AND.IASM.LT.IBSM) INCLUDED = -1
        IF(INCLUDED.EQ.1) THEN
          IOFF(IASM) = IB
          IB = IB + NA(IASM)*NB(IBSM)
        ELSE
          IOFF(IASM) = 0
        END IF
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE IOFF_SYMBLK_MAT_PAR(NSMST,NA,NB,ITOTSM,IOFF,IRESTRICT,
     &                               IACTLIST,IOFF_ACT,IVERS)
*
* Offset for symmetry blocked matrix. 
* Row-symmetry used as primary index
*
* Jeppe Olsen, Summer of 99
*               July 2000 : IRESTRICT added
* Modified for relativistic case, Timo Fleig, July 2001
* adaption of sequential routine, Stefan Knecht May 2007
*
* IRESTRICT = 0 : No restrictions of symmetry blocks
*           = 1 : Only blocks with ASM .GE. BSM are included
*           =-1 : Only blocks with ASM .LT. BSM are included
* IVERS     = 0 : c-block routine
*           = 1 : s-block routine
#include "implicit.inc"
#include "ipoist8.inc"
*. General input  
#include "mxpdim.inc"
#include "symm.inc"
#include "parluci.h"
*. Specific input 
      INTEGER NA(NSMST),NB(NSMST), IACTLIST(*)
*. Output
      INTEGER IOFF(NSMST)
*
      IB = 1
      IOFF_SCR = 0
      IOFF_SCR = IOFF_ACT
      IF( IVERS .eq. 0 ) THEN
        DO IASM = 1, NSMST
          IF( IACTLIST( IOFF_SCR ) .eq. 0 ) GOTO 100
          IBSM = IDBGMULT(INVELM(IASM),ITOTSM)
          INCLUDED = 0
          IF(IRESTRICT.EQ.0) INCLUDED = 1
          IF(IRESTRICT.EQ.1.AND.IASM.GE.IBSM) INCLUDED = 1
          IF(IRESTRICT.EQ.-1.AND.IASM.LT.IBSM) INCLUDED = -1
          IF(INCLUDED.EQ.1) THEN
            IOFF(IASM) = IB
            IB = IB + NA(IASM)*NB(IBSM)
          ELSE
            IOFF(IASM) = 0
          END IF
 100      CONTINUE
          IOFF_SCR = IOFF_SCR + 1
        END DO
      ELSE
        DO IASM = 1, NSMST
          IF( IACTLIST( IOFF_SCR ) .ne. MYPROC ) GOTO 200
          IBSM = IDBGMULT(INVELM(IASM),ITOTSM)
          INCLUDED = 0
          IF(IRESTRICT.EQ.0) INCLUDED = 1
          IF(IRESTRICT.EQ.1.AND.IASM.GE.IBSM) INCLUDED = 1
          IF(IRESTRICT.EQ.-1.AND.IASM.LT.IBSM) INCLUDED = -1
          IF(INCLUDED.EQ.1) THEN
            IOFF(IASM) = IB
            IB = IB + NA(IASM)*NB(IBSM)
          ELSE
            IOFF(IASM) = 0
          END IF
 200      CONTINUE
          IOFF_SCR = IOFF_SCR + 1
        END DO
      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
*   Routine for multiplying matrix (or vector) with scalar
*
      subroutine matfac(A,NROWS,NCOLUMNS,IFAC)
      implicit real*8 (A-H,O-Z)
      dimension A(NROWS,NCOLUMNS)
*
*  multiply
*
      do J=1,NROWS,1
         do I=1,NCOLUMNS,1
            A(I,J) = IFAC * A(I,J)
         end do
      end do
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*
*   Routine for multiplying two rectangular matrices A and B
*
      subroutine matmult(A,NROWSA,B,NCOLSB,C,IDIM)
      implicit real*8 (A-H,O-Z)
      dimension C(NROWSA,NCOLSB),A(NROWSA,IDIM),
     &          B(IDIM,NCOLSB)
*
*   put zeroes in matrix C
*
      do 5 M=1,NCOLSB,1
         do 7 I=1,NROWSA,1
            C(I,M) = 0.
7        continue
5     continue
      do 10 L=1,NCOLSB,1
         do 20 K=1,NROWSA,1
            do 30 J=1,IDIM,1
               C(K,L) = C(K,L) + A(K,J) * B(J,L)
30          continue
20       continue
10    continue
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Sort matrix elements according to new index pattern
*  Output matrix is resorted input matrix
*
      subroutine matsort(ZMAT,ZMATHELP,IINOUT,JINOUT,IDIM,JDIM)
*
      implicit real*8 (A-H,O-Z)
*
      dimension IINOUT(IDIM),JINOUT(JDIM)
      dimension ZMAT(IDIM,JDIM),ZMATHELP(IDIM,JDIM)
*
* IINOUT(INPUT INDEX) = OUTPUT INDEX
* JINOUT(INPUT INDEX) = OUTPUT INDEX
*
      do I=1,IDIM,1
        do J=1,JDIM,1
          if (IINOUT(I).gt.IDIM.or.JINOUT(J).gt.JDIM.
     &        or.IINOUT(I).lt.1.or.JINOUT(J).lt.1) then
            write(6,*) 'Problem in matsort'
            write(6,'(A7,I6,A3,I6)') 'IINOUT(',I,') =',IINOUT(I)
            write(6,'(A7,I6,A3,I6)') 'JINOUT(',J,') =',JINOUT(J)
            Call Abend2( 'Quitting.' )
          end if
          ZMATHELP(IINOUT(I),JINOUT(J)) = ZMAT(I,J)
        end do
      end do
*
      call copvec(ZMATHELP,ZMAT,IDIM*JDIM)
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE MERGE_VEC(LU,VEC)
*
* A vector is defined as a multiple record file on LU
* Assemble vector
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION VEC(*)
*
      NTEST = 000
*
*. Loop over records
      IBASE = 1
      IBLK = 0
 1000 CONTINUE
      CALL IFRMDS(LBLK,1,-1,LU)
      IBLK = IBLK + 1
       IF(LBLK .GE. 0 ) THEN
         CALL FRMDSC(VEC(IBASE),LBLK,-1,LU,IMZERO,IAMPACK)
         IBASE = IBASE + LBLK
       END IF
      if (NTEST.ge.100) then
        WRITE(6,*) ' IBLK,LBLK,IBASE',IBLK,LBLK,IBASE
      end if
      IF(LBLK.GE.0) GOTO 1000
*
      LENGTH = IBASE - 1
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Vector assembled from file LU = ', LU
        WRITE(6,*) ' Length = LENGTH '
        CALL WRTMAT(VEC,1,LENGTH,1,LENGTH)
      END IF
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* General routine for picking blocks of elements over
* (e.g. active) indices from a full lower triangle matrix over
* all AO symmetries. (Optional: square matrix)
*
      subroutine pick_active(NIRR,NORBTOT,NORB,TRIAIN,TRIAOUT,
     &                       RMATIN,RMATOUT,NTRIAIN,NTRIAOUT,
     &                       IFACIN,IFACOUT,IDIMIN,IDIMOUT,
     &                       IINFORM,IOUTFORM,IPRINT)
*
      implicit real*8 (A-H,O-Z)
*
      dimension NORBTOT(NIRR),NORB(NIRR)
      dimension TRIAIN(NTRIAIN),TRIAOUT(NTRIAOUT),
     &          RMATIN(IDIMIN,IDIMIN),RMATOUT(IDIMOUT,IDIMOUT)
*
*  IINFORM = 1: Input lower triangle
*            2: Input full square
*  IOUTFORM = 1: Output lower triangle
*             2: Output full square
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.100) then
        write(6,*) 'Checking input matrices in pick_active'
        if (IINFORM.eq.1) then
          write(6,*) 'Input array triangle'
          call wrtmat(TRIAIN,1,NTRIAIN,1,NTRIAIN)
        else if (IINFORM.eq.2) then
          write(6,*) 'Input array block'
          call wrtmat(RMATIN,IDIMIN,IDIMIN,IDIMIN,IDIMIN)
        else
          write(6,*) 'Problem in pick_active.'
          write(6,*) 'IINFORM has disallowed value ',IINFORM
          Call Abend2( 'Quitting.' )
        end if
      end if
*
      if (IINFORM.eq.1) then
        call unfold(RMATIN,TRIAIN,IDIMIN,1,NTRIAIN,IFACIN)
      end if
*
*  J refers to columns, I to rows
      do JRR=1,NIRR,1
        do IRR=1,NIRR,1
*
*  Input start row
*  Output start row
          IINSTA = 1
          IOUTSTA = 1
          do KRR=1,IRR-1,1
            IINSTA = IINSTA + NORBTOT(KRR)
            IOUTSTA = IOUTSTA + NORB(KRR)
          end do
          IINSTA = IINSTA + NORBTOT(IRR) - NORB(IRR)
          IINEND = IINSTA + NORB(IRR) - 1
*
*  Input start column
*  Output start column
          JINSTA = 1
          JOUTSTA = 1
          do LRR=1,JRR-1,1
            JINSTA = JINSTA + NORBTOT(LRR)
            JOUTSTA = JOUTSTA + NORB(LRR)
          end do
          JINSTA = JINSTA + NORBTOT(JRR) - NORB(JRR)
          JINEND = JINSTA + NORB(JRR) - 1
*
*  Do the job ...
          if (NTEST.ge.10) then
            write(6,*)
            write(6,*) 'IRR  ......   ',IRR
            write(6,*) 'JRR  ......   ',JRR
            write(6,*) 'IINSTA ....    ',IINSTA
            write(6,*) 'IINEND ....    ',IINEND
            write(6,*) 'JINSTA ....    ',JINSTA
            write(6,*) 'JINEND ....    ',JINEND
            write(6,*) 'IOUTSTA ...    ',IOUTSTA
            write(6,*) 'JOUTSTA ...    ',JOUTSTA
          end if
          ICOUNT = 0
          JSTART = JOUTSTA
          ISTART = IOUTSTA
          do JCOL=JINSTA,JINEND,1
            do I=IINSTA,IINEND,1
              RMATOUT(ISTART,JSTART) = RMATIN(I,JCOL)
              ISTART = ISTART + 1
              ICOUNT = ICOUNT + 1
            end do
            ISTART = IOUTSTA
            JSTART = JSTART + 1
          end do
          if (NTEST.ge.10) then
            write(6,*) 'No. of integrals transferred: '
            write(6,*) ICOUNT
          end if
        end do
      end do
*
      if (IOUTFORM.eq.1) then
        call fold(TRIAOUT,RMATOUT,IDIMOUT,1,NTRIAOUT)
        if (NTEST.ge.3) then
          write(6,*) 'Checking output from pick_active '
          do INT=1,NTRIAOUT,1
            write(6,'(I6,1F18.14)') INT,TRIAOUT(INT)
          end do
        end if
      else if (IOUTFORM.eq.2) then
        if (NTEST.ge.3) then
          write(6,*) 'Checking output from pick_active '
          call wrtmat(RMATOUT,IDIMOUT,IDIMOUT,IDIMOUT,IDIMOUT)
        end if
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
* General routine for picking blocks of elements over
* (e.g. active) indices from a full lower triangle matrix over
* all AO symmetries. (Optional: square matrix)
*
*  Array structure is given by NFCT(NINNER,NOUTER) for rows
*                       and by NFCT2(NINNER2,NOUTER2) for columns
*  IINNER gives the index of the desired block(s) for row index
*  IINNER2 for column index
*   in the array RMATIN.
*  It is assumed that the outer index runs over all present blocks.
*
      subroutine pick_active2(NOUTER,NINNER,NOUTER2,NINNER2,
     &                        NFCT,NFCT2,
     &                        IINNER,IINNER2,
     &                        TRIAIN,TRIAOUT,
     &                        RMATIN,RMATOUT,NTRIAIN,NTRIAOUT,
     &                        IFACIN,IFACOUT,IDIMIN,
     &                        IDIMOUTR,IDIMOUTC,
     &                        IINFORM,IOUTFORM,IPRINT)
*
*  Modified from pick_active for substructure arrays
*  Timo Fleig, June 2000 (it's a hot summer!)
*
      implicit real*8 (A-H,O-Z)
*
      dimension NFCT(NINNER,NOUTER),NFCT2(NINNER2,NOUTER2)
      dimension TRIAIN(NTRIAIN),TRIAOUT(NTRIAOUT),
     &          RMATIN(IDIMIN,IDIMIN),
     &          RMATOUT(IDIMOUTR,IDIMOUTC)
*
*  IINFORM = 1: Input lower triangle
*            2: Input full square
*  IOUTFORM = 1: Output lower triangle
*             2: Output full square
*
      NTESTL = 05
      NTEST = max(NTESTL,IPRINT)
*
      if (NTEST.ge.100) then
        write(6,*) 'Checking input matrices in pick_active'
        if (IINFORM.eq.1) then
          write(6,*) 'Input array triangle'
          call wrtmat(TRIAIN,1,NTRIAIN,1,NTRIAIN)
        else if (IINFORM.eq.2) then
          write(6,*) 'Input array block'
          call wrtmat(RMATIN,IDIMIN,IDIMIN,IDIMIN,IDIMIN)
        else
          write(6,*) 'Problem in pick_active.'
          write(6,*) 'IINFORM has disallowed value ',IINFORM
          Call Abend2( 'Quitting.' )
        end if
      end if
*
      if (IINFORM.eq.1) then
        call unfold(RMATIN,TRIAIN,IDIMIN,1,NTRIAIN,IFACIN)
      end if
*
*  J refers to columns, I to rows
      do JRR=1,NOUTER2,1
        do IRR=1,NOUTER,1
*
*  Input start row
*  Output start row
          IINSTA = 1
          IOUTSTA = 1
          do KRR=1,IRR-1,1
            do INN=1,NINNER,1
              IINSTA = IINSTA + NFCT(INN,KRR)
            end do
            IOUTSTA = IOUTSTA + NFCT(IINNER,KRR)
          end do
          do INN=1,IINNER-1,1
            IINSTA = IINSTA + NFCT(IINNER,IRR)
          end do
          IINEND = IINSTA + NFCT(IINNER,IRR) - 1
*
*  Input start column
*  Output start column
          JINSTA = 1
          JOUTSTA = 1
          do LRR=1,JRR-1,1
            do JNN=1,NINNER2,1
              JINSTA = JINSTA + NFCT2(JNN,LRR)
            end do
            JOUTSTA = JOUTSTA + NFCT2(IINNER2,LRR)
          end do
          do JNN=1,IINNER2-1,1
            JINSTA = JINSTA + NFCT2(IINNER2,JRR)
          end do
          JINEND = JINSTA + NFCT2(IINNER2,JRR) - 1
*
*  Do the job ...
          if (NTEST.ge.10) then
            write(6,*)
            write(6,*) 'IRR  ......   ',IRR
            write(6,*) 'JRR  ......   ',JRR
            write(6,*) 'IINSTA ....    ',IINSTA
            write(6,*) 'IINEND ....    ',IINEND
            write(6,*) 'JINSTA ....    ',JINSTA
            write(6,*) 'JINEND ....    ',JINEND
            write(6,*) 'IOUTSTA ...    ',IOUTSTA
            write(6,*) 'JOUTSTA ...    ',JOUTSTA
          end if
          ICOUNT = 0
          JSTART = JOUTSTA
          ISTART = IOUTSTA
          do JCOL=JINSTA,JINEND,1
            do I=IINSTA,IINEND,1
              RMATOUT(ISTART,JSTART) = RMATIN(I,JCOL)
              ISTART = ISTART + 1
              ICOUNT = ICOUNT + 1
            end do
            ISTART = IOUTSTA
            JSTART = JSTART + 1
          end do
          if (NTEST.ge.10) then
            write(6,*) 'No. of integrals transferred: '
            write(6,*) ICOUNT
          end if
        end do
      end do
*
      if (IOUTFORM.eq.1) then
        call fold(TRIAOUT,RMATOUT,IDIMOUT,1,NTRIAOUT)
        if (NTEST.ge.3) then
          write(6,*) 'Checking output from pick_active '
          do INT=1,NTRIAOUT,1
            write(6,'(I6,1F18.14)') INT,TRIAOUT(INT)
          end do
        end if
      else if (IOUTFORM.eq.2) then
        if (NTEST.ge.3) then
          write(6,*) 'Checking output from pick_active '
          call wrtmat(RMATOUT,IDIMOUTR,IDIMOUTC,IDIMOUTR,IDIMOUTC)
        end if
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      subroutine rdvcdc(VE,SCR,LU,LBLK,IRILP,IPRT)
*
* Read vector from disc to array
*
* Complex version      Timo Fleig, Dec 04, 2001
*
#include "implicit.inc"
*. Output
      DIMENSION VE(*),SCR(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRT)
*
      CALL REWINE(LU,LBLK)
*
*
*. Loop over blocks of file
*
      if (LBLK.gt.0) then
        LENGTH = LBLK
        call frmdsc(VE,LENGTH,-1,LU,IMZERO,IAMPACK)
        stop 'Take care: Check me first!!.'
      else
        NASS = 0
        IOFF = 1
        do IRI = 1,IRILP,1
 1000     continue
          call ifrmds(LENGTH,1,LBLK,LU)
          if (LENGTH.ge.0) then
            call frmdsc(SCR,LENGTH,LBLK,LU,IMZERO,IAMPACK)
            if (NTEST.ge.10) then
              write(6,*) ' Record read in: '
              call wrtmat(SCR,1,LENGTH,1,LENGTH)
            end if
            call copvec(SCR,VE(IOFF),LENGTH)
            IOFF = IOFF + LENGTH
          end if
          if (LENGTH.ge.0) then
            NASS = NASS + LENGTH
            goto 1000
          end if
*
          if (NTEST.ge.5) then
            if (IRI.eq.1) then
              write(6,*) ' Read vector in rdvcdc (real):  '
              call wrtmat(VE(1),1,NASS,1,NASS)
            else
              write(6,*) ' Read vector in rdvcdc (imag):  '
              call wrtmat(VE(1+NASS/2),1,NASS/2,1,NASS/2)
            end if
          end if
*
        end do
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE REO4MT(XIN,XOUT,LIN,IREO,FACIN,FACOUT)
*
* Reorder a 4 full index matrix and add to output matrix
*
* XOUT(IREO,JREO,KREO,LREO) = FACOUT*XOUT(IREO,JREO,KREO,LREO)
*                           + FACIN *XIN(I,J,K,L)
*
* Ordering of indeces is XIN(1,2,3,4)
*
* IREO : Input => output order
*
* fixed bug. IREO wasn't dimensioned ...    Timo Fleig, Aug. '98
*
      IMPLICIT REAL*8(A-H,O-Z)
*.input
      DIMENSION XIN(*)
      dimension LIN(*),IREO(*)
*. output
      DIMENSION XOUT(*)
*.local scratch
      DIMENSION ISCR(4),LOUT(4)
*. Dimension of output matrix
      DO IND = 1, 4
        LOUT(IREO(IND)) = LIN(IND)
      END DO
*. Loop over element in input order
      IELMNT_IN = 0
*
* I4
*  I3
*   I2
*    I1
*
      DO I4 = 1, LIN(4)
       DO I3 = 1, LIN(3)
        DO I2 = 1, LIN(2)
         DO I1 = 1, LIN(1)
           IELMNT_IN = IELMNT_IN + 1
*. output order
           ISCR(IREO(1)) = I1
           ISCR(IREO(2)) = I2
           ISCR(IREO(3)) = I3
           ISCR(IREO(4)) = I4
**.output address
           IELMNT_OUT = (ISCR(4)-1)*LOUT(3)*LOUT(2)*LOUT(1)
     &                + (ISCR(3)-1)        *LOUT(2)*LOUT(1)
     &                + (ISCR(2)-1)                *LOUT(1)
     &                +  ISCR(1)
           IF(FACOUT.EQ.0.0) THEN
            XOUT(IELMNT_OUT) = FACIN * XIN(IELMNT_IN)
           ELSE
            XOUT(IELMNT_OUT) = FACIN * XIN(IELMNT_IN)
     &                       + FACOUT*XOUT(IELMNT_OUT)
           END IF
         END DO
        END DO
       END DO
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE VECSUMC(C,A,B,FACAR,FACAI,FACBR,FACBI,NDIM)
*
* Multiply two complex vectors with scalars and add
*
*C(I) = FACA*A(I) + FACB*B(I)
*
* NDIM is length of real/imag part
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION C(NDIM,2)
*. Input
      DIMENSION A(NDIM,2),B(NDIM,2)
*
*
      DO I = 1, NDIM
*. Real part of C
        C(I,1) = FACAR * A(I,1)  + FACBR * B(I,1)
     &         - FACAI * A(I,2) -  FACBI * B(I,2)
*. Imaginary part of C
        C(I,2) = FACAR * A(I,2) + FACBR * B(I,2)
     &         + FACAI * A(I,1) + FACBI * B(I,1)
      END DO
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
*  Unfold lower triangle square matrix to full size;
*  I suppose that the input matrix is a one-dim. field with
*  N(N/2 + 1) elements
*  1) write vector to matrix
*  2) unfold considering symmetry (e.g. hermitian or not)
*  output matrix id partition of a fully dimensioned matrix!
*
      subroutine vectopart(ARRAY,VECTOR,ITOTDIM,ISYM,NTRIAEL,
     &                     IFAC,IPARTDIM,IROWOFF,ICOLOFF)
      implicit real*8 (A-H,O-Z)
      dimension ARRAY(ITOTDIM,ITOTDIM),VECTOR(NTRIAEL)
*
      NTEST = 000
*
*  copy to array
*
      ICOUNT = 0
      do N = 1,IPARTDIM,1
         do M = 1,N,1
            ARRAY(N+IROWOFF-1,M+ICOLOFF-1)
     &            = IFAC * VECTOR(ICOUNT+M)
         end do
         ICOUNT = ICOUNT + N
      end do
*
*  multiply with sign and generate second entry (unfold)
*
      do I=1,IPARTDIM,1
         do J=I+1,IPARTDIM,1
            ARRAY(I+IROWOFF-1,J+ICOLOFF-1)
     &            = ISYM * ARRAY(J+IROWOFF-1,I+ICOLOFF-1)
         end do
      end do
*
      if (NTEST.ge.100) then
        write(6,*) 'Info from vectopart:'
        write(6,*) 'Dimension of partition is: ',IPARTDIM
        write(6,*) 'Row offset is ',IROWOFF
        write(6,*) 'Column offset is ',ICOLOFF
        call wrtmat(ARRAY,ITOTDIM,ITOTDIM,ITOTDIM,ITOTDIM)
      end if
*
      return
      end
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE WRSVCD_REL(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,
     &                      NBLOCK,IBLOCKL,IBLOCKD,IAM_REAL,JROOT
#if defined (VAR_MPI2)
     &                      ,LU1LIST
#endif
     &                      )
*
* Write scattered vector to disc,
* every node writes only it's part.
* Information about the length of a blcok is stored in IBLOCKL.
* Information about which block belongs to which node is stored in
* IBLOCKD.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*.Input
      DIMENSION IPLAC(*),VAL(*), IBLOCKL(*), IBLOCKD(*)
*.Scratch
      DIMENSION VEC1(*)
      LOGICAL FOUND_ONE
#if defined (VAR_MPI2)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1, IOFF_SCR2
      DIMENSION LU1LIST(*)
#endif
C
      FOUND_ONE = .FALSE.
      IONE = 1
C
      IF(LBLK.GT.0) THEN
C        Write the vector without markers in one block
         CALL DZERO(VEC1,NDIM)
         DO 100 IEFF = 1, NSCAT
           VEC1(IPLAC(IEFF)) = VAL(IEFF)
  100    CONTINUE
         CALL TODSC(VEC1,NDIM,-1,LU)
      ELSE
C       Write the vector with the block format of file LUFORM
        CALL REWINE(LUFORM,-1)
        IBOT = 1
C       Loop over records
 1000   CONTINUE
#if defined (VAR_MPI2)
C       initialize offset
        IOFF_SCR1 = 0
        IOFF_SCR2 = 0
        INT_IOFF1  = 0
        INT_IOFF2  = 0
C       
        IOFF_SCR1 = MY_LU1_OFF  
     &            + MY_VEC2_IOFF * ( JROOT    - 1 ) 
     &            + MY_VEC1_IOFF * ( IAM_REAL - 1 )
        INT_IOFF1  = 1 + MY_ACT_BLK2 * ( JROOT    - 1 ) 
     &                 + MY_ACT_BLK1 * ( IAM_REAL - 1 )
CSK        WRITE(LUWRT,*) 'OFFSET at START',IOFF_SCR1
CSK        WRITE(LUWRT,*) 'INT_IOFF1 at START ',INT_IOFF1
C
        DO IBLK = 1, NBLOCK
          IF( IBLOCKD( IBLK ) .ne. MYPROC ) GOTO 300
CSK          WRITE(LUWRT,*) 'OFFSET new for block',IOFF_SCR1,IBLK -
CSK     &                                          INT_IOFF2, INT_IOFF2
CSK          WRITE(LUWRT,*) 'INT_IOFF1 new for block ',INT_IOFF1,IBLK -
CSK     &                                              INT_IOFF2, INT_IOFF2
          LBL = IBLOCKL( IBLK ) 
C         
#else
        CALL IFRMDS(LBL,1,-1,LUFORM)
        CALL ITODS(LBL,1,-1,LU)
        IF(LBL.GE.0) THEN
#endif
C?        write(6,*) ' IBOT,IBOT + LBL -1 ',IBOT,
C?   &    IBOT + LBL -1
          CALL DZERO(VEC1,LBL)
C         *******************************************************
C         *******************************************************
C                Find the places where a 1.0 should be set
C         *******************************************************
C         *******************************************************
          DO 200 IEFF = 1, NSCAT
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1) THEN
               VEC1(IPLAC(IEFF)-IBOT+1) = VAL(IEFF)
               FOUND_ONE = .TRUE.
            END IF
CSK          IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
CSK     &      write(luwrt,*) ' Catch : IPLAC(IEFF) VAL(IEFF) ',
CSK     &      IPLAC(IEFF),VAL(IEFF)
  200     CONTINUE
C         
#if defined (VAR_MPI2)
C
C         new offset
C
          IOFF_SCR1 = IOFF_SCR1 + IOFF_SCR2
          INT_IOFF1 = INT_IOFF1 + INT_IOFF2
C
C
          IF( FOUND_ONE ) THEN
CSK          WRITE(LUWRT,*) 'THIS IS WHAT I WILL WRITE',MYPROC
CSK          CALL WRTMATMN(VEC1,1,LBL,1,LBL,LUWRT)
CSK          WRITE(LUWRT,*) 'OFFSET',IOFF_SCR1
CSK          WRITE(LUWRT,*) 'INT_IOFF1 is ',INT_IOFF1
C
            LU1LIST( INT_IOFF1 ) = IONE
            CALL MPI_FILE_WRITE_AT(ILU1,IOFF_SCR1,VEC1,LBL,
     &                             MPI_REAL8,ISTAT,IERR)
C
          END IF
C
          FOUND_ONE = .FALSE.
#else
          CALL TODSC(VEC1,LBL,-1,LU)
#endif
C         Skip record on LUFORM
#if !defined (VAR_MPI2)
          CALL FRMDSC(VEC1,LBL,-1,LUFORM,IMZERO,IAMPACK)
          IBOT = IBOT + LBL
          GOTO 1000
        END IF
#else
C         keep track of correct offset
          IOFF_SCR2 = LBL
          INT_IOFF2 = IONE
  300     CONTINUE
C         keep IBOT up-to-date for all nodes
          LBL = IBLOCKL( IBLK )
          IBOT = IBOT + LBL
        END DO
#endif
C
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE WRTVCDC(SEGMNT,LU,IREW,IRC,LBLK)
C
C PRINT COMPLEX VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU
        ELSE
          CALL REWINE(LU,LBLK)
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
      DO IRI = 1,IRC,1
        IF(IRI.EQ.1) THEN
           WRITE(6,*) ' Real Part '
        ELSE
           WRITE(6,*) ' Imaginary part '
        END IF
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU) LBL
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,I_AM_PACKED)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,6)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      END DO
*. End of loop over real/imaginary parts
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE ZERORC_REL(MBLOCK,IFIL)
*
* A record was known to be identical  zero
*
*
* Write corresponding info to file IFIL
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISCR(2)
* Zero record
      ISCR(1) = 1
!  no packing
      ISCR(2) = 0
*
      CALL ITODS(ISCR,2,2,IFIL)
*
      END
