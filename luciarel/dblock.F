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
*
      integer ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      common/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      common/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
      common/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      IDUM = 0
      call memman(IDUM,IDUM,'MARK  ',IDUM,'DBLOCK')
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
      call memman(KCONSPA_S,NOCTPA_S**2,'ADDL  ',1,'CNSPAC')
      call memman(KCONSPB_S,NOCTPB_S**2,'ADDL  ',1,'CNSPBC')
C     SPGRPCON(IOFSPGRP,NSPGRP,NGAS,MXPNGAS,IELFSPGRP,ISPGRPCON,IPRNT)
*
      call SPGRPCON(IOCTPA_S,NOCTPA_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA_S),IPRCIX)
      call SPGRPCON(IOCTPB_S,NOCTPB_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB_S),IPRCIX)
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      call memman(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      call memman(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      call STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
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
        LSCR1 = MXSB
      else if (ICISTR.eq.3) then
        LSCR1 = MXSOOB
      end if
      LSCR1 = max(LSCR1,LCSBLK)
      if (IPRCIX.ge.2)
     &write(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
*. Arrays giving allowed type combinations
      call memman(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',2,'CIOIO ')
      call memman(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
*. sigma needed for MXRESC
      call IAIBCM(ISSPC,IATP,IBTP,WORK(KSIOIO))
      call IAIBCM(ICSPC,JATP,JBTP,WORK(KCIOIO))
*. Arrays giving block type
      call memman(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
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
      call memman(KRHO1S,LRHO1,'ADDL  ',2,'RHO1S ')
CTF
*. Allocation of memory for annihilation/creation from supergroups
CTF setting this undefined variable to two:
      MXSXST = 2
      MXADKBLK = 2
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0,
     &            MXSTBL*MXTSOB_L)
      call memman(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      call memman(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      call memman(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      call memman(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      call memman(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      call memman(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      call memman(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      call memman(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
*. Type of blocks (all types set to 1 for now )
      call ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
*.Some TTS arrays
      NTTS = MXNTTS
*. for partitioning of vector
      call memman(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      call memman(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      call memman(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      call memman(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      call memman(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
*. Space for four blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (max(NAEL_S,NAEL_C,NBEL_S,NBEL_C)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (max(NAEL_S,NBEL_S,NAEL_C,NBEL_C)) * NOCOB
      do I1234=1,4,1
        call memman(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
        call memman(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        call memman(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      end do
      call memman(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
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
      call memman(IDUM,IDUM,'FLUSM ',2,'DBLOCK')
      return
      end
