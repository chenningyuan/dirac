      SUBROUTINE SBLOCK_REL(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
     &                      LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END,
     &                      IRIS,IRIC,I_DO_TRNSP)
*
* Generate a set of sigma blocks,
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order
* If LUCBLK.GT.0, the blocks of C corresponding to IBLOCK
* are stored on LUCBLK
*
* ICONSPA,ICONSPB  added October 1996
* ICBAT_RES, ICBAT_INI, IBBAT_END added august 1997
*
* If ICBAT_RES .eq.1 then it as assumed that only
* Cbatches ICBAT_INI to ICBAT_END are stored on  LUC
*
* Generalized to NOPAIR Hamiltonian, T.F + J.O., March21, 7.22
*                                    (not bad for a saturday)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
* IA/BS/CTP : Type (number of elecs in A/B,S/C)
#include "cands.inc"
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
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
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRDR/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCRR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
CTF
* Scratch space for rsbb2a, rsbb2bn, rsbb2bvn, skickj
      dimension SOMESCR(MXTSOB*MXTSOB*MXTSOB*MXTSOB)
      dimension XIJILS(MXTSOB)
*
      IDUM = 0
      call memmar(IDUM,IDUM,'MARK  ',IDUM,'SBLOCK')
*
      NTEST = 0
      IF(NTEST.GE.5)
     &WRITE(6,*) ' SBLOCK : ISSPC,ICSPC ', ISSPC,ICSPC
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
      IF(LUCBLK.GT.0) THEN
        CALL REWINO(LUCBLK)
      END IF
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
      CALL SPGRPCON(IOCTPA_S,NOCTPA_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA_S),IPRCIX)
      CALL SPGRPCON(IOCTPB_S,NOCTPB_S,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB_S),IPRCIX)
*
*
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      call memmar(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      call memmar(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM_REL(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MXSTBL = MXNSTR
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN(MXINKA,MXSTBL)
      MAXK = MIN(MXINKA,MXSTBL)
      MAXIK = MAX(MAXI,MAXK)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
      DO IOBSM = 1, NSMOB
       MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
      END DO
      END DO
      MAXIJ = MXTSOB_L ** 2
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      call memmar(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*.Local scratch arrays for blocks of C and sigma
      IF(ICISTR.LE.2) THEN
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      if (IPRCIX.ge.2) then
        WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      end if
*. Arrays giving allowed type combinations
      call memmar(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',2,'CIOIO ')
      call memmar(KSIOIO,NOCTPA_S*NOCTPB_S,'ADDL  ',2,'SIOIO ')
*. sigma needed for MXRESC
C          IAIBCM_REL(ICISPC,IAIB)
      CALL IAIBCM_REL(ISSPC,IATP,IBTP,WORK(KSIOIO))
      CALL IAIBCM_REL(ICSPC,JATP,JBTP,WORK(KCIOIO))
*. Arrays giving block type
      call memmar(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. for SIgmav reflection - not active
CTF   KSVST = 1
*  Why the heck is WORK(1) used here??
*  This is damn dangerous. ISTRFL allocated statically in sblocks
*
*.scratch space for projected matrices and a CI block
      LSCR2 = LCJRES
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
C  memory for blocks already allocated, so
      KC2 = KVEC3
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*. Allocation of memory for annihilation/creation from supergroups
CTF setting this undefined variable to two:
      MXSXST = 2
      MXADKBLK = 2
*. The size of XI* should be obtained from MXRESCP, but for now ..
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0,
     &             MXSTBL*MXTSOB_L)
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
      CALL ZBLTP_REL(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
*.Some TTS arrays
      NTTS = MXNTTS
*. for partitioning of vector
      call memmar(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      call memmar(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      call memmar(KLI1BT,2*NTTS  ,'ADDL  ',1,'I1BTC ')
      call memmar(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      call memmar(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
*. Space for four blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL_S,NAEL_C,NBEL_S,NBEL_C)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL_S,NBEL_S,NAEL_C,NBEL_C)) * NOCOB
      DO I1234 = 1, 2
      call memmar(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      call memmar(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
      call memmar(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      call memmar(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
* Fed in directly from common block oper.inc.
* Set in readin !!
      IF(IHAM12.EQ.2) THEN
* One- and two-electron Hamiltonian active
        IDOH2 = 1
      ELSE
* Only one-electron Hamiltonian active
        IDOH2 = 0
      END IF
*
      CALL SBLOCKS_REL(WORK(KPICI1),
     &            NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KC2),
     &            WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
     &            NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            WORK(KNSTSO(JATP)),WORK(KNSTSO2(JBTP)),
     &            NAEL_S,IATP,NBEL_S,IBTP,
     &            NAEL_C,JATP,NBEL_C,JBTP,
     &            IOCTPA_S,IOCTPB_S,NOCTPA_S,NOCTPB_S,
     &            IOCTPA_C,IOCTPB_C,NOCTPA_C,NOCTPB_C,
     &            NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,
     &            NOBPTS2,IOBPTS2,MXPNGAS,
     &            ITSOB,MAXIJ,MAXK,MAXI,INSCR,LSCR1,
     &            LSCR1,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
     &            SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
     &            ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
     &            MXTSOB,SOMESCR,XIJILS,
     &            WORK(KI1),WORK(KXI1S),
     &            WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,
     &            PSSIGN,IPRSIG,LUC,ICJKAIB,WORK(KCJRES),
     &            WORK(KSIRES),WORK(KI3),WORK(KXI3S),
     &            WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
     &            MOCAA,MOCAB,IAPR,
     &            WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &            IRESTRICT,WORK(KCONSPA_S),WORK(KCONSPB_S),
     &            WORK(KLSCLFAC),
     &            LUCBLK,IPERTOP,IH0INSPC,ISPINFREE,
CTI  &            WORK(KLH0SPC),
     &            ICBAT_RES,ICBAT_INI,ICBAT_END,
     &            IRELACI,IRIS,IRIC,I_DO_TRNSP,IELCI)
*
C?    write(6,*) 'Checking IDC  ',IDC
      IF(IDC.EQ.2) THEN
*. reform
        CALL RFTTS(HCB,CB,IBLOCK(1,IBOFF),NBLOCK,
     &             1,NSMST,NOCTPA_S,NOCTPB_S,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO2(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,IBLOCK(1,IBOFF),NBLOCK,NSMST,NOCTPA_S,NOCTPB_S,
     &              WORK(KNSTSO(IATP_S)), WORK(KNSTSO2(IBTP_S)),
     &              IDC,1,NTEST)
      END IF
*
      IF(LUCBLK.GT.0) THEN
        CALL ITODS(-1,1,-1,LUCBLK)
      END IF
*. Eliminate local memory
      IDUM = 0
      call memmar(IDUM ,IDUM,'FLUSM ',2,'SBLOCK')
      RETURN
      END
