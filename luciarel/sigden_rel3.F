      SUBROUTINE SIGDEN_REL3(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,T,ISIGDEN,
     &                       IRIS,IRIC)
*
* Generate a set of sigma blocks, 
* The NBLOCK specified in IBLOCK starting from IBOFF,
* be more specific.
*
* The blocks are delivered in HCB
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*. Sigma blocks require
      INTEGER IBLOCK(8,*)
*
#include "cands.inc"
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
#include "ctcc.inc"
      COMMON/CMXCJ/MXCJ
*
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
#include "cintfo.inc"
*
      DIMENSION CB(*), HCB(*)
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SIGDE3')
*
      NTESTL = 00
      IF(NTESTL.GE.5)
     &WRITE(6,*) ' SIGDE3 : ISSPC,ICSPC ', ISSPC,ICSPC
*. MXINKA is pt used for length of batch of strings and Toccs
*
*.Local scratch arrays for blocks of C and sigma
      LSCR1 = MXSOOB_AS
      LSCR1 = MAX(LSCR1,LCSBLK)
*.T-elements                       
      INTSCR = MIN(MXINKA**4,MX_TBLK_MX)
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*. Offsets for alpha and beta supergroups
*
      NOCTPA_C = NOCTYP(IACTP)
      NOCTPB_C = NOCTYP(IBCTP)
*
      NAEL_C = NELEC(IACTP)
      NBEL_C = NELEC(IBCTP)
*
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KCIOIO,NOCTPA_C*NOCTPB_C,'ADDL  ',2,'CIOIO ')
      CALL IAIBCM(ICSPC,IACTP,IBCTP,WORK(KCIOIO))
*
      KC2 = KVEC3
*
      KCJRES = KC2
      KSIRES = KC2 + MXINKA**4
*
*. Arrays for storing NEL consecutive annihilations/creations
*  (resolution matrices)
* Number of alpha/beta excitation types
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                      1,NAOBEX_TP,NBOBEX_TP,IDUMMY,IDUMMY)
*. And the alpha-and beta-excitations
      LENA = 2*NGAS*NAOBEX_TP
      LENB = 2*NGAS*NBOBEX_TP
      CALL MEMMAN(KLAOBEX,LENA,'ADDL  ',2,'IAOBEX')
      CALL MEMMAN(KLBOBEX,LENB,'ADDL  ',2,'IAOBEX')
      CALL SPOBEX_TO_ABOBEX(WORK(KLSOBEX),NSPOBEX_TP,NGAS,
     &                      0,NAOBEX_TP,NBOBEX_TP,
     &                      WORK(KLAOBEX),WORK(KLBOBEX))
*
*. For alpha excitations, C and S vector (might be flipped around!)
      IOCTPA_C = IBSPGPFTP(IACTP)
      CALL LEN_GENOP_STR_MAP(NAOBEX_TP,WORK(KLAOBEX),NOCTPA_C,
     &                       NELFSPGP(1,IOCTPA_C),NOBPT,NGAS,
     &                       MAXLENAC)
      IOCTPA_S = IBSPGPFTP(IASTP)
      NOCTPA_S = NOCTYP(IASTP)
      CALL LEN_GENOP_STR_MAP(NAOBEX_TP,WORK(KLAOBEX),NOCTPA_S,
     &                       NELFSPGP(1,IOCTPA_S),NOBPT,NGAS,
     &                       MAXLENAS)
*
*. For beta excitations
      IOCTPB_C = IBSPGPFTP(IBCTP)
      CALL LEN_GENOP_STR_MAP(NBOBEX_TP,WORK(KLBOBEX),NOCTPB_C,
     &                       NELFSPGP(1,IOCTPB_C),NOBPT,NGAS,
     &                       MAXLENBC)
      IOCTPB_S = IBSPGPFTP(IBSTP)
      NOCTPB_S = NOCTYP(IBSTP)
      CALL LEN_GENOP_STR_MAP(NBOBEX_TP,WORK(KLBOBEX),NOCTPB_S,
     &                       NELFSPGP(1,IOCTPB_S),NOBPT,NGAS,
     &                       MAXLENBS)
*
      MAXLEN_I1 = MAX(MAXLENAC,MAXLENBC,MAXLENAS,MAXLENBS)
C?    WRITE(6,*) ' MAXLEN_I1 = ', MAXLEN_I1
*
      LSCR3 = MAXLEN_I1
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
*
*. Arrays for storing occupations of T-operators
      CALL MEMMAN(KTOCC1,MX_ST_TSOSO_BLK_MX,'ADDL  ',2,'TOCC1 ')
      CALL MEMMAN(KTOCC2,MX_ST_TSOSO_BLK_MX,'ADDL  ',2,'TOCC2 ')
      CALL MEMMAN(KTOCC3,MX_ST_TSOSO_BLK_MX,'ADDL  ',2,'TOCC3 ')
      CALL MEMMAN(KTOCC4,MX_ST_TSOSO_BLK_MX,'ADDL  ',2,'TOCC4 ')
*
*. Arrays for storing NEL consecutive annihilations/creations
      LSCR4 = MXINKA*MXINKA*MXINKA
      CALL MEMMAN(KI1G  ,LSCR4,'ADDL  ',1,'I1G   ')
      CALL MEMMAN(KXI1G,LSCR4,'ADDL  ',2,'XIG  ')
*
      CALL MEMMAN(KI2G  ,LSCR4,'ADDL  ',1,'I2G   ')
      CALL MEMMAN(KXI2G,LSCR4,'ADDL  ',2,'XI2G  ')
*
      CALL MEMMAN(KI3G  ,LSCR4,'ADDL  ',1,'I3G   ')
      CALL MEMMAN(KXI3G,LSCR4,'ADDL  ',2,'XI3G  ')
*
      CALL MEMMAN(KI4G  ,LSCR4,'ADDL  ',1,'I4G   ')
      CALL MEMMAN(KXI4G,LSCR4,'ADDL  ',2,'XI4G  ')
*
      CALL MEMMAN(KI1GE  ,LSCR4,'ADDL  ',1,'I1GE  ')
      CALL MEMMAN(KXI1GE,LSCR4,'ADDL  ',2,'XIGE  ')
*
      CALL MEMMAN(KI2GE  ,LSCR4,'ADDL  ',1,'I2GE  ')
      CALL MEMMAN(KXI2GE,LSCR4,'ADDL  ',2,'XI2GE ')
*
      CALL MEMMAN(KI3GE  ,LSCR4,'ADDL  ',1,'I3GE  ')
      CALL MEMMAN(KXI3GE,LSCR4,'ADDL  ',2,'XI3GE ')
*
      CALL MEMMAN(KI4GE  ,LSCR4,'ADDL  ',1,'I4GE  ')
      CALL MEMMAN(KXI4GE,LSCR4,'ADDL  ',2,'XI4GE ')
*.TTS arrays  for partitioning of vector 
      NTTS = MXNTTS
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,2*NTTS,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
*. For scaling for each TTS block
      CALL MEMMAN(KLSCLFAC ,NTTS,'ADDL  ',2,'SCLFAC')
*. Arrays giving block type
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Find batches of C - strings
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP))
*. Use all symmetryblocks of given TT
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IACTP)),
     &     WORK(KNSTSO2(IBCTP)),NOCTPA_C,NOCTPB_C,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),NCBATCH,WORK(KLLBT),
     &     WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     0,ITTSS_ORD)
*. Space for four blocks of string occupations and arrays of 
*. reordering arrays
*. Also used to hold an NORB*NORB matrix  
      LZSCR = (MAX(NAEL_C,NBEL_C)+3)*(NOCOB+1) 
     &         + 2 * NOCOB + NOCOB*NOCOB
      LZ    = (MAX(NAEL_C,NBEL_C)+2) * NOCOB
*. Set up to two blocks for orbital conserving operator 
      DO I1234 = 1, 2
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
*
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',2,'KLZSCR')
      call memman(KLSPOBEX_FRZ,NSPOBEX_TP,'ADDL  ',1,'SPOBFR')
*
      CALL SIGDEN_REL4(NBLOCK,IBLOCK(1,IBOFF),CB,HCB,WORK(KVEC3),
     &     WORK(KNSTSO(IASTP)),WORK(KNSTSO2(IBSTP)),
     &     WORK(KNSTSO(IACTP)),WORK(KNSTSO2(IBCTP)),
     &     IASTP,IBSTP,IACTP,IBCTP,MXINKA,WORK(KINSCR),
     &     WORK(KI1),WORK(KXI1S),WORK(KI2),WORK(KXI2S),
     &     WORK(KI3),WORK(KXI3S),WORK(KI4),WORK(KXI4S),
     &     IPRSIG,IPRDEN,LUC,WORK(KCJRES),WORK(KSIRES),
     &     WORK(KLLBT),WORK(KLLEBT),WORK(KLI1BT),WORK(KLIBT),
     &     ISSM,ICSM,ISIGDEN,NCBATCH,T,
     &     IDIAG,IRIS,IRIC,
     &     NSPOBEX_TP,WORK(KLSOBEX),WORK(KLIBSOBEX),
     &     WORK(KLSPOBEX_AC),
     &     WORK(KTOCC1),WORK(KTOCC2), WORK(KTOCC3),WORK(KTOCC4),
     &     WORK(KI1G),WORK(KXI1G), WORK(KI2G),WORK(KXI2G),
     &     WORK(KI3G),WORK(KXI3G), WORK(KI4G),WORK(KXI4G),
     &     WORK(KI1GE),WORK(KXI1GE), WORK(KI2GE),WORK(KXI2GE),
     &     WORK(KI3GE),WORK(KXI3GE), WORK(KI4GE),WORK(KXI4GE),
     &     WORK(KSIGNNHX),WORK(KLSPOBEX_FRZ))
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SIGDE3')
      RETURN
      END
