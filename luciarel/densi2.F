      SUBROUTINE DENSI2(I12,RHO1,
CTF  &                  RHO2,
     &                  L,R,LUL,LUR)
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
C     COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC,IASTP,IBSTP,IACTP,IBCTP
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
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
*./LUCINP/
*. Used : NSMOB
#include "lucinp.inc"
#include "clunit.inc"
*. Scratch for string information
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*. Specific input
      REAL*8 L
      DIMENSION L(*),R(*)
*.Output
      DIMENSION RHO1(*)
CTF   dimension RHO2(*)
*. Before I forget it :
      CALL QENTER('DENSI')
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DENSI ')
      ZERO = 0.0D0
      CALL SETVEC(RHO1,ZERO,NACOB ** 2)
      IF(I12.EQ.2)
     &CALL SETVEC(RHO2,ZERO,NACOB ** 2 *(NACOB**2+1)/2)
*
      WRITE(6,*) ' ISSPC ICSPC in DENSI2 ',ISSPC,ICSPC
      WRITE(6,*) ' IASTP IBSTP in DENSI2 ',IASTP,IBSTP
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
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
      CALL MEMMAN(KCONSPA,NOCTPA**2,'ADDL  ',1,'CONSPA')
      CALL MEMMAN(KCONSPB,NOCTPB**2,'ADDL  ',1,'CONSPB')
      CALL SPGRPCON(IOCTPA,NOCTPA,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPA),IPRCIX)
      CALL SPGRPCON(IOCTPB,NOCTPB,NGAS,MXPNGAS,NELFSPGP,
     &              WORK(KCONSPB),IPRCIX)
*. Largest block of strings in zero order space
      MAXA0 = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTYP(IATP),2)
      MAXB0 = IMNMX(WORK(KNSTSO(IBTP)),NSMST*NOCTYP(IBTP),2)
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
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
        LSCR1 = MXSB
      ELSE IF(ICISTR.EQ.3) THEN
        LSCR1 = MXSOOB
      END IF
      LSCR1 = MAX(LSCR1,LCSBLK)
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' ICISTR,LSCR1 ',ICISTR,LSCR1
      IF(ICISTR.EQ.1) THEN
        CALL MEMMAN(KCB,LSCR1,'ADDL  ',2,'KCB   ')
        CALL MEMMAN(KSB,LSCR1,'ADDL  ',2,'KSB   ')
      END IF
*.SCRATCH space for block of two-electron density matrix
* A 4 index block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRDEN.GE.2)
     &WRITE(6,*) ' Density scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
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
      CALL MXRESC(WORK(KCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
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
COLD  CALL MEMMAN(KC2,LSCR12,'ADDL  ',2,'KC2   ')
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR12
      KSSCR = KC2
      KCSCR = KC2 + LSCR2
*
*. Space for annihilation/creation mappings
      MAXIK = MAX(MAXI,MAXK)
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      CALL MEMMAN(KI1,  LSCR3       ,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KI2,  LSCR3       ,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KI3,  LSCR3       ,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KI4,  LSCR3       ,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
      CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
      CALL MEMMAN(KXI3S,LSCR3       ,'ADDL  ',2,'XI3S  ')
      CALL MEMMAN(KXI4S,LSCR3       ,'ADDL  ',2,'XI4S  ')
*. Arrays giving block type
      CALL MEMMAN(KSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
      CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KSBLTP),WORK(KSVST))
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.0 OOS arrayy
      NOOS = NOCTPA*NOCTPB*NSMST
* scratch space containing active one body
      CALL MEMMAN(KRHO1S,NACOB ** 2,'ADDL  ',2,'RHO1S ')
*. For natural orbitals
      CALL MEMMAN(KRHO1P,NACOB*(NACOB+1)/2,'ADDL  ',2,'RHO1P ')
      CALL MEMMAN(KXNATO,NACOB **2,'ADDL  ',2,'XNATO ')
*. Natural orbitals in symmetry blocks
      CALL MEMMAN(KRHO1SM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KXNATSM,NACOB ** 2,'ADDL  ',2,'RHO1S ')
      CALL MEMMAN(KOCCSM,NACOB ,'ADDL  ',2,'RHO1S ')
*
*. Space for two  blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*. Space for two blocks of string occupations and four arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      DO K12 = 1, 2
        CALL MEMMAN(KLOCSTR(K12),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
      END DO
      DO I1234 = 1, 4
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
*. Arrays for partitioning of Left vector = sigma
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTL ,NTTS  ,'ADDL  ',1,'LBT_L  ')
      CALL MEMMAN(KLLEBTL,NTTS  ,'ADDL  ',1,'LEBT_L ')
      CALL MEMMAN(KLI1BTL,NTTS  ,'ADDL  ',1,'I1BT_L ')
      CALL MEMMAN(KLIBTL ,8*NTTS,'ADDL  ',1,'IBT_L  ')
      CALL MEMMAN(KLSCLFCL,NTTS, 'ADDL  ',2,'SCLF_L')
      CALL PART_CIV2(IDC,WORK(KSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KSIOIO),ISMOST(1,ISSM),
     &     NBATCHL,WORK(KLLBTL),WORK(KLLEBTL),
     &     WORK(KLI1BTL),WORK(KLIBTL),0)
*. Number of BLOCKS
        NBLOCKL = IFRMR(WORK(KLI1BTL),1,NBATCHL)
     &         + IFRMR(WORK(KLLBTL),1,NBATCHL) - 1
*. Arrays for partitioning of Right  vector = C
      NTTS = MXNTTS
      CALL MEMMAN(KLLBTR ,NTTS  ,'ADDL  ',1,'LBT_R  ')
      CALL MEMMAN(KLLEBTR,NTTS  ,'ADDL  ',1,'LEBT_R ')
      CALL MEMMAN(KLI1BTR,NTTS  ,'ADDL  ',1,'I1BT_R ')
      CALL MEMMAN(KLIBTR ,8*NTTS,'ADDL  ',1,'IBT_R  ')
      CALL MEMMAN(KLSCLFCR,NTTS, 'ADDL  ',2,'SCLF_R')
      CALL PART_CIV2(IDC,WORK(KCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,NSMST,LSCR1,
     &     WORK(KCIOIO),ISMOST(1,ICSM),
     &     NBATCHR,WORK(KLLBTR),WORK(KLLEBTR),
     &     WORK(KLI1BTR),WORK(KLIBTR),0)
*. Number of BLOCKS
        NBLOCKR = IFRMR(WORK(KLI1BTR),1,NBATCHR)
     &         + IFRMR(WORK(KLLBTR),1,NBATCHR) - 1
        WRITE(6,*) ' DENSI2T :NBLOCKR =',NBLOCKR
      WRITE(6,*) ' Memcheck in densi2 before GASDN2'
      CALL LMEMCHK

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
      CALL LMEMCHK
*
* Natural Orbitals etc.
*
C     CALL NATORB(RHO1,NSMOB,NTOOBS,NACOBS,NINOBS,
C    &            IREOST,WORK(KXNATO),
C    &            WORK(KRHO1SM),WORK(KOCCSM),NACOB,
C    &            WORK(KRHO1P),IPRDEN)
      IF(IPRDEN.GE.5.AND.I12.EQ.2) THEN
        WRITE(6,*) ' Two-electron density '
        CALL PRSYM(RHO2,NACOB**2)
      END IF
*. Eliminate local memory
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'DENSI ')
      CALL QEXIT('DENSI')
C     WRITE(6,*) ' Leaving DENSI '
      RETURN
      END
