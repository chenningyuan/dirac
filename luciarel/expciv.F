      SUBROUTINE EXPCIV(ISM,ISPCIN,LUIN,
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
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'EXPCIV')
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
      CALL MEMMAN(KLIABI,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABI')
      CALL MEMMAN(KLIABU,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABU')
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
        CALL MEMMAN(KLBLIN,NSMST,'ADDL  ',1,'KLBLIN')
        CALL MEMMAN(KLBLUT,NSMST,'ADDL  ',1,'KLBLUT')
        CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLIN),IDUMMY)
        CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLUT),IDUMMY)

*
*. Number of dets etc per TTS block
*
        CALL MEMMAN(KLTTSII,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
        CALL MEMMAN(KLTTSNI,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLIN),NSMST,WORK(KLIABI),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSII),WORK(KLTTSNI),
     &            NCOMBI,0)
*
        CALL MEMMAN(KLTTSIU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
        CALL MEMMAN(KLTTSNU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
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
      CALL MEMMAN(KLBLI,LENGTHI,'ADDL  ',2,'KLBLI ')
      CALL MEMMAN(KLBLU,LENGTHU,'ADDL  ',2,'KLBLU ')
*
*. and now : Let another subroutine complete the taks
*
      CALL REWINO(LUIN)
      CALL REWINO(LUUT)
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
C?    CALL REWINO(LUIN)
*
      write(6,*) ' ICISTR = ', ICISTR
      DO IROOT = 1, NROOT
*. Input vector should be first vector on file so
        IF(IROOT.EQ.1) THEN
          LLUIN = LUIN
        ELSE
          IF(ICISTR.EQ.1) THEN
            CALL REWINO(LUSCR)
            CALL FRMDSC(WORK(KLBLI),NCOMBI,-1,LUIN,IMZERO,IAMPACK)
            CALL  TODSC(WORK(KLBLI),NCOMBI,-1,LUSCR)
            CALL REWINO(LUSCR)
            LLUIN = LUSCR
          ELSE
*. With the elegance of an elephant
            CALL REWINO(LUSCR)
            CALL REWINO(LUIN)
            DO JROOT = 1, IROOT
              CALL REWINO(LUSCR)
              CALL COPVCD(LUIN,LUSCR,WORK(KLBLI),0,-1)
            END DO
            CALL REWINO(LUSCR)
C?          write(6,*) ' Vector copied to LUSCR '
C?          CALL WRTVCD(WORK(KLBLI),LUSCR,1,-1)
C?          CALL REWINO(LUSCR)
            LLUIN = LUSCR
          END IF
        END IF
*
        CALL EXPCIVS(LLUIN,WORK(KLBLI),NCOMBI,WORK(KLTTSNI),
     &       WORK(KLTTSII),
     &       WORK(KLIABI),NOCTPA,NOCTPB,WORK(KLBLIN),
     &       LUUT,WORK(KLBLU),NCOMBU,WORK(KLTTSNU),
     &       WORK(KLTTSIU),WORK(KLIABU),
     &       WORK(KLBLUT),
     &       ICISTR,IDC,NSMST,
     &       LBLK)
*
      END DO
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINO(LUIN)
        CALL REWINO(LUUT)
        DO IROOT = 1, NROOT
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC(WORK(KLBLU),NCOMBU,-1,LUUT,IMZERO,IAMPACK)
            CALL  TODSC(WORK(KLBLU),NCOMBU,-1,LUIN)
          ELSE
            CALL COPVCD(LUUT,LUIN,WORK(KLBLU),0,-1)
          END IF
        END DO
      END IF
*
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
*
C?    CALL REWINO(LUUT)
C?    DO IROOT = 1, NROOT
C?    WRITE(6,*) ' Root number ', IROOT
*
C?      IF(ICISTR.EQ.1) THEN
C?          CALL WRTMAT(WORK(KLBLU),1,NCOMBUT,1,NCOMBUT)
C?      ELSE
C?        CALL WRTVCD(WORK(KLBLU),LUUT,0,-1)
C?      END IF
C?    END DO
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'EXPCV2')
*
      RETURN
      END
