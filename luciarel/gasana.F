      SUBROUTINE GASANA(C,NBLOCK,IBLOCK,IBLTP,LUC,LUSCR,
     &                  ICISTR)
*
*
*
* Analyze CI vector
*
* Jeppe Olsen, August 1995
* Driven By IBLOCK, May 1997
*           String occupations added, Feb. 98
*
*
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
* Modified for relativistic CI calculations
*
*        Timo Fleig, Feb. 2000
*´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´
*
      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
#include "mxpdim.inc"
#include "orbinp.inc"
#include "strbas.inc"
#include "cicisp.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "cprnt.inc"
#include "cands.inc"
*
      DIMENSION IBLOCK(8,NBLOCK),IBLTP(*)
      dimension C(*)
*
      CALL QENTER('ANACI')
      CALL MEMMAN(KLOFF,DUMMY,'MARK  ',DUMMY,'GASANA')
*
** Specifications of internal space
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRWFC)
*
      LBLK = -1
      IBOFF = 0
      call rewino(LUSCR)
      call itods(-1,1,-1,LUSCR)
* Type of alpha and beta strings
      do IRIR = 1,2,1
        if (NTEST.ge.50) write(6,*) 'right r/i is  ',IRIR
        if (IRIR.EQ.2) then
*. Skip EOV mark between real and imaginary part
          call IFRMDS(IONEM,1,-1,LUC)
        end if
        do IMK2_R = 1,NMS2VAL,1
          MK2_R = MS2VAL(IMK2_R)
          IACTP = IST_FOR_DT(1,IMK2_R)
          IBCTP = IST_FOR_DT(2,IMK2_R)
*
*  Offset for IBLOCK(sym,type) information
          if (IMK2_R.eq.1) then
            IBOFF = 1
          else
            IBOFF = IBOFF + NBLK_MS2(IMK2_R-1)
          end if
*
          call rewino(LUSCR)
          NBL_C = NBLK_MS2(IMK2_R)
          CALL COPNBLKD(LUC,LUSCR,C,NBL_C,0,LBLK)
          CALL ITODS(-1,1,-1,LUSCR)
          if (NTEST.ge.50) then
            WRITE(6,*)
     &      ' ********************************************* '
            WRITE(6,*)
     &      ' *                                           * '
            WRITE(6,*)
     &      ' * Start of contribution from given MK_2     * '
            WRITE(6,*)
     &      ' *                                           * '
            WRITE(6,*)
     &      ' ********************************************* '
            WRITE(6,*)
            WRITE(6,'(4X,A,2I3)')
     &      ' MK2/IRI vector of C :  MK2_R, IRIR =', MK2_R,IRIR
            WRITE(6,*)
          end if
          write(6,*)
          write(6,*) '´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´ '
          write(6,'(A,X,I3,X,I2)')
     &              ' Contrib. Kramers proj., real/im  ',MK2_R,IRIR
          write(6,*) '´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´´ '
          NAEL = NELEC(IACTP)
          NBEL = NELEC(IBCTP)
*
          NOCTPA = NOCTYP(IACTP)
          NOCTPB = NOCTYP(IBCTP)
*
          IOCTPA = IBSPGPFTP(IACTP)
          IOCTPB = IBSPGPFTP(IBCTP)
*
          IASPGPTP = IACTP
          IBSPGPTP = IBCTP
*
          IF(NTEST.GE.10) THEN
            WRITE(6,*) ' ================'
            WRITE(6,*) ' GASANA speaking '
            WRITE(6,*) ' ================'
            WRITE(6,'(A26)') ' IACTP   IBCTP   NAEL   NBEL'
            WRITE(6,'(3X,4I6)')  IACTP,IBCTP,NAEL,NBEL
            WRITE(6,'(A14,2I6)') ' NOCTPA NOCTPB ', NOCTPA,NOCTPB
          END IF
*
**. Info on block structure of space
*
*
CTF       IF( ICISTR .NE. 1 ) CALL REWINO(LUC)
*. Number of terms to be printed
C?        WRITE(6,*) ' GASANA : IPRNCIV = ', IPRNCIV
          IF(IPRNCIV.EQ.0) THEN
            THRES = 0.02
            MAXTRM = 200
          ELSE
            THRES = 0.0D0
*. Well atmost 100000 coefs - to save disk ..
            MAXTRM = 100000
          END IF
*
          IUSLAB = 0
*. Number of occupation classes
          IWAY = 1
          NEL = NAEL + NBEL
          CALL OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,
     &               IGSOCC(1,1),IGSOCC(1,2))
*. and then the occupation classes
          CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')
          IWAY = 2
          CALL OCCLS(IWAY,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
     &               IGSOCC(1,1),IGSOCC(1,2))
*
          CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
          CALL MEMMAN(KLBSTR,MXNSTR*NBEL,'ADDL  ',1,'KLBSTR')
*
          LENGTH = NOCCLS*10
          CALL MEMMAN(KNCPMT,LENGTH,'ADDL  ',1,'KNCPMT')
          CALL MEMMAN(KWCPMT,LENGTH,'ADDL  ',2,'KWCPMT')
*
*. Occupation of strings of given sym and supergroup
          CALL GASANAS(C,LUSCR,NBL_C,WORK(KNSTSO(IACTP)),
     &                WORK(KNSTSO2(IBCTP)),
     &                NOCTPA,NOCTPB,MXPNGAS,IOCTPA,IOCTPB,
     &                NBLOCK,IBLOCK(1,IBOFF),
     &                THRES,MAXTRM,NAEL,NBEL,
     &                WORK(KLASTR),WORK(KLBSTR),
     &                IBLTP,NSMST,IUSLAB,
     &                IDUMMY,WORK(KNCPMT),WORK(KWCPMT),
     &                NELFSPGP,
     &                NOCCLS,NGAS,WORK(KLOCCLS),
     &                IASPGPTP,IBSPGPTP,
     &                ICISTR,NTOOB,IPRNCIV,IPRWFC)
*
        end do
      end do
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM',IDUM,'GASANA')
      CALL QEXIT('ANACI')
*
C?    Call Abend2( ' Jeppe forced me to stop after INTANA ' )
      RETURN
      END
