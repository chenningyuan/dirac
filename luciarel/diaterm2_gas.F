      SUBROUTINE DIATERM2_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
     &                       JPERT,J12,JDC)
* = DIATERM_GAS, just J12 added !
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
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
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*
      INTEGER IBLOCK(8,*)
*
      DIMENSION VEC(*)
*
      CALL QENTER('DIATR')
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIATRM')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
C     IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
C       I12 = 2
C       IPERTOP = 0
C     ELSE
*. Use perturbation operator
C       IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
C         I12 = 1
C         IPERTOP = 1
C       ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
C         I12 = 2
C         IPERTOP = 0
C       END IF
C     END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   DIATERM2_GAS speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,J12,IPERTOP',JPERT,IPART,J12,IPERTOP
      END IF
*. A bit of scracth
      CALL MEMMAN(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      CALL MEMMAN(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      CALL MEMMAN(KLSCR2,2*NTOOB**2,'ADDL  ',2,'KLSC2 ')
      CALL MEMMAN(KLXA  ,NACOB,   'ADDL  ',2,'KLXA  ')
      CALL MEMMAN(KLXB  ,NACOB,   'ADDL  ',2,'KLXB  ')
      CALL MEMMAN(KLSCR ,2*NACOB, 'ADDL  ',2,'KLSCR ')
      CALL MEMMAN(KLH1D ,NACOB,   'ADDL  ',2,'KLH1D ')
*. Space for blocks of strings
C     WRITE(6,*) ' MXNSTR in DIATERM', MXNSTR
      CALL MEMMAN(KLASTR,MXNSTR*NAEL,'ADDL  ',1,'KLASTR')
      CALL MEMMAN(KLBSTR,MXNSTR*NAEL,'ADDL  ',1,'KLBSTR')
      MAXA = IMNMX(WORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL MEMMAN(KLRJKA,MAXA,'ADDL  ',2,'KLRJKA')
*. Diagonal of one-body integrals and coulomb and exchange integrals
*. Integrals assumed in place so :
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
C     CALL GT1DIA(WORK(KLH1D))
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      IF(J12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*. Core energy not included
      ECOREP = 0.0D0
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS)
*
      CALL DIATERMS_GAS(NAEL,WORK(KLASTR),NBEL,WORK(KLBSTR),
     &             NACOB,VEC,NSMST,WORK(KLH1D),
     &             JDC,WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &             ECOREP,0,0,
     &             IPRDIA,NTOOB,
     &             WORK(KLRJKA),
     &             J12,IBLOCK(1,IOFF),NBLOCK,ITASK,FACTOR,0,0)
C    &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
*.Flush local memory
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ',IDUMMY,'DIATRM')
      CALL QEXIT('DIATR')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTS(VEC,IBLOCK(1,IOFF),NBLOCK,
     &              NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC)
      END IF
*
      RETURN
      END
