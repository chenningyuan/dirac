      SUBROUTINE GET_3BLKS(KVEC1,KVEC2,KC2,LSCR2)
*
* Allocate the three blocks VEC1, VEC2, C2 used in sigma, densi etc
*
* One of the blocks (KC2) is dimensioned so two partitioned blocks of
* length LSCR2 can be stored
*
* Jeppe Olsen, Jan 1998
*               Relativistic version, Feb. 98
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "strbas.inc"
#include "glbbas.inc"
#include "cprnt.inc"
#include "oper.inc"
#include "gasstr.inc"
*
#include "cgas.inc"
#include "lucinp.inc"
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*
      IDUM = 1
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_3B')
*
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRDIA)
*
      ISPC = ISSPC
      ISM  = ISSM
*
      LSCR2 = 0
*
      DO IMS2 = 1, NMS2VAL
*
        if (NTEST.ge.500) write(6,*) 'Doing IMS2 = ',IMS2
        IATP = IST_FOR_DT(1,IMS2)
        IBTP = IST_FOR_DT(2,IMS2)
*
        NAEL = NELEC(IATP)
        NBEL = NELEC(IBTP)
*
        NOCTPA = NOCTYP(IATP)
        NOCTPB = NOCTYP(IBTP)
*
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
*
        CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
        CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*
        if (NTEST.ge.500) then
          write(6,*) 'KLCIOIO at ',KLCIOIO
          write(6,*) 'Length is ',NOCTPA*NOCTPB
          write(6,*) 'KLCBLTP at ',KLCBLTP
          write(6,*) 'Length is ',NSMST
        end if
*
        CALL IAIBCM(ISPC,IATP,IBTP,WORK(KLCIOIO))
CTI        KSVST = 1
CTI        CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTP),WORK(KSVST))
        CALL ZBLTP(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTP))
*. Largest block of strings in zero order space
        MXSTBL0 = MXNSTR
*. This must be generalized !
*. alpha and beta strings with an electron removed
        IATPM1 = 3
        IBTPM1 = 4
*. alpha and beta strings with two electrons removed
        IATPM2 = 5
        IBTPM2 = 6
*. Largest number of strings of given symmetry and type
        MAXA = MXNSTR
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
          MAXB1 = IMNMX(WORK(KNSTSO2(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
          MAXB = MAX(MAXB,MAXB1)
        END IF
        IF(NBEL.GE.2) THEN
          MAXB1 = IMNMX(WORK(KNSTSO2(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
          MAXB = MAX(MAXB,MAXB1)
        END IF
        MXSTBL = MAX(MAXA,MAXB)
        IF(IPRCIX.GE.2 ) WRITE(6,*)
     &  ' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
        MAXI = MIN( MXINKA,MXSTBL)
        MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
        CALL MXRESC(WORK(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                NSMST,NSTFSMSPGP,MXPNSMST,
     &                NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &                NELFSPGP,
     &                MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
        LSCR2 = MAX(LSCR2,MXCJ,MXCIJA,MXCIJB)
        LCJRES = LSCR2
*. Note : LCJRES and MXADKBLK are being save in /CICISP/
        IF(IPRCIX.GE.2) THEN
          WRITE(6,*) 'GET_3BL: MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
          WRITE(6,*) 'GET_3BL: MXADKBLK ', MXADKBLK
          WRITE(6,*) ' Space for resolution matrices ',LSCR2
        END IF
      END DO
*     ^ End loop over MK2 values
*
      IF(ICISTR.EQ.1) THEN
        WRITE(6,*) ' ICISTR = 0 not in action ANYMORE'
        Call Abend2( ' ICISTR = 0 not in action ANYMORE' )
      ELSE IF (ICISTR.EQ.2) THEN
        LBLOCK = MXSB
      ELSE IF (ICISTR.EQ.3) THEN
        LBLOCK = MXSOOB
      END IF
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LSCR12 = MAX(LBLOCK,2*LSCR2)
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_3B')
*
      CALL MEMMAN(KVEC1,LBLOCK,'ADDS  ',2,'VEC1  ')
      CALL MEMMAN(KVEC2,LBLOCK,'ADDS  ',2,'VEC2  ')
      CALL MEMMAN(KC2,LSCR12,'ADDS  ',2,'C2    ')
*
      RETURN
      END
