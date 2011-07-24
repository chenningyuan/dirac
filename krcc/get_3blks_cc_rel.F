***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_3BLKS_CC_REL(LBLOCK,LSCR12,
     &                            WORK,KFREE,LFREE)
*
* Allocate the three blocks VEC1, VEC2, C2 used in sigma, densi etc
*
* One of the blocks (KC2) is dimensioned so two partitioned blocks of
* length LSCR2 can be stored
*
* Jeppe Olsen, Jan 1998
*               Relativistic version, Feb. 98
*
* IFLAG = 0 : Determine values only, do not allocate
* IFLAG = 1 : Determine values and allocate
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*
#include "mxpdim.inc"
#include "cicisp.inc"
!#include "wrkspc.inc"
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
#include "cands.inc"
*
      DIMENSION WORK(*)
*
      IDUM = 1
C     call memmar(KDUM,IDUM,'MARK  ',IDUM,'GET_3B')
*
      NTESTL = 1000
      NTEST = MAX(NTESTL,IPRDIA)
*
      ISPC = ISSPC
      ISM  = ISSM
*
      LSCR2 = 0
*
      DO IMS2 = 1, NMS2VAL
*
        if (NTEST.ge.100) write(6,*) 'Doing MS2 = ',MS2VAL(IMS2)
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
        CALL MEMGET('REAL',KLCIOIO,NOCTPA*NOCTPB,WORK,KFREE,LFREE)
        CALL MEMGET('REAL',KLCBLTP,NSMST,WORK,KFREE,LFREE)
*
        if (NTEST.ge.500) then
          write(6,*) 'KLCIOIO at ',KLCIOIO
          write(6,*) 'Length is ',NOCTPA*NOCTPB
          write(6,*) 'KLCBLTP at ',KLCBLTP
          write(6,*) 'Length is ',NSMST
        end if
*
        print*,'ISPC,IATP,IBTP',ISPC,IATP,IBTP
        CALL IAIBCM_KRCC(ISPC,IATP,IBTP,WORK(KLCIOIO))
        CALL ZBLTP_KRCC(ISMOST(1,ISSM),NSMST,IDC,WORK(KLCBLTP))
*. Largest block of strings in zero order space
        MXSTBL0 = MXNSTR
*. This must be generalized !
*. alpha and beta strings with an electron removed
        IATPM1 = ITPM1(1,IMS2)
        IBTPM1 = ITPM1(2,IMS2)
*. alpha and beta strings with two electrons removed
        IATPM2 = ITPM2(1,IMS2)
        IBTPM2 = ITPM2(2,IMS2)
*. Largest number of strings of given symmetry and type
        MAXA = MXNSTR
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
          MAXB1 = IMNMX_REL(WORK(KNSTSO2(IBTPM1)),NSMST*NOCTYP(IBTPM1),
     &                      2)
          MAXB = MAX(MAXB,MAXB1)
        END IF
        IF(NBEL.GE.2) THEN
          MAXB1 = IMNMX_REL(WORK(KNSTSO2(IBTPM2)),NSMST*NOCTYP(IBTPM2),
     &                      2)
          MAXB = MAX(MAXB,MAXB1)
        END IF
        MXSTBL = MAX(MAXA,MAXB)
        IF(IPRCIX.GE.2 ) WRITE(6,*)
     &  ' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
        MAXK = MIN(MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
        CALL MXRESC_REL(WORK(KLCIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &                  NSMST,NSTFSMSPGP,NSTFSMSPGP2,MXPNSMST,
     &                  NSMOB,MXPNGAS,NGAS,NOBPTS,NOBPTS2,IPRCIX,MAXK,
     &                  NELFSPGP,
     &                  MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
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
C     call memmar(KDUM,IDUM,'FLUSM ',IDUM,'GET_3B')
      CALL MEMREL('3BLKS',WORK,KLCIOIO,KLCIOIO,KFREE,LFREE)
*
      RETURN
      END
