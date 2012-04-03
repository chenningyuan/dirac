      SUBROUTINE GET_3BLKS_GCC_REL(KVEC1,KVEC2,KC2,MXCJ,
     &                             WORK,KFREE,LFREE)
*
* Allocate the three blocks VEC1, VEC2, C2 used in sigma, densi etc
*
* Version for general coupled cluster
*
* Size of Intermediate arrays are simply obtained from allowed max
* on resolution strings
*
* Jeppe Olsen, Summer of 99
*
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cicisp.inc"
#include "orbinp.inc"
#include "clunit.inc"
#include "csm.inc"
#include "cstate.inc"
#include "crun.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "glbbas.inc"
#include "gasstr.inc"
#include "cgas.inc"
#include "lucinp.inc"
*
      DIMENSION WORK(*)
*
      IDUM = 1
*
      NTEST = 000
*
* 1 : Memory required for CC
*
      CALL GET_3BLKS_CC_REL(LVEC1,LC2,WORK,KFREE,LFREE)
      WRITE(6,*) ' Length of vectors for CI  : ', LVEC1, LC2
      print*,'LC2,MXSOOB_AS,MXINKA_CC',LC2,MXSOOB_AS,MXINKA_CC
*
* 2 : Length of vectors used in CC
      MXCJ = MXINKA_CC ** 4
      LBLOCK = MXSOOB_AS
      print*,' BEFORE COMPARISON.'
      print*,' MXINKA_CC,LBLOCK,LVEC1,LC2 ',MXINKA_CC,LBLOCK,LVEC1,LC2
*
      LBLOCK = MAX(LBLOCK,LCSBLK)
      LSCR12 = MAX(LBLOCK,2*MXCJ)  
*
      WRITE(6,*) ' Length of vectors for CC : ', LBLOCK,LSCR12
*
* 3 : Max lengths 
*
      LBLOCK = MAX(LBLOCK,LVEC1)
      LSCR12 = MAX(LSCR12,LC2)
*
      print*,'lblock,lscr12  in GET_3BLKS_GCC_REL: ',LBLOCK,LSCR12
      CALL MEMGET('REAL',KVEC1,IRECOM*LBLOCK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KVEC2,IRECOM*LBLOCK,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KC2,IRECOM*LSCR12,WORK,KFREE,LFREE)
*
      RETURN
      END
