      SUBROUTINE GASDIAT_REL(DIAG,LUDIA,ECORE,ICISTR,IHAM12,
     &           IBLTP,NBLOCK,IBLKFO,KPICI1,INTSPACE,
     &           IDENSI,NIRREP,ISPINFREE,ENVIRO,SIGDEN_ROUTE,IRECO)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
* For DBG-RELCI, March21 98
* Modifed fetch of 1 electron integrals, May 98
*
      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
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
#include "symm.inc"
*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('CIDIA')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*. Largest Number of electrons around
      MAXNEL = 0
      DO IMK2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMK2)
        IBTP = IST_FOR_DT(2,IMK2)
        MAXNEL = MAX(MAXNEL,NELEC(IATP),NELEC(IBTP))
      END DO
*. Local memory
      IDUM = 0
      call memmar(IDUM,  IDUM,    'MARK  ',IDUM,'GASDIA')
*
CTF  Total number of orbitals and active orbitals (spinors, actually)
CTF  from BOTH sets, barred and unbarred !
      NTOOBTOT = NTOOB + NTOOB2
      NACOBTOT = NACOB + NACOB2
*
*
      call memmar(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      call memmar(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLK2  ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLSCR2,2*NTOOBTOT**2,'ADDL  ',2,'KLSC2 ')
      call memmar(KLXA  ,NACOBTOT,   'ADDL  ',2,'KLXA  ')
      call memmar(KLXB  ,NACOBTOT,   'ADDL  ',2,'KLXB  ')
      call memmar(KLSCR ,2*NACOBTOT, 'ADDL  ',2,'KLSCR ')
C     CALL MEMMAR(KLH1D ,NACOBTOT, 'ADDL  ',2,'KLH1D ')
      call memmar(KLH1DA,NTOOB, 'ADDL  ',2,'KLH1DA')
      call memmar(KLH1DB,NTOOB, 'ADDL  ',2,'KLH1DB')
*. Space for blocks of strings
      call memmar(KLASTR,MXNSTR*MAXNEL,'ADDL  ',1,'KLASTR')
      call memmar(KLBSTR,MXNSTR*MAXNEL,'ADDL  ',1,'KLBSTR')
      KLSVST = 1
      call memmar(KLRJKA,MXNSTR, 'ADDL  ',2,'KLRJKA')
*
      if (IRECO.eq.1) then
        if (IHAM12.eq.1.or.IHAM12.eq.2) then
          call geth1dia4(WORK(KLH1DA),WORK(KLH1DB),WORK(KPICI1),
     &                   INTSPACE,ENVIRO,SIGDEN_ROUTE)
        end if
*. Obtain 2e integrals needed for diagonal
        if (IHAM12.EQ.2.or.IHAM12.eq.3) then
          call get_dia_2e(WORK(KLJ),WORK(KLK),WORK(KLK2),IDENSI,
     &                    ISPINFREE,
     &                    ISPSF_MO,ENVIRO,SIGDEN_ROUTE,NIRR_PN,IPRDIA)
        end if
      end if
*
      DO IMK2 = 1, NMS2VAL
        if (NTEST.ge.10) then
          write(6,*) '=============================='
          write(6,*) 'Treating MK space in diagonal '
          write(6,'(A15,I4)') '         MK2 = ',MS2VAL(IMK2)
          write(6,*) '=============================='
        end if
        IF(IMK2.EQ.1) THEN
          IBOFF = 1
        ELSE
          IBOFF = IBOFF + NBLK_MS2(IMK2-1)
        END IF
        NNBLK = NBLK_MS2(IMK2)
*
        IATP = IST_FOR_DT(1,IMK2)
        IBTP = IST_FOR_DT(2,IMK2)
*
        NAEL = NELEC(IATP)
        NBEL = NELEC(IBTP)
*
        NOCTPA = NOCTYP(IATP)
        NOCTPB = NOCTYP(IBTP)
*
        IASPGPTP = IATP
        IBSPGPTP = IBTP
*
*. Offsets for alpha and beta supergroups
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
        IF(NTEST.GE.10) THEN
          WRITE(6,*) ' ================'
          WRITE(6,*) ' GASDIA speaking '
          WRITE(6,*) ' ================'
          WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
          write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
          write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        END IF
*
** Info on block structure of internal state
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
*
        CALL GASDIAS_REL(IMK2,NMS2VAL,MS2VAL,NAEL,WORK(KLASTR),NBEL,
     &         WORK(KLBSTR),NTOOB,DIAG,NSMST,NIRREP,
     &         WORK(KLH1DA),WORK(KLH1DB),
     &         WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &         WORK(KLK),WORK(KLK2),
     &         WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &         LUDIA,ECORE,PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &         WORK(KLRJKA),IHAM12,IBLTP,NNBLK,IBLKFO(1,IBOFF),
     &         IASPGPTP,IBSPGPTP,IREOTS,NGAS,NIRR_PN,MXTSOB,ISPINFREE)
*
      END DO
*     ^End of loop over MK2 spaces
c     stop 'AFTER CONSTRUCTION OF DIAGONAL.'
*.Flush local memory
      call memmar(IDUM,  IDUM,    'FLUSM ',IDUM,'GASDIA')
      CALL QEXIT('CIDIA')
*
      RETURN
      END
