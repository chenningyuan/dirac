*  Routine connects LUCIAREL info with MOLECULE SWEDEN info
*  with respect to basis functions and related things.
*  (at least is planned to do so ...)

      subroutine gtmoswe(XINT_L,NTOTINT)
*
      implicit real*8 (A-H,O-Z)
*
#include "mxpdim.inc"
#include "cgas.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "wrkspc.inc"
*. Input
      dimension NORBMOL(MXNSPIR),NBASMOL(MXNSPIR)
      dimension NFRO(MXNSPIR)
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*. Output
      DIMENSION XINT_L(*)
*
      IDUMMY = 1
      CALL MEMMAN(IDUMMY,IDUMMY,'MARK  ', IDUMMY,'GTMOSW')
      CALL MEMMAN(KLSWINT,NINT2,'ADDL  ',2,'SWEINT')
      write(6,*) 'SWEINT length NINT2  ',NINT2
*
*  Call integral fetching routine ...
*
      do ISYM = 1,MXNSPIR,1
         NORBCT = 0
         NBASCT = 0
         do IGAS = 1,NGAS,1
            NORBCT = NORBCT + NGSSH(ISYM,IGAS)
            NBASCT = NBASCT + NGSSH(ISYM,IGAS)
         end do
         NORBMOL(ISYM) = NORBCT
         NBASMOL(ISYM) = NBASCT
         NFRO(ISYM) = 0
      end do
*
*  The following info is superfluous, as it is provided by LUCIAREL's
*  input. Still, we retain the possibility of using a frozen core in
*  rdmole2 by keeping the arrays NBASMOL, NORBMOL, and NFRO for
*  future purposes.
*
C     LUONEM =
C     REWIND LUONEM
C     CALL SEARCH(LUONEM,'INDATACI')
C     READ (LUONEM) NIRR_PN,(NBASMOL(I),I=1,NIRR_PN),
C    &              (NORBMOL(I),I=1,NIRR_PN),
C    &              PNUCCO,(NFRO(I),I=1,NIRR_PN)
C     IF(IPRINT.GT.1)
C    &   WRITE(6,611) NIRR_PN,(NBASMOL(I),I=1,NIRR_PN),
C    &                (NORBMOL(I),I=1,NIRR_PN),
C    &                (NFRO(I),I=1,NIRR_PN),PNUCCO
611   FORMAT(' NIRR_PN :',I5/' NBASMOL(1:NIRR_PN),NORBMOL(1:NIRR_PN),
     &         NFRO(1:NIRR_PN) '/24I5/' MODIFIED CORE-ENERGY',D13.6)
*
      call rdmole2(WORK(KLSWINT),NORBMOL,NBASMOL,NFRO,NTOTINT,
     &             XINT_L)
CTF
*
      CALL MEMMAN(IDUMMY,IDUMMY,'FLUSM ', IDUMMY,'GTMOSW')
      return
      end
