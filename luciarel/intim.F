      SUBROUTINE INTIM(IPRTRA)
*
* Interface to external integrals
*
* If NOINT .ne. 0, only pointers are constructed
* Jeppe Olsen, Winter of 1991
*              Last revision : Jan 98, double group integrals added
* Note : In this version only pointers are constructed, the rest is 
*        realized in picasso.
*    Timo Fleig
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "glbbas.inc"
#include "clunit.inc"
#include "lucinp.inc"
#include "csm.inc"
#include "cintfo.inc"
#include "orbinp.inc"
#include "symm.inc"
#include "input.inc"
      COMMON/CECORE/ECORE
*
      NTEST = 0000
*
      CALL LMEMCHK
*
*. : Pointers for symmetry blocks of integrals
*
      if (SIGDEN_ROUTE.eq.'TRADIT') then
        CALL INTPNT_PNG(WORK(KPINT1),WORK(KPINT2))
*
* Pointers for double group integrals
        IF(IDBG.EQ.1) CALL INTPNT_DBG
      end if
*
*. Pointer for orbital indeces for symmetry blocked matrices
      CALL ORBINH1(WORK(KINH1),NTOOBS,NTOOB,NSMOB)
      I_AM_TIMO_WORLD = 0
      if (INTIMP.eq.8) I_AM_TIMO_WORLD = 1
CTF
*  We read integrals in routine DIST_2EINTS coming up.
*  There, the individual environment flags are relevant.
CTF
      I_DO_2EINTS = 1
      IF(I_DO_2EINTS .EQ. 1 ) THEN
*. Skip 1-electron integrals
c       CALL REWINO(LU2INT)
c       DO INT1 = 1, NINT1
c         READ(LU2INT,'(E22.15)') XDUMMY
c       END DO
*. 2-electron integrals
c       READ(LU2INT,'(E22.15)') (WORK(KINT2-1+INT2),INT2=1,NINT2)
      END IF
c     WRITE(6,*) ' Integrals read in '
c     CALL WRTMAT(WORK(KINT2),NINT2,1,NINT2,1)
*
      IF (I_AM_TIMO_WORLD.EQ.0.AND.
     &       INTIMP.EQ.1.AND.NOINT.EQ.0) THEN
*
* ===============
* . MOLCAS format
* ===============
*
        WRITE(6,*) ' Integrals imported from MOLCAS files'
*.Initialize buffers, open
        CALL MKLUNDIO
*. Load one-electron integrals
CTF     CALL GETH0(WORK(KINT1))
*. And two-electron integrals if desired
        IF(INCORE.EQ.1) THEN
          CALL INTIMM(WORK(KINT2),NIRR_PN,SYDI,IPRTRA)
        END IF
*
* Check read integrals
        if (NTEST.ge.1000) then
          do IINT = 0,NINT2-1,1
            write(6,'(A4,I3,X,F18.14)') 'No. ',IINT+1,WORK(KINT2+IINT)
          end do
        end if
*
* Formatted input of symmetry non-vanishing integrals
*
      ELSE IF (INTIMP.EQ.3.AND.NOINT.EQ.0 ) THEN
        WRITE(6,*) ' Integrals imported formatted (E22.15) '
*.1 : One-electron integrals
        REWIND LU2INT
        READ(LU2INT,'(E22.15)') (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
        IF(INCORE.EQ.1.OR.EXTSPC.EQ.0) THEN
          READ(LU2INT,'(E22.15)') (WORK(KINT2-1+INT2),INT2=1,NINT2)
        END IF
*.3 : Core energy
        WRITE(6,*) ' Core energy will be read '
        READ(LU2INT,'(E22.15)') ECORE
*
*. Proceed in Normal MOLCAS way so
        INTIMP = 1
*
* ===============
* . SIRIUS format
* ===============
*
      ELSE IF (INTIMP .EQ. 5.AND.NOINT.EQ.0  ) THEN
        WRITE(6,*) ' Integrals imported from SIRIUS files'
*. Load one-electron integrals
CTF     CALL GETH0S(WORK(KINT1),NOCOB)
*. And two-electron integrals if desired
        IF(INCORE.EQ.1) THEN
          LSCR = 100000
          CALL MEMMAN(KLSCR,LSCR,'ADDL  ',2,'KLSCR ')
CTF       CALL READMO(WORK(KLSCR),LSCR,WORK(KINT2))
        END IF
      else if (I_AM_TIMO_WORLD.eq.1.and.INTIMP.eq.8
     &         .and.NOINT.eq.0) then
*
* ===================
*  MOLECULE format
* ===================
*
         if (ARCHI.eq.'RS6') then
           write(6,*) 'Reading from file 13 in MOLECULE format'
           write(6,*) 'Reading from file 13 in MOLECULE format'
           write(6,*) 'Reading from file 13 in MOLECULE format'
           write(6,*) 'Reading from file 13 in MOLECULE format'
         else if (ARCHI.eq.'SGO') then
           write(6,*) 'Reading from ASCII files'
           write(6,*) 'Reading from ASCII files'
           write(6,*) 'Reading from ASCII files'
           write(6,*) 'Reading from ASCII files'
         end if
CTF
*  Interface to MOLECULE-SWEDEN integrals (Jan Almloef, 1986)
*  Patch from program tradpt written by Christel Marian in 1987.
CTF
         call gtmoswe(WORK(KINT2),NTOTINT)
         write(6,*) 'Total number of 2-el. ints
     &               that have been read:         ',NTOTINT
      END IF
*
CTF
*
C     write(6,*) ' Memory check   INTIM 2 '
*. integrals (ij!kk) -(ik!kj)
C       ZIJPNT(IJPNT,IJSM,IGEJ,NOBSM,IOBSM,LOBSM,ISTOB,NTOTOB)
*... Here I am
CTF      CALL ZIJPNT(WORK(KPNIJ),1,1,NSMOB,IBSO,NTOOBS,IREOST,NTOOB)
C?    write(6,*) ' Home from ZIJPNT '
C    IJKK(XIJKK,IGEJ,IJPNT,NIJ,IJSM,NOBSM,LOBSM,IOBSM,NTOTOB)
CTF      IF(NOINT.EQ.0)
CTF     &CALL IJKK(WORK(KIJKK),1,WORK(KPNIJ),1,NSMOB,NTOOBS,
CTF     &          IBSO,IREOST,NTOOB)
CTF      WRITE(6,*) ' Memcheck after IJKK '
CTF      CALL LMEMCHK
C?    write(6,*) ' home from IJKK '
*. Change one-electron integrals to inactive fock matrix
CTF      IF(NOINT.EQ.0) CALL FI(WORK(KINT1),ECC)
CTF      ECORE = ECORE + ECC
CTF      WRITE(6,*) ' Updated core energy ',ECORE
CTF      CALL LMEMCHK
*
C?    WRITE(6,*) ' IDMPIN ', IDMPIN
CTF      IF (IDMPIN.EQ.1 ) THEN
CTF        WRITE(6,*)
CTF     &   ' Integrals written formatted (E22.15) on unit 90'
CTF        LU90 =
CTF        REWIND LU90
*.1 : One-electron integrals
CTF        WRITE(LU90,'(E22.15)')
CTF     &   (WORK(KINT1-1+INT1),INT1=1,NINT1)
*.2 : Two-electron integrals
CTF        WRITE(LU90,'(E22.15)')
CTF     &   (WORK(KINT2-1+INT2),INT2=1,NINT2)
*.3. Core energy
CTF        WRITE(LU90,'(E22.15)')ECORE
*.4  Rewind to empty buffer
CTF        REWIND LU90

C!    stop ' Jeppe forced my to stop in INTIM '
      RETURN
      END
