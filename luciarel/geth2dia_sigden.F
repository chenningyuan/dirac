* Get diagonal two-electron integrals in SIGDEN formalism.
*   Timo Fleig, August 2001
* Modified to account for integrals over true spinors, i.e.
* explicit reference to unbarred (+) and barred (-) indices (IKFAC).
* Revised routine.
*   Timo Fleig, February 2002
*
      subroutine geth2dia_sigden(XOUT,XINT,INDEX,ISEA,JSEA,TYPE,
     &                           IOFF2,N12,I2DIATP,IKFAC,IKPERM,IPRNT)
*
#include "implicit.inc"
*
      character*4 TYPE
*
      dimension INDEX(4,N12),IAC(2),JAC(2),IOFF2(*),IKPERM(4)
      dimension XINT(N12)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      IERROR = 1
      N2B = IOFF2(I2DIATP)
      do INT = N2B,N12,1
        if (TYPE.eq.'COUL') then
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(2),INT)
          K = INDEX(IKPERM(3),INT) * IKFAC
          L = INDEX(IKPERM(4),INT) * IKFAC
        else if (TYPE.eq.'EXC1') then
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(4),INT) * IKFAC
          K = INDEX(IKPERM(2),INT)
          L = INDEX(IKPERM(3),INT) * IKFAC
        else if (TYPE.eq.'EXC2') then
          I = INDEX(IKPERM(1),INT)
          J = INDEX(IKPERM(4),INT)
          K = INDEX(IKPERM(2),INT) * IKFAC
          L = INDEX(IKPERM(3),INT) * IKFAC
        end if
*
        if (I.eq.ISEA.and.
     &      J.eq.ISEA.and.
     &      K.eq.JSEA.and.
     &      L.eq.JSEA) then
          if (NTEST.ge.100) then
            write(6,'(A,X,A,I5,3X,4I3,1F18.10)') 
     &            'Diagonal int type ',TYPE,
     &             INT,(INDEX(II,INT),II=1,4),XINT(INT)
            write(6,*)
          end if
          XOUT = XINT(INT)
          IERROR = 0
* and exit.
          goto 10
        end if
      end do
      if (IERROR.eq.1) then
        write(6,*) 'No integral found in geth2dia_sigden.'
        write(6,*) 'N2B,N12 = ',N2B,N12
        write(6,*) 'Index pair ',ISEA,JSEA
        write(6,*) 'Indices ',I,J,K,L
        stop 'Quitting.'
      end if
*
  10  continue
      return
      end
