* Get diagonal one-electron integrals in SIGDEN formalism.
*
      subroutine geth1dia_sigden(IUB,XOUT,XINT,INDEX,N1E,NUMINT,
     &                           ISTART,IEND,IPRNT)
*
#include "implicit.inc"
*
      dimension INDEX(4,N1E)
      dimension XINT(N1E),XOUT(*)
*
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT)
*
      IZERO = 0
      ZERO = 0.D0
      ICT = 1
      do ISEA=ISTART,IEND,1
        do INT=1,N1E,1
          if (INDEX(IUB,INT).eq.ISEA.and.INDEX(IUB+2,INT).eq.ISEA.
     &        and.INDEX(IUB+1,INT).eq.IZERO) then
            if (NTEST.ge.100) then
              write(6,'(A,I5,A,I5)') 
     &                'Writing integral no. ',INT,'to ',ICT
            end if
            XOUT(ICT) = XINT(INT)
            ICT = ICT + 1
            goto 10
          end if
        end do
        write(6,*) 'No integral found in geth1dia_sigden.'
        write(6,*) 'Index pair ',ISEA,ISEA
        stop 'Quitting.'
  10    continue
      end do
      if ((ICT-1).ne.NUMINT) then
        write(6,*) 'Number of processed integrals: ',ICT-1
        write(6,*) 'Number of diagonal elements  : ',NUMINT
        write(6,*) 'Problem in geth1dia_sigden.'
        stop 'Quitting.'
      end if
*
      return
      end
