* Write density matrix for current eigenvector in 
* SIGDEN ordering.
*
*   October 18, 2001, Timo Fleig
*
*   Revised for operator ordering, Mar 20, 2002
*
      subroutine wrtdens(DENS,INDEX,IHINDX,SIGN,IOP_REO,ISF)
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
*
      dimension DENS(NINT_CC),SIGN(NSPOBEX_TP)
      dimension INDEX(4,NINT_CC),IND_R(4),IHINDX(4,NSPOBEX_TP),
     &          IOP_REO(4,NSPOBEX_TP)
*
      NINT_TOT = N1ELINT + N2ELINT
      if (ISF.eq.1) then
        IRILP = 1
      else
        IRILP = 2
      end if
*
      do IRI = 1,IRILP,1
        IRUN = 1
        IRIOFF = 0
        do JRI = 1,IRI-1,1
          IRIOFF = IRIOFF + NINT_TOT
        end do
        write(6,*)
        write(6,*)
        if (IRI.eq.1) write(6,*) ' Real density elements     : '
        if (IRI.eq.2) write(6,*) ' Imaginary density elements: '
        write(6,*) ' -------------------------------'
        do IXTP = 1,NSPOBEX_TP,1
          SIGN_REO = 1.D0
          write(6,*)
          write(6,'(A,I4)') ' Excitation type ... ',IXTP
          write(6,*)        ' -------------------------'
          NINT_TP = NELM_CC(WORK(KLIBSOBEX),IXTP,NSPOBEX_TP,NINT_TOT)
          do INT = 1,NINT_TP,1
*
            if (INDEX(1,IRUN).ne.0.and.INDEX(2,IRUN).ne.0) then
*  Two-particle density element.
*
*  INDEX array refers to integral indexing.
*  The corresponding second quantization operators then result
*  from the following resorting to CACA ordering:
*  CA AA CB AB for printing:
              IND_R(IOP_REO(1,IXTP)) = INDEX(IHINDX(1,IXTP),IRUN)
              IND_R(IOP_REO(2,IXTP)) = INDEX(IHINDX(2,IXTP),IRUN)
              IND_R(IOP_REO(3,IXTP)) = INDEX(IHINDX(3,IXTP),IRUN)
              IND_R(IOP_REO(4,IXTP)) = INDEX(IHINDX(4,IXTP),IRUN)
              SIGN_REO = -1.D0
            else
*  One-particle density element. No modification.
              IND_R(1) = INDEX(1,IRUN)
              IND_R(2) = INDEX(2,IRUN)
              IND_R(3) = INDEX(3,IRUN)
              IND_R(4) = INDEX(4,IRUN)
              SIGN_REO = 1.D0
            end if
            write(6,'(X,4I3,2X,1F18.10)') 
     &           (IND_R(I),I=1,4,1),DENS(IRUN+IRIOFF)
     &                              *SIGN_REO*SIGN(IXTP)
            IRUN = IRUN + 1
          end do
        end do
      end do
c     stop 'AFTER WRTDENS.'
*
      return
      end
