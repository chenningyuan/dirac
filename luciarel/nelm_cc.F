* Find number of elements in a given excitation class.
*
      function nelm_cc(IOFF,ITP,NXTP,NT)
*
#include "implicit.inc"
*
      dimension IOFF(NXTP)
*
      if (ITP.lt.NXTP) then
        NELM_CC = IOFF(ITP+1) - IOFF(ITP)
      else if (ITP.eq.NXTP) then
        NELM_CC = NT - IOFF(ITP) + 1
      else
        write(6,*) 'Requested excitation type ',ITP
        write(6,*) 'exceeds total number of types ',NXTP
        stop 'QUITTING in function nelm_cc.'
      end if
*
      return
      end
