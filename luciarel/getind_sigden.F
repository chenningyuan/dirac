* Slave (yes, slave. We are NOT politically correct in 
* FORTRAN programs!) routine for fetching desired integral
* from DIRAC list; both 1- and 2-electron
*
      subroutine getind_sigden(YINT,INDEX,IRI,ITOTDIM,
     &                         NUMINT,KR,INDK,INDL,INVREO,
     &                         IOFF12,XIN,ONEINT,IPRNT)
*
#include "implicit.inc"
*
      dimension INDEX(4),ICSAVE(4),INVREO(-ITOTDIM:ITOTDIM)
      dimension ONEINT(ITOTDIM*2,ITOTDIM*2)
*
      NTESTL = 00
      NTEST = max(IPRNT,NTESTL)
*
      if (NTEST.ge.20) then
        write(6,*) '1-el. integral array in getind_sigden:'
        call wrtmat(ONEINT,ITOTDIM*2,ITOTDIM*2,ITOTDIM*2,ITOTDIM*2)
      end if
*
* 1- or 2-electron integral?
      NCT = 0
      do IC=1,4,1
        if (INDEX(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do 
      if (NCT.eq.2) then
        if (NTEST.ge.10) then
          write(6,*) '1-electron integral.'
          write(6,'(A,2I4)') 'Index  pair  : ',INDEX(ICSAVE(1)),
     &                                         INDEX(ICSAVE(2))
          write(6,'(A,2I4)') 'reordered to : ',INVREO(INDEX(ICSAVE(1))),
     &                                         INVREO(INDEX(ICSAVE(2)))
        end if
        I = INVREO(INDEX(ICSAVE(1)))
        J = INVREO(INDEX(ICSAVE(2)))
        YINT = ONEINT(I,J)
*
* Sign for imaginary 1-electron integrals:
* (According to h_(IJ) = -h_(I_J_) )
        if (IRI.eq.2.and.INDEX(2).eq.0) YINT = - YINT
* FIXME : Why is it the unbarred integrals have to be multiplied??
*
      else if (NCT.eq.4) then
        if (NTEST.ge.10) then
          write(6,*) '2-electron integral.'
          write(6,'(A,4I4)') 'Index quadr. : ',INDEX(1),
     &                                         INDEX(2),
     &                                         INDEX(3),
     &                                         INDEX(4)
          write(6,'(A,4I4)') 'reordered to : ',INVREO(INDEX(1)),
     &                                         INVREO(INDEX(2)),
     &                                         INVREO(INDEX(3)),
     &                                         INVREO(INDEX(4))
        end if
        I = INVREO(INDEX(1))
        J = INVREO(INDEX(2))
        K = INVREO(INDEX(3))
        L = INVREO(INDEX(4))
*
        call getind_index(YINT,I,J,K,L,NUMINT,KR,INDK,INDL,
     &                    IOFF12,IRI,XIN,IPRNT)
*	
      else
        write(6,*) 'No QFT implemented YET in LUCIAREL.'
        write(6,*) 'you specified ',4-NCT,' indices.'
        write(6,*) 'INDEX array: ',(INDEX(I),I=1,4,1)
        call abend2('Quitting.')
      end if
*
      return
      end
