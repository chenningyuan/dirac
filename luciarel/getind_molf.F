* Determine 1- or 2- particle integral
* Fetch integral.
*
      subroutine getind_molf(YINT,INDEX,IRI,IDIM,INVREO,
     &                       KR,NASHT,DMFC,VMUUUU,IPRNT)
*
#include "implicit.inc"
*
      dimension VMUUUU(*),DMFC(*)
      dimension ICSAVE(4),INDEX(4),INVREO(-IDIM:IDIM),KR(-IDIM:IDIM)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRNT)
*
      if (NTEST.ge.5) then
        write(6,*) 
        write(6,*) ' GETIND_MOLF : INDEX array'
        write(6,*) (INDEX(I),I=1,4,1)
      end if
*
      if (NTEST.ge.20) then
        write(6,*) '1-el. integral array in getind_sigden:'
        call wrtmat(DMFC,IDIM*2,IDIM*2,IDIM*2,IDIM*2)
      end if
*
      NCT = 0
      do IC=1,4,1
        if (INDEX(IC).ne.0) then
          NCT = NCT + 1
          ICSAVE(NCT) = IC
        end if
      end do
*
      if (NCT.eq.2) then
        if (NTEST.ge.10) write(6,*) '1-electron integral.'
        I = KR(INVREO(INDEX(ICSAVE(1))))
        J = KR(INVREO(INDEX(ICSAVE(2))))
* Compute input address
        IIN =  (IRI-1)*NASHT*NASHT +
     &          (J-1) *NASHT +
     &           I
        if (NTEST.ge.5) then
          write(6,*) 'Index pair ',I,J,' address :',IIN
        end if
        YINT = DMFC(IIN)
*
* Sign for imaginary 1-electron integrals:
* (According to h_(IJ) = -h_(I_J_). 
*  Reference are the unbarred spinors.)
        if (IRI.eq.2.and.INDEX(1).eq.0) YINT = - YINT
*
        if (NTEST.ge.2) write(6,'(2I3,2X,1F18.10)') I,J,YINT
*
      else if (NCT.eq.4) then
        I = KR(INVREO(INDEX(1)))
        J = KR(INVREO(INDEX(2)))
        K = KR(INVREO(INDEX(3)))
        L = KR(INVREO(INDEX(4)))
*
* Compute input address
        IIN = (IRI-1)*NASHT*NASHT*NASHT*NASHT +
     &         (L-1) *NASHT*NASHT*NASHT +
     &         (K-1) *NASHT*NASHT +
     &         (J-1) *NASHT +
     &          I
        if (NTEST.ge.5) then
          write(6,*) 'Index quad ',I,J,K,L,' address :',IIN
        end if
        YINT = VMUUUU(IIN)
        if (NTEST.ge.2) write(6,'(4I3,2X,1F18.10)') I,J,K,L,YINT
*
      else
        write(6,*) 'No QFT implemented YET in LUCIAREL.'
        write(6,*) 'you specified ',4-NCT,' indices.'
        call abend2('Quitting.')
      end if
*
      return
      end
