*     calculate integer A over integer B
*
      subroutine over(AOVERB,A,B)
      implicit none
      real*8 AOVERB
      integer A,B,C,FACULA,FACULB,FACULC
      if (B.gt.A) then
         write(*,*) 'CRAP! Upper < Lower!'
         Call Abend()
      end if
      C = A-B
*
*  Calculate faculties
*
      call faculty(FACULA,A)
      call faculty(FACULB,B)
      call faculty(FACULC,C)
*
      AOVERB = FACULA / (FACULB*FACULC)
      return
      end
