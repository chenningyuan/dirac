      subroutine findhigh(IARRAY,LENGTH,MVALUE)
*
*  Find highest integer on IARRAY
*  Positive integers assumed
*
      implicit real*8(A-H,O-Z)
*
      dimension IARRAY(LENGTH)
*
      MVALUE = 0
      do I=1,LENGTH,1
         if (IARRAY(I).gt.MVALUE) MVALUE = IARRAY(I)
      end do
      return
      end
