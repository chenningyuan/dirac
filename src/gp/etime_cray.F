!      integer function lnblnk(string)
!         character string*(*)
!         lnblnk = len_trim(string)
!      end

!      subroutine etime(tarray,result)
!         real, intent(inout), dimension(2) :: tarray
!         real, intent(inout) :: result
! 
!         call myetime(tarray,result)
!      end

      function etime(tarray)
         real, intent(inout), dimension(2) :: tarray
         real :: etime
               
         call myetime(tarray,result)
      end 

      subroutine myetime(tarray,result)
      real,intent(out), dimension(2) :: tarray
      real,intent(out) :: result
      INTEGER(8),save :: count_prev = 0
      INTEGER(8) :: count, count_rate, count_max

      if ( count_prev == 0  ) then
       CALL SYSTEM_CLOCK(count, count_rate, count_max)
      end if

      CALL SYSTEM_CLOCK(count, count_rate, count_max)

      tarray(1) = 0.99*real(count - count_prev) / count_rate
      tarray(2) = 0.01*real(count - count_prev) / count_rate
      result = real(count - count_prev) / count_rate
      count_prev = count
      end subroutine

