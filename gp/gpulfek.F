

      function ifind(i,arr,len)
C     Return the index of the first occurance of i in arr, or 0 if not
C     found.
      implicit none
      integer ifind, k,i,arr,len
      dimension arr(len)
      do k=1,len
         if (arr(k).eq.i) then
            ifind = k
            return
         endif
      enddo
      ifind = 0
      end
