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

#ifdef MOD_AUTODIFF
c     Return a derivative of the taylor series stored in coeff,
c     with exponents as the derivative order for each variable.
      function taylor_deriv(coeff, nvar, exponents) result(d)

      implicit none

      integer :: i, nvar, exponents(nvar), taylor_index
      real(8) :: d, coeff(*), fact

      common /factrl/ fact(30)

      d = coeff(taylor_index(nvar, exponents))

      do i = 1, nvar
        if (exponents(i) > 0) then
          d = d*fact(exponents(i))
        end if
      end do


      end function
#endif /* ifdef MOD_AUTODIFF */
