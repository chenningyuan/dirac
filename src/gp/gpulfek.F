!dirac_copyright_start
!      Copyright (c) by the authors of DIRAC.
!
!      This program is free software; you can redistribute it and/or
!      modify it under the terms of the GNU Lesser General Public
!      License version 2.1 as published by the Free Software Foundation.
!
!      This program is distributed in the hope that it will be useful,
!      but WITHOUT ANY WARRANTY; without even the implied warranty of
!      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!      Lesser General Public License for more details.
!
!      If a copy of the GNU LGPL v2.1 was not distributed with this
!      code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!dirac_copyright_end

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

      LOGICAL FUNCTION ILL_NUMBER(x)
! Returns true if the input number is "ill", NaN,+Inf or -Inf.
! Calls internal C-function "is_ill_number(x)" in gp/gpc.c
      real*8, intent(in) :: x
      external :: is_ill_number
      if (is_ill_number(x) .ne. 0) then
          ILL_NUMBER = .True.
      else
          ILL_NUMBER = .False.
      endif
      return
      end 
