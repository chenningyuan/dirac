!
!     Second-quantized dble(1)-body and two-body matrix elements of total spin operator S^2
!
!          S^2 = S2_1 + S2_2
!
!     dble(1)-body part of S^2 = S2_1
!
!          dble(1)-body part of Sz    +
!          dble(1)-body part of Sz*Sz +
!          dble(1)-body part of S-*S+
!
!     S2_1 = + sum(i,j) [ <i|s2_1|j> i+ j ]
!
!     <i|s2_1|j> = + <i|sz|j>
!                  + sum (k) [ <i|s-|k>*<k|s+|j> ]
!                  + sum (k) [ <i|sz|k>*<k|sz|j> ]
!
!     two-body part of S^2 = S2_2
!
!          two-body part of Sz*Sz +
!          two-body part of S-*S+
!
!     S2_2 = + sum(i>j,k>l) [ <i,j|s2_2|k,l> i+ j+ l k ]
!
!     where, 
!
!     <i,j|s2_2|k,l> is anti-symmetrized element defined as follows.
!
!     <i,j|s2_2|k,l> = + <i|sz|k>*<j|sz|l> - <i|sz|l>*<j|sz|k>
!                      - <j|sz|k>*<i|sz|l> + <j|sz|l>*<i|sz|k>
!                      + <i|s-|k>*<j|s+|l> - <i|s-|l>*<j|s+|k>
!                      - <j|s-|k>*<i|s+|l> + <j|s-|l>*<i|s+|k>
!
!
      module S2MatEle

      use Parameters
      use IOModule
      use GMSMod
      implicit none
      private

      logical :: Init = .false.
      double precision, allocatable :: s2_1(:,:), s2_2(:,:)
      public :: s2_1, s2_2
      public :: two_s_plus1
      public :: InitS2MatEle
      public :: EndS2MatEle

      contains

      subroutine InitS2MatEle

      implicit none

      if (Init) then
         write(UOUT,*) 'attempt to reinitialize init_s2matele'
         call flushout
         call quit('InitS2MatEle: attempt to reinitialize init_s2matele')
      endif

      Init  = .true.

      if (get_iprint() >= 7) then
        write(UOUT,'(/,2x,a,i4)') 'InitS2MatEle:   NVMO=',NVMO
        write(UOUT,'(2x,a,i8)')   'InitS2MatEle: NVMOSQ=',NVMOSQ
      endif

      allocate(s2_1(NVMO,NVMO))
      allocate(s2_2(NVMOSQ,NVMOSQ))

      call CalcS2MatEle

      return

      end subroutine InitS2MatEle

      subroutine EndS2MatEle

      implicit none

      if (.not.Init) then
         write(UOUT,*) 'error in end_s2matele'
         call flushout
         call quit('error in end_s2matele')
      endif

      Init = .false.
      deallocate(s2_1)
      deallocate(s2_2)

      return

      end subroutine EndS2MatEle

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine CalcS2MatEle

      implicit none
      integer i, j, k, l, ij, kl

!     dble(1)-body part of S^2 =
!
!          dble(1)-body part of Sz    +
!          dble(1)-body part of Sz*Sz +
!          dble(1)-body part of S-*S+
!
!     +sum(i,j) [ <i|s2_1|j> i+ j ]
!
!     <i|s2_1|j> = + <i|sz|j>
!                  + sum (k) [ <i|s-|k>*<k|s+|j> ]
!                  + sum (k) [ <i|sz|k>*<k|sz|j> ]
!

      if (get_iprint() >= 7) then
        write(UOUT,'(/,2x,a,i4)') 'in CalcS2MatEle: NVMO=',NVMO
      endif

      do j = 1, NVMO
      do i = 1, NVMO

         s2_1(i,j) = dble(0)
         s2_1(i,j) = s2_1(i,j) + opsz(i,j)

         do k = 1, NVMO

            s2_1(i,j) = s2_1(i,j) + opsminus(i,k)*opsplus(k,j)
            s2_1(i,j) = s2_1(i,j) + opsz(i,k)*opsz(k,j)

         enddo

         if (abs(s2_1(i,j)).le.1.0d-10) s2_1(i,j) = dble(0)

      enddo
      enddo

!     two-body part of S^2 =
!
!          two-body part of Sz*Sz +
!          two-body part of S-*S+
!
!     +sum(i>j,k>l) [ <i,j|s2_2|k,l> i+ j+ l k ]
!
!     <i,j|s2_2|k,l> = + <i|sz|k>*<j|sz|l> - <i|sz|l>*<j|sz|k>
!                      - <j|sz|k>*<i|sz|l> + <j|sz|l>*<i|sz|k>
!                      + <i|s-|k>*<j|s+|l> - <i|s-|l>*<j|s+|k>
!                      - <j|s-|k>*<i|s+|l> + <j|s-|l>*<i|s+|k>

      do i = 2, NVMO
      do j = 1, (i-1)

         ij = MAXMIN(i,j)

         do k = 2, NVMO
         do l = 1, (k-1)

            kl = MAXMIN(k,l)

            s2_2(kl,ij) = dble(0)
            s2_2(kl,ij) = s2_2(kl,ij) + opsz(i,k)*opsz(j,l)
            s2_2(kl,ij) = s2_2(kl,ij) - opsz(i,l)*opsz(j,k)
            s2_2(kl,ij) = s2_2(kl,ij) - opsz(j,k)*opsz(i,l)
            s2_2(kl,ij) = s2_2(kl,ij) + opsz(j,l)*opsz(i,k)
            s2_2(kl,ij) = s2_2(kl,ij) + opsminus(i,k)*opsplus(j,l)
            s2_2(kl,ij) = s2_2(kl,ij) - opsminus(i,l)*opsplus(j,k)
            s2_2(kl,ij) = s2_2(kl,ij) - opsminus(j,k)*opsplus(i,l)
            s2_2(kl,ij) = s2_2(kl,ij) + opsminus(j,l)*opsplus(i,k)

            if (abs(s2_2(kl,ij)).le.1.0d-10) s2_2(kl,ij) = dble(0)

         enddo
         enddo

      enddo
      enddo

      return
      
      end subroutine CalcS2MatEle

      function opsplus(i,j)
!      ...  returns <i|s+|j>

     ! integer i, j
      integer, intent(in) :: i, j
      double precision opsplus

      opsplus = dble(0)
      if (SpinOfActvOrb(i).eq.BEETA_SPIN) return
      if (SpinOfActvOrb(j).eq.ALPHA_SPIN) return
      if ((j-i).ne.1) return
      opsplus = dble(1)

      return

      end function opsplus

      function opsminus(i,j)
!      ... returns <i|s-|j>

     ! integer i, j
      integer, intent(in) :: i, j
      double precision opsminus

      opsminus = dble(0)

      if (SpinOfActvOrb(i).eq.ALPHA_SPIN) return
      if (SpinOfActvOrb(j).eq.BEETA_SPIN) return
      if ((i-j).ne.1) return
      opsminus = dble(1)

      return

      end function opsminus

      function opsz(i,j) result(opsz_out)
!      ... returns <i|sz|j>
!     MI: modified to this form when catching some bugs...
!     ... to fix it it's sufficient switch off pgdbg 5.2.4 -pg flag ...
      integer, intent(in) :: i, j
      double precision opsz_out

      opsz_out = dble(0)

      if (i.ne.j) return
      if (SpinOfActvOrb(i).eq.ALPHA_SPIN) opsz_out = +half
      if (SpinOfActvOrb(i).eq.BEETA_SPIN) opsz_out = -half

      return

      end function opsz

      integer function two_s_plus1(ssplus1)

      implicit none
      double precision ssplus1
      double precision d
       
      d = one + dble(4)*ssplus1
      d = sqrt(d)
      two_s_plus1 = nint(d)

      return

      end function two_s_plus1

#include "maxmin.h"

      end module S2MatEle
