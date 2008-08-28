!     n_choose_r(n,r) returns nCr, the number of ways
!     of choosing r objects from n numbered objects
!     without regard to the order of choice.
!
!     direct calculation using factorials runs into
!     integer overflow problems. so use pascal triangle
!     with recursion formula (n,r) = (n-1,r) + (n-1,r-1)
!     
      function n_choose_r(n,r)

      implicit none
      integer n, r, n_choose_r
      integer n_choose_r_internal
      logical lerr

      lerr = .false.
      if ((n.lt.0).or.(r.lt.0)) lerr = .true.
      if (n.lt.r) lerr = .true.

   ! if (lerr) stop "invalid parameters to n_choose_r"
      if (lerr) then
        write(*,*) 'n=',n,' r=',r     
        call quit("invalid parameters to n_choose_r")
      endif

      n_choose_r = n_choose_r_internal(n,min(r,n-r))

      return

      end

!     n_choose_r_internal(n,r) returns nCr

      function n_choose_r_internal(n,r)

      implicit none
      integer n, r, n_choose_r_internal
      integer i, ir
      integer col(r+1), colnext(r+1)

      col(1) = 1
      do ir = 1, r
         col(ir+1) = 0
      enddo
      do i = 1, n
         colnext(1) = 1
         do ir = 1, r
            colnext(ir+1) = col(ir) + col(ir+1)
         enddo
         do ir = 0, r
            col(ir+1) = colnext(ir+1)
         enddo
      enddo
      n_choose_r_internal = col(r+1)

      return

      end
