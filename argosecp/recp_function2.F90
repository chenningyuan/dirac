MODULE RECP_FUNCTION2

CONTAINS

FUNCTION FN_IODFAC(l1,l2,m1,m2,n1,n2)
!  compute the product of odd number factorials which gives the
!  one-center overlaps of cartesian gaussian aos.

  IMPLICIT NONE

! vaiables : global
  INTEGER :: l1, l2, m1, m2, n1, n2
  INTEGER :: FN_IODFAC
! vaiables : local
  INTEGER :: i, iprd

  if(mod(l1+l2,2).eq.1 .or.  mod(m1+m2,2).eq.1 .or.  mod(n1+n2,2).eq.1) then
    fn_iodfac = 0
  else
    iprd = 1
    do i=3,l1+l2-1,2
       iprd = i*iprd
    enddo
    do i=3,m1+m2-1,2
      iprd = i*iprd
    enddo
    do i=3,n1+n2-1,2
      iprd = i*iprd
    enddo
    fn_iodfac = iprd
  endif

END FUNCTION FN_IODFAC


SUBROUTINE RECP_FUNCTION_FACAB (binom,An,Bn,crda,crdb,xab)

END SUBROUTINE RECP_FUNCTION_FACAB


END MODULE RECP_FUNCTION2
