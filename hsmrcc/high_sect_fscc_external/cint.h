!
! arrays for integrals (vint) and indexes (iint)
! nint - number of stored integrals & indexes
!

      integer nint,idum1,            iint(ndint)
      real*8             vint(ndint)

      common /cint/ nint,idum1,vint,iint

!       common /cint/ nint,idum1,vint(ndint),iint(ndint)
