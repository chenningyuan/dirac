MODULE RECP_FUNCTION

CONTAINS

!SUBROUTINE RECP_ONEINT_C_CORE
!END SUBROUTINE RECP_ONEINT_C_CORE

!SUBROUTINE RECP_ONEINT_C_PSEUDO1 &
!          (a,ang,ccr,crda,crdb,gout,ipt,lmnv,ltot1,ncr,nkcrl,nkcru,q,qsum,xab,yab,zab,zcr)
!END SUBROUTINE RECP_ONEINT_C_PSEUDO1

!SUBROUTINE RECP_ONEINT_C_PSEUDO1_RAD
!SUBROUTINE RECP_ONEINT_C_PSEUDO1_ANG

SUBROUTINE RECP_FUNCTION_FACAB (binom,An,Bn,crda,crdb,xab)

  IMPLICIT NONE

! variables : global
  INTEGER :: An, Bn
! dimension binom(*), crda(*), crdb(*), xab(*)
  REAL(8) :: binom
  REAL(8) :: crda
  REAL(8) :: crdb
  REAL(8) :: xab
  DIMENSION binom(*), crda(*), crdb(*), xab(*)

! variables : local
  INTEGER :: Ai, SUM_A
  INTEGER :: Bi, SUM_B

  CALL RECP_FUNCTION_WZERO(An+Bn-1, xab, 1)

  SUM_A=(An*(An-1))/2  ! sum from 1 to An
  SUM_B=(Bn*(Bn-1))/2
  do Ai=1,An
    do Bi=1,Bn
      xab((Ai-1)+Bi) = xab((Ai-1)+Bi) &
                       + binom(SUM_A+Ai) * crda((An+1)-Ai) * binom(SUM_B+Bi) * crdb((Bn+1)-Bi)
    enddo
  enddo
END SUBROUTINE RECP_FUNCTION_FACAB


SUBROUTINE RECP_FUNCTION_WZERO (n,a,inca)
! write zero to a(1+inca*i), i=0~n-1

  IMPLICIT NONE
! variables : global
  INTEGER :: n, inca
  REAL(8) :: a(*)
! variables : local
  INTEGER :: ia, i

  ia = 1
  do i = 1, n
     a(ia) = 0.d0
     ia    = ia + inca
  enddo
END SUBROUTINE RECP_FUNCTION_WZERO

!------------------------------------------------------
SUBROUTINE RECP_FUNCTION_QCOMP(alpha,dfac,n,l,t,xk)

!  compute q(n,l) scaled by sqrt(pi)*exp(-t) to prevent overflows
!  arguments are alpha, xk, and t=xk**2/(4*alpha)
!  no restriction on the magnitude of t
!  increase dfac array to raise n, l restrictions
 
      implicit real*8 (a-h,o-z)
      parameter (accpow=1.0d-14, accasy=1.0d-10)
      dimension dfac(*)
      dimension tmin(0:8)
      data tmin/31.0d0,28.0d0,25.0d0,23.0d0,22.0d0,20.0d0,19.0d0,18.0d0,15.0d0/
 
      if(mod(n+l,2).ne.0.or.n.le.l) go to 30

!  ----------------------------------------------
!  use alternating series (n+l.le.22.and.l.le.10)
!  ----------------------------------------------
      if(l.eq.0) then
        xkp=1.0d0
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=xkp * dfac(n+l+1)/ ((alpha+alpha)**((n-l)/2) * sqrt(4.0d0*alpha) * dfac(l+l+3))
      num=l-n+2
      xden=(l+l+3)
      term=1.0d0
      sum=term
      xc=-1.0d0
   10 if(num.ne.0) then
        fnum=num
        term=term*fnum*t/(xden*xc)
        xc=xc-1.0d0
        sum=sum+term
        num=num+2
        xden=xden+2.0d0
        go to 10
      endif
      qcomp=prefac*sum
      return
 
   30 if(t.lt.tmin(min(n,8))) go to 60
 
!  -------------------------------------
!  use asymptotic series (arbitrary n,l)
!  -------------------------------------
      xkp=(xk/(alpha+alpha))**(n-2)
      prefac=xkp/((alpha+alpha)*sqrt(4.0d0*alpha))
      sum=1.0d0
      term=1.0d0
      fac1=(l-n+2)
      fac2=(1-l-n)
      xc=1.0d0
   40 term=term*fac1*fac2/(4.0d0*xc*t)
      if(term.eq.0.0d0) go to 50
      sum=sum+term
      if(abs(term/sum).lt.accasy) go to 50
      fac1=fac1+2.0d0
      fac2=fac2+2.0d0
      xc=xc+1.0d0
      go to 40
   50 qcomp=prefac*sum
      return

!  ----------------------------------------
!  use power series (n+l.le.22.and.l.le.10)
!  ----------------------------------------
   60 if(l.eq.0) then
        xkp=1.0d0
      else
        xkp=(xk/(alpha+alpha))**l
      endif
      prefac=exp(-t)*xkp/(alpha+alpha)**((n-l+1)/2)
      if(mod(n+l,2).eq.0) then
        prefac=prefac/sqrt(4.0d0*alpha)
      else
        prefac=0.56418958354775629d0*prefac
      endif
      xnum=(l+n-1)
      xden=(l+l+1)
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=0.0d0
   70 xnum=xnum+2.0d0
      xden=xden+2.0d0
      xj=xj+1.0d0
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if((term/sum).gt.accpow) go to 70
      qcomp=prefac*sum
      return

END SUBROUTINE RECP_FUNCTION_QCOMP

FUNCTION RECP_FNIPQ(VALUE_IN)
! From ipq
!  integer        ipq
!  common /dim21/ ipq(1000) !256->1000 
  IMPLICIT NONE
  INTEGER  :: VALUE_IN 
  INTEGER  :: RECP_FNIPQ 
  RECP_FNIPQ = (VALUE_IN * (VALUE_IN - 1) ) / 2
  RETURN
END FUNCTION RECP_FNIPQ

END MODULE RECP_FUNCTION
