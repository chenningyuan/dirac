MODULE RECP_RAD1

CONTAINS

SUBROUTINE RECP_RAD1_MAIN(aa,aarr1,alpt,arp2,ccr,fctr2,  &
           kcrl,kcru,lamu,ltot1,ncr,rk,tol,zcr)
! compute type 1 radial integrals
  USE RECP_IPT 
  USE RECP_OUTPUT 
  implicit real(8) (a-h,o-z)
#include "inc_print.h"
  dimension ccr(*), ncr(*), zcr(*)
! global variables
  REAL(8) :: aa
! local variables
  INTEGER :: AREP_INDEX, AREP_N
  REAL(8) :: alpha

  CALL OUTPUT_LOC('RECP_RAD1_MAIN','E')
! ---------------------  
! AREP block(kcrl~kcru) 
! ---------------------  
  do AREP_INDEX =kcrl,kcru
!    AREP : AREP_C * r^(AREP_N-2) * exp(-AREP_Z*r*r)
     AREP_N = ncr(AREP_INDEX)     
!    alpha_j= alpha_A + alpha_B + AREP_Z ?? 
     alpha  = aa+zcr(AREP_INDEX)  

!    # exponential factor from q functions included in dum
     dum=aarr1+zcr(AREP_INDEX)*arp2/alpha
     if (dum.le.tol) then
 
        prd=fctr2*ccr(AREP_INDEX)*exp(-dum)
 
        if(rk.eq.0.0d0) then
          t=0.0d0
          do n=1,ltot1-mod(ltot1-1,2),2
            IPT_P1_QSUM(n,1) = IPT_P1_QSUM(n,1)  &
             + prd * qcomp(alpha,IPT_SH_DFAC,AREP_N+n-1,0,t,rk)
          enddo
        else
          t=alpt/alpha
          do lam=1,lamu
            IPT_P1_QSUM(lam,lam) = IPT_P1_QSUM(lam,lam)  &
             + prd * qcomp(alpha,IPT_SH_DFAC,AREP_N+lam-1,lam-1,t,rk)
          enddo
        endif
     endif
  enddo
 
! ---------------------  

  if(rk.ne.0.0d0) then
    f2lam3=(lamu+lamu-3)
    do lam=lamu-2,1,-1
      do n=lam+2,lamu-mod(lamu-lam,2),2
!              sum(d*Q)   =       sum(d*Q)     + ((2l+3)/k)* sum(d*Q)
        IPT_P1_QSUM(n,lam)=IPT_P1_QSUM(n,lam+2)+(f2lam3/rk)*IPT_P1_QSUM(n-1,lam+1)
      enddo
      f2lam3=f2lam3-2.0d0
    enddo
  endif

  CALL OUTPUT_LOC('RECP_RAD1_MAIN','X')

END SUBROUTINE RECP_RAD1_MAIN


END MODULE RECP_RAD1
