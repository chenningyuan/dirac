MODULE RECP_INT_C2

CONTAINS

SUBROUTINE RECP_INT_C2_RAD1(aa,aarr1,alpt,arp2,ccr,dfac,fctr2,kcrl,kcru,lamu,ltot1,ncr,qsum,rk,tol,zcr)
 
! compute type 1 radial integrals
 
  implicit real(8) (a-h,o-z)
  dimension ccr(*), dfac(*), ncr(*), qsum(ltot1,*), zcr(*)

! variables : global
  REAL(8) :: aa

! variables : local
  INTEGER :: AREP_INDEX, AREP_N
  REAL(8) :: alpha

! IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : entering... RECP_INT_C2_RAD1'
  PRINT *,' * ECPINT : entering... RECP_INT_C2_RAD1'

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
            qsum(n,1)=qsum(n,1) + prd * qcomp(alpha,dfac,AREP_N+n-1,0,t,rk)
          enddo
        else
          t=alpt/alpha
          do lam=1,lamu
            qsum(lam,lam) = qsum(lam,lam) + prd * qcomp(alpha,dfac,AREP_N+lam-1,lam-1,t,rk)
          enddo
        endif
     endif
  enddo
 
! ---------------------  

  if(rk.ne.0.0d0) then
    f2lam3=(lamu+lamu-3)
    do lam=lamu-2,1,-1
      do n=lam+2,lamu-mod(lamu-lam,2),2
!       sum(d*Q)    = sum(d*Q)      + ((2l+3)/k)* sum(d*Q)
        qsum(n,lam) = qsum(n,lam+2) + (f2lam3/rk)*qsum(n-1,lam+1)
      enddo
      f2lam3=f2lam3-2.0d0
    enddo
  endif

  IF (PRINT_LEVEL .GE. 5) PRINT *,' * ECPINT : exiting...  RECP_INT_C2_RAD1'

END SUBROUTINE RECP_INT_C2_RAD1

! -----------------------------

SUBROUTINE RECP_INT_C2_ANG1(ang,dfac,nanb,lalb,mamb,lamu,lmf,lml,lmx,lmy,lmz,ltot1,xab,yab,zab,xk,yk,zk,zlm)
!
!  compute type 1 angular integrals
!
      implicit real(8) (a-h,o-z)
      dimension ang(ltot1,*),dfac(*),lmf(*),lml(*),lmx(*),lmy(*),lmz(*),xab(*),yab(*),zab(*),zlm(*)
 
      call wzero(ltot1*lamu,ang,1)
      do 96 n=1,nanb
      if(xab(n).eq.0.0d0) go to 96
      do 94 l=1,lalb
      if(yab(l).eq.0.0d0) go to 94
      do 92 m=1,mamb
      if(zab(m).eq.0.0d0) go to 92
      nlm=((n-2)+l)+m
      lamlo=mod(nlm-1,2)+1
      lamhi=min(nlm,lamu)
      if(lamlo.gt.lamhi) go to 92
      do 90 lam=lamlo,lamhi,2
      l2=lam+lam-1
      angt=0.0d0
      loc=(lam-1)**2
      do mu1=1,l2
         istart=lmf(loc+mu1)
         if(mod(n,2).ne.mod(lmx(istart),2) .and. &
            mod(l,2).ne.mod(lmy(istart),2) .and. &
            mod(m,2).ne.mod(lmz(istart),2)) then
            pre=0.0d0
            aint=0.0d0
            iend=lml(loc+mu1)
            do i=istart,iend
               indx=lmx(i)
               indy=lmy(i)
               indz=lmz(i)
               if(indx.eq.0) then
                 xkp=1.0d0
               else
                 xkp=xk**indx
               endif
               if(indy.eq.0) then
                 ykp=1.0d0
               else
                 ykp=yk**indy
               endif
               if(indz.eq.0) then
                 zkp=1.0d0
               else
                 zkp=zk**indz
               endif
               pre=pre+zlm(i)*xkp*ykp*zkp
               aint = aint + zlm(i)*dfac(n+indx)*dfac(l+indy)*dfac(m+indz) &
                             / dfac((n+indx)+(l+indy)+(m+indz))
            enddo
            angt=angt+pre*aint
         endif
      enddo
      ang(nlm,lam)=ang(nlm,lam)+((xab(n)*yab(l))*zab(m))*angt
   90 continue
   92 continue
   94 continue
   96 continue
      return
END SUBROUTINE RECP_INT_C2_ANG1

END MODULE RECP_INT_C2
