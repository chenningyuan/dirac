!-------------------------------------------------------------------------------------
!
! Library:      InteRest
!
! Module:       module_interest_one.f90 
!
! Description:  InteRest library -- one-electron integral routines 
!
! Licensing:    
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revision:     2.0   
!
!-------------------------------------------------------------------------------------
MODULE INTEREST

  implicit none

  private
 
  public :: interest_overlap

CONTAINS


  SUBROUTINE interest_overlap( gout,na,nb,nrkb,         &
                               la,alpha,ax,ay,az,anorm, &
                               lb,beta, bx,by,bz,bnorm  )

    implicit none
  
    !input
    integer, intent(in) :: la,lb 
    real(8), intent(in) :: alpha, anorm, ax, ay, az 
    real(8), intent(in) :: beta,  bnorm, bx, by, bz 
   
    !output
    integer, intent(out) :: na
    integer, intent(out) :: nb
    integer, intent(out) :: nrkb
    real(8), intent(out) :: gout(*)
  
    !local
    logical :: diag
    integer :: ij, losr, m, n, mlow
    integer :: i, li, ix, iy, iz, ic, istart, istop
    integer :: j, lj, jx, jy, jz, jc, jstart, jstop
    real(8) :: pi, pi34, pexp, rexp, o2p, piop, cntr, rxab, ryab, rzab
    real(8) :: pexp2, bop, bxop, byop, bzop, sxtwo, sytwo, sztwo, sxone, syone, szone 
        
    integer, parameter :: lmx   = 8
    integer, parameter :: lmxcc = 120
    integer            :: ncc (lmx) 
    integer            :: last(lmx) 
    integer            :: lin(lmxcc) 
    integer            :: min(lmxcc) 
    integer            :: nin(lmxcc) 
    real(8)            :: sx(2*lmx-1,lmx)
    real(8)            :: sy(2*lmx-1,lmx)
    real(8)            :: sz(2*lmx-1,lmx)
   
  
    !sanity checks
    if( max(la,lb)>lmx )then
      write(6,'(//2x,a//)') 'ERROR: the max. allowed angular quantum number exceeded => stop'
      call flush(6)
      stop
    endif
  
    !init
    call init()
    nrkb = 1
    na   = ncc(la)
    nb   = ncc(lb)
    pi   = 4.0d0*datan(1.0d0)    
    pi34 = (0.5d0/pi)**0.75d0    

    !bra/ket
    li     = la 
    lj     = lb 
    istart = last(li) + 1    
    jstart = last(lj) + 1 
    istop  = last(li) + ncc(li)
    jstop  = last(lj) + ncc(lj)

    !pre-calculations
    pexp = alpha + beta
    rexp = alpha*beta/pexp
    piop = dsqrt(pi/pexp)
    cntr = anorm*bnorm

    rxab = ax - bx 
    ryab = ay - by 
    rzab = az - bz 
    if( max(dabs(rxab),dabs(ryab),dabs(rzab))<1.d-12 )then
      rxab = 0.0d0 
      ryab = 0.0d0
      rzab = 0.0d0
      sx(1,1) = piop
      sy(1,1) = piop
      sz(1,1) = piop
    else
      sx(1,1) = dexp(-rexp*rxab*rxab)*piop
      sy(1,1) = dexp(-rexp*ryab*ryab)*piop
      sz(1,1) = dexp(-rexp*rzab*rzab)*piop
    endif


    losr = li + lj - 1
    if( losr.ge.2 )then

      pexp2 = 2.0d0*pexp
      bop   =-beta/pexp
      bxop  = bop*rxab
      byop  = bop*ryab
      bzop  = bop*rzab
      sxtwo = 0.0d0
      sytwo = 0.0d0
      sztwo = 0.0d0
      do m=2,losr,1  
        mlow = m - 1 
        o2p  = dble(mlow-1)/pexp2
        sxone   = sxtwo
        syone   = sytwo
        szone   = sztwo
        sxtwo   = sx(mlow,1)
        sytwo   = sy(mlow,1)
        sztwo   = sz(mlow,1)
        sx(m,1) = bxop*sxtwo + o2p*sxone
        sy(m,1) = byop*sytwo + o2p*syone
        sz(m,1) = bzop*sztwo + o2p*szone
      enddo         

      do m=2,lj,1
        mlow = m-1
        losr = losr-1 
        do n=1,losr,1   
          sx(n,m) = rxab*sx(n,mlow) + sx(n+1,mlow)
          sy(n,m) = ryab*sy(n,mlow) + sy(n+1,mlow)
          sz(n,m) = rzab*sz(n,mlow) + sz(n+1,mlow)
        enddo
      enddo
  
    endif 


    ij = 0
    do jc=jstart,jstop 
      jx = lin(jc)
      jy = min(jc)
      jz = nin(jc)
      do ic=istart,istop 
        ix = lin(ic)
        iy = min(ic)
        iz = nin(ic)
        ij = ij + 1
        gout(ij) = cntr*sx(ix,jx)*sy(iy,jy)*sz(iz,jz)
      enddo
    enddo

 CONTAINS
  
    subroutine init()
  
      !local
      integer :: l, lx, ly
      integer :: lval, mval, nval
      integer :: nc, nccgod, madr 
  
      nc     =  0
      nccgod =  0
      do l=1,lmx
        nc      = nc + l
        ncc(l)  = nc
        last(l) = nccgod
        nccgod  = nccgod + nc
      enddo
  
  
      madr = 0
      do l=0,lmx-1
        do lx=0,l
          lval = l - lx
          do ly=0,lx
            mval = lx - ly
            nval = ly
            madr      = madr + 1
            lin(madr) = lval + 1
            min(madr) = mval + 1
            nin(madr) = nval + 1
          enddo
        enddo
      enddo
  
    end subroutine
  
  END SUBROUTINE

END MODULE
