! -------------------------------------------------------------------------
!
! Program:      Interest
!
! File:         module_interest.f90 
!
! Description:  Performing integral calculations over normalized cartesian/spherical
!               unbalanced, or kinetically/magnetically balanced gaussian type functions
!
! Licensing:    This code is distributed under the GNU LGPL license
!
! Author:       Michal Repisky (michal.repisky@uit.no)
!
! Revisions:    
!
! -------------------------------------------------------------------------
  subroutine interest_initialize()

    use module_interest_osr
    implicit none

    if( is_interest_initialized )then
      write(6,*)
      write(6,*)'Module "interest" has already been initialized ...'
      write(6,*)
    else
      write(6,*)
      write(6,*)'Module "interest" is initializing now ...'
      write(6,*)
      call interest_osr_initialize()
      is_interest_initialized = .true.
    endif

  end subroutine

! -------------------------------------------------------------------------
  subroutine interest_overlap(fint,nint,              &
                              la,alpha,ax,ay,az,anorm,&
                              lb,beta, bx,by,bz,bnorm,&
                              rkb)
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: nint
    real(8), intent(out) :: fint(*)
    !-- input --!
    logical, intent(in) :: rkb 
    integer, intent(in) :: la,lb
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    !-- local --!
    integer :: n1
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: pexp,rexp,rxpa,rypa,rzpa,cntr


    if( .not.rkb )then
      nint = ncc(la)*ncc(lb)
    else
      nint = ncc(la)*ncc(lb)*10
    endif
    fint(1:nint) = 0.0d0

    rxab = ax - bx
    ryab = ay - by
    rzab = az - bz
    pexp = alpha + beta
    rrab = rxab*rxab + ryab*ryab + rzab*rzab

! todo: set the prescreening here

    if( rrab.lt.1.d-12 )then
      rxpa = 0.0d0 
      rypa = 0.0d0
      rzpa = 0.0d0
      rxab = 0.0d0  
      ryab = 0.0d0  
      rzab = 0.0d0  
      cntr = anorm*bnorm 
    else
      rxpa = ( alpha*ax + beta*bx )/pexp - ax
      rypa = ( alpha*ay + beta*by )/pexp - ay
      rzpa = ( alpha*az + beta*bz )/pexp - az
      rexp = alpha*beta/pexp
      cntr = anorm*bnorm*dexp(-rexp*rrab)
    endif

    if( .not.rkb )then
      call interest_osr_class_overlap(n1,fint,cntr,(la  ),(la+lb-1),pexp,rxpa,rypa,rzpa)
      call interest_hrr_bra(fint,nint,1,la,lb,rxab,ryab,rzab)
    else
      call interest_osr_class_overlap(n1,fint,cntr,(la-1),(la+lb+1),pexp,rxpa,rypa,rzpa)
      call interest_hrr_bra_rkb(fint,nint,1,la,lb,rxab,ryab,rzab,alpha,beta)
    endif

  end subroutine

! -------------------------------------------------------------------------
  subroutine interest_nuclear_attraction_point_nucleus(fint,nint,                     &
                                                       la,alpha,ax,ay,az,anorm,       &
                                                       lb,beta, bx,by,bz,bnorm,       &
                                                       ncentr,nuc_charge,coord_of_nuc,&
                                                       rkb)
    use module_interest_osr
    use module_interest_hrr
    implicit none

    !-- output --!
    integer, intent(out) :: nint
    real(8), intent(out) :: fint(*)
    !-- input --!
    logical, intent(in) :: rkb 
    integer, intent(in) :: la,lb
    real(8), intent(in) :: alpha,ax,ay,az,anorm
    real(8), intent(in) :: beta, bx,by,bz,bnorm
    integer, intent(in) :: ncentr 
    integer, intent(in) :: nuc_charge(ncentr) 
    real(8), intent(in) :: coord_of_nuc(3,ncentr) 
    !-- local --!
    integer :: lab,n1,n2,iat
    real(8) :: rxab,ryab,rzab,rrab
    real(8) :: pexp,rexp,rxpa,rypa,rzpa
    real(8) :: cx,cy,cz,px,py,pz,rxcp,rycp,rzcp
    real(8) :: cntr1,cntr,tval,fnbra,fdbra
    real(8), dimension(ncc2) :: fintl


    if( .not.rkb )then
      nint = ncc(la)*ncc(lb)
    else
      nint = ncc(la)*ncc(lb)*10
    endif
    fint(1:nint) = 0.0d0

    rxab = ax - bx
    ryab = ay - by
    rzab = az - bz
    pexp = alpha + beta
    rrab = rxab*rxab + ryab*ryab + rzab*rzab

! todo: set the prescreening here

    if( rrab.lt.1.d-12 )then
      px    = ax 
      py    = ay
      pz    = az
      rxpa  = 0.0d0 
      rypa  = 0.0d0
      rzpa  = 0.0d0
      rxab  = 0.0d0  
      ryab  = 0.0d0  
      rzab  = 0.0d0  
      cntr1 = 2.0d0*(pi/pexp)*anorm*bnorm 
    else
      px    = ( alpha*ax + beta*bx )/pexp
      py    = ( alpha*ay + beta*by )/pexp
      pz    = ( alpha*az + beta*bz )/pexp
      rxpa  = px - ax
      rypa  = py - ay
      rzpa  = pz - az
      rexp  = alpha*beta/pexp
      cntr1 = 2.0d0*(pi/pexp)*anorm*bnorm*dexp(-rexp*rrab)
    endif

    fnbra = 1.0d0 
    fdbra = 1.0d0/(2.0d0*pexp)  

    do iat=1,ncentr
      if( nuc_charge(iat) <= 0 )cycle
      cx = coord_of_nuc(1,iat)
      cy = coord_of_nuc(2,iat)
      cz = coord_of_nuc(3,iat)
      rxcp = cx - px
      rycp = cy - py
      rzcp = cz - pz
      cntr = dfloat(nuc_charge(iat))*cntr1 
      tval = pexp*(rxcp*rxcp + rycp*rycp + rzcp*rzcp)
      if( .not.rkb )then
        call interest_osr_class_nuclear(n1,n2,fintl,cntr,tval,(la  ),(la+lb-1),       &
                                        1,fnbra,fdbra,rxpa,rypa,rzpa,rxcp,rycp,rzcp,  &
                                        1,1,1,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      else
        call interest_osr_class_nuclear(n1,n2,fintl,cntr,tval,(la-1),(la+lb+1),       &
                                        1,fnbra,fdbra,rxpa,rypa,rzpa,rxcp,rycp,rzcp,  &
                                        1,1,1,0.0d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0)
      endif
      fint(1:n1) = fint(1:n1) - fintl(1:n1)                 
    enddo  

    if( .not.rkb )then
      call interest_hrr_bra(fint,nint,1,la,lb,rxab,ryab,rzab)
    else
      call interest_hrr_bra_rkb(fint,nint,1,la,lb,rxab,ryab,rzab,alpha,beta)
    endif
  end subroutine

! --------------------------------------------------------------------------
