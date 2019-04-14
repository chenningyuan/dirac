module exacorr_cholesky

!routines to perform cholesky decomposition
!Johann Pototschnig, 04-2019

use, intrinsic:: ISO_C_BINDING

implicit none

complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

private

public cholesky_decompose

interface cholesky_decompose
   module procedure cholesky_decompose_real
end interface

contains

subroutine cholesky_decompose_real (matrix, mdim, CD_out, threshold)

  use exacorr_datatypes
  use exacorr_global

  implicit none

  !in: ao information
  real(8), intent(in)                  :: threshold
  real(8), intent(in)                  :: matrix(mdim,mdim)

  !out: 3 component tensors
  real(8), intent(out), allocatable :: CD_out(:,:)

  !work
  real(8), allocatable :: Dm(:)
  integer(8)           :: m, i, j, nonzero
  integer(8)           :: i_max
  real(8)              :: D_max,D_max_denom
  real(8), allocatable :: Gkl(:)
  real(8), allocatable :: CD(:,:)
  integer              :: mdim, log_conv

  !allocate tensors for cholesky decompostion
  allocate(CD(mdim,mdim))
  allocate(Dm(mdim))
  allocate(Gkl(mdim))

  nonzero = 0
  do j = 1, mdim
     do i = 1, mdim
        if (matrix(i,j).gt.threshold) nonzero = nonzero+1
     end do
  end do

  do i = 1, mdim
     Dm(i) = matrix(i,i)
  end do

  !find largest value
  i_max=maxloc(abs(Dm),1)
  D_max=Dm(i_max)
  write (*,'(a,a,i6,a,e12.2)') 'start iter',' largest element ind = ', i_max, ' value = ', abs(D_max)
  
  !loop until threshold is reached

  log_conv = ceiling(log10(abs(D_max)))
  CD=ZERO
  m=0
  do while (abs(D_max)>threshold)
    
    !compute next element
    do i = 1,m
       CD(:,m+1)=CD(:,m+1)+CD(:,i)*CD(i_max,i)
    end do
    CD(:,m+1)=matrix(:,i_max)-CD(:,m+1)
    D_max_denom = ONE / sqrt(D_max)
    CD(:,m+1)=CD(:,m+1) * D_max_denom

    m=m+1

    !update diagonal
    Gkl(:)=CD(:,m)*CD(:,m)
    Dm=Dm-Gkl

    !find largest element
    i_max=maxloc(abs(Dm),1)
    D_max=Dm(i_max)

    if (ceiling(log10(abs(D_max))) /= log_conv) then
       write (*,'(a,i6,a,i6,a,e12.2)') 'iter',m,' largest element ind = ', i_max, ' value = ', abs(D_max)
       log_conv = ceiling(log10(abs(D_max)))
    end if

    if (m.ge.mdim) then
      stop 'error: cholesky not converged'
    end if
  end do

  print *, 'Cholesky decomposition done, used ',m,' elements'
  print *, 'Matrix compressed to ',real(m)/real(mdim)*100,'% of the original size'
  print *, 'Original number of non-zero elements was ', real(nonzero)/real(mdim*mdim)*100,'%'

  deallocate(Gkl)
  deallocate(Dm)

  !hand tensor over
  allocate(CD_out(mdim, m))
  CD_out=CD(:,1:m)
  deallocate(CD)
end subroutine cholesky_decompose_real

subroutine cholesky_decompose_aomat (CD_out, threshold)

  use exacorr_datatypes
  use exacorr_global

  implicit none

  !in: ao information
  real(8), intent(in)                  :: threshold
  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells

  !out: 3 component tensors
  complex(8), intent(out), allocatable :: CD_out(:,:,:)

  !work
  complex(8), allocatable :: Dm(:,:)
  integer(8)              :: m, i, j
  integer(8)              :: i_max(1:2)
  complex(8)              :: D_max,D_max_denom
  complex(8), allocatable :: Gkl(:,:)
  complex(8), allocatable :: CD(:,:,:)

  !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
  call interest_initialize(.false.)

  !get all information about the Gauss type orbitals, allocate retrieve data
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)

  !allocate tensors for cholesky decompostion
  allocate(CD(nao*nao, nao, nao))
  allocate(Dm(nao,nao))
  allocate(Gkl(nao,nao))

  call get_2el_diag(nao, gto, nshells, Dm)

  !find largest value
  i_max=maxloc(abs(Dm))
  D_max=Dm(i_max(1), i_max(2))
  write (*,'(a,a,2i6,a,e12.2)') 'start iter',' largest element ind = ', i_max, ' value = ', abs(D_max)
  
  !loop until threshold is reached

  CD=ZERO
  m=0
  do while (abs(D_max)>threshold)
    
    !compute next element
    do i = 1,m
       CD(m+1,:,:)=CD(m+1,:,:)+CD(i,:,:)*CD(i,i_max(1),i_max(2))
    end do
    call get_2el_int(nao, gto, nshells, i_max(1), i_max(2), Gkl)
    CD(m+1,:,:)=Gkl-CD(m+1,:,:)
    D_max_denom = ONE / zsqrt(D_max)
    CD(m+1,:,:)=CD(m+1,:,:) * D_max_denom

    m=m+1

    !update diagonal
    Gkl(:,:)=CD(m,:,:)*CD(m,:,:)
    Dm=Dm-Gkl

    !find largest element
    i_max=maxloc(abs(Dm))
    D_max=Dm(i_max(1), i_max(2))
    write (*,'(a,i6,a,2i6,a,e12.2)') 'iter',m,' largest element ind = ', i_max, ' value = ', abs(D_max)

    if (m.ge.nao*nao) then
      stop 'error: cholesky not converged'
    end if
  end do

  print *, 'Cholesky dcomposition done '

  deallocate(Gkl)
  deallocate(Dm)
  deallocate(gto)

  !hand tensor over
  allocate(CD_out(m, nao, nao))
  CD_out=CD(1:m,:,:)
  deallocate(CD)
end subroutine cholesky_decompose_aomat

subroutine get_2el_diag(nao, gto, nshells, Dm)

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  !returns two electron integrals with i=k and j=l
  integer, intent(in) :: nao, nshells
  type(basis_func_info_t), intent(in) :: gto(:)  ! arrays with basis function information

  complex(8), intent(inout) :: Dm(:,:)

  integer :: basis_angular
  integer :: lq, lp, nq, np, q, p
  integer :: i, j, ijij
  integer :: ioff, joff
  real(8) :: eq,ep,cq,cp
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0

  basis_angular = get_basis_angular()
  joff = 0
  do q = 1, nshells
     lq   =  gto(q)%orb_momentum
     eq   =  gto(q)%exponent(1)
     xq   =  gto(q)%coord(1)
     yq   =  gto(q)%coord(2)
     zq   =  gto(q)%coord(3)
     cq   =  gto(q)%coefficient(1)
     nq   =  nfunctions(lq,basis_angular)

     ioff = 0
     do p = 1, nshells
        lp   =  gto(p)%orb_momentum
        ep   =  gto(p)%exponent(1)
        xp   =  gto(p)%coord(1)
        yp   =  gto(p)%coord(2)
        zp   =  gto(p)%coord(3)
        cp   =  gto(p)%coefficient(1)
        np   =  nfunctions(lp,basis_angular)

        !get eri values
        call interest_eri('llll',FIJKL,gout,&
                          lp,ep,xp,yp,zp,cp,&
                          lq,eq,xq,yq,zq,cq,&
                          lp,ep,xp,yp,zp,cp,&
                          lq,eq,xq,yq,zq,cq )

        do j = 1, nq
          do i = 1, np
            ijij = (j-1)*np*nq*np+(i-1)*nq*np+(j-1)*np+i
            Dm(i+ioff,j+joff) = dcmplx(gout(ijij),0.D0)
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine get_2el_diag

subroutine get_2el_int(nao, gto, nshells, k_fix, l_fix, Gkl)
  !returns two electron integrals with two fixed indices

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  integer, intent(in)                 :: nao, nshells
  integer(8), intent(in)              :: k_fix, l_fix
  type(basis_func_info_t), intent(in) :: gto(:)  ! arrays with basis function information

  complex(8), intent(inout) :: Gkl(:,:)

  integer      :: basis_angular
  integer      :: i, j, k, l, ioff,joff, koff, loff
  integer      :: ijkl
  integer      :: ish, jsh, ksh, lsh
  integer      :: li,lj,lk,ll, ni,nj,nk,nl
  real(8)      :: ei,ej,ek,el,ci,cj,ck,cl
  real(8)      :: xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0

  basis_angular = get_basis_angular()
  loff = 0
! This search could be avoided by generalizing the routine that points from a function to the shell
  do lsh = 1, nshells
    ll   =  gto(lsh)%orb_momentum
    el   =  gto(lsh)%exponent(1)
    xl   =  gto(lsh)%coord(1)
    yl   =  gto(lsh)%coord(2)
    zl   =  gto(lsh)%coord(3)
    cl   =  gto(lsh)%coefficient(1)
    nl   =  nfunctions(ll,basis_angular)
    if ( l_fix>loff .and. l_fix.le.loff+nl) then
        l = l_fix - loff
        EXIT
    else
       loff = loff + nl
    end if
  end do

  koff = 0
  do ksh = 1, nshells
    lk   =  gto(ksh)%orb_momentum
    ek   =  gto(ksh)%exponent(1)
    xk   =  gto(ksh)%coord(1)
    yk   =  gto(ksh)%coord(2)
    zk   =  gto(ksh)%coord(3)
    ck   =  gto(ksh)%coefficient(1)
    nk   =  nfunctions(lk,basis_angular)
    if ( k_fix>koff .and. k_fix.le.koff+nk) then
        k = k_fix - koff
        EXIT
    else
       koff = koff + nk
    end if
  end do

  joff = 0
  do jsh = 1, nshells
    lj   =  gto(jsh)%orb_momentum
    ej   =  gto(jsh)%exponent(1)
    xj   =  gto(jsh)%coord(1)
    yj   =  gto(jsh)%coord(2)
    zj   =  gto(jsh)%coord(3)
    cj   =  gto(jsh)%coefficient(1)
    nj   =  nfunctions(lj,basis_angular)
    ioff = 0
    do ish = 1, nshells
      li   =  gto(ish)%orb_momentum
      ei   =  gto(ish)%exponent(1)
      xi   =  gto(ish)%coord(1)
      yi   =  gto(ish)%coord(2)
      zi   =  gto(ish)%coord(3)
      ci   =  gto(ish)%coefficient(1)
      ni   =  nfunctions(li,basis_angular)

!     output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
      call interest_eri('llll',fijkl,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )
      do j = 1, nj
        do i = 1, ni
           ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
           Gkl(i+ioff,j+joff) = dcmplx(gout(ijkl),0.D0)
        end do
      end do
      
      ioff = ioff + ni
    end do
    joff = joff + nj
  end do
  
end subroutine get_2el_int


end module exacorr_cholesky
