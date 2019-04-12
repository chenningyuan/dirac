module exacorr_cholesky

!routines to perform cholesky decomposition
!Johann Pototschnig, 04-2019

use, intrinsic:: ISO_C_BINDING

implicit none

complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

private

public cholesky_decompose

contains

subroutine cholesky_decompose (CD_out, threshold)

  use exacorr_datatypes
  use exacorr_global

  implicit none

  !in: ao information
  real(8), intent(in)                  :: threshold
  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular ! 1=cartesian, 2=spherical

  !out: 3 component tensors
  complex(8), intent(out), allocatable :: CD_out(:,:,:)

  !work
  complex(8), allocatable :: Dm(:,:)
  integer(8)              :: m, i, j
  integer(8)              :: i_max(1:2)
  complex(8)              :: D_max
  complex(8), allocatable :: Res(:,:), Gkl(:,:)
  complex(8), allocatable :: CD(:,:,:)

  !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
  call interest_initialize(.false.)

  !get all information about the Gauss type orbitals, allocate retrieve data
  nao = get_nao()
  basis_angular = get_basis_angular()
  call get_gtos(1,nao,gto,nshells)

  !allocate tensors for cholesky decompostion
  allocate(CD(nao*nao, nao, nao))
  CD=ZERO
  allocate(Dm(nao,nao))
  allocate(Res(nao,nao))
  allocate(Gkl(nao,nao))

  call get_2el_diag(nao, gto, basis_angular, nshells, Dm)

  !find largest value
  m=0
! print *, ' -------- iteration m = ', m, ' ------------ '
  i_max=maxloc(abs(Dm))
  D_max=Dm(i_max(1), i_max(2))
  print *, ' largest element ind = ', i_max, ' value = ', D_max
  write (*,'(a,a,2i6,a,e12.2)') 'start iter',' largest element ind = ', i_max, ' value = ', abs(D_max)
  
  !loop until threshold is reached

  do while (abs(D_max)>threshold)
    
    !compute next element
    Res=ZERO
    do i = 1,m
       Res=Res+CD(i,:,:)*CD(i,i_max(1),i_max(2))
    end do
    call get_2el_int(nao, gto, basis_angular, nshells, i_max(1), i_max(2), Gkl)
    Res=Gkl-Res
    Res=Res/zsqrt(D_max)

    m=m+1
!   print *, ' -------- iteration m = ', m, ' ------------ '

    !save vector
    CD(m,:,:)=Res

    !update diagonal
    Gkl=ZERO
    do i = 1,nao
      do j= 1, nao
        Gkl(i,j)=Res(i,j)*Res(i,j)
      end do
    end do
    Dm=Dm-Gkl

    !find largest element
!   print *, 'old one =', Dm(i_max(1), i_max(2))
    i_max=maxloc(abs(Dm))
    D_max=Dm(i_max(1), i_max(2))
    write (*,'(a,i6,a,2i6,a,e12.2)') 'iter',m,' largest element ind = ', i_max, ' value = ', abs(D_max)

    if (m.ge.nao*nao) then
      stop 'error: cholesky not converged'
    end if
  end do

  print *, 'Cholesky dcomposition done '

  deallocate(Gkl)
  deallocate(Res)
  deallocate(Dm)
  deallocate(gto)

  !hand tensor over
  allocate(CD_out(m, nao, nao))
  do i = 1,m
       CD_out(i,:,:)=CD(i,:,:)
  end do
  deallocate(CD)
end subroutine cholesky_decompose

subroutine get_2el_diag(nao, gto, basis_angular, nshells, Dm)

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  !returns two electron integrals with i=k and j=l
  integer, intent(in) :: nao, basis_angular, nshells
  type(basis_func_info_t), intent(in) :: gto(:)  ! arrays with basis function information

  complex(8), intent(inout) :: Dm(:,:)

  integer :: lq, lp, nq, np, q, p
  integer :: i, j, k, l, ijkl
  integer :: ioff, joff, koff, loff
  real(8) :: eq,ep,cq,cp
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0

  joff = 0
  loff = 0
  do q = 1, nshells
     lq   =  gto(q)%orb_momentum
     eq   =  gto(q)%exponent(1)
     xq   =  gto(q)%coord(1)
     yq   =  gto(q)%coord(2)
     zq   =  gto(q)%coord(3)
     cq   =  gto(q)%coefficient(1)
     nq   =  nfunctions(lq,basis_angular)

     ioff = 0
     koff = 0
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

        ijkl=0
        do  l = 1, nq
          do k = 1, np
            do j = 1, nq
              do i = 1, np
                ijkl = ijkl+1

                if (i==k) then
                  if (j==l) Dm(i+ioff,j+joff) = dcmplx(gout(ijkl),0.D0)
                end if

              end do
            end do
          end do
        end do

        ioff = ioff + np
        koff = koff + np
     end do
     joff = joff + nq
     loff = loff + nq
  end do

end subroutine get_2el_diag

subroutine get_2el_int(nao, gto, basis_angular, nshells, l_fix, k_fix, Gkl)
  !returns two electron integrals with two fixed indeces

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  integer, intent(in)                 :: nao, basis_angular, nshells
  integer(8), intent(in)              :: l_fix, k_fix
  type(basis_func_info_t), intent(in) :: gto(:)  ! arrays with basis function information

  complex(8), intent(inout) :: Gkl(:,:)



  integer      :: i, j, k, l, ioff,joff, koff, loff
  integer      :: ijkl
  integer      :: ish, jsh, ksh, lsh
  integer      :: li,lj,lk,ll, ni,nj,nk,nl
  real(8)      :: ei,ej,ek,el,ci,cj,ck,cl
  real(8)      :: xi,yi,zi,xj,yj,zj,xk,yk,zk,xl,yl,zl
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0

  loff = 0
  do lsh = 1, nshells
    ll   =  gto(lsh)%orb_momentum
    el   =  gto(lsh)%exponent(1)
    xl   =  gto(lsh)%coord(1)
    yl   =  gto(lsh)%coord(2)
    zl   =  gto(lsh)%coord(3)
    cl   =  gto(lsh)%coefficient(1)
    nl   =  nfunctions(ll,basis_angular)
    if ( l_fix>loff .and. l_fix.le.loff+nl) then
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
      ijkl = 0
      do l = 1, nl
        do k = 1, nk
          do j = 1, nj
            do i = 1, ni
              ijkl = ijkl + 1
              if (koff+k==k_fix) then
                if (loff+l==l_fix) then
                  Gkl(i+ioff,j+joff) = dcmplx(gout(ijkl),0.D0)
                end if
              end if
            end do
          end do
        end do
      end do
      
      ioff = ioff + ni
    end do
    joff = joff + nj
  end do
  
end subroutine get_2el_int


end module exacorr_cholesky
