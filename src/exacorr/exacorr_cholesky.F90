module exacorr_cholesky

!routines to perform cholesky decomposition
!Johann Pototschnig, 04-2019

use, intrinsic:: ISO_C_BINDING

implicit none

complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

private

public talsh_get_diag_real
public talsh_get_diag_comp
public talsh_find_shells_real
public talsh_find_shells_comp
public talsh_get_shell_real
public talsh_get_shell_comp
public talsh_find_el_real
public talsh_find_el_comp
public talsh_square_el
public talsh_set_mask_done

public get_offset
public get_shellsize

contains

subroutine talsh_get_diag_real(Diag)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag

  real(8), pointer     :: Dm(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2)
  integer              :: ierr

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer :: lq, lp, nq, np, q, p
  integer :: i, j, ijij
  integer :: ioff, joff
  real(8) :: eq,ep,cq,cp
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0
  
  ! get acces to the tensor
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Diag,body_p,R8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  ! compute diagonal elements
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
            Dm(i+ioff,j+joff) = gout(ijij)
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine talsh_get_diag_real

subroutine talsh_get_diag_comp(Diag)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag

  complex(8), pointer  :: Dm(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2)
  integer              :: ierr

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer :: lq, lp, nq, np, q, p
  integer :: i, j, ijij
  integer :: ioff, joff
  real(8) :: eq,ep,cq,cp
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0
  
  ! get acces to the tensor
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Diag,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  ! compute diagonal elements
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

end subroutine talsh_get_diag_comp


subroutine talsh_find_shells_real(Diag, Filter, lshell, kshell, lsize, ksize, D_max)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag
  type(talsh_tens_t), intent(inout) :: Filter
  integer, intent(out)              :: lshell, kshell, lsize, ksize
  real(8), intent(out)              :: D_max


  real(8), pointer     :: Dm(:,:)
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), i_max(1:2)
  real(8)              :: D_tot
  integer              :: ierr

  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Diag,body_p,R8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=rfilt.gt.0

  !find largest element
  i_max=maxloc(abs(Dm))
  D_tot=Dm(i_max(1), i_max(2))

  !find largest non treated element
  i_max=maxloc(abs(Dm), mask=lfilt)
  D_max=Dm(i_max(1), i_max(2))

  if (abs(D_tot).gt.abs(D_max)) then
    print *, 'Warning: error in cholesky will not be smaller than ', D_tot
  end if

  lshell=get_shell_index(i_max(1))
  kshell=get_shell_index(i_max(2))

  lsize=get_shell_size(lshell)
  ksize=get_shell_size(kshell)

  deallocate(lfilt)

end subroutine talsh_find_shells_real

subroutine talsh_find_shells_comp(Diag, Filter, lshell, kshell, D_max)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag
  type(talsh_tens_t), intent(inout) :: Filter
  integer, intent(out)              :: lshell, kshell
  complex(8), intent(out)           :: D_max


  complex(8), pointer  :: Dm(:,:)
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), i_max(1:2)
  complex(8)           :: D_tot
  integer              :: ierr

  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Diag,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=rfilt.gt.0

  !find largest element
  i_max=maxloc(abs(Dm))
  D_tot=Dm(i_max(1), i_max(2))

  !find largest non treated element
  i_max=maxloc(abs(Dm), mask=lfilt)
  D_max=Dm(i_max(1), i_max(2))

  if (abs(D_tot).gt.abs(D_max)) then
    print *, 'Warning: error in cholesky will not be smaller than ', D_tot
  end if

  lshell=get_shell_index(i_max(1))
  kshell=get_shell_index(i_max(2))

  deallocate(lfilt)

end subroutine talsh_find_shells_comp


subroutine talsh_get_shell_real(Dshell, lsh, ksh)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none 

  integer, intent(in)            :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Dshell

  real(8), pointer     :: Dlk(:,:,:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims4(1:4)
  integer              :: ierr


  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  !variables required to compute 2el integrals
  integer(8)         :: ijkl
  integer(8)         :: ioff,joff,koff,loff
  integer(8)         :: ish,jsh
  integer(8)         :: kshell,lshell
  integer            :: i, ni, li
  real(8)            :: ei, ci, xi, yi, zi
  integer            :: j, nj, lj
  real(8)            :: ej, cj, xj, yj, zj
  integer            :: k, nk, lk
  real(8)            :: ek, ck, xk, yk, zk
  integer            :: l, nl, ll
  real(8)            :: el, cl, xl, yl, zl
  real(8), parameter :: FIJKL = 1.0
  integer, parameter :: LMAX = 21 ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), save      :: gout(LMAX*LMAX*LMAX*LMAX) 
  real(8)            :: limit

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !set parameters for the fixed shell
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  el   =  gto(lsh)%exponent(1)
  xl   =  gto(lsh)%coord(1)
  yl   =  gto(lsh)%coord(2)
  zl   =  gto(lsh)%coord(3)
  cl   =  gto(lsh)%coefficient(1)
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  ek   =  gto(ksh)%exponent(1)
  xk   =  gto(ksh)%coord(1)
  yk   =  gto(ksh)%coord(2)
  zk   =  gto(ksh)%coord(3)
  ck   =  gto(ksh)%coefficient(1)
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  rank=talsh_tensor_rank(Dshell)
  if (rank.ne.4) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Dshell,rank,dims4)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Dshell,body_p,R8,0,DEV_HOST)
  call c_f_pointer(body_p,Dlk,dims4)

  !loop over all the shells for the other indices
  joff = 0
  do jsh = 1, nshells
    lj   =  gto(jsh)%orb_momentum
    ej   =  gto(jsh)%exponent(1)
    xj   =  gto(jsh)%coord(1)
    yj   =  gto(jsh)%coord(2)
    zj   =  gto(jsh)%coord(3)
    cj   =  gto(jsh)%coefficient(1)
    nj   =  nfunctions(lj,basis_angular)

    !neglect small elements
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
      call interest_eri('llll',FIJKL,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )

      !save all elements for this shell
      do k = 1, nk
        do l = 1, nl
          do j = 1, nj
            do i = 1, ni
              ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
              Dlk(i+ioff,j+joff,l,k) = gout(ijkl)
            end do
          end do
        end do
      end do

      ioff = ioff + ni
    end do

    joff = joff + nj
  end do

end subroutine talsh_get_shell_real

subroutine talsh_get_shell_comp(Dshell, lsh, ksh)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none 

  integer, intent(in)            :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Dshell

  complex(8), pointer  :: Dlk(:,:,:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims4(1:4)
  integer              :: ierr


  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  !variables required to compute 2el integrals
  integer(8)         :: ijkl
  integer(8)         :: ioff,joff,koff,loff
  integer(8)         :: ish,jsh
  integer(8)         :: kshell,lshell
  integer            :: i, ni, li
  real(8)            :: ei, ci, xi, yi, zi
  integer            :: j, nj, lj
  real(8)            :: ej, cj, xj, yj, zj
  integer            :: k, nk, lk
  real(8)            :: ek, ck, xk, yk, zk
  integer            :: l, nl, ll
  real(8)            :: el, cl, xl, yl, zl
  real(8), parameter :: FIJKL = 1.0
  integer, parameter :: LMAX = 21 ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), save      :: gout(LMAX*LMAX*LMAX*LMAX) 
  real(8)            :: limit

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !set parameters for the fixed shell
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  el   =  gto(lsh)%exponent(1)
  xl   =  gto(lsh)%coord(1)
  yl   =  gto(lsh)%coord(2)
  zl   =  gto(lsh)%coord(3)
  cl   =  gto(lsh)%coefficient(1)
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  ek   =  gto(ksh)%exponent(1)
  xk   =  gto(ksh)%coord(1)
  yk   =  gto(ksh)%coord(2)
  zk   =  gto(ksh)%coord(3)
  ck   =  gto(ksh)%coefficient(1)
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  rank=talsh_tensor_rank(Dshell)
  if (rank.ne.4) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Dshell,rank,dims4)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Dshell,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,Dlk,dims4)

  !loop over all the shells for the other indices
  joff = 0
  do jsh = 1, nshells
    lj   =  gto(jsh)%orb_momentum
    ej   =  gto(jsh)%exponent(1)
    xj   =  gto(jsh)%coord(1)
    yj   =  gto(jsh)%coord(2)
    zj   =  gto(jsh)%coord(3)
    cj   =  gto(jsh)%coefficient(1)
    nj   =  nfunctions(lj,basis_angular)

    !neglect small elements
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
      call interest_eri('llll',FIJKL,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )

      !save all elements for this shell
      do k = 1, nk
        do l = 1, nl
          do j = 1, nj
            do i = 1, ni
              ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
              Dlk(i+ioff,j+joff,l,k) = dcmplx(gout(ijkl),0.D0)
            end do
          end do
        end do
      end do

      ioff = ioff + ni
    end do

    joff = joff + nj
  end do

end subroutine talsh_get_shell_comp


subroutine talsh_find_el_real(Diag, Filter, lsh, ksh, l, k, D_el)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none
  
  integer, intent(in)               :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Filter
  type(talsh_tens_t), intent(inout) :: Diag
  integer, intent(out)              :: l, k
  real(8), intent(out)              :: D_el

  real(8), pointer     :: Dm(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), mind(1:2)
  integer              :: ierr
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer(8)         :: koff, loff
  integer            :: kshell, lshell
  integer            :: nk, nl
  integer            :: ll, lk

  logical, allocatable :: Dmask(:,:)

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !get shell information
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Diag,body_p,R8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  ! set up Mask
  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=.false.

  allocate(Dmask(dims2(1),dims2(2)))
  Dmask=.false.
  Dmask(koff+1:koff+nk,loff+1:loff+nl)=.true.
  where (Dmask) lfilt = rfilt.gt.0

  mind=maxloc(abs(Dm), mask=lfilt)

  k=mind(1)
  l=mind(2)
  D_el=Dm(l,k)
  
  deallocate(lfilt)
  deallocate(Dmask)

end subroutine talsh_find_el_real

subroutine talsh_find_el_comp(Diag, Filter, lsh, ksh, l, k, D_el)

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global

  implicit none
  
  integer, intent(in)               :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Filter
  type(talsh_tens_t), intent(inout) :: Diag
  integer, intent(out)              :: l, k
  complex(8), intent(out)           :: D_el

  complex(8), pointer  :: Dm(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), mind(1:2)
  integer              :: ierr
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer(8)         :: koff, loff
  integer            :: kshell, lshell
  integer            :: nk, nl
  integer            :: ll, lk

  logical, allocatable :: Dmask(:,:)

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !get shell information
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Diag,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  ! set up Mask
  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=.false.

  allocate(Dmask(dims2(1),dims2(2)))
  Dmask=.false.
  Dmask(koff+1:koff+nk,loff+1:loff+nl)=.true.
  where (Dmask) lfilt = rfilt.gt.0

  if (any(lfilt)) then
    mind=maxloc(abs(Dm), mask=lfilt)
    k=mind(1)
    l=mind(2)
    D_el=Dm(l,k)
  else
    k=0
    l=0
    D_el=ZERO
  end if
  
  deallocate(lfilt)
  deallocate(Dmask)

end subroutine talsh_find_el_comp

subroutine talsh_square_el(h_tensor)

  use tensor_algebra
  use talsh

  implicit none
  
  type(talsh_tens_t), intent(inout)    :: h_tensor

  integer  ::  dims2(1:2)
  integer  ::  dims3(1:3)
  integer  ::  dims4(1:4)
  integer  ::  ierr
  integer  ::  rank, rrank
  integer  ::  i, j, k, l
  type(C_PTR)         :: body_p
  complex(8), pointer :: tens2(:, :)
  complex(8), pointer :: tens3(:, :, :)
  complex(8), pointer :: tens4(:, :, :, :)

  rank=talsh_tensor_rank(h_tensor)
        
  if (rank.eq.2) then
    ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
    if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'

    ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
    call c_f_pointer(body_p,tens2,dims2)
    
    do i = 1, dims2(1)
      do j = 1, dims2(2)
        tens2(i,j)=tens2(i,j)*tens2(i,j)
      end do
    end do

  else if (rank.eq.3) then
    ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
    if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'

    ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
    call c_f_pointer(body_p,tens3,dims3)

    do i = 1, dims3(1)
      do j = 1, dims3(2)
        do k = 1, dims3(3)
          tens3(i, j, k)=tens3(i, j, k)*tens3(i, j, k)
        end do
      end do
    end do

  else if (rank.eq.4) then
    ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
    if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'

    ierr=talsh_tensor_get_body_access(h_tensor,body_p,C8,0,DEV_HOST)
    call c_f_pointer(body_p,tens4,dims4)

    do i = 1, dims4(1)
      do j = 1, dims4(2)
        do k = 1, dims4(3)
          do l = 1, dims4(4)
            tens4(i, j, k, l)=tens4(i, j, k, l)*tens4(i, j, k, l)
          end do
        end do
      end do
    end do

  else
    stop 'error in print_tensor: only ranks 2,3 and 4 implemented'
  end if

end subroutine talsh_square_el

subroutine talsh_set_mask_done(Filter,l,k)

  use tensor_algebra
  use talsh
  
  implicit none
  
  integer, intent(in)               :: l, k
  type(talsh_tens_t), intent(inout) :: Filter

  type(C_PTR)          :: body_p
  integer              :: dims2(1:2)
  integer              :: ierr
  real(4), pointer     :: rfilt(:,:)
  integer              :: rank

  ! set up Mask
  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  
  rfilt(l,k)=-1.0

end subroutine talsh_set_mask_done

subroutine get_offset(l,l0)

  use exacorr_global

  implicit none

  integer, intent(in)  :: l
  integer, intent(out) :: l0
  
  integer :: i, lshells

  lshells=get_shell_index(l)
  
  l0=0
  do i = 1,lshells-1
    l0=l0+get_shell_size(i)
  end do

end subroutine get_offset

integer function get_shellsize(lsh)

  use exacorr_global

  implicit none

  integer, intent(in) :: lsh

  get_shellsize = get_shell_size(lsh)

end function get_shellsize

end module exacorr_cholesky
