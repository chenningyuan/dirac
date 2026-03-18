module eti_CD_routines

!routines to perform cholesky decomposition

!       written by Johann Pototschnig, Summer 2019
!       some refactoring, Lucas Visscher, april 2020


use tensor_algebra
use exacorr_tensor_interface
use exacorr_scaling_parameters
use eti_common_routines
use exacorr_utils
!use, intrinsic:: ISO_C_BINDING

implicit none

private

!talsh routines:
public decompose_cholesky_eti

contains

    subroutine decompose_cholesky_eti (CD, t_cholesky, print_level)

!       Cholesky decomposition based on treating the shells separately using tensor operations
!       possible additions: 
!                          SVD (How to combine with tensor library)
!                          Additional screening in the computation of integrals

!       written by Johann Pototschnig, Summer 2019
!       refactoring to reduce amount of memory copies, LV 2023

        use exacorr_global

        implicit none

        real(8), intent(in)                  :: t_cholesky ! cholesky threshold
        type(tensor_t), intent(inout)    :: CD
        integer, intent(in)                  :: print_level

        integer, parameter :: CD_blocksize = 100
        integer(C_INT)     :: ierr
        integer            :: nao, lshell, kshell, lsize, ksize
        integer            :: l, k, l0, k0
        real(8)            :: D_max, D_el
        type(tensor_t) :: Diag, Res, Temp
        type(tensor_t) :: Dlk
        type(tensor_t) :: CD_row, CD_column
        type(tensor_t),allocatable :: CD_blocks(:)
        integer(C_INT)     :: dims2(1:2)
        integer(C_INT)     :: dims3(1:3), off3(1:3)
        integer(C_INT)     :: dims4(1:4), off4(1:4)
        integer(C_INT)     :: rank
        integer(C_INT)     :: i, m, mblock, iblock
        type(tensor_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)
        real(8), pointer     :: CDR(:,:,:),DiagR(:,:)
        complex(8),pointer   :: CDC(:,:,:)
        type(C_PTR)          :: body_p
        logical, allocatable :: Mask(:,:)

        ! accuracy for tensor printing
        real(8), parameter :: ACC=1.0E-15

        !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
        call interest_initialize(.false.)

        if (print_level.gt.2) write(*,*) "========= entered decompose_cholesky_eti ========= "

        !Start by setting the maximum dimensions of tensors and initializing the offsets
        nao = get_nao()
        dims2       = nao
        dims3(1:2)  = nao
        dims3(3)    = CD_blocksize
        dims4       = nao
        off3        = 0
        off4        = 0
        one_dims(1) = 1

        ! allocate the stupid one tensor that is needed for multiplying with a scalar in TALSH
        ierr=construct_tensor_t(one_tensor,R8,one_dims(1:0),init_val=ONE)

        ! tensor to hold the diagonal of the integral matrix and also a pointer for direct use thereof in Fortran
        ierr=construct_tensor_t(Diag,R8,dims2,init_val=ZERO)
        ierr=direct_data_pointer(Diag,body_p)
        call c_f_pointer(body_p,DiagR,dims2)

        ! the following allocates an array to hold Cholesky blocks, but no tensor memory is allocated yet
        allocate(CD_blocks(max(1,nao*nao/CD_blocksize)))

        ! tensors to hold one column or one row of a block of the Cholesky tensor
        ierr=construct_tensor_t(Res,R8,dims2,init_val=ZERO)
        ierr=construct_tensor_t(Temp,R8,dims2,init_val=ZERO)
        ierr=construct_tensor_t(CD_row,R8,(/ CD_blocksize /),init_val=ZERO)
        ierr=construct_tensor_t(CD_column,R8,(/ nao, nao /),init_val=ZERO)

        ! a mask array to mask out columns that are already considered
        allocate (Mask(nao,nao))
        Mask = .true.

        call get_diag(DiagR)
        if (print_level.gt.2) write(*,*) "--- diagonal computed "

        !find the shell pair that has the largest diagonal element
        call find_shells(DiagR, Mask, lshell, kshell, D_max)
        if (print_level.gt.6) write(*,*) "--- shell pair found: ", lshell, kshell, D_max

        m=0

        !outermost loop over shell pairs that runs until all significant shell pairs have been added
        do while (abs(D_max)>t_cholesky)
          lsize=get_shell_size(lshell)
          ksize=get_shell_size(kshell)
          if (print_level.gt.2) then 
            write (*,'(a,i6,a,a,i6,i6,a,e12.2)') 'cholesky length =',m,' largest element:', &
            ' shell = ', lshell, kshell, ' value = ', abs(D_max)
          end if

          dims4(3)=lsize
          dims4(4)=ksize

          ! allocate the tensor that will hold all generated integrals for this shell pair
          ierr=construct_tensor_t(Dlk,R8,dims4,init_val=ZERO)

          !get the integral block for the current shell pair
          call get_integral_block(Dlk, lshell, kshell)
          if (print_level.gt.15) call print_tensor(Dlk, ACC, "Dlk")

          !find the initially largest diagonal element for this shell pair
          call eti_find_el(DiagR, Mask, lshell, kshell, l, k, D_el)
          if (print_level.gt.6) write(*,*) '--- largest value found: ', l, k, D_el
          call get_shell_offset(l,l0)
          call get_shell_offset(k,k0)

          !innermost loop over significant index pairs within the shell pair
          do while(abs(D_el)>abs(D_max/1000))

            off3(1)=l-1
            off3(2)=k-1
            off3(3)=0
            !sum up previous CD vectors for the current index pair
            ierr=initialize_tensor_t(Res, ZERO)
            mblock = m / CD_blocksize
            do iblock = 1, mblock
               !for the filled blocks we can use full contractions, but we need to extract the current row first
               ierr=reshape_static_tensor_t(CD_row,(/ 1,1,CD_blocksize /))
               ierr=slice_tensor_t(CD_row,CD_blocks(iblock),off3)
               ierr=reshape_static_tensor_t(CD_row,(/ CD_blocksize /))
               ierr=t_einsum("L(p,q)+=I(p,q,u)*R(u)",Res,CD_blocks(iblock),CD_row)
            enddo
            do i = 1, mod(m,CD_blocksize)
               !for the last incomplete block we slice out and multiply each column separately
               ierr=reshape_static_tensor_t(CD_column,(/ nao, nao, 1 /))
               ierr=slice_tensor_t(CD_column,CD_blocks(mblock+1),(/0,0,i-1/))
               ierr=reshape_static_tensor_t(CD_column,(/ nao, nao /))
               ! the row is now a single number that we extract via the Fortran pointer array
               ierr=t_einsum("L(p,q)+=I(p,q)",Res,CD_column,scale=ONE*CDR(l,k,i))
            end do

            !get Dlk for the current l and k
            ierr=initialize_tensor_t(Temp, ZERO)
            dims4(3)=1
            dims4(4)=1
            ierr=reshape_static_tensor_t(Temp,dims4)
            off4=0
            off4(3)=l-l0-1
            off4(4)=k-k0-1
            ierr=slice_tensor_t(Temp,Dlk,off4)
            ierr=reshape_static_tensor_t(Temp,dims2)

            !sum and denominate
            ierr=t_einsum("T(p,q)+=V(p,q)",Temp,Res,scale=MINUS_ONE)
            ierr=initialize_tensor_t(Res, ZERO)
            ierr=t_einsum("W(p,q)=V(p,q)",Res,Temp,scale=ONE/sqrt(D_el))
            if (print_level.gt.15) call print_tensor(Res, ACC, "CD entry")

            !insert the newly obtained Cholesky vector at the end of the last cholesky block
            dims3(1:2)=nao
            dims3(3)=1
            ierr=reshape_static_tensor_t(Res,dims3)
            off3(1:2)=0
            off3(3)=mod(m,CD_blocksize)
            if (off3(3) == 0) then
               ! first vector of a new block, allocate this new block
               dims3(1:2)=nao
               dims3(3)=CD_blocksize
               ierr=construct_tensor_t(CD_blocks(mblock+1),R8,dims3,init_val=ZERO)
               ! the fortran array pointer should refer to this final, incomplete, block
               ierr=direct_data_pointer(CD_blocks(mblock+1),body_p)
               call c_f_pointer(body_p,CDR,dims3)
            end if
            ierr=insert_static_tensor_t(CD_blocks(mblock+1),Res,off3)

            !update diagonal
            ierr=reshape_static_tensor_t(Res,dims2)
            call eti_square_el(Res)
            ierr=t_einsum("W(p,q)+=V(p,q)",Diag,Res,scale=MINUS_ONE)
            if (print_level.gt.15) call print_tensor(Diag, ACC, "Diag")

            !update mask
            Mask(l,k) = .false.

            m=m+1

            if (m>nao*nao) write(*,*) 'Cholesky decomposition did not work'

            !find next element in this shell
            call eti_find_el(DiagR, Mask, lshell, kshell, l, k, D_el)
            if (print_level.gt.6) write(*,*)  '--- largest value found: ', l, k, D_el

          end do

          !find next shell combination that contains the largest element
          call find_shells(DiagR, Mask, lshell, kshell, D_max)

          ierr=destroy_tensor_t(Dlk)

        end do

        if (print_level.gt.-1) write(*,*) 'Cholesky decomposition done, used ',m,' elements'

        ! allocate the complex output tensor
        dims3(1:2) = nao
        dims3(3)   = m
        ierr=construct_tensor_t(CD,C8,dims3)
        ierr=direct_data_pointer(CD,body_p)
        call c_f_pointer(body_p,CDC,dims3)

        ! copy the fully filled blocks into the complex output tensor and then delete them
        dims3(1:2) = nao
        dims3(3)   = CD_blocksize
        i = 0
        do iblock = 1, mblock
           ierr=direct_data_pointer(CD_blocks(iblock),body_p)
           call c_f_pointer(body_p,CDR,dims3)
           CDC(:,:,i+1:i+CD_blocksize) = dcmplx(CDR,0.D0)
           i = i + CD_blocksize
           ierr=destroy_tensor_t(CD_blocks(iblock))
        end do
        ! copy the final, partially filled, block and free also this memory
        if (m > i) then
           ierr=direct_data_pointer(CD_blocks(mblock+1),body_p)
           call c_f_pointer(body_p,CDR,dims3)
           CDC(:,:,i+1:m) = dcmplx(CDR(:,:,1:m-i),0.D0)
           ierr=destroy_tensor_t(CD_blocks(mblock+1))
        end if

        !clean up
        ierr=destroy_tensor_t(CD_row)
        ierr=destroy_tensor_t(CD_column)
        ierr=destroy_tensor_t(Diag)
        ierr=destroy_tensor_t(Res)
        ierr=destroy_tensor_t(Temp)
        ierr=destroy_tensor_t(one_tensor)

        if (print_level.gt.2) write(*,*)  "========= leaving decompose_cholesky_eti ========= "

    end subroutine decompose_cholesky_eti

subroutine get_diag(DiagR)

! routine to compute diagonal 2D tensor

  use exacorr_ao
  use exacorr_global
  use exacorr_eri

  implicit none

  real(8), pointer     :: DiagR(:,:)
  integer              :: ierr

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer :: lq, lp, nq, np, q, p
  integer :: i, j, ijij
  integer :: ioff, joff
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0
  
  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  ! compute diagonal elements
  joff = 0
  do q = 1, nshells
     lq   =  gto(q)%orb_momentum
     xq   =  gto(q)%coord(1)
     yq   =  gto(q)%coord(2)
     zq   =  gto(q)%coord(3)
     nq   =  nfunctions(lq,basis_angular)

     ioff = 0
     do p = 1, nshells
        lp   =  gto(p)%orb_momentum
        xp   =  gto(p)%coord(1)
        yp   =  gto(p)%coord(2)
        zp   =  gto(p)%coord(3)
        np   =  nfunctions(lp,basis_angular)

        !get eri values
        call compute_eri('llll',FIJKL,gout(1:np*nq*np*nq),&
                         lp,gto(p)%exponent,xp,yp,zp,gto(p)%coefficient,&
                         lq,gto(q)%exponent,xq,yq,zq,gto(q)%coefficient,&
                         lp,gto(p)%exponent,xp,yp,zp,gto(p)%coefficient,&
                         lq,gto(q)%exponent,xq,yq,zq,gto(q)%coefficient )

        do j = 1, nq
          do i = 1, np
            ijij = (j-1)*np*nq*np+(i-1)*nq*np+(j-1)*np+i
            DiagR(i+ioff,j+joff) = gout(ijij)
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine get_diag

subroutine find_shells(DiagR, Mask, lshell, kshell, D_max)

! routine to find shell pair with the largest element

  use exacorr_global, only: get_shell_index

  implicit none

  real(8), pointer, intent(inout)   :: DiagR(:,:)
  integer, intent(out)              :: lshell, kshell
  real(8), intent(out)              :: D_max
  logical, intent(in)               :: Mask(:,:)


  integer              :: i_max(1:2)
  real(8)              :: D_tot

  !find largest element
  i_max=maxloc(abs(DiagR))
  D_tot=DiagR(i_max(1), i_max(2))

  !find largest non treated element
  i_max=maxloc(abs(DiagR), mask=Mask)
  D_max=DiagR(i_max(1), i_max(2))

  if (abs(D_tot).gt.abs(D_max)) then
    print *, 'Warning: error in cholesky will not be smaller than ', D_tot
  end if

  lshell=get_shell_index(i_max(1))
  kshell=get_shell_index(i_max(2))

end subroutine find_shells

subroutine get_integral_block(Dshell, lsh, ksh)

! routine to compute all ao integrals for a shell pair

  use tensor_algebra
  use exacorr_tensor_interface
  use exacorr_ao
  use exacorr_global
  use exacorr_eri

  implicit none 

  integer, intent(in)            :: lsh, ksh
  type(tensor_t), intent(inout) :: Dshell

  real(8), pointer     :: DlkR(:,:,:,:)
  complex(8), pointer  :: DlkC(:,:,:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims4(1:4)
  integer              :: ierr


  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular
  integer                              :: nData,DataKind(1)

  !variables required to compute 2el integrals
  integer(8)         :: ijkl
  integer(8)         :: ioff,joff,koff,loff
  integer(8)         :: ish,jsh
  integer(8)         :: kshell,lshell
  integer            :: i, ni, li
  real(8)            :: xi, yi, zi
  integer            :: j, nj, lj
  real(8)            :: xj, yj, zj
  integer            :: k, nk, lk
  real(8)            :: xk, yk, zk
  integer            :: l, nl, ll
  real(8)            :: xl, yl, zl
  real(8), parameter :: FIJKL = 1.0
  integer, parameter :: LMAX = 21 ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), save      :: gout(LMAX*LMAX*LMAX*LMAX) 
  real(8)            :: limit

  ! get basis information
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)
  basis_angular = get_basis_angular()

  !set parameters for the fixed shells
  loff=0
  do lshell = 1, lsh-1
    ll   =  gto(lshell)%orb_momentum
    nl   =  nfunctions(ll,basis_angular)
    loff = loff + nl
  end do
  ll   =  gto(lsh)%orb_momentum
  xl   =  gto(lsh)%coord(1)
  yl   =  gto(lsh)%coord(2)
  zl   =  gto(lsh)%coord(3)
  nl   =  nfunctions(ll,basis_angular)

  koff=0
  do kshell = 1, ksh-1
    lk   =  gto(kshell)%orb_momentum
    nk   =  nfunctions(lk,basis_angular)
    koff = koff + nk
  end do
  lk   =  gto(ksh)%orb_momentum
  xk   =  gto(ksh)%coord(1)
  yk   =  gto(ksh)%coord(2)
  zk   =  gto(ksh)%coord(3)
  nk   =  nfunctions(lk,basis_angular)

  ! get acces to the tensor
  ierr = get_nmodes(Dshell,rank)
  if (rank.ne.4) stop 'error: wrong rank in get_shell_'
  ierr = get_extents_static(Dshell,dims4)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_h'
  ierr=get_data_kind_tensor_t(Dshell,DataKind(1))
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=direct_data_pointer(Dshell,body_p) ! ,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DlkR,dims4)
  case (C8)
     call c_f_pointer(body_p,DlkC,dims4)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

  !loop over all the shells for the other indices
  joff = 0
  do jsh = 1, nshells
    lj   =  gto(jsh)%orb_momentum
    xj   =  gto(jsh)%coord(1)
    yj   =  gto(jsh)%coord(2)
    zj   =  gto(jsh)%coord(3)
    nj   =  nfunctions(lj,basis_angular)

    ioff = 0
    do ish = 1, nshells
      li   =  gto(ish)%orb_momentum
      xi   =  gto(ish)%coord(1)
      yi   =  gto(ish)%coord(2)
      zi   =  gto(ish)%coord(3)
      ni   =  nfunctions(li,basis_angular)

!     output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
      call compute_eri('llll',FIJKL,gout(1:ni*nj*nk*nl),&
                       lk,gto(ksh)%exponent,xk,yk,zk,gto(ksh)%coefficient,&
                       ll,gto(lsh)%exponent,xl,yl,zl,gto(lsh)%coefficient,&
                       li,gto(ish)%exponent,xi,yi,zi,gto(ish)%coefficient,&
                       lj,gto(jsh)%exponent,xj,yj,zj,gto(jsh)%coefficient )

      !save all elements for this shell pair
      do k = 1, nk
        do l = 1, nl
          do j = 1, nj
            do i = 1, ni
              ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
              if (DataKind(1)==R8) then
                 DlkR(i+ioff,j+joff,l,k) = gout(ijkl)
              else
                 DlkC(i+ioff,j+joff,l,k) = dcmplx(gout(ijkl),0.D0)
              end if
            end do
          end do
        end do
      end do

      ioff = ioff + ni
    end do

    joff = joff + nj
  end do

end subroutine get_integral_block

subroutine eti_find_el(DiagR, Mask, lsh, ksh, l, k, D_el)

! routine to locate largest element in a shell pair

  use exacorr_ao
  use exacorr_global

  implicit none
  
  integer, intent(in)               :: lsh, ksh
  logical, intent(in)               :: Mask(:,:)
  real(8), pointer, intent(inout)   :: DiagR(:,:)
  integer, intent(out)              :: l, k
  real(8), intent(out)              :: D_el

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer(8)         :: koff, loff
  integer            :: kshell, lshell
  integer            :: nk, nl
  integer            :: ll, lk
  integer            :: mind(2)

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

  ! set up DMask for indices inside this shell pair
  allocate(Dmask(nao,nao))
  Dmask=.false.
  Dmask(koff+1:koff+nk,loff+1:loff+nl)=.true.
  Dmask = DMask .and. Mask

  if (any(Dmask)) then
     mind=maxloc(abs(DiagR), mask=Dmask)
     k=mind(1)
     l=mind(2)
     D_el=DiagR(l,k)
  else
    k=0
    l=0
    D_el=0.D0
  end if
  
  deallocate(Dmask)

end subroutine eti_find_el

subroutine eti_square_el(h_tensor)
! squaring talsh tensor
  implicit none
  
  type(tensor_t), intent(inout)    :: h_tensor

  integer, allocatable  ::  dims(:)
  integer  ::  ierr
  integer  ::  rank, rrank
  integer  ::  i, j, k, l
  type(C_PTR)         :: body_p
  real(8), pointer :: tens2(:, :)
  real(8), pointer :: tens3(:, :, :)
  real(8), pointer :: tens4(:, :, :, :)

  ierr = get_nmodes(h_tensor,rank)
  allocate(dims(rank))
  ierr = get_extents(h_tensor,dims)
  ierr = direct_data_pointer(h_tensor,body_p)
        
  if (rank.eq.2) then
! this makes no sense ?
!   ierr = talsh_tensor_dimensions(h_tensor,rrank,dims2)
!   if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'
    call c_f_pointer(body_p,tens2,dims)
    
    do i = 1, dims(1)
      do j = 1, dims(2)
        tens2(i,j)=tens2(i,j)*tens2(i,j)
      end do
    end do

  else if (rank.eq.3) then
! this makes no sense ?
!   ierr = talsh_tensor_dimensions(h_tensor,rrank,dims3)
!   if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'
    call c_f_pointer(body_p,tens3,dims)

    do i = 1, dims(1)
      do j = 1, dims(2)
        do k = 1, dims(3)
          tens3(i, j, k)=tens3(i, j, k)*tens3(i, j, k)
        end do
      end do
    end do

  else if (rank.eq.4) then
! this makes no sense ?
!   ierr = talsh_tensor_dimensions(h_tensor,rrank,dims4)
!   if (ierr.ne.0 .or. rrank.ne.rank) stop 'error in square_el: wrong rank'
    call c_f_pointer(body_p,tens4,dims)

    do i = 1, dims(1)
      do j = 1, dims(2)
        do k = 1, dims(3)
          do l = 1, dims(4)
            tens4(i, j, k, l)=tens4(i, j, k, l)*tens4(i, j, k, l)
          end do
        end do
      end do
    end do

  else
    stop 'error in print_tensor: only ranks 2,3 and 4 implemented'
  end if
  deallocate(dims)

end subroutine eti_square_el

end module eti_CD_routines
