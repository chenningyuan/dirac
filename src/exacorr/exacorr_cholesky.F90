module exacorr_cholesky

!routines to perform cholesky decomposition
!Johann Pototschnig, 04-2019

use, intrinsic:: ISO_C_BINDING

implicit none

complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

private

public cholesky_decompose
public cholesky_decompose_aomat
public cholesky_decompose_aomat_real

public find_shells_talsh
public get_diag_talsh
public get_shell_talsh
public find_el_talsh


interface cholesky_decompose
   module procedure cholesky_decompose_real_matrix
   module procedure cholesky_decompose_real_tensor
end interface

contains

subroutine cholesky_decompose_real_matrix (matrix, mdim, CD_out, threshold)

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
    
    if (m.ge.mdim) then
      stop 'error: cholesky not converged'
    end if

    !compute next element
    do i = 1,m
       CD(:,m+1)=CD(:,m+1)+CD(:,i)*CD(i_max,i)
    end do
    CD(:,m+1)=matrix(:,i_max)-CD(:,m+1)
    if (D_max < threshold) then
!      this should not happen: stop (alternative would be to have an error output)
       stop 'input matrix is not positive definite: Cholesky decomposition not possible'
    end if
    D_max_denom = ONE / dsqrt(D_max)
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

  end do

  print *, 'Cholesky decomposition done, used ',m,' elements'
  print *, 'Matrix compressed to ',real(m)/real(mdim)*100,'% of the original size'

  deallocate(Gkl)
  deallocate(Dm)

  !hand tensor over
  allocate(CD_out(mdim, m))
  CD_out=CD(:,1:m)
  deallocate(CD)
end subroutine cholesky_decompose_real_matrix

subroutine cholesky_decompose_real_tensor (tensor, ndim, CD_out, threshold)

  use exacorr_datatypes
  use exacorr_global

  implicit none

  !in: real 4-d tensor
  real(8), intent(in)                  :: threshold
  real(8), intent(in)                  :: tensor(ndim,ndim,ndim,ndim)

  !out: 3 component tensors
  real(8), intent(out), allocatable :: CD_out(:,:,:)

  !work
  real(8), allocatable :: Dm(:,:)
  integer(8)           :: m, i, j, nonzero, nkeep
  integer(8)           :: i_max(2)
  real(8)              :: D_max,D_max_denom
  real(8), allocatable :: Gkl(:,:)
  real(8), allocatable :: CD(:,:,:)
  integer              :: ndim, log_conv
  !work for SVD via lapack
  real(8), allocatable :: RES(:,:),ev(:),work(:)
  integer, allocatable :: iwork(:)
  integer              :: lwork,liwork,info

  !allocate tensors for cholesky decomposition
  allocate(CD(ndim,ndim,ndim*ndim))
  allocate(Dm(ndim,ndim))
  allocate(Gkl(ndim,ndim))

  ! allocate matrices for singular value decomposition
  allocate(RES(ndim,ndim))
  allocate(ev(ndim))
  lwork = 1 + 6*ndim + 2*ndim**2 ! this is the minimum, maybe should scale this up
  allocate(work(lwork))
  liwork = 3 + 5*ndim ! again the minimum, for now
  allocate(iwork(8*ndim))

  ! set diagonal matrix
  do j = 1, ndim
     do i = 1, ndim
        Dm(i,j) = tensor(i,j,i,j)
     end do
  end do

  !find largest value
  i_max=maxloc(abs(Dm))
  D_max=Dm(i_max(1), i_max(2))
  write (*,'(a,a,2i6,a,e12.2)') 'start iter',' largest element ind = ', i_max, ' value = ', abs(D_max)

  !loop until threshold is reached

  log_conv = ceiling(log10(abs(D_max)))
  CD=ZERO
  m=0
  nonzero = 0
  do while (abs(D_max)>threshold)

    !compute next element
    do i = 1,m
       CD(:,:,m+1)=CD(:,:,m+1)+CD(:,:,i)*CD(i_max(1),i_max(2),i)
    end do
    CD(:,:,m+1)=tensor(:,:,i_max(1),i_max(2))-CD(:,:,m+1)
    D_max_denom = ONE / sqrt(D_max)
    CD(:,:,m+1)=CD(:,:,m+1) * D_max_denom

    m=m+1
    !compress CD vector
    RES = CD(:,:,m)
    call dsyevd ('V','L',ndim,RES,ndim,ev,work,lwork,iwork,liwork,info)
    if (info /= 0) print*,'iteration ',m,' dsyevd info =',info
    nkeep = ndim
    do i = 1, ndim
       if (ev(i) > threshold) then
          nkeep = i
          exit
       end if
    end do
    nonzero = nonzero + (nkeep+1) * ndim

    !update diagonal
    Gkl(:,:)=CD(:,:,m)*CD(:,:,m)
    Dm=Dm-Gkl

    !find largest element
    i_max=maxloc(abs(Dm))
    D_max=Dm(i_max(1),i_max(2))

    if (ceiling(log10(abs(D_max))) /= log_conv) then
       write (*,'(a,i6,a,2i6,a,e12.2)') 'iter',m,' largest element ind = ', i_max, ' value = ', abs(D_max)
       log_conv = ceiling(log10(abs(D_max)))
    end if

    if (m.ge.ndim*ndim) then
      stop 'error: cholesky not converged'
    end if
  end do

  print *, 'Cholesky decomposition done, used ',m,' elements'
  print *, 'Tensor compressed to ',real(m)/real(ndim*ndim)*100,'% of the original size'
  print *, 'Additional compression possible to ', real(nonzero)/real(ndim)**4*100,'%'

  deallocate(Gkl)
  deallocate(Dm)

  !hand tensor over
  allocate(CD_out(ndim,ndim, m))
  CD_out=CD(:,:,1:m)
  deallocate(CD)
end subroutine cholesky_decompose_real_tensor

subroutine cholesky_decompose_aomat (CD_out, threshold)

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri
  
  implicit none

  !in: ao information
  real(8), intent(in)                  :: threshold
  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells

  !out: 3 component tensors
  complex(8), intent(out), allocatable :: CD_out(:,:,:)

  !quanities for doing the main work
  integer(8)              :: m
  integer(8)              :: i_max(1:2)
  complex(8)              :: D_max,D_max_denom, D_tot
  complex(8), allocatable :: Dm(:,:)
  complex(8), allocatable :: Gkl(:,:,:)
  complex(8), allocatable :: CD(:,:,:)
  logical, allocatable    :: M_do(:,:), M_ncomp(:,:)

  !variables required to compute 2el integrals
  integer            :: basis_angular
  integer(8)         :: ijkl
  integer(8)         :: ioff,joff,koff,loff
  integer(8)         :: ish,jsh,ksh,lsh
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

  ! variables required to treat all significant elements of the shell
  integer            :: max_ind(1:LMAX*LMAX,1:2)
  integer            :: max_m, q
  real(8)            :: deltaD = 1.0E-3

  !get all information about the Gauss type orbitals, allocate retrieve data
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)

  !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
  call interest_initialize(.false.)
  !set threshold for which integrals are include
  limit=threshold*threshold*1.0E-6
  !initialization for 2el integrals
  basis_angular = get_basis_angular()

  !allocate tensors for cholesky decompostion
  allocate(CD(nao, nao, nao*nao))
  allocate(Dm(nao,nao))
  allocate(M_do(nao,nao))
  M_do=.true.
  allocate(M_ncomp(nao,nao))

  call get_2el_diag(nao, gto, nshells, Dm)

  !find largest value
  i_max=maxloc(abs(Dm))
  D_max=Dm(i_max(1), i_max(2))
  write (*,'(a,a,2i6,a,e12.2)') 'start iter',' largest element ind = ', i_max, ' value = ', abs(D_max)
  
  !loop until threshold is reached
  CD=ZERO
  m=0
  do while (abs(D_max)>threshold)
    
    !find shell and save values for integral computation to variables
    loff = 0
    ! This search could be avoided by generalizing the routine that points from a function to the shell
    do lsh = 1, nshells
      ll   =  gto(lsh)%orb_momentum
      nl   =  nfunctions(ll,basis_angular)
      if ( i_max(2)>loff .and. i_max(2).le.loff+nl) then
        el   =  gto(lsh)%exponent(1)
        xl   =  gto(lsh)%coord(1)
        yl   =  gto(lsh)%coord(2)
        zl   =  gto(lsh)%coord(3)
        cl   =  gto(lsh)%coefficient(1)
        EXIT
      else
       loff = loff + nl
      end if
    end do

    koff = 0
    do ksh = 1, nshells
      lk   =  gto(ksh)%orb_momentum
      nk   =  nfunctions(lk,basis_angular)
      if ( i_max(1)>koff .and. i_max(1).le.koff+nk) then
        ek   =  gto(ksh)%exponent(1)
        ck   =  gto(ksh)%coefficient(1)
        xk   =  gto(ksh)%coord(1)
        yk   =  gto(ksh)%coord(2)
        zk   =  gto(ksh)%coord(3)
        EXIT
      else
       koff = koff + nk
      end if
    end do

    !find all large elements in the shell, allocate accordingly
    M_ncomp=.false.
    M_ncomp(koff+1:koff+nk,loff+1:loff+nl)=.true.
    max_ind(1,:)=i_max
    M_ncomp(max_ind(1,1),max_ind(1,2))=.false.
    
    do q =2,LMAX*LMAX
      if (count(M_ncomp)>0) then
        max_ind(q,:)=maxloc(abs(Dm), mask=M_ncomp)
        if (abs(Dm(max_ind(q,1),max_ind(q,2)))<abs(D_max)*deltaD) then
          max_m=q-1
          EXIT
        end if
        M_ncomp(max_ind(q,1),max_ind(q,2))=.false.
      else
        max_m=q-1
        EXIT
      end if
    end do
    
    allocate(Gkl(nao,nao,max_m))

    !compute nessecary two electron integrals
    M_ncomp=.true.

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
      if (any(abs(Dm(:,joff+1:joff+nj)*D_max)>limit)) then
        ioff = 0
        do ish = 1, nshells
          li   =  gto(ish)%orb_momentum
          ei   =  gto(ish)%exponent(1)
          xi   =  gto(ish)%coord(1)
          yi   =  gto(ish)%coord(2)
          zi   =  gto(ish)%coord(3)
          ci   =  gto(ish)%coefficient(1)
          ni   =  nfunctions(li,basis_angular)

          !neglect small elements
          if (any(abs(Dm(ioff+1:ioff+ni,joff+1:joff+nj)*D_max)>limit)) then
!           output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
            call interest_eri('llll',FIJKL,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )

            !save all significant elements for this shell
            do q = 1,max_m
              k = max_ind(q,1) - koff
              l = max_ind(q,2) - loff
              do j = 1, nj
                do i = 1, ni
                  ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
                  Gkl(i+ioff,j+joff,q) = dcmplx(gout(ijkl),0.D0)
                  M_ncomp(i+ioff,j+joff)=.false.
                end do
              end do
            end do

          end if
          ioff = ioff + ni
        end do
      end if
      joff = joff + nj
    end do

    !add next elements
    do q = 1,max_m
      !check if previous vectors already took care of this element
      if (abs(Dm(max_ind(q,1),max_ind(q,2)))>abs(D_max)*deltaD) then

        !compute new element 
        do i = 1,m
          CD(:,:,m+1)=CD(:,:,m+1)+CD(:,:,i)*CD(max_ind(q,1),max_ind(q,2),i)
        end do
        CD(:,:,m+1)=Gkl(:,:,q)-CD(:,:,m+1)
        !zero out elements that are not treated but should be small
        where (M_ncomp) CD(:,:,m+1)=ZERO
        D_max_denom = ONE / zsqrt(Dm(max_ind(q,1),max_ind(q,2)))
        CD(:,:,m+1)=CD(:,:,m+1) * D_max_denom

        ! next element
        m=m+1   

        !update diagonal
        Gkl(:,:,1)=CD(:,:,m)*CD(:,:,m)
        Dm=Dm-Gkl(:,:,1)

        !make mask for searching new elements
        M_do(max_ind(q,1),max_ind(q,2))=.false.

      end if
    end do

    !find largest element
    i_max=maxloc(abs(Dm))
    D_tot=Dm(i_max(1), i_max(2))

    !find largest non treated element
    i_max=maxloc(abs(Dm), mask=M_do)
    D_max=Dm(i_max(1), i_max(2))

    if (abs(D_tot).gt.abs(D_max)) then
      print *, 'Warning: error in cholesky will not be smaller than ', D_tot
    end if

    write (*,'(a,i6,a,2i6,a,e12.2)') 'cholesky length =',m,' largest element ind = ', i_max, ' value = ', abs(D_max)

    if (m.ge.nao*nao) then
      stop 'error: cholesky not converged'
    end if

    deallocate(Gkl)
  end do

  print *, 'Cholesky dcomposition done '

  ! clean up
  deallocate(Dm)
  deallocate(M_do)
  deallocate(M_ncomp)
  deallocate(gto)

  !hand tensor over
  allocate(CD_out(nao, nao, m))
  CD_out=CD(:,:,1:m)
  deallocate(CD)
end subroutine cholesky_decompose_aomat

subroutine cholesky_decompose_aomat_real(CD_out, threshold)

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri
  
  implicit none

  !in: ao information
  real(8), intent(in)                  :: threshold
  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells

  !out: 3 component tensors
  real(8), intent(out), allocatable :: CD_out(:,:,:)

  !quanities for doing the main work
  integer(8)           :: m
  integer(8)           :: i_max(1:2)
  real(8)              :: D_max,D_max_denom, D_tot
  real(8), allocatable :: Dm(:,:)
  real(8), allocatable :: Gkl(:,:,:)
  real(8), allocatable :: CD(:,:,:)
  logical, allocatable :: M_do(:,:), M_ncomp(:,:)

  !variables required to compute 2el integrals
  integer            :: basis_angular
  integer(8)         :: ijkl
  integer(8)         :: ioff,joff,koff,loff
  integer(8)         :: ish,jsh,ksh,lsh
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

  ! variables required to treat all significant elements of the shell
  integer            :: max_ind(1:LMAX*LMAX,1:2)
  integer            :: max_m, q
  real(8)            :: deltaD = 1.0E-3

  !get all information about the Gauss type orbitals, allocate retrieve data
  nao = get_nao()
  call get_gtos(1,nao,gto,nshells)

  !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
  call interest_initialize(.false.)
  !set threshold for which integrals are include
  limit=threshold*threshold*1.0E-6
  !initialization for 2el integrals
  basis_angular = get_basis_angular()

  !allocate tensors for cholesky decompostion
  allocate(CD(nao, nao, nao*nao))
  allocate(Dm(nao,nao))
  allocate(M_do(nao,nao))
  M_do=.true.
  allocate(M_ncomp(nao,nao))

  call get_2el_diag_real(nao, gto, nshells, Dm)

  !find largest value
  i_max=maxloc(abs(Dm))
  D_max=Dm(i_max(1), i_max(2))
  write (*,'(a,a,2i6,a,e12.2)') 'start iter',' largest element ind = ', i_max, ' value = ', abs(D_max)
  
  !loop until threshold is reached
  CD=0.0
  m=0
  do while (abs(D_max)>threshold)
    
    !find shell and save values for integral computation to variables
    loff = 0
    ! This search could be avoided by generalizing the routine that points from a function to the shell
    do lsh = 1, nshells
      ll   =  gto(lsh)%orb_momentum
      nl   =  nfunctions(ll,basis_angular)
      if ( i_max(2)>loff .and. i_max(2).le.loff+nl) then
        el   =  gto(lsh)%exponent(1)
        xl   =  gto(lsh)%coord(1)
        yl   =  gto(lsh)%coord(2)
        zl   =  gto(lsh)%coord(3)
        cl   =  gto(lsh)%coefficient(1)
        EXIT
      else
       loff = loff + nl
      end if
    end do

    koff = 0
    do ksh = 1, nshells
      lk   =  gto(ksh)%orb_momentum
      nk   =  nfunctions(lk,basis_angular)
      if ( i_max(1)>koff .and. i_max(1).le.koff+nk) then
        ek   =  gto(ksh)%exponent(1)
        ck   =  gto(ksh)%coefficient(1)
        xk   =  gto(ksh)%coord(1)
        yk   =  gto(ksh)%coord(2)
        zk   =  gto(ksh)%coord(3)
        EXIT
      else
       koff = koff + nk
      end if
    end do

    !find all large elements in the shell, allocate accordingly
    M_ncomp=.false.
    M_ncomp(koff+1:koff+nk,loff+1:loff+nl)=.true.
    max_ind(1,:)=i_max
    M_ncomp(max_ind(1,1),max_ind(1,2))=.false.
    
    do q =2,LMAX*LMAX
      if (count(M_ncomp)>0) then
        max_ind(q,:)=maxloc(abs(Dm), mask=M_ncomp)
        if (abs(Dm(max_ind(q,1),max_ind(q,2)))<abs(D_max)*deltaD) then
          max_m=q-1
          EXIT
        end if
        M_ncomp(max_ind(q,1),max_ind(q,2))=.false.
      else
        max_m=q-1
        EXIT
      end if
    end do
    
    allocate(Gkl(nao,nao,max_m))

    !compute nessecary two electron integrals
    M_ncomp=.true.

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
      if (any(abs(Dm(:,joff+1:joff+nj)*D_max)>limit)) then
        ioff = 0
        do ish = 1, nshells
          li   =  gto(ish)%orb_momentum
          ei   =  gto(ish)%exponent(1)
          xi   =  gto(ish)%coord(1)
          yi   =  gto(ish)%coord(2)
          zi   =  gto(ish)%coord(3)
          ci   =  gto(ish)%coefficient(1)
          ni   =  nfunctions(li,basis_angular)

          !neglect small elements
          if (any(abs(Dm(ioff+1:ioff+ni,joff+1:joff+nj)*D_max)>limit)) then
!           output order of eri is (5,c,d,5,a,b), so input k,l first and then i,j to get the order that we want
            call interest_eri('llll',FIJKL,gout,&
                        lk,ek,xk,yk,zk,ck,&
                        ll,el,xl,yl,zl,cl,&
                        li,ei,xi,yi,zi,ci,&
                        lj,ej,xj,yj,zj,cj )

            !save all significant elements for this shell
            do q = 1,max_m
              k = max_ind(q,1) - koff
              l = max_ind(q,2) - loff
              do j = 1, nj
                do i = 1, ni
                  ijkl = (l-1)*nk*nj*ni+(k-1)*nj*ni+(j-1)*ni+i
                  Gkl(i+ioff,j+joff,q) = gout(ijkl)
                  M_ncomp(i+ioff,j+joff)=.false.
                end do
              end do
            end do

          end if
          ioff = ioff + ni
        end do
      end if
      joff = joff + nj
    end do

    !add next elements
    do q = 1,max_m
      !check if previous vectors already took care of this element
      if (abs(Dm(max_ind(q,1),max_ind(q,2)))>abs(D_max)*deltaD) then

        !compute new element 
        do i = 1,m
          CD(:,:,m+1)=CD(:,:,m+1)+CD(:,:,i)*CD(max_ind(q,1),max_ind(q,2),i)
        end do
        CD(:,:,m+1)=Gkl(:,:,q)-CD(:,:,m+1)
        !zero out elements that are not treated but should be small
        where (M_ncomp) CD(:,:,m+1)=0.0
        D_max_denom = ONE / sqrt(Dm(max_ind(q,1),max_ind(q,2)))
        CD(:,:,m+1)=CD(:,:,m+1) * D_max_denom

        ! next element
        m=m+1   

        !update diagonal
        Gkl(:,:,1)=CD(:,:,m)*CD(:,:,m)
        Dm=Dm-Gkl(:,:,1)

        !make mask for searching new elements
        M_do(max_ind(q,1),max_ind(q,2))=.false.

      end if
    end do

    !find largest element
    i_max=maxloc(abs(Dm))
    D_tot=Dm(i_max(1), i_max(2))

    !find largest non treated element
    i_max=maxloc(abs(Dm), mask=M_do)
    D_max=Dm(i_max(1), i_max(2))

    if (abs(D_tot).gt.abs(D_max)) then
      print *, 'Warning: error in cholesky will not be smaller than ', D_tot
    end if

    write (*,'(a,i6,a,2i6,a,e12.2)') 'cholesky length =',m,' largest element ind = ', i_max, ' value = ', abs(D_max)

    if (m.ge.nao*nao) then
      stop 'error: cholesky not converged'
    end if

    deallocate(Gkl)
  end do

  print *, 'Cholesky dcomposition done '

  ! clean up
  deallocate(Dm)
  deallocate(M_do)
  deallocate(M_ncomp)
  deallocate(gto)

  !hand tensor over
  allocate(CD_out(nao, nao, m))
  CD_out=CD(:,:,1:m)
  deallocate(CD)
end subroutine cholesky_decompose_aomat_real

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

subroutine get_2el_diag_real(nao, gto, nshells, Dm)

  use exacorr_datatypes
  use exacorr_global
  use module_interest_eri

  implicit none
  
  !returns two electron integrals with i=k and j=l
  integer, intent(in) :: nao, nshells
  type(basis_func_info_t), intent(in) :: gto(:)  ! arrays with basis function information

  real(8), intent(inout) :: Dm(:,:)

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
            Dm(i+ioff,j+joff) = gout(ijij)
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine get_2el_diag_real

subroutine find_shells_talsh(Diag, Filter, lshell, kshell, lsize, ksize, D_max)

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
  ierr=talsh_tensor_get_body_access(Diag,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  
  !find largest element
  i_max=maxloc(abs(Dm))
  D_tot=Dm(i_max(1), i_max(2))

  !find largest non treated element
  lfilt=rfilt.lt.0
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

end subroutine find_shells_talsh

subroutine get_diag_talsh(Diag)

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
            Dm(i+ioff,j+joff) = gout(ijij)
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine get_diag_talsh

subroutine get_shell_talsh(Dshell, lsh, ksh)

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
              Dlk(i+ioff,j+joff,l,k) = gout(ijkl)
            end do
          end do
        end do
      end do

      ioff = ioff + ni
    end do

    joff = joff + nj
  end do

end subroutine get_shell_talsh


subroutine find_el_talsh(Diag, Filter, lsh, ksh, l, k, D_el)

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
  ierr=talsh_tensor_get_body_access(Diag,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,Dm,dims2)

  ! set up Mask
  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,C8,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=.false.

  allocate(Dmask(dims2(1),dims2(2)))
  Dmask=.false.
  Dmask(koff+1:koff+nk,loff+1:loff+nl)=.true.
  where (Dmask) lfilt = rfilt.lt.0

  mind=maxloc(abs(Dm), mask=lfilt)

  k=mind(1)
  l=mind(2)
  D_el=Dm(l,k)
  
  deallocate(lfilt)
  deallocate(Dmask)

end subroutine find_el_talsh

end module exacorr_cholesky
