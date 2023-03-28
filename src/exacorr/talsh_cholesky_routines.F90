module talsh_CD_routines

!routines to perform cholesky decomposition

!       written by Johann Pototschnig, Summer 2019
!       some refactoring, Lucas Visscher, april 2020


use tensor_algebra
use talsh
use talsh_common_routines
use exacorr_utils
!use, intrinsic:: ISO_C_BINDING

implicit none

complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)

private

!talsh routines:
public decompose_cholesky_talsh_vec
public decompose_cholesky_talsh_all

contains

subroutine decompose_cholesky_talsh_vec (CD, m, t_cholesky, nao, print_level)

!       Cholesky decomposition based on treating the shells seperately using talsh tensor operations
!       possible additions: 
!                          SVD (How to combine with tensor library)
!                          Additional screening in the computation of integrals

!       written by Johann Pototschnig, Summer 2019

        use exacorr_global

        implicit none

        real(8), intent(in)                            :: t_cholesky ! cholesky threshold
        type(talsh_tens_t), allocatable, intent(inout) :: CD(:)
        integer, intent(in)                            :: nao, print_level
        integer(C_INT), intent(out)                    :: m

        integer(C_INT)     :: ierr
        integer            :: lshell, kshell, lsize, ksize
        integer            :: l, k, l0, k0, i
        complex(8)         :: D_max, D_el, factor
        type(talsh_tens_t) :: Mask, Diag, Res, Temp
        type(talsh_tens_t) :: Dlk !, CD_v
        integer(C_INT)     :: dims2(1:2),dims4(1:4), hdims(1:2)
        integer(C_INT)     :: off2(1:2),off4(1:4)
        integer(C_INT)     :: rank
        type(talsh_tens_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        ! accuracy for tensor printing
        real(8), parameter :: ACC=1.0E-9

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

        if (print_level.gt.2) write(*,*) "========= entered decompose_cholesky_talsh ========= "

        hdims=1
        dims2=nao
        dims4=nao
        off2=0
        off4=0

        !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
        call interest_initialize(.false.)

        !alllocate cholesky tensor
        allocate(CD(nao*nao))

        ierr=talsh_tensor_construct(Mask,R4,dims2,init_val=ONE)
        ierr=talsh_tensor_construct(Diag,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Res,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Temp,C8,dims2,init_val=ZERO)
        !ierr=talsh_tensor_construct(CD_v,C8,hdims,init_val=ZERO)

        call talsh_get_diag(Diag)
        if (print_level.gt.2) write(*,*) "--- diagonal computed "

        call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)
        if (print_level.gt.6) write(*,*) "--- shell found: ", lshell, kshell, lsize, ksize, D_max

        m=0

        do while (abs(D_max)>t_cholesky)
          lsize=get_shell_size(lshell)
          ksize=get_shell_size(kshell)
          if (print_level.gt.2) then
            write (*,'(a,i6,a,a,i6,i6,a,e12.2)') ' cholesky length =',m,' largest element:', &
            ' shell = ', lshell, kshell, ' value = ', abs(D_max)
          end if

          dims4(3)=lsize
          dims4(4)=ksize

          ierr=talsh_tensor_construct(Dlk,C8,dims4,init_val=ZERO)

          !get the values for the current shell
          call talsh_get_shell(Dlk, lshell, kshell)
          if (print_level.gt.15) call print_tensor(Dlk, ACC, "Dlk")

          call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
          if (print_level.gt.6) write(*,*) '--- largest value found: ', l, k, D_el

          do while(abs(D_el)>abs(D_max/1000))
            call get_shell_offset(l,l0)
            call get_shell_offset(k,k0)

            !create current elment of the cholesky vector
            ierr=talsh_tensor_construct(CD(m+1),C8,dims2,init_val=ZERO)
            ierr=talsh_tensor_init(Res, ZERO)

            !sum up previous CD vectors for the indices
            off2(1)=l
            off2(2)=k
            do i = 1, m
                !ierr=talsh_tensor_init(CD_v, ZERO)
                !ierr=talsh_tensor_slice(CD_v,CD(i),off2,0,DEV_HOST)
                !ierr=talsh_tensor_reshape(CD_v,hdims(1:0))
                !if (print_level.gt.6) call print_tensor(CD_v, ACC, "CD_v")
                !ierr=talsh_tensor_contract("L(p,q)+=X(p,q)*Y()",Res,CD(i),CD_v)
                !ierr=talsh_tensor_reshape(CD_v,hdims)
                factor=get_element(CD(i),off2)
                if (print_level.gt.10) write(*,*) 'factor=',factor
                ierr=talsh_tensor_add("L(p,q)+=X(p,q)",Res,CD(i),scale=factor)
                if (print_level.gt.15) call print_tensor(Res, ACC, "Res")
            end do

            !get Dlk for the current l and k
            ierr=talsh_tensor_init(Temp, ZERO)
            dims4(3)=1
            dims4(4)=1
            ierr=talsh_tensor_reshape(Temp,dims4)
            off4=0
            off4(3)=l-l0-1
            off4(4)=k-k0-1
            ierr=talsh_tensor_slice(Temp,Dlk,off4,0,DEV_HOST)
            ierr=talsh_tensor_reshape(Temp,dims2)

            !sum and denominate
            ierr=talsh_tensor_contract("T(p,q)+=V(p,q)",Temp,Res,one_tensor,scale=MINUS_ONE)
            !factor=ONE / zsqrt(D_el)
            factor=ONE / sqrt(D_el)
            ierr=talsh_tensor_add("W(p,q)+=V(p,q)",CD(m+1),Temp,scale=factor)
            
            !update diagonal
            ierr=talsh_tensor_init(Res, ZERO)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Res,CD(m+1),one_tensor)
            if (print_level.gt.15) call print_tensor(Res, ACC, "CD entry")
            call talsh_square_el(Res)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Diag,Res,one_tensor,scale=MINUS_ONE)
            if (print_level.gt.15) call print_tensor(Diag, ACC, "Diag")

            !update mask
            call talsh_set_mask_done(Mask,l,k)

            m=m+1

            if (m+1>nao*nao) then
                write(*,*) 'WARNING: Cholesky decomposition did not work'
                write(*,*) 'WARNING: Final error:',D_el
                exit
            end if

            !find next element in this shell
            call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
            if (print_level.gt.6) write(*,*)  '--- largest value found: ', l, k, D_el

          end do

          !find shell with next largest element
          call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)

          ierr=talsh_tensor_destruct(Dlk)

          if (m+1>nao*nao) then
                exit
          end if

        end do
        if (print_level.gt.-1) write(*,*) ' Cholesky decomposition done, used ',m,' elements'

        !clean up
        ierr=talsh_tensor_destruct(Mask)
        ierr=talsh_tensor_destruct(Diag)
        ierr=talsh_tensor_destruct(Res)
        ierr=talsh_tensor_destruct(Temp)
        ierr=talsh_tensor_destruct(one_tensor)

        if (print_level.gt.2) write(*,*)  "========= leaving decompose_cholesky_talsh ========= "

    end subroutine decompose_cholesky_talsh_vec

    subroutine decompose_cholesky_talsh_all (CD, t_cholesky, nao, print_level)

!       Cholesky decomposition based on treating the shells seperately using talsh tensor operations
!       possible additions: 
!                          SVD (How to combine with tensor library)
!                          Additional screening in the computation of integrals

!       written by Johann Pototschnig, Summer 2019

        use exacorr_global

        implicit none

        real(8), intent(in)                  :: t_cholesky ! cholesky threshold
        type(talsh_tens_t), intent(inout)    :: CD
        integer, intent(in)                  :: nao,print_level

        integer(C_INT)     :: ierr
        integer            :: lshell, kshell, lsize, ksize
        integer            :: l, k, l0, k0
        complex(8)         :: D_max, D_el, denom_fac
        type(talsh_tens_t) :: Mask, Diag, Res, Temp
        type(talsh_tens_t) :: Dlk, CD_v, CD_m, CD_int
        integer(C_INT)     :: dims2(1:2)
        integer(C_INT)     :: dims3(1:3), off3(1:3)
        integer(C_INT)     :: dims4(1:4), off4(1:4)
        integer(C_INT)     :: rank
        integer(C_INT)     :: m
        type(talsh_tens_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)

        ! accuracy for tensor printing
        real(8), parameter :: ACC=1.0E-15

        one_dims(1) = 1
        ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)

        if (print_level.gt.2) write(*,*) "========= entered decompose_cholesky_talsh ========= "

        dims2=nao
        dims3=nao
        dims4=nao
        off4=0
        off3=0

        !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
        call interest_initialize(.false.)

        ierr=talsh_tensor_construct(Mask,R4,dims2,init_val=ONE)
        ierr=talsh_tensor_construct(Diag,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Res,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Temp,C8,dims2,init_val=ZERO)


        call talsh_get_diag(Diag)
        if (print_level.gt.2) write(*,*) "--- diagonal computed "

        call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)
        if (print_level.gt.6) write(*,*) "--- shell found: ", lshell, kshell, lsize, ksize, D_max

        m=0
        dims3(3)=1
        ierr=talsh_tensor_construct(CD,C8,dims3,init_val=ZERO)
        ierr=talsh_tensor_construct(CD_int,C8,dims3,init_val=ZERO)

        do while (abs(D_max)>t_cholesky)
          lsize=get_shell_size(lshell)
          ksize=get_shell_size(kshell)
          if (print_level.gt.2) then 
            write (*,'(a,i6,a,a,i6,i6,a,e12.2)') 'cholesky length =',m,' largest element:', &
            ' shell = ', lshell, kshell, ' value = ', abs(D_max)
          end if

          dims4(3)=lsize
          dims4(4)=ksize

          ierr=talsh_tensor_construct(Dlk,C8,dims4,init_val=ZERO)

          ierr=talsh_tensor_destruct(CD_int)
          dims3(3)=m+lsize*ksize
          ierr=talsh_tensor_construct(CD_int,C8,dims3,init_val=ZERO)
          if (m>0) then
            ierr = talsh_tensor_dimensions(CD,rank,dims3)
            off3=0
            ierr=talsh_tensor_insert(CD_int,CD,off3)
          end if

          !get the values for the current shell
          call talsh_get_shell(Dlk, lshell, kshell)
          if (print_level.gt.15) call print_tensor(Dlk, ACC, "Dlk")

          call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
          if (print_level.gt.6) write(*,*) '--- largest value found: ', l, k, D_el

          do while(abs(D_el)>abs(D_max/1000))
            call get_shell_offset(l,l0)
            call get_shell_offset(k,k0)

            !sum up previous CD vectors for the indices
            ierr=talsh_tensor_init(Res, ZERO)
            if (m>0) then
              dims3(1:2)=1
              dims3(3)=m
              ierr=talsh_tensor_construct(CD_v,C8,dims3,init_val=ZERO)
              dims3(1:2)=nao
              off3(1)=l-1
              off3(2)=k-1
              off3(3)=0
              ierr=talsh_tensor_slice(CD_v,CD_int,off3,0,DEV_HOST)
              ierr=talsh_tensor_reshape(CD_v,(/ m /))

              dims3(3)=m
              ierr=talsh_tensor_construct(CD_m,C8,dims3,init_val=ZERO)
              off3=0
              ierr=talsh_tensor_slice(CD_m,CD_int,off3,0,DEV_HOST)
              if (print_level.gt.15) call print_tensor(CD_m, ACC, "CD_m")

              ierr=talsh_tensor_contract("T(p,q)+=C(p,q,u)*D(u)",Res,CD_m,CD_v)

              ierr=talsh_tensor_destruct(CD_m)
              ierr=talsh_tensor_destruct(CD_v)
            end if

            !get Dlk for the current l and k
            ierr=talsh_tensor_init(Temp, ZERO)
            dims4(3)=1
            dims4(4)=1
            ierr=talsh_tensor_reshape(Temp,dims4)
            off4=0
            off4(3)=l-l0-1
            off4(4)=k-k0-1
            ierr=talsh_tensor_slice(Temp,Dlk,off4,0,DEV_HOST)
            ierr=talsh_tensor_reshape(Temp,dims2)

            !sum and denominate
            ierr=talsh_tensor_contract("T(p,q)+=V(p,q)",Temp,Res,one_tensor,scale=MINUS_ONE)
            ierr=talsh_tensor_init(Res, ZERO)
            ! aspg, have to replace zsqrt by a general routine, this is a gnu extension
            !denom_fac=ONE / zsqrt(D_el)
            denom_fac=ONE / sqrt(D_el)
            ierr=talsh_tensor_add("W(p,q)=V(p,q)",Res,Temp,scale=denom_fac)
            if (print_level.gt.15) call print_tensor(Res, ACC, "CD entry")

            !add Cholesky vector
            dims3(3)=1
            ierr=talsh_tensor_reshape(Res,dims3)
            off3(3)=m
            ierr=talsh_tensor_insert(CD_int,Res,off3)
            ierr=talsh_tensor_reshape(Res,dims2)

            !update diagonal
            call talsh_square_el(Res)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Diag,Res,one_tensor,scale=MINUS_ONE)
            if (print_level.gt.15) call print_tensor(Diag, ACC, "Diag")

            !update mask
            call talsh_set_mask_done(Mask,l,k)

            m=m+1

            !update CD
            ierr=talsh_tensor_destruct(CD)
            dims3(3)=m
            ierr=talsh_tensor_construct(CD,C8,dims3,init_val=ZERO)
            off3=0
            ierr=talsh_tensor_slice(CD,CD_int,off3,0,DEV_HOST)
            if (print_level.gt.15) call print_tensor(CD, ACC, "CD")

            if (m>nao*nao) write(*,*) 'Cholesky decomposition did not work'

            !find next element in this shell
            call talsh_find_el(Diag, Mask, lshell, kshell, l, k, D_el)
            if (print_level.gt.6) write(*,*)  '--- largest value found: ', l, k, D_el

          end do

          !find shell with next largest element
          call talsh_find_shells(Diag, Mask, lshell, kshell, D_max)

          ierr=talsh_tensor_destruct(Dlk)

        end do

        if (print_level.gt.-1) write(*,*) 'Cholesky decomposition done, used ',m,' elements'

        !clean up
        ierr=talsh_tensor_destruct(CD_int)

        ierr=talsh_tensor_destruct(Mask)
        ierr=talsh_tensor_destruct(Diag)
        ierr=talsh_tensor_destruct(Res)
        ierr=talsh_tensor_destruct(Temp)
        ierr=talsh_tensor_destruct(one_tensor)

        if (print_level.gt.2) write(*,*)  "========= leaving decompose_cholesky_talsh ========= "

    end subroutine decompose_cholesky_talsh_all

subroutine talsh_get_diag(Diag)

! routine to compute diagonal 2D tensor

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use exacorr_eri

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag

  real(8), pointer     :: DmR(:,:)
  complex(8),pointer   :: DmC(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2)
  integer              :: ierr
  integer              :: DataKind(1),nData

  type(basis_func_info_t), allocatable :: gto(:)  ! arrays with basis function information
  integer                              :: nao, nshells
  integer                              :: basis_angular

  integer :: lq, lp, nq, np, q, p
  integer :: i, j, ijij
  integer :: ioff, joff
  real(8) :: xq,yq,zq,xp,yp,zp
  real(8),save :: gout(21*21*21*21) ! hardwired for maximum  l-value of 5 ? (s=1,p=3,d=6,f=10,g=15,h=21)
  real(8), parameter :: FIJKL = 1.0
  
  ! get acces to the tensor after determining its type and dimensions
  rank=talsh_tensor_rank(Diag)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Diag,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_data_kind(Diag,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Diag,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DmR,dims2)
  case (C8)
     call c_f_pointer(body_p,DmC,dims2)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

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
            if (DataKind(1)==R8) then
               DmR(i+ioff,j+joff) = gout(ijij)
            else
               DmC(i+ioff,j+joff) = dcmplx(gout(ijij),0.D0)
            end if
          end do
        end do

        ioff = ioff + np
     end do
     joff = joff + nq
  end do

end subroutine talsh_get_diag

subroutine talsh_find_shells(Diag, Filter, lshell, kshell, D_max)

! routine to find shell with largest element

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global

  implicit none

  type(talsh_tens_t), intent(inout) :: Diag
  type(talsh_tens_t), intent(inout) :: Filter
  integer, intent(out)              :: lshell, kshell
  complex(8), intent(out)           :: D_max


  integer              :: nData,DataKind(1)
  real(8), pointer     :: DmR(:,:)
  complex(8), pointer  :: DmC(:,:)
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
  ierr=talsh_tensor_data_kind(Diag,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Diag,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DmR,dims2)
  case (C8)
     call c_f_pointer(body_p,DmC,dims2)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

  rank=talsh_tensor_rank(Filter)
  if (rank.ne.2) stop 'error: wrong rank in find_shells_talsh'
  ierr = talsh_tensor_dimensions(Filter,rank,dims2)
  if (ierr.ne.0) stop 'error: wrong dimension in find_shells_talsh'
  ierr=talsh_tensor_get_body_access(Filter,body_p,R4,0,DEV_HOST)
  call c_f_pointer(body_p,rfilt,dims2)
  allocate(lfilt(dims2(1),dims2(2)))
  lfilt=rfilt.gt.0

  if (DataKind(1)==R8) then
     !find largest element
     i_max=maxloc(abs(DmR))
     D_tot=DmR(i_max(1), i_max(2))

     !find largest non treated element
     i_max=maxloc(abs(DmR), mask=lfilt)
     D_max=dcmplx(DmR(i_max(1), i_max(2)),0.D0)
  else
     !find largest element
     i_max=maxloc(abs(DmC))
     D_tot=DmC(i_max(1), i_max(2))

     !find largest non treated element
     i_max=maxloc(abs(DmC), mask=lfilt)
     D_max=DmC(i_max(1), i_max(2))
  end if

  if (abs(D_tot).gt.abs(D_max)) then
    print *, 'Warning: error in cholesky will not be smaller than ', D_tot
  end if

  lshell=get_shell_index(i_max(1))
  kshell=get_shell_index(i_max(2))

  deallocate(lfilt)

end subroutine talsh_find_shells

subroutine talsh_get_shell(Dshell, lsh, ksh)

! routine to compute all ao integrals for a shell pair

  use tensor_algebra
  use talsh
  use exacorr_datatypes
  use exacorr_global
  use exacorr_eri

  implicit none 

  integer, intent(in)            :: lsh, ksh
  type(talsh_tens_t), intent(inout) :: Dshell

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
  rank=talsh_tensor_rank(Dshell)
  if (rank.ne.4) stop 'error: wrong rank in get_shell_talsh'
  ierr = talsh_tensor_dimensions(Dshell,rank,dims4)
  if (ierr.ne.0) stop 'error: wrong dimension in get_shell_talsh'
  ierr=talsh_tensor_data_kind(Dshell,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Dshell,body_p,DataKind(1),0,DEV_HOST)
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

end subroutine talsh_get_shell

subroutine talsh_find_el(Diag, Filter, lsh, ksh, l, k, D_el)

! routine to locate largest element in a shell

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

  real(8), pointer     :: DmR(:,:)
  complex(8), pointer  :: DmC(:,:)
  type(C_PTR)          :: body_p
  integer              :: rank
  integer              :: dims2(1:2), mind(1:2)
  integer              :: ierr
  real(4), pointer     :: rfilt(:,:)
  logical, allocatable :: lfilt(:,:)
  integer              :: nData,DataKind(1)

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
  ierr=talsh_tensor_data_kind(Diag,nData,DataKind)
  if (ierr.ne.0) stop 'error in getting DataKind'
  ierr=talsh_tensor_get_body_access(Diag,body_p,DataKind(1),0,DEV_HOST)
  select case (DataKind(1))
  case (R8)
     call c_f_pointer(body_p,DmR,dims2)
  case (C8)
     call c_f_pointer(body_p,DmC,dims2)
  case default
     print*, "wrong datakind :",DataKind(1)
  end select

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
    if (DataKind(1)==R8) then
       mind=maxloc(abs(DmR), mask=lfilt)
       k=mind(1)
       l=mind(2)
       D_el=dcmplx(DmR(l,k),0.D0)
     else
       mind=maxloc(abs(DmC), mask=lfilt)
       k=mind(1)
       l=mind(2)
       D_el=DmC(l,k)
     end if
  else
    k=0
    l=0
    D_el=ZERO
  end if
  
  deallocate(lfilt)
  deallocate(Dmask)

end subroutine talsh_find_el

subroutine talsh_square_el(h_tensor)

! squaring talsh tensor

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

! updating mask

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

end module talsh_CD_routines
