module eti_cis_density
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy
    use tensor_algebra
    use exacorr_tensor_interface
    use exacorr_datatypes
    use exacorr_cfg, only: exacc_input
    use exacorr_utils
    use eti_common_routines
    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(tensor_t) :: one_tensor

    private

    public get_cis_density_exe_exe, get_cis_density_exe_grd, reshape_cis

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine get_cis_density_exe_exe(exa_input,eVectors,tar_state,ini_state,Tra_Den_exe_exe_tensor)  
                 
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: tar_state,ini_state

    type(eti_dm_tens_t)         :: Tra_Den_exe_exe_tensor
    type(tensor_t)            :: Diagonal_tensor

    !complex(kind=8), pointer     :: Den_eff(:,:)
    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    integer(C_INT)     :: one_dims(1)
    type(tensor_t) :: CIS_eigen_tensor_left, CIS_eigen_tensor_right
    complex(kind=8), pointer, contiguous :: h_tens(:,:)
    complex(kind=8), pointer, contiguous :: diagonal_tens(:)


    integer              :: i, j
    integer              :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
    integer(INTD), dimension(2) :: vv_dims, oo_dims
    integer(INTD), dimension(2) :: cis_eig_dims

    nvir = exa_input%nvir
    nocc = exa_input%nocc
    cis_eig_dims(1) = nvir
    cis_eig_dims(2) = nocc
    vv_dims = nvir
    oo_dims = nocc
    one_dims(1) = 1

    ierr=construct_tensor_t(Tra_Den_exe_exe_tensor%prime_vv,C8,vv_dims,init_val=ZERO)
    ierr=construct_tensor_t(Tra_Den_exe_exe_tensor%prime_oo,C8,oo_dims,init_val=ZERO)
    ierr=construct_tensor_t(CIS_eigen_tensor_left,C8,cis_eig_dims,init_val=ZERO)
    ierr=construct_tensor_t(CIS_eigen_tensor_right,C8,cis_eig_dims,init_val=ZERO)
    ierr=construct_tensor_t(Diagonal_tensor,C8,one_dims(1:0),init_val=ZERO)



    !       Reshape CIS eigenvector matrix
    call reshape_cis(eVectors,CIS_eigen_tensor_left,cis_eig_dims,ini_state)
    call reshape_cis(eVectors,CIS_eigen_tensor_right,cis_eig_dims,tar_state)

    ierr=t_einsum("TDM(a,b)+=CISL(a,i)*CISR(b,i)",Tra_Den_exe_exe_tensor%prime_vv,&
                                CIS_eigen_tensor_left,CIS_eigen_tensor_right,scale=ONE)
    ierr=t_einsum("TDM(i,j)+=CISL(a,i)*CISR(a,j)",Tra_Den_exe_exe_tensor%prime_oo,&
                                CIS_eigen_tensor_left,CIS_eigen_tensor_right,scale=ONE)
    ierr=t_einsum("Dia()+=CISL(a,i)*CISR(a,i)",Diagonal_tensor,CIS_eigen_tensor_left,&
                                CIS_eigen_tensor_right,scale=ONE)

    ierr=direct_data_pointer(Diagonal_tensor,body_p)
    call c_f_pointer(body_p,diagonal_tens,one_dims) ! to use as a regular Fortran array

    ! print *, "error befre get_body_access_oo"
    ierr=direct_data_pointer(Tra_Den_exe_exe_tensor%prime_oo,body_p)
    call c_f_pointer(body_p,h_tens,oo_dims) ! to use as a regular Fortran array

    do i=1, nocc 
        h_tens(i,i) = diagonal_tens(1)-h_tens(i,i)
    end do

    ierr=destroy_tensor_t(CIS_eigen_tensor_left)
    ierr=destroy_tensor_t(CIS_eigen_tensor_right)
    ierr=destroy_tensor_t(Diagonal_tensor)

    ! call print_date('Leave CIS effective excited states density routine')
    
    
    end subroutine get_cis_density_exe_exe
!---------------------------------------------------------------------------------------------------------------

subroutine get_cis_density_exe_grd(exa_input,eVectors,N_state,Tra_Den_exe_grd_tensor)
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: N_state
    type(eti_dm_tens_t)         :: Tra_Den_exe_grd_tensor


    !complex(kind=8), pointer      :: Den_eff(:,:)
    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    type(tensor_t) :: Iden_tensor,CIS_eigen_tensor
    complex(kind=8), pointer, contiguous :: iden_tens(:,:)

    integer :: i, j, k, q
    integer              :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
    real(kind=8), allocatable :: norm_cis(:)
    integer(INTD), dimension(2) :: vo_dims, oo_dims
    integer(INTD), dimension(2) :: cis_eig_dims
    

    nvir = exa_input%nvir
    nocc = exa_input%nocc
    cis_eig_dims(1) = nvir
    cis_eig_dims(2) = nocc
    vo_dims(1) = nvir
    vo_dims(2) = nocc
    oo_dims = nocc

    ierr=construct_tensor_t(Tra_Den_exe_grd_tensor%prime_vo,C8,vo_dims,init_val=ZERO)
    ierr=construct_tensor_t(Iden_tensor,C8,oo_dims,init_val=ZERO)
    ierr=construct_tensor_t(CIS_eigen_tensor,C8,cis_eig_dims,init_val=ZERO)
    

!   Reshape CIS eigenvector matrix
    call reshape_cis(eVectors,CIS_eigen_tensor,cis_eig_dims,N_state)

!   Identity matrix
    ierr=direct_data_pointer(Iden_tensor,body_p)
    call c_f_pointer(body_p,iden_tens,oo_dims) ! to use as a regular Fortran array
    iden_tens = ZERO
    do i = 1, nocc
        iden_tens(i,i) = 1
    end do

    ierr=t_einsum("TDM(a,j)+=CISR(a,i)*CISR(i,j)",Tra_Den_exe_grd_tensor%prime_vo,&
                                CIS_eigen_tensor,Iden_tensor,scale=ONE)

    ierr=destroy_tensor_t(CIS_eigen_tensor)

end subroutine get_cis_density_exe_grd

!------------------------------------------------------------------------------------------------------------

subroutine reshape_cis(eVectors,CIS_eigen_tensor,cis_eig_dims,N_state)

    complex(kind=8), pointer    :: eVectors(:,:)
    type(tensor_t)          :: CIS_eigen_tensor
    integer                     :: N_state
    integer(INTD), dimension(2) :: cis_eig_dims

    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    complex(kind=8), pointer, contiguous :: cis_eig_tens(:,:)

    integer :: i, j, q
    
    ierr=direct_data_pointer(CIS_eigen_tensor,body_p)
    call c_f_pointer(body_p,cis_eig_tens,cis_eig_dims) ! to use as a regular Fortran array
    cis_eig_tens = ZERO
        do j = 1, cis_eig_dims(2)
            do i = 1, cis_eig_dims(1)
                q = (j-1)*cis_eig_dims(1) + i
                cis_eig_tens(i,j) = eVectors(q,N_state)
            end do
        end do

end subroutine reshape_cis    

end module



