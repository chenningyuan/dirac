module talsh_cis_density
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public get_cis_density_exe_exe, get_cis_density_exe_grd

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine get_cis_density_exe_exe(exa_input,eVectors,tar_state,N_state,Tra_Den_exe_exe_tensor)  
                 
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: tar_state,N_state

    type(talsh_tens_t)            :: Tra_Den_exe_exe_tensor

    complex(kind=8), pointer      :: Den_eff(:,:)
    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    type(talsh_tens_t) :: CIS_eigen_tensor
    complex(8), pointer, contiguous :: h_tens(:,:)

    integer              :: i, j
    integer              :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
    integer(INTD), dimension(2) :: vv_dims
    integer(INTD), dimension(2) :: cis_eig_dims

    nvir = exa_input%nvir
    nocc = exa_input%nocc
    cis_eig_dims(1) = nvir
    cis_eig_dims(2) = nocc
    vv_dims = nvir

    ierr=talsh_tensor_construct(Tra_Den_exe_exe_tensor,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(CIS_eigen_tensor,C8,cis_eig_dims,init_val=ZERO)


    !       Reshape CIS eigenvector matrix
    call reshape_cis(eVectors,CIS_eigen_tensor,cis_eig_dims,N_state)

    ierr=talsh_tensor_contract("TDM(a,b)+=CISR(a,i)*CISR(b,i)",Tra_Den_exe_exe_tensor,CIS_eigen_tensor,CIS_eigen_tensor,scale=ONE)

    !        ierr=talsh_tensor_get_body_access(Tra_Den_exe_exe_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    !        call c_f_pointer(body_p,h_tens,vv_dims) ! to use as a regular Fortran array

    !        do i=1, nvir 
    !            do j =1, nvir
    !                print *, h_tens(i,j)
    !            end do
    !        end do

    ierr=talsh_tensor_destruct(CIS_eigen_tensor)

    call print_date('Leave CIS effective excited states density routine')
    
    
        end subroutine get_cis_density_exe_exe
!---------------------------------------------------------------------------------------------------------------

subroutine get_cis_density_exe_grd(exa_input,eVectors,N_state,Tra_Den_exe_grd_tensor)
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: N_state
    type(talsh_tens_t)            :: Tra_Den_exe_grd_tensor


    complex(kind=8), pointer      :: Den_eff(:,:)
    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    type(talsh_tens_t) :: Iden_tensor,CIS_eigen_tensor
    complex(8), pointer, contiguous :: h_tens(:,:), iden_tens(:,:)

    integer :: i, j, k, q
    integer              :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
    real(8), allocatable :: norm_cis(:)
    integer(INTD), dimension(2) :: oo_dims, vv_dims
    integer(INTD), dimension(2) :: cis_eig_dims
    

    ierr=talsh_tensor_construct(Tra_Den_exe_grd_tensor,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(Iden_tensor,C8,vv_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(CIS_eigen_tensor,C8,cis_eig_dims,init_val=ZERO)
    

!   Reshape CIS eigenvector matrix
    call reshape_cis(eVectors,CIS_eigen_tensor)

!   Identity matrix
    oo_dims = nocc
    ierr=talsh_tensor_get_body_access(Iden_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,iden_tens,oo_dims) ! to use as a regular Fortran array
    iden_tens = ZERO
    do i = 1, nocc
        iden_tens(i,i) = 1
    end do

    ierr=talsh_tensor_contract("TDM(a,j)+=CISR(a,i)*CISR(i,j)",Tra_Den_exe_grd_tensor,CIS_eigen_tensor,Iden_tensor,scale=ONE)

!    ierr=talsh_tensor_get_body_access(Tra_Den_exe_exe_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
!    call c_f_pointer(body_p,h_tens,vv_dims) ! to use as a regular Fortran array

!    do i=1, nvir 
!        do j =1, nvir
!            print *, h_tens(i,j)
!        end do
!    end do
    ierr=talsh_tensor_destruct(CIS_eigen_tensor)


end subroutine get_cis_density_exe_grd

!------------------------------------------------------------------------------------------------------------

subroutine reshape_cis(eVectors,CIS_eigen_tensor,cis_eig_dims,N_state)

    complex(kind=8), pointer    :: eVectors(:,:)
    type(talsh_tens_t)          :: CIS_eigen_tensor
    integer                     :: N_state
    integer(INTD), dimension(2) :: cis_eig_dims

    type(C_PTR):: body_p
    integer(C_INT)       :: ierr
    complex(8), pointer, contiguous :: cis_eig_tens(:,:)

    integer :: i, j, q
    
    ierr=talsh_tensor_get_body_access(CIS_eigen_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,cis_eig_tens,cis_eig_dims) ! to use as a regular Fortran array
    cis_eig_tens = ZERO
        do j = 1, nocc
            do i = 1, nvir
                q = (j-1)*nvir + i
                cis_eig_tens(i,j) = eVectors(q,N_state)
            end do
        end do

end subroutine reshape_cis    

!---------------------------------------------------------------------------------------------------------


