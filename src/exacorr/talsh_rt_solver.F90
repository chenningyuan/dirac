module talsh_rt_solver
!   This module contains the routines needed to solve real-time equations
! Written by Xiang Yuan, June 2024
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines



    contains

!---------------------------------------------------------------------------------------------------------------

subroutine talsh_cis_rt_solver_driver(exa_input)

    use exacorr_global
    use talsh_ao_to_mo


    type(exacc_input), intent(in) :: exa_input


    integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
    integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
    integer(C_INT)       :: ierr
    complex(8), pointer :: result_tens(:,:,:,:),sigma1(:,:)
    integer(INTD)       :: fo_dims(1:2),fv_dims(1:2)
    type(C_PTR):: body_p
    integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
    integer(C_INT):: host_arg_max
    complex(8), pointer, contiguous :: h_tens(:,:)
    integer :: i,j,a,b,N_step
    integer :: bj,ai


!   fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
    type(talsh_intg_tens_t) :: int_t
!   single-body and body_two tensors
    type(talsh_tens_t) :: body_one_tensor, diag_foo_tensor, diag_fvv_tensor
    type(talsh_tens_t) :: id_mat_occ, id_mat_vir
!   intermediate tensors (same notation as in relccsd)
    integer(INTD), dimension(2) :: oo_dims, vv_dims
    integer(INTD), dimension(4) :: oooo_dims, vovo_dims, voov_dims
!   scalars (need to be defined as tensor types)
    type(talsh_tens_t) :: result_tensor,sigma_tensor
    integer(C_INT)     :: one_dims(1)
!   orbital energies
    real(8), allocatable :: eps_occ(:),eps_vir(:)


    call get_integrals()
    call construct_initial_cis_tensors()

!        Get orbital energy
    call get_orbital_energies(nocc,mo_occ,eps_occ)
    call get_orbital_energies(nvir,mo_vir,eps_vir)

    do i=1, N_step
        call get_F()
        call get_cis_tensors()
        call get_cis_0_tensor()
    end do

!---------------------------------------------------------------------------------------------------------------
    subroutine get_integrals()


    end subroutine
!---------------------------------------------------------------------------------------------------------------

    subroutine construct_initial_cis_tensors()


    end subroutine
!---------------------------------------------------------------------------------------------------------------

    subroutine get_F()


    end subroutine 
!---------------------------------------------------------------------------------------------------------------

    subroutine get_cis_tensors()


    end subroutine 
!---------------------------------------------------------------------------------------------------------------


    subroutine get_cis_0_tensor()


    end subroutine 
!---------------------------------------------------------------------------------------------------------------



end subroutine 
