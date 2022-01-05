module talsh_eom_QR
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use exacorr_global
    use intermediates_eom_rsp
    use talsh_eom_rspeq_solver
    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public talsh_eom_QR_driver

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine talsh_eom_QR_driver(exa_input,eValues,eVectors)
!   Written by Xiang Yuan July 2021, Amsterdam, Netherland 

    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    complex(kind=8), pointer      :: eValues(:) 
    integer                       :: All_final_states
    integer                       :: All_opera
    integer                       :: f_state
    complex(kind=8), pointer      :: onetdmee(:,:,:)
    complex(kind=8), pointer      :: onetdmeg(:,:)
    complex(kind=8)               :: TPA_one_state
    complex(kind=8), pointer      :: TPA(:)

    
    call print_date('Enter CC-CI quadratic response')
    
    call form_eom_quarsp_intermediates()
    call solve_eom_linrsp_left()
    call solve_eom_linrsp_right()
    call construct_eom_quarsp()
    call print_date('Leave CC-CI quadratic response')    

end subroutine talsh_eom_QR_driver


!------------------------------------------------------------------------------------------------------------------

subroutine get_prop_integral(exa_input,orb_int_tensor,nmo1,mo1_list,nmo2,mo2_list,i_prop)
    ! Read proporty integrals from  MDPROP file 
    type(exacc_input), intent(in) :: exa_input
    type(talsh_tens_t)        :: orb_int_tensor
    

    type(C_PTR) :: body_p
    integer(C_INT) :: ierr
    integer                      :: nocc, nvir
    integer                      :: i_prop
    complex(kind=8), pointer, contiguous :: orbital_integral_tens(:,:)

    integer                :: nmo1          ! the length of the mo basis
    integer                :: mo1_list(:)   ! and their indices
    integer                :: nmo2          ! the length of the mo basis
    integer                :: mo2_list(:)   ! and their indices
    integer                :: min_occ       ! smallest index, assumes same closed shell as MOLTRA

    integer(INTD),dimension(2)   :: mo_dims


    mo_dims(1) = nmo1
    mo_dims(2) = nmo2


    ierr=talsh_tensor_construct(orb_int_tensor,C8,mo_dims,init_val=ZERO)
    ierr=talsh_tensor_get_body_access(orb_int_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,orbital_integral_tens,mo_dims)


    min_occ = minval(exa_input%mo_occ)
    !call add_finitefield(orbital_integral_tens,exa_input,iff,mo1_list,nmo1,mo2_list,nmo2,min_occ)
    call get_ff_mat (orbital_integral_tens,i_prop,mo1_list,nmo1,mo2_list,nmo2,min_occ)
    ! call print_date('Get property integrals')    



end subroutine get_prop_integral
!-------------------------------------------------------------------------------------------------------------------

end module
