module talsh_cis_property
!   This module contains the routines needed to compute CI-single wave-function and CIS excitation energy
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use talsh_common_routines
    use talsh_cis_density
    use exacorr_global
    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                ONE_HALF=(0.5D0,0.D0), &
                                MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor

    private

    public talsh_cis_pro_driver

    contains

!-----------------------------------------------------------------------------------------------------------------

subroutine talsh_cis_pro_driver(exa_input,eValues,eVectors,tar_state)
!   Written by Xiang Yuan July 2021, Amsterdam, Netherland  
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    complex(kind=8), pointer      :: eValues(:) 
    integer                       :: tar_state, Total_state, N_state
    integer                       :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
 

    type(talsh_tens_t)        :: Tra_Den_exe_exe_tensor, Tra_Den_exe_grd_tensor
    type(talsh_tens_t)        :: pro1_orbital_integral_exe_grd_tensor
    type(talsh_tens_t)        :: pro1_orbital_integral_exe_exe_tensor
    type(talsh_tens_t)        :: pro2_orbital_integral_exe_grd_tensor
    type(talsh_tens_t)        :: pro2_orbital_integral_exe_exe_tensor
    complex(kind=8)           :: term1_N_state, term2_N_state, property_All_state
    integer                   :: iff

    nvir = exa_input%nvir
    nocc = exa_input%nocc
    
    call print_date('Enter CIS transition property routine')
    call get_prop_integral(exa_input,pro1_orbital_integral_exe_grd_tensor,& 
                            nvir,exa_input%mo_vir,nocc,exa_input%mo_occ,iff)

    call get_prop_integral(exa_input,pro1_orbital_integral_exe_exe_tensor,&
                            nvir,exa_input%mo_vir,nvir,exa_input%mo_vir,iff)

    call get_prop_integral(exa_input,pro2_orbital_integral_exe_grd_tensor,& 
                            nvir,exa_input%mo_vir,nocc,exa_input%mo_occ,iff)

    call get_prop_integral(exa_input,pro2_orbital_integral_exe_exe_tensor,&
                            nvir,exa_input%mo_vir,nvir,exa_input%mo_vir,iff)

    Total_state = size(eVectors,2)
    property_All_state = (0.0,0.0)

    do N_state =1, Total_state
        call get_cis_density_exe_exe(exa_input,eVectors,tar_state,N_state,Tra_Den_exe_exe_tensor)
        call get_cis_density_exe_grd(exa_input,eVectors,N_state,Tra_Den_exe_grd_tensor)
        call get_property(exa_input,eValues,N_state,&
                        pro1_orbital_integral_exe_grd_tensor,pro2_orbital_integral_exe_exe_tensor,& 
                        Tra_Den_exe_exe_tensor,Tra_Den_exe_grd_tensor,term1_N_state)
        call get_property(exa_input,eValues,N_state,&
                        pro2_orbital_integral_exe_grd_tensor,pro1_orbital_integral_exe_exe_tensor,& 
                        Tra_Den_exe_exe_tensor,Tra_Den_exe_grd_tensor,term2_N_state)
        property_All_state = property_All_state + term1_N_state + term2_N_state
    end do 

    call print_date('Leave CIS transition property routine')    

end subroutine talsh_cis_pro_driver

!-----------------------------------------------------------------------------------------------------------------        

subroutine get_prop_integral(exa_input,orb_int_tensor,nmo1,mo1_list,nmo2,mo2_list,iff)
    ! Read proporty integrals from  MDPROP file 
    type(exacc_input), intent(in) :: exa_input
    type(talsh_tens_t)        :: orb_int_tensor
    

    type(C_PTR) :: body_p
    integer(C_INT) :: ierr
    integer                      :: nocc, nvir
    integer                      :: iff
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
    call add_finitefield(orbital_integral_tens,exa_input,iff,mo1_list,nmo1,mo2_list,nmo2,min_occ)


end subroutine get_prop_integral
!-------------------------------------------------------------------------------------------------------------------

subroutine get_property(exa_input,eValues,N_state,pro_orbital_integral_exe_grd_tensor,& 
                pro_orbital_integral_exe_exe_tensor, Tra_Den_exe_exe_tensor,&
                Tra_Den_exe_grd_tensor,property_N_state)

    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eValues(:) 
    integer(C_INT) :: ierr
    integer(C_INT)     :: one_dims(1)
    integer :: N_state


    type(talsh_tens_t)        :: Tra_Den_exe_exe_tensor, Tra_Den_exe_grd_tensor
    type(talsh_tens_t)        :: pro_orbital_integral_exe_grd_tensor
    type(talsh_tens_t)        :: pro_orbital_integral_exe_exe_tensor
    type(talsh_tens_t)        :: TDM_exe_grd_tensor, TDM_exe_exe_tensor

    complex(kind=8), pointer      :: TDM_exe_grd_tens(:), TDM_exe_exe_tens(:)
    type(C_PTR)                   :: body_p
    complex(kind=8)               :: property_N_state
    complex(kind=8)               :: ene_exc



    one_dims(1) = 1
    ene_exc = eValues(N_state)

    ierr=talsh_tensor_get_body_access(TDM_exe_grd_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,TDM_exe_grd_tens,one_dims)
    ierr = talsh_tensor_contract("TDM+=OL(a,i)*DM(a,i)",TDM_exe_grd_tensor,pro_orbital_integral_exe_grd_tensor,&
                                Tra_Den_exe_grd_tensor,scale=ONE)


    ierr=talsh_tensor_get_body_access(TDM_exe_exe_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,TDM_exe_exe_tens,one_dims)
    ierr = talsh_tensor_contract("TDM+=OL(a,b)*DM(a,b)",TDM_exe_exe_tensor,pro_orbital_integral_exe_exe_tensor,&
                                Tra_Den_exe_exe_tensor,scale=ONE)


    property_N_state = TDM_exe_grd_tens(1)*TDM_exe_exe_tens(1)/ene_exc

end subroutine get_property

end module