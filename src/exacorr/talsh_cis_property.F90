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

subroutine talsh_cis_pro_driver(exa_input,eValues,eVectors)
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

    
    call print_date('Enter CIS transition property routine')
    All_final_states = size(eValues)
    ! All_final_states = 1
    ! the number of operator component should be determined by input or here?
    All_opera = 3 
    allocate(TPA(All_final_states))
    
    call get_all_one_photon_tdm(exa_input, onetdmeg, onetdmee, All_opera, eVectors)
    do f_state = 1, All_final_states
        call get_TPA(exa_input, TPA_one_state, onetdmeg, onetdmee, All_opera, f_state, eValues)
        TPA(f_state) = TPA_one_state
    end do

    call print_two_photon(TPA)
    call print_date('Leave CIS transition property routine')    

end subroutine talsh_cis_pro_driver

!-----------------------------------------------------------------------------------------------------------------  

subroutine get_TPA(exa_input, TPA_one_state, onetdmeg, onetdmee, All_opera, f_state, eValues)
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eValues(:) 
    complex(kind=8), pointer      :: onetdmeg(:,:)
    complex(kind=8), pointer      :: onetdmee(:,:,:)
    integer                       :: All_opera
    integer                       :: f_state
    integer                       :: ope1, ope2
    integer                       :: F, G, H
    integer                       :: S_ty
    complex(kind=8)               :: S_aabb, S_abab, S_abba
    complex(kind=8)               :: TPA_one_state

    ! call print_date('Start the CIS Two photon absorption routine')
! now we use parallel linearly polarized light, but inprinciple we could use other parameter.
    F=2
    G=2
    H=2
    TPA_one_state = (0.0,0.0)
    do ope1 = 1, All_opera
        do ope2 = 1, All_opera
            S_ty = 1
            call get_trans_streng_matx_S(S_aabb, ope1, ope2, f_state, eValues, onetdmeg, onetdmee, S_ty)
            S_ty = 2
            call get_trans_streng_matx_S(S_abab, ope1, ope2, f_state, eValues, onetdmeg, onetdmee, S_ty)
            S_ty = 3          
            call get_trans_streng_matx_S(S_abba, ope1, ope2, f_state, eValues, onetdmeg, onetdmee, S_ty)
            TPA_one_state = TPA_one_state+(F/30.0)*S_aabb + (G/30.0)*S_abab + (H/30.0)*S_abba
!             print *, "S_aabb", S_aabb
!             print *, "S_abab", S_abab
!             print *, "S_abba", S_abba
!             print *, "TPA_one_state", TPA_one_state 
        end do
    end do
    ! call print_date('Leave CIS TPA routine')

end subroutine get_TPA

!------------------------------------------------------------------------------------------------------------------

subroutine get_trans_streng_matx_S(S, ope1, ope2, f_state, eValues, onetdmeg, onetdmee, S_ty)
    complex(kind=8), pointer      :: eValues(:) 
    complex(kind=8), pointer      :: onetdmeg(:,:)
    complex(kind=8), pointer      :: onetdmee(:,:,:)
    integer                       :: f_state
    integer                       :: ope1, ope2
    integer                       :: S_ty
    complex(kind=8)               :: S
    real(kind=8)                  :: omega1, omega2

    integer                       :: inte_state
    integer                       :: all_inte_states
    integer                       :: ope1_obj, ope2_obj, ope3_obj, ope4_obj
    complex(kind=8)               :: term1, term2
    complex(kind=8)               :: M0f, Mf0


    ! call print_date('Start the Strength matrix element routine')
    all_inte_states = size(eValues)
    ! omega should be in input
    omega1= 0.0
    omega2 =0.0

    select case (S_ty)
        case (1)
            ope1_obj = ope1
            ope2_obj = ope1
            ope3_obj = ope2
            ope4_obj = ope2
        case (2)
            ope1_obj = ope1
            ope2_obj = ope2
            ope3_obj = ope1
            ope4_obj = ope2
        case (3)
            ope1_obj = ope1
            ope2_obj = ope2
            ope3_obj = ope2
            ope4_obj = ope1
    end select

    M0f = (0.0, 0.0)
    Mf0 = (0.0, 0.0)
    do inte_state = 1, all_inte_states
        term1 = -onetdmee(f_state, inte_state, ope2_obj)*onetdmeg(inte_state, ope1_obj)/(eValues(inte_state)-omega1)
        term2 = -onetdmee(f_state, inte_state, ope1_obj)*onetdmeg(inte_state, ope2_obj)/(eValues(inte_state)-omega2)
        M0f = M0f + term1 + term2
    end do

    do inte_state = 1, all_inte_states
        term1 = -onetdmee(inte_state, f_state, ope4_obj)*onetdmeg(inte_state, ope3_obj)/(eValues(inte_state)-omega1)
        term2 = -onetdmee(inte_state, f_state, ope3_obj)*onetdmeg(inte_state, ope4_obj)/(eValues(inte_state)-omega2)
        Mf0 = Mf0 + term1 + term2
    end do

    S = M0f * Mf0

    ! call print_date('Leave the Strength matrix element routine')

end subroutine get_trans_streng_matx_S

!---------------------------------------------------------------------------------------------------------------------

subroutine get_all_one_photon_tdm(exa_input, onetdmeg, onetdmee, All_opera, eVectors)

    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: All_opera
    complex(kind=8), pointer      :: onetdmee(:,:,:)
    complex(kind=8), pointer      :: onetdmeg(:,:)


    call print_date('Enter One Photon transition moments routine')
    call get_exc_grd(exa_input, onetdmeg, All_opera, eVectors)
    call get_exc_exc(exa_input, onetdmee, All_opera, eVectors)
    call print_onetdm(onetdmeg, onetdmee)
    call print_date('Leave One Photon transition moments routine')

end subroutine get_all_one_photon_tdm

!------------------------------------------------------------------------------------------------------------------

subroutine get_exc_grd(exa_input, onetdmeg, All_opera, eVectors)
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: All_opera
    complex(kind=8), pointer      :: onetdmeg(:,:)

    integer                       :: ope 
    integer                       :: f_state
    integer                       :: all_state
    integer                       :: nvir, nocc
    complex(kind=8)               :: onetdmeg_tem

    type(talsh_dm_tens_t)       :: Tra_Den_exe_grd_tensor
    type(talsh_tens_t)        :: orb_inte_exe_grd_tensor


    nvir = exa_input%nvir
    nocc = exa_input%nocc
    all_state = size(eVectors,2)
    ! call print_date('Enter exc-grd routine')
    allocate(onetdmeg(all_state,All_opera))
    do ope=1, All_opera
        do f_state =1, all_state 
            call get_prop_integral(exa_input, orb_inte_exe_grd_tensor,& 
                        nvir,exa_input%mo_vir,nocc,exa_input%mo_occ,ope)
            ! print *, "get integral"
            call get_cis_density_exe_grd(exa_input,eVectors,f_state,Tra_Den_exe_grd_tensor)

            ! print *, "get density mat"
            call get_propeg(onetdmeg_tem, orb_inte_exe_grd_tensor, Tra_Den_exe_grd_tensor)

            ! print *, "onetdmeg", onetdmeg_tem
            onetdmeg(f_state,ope) = onetdmeg_tem
        end do
    end do
    ! call print_date('Leave exc-grd routine')

end subroutine get_exc_grd

!------------------------------------------------------------------------------------------------------------------

subroutine get_exc_exc(exa_input, onetdmee, All_opera, eVectors)
    type(exacc_input), intent(in) :: exa_input
    complex(kind=8), pointer      :: eVectors(:,:)
    integer                       :: All_opera
    complex(kind=8), pointer      :: onetdmee(:,:,:)

    integer                       :: ope 
    integer                       :: f_state, ini_state
    integer                       :: all_state
    integer                       :: nvir, nocc
    complex(kind=8)               :: onetdmee_tem
    type(talsh_dm_tens_t)           :: Tra_Den_exe_exe_tensor
    type(talsh_tens_t)            :: orb_inte_exe_exe_vv_tensor
    type(talsh_tens_t)            :: orb_inte_exe_exe_oo_tensor


    nvir = exa_input%nvir
    nocc = exa_input%nocc
    all_state = size(eVectors,2)
    allocate(onetdmee(all_state,all_state,All_opera))
    do ope=1, All_opera
        do ini_state =1, all_state
            do f_state =1, all_state 
                ! print *, "error before get integ vv?"
                call get_prop_integral(exa_input,orb_inte_exe_exe_vv_tensor,&
                nvir,exa_input%mo_vir,nvir,exa_input%mo_vir,ope)
                
                ! print *, "error before get integ oo?"
                call get_prop_integral(exa_input,orb_inte_exe_exe_oo_tensor,&
                nocc,exa_input%mo_occ,nocc,exa_input%mo_occ,ope)

                ! print *, "error before getcisden?"
                call get_cis_density_exe_exe(exa_input,eVectors,f_state,ini_state,Tra_Den_exe_exe_tensor)

                ! print *, "error before gpe?"
                call get_propee(onetdmee_tem, orb_inte_exe_exe_vv_tensor, orb_inte_exe_exe_oo_tensor, Tra_Den_exe_exe_tensor)
                
                ! print *, "onetdmee", onetdmee_tem
                onetdmee(f_state,ini_state,ope) = onetdmee_tem
            end do
        end do
    end do

end subroutine get_exc_exc

!-----------------------------------------------------------------------------------------------------------------

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

subroutine get_propeg(onetdmeg_tem, orb_inte_exe_grd_tensor, Tra_Den_exe_grd_tensor)

    integer(C_INT) :: ierr
    integer(C_INT)     :: one_dims(1)
    integer :: N_state
    complex(kind=8)               :: onetdmeg_tem


    type(talsh_dm_tens_t)           :: Tra_Den_exe_grd_tensor
    type(talsh_tens_t)            :: orb_inte_exe_grd_tensor
    type(talsh_tens_t)            :: TDM_exe_grd_tensor

    complex(kind=8), pointer      :: TDM_exe_grd_tens(:)
    type(C_PTR)                   :: body_p

    one_dims(1) = 1

    ierr=talsh_tensor_construct(TDM_exe_grd_tensor,C8,one_dims(1:0),init_val=ZERO)
    ierr=talsh_tensor_get_body_access(TDM_exe_grd_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,TDM_exe_grd_tens,one_dims)

    ierr = talsh_tensor_contract("TDM()+=OL(a,i)*DM(a,i)",TDM_exe_grd_tensor,orb_inte_exe_grd_tensor,&
                                Tra_Den_exe_grd_tensor%prime_vo,scale=ONE)
    onetdmeg_tem = TDM_exe_grd_tens(1)
    ierr=talsh_tensor_destruct(TDM_exe_grd_tensor)
    ierr=talsh_tensor_destruct(orb_inte_exe_grd_tensor)
    ierr=talsh_tensor_destruct(Tra_Den_exe_grd_tensor%prime_vo)

    
end subroutine get_propeg

!------------------------------------------------------------------------------------------------------------------------

subroutine get_propee(onetdmee_tem, orb_inte_exe_exe_vv_tensor, orb_inte_exe_exe_oo_tensor, Tra_Den_exe_exe_tensor)

    integer(C_INT) :: ierr
    integer(C_INT)     :: one_dims(1)
    integer :: N_state


    complex(kind=8)           :: onetdmee_tem
    type(talsh_dm_tens_t)       :: Tra_Den_exe_exe_tensor
    type(talsh_tens_t)        :: orb_inte_exe_exe_vv_tensor
    type(talsh_tens_t)        :: orb_inte_exe_exe_oo_tensor
    type(talsh_tens_t)        :: TDM_exe_exe_vv_tensor, TDM_exe_exe_oo_tensor

    complex(kind=8), pointer      :: TDM_exe_exe_vv_tens(:), TDM_exe_exe_oo_tens(:)
    type(C_PTR)                   :: body_p

    one_dims(1) = 1


    ierr=talsh_tensor_construct(TDM_exe_exe_vv_tensor,C8,one_dims(1:0),init_val=ZERO) 
    ierr=talsh_tensor_construct(TDM_exe_exe_oo_tensor,C8,one_dims(1:0),init_val=ZERO) 


    ierr=talsh_tensor_get_body_access(TDM_exe_exe_vv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,TDM_exe_exe_vv_tens,one_dims)
    ierr = talsh_tensor_contract("TDM()+=OL(a,b)*DM(a,b)",TDM_exe_exe_vv_tensor,orb_inte_exe_exe_vv_tensor,&
                    Tra_Den_exe_exe_tensor%prime_vv,scale=ONE)

               
    ierr=talsh_tensor_get_body_access(TDM_exe_exe_oo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_pointer(body_p,TDM_exe_exe_oo_tens,one_dims)
    ierr = talsh_tensor_contract("TDM()+=OL(i,j)*DM(i,j)",TDM_exe_exe_oo_tensor,orb_inte_exe_exe_oo_tensor,&
                    Tra_Den_exe_exe_tensor%prime_oo,scale=ONE)

    onetdmee_tem = TDM_exe_exe_vv_tens(1)+TDM_exe_exe_oo_tens(1)

    ierr=talsh_tensor_destruct(TDM_exe_exe_vv_tensor)
    ierr=talsh_tensor_destruct(TDM_exe_exe_oo_tensor)
    ierr=talsh_tensor_destruct(orb_inte_exe_exe_vv_tensor)
    ierr=talsh_tensor_destruct(orb_inte_exe_exe_oo_tensor)
    ierr=talsh_tensor_destruct(Tra_Den_exe_exe_tensor%prime_vv)
    ierr=talsh_tensor_destruct(Tra_Den_exe_exe_tensor%prime_oo)
end subroutine get_propee

!-------------------------------------------------------------------------------------------------------------------------

subroutine print_onetdm(onetdmeg, onetdmee)
    complex(kind=8), pointer      :: onetdmee(:,:,:)
    complex(kind=8), pointer      :: onetdmeg(:,:)
    integer                       :: i,j,k

    call print_date("Print One photon Transition moment (Ground to Excited)")

    do i=1, size(onetdmeg,1)
        print *, "Groud to ", i, "Excited Transition moment:"
        do j=1, size(onetdmeg,2)
            print "(I5,(F27.5))",j, onetdmeg(i,j)
        end do
    end do

    ! do j=1, size(onetdmeg,2)
    !     print *, "Component:", j
    !     do i =1, size(onetdmeg,1)
    !         print *, onetdmeg(i,j)
    !     end do
    ! end do
    call print_date("Print One photon Transition moment (Excited to Excited)")
    do i=1, size(onetdmee,1)
        do j=1, size(onetdmee,2)
            print *, "Excited", i, "to", j, "Excited Transition moment:"
            do k=1, size(onetdmee,3)
                print "(I5,(F27.5))",k, onetdmee(i,j,k)
            end do
        end do
    end do

    ! do k=1, size(onetdmee,3)
    !     print *, "Component:", k
    !     do j =1, size(onetdmee,2)
    !         do i=1, size(onetdmee,1)
    !             print *, onetdmee(i,j,k)
    !         end do 
    !     end do
    ! end do


end subroutine print_onetdm
!-----------------------------------------------------------------------------------------------------------------------

subroutine print_two_photon(TPA)
        
    complex(kind=8), pointer      :: TPA(:)
    integer                       :: i

    call print_date("Print Two photon cross section") 
    do i=1, size(TPA,1)
        print *, i, "Excited Transition moment:"
            print "(I5,(F27.5))",i, abs(TPA(i))
    end do


end subroutine print_two_photon


end module
