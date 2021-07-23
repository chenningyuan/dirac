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


    type(talsh_tens_t)        :: Tra_Den_exe_exe_tensor, Tra_Den_exe_grd_tensor
    type(talsh_tens_t)        :: pro_orbital_integral_exe_grd_tensor
    type(talsh_tens_t)        :: pro_orbital_integral_exe_exe_tensor
    complex(16)               :: property_N_state, property_All_state

    call print_date('Enter CIS transition property routine')
    call get_prop_integral(exa_input,pro_orbital_integral_exe_grd_tensor,&
                            pro_orbital_integral_exe_exe_tensor)

    Total_state = size(eVectors,2)
    property_All_state = (0.0,0.0)

    do N_state =1, Total_state
        call get_cis_density_exe_exe(exa_input,eVectors,tar_state,N_state,Tra_Den_exe_exe_tensor)
        call get_cis_density_exe_grd(exa_input,eVectors,N_state,Tra_Den_exe_grd_tensor)
        call get_property(exa_input,eValues,pro_orbital_integral_exe_grd_tensor,& 
                        pro_orbital_integral_exe_exe_tensor, Tra_Den_exe_exe_tensor,&
                        Tra_Den_exe_grd_tensor,property_N_state)
        property_All_state = property_All_state + property_N_state
    end do 

    call print_date('Leave CIS transition property routine')    

end subroutine talsh_cis_pro_driver

!-----------------------------------------------------------------------------------------------------------------        

subroutine get_prop_integral(exa_input,orb_int_exe_grd_tensor,orb_int_exe_exe_tensor)
    ! Read proporty integrals from  MDPROP file 
    type(exacc_input), intent(in) :: exa_input
    type(talsh_tens_t)        :: orb_int_exe_grd_tensor
    type(talsh_tens_t)        :: orb_int_exe_exe_tensor
    

    type(C_PTR) :: body_p
    integer(C_INT) :: ierr
    integer(INTD), dimension(2)  :: vo_dims, vv_dims
    integer                      :: nocc, nvir
    complex(8), pointer, contiguous :: orbital_integral_exe_grd_tens(:,:)
    complex(8), pointer, contiguous :: orbital_integral_exe_exe_tens(:,:)

    nocc = exa_input%nocc
    nvir = exa_input%nvir
    vv_dims =nvir
    vo_dims(1) = nvir
    vo_dims(2) = nocc

    ierr=talsh_tensor_construct(orb_int_exe_grd_tensor,C8,vo_dims,init_val=ZERO)
    ierr=talsh_tensor_construct(orb_int_exe_exe_tensor,C8,vv_dims,init_val=ZERO)


    ierr=talsh_tensor_get_body_access(orb_int_exe_grd_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_point(body_p,orbital_integral_exe_grd_tens,vo_dims)
    call read_integral(exa_input,orbital_integral_exe_grd_tens,vo_dims)


    ierr=talsh_tensor_get_body_access(orb_int_exe_exe_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
    call c_f_point(body_p,orbital_integral_exe_exe_tens,vv_dims)
    call read_integral(exa_input,orbital_integral_exe_exe_tens,vv_dims)

end subroutine get_prop_integral

!------------------------------------------------------------------------------------------------------------------

subroutine read_integral(exa_input,orbital_integral_tens,mo_dims)

    implicit none

!      routine to read property integrals and add these to the int_1el data

    type(exacc_input), intent(in) :: exa_input
    complex(8), pointer, contiguous :: orbital_integral_tens(:,:)
    integer(INTD), dimension(2)  :: mo_dims
    integer, parameter     :: NOERROR  = 0 ! no error
    integer, parameter     :: NOTFOUND = 1 ! file not found
    integer, parameter     :: CORRUPT  = 2 ! file unreadable

    integer :: nkr, n_prop, i_prop
    complex(8), allocatable :: h_prop(:,:)

    character*32           :: label
    integer                :: ierr
    logical                :: tobe
    integer                :: luone
    integer                :: total_prop, max_prop
    integer                :: nocc, nvir, orbital_index1, orbital_index2

    call read_mdprop()
    call add_finitefield()

    ierr=NOERROR

    nocc = exa_input%nocc
    nvir = exa_input%nvir
    total_prop = nocc + nvir 
    max_prop = 100

    if (mo_dims(1)==nocc) then
        orbital_index1 = 1
    else 
        orbital_index1 = nocc/2 + 1
    end if 

    if (mo_dims(2)==nocc) then
        orbital_index2 = 1
    else 
        orbital_index2 = nocc/2 + 1
    end if 


!      open file and read setup infromation
    inquire (file='MDPROP',exist=tobe)
    if (.not.tobe) then
       print*, "error: file MDPROP not found"
       ierr = NOTFOUND
       return
    else
       call get_free_fileunit(luone)
       open (luone,file='MDPROP',status='OLD',FORM='UNFORMATTED',
  &     access='SEQUENTIAL')
    end if

!      allocate auxilliary array to read in
    allocate(h_prop(total_prop, total_prop))

!      count the number of matrices that we need to store
!      MDROP has a very simple format: record with label followed by record containing the data
    n_prop = 0
    do while (n_prop <= max_prop)
       read (luone,err=10,end=10) label
       if (label(25:32) == 'EOFLABEL') exit
       read (luone,err=10,end=10) h_prop
       n_prop = n_prop + 1
    end do
    rewind (luone)


!      store property integrals
    if (n_prop == exa_input%nocc + exa_input%nvir) then
        allocate(orbital_integral_tens(mo_dims(1),mo_dims(2)))
        nkr = total_prop / 2 ! number of Kramers pairs is always half the number of spinors

        i_prop = orbital_index1 - 1
        do while (i_prop < n_prop)
            read (luone,err=10,end=10) label
            i_prop = i_prop + 1
            read (luone,err=10,end=10) h_prop
            ! reorder and store in the data type
            orbital_integral_tens(1:mo_dims(1):2,1:mo_dims(2):2)
        &         = h_prop(orbital_index1:nkr,orbital_index2:nkr)     ! unbar-unbar
            orbital_integral_tens(2:mo_dims(1):2,1:mo_dims(2):2)
        &         = h_prop(nkr+orbital_index1:total_prop,orbital_index2:nkr)   ! bar-unbar
            orbital_integral_tens(1:mo_dims(1):2,2:mo_dims(2):2)
        &         = h_prop(orbital_index1:nkr,nkr+orbital_index2:total_prop)   ! unbar-bar
            orbital_integral_tens(2:mo_dims(1):2,2:mo_dims(2):2)
        &         = h_prop(nkr+orbital_index1:total_prop,nkr+orbital_index2:total_prop) ! bar-bar
        end do

        close (luone,status='KEEP')
        deallocate (h_prop)
    else
        print *, "Error: The number of correlated orbitals in MOLTRA is consisten with Exacorr!!"
    end if 



end subroutine read_integral  
!-------------------------------------------------------------------------------------------------------------------

subroutine get_property(exa_input,eValues,pro_orbital_integral_exe_grd_tensor,& 
                pro_orbital_integral_exe_exe_tensor, Tra_Den_exe_exe_tensor,&
                Tra_Den_exe_grd_tensor,property_N_state)



end subroutine get_property