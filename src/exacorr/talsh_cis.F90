module talsh_cis
!   This module contains the routines needed to compute CI-single wave-function 
    use tensor_algebra
    use talsh
    use exacorr_datatypes
    use exacorr_utils
    use intermediates
    use talsh_common_routines
    use talsh_cc
    implicit none
    complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER_C=(0.25D0,0.D0), MINUS_ONE_QUARTER_C=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
    real(8), parameter    :: ONE_QUARTER=0.25D0

    type(talsh_tens_t) :: one_tensor
    private

    public talsh_cis_driver
    interface diis
    module procedure diis
    end interface

    interface get_subtensor
    module procedure get_subtensor
    end interface

    contains
!---------------------------------------------------------------------------------------------------------------

    subroutine talsh_cis_driver (exa_input)
    ! This routine drives CIS calculation including construction of CIS Hamiltonian and diagnolization
        use exacorr_global
        use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
        use talsh_gradient


        type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer(C_INT)       :: ierr
         complex(8), pointer :: result_tens(:,:,:,:)
         complex(8), pointer :: hf_tens(:,:,:,:)
         type(C_PTR):: body_p
         type(C_PTR):: body_hf_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max
         integer :: i,j,a,b

!        Common tensors
         type(talsh_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens) :: int_t
!        single-body and two-body tensors
         type(talsh_tens_t) :: one-body_tensor, two-body_tensor
!        intermediate tensors (same notation as in relccsd)
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims, vovo_dims
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor
         integer(C_INT)     :: one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CIS control and result variables
         real(8) :: scf_energy, cis_energy
         real(8) :: t_target_precision
         integer :: ncycles

!        Initialize libraries
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)

        !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir

!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)

         !Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc

!         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(result_tensor,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,one_dims)

         !Constrcut one-body tensor
         ierr=talsh_tensor_construct(one-body_tensor,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(one-body_tensor,body_hf_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_hf_p,hf_tens,vovo_dims) ! to use as a regular Fortran array
         !compute HF energy and fock matrix 
         call get_CC_fockmatrix (nocc,nvir,mo_occ,mo_vir,eps_occ,eps_vir,comm_t%foo,comm_t%fov,comm_t%fvv, &
                                  int_t, scf_energy, exa_input%level_shift, exa_input%print_level)
!        set one-body tensor element 
            hf_tens = ZERO
             do b = 1, nvir
              do i = 1, nocc
               do a = 1, nvir
                do j = 1, nocc
                if (i==j .and. a==b) then
                    hf_tens(b,i,a,j) = scf_energy+eps_vir(a)-eps_occ(j)
                end if 
               end do
              end do
             end do 
            end do

         !Construct two-body tensor 
         call get_CC_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%print_level)

         ! create CIS hamiltonian tensor 
         ierr=talsh_tensor_add("HCIS(b,i,a,j)+=OB(b,i,a,j)",result_tensor,one-body_tensor)
         ierr=talsh_tensor_add("HCIS(b,i,a,j)+=TB(b,i,a,j)",result_tensor,int_t%vovo)

         call diagoanl_cis_matrix(result_tens)
         ! print CIS matrix
         !do b = 1, nvir
         !   do i = 1, nocc
         !    do a = 1, nvir
         !     do j = 1, nocc
         !       write(*,*) result_tens(b,i,a,j)
         !    end do
         !   end do
         !  end do 
         ! end do
    end subroutine talsh_cis_driver 


    
!----------------------------------------------------------------------------------------------------------------------    
    subroutine diagoanl_cis_matrix


    end subroutine diagoanl_cis_matrix


