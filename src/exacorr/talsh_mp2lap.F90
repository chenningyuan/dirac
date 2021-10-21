module talsh_mp2lap
!This module contains the routines needed to calculate the frozen mp2 natural orbitals 
!using Cholesky decompostion and Laplace transformation

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils

        use, intrinsic:: ISO_C_BINDING

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0)
        real(8), parameter    :: ONE_QUARTER_R=0.25D0

        private

        public talsh_mp2lap_driver

       contains

!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine talsh_mp2lap_driver (exa_input)
!       This module is used to construct the MP2 density matrix by talsh 
!       Written by Johann Pototschnig

        use talsh_cholesky_integrals
        use type_laplace, only : laplace_info, deallocate_laplace

        implicit none
   
        type(exacc_input), intent(in) :: exa_input

        !variables to initlize library
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max

        integer(C_INT)      :: ierr, rank

        !variables to extract MP2 energy
        type(C_PTR)         :: body_p
        type(talsh_tens_t)  :: mp2_tensor
        integer(INTD)       :: mp2_dims(1)
        complex(8), pointer :: mp2_tens(:)
        real(8)             :: mp2_energy

        !cholesky vector
        integer, dimension(3)           :: uov_dims, uvo_dims
        type(talsh_tens_t)              :: uov_tensor, uvo_tensor

        !Laplace denominator
        integer                                :: nocc, nvir, nchol, i
        real(8), allocatable                   :: eps_occ(:),eps_vir(:),eps_occ_neg(:), zero_chol(:)
        type(laplace_info)                     :: laplace
        integer, dimension(2)                  :: dims2
        type(talsh_tens_t)                     :: L1,L2,L3

        call print_date("Entering MP2_LAP Module")

        !Initialize TALSH
        buf_size=exa_input%talsh_buff*buf_size
        ierr=talsh_init(buf_size,host_arg_max)
        call print_date('Initialized talsh library')
        write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

        call interest_initialize(.true.)
        call print_date('Initialized interest library')

        !Initialize scalars to extract values
        mp2_dims(1) = 1
        ierr=talsh_tensor_construct(mp2_tensor,C8,mp2_dims(1:0),init_val=ZERO)
        ierr=talsh_tensor_get_body_access(mp2_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,mp2_tens,mp2_dims)

        nocc=exa_input%nocc
        nvir=exa_input%nvir

        !perform cholesk decompostion / integral transformation
        call get_cholesky_vec (nocc,nvir,exa_input%mo_occ,exa_input%mo_vir,uov_tensor, uvo_tensor, &
                               exa_input%t_cholesky,exa_input%print_level)
        rank=talsh_tensor_rank(uov_tensor)
        if (rank.eq.3) then
           ierr = talsh_tensor_dimensions(uov_tensor,rank,uov_dims)
           ierr = talsh_tensor_dimensions(uvo_tensor,rank,uvo_dims)
        else
           stop 'error in talsh_mp2lap_driver: tensor rank wrong'
        endif
        write(*,*) 'uov_dims', uov_dims
        nchol=uov_dims(0)
        allocate (zero_chol(uov_dims(3)))
        call print_date("Cholesky vector computed / transformed")

        !Get orbital energies from the DIRAC restart file
        allocate (eps_occ(nocc))
        allocate (eps_vir(nvir))
        call get_orbital_energies(nocc,exa_input%mo_occ,eps_occ)
        call get_orbital_energies(nvir,exa_input%mo_vir,eps_vir)
        call print_date("retrieved orbital energies")

        !Set up Laplace points
        allocate(eps_occ_neg(nocc))
        eps_occ_neg(1:nocc)=-eps_occ(1:nocc)
!       determine Laplace quadrature parameters with the minimax algorithm
        laplace%bounds(1) = 2.0d0*(eps_vir(1)-eps_occ(nocc))
        laplace%bounds(2) = 2.0d0*(eps_vir(nvir)-eps_occ(1))
        if (exa_input%print_level>0) then
          write(*,*) '--- Laplace bounds :', laplace%bounds
          write(*,*) '--- vir orb. :', eps_vir
          write(*,*) '--- occ orb. :', eps_occ
        end if
        
!       get laplace points
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

!       start computation
        mp2_energy=0

        do i=1,laplace%num

          ierr=talsh_tensor_construct(L1,C8,uov_dims,init_val=ZERO)
          ierr=talsh_tensor_add("R(k,i,a)+=O(k,i,a)",L1,uov_tensor)
          call laplace_mult(L1,laplace%val(i)/2.0d0,zero_chol,eps_occ_neg,eps_vir)
          

          ierr=talsh_tensor_construct(L2,C8,uvo_dims,init_val=ZERO)   
          ierr=talsh_tensor_add("R(k,a,i)+=O(k,a,i)",L2,uvo_tensor)
          call laplace_mult(L2,laplace%val(i)/2.0d0,zero_chol,eps_vir,eps_occ_neg)

          mp2_tens(1)=ZERO
          ierr=talsh_tensor_contract("E()+=T(k,i,a)*V(k,a,i)",mp2_tensor,L1,L2)
          if (exa_input%print_level>0) write(*,*) ' MP2 contribution: ',mp2_tens(1)

          ierr=talsh_tensor_destruct(L1)
          ierr=talsh_tensor_destruct(L2)

          mp2_energy=mp2_energy+laplace%wei(i)*mp2_tens(1)
           
        end do

        mp2_energy=-mp2_energy*ONE_QUARTER_R
        !minus because the regular denominator is negative
        !which is of course not possible with Laplace
        write(*,*) '    MP2 energy = ',mp2_energy,'(Laplace)'


        !clean up
        ierr=talsh_tensor_destruct(uvo_tensor)
        ierr=talsh_tensor_destruct(uov_tensor)
        call deallocate_laplace(laplace)
        deallocate(eps_occ, eps_vir)


        !shutdown TALSH
        ierr=talsh_stats()
        ierr = talsh_shutdown()

        call print_date('Leaving MP2_LAP Module')
      end subroutine talsh_mp2lap_driver


end module talsh_mp2lap
