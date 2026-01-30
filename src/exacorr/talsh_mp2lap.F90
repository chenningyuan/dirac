module eti_mp2lap
!This module contains the routines needed to calculate the frozen mp2 natural orbitals 
!using Cholesky decompostion and Laplace transformation

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_utils
        use exacorr_global
        use eti_cholesky_integrals
        use eti_common_routines
        use type_laplace, only : laplace_info, deallocate_laplace
        use, intrinsic:: ISO_C_BINDING

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0)
        real(8), parameter    :: ONE_QUARTER_R=0.25D0

        private

        public eti_mp2lap_driver
        public mp2lap_get_DM

       contains

!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine eti_mp2lap_driver (exa_input)
!       This module is used to construct the MP2 density matrix by talsh 
!       Written by Johann Pototschnig

        use exacorr_cfg, only: exacc_input
        implicit none
   
        type(exacc_input), intent(in) :: exa_input

        !variables to initlize library
        integer(C_INT)      :: ierr, rank

        !variables to extract MP2 energy
        type(C_PTR)         :: body_p
        type(tensor_t)  :: mp2_tensor
        integer(INTD)       :: mp2_dims(1)
        complex(8), pointer :: mp2_tens(:)
        real(8)             :: mp2_energy

        !cholesky vector
        integer, dimension(3)           :: uov_dims, uvo_dims
        type(tensor_t)              :: uov_tensor, uvo_tensor

        !Laplace denominator
        integer                                :: nocc, nvir, nchol, i, i1
        real(8), allocatable                   :: eps_occ(:),eps_vir(:),eps_occ_neg(:), zero_chol(:)
        type(laplace_info)                     :: laplace
        integer(INTD)                          :: dims4(1:4), dims2(1:2)
        type(tensor_t)                     :: L1,L2

        !intermediates
        type(tensor_t)                     :: temp

        !avoiding nvir
        type(tensor_t)                     :: T1,T2,R1,R2


        call print_date("Entering MP2_LAP Module")

        call eti_executor_initialize(exa_input%talsh_buff,get_default_executor())

        call interest_initialize(.true.)
        call print_date('Initialized interest library')

        nocc=exa_input%nocc
        nvir=exa_input%nvir

        !Initialize scalars to extract values
        mp2_dims(1) = 1
        ierr=construct_tensor_t(mp2_tensor,C8,mp2_dims(1:0),init_val=ZERO)
        ierr=direct_data_pointer(mp2_tensor,body_p)
        call c_f_pointer(body_p,mp2_tens,mp2_dims)

        !perform cholesk decompostion / integral transformation
        call get_cholesky_vec (nocc,nvir,exa_input%mo_occ,exa_input%mo_vir,uov_tensor, uvo_tensor, &
                               exa_input%t_cholesky,exa_input%print_level)

        ierr = get_nmodes(uov_tensor,rank)
        ierr = get_nmodes(uvo_tensor,rank)
        if (rank.eq.3) then
           ierr = get_extents_static(uov_tensor,uov_dims)
           ierr = get_extents_static(uvo_tensor,uvo_dims)
        else
           stop 'error in talsh_mp2lap_driver: tensor rank wrong'
        endif
        nchol=uov_dims(3)

        allocate (zero_chol(nchol))
        zero_chol=0.0
        
        call print_date("Cholesky vector computed / transformed")

        !Get orbital energies from the DIRAC restart file
        allocate (eps_occ(nocc))
        allocate (eps_vir(nvir))
        call get_orbital_energies(nocc,exa_input%mo_occ,eps_occ)
        call get_orbital_energies(nvir,exa_input%mo_vir,eps_vir)
        !using level shifts as recomputation requires all integrals, introduces deviations
        call shift_orbital_energy(eps_vir,eps_occ,exa_input%level_shift)
        call print_date("retrieved orbital energies")

        !Set up Laplace points
        allocate(eps_occ_neg(nocc))
        eps_occ_neg(1:nocc)=-eps_occ(1:nocc)
!       determine Laplace quadrature parameters with the minimax algorithm
        laplace%bounds(1) = 2.0d0*(MINVAL(eps_vir)-MAXVAL(eps_occ))
        laplace%bounds(2) = 2.0d0*(MAXVAL(eps_vir)-MINVAL(eps_occ))

        if (exa_input%print_level>5) then
          write(*,*) '--- Laplace bounds :', laplace%bounds
          write(*,*) '--- occ orb. :', eps_occ
          write(*,*) '--- vir orb. :', eps_vir
          write(*,*) '--- limits :',MINVAL(eps_vir),MAXVAL(eps_occ),MAXVAL(eps_vir),MINVAL(eps_occ)
        end if
        
!       get laplace points
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

        !create variables needed for loops
        ierr=construct_tensor_t(L1,C8,uov_dims,init_val=ZERO)
        ierr=construct_tensor_t(L2,C8,uvo_dims,init_val=ZERO)

        if (nvir<nchol) then
            dims4=nocc
            dims4(1:2)=nvir
            ierr=construct_tensor_t(temp,C8,dims4,init_val=ZERO)
        else
            dims4=nchol
            dims4(1:2)=nocc
            ierr=construct_tensor_t(R1,C8,dims4,init_val=ZERO)
            ierr=construct_tensor_t(R2,C8,dims4,init_val=ZERO)
            dims2=nchol
            ierr=construct_tensor_t(T1,C8,dims2,init_val=ZERO)
            ierr=construct_tensor_t(T2,C8,dims2,init_val=ZERO)
        end if

        mp2_energy=0

        do i=1,laplace%num

          ierr=initialize_tensor_t(L1,ZERO)
          ierr=t_einsum("R(i,a,k)+=O(i,a,k)",L1,uov_tensor)
          call laplace_mult(L1,laplace%val(i)/2.0d0,eps_occ_neg,eps_vir,zero_chol)

          ierr=initialize_tensor_t(L2,ZERO)
          ierr=t_einsum("R(a,i,k)+=O(a,i,k)",L2,uvo_tensor)
          call laplace_mult(L2,laplace%val(i)/2.0d0,eps_vir,eps_occ_neg,zero_chol)

          if (nvir<nchol) then
              !using nooc
              ierr=initialize_tensor_t(temp,ZERO)
              ierr=t_einsum("D(a,b,i,j)+=T+(i,a,k)*V(b,j,k)",temp,L1,L2)
              ierr=t_einsum("D(a,b,i,j)+=T+(i,b,k)*V(a,j,k)",temp,L1,L2,scale=MINUS_ONE)
              mp2_tens(1)=ZERO
              ierr=t_einsum("E()+=T+(a,b,i,j)*V(a,b,i,j)",mp2_tensor,temp,temp)
          else
              !using nchol
              mp2_tens(1)=ZERO
              ierr=initialize_tensor_t(T1,ZERO)
              ierr=t_einsum("M(k,n)+=V+(i,a,n)*T(i,a,k)",T1,L1,L1)
              ierr=initialize_tensor_t(T2,ZERO)
              ierr=t_einsum("N(k,n)+=T+(b,j,k)*V(b,j,n)",T2,L2,L2)
              ierr=t_einsum("E()+=M(k,n)*N(k,n)",mp2_tensor,T1,T2)

              ierr=initialize_tensor_t(R1,ZERO)
              ierr=t_einsum("M(i,j,k,n)+=T(b,j,k)*V(i,b,n)",R1,L2,L1)
              ierr=initialize_tensor_t(R2,ZERO)
              ierr=t_einsum("N(i,j,k,n)+=T(i,a,k)*V(a,j,n)",R2,L1,L2)
              ierr=t_einsum("E()+=M+(i,j,k,n)*N(i,j,k,n)",mp2_tensor,R1,R2,scale=MINUS_ONE)

              mp2_tens(1)=mp2_tens(1)*2

          end if

          if (exa_input%print_level>5) then 
              write(*,*) ' MP2 contribution : ',mp2_tens(1)
          end if

          mp2_energy=mp2_energy+mp2_tens(1)*laplace%wei(i)

        end do

        mp2_energy=-mp2_energy*ONE_QUARTER_R
        write(*,*) '    MP2 energy = ',mp2_energy,'(Laplace)'

        !clean up
        if (nvir<nchol) then
            ierr=destroy_tensor_t(temp)
        else
            ierr=destroy_tensor_t(T1)
            ierr=destroy_tensor_t(T2)
            ierr=destroy_tensor_t(R1)
            ierr=destroy_tensor_t(R2)
        end if
        ierr=destroy_tensor_t(mp2_tensor)
        !cholesky
        ierr=destroy_tensor_t(uvo_tensor)
        ierr=destroy_tensor_t(uov_tensor)
        ierr=destroy_tensor_t(L1)
        ierr=destroy_tensor_t(L2)
        !laplace
        call deallocate_laplace(laplace)
        deallocate(eps_occ)
        deallocate(eps_vir)
        deallocate(zero_chol)

        call eti_executor_shutdown

        call print_date('Leaving MP2_LAP Module')
      end subroutine eti_mp2lap_driver

      !-------------------------------------------------------------------------------------------------------------------------------------
      subroutine mp2lap_get_DM (exa_input)
!       This module is used to construct the MP2 density matrix by talsh 
!       Written by Johann Pototschnig

        use exacorr_cfg, only: exacc_input
        implicit none
   
        type(exacc_input), intent(in) :: exa_input

        !variables to initlize library
        integer(C_INT)      :: ierr, rank

        !cholesky vector
        integer, dimension(3)           :: uov_dims, uvo_dims
        type(tensor_t)              :: uov_tensor, uvo_tensor

        !Laplace denominator
        integer                                :: nocc, nvir, nchol, i, i1
        real(8), allocatable                   :: eps_occ(:),eps_vir(:),eps_occ_neg(:), zero_chol(:)
        type(laplace_info)                     :: laplace
        integer(INTD)                          :: dims4(1:4), dims2(1:2), dims1(1)
        type(tensor_t)                     :: L1p,L2p,L1q,L2q

        !intermediates
        type(tensor_t)                     :: temp, temp1, mp2_dens

        !avoiding nvir
        type(tensor_t)                     :: T1,R1,R2,S1,one_mp2lap

        call print_date("Entering MP2_LAP Module")
        
        call eti_executor_initialize(exa_input%talsh_buff,get_default_executor())

        call interest_initialize(.true.)
        call print_date('Initialized interest library')

        nocc=exa_input%nocc
        nvir=exa_input%nvir

        !Create MP2 density matrix
        dims2=nvir
        ierr=construct_tensor_t(mp2_dens,C8,dims2,init_val=ZERO)

        !perform cholesk decompostion / integral transformation
        call get_cholesky_vec (nocc,nvir,exa_input%mo_occ,exa_input%mo_vir,uov_tensor, uvo_tensor, &
                               exa_input%t_cholesky,exa_input%print_level)
        ierr = get_nmodes(uov_tensor,rank)
        ierr = get_nmodes(uvo_tensor,rank)
        if (rank.eq.3) then
           ierr = get_extents_static(uov_tensor,uov_dims)
           ierr = get_extents_static(uvo_tensor,uvo_dims)
        else
           stop 'error in eti_mp2lap_driver: tensor rank wrong'
        endif
        nchol=uov_dims(3)

        allocate (zero_chol(nchol))
        zero_chol=0.0
        
        call print_date("Cholesky vector computed / transformed")

        !Get orbital energies from the DIRAC restart file
        allocate (eps_occ(nocc))
        allocate (eps_vir(nvir))
        call get_orbital_energies(nocc,exa_input%mo_occ,eps_occ)
        call get_orbital_energies(nvir,exa_input%mo_vir,eps_vir)
        !using level shifts as recomputation requires all integrals, introduces deviations
        call shift_orbital_energy(eps_vir,eps_occ,exa_input%level_shift)
        call print_date("retrieved orbital energies")

        !Set up Laplace points
        allocate(eps_occ_neg(nocc))
        eps_occ_neg(1:nocc)=-eps_occ(1:nocc)
!       determine Laplace quadrature parameters with the minimax algorithm
        laplace%bounds(1) = 2.0d0*(MINVAL(eps_vir)-MAXVAL(eps_occ))
        laplace%bounds(2) = 2.0d0*(MAXVAL(eps_vir)-MINVAL(eps_occ))

        if (exa_input%print_level>5) then
          write(*,*) '--- Laplace bounds :', laplace%bounds
          write(*,*) '--- occ orb. :', eps_occ
          write(*,*) '--- vir orb. :', eps_vir
          write(*,*) '--- limits :',MINVAL(eps_vir),MAXVAL(eps_occ),MAXVAL(eps_vir),MINVAL(eps_occ)
        end if
        
!       get laplace points
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

        !create variables needed for loops
        ierr=construct_tensor_t(L1p,C8,uov_dims,init_val=ZERO)
        ierr=construct_tensor_t(L2p,C8,uvo_dims,init_val=ZERO)
        ierr=construct_tensor_t(L1q,C8,uov_dims,init_val=ZERO)
        ierr=construct_tensor_t(L2q,C8,uvo_dims,init_val=ZERO)

        if (nocc<nchol) then
            dims4=nocc
            dims4(1:2)=nvir
            ierr=construct_tensor_t(temp,C8,dims4,init_val=ZERO)
            ierr=construct_tensor_t(temp1,C8,dims4,init_val=ZERO)  

        else
            dims4=nchol
            dims4(1:2)=nvir
            ierr=construct_tensor_t(R1,C8,dims4,init_val=ZERO)
            ierr=construct_tensor_t(R2,C8,dims4,init_val=ZERO)
            dims2=nchol
            ierr=construct_tensor_t(T1,C8,dims2,init_val=ZERO)
            dims2=nvir
            ierr=construct_tensor_t(S1,C8,dims2,init_val=ZERO)
            dims1(1) = 1
            ierr=construct_tensor_t(one_mp2lap,C8,dims1(1:0),init_val=ONE)
        end if

        do i=1,laplace%num

          ierr=initialize_tensor_t(L1p,ZERO)
          ierr=t_einsum("R(i,a,k)+=O(i,a,k)",L1p,uov_tensor)
          call laplace_mult(L1p,laplace%val(i)/1.0d0,eps_occ_neg,eps_vir,zero_chol)

          ierr=initialize_tensor_t(L2p,ZERO)
          ierr=t_einsum("R(a,i,k)+=O(a,i,k)",L2p,uvo_tensor)
          call laplace_mult(L2p,laplace%val(i)/1.0d0,eps_vir,eps_occ_neg,zero_chol)

          if (nocc<nchol) then
              !using nooc
              ierr=initialize_tensor_t(temp,ZERO)
              ierr=t_einsum("D(a,b,i,j)+=T+(i,a,k)*V(b,j,k)",temp,L1p,L2p)
              ierr=t_einsum("D(a,b,i,j)+=T+(i,b,k)*V(a,j,k)",temp,L1p,L2p,scale=MINUS_ONE)
          end if

          ! ---------------- MP2 density matrix
          do i1=1,laplace%num

            ierr=initialize_tensor_t(L1q,ZERO)
            ierr=t_einsum("R(i,a,k)+=O(i,a,k)",L1q,uov_tensor)
            call laplace_mult(L1q,laplace%val(i1)/1.0d0,eps_occ_neg,eps_vir,zero_chol)
                
            ierr=initialize_tensor_t(L2q,ZERO)
            ierr=t_einsum("R(a,i,k)+=O(a,i,k)",L2q,uvo_tensor)
            call laplace_mult(L2q,laplace%val(i1)/1.0d0,eps_vir,eps_occ_neg,zero_chol)

            if (nocc<nchol) then
              ierr=initialize_tensor_t(temp1,ZERO)
              ierr=t_einsum("D(a,b,i,j)+=T+(i,a,k)*V(b,j,k)",temp1,L1q,L2q)
              ierr=t_einsum("D(a,b,i,j)+=T+(i,b,k)*V(a,j,k)",temp1,L1q,L2q,scale=MINUS_ONE)
                
              ierr=t_einsum("E(a,b)+=T+(a,c,i,j)*V(b,c,i,j)",mp2_dens,temp,temp1, &
                                          scale=ONE_HALF*laplace%wei(i)*laplace%wei(i1))
            else
              !using nchol
              ierr=initialize_tensor_t(R2,ZERO)
              ierr=t_einsum("N(b,a,k,n)+=T+(b,j,k)*V(a,j,n)",R2,L2p,L2q)
              ierr=initialize_tensor_t(T1,ZERO)
              ierr=t_einsum("M(k,n)+=V+(i,c,n)*T(i,c,k)",T1,L1q,L1p)
              ierr=initialize_tensor_t(S1,ZERO)
              ierr=t_einsum("E(a,b)+=M(k,n)*N(b,a,k,n)",S1,T1,R2, &
                                          scale=ONE_HALF*laplace%wei(i)*laplace%wei(i1))
              ierr=t_einsum("D(a,b)+=E(a,b)",mp2_dens,S1)
              ierr=t_einsum("D(a,b)+=E+(b,a)*K()",mp2_dens,S1,one_mp2lap)
          
              ierr=initialize_tensor_t(R1,ZERO)
              ierr=t_einsum("M(c,a,k,n)+=T+(i,a,n)*V(i,c,k)",R1,L1q,L1p)
              ierr=initialize_tensor_t(R2,ZERO)
              ierr=t_einsum("N(b,c,k,n)+=T+(b,j,k)*V(c,j,n)",R2,L2p,L2q)
              ierr=initialize_tensor_t(S1,ZERO)
              ierr=t_einsum("E(a,b)+=M(c,a,k,n)*N(b,c,k,n)",S1,R1,R2, &
                                          scale=MINUS_ONE_HALF*laplace%wei(i)*laplace%wei(i1))
              ierr=t_einsum("D(a,b)+=E(a,b)",mp2_dens,S1)
              ierr=t_einsum("D(a,b)+=E+(b,a)*K()",mp2_dens,S1,one_mp2lap)

            end if
          
          end do

        end do

        call write_2dtensor_as_matrix(nvir,mp2_dens,'DM_complex')

        !clean up
        if (nocc<nchol) then
            ierr=destroy_tensor_t(temp)
            ierr=destroy_tensor_t(temp1)
        else
            ierr=destroy_tensor_t(T1)
            ierr=destroy_tensor_t(R1)
            ierr=destroy_tensor_t(R2)
        end if
        ierr=destroy_tensor_t(mp2_dens)
        !cholesky
        ierr=destroy_tensor_t(uvo_tensor)
        ierr=destroy_tensor_t(uov_tensor)
        ierr=destroy_tensor_t(L1p)
        ierr=destroy_tensor_t(L2p)
        ierr=destroy_tensor_t(L1q)
        ierr=destroy_tensor_t(L2q)
        ierr=destroy_tensor_t(one_mp2lap)
        !laplace
        call deallocate_laplace(laplace)
        deallocate(eps_occ)
        deallocate(eps_vir)
        deallocate(zero_chol)

        call eti_executor_shutdown

        call print_date('Leaving MP2_LAP Module')
      end subroutine mp2lap_get_DM


end module eti_mp2lap
