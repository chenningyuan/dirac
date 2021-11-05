module talsh_mp2lap
!This module contains the routines needed to calculate the frozen mp2 natural orbitals 
!using Cholesky decompostion and Laplace transformation

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils
        use exacorr_global
        use talsh_cholesky_integrals
        use talsh_common_routines
        use type_laplace, only : laplace_info, deallocate_laplace

        use, intrinsic:: ISO_C_BINDING

        implicit none

        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0),MINUS_ONE_HALF=(-0.5D0,0.D0)
        real(8), parameter    :: ONE_QUARTER_R=0.25D0

        private

        public talsh_mp2lap_driver
        public mp2lap_get_DM

       contains

!-------------------------------------------------------------------------------------------------------------------------------------
       subroutine talsh_mp2lap_driver (exa_input)
!       This module is used to construct the MP2 density matrix by talsh 
!       Written by Johann Pototschnig

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
        integer                                :: nocc, nvir, nchol, i, i1
        real(8), allocatable                   :: eps_occ(:),eps_vir(:),eps_occ_neg(:), zero_chol(:)
        type(laplace_info)                     :: laplace
        integer(INTD)                          :: dims4(1:4), dims2(1:2)
        type(talsh_tens_t)                     :: L1,L2

        !intermediates
        type(talsh_tens_t)                     :: temp

        !avoiding nvir
        type(talsh_tens_t)                     :: T1,T2,R1,R2

        call print_date("Entering MP2_LAP Module")

        !Initialize TALSH
        buf_size=exa_input%talsh_buff*buf_size
        ierr=talsh_init(buf_size,host_arg_max)
        call print_date('Initialized talsh library')
        write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

        call interest_initialize(.true.)
        call print_date('Initialized interest library')

        nocc=exa_input%nocc
        nvir=exa_input%nvir

        !Initialize scalars to extract values
        mp2_dims(1) = 1
        ierr=talsh_tensor_construct(mp2_tensor,C8,mp2_dims(1:0),init_val=ZERO)
        ierr=talsh_tensor_get_body_access(mp2_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
        call c_f_pointer(body_p,mp2_tens,mp2_dims)

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
        laplace%bounds(1) = 2.0d0*(eps_vir(1)-eps_occ(nocc))
        laplace%bounds(2) = 2.0d0*(eps_vir(nvir)-eps_occ(1))
        if (exa_input%print_level>10) then
          write(*,*) '--- Laplace bounds :', laplace%bounds
          write(*,*) '--- occ orb. :', eps_occ
          write(*,*) '--- vir orb. :', eps_vir
        end if
        
!       get laplace points
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

        !create variables needed for loops
        ierr=talsh_tensor_construct(L1,C8,uov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(L2,C8,uvo_dims,init_val=ZERO)

        if (nvir<nchol) then
            dims4=nocc
            dims4(1:2)=nvir
            ierr=talsh_tensor_construct(temp,C8,dims4,init_val=ZERO)
        else
            dims4=nchol
            dims4(1:2)=nocc
            ierr=talsh_tensor_construct(R1,C8,dims4,init_val=ZERO)
            ierr=talsh_tensor_construct(R2,C8,dims4,init_val=ZERO)
            dims2=nchol
            ierr=talsh_tensor_construct(T1,C8,dims2,init_val=ZERO)
            ierr=talsh_tensor_construct(T2,C8,dims2,init_val=ZERO)
        end if

        mp2_energy=0

        do i=1,laplace%num

          ierr=talsh_tensor_init(L1,ZERO)
          ierr=talsh_tensor_add("R(i,a,k)+=O(i,a,k)",L1,uov_tensor)
          call laplace_mult(L1,laplace%val(i)/2.0d0,eps_occ_neg,eps_vir,zero_chol)

          ierr=talsh_tensor_init(L2,ZERO)
          ierr=talsh_tensor_add("R(a,i,k)+=O(a,i,k)",L2,uvo_tensor)
          call laplace_mult(L2,laplace%val(i)/2.0d0,eps_vir,eps_occ_neg,zero_chol)

          if (nvir<nchol) then
              !using nooc
              ierr=talsh_tensor_init(temp,ZERO)
              ierr=talsh_tensor_contract("D(a,b,i,j)+=T+(i,a,k)*V(b,j,k)",temp,L1,L2)
              ierr=talsh_tensor_contract("D(a,b,i,j)+=T+(i,b,k)*V(a,j,k)",temp,L1,L2,scale=MINUS_ONE)
              mp2_tens(1)=ZERO
              ierr=talsh_tensor_contract("E()+=T+(a,b,i,j)*V(a,b,i,j)",mp2_tensor,temp,temp)
          else
              !using nchol
              mp2_tens(1)=ZERO
              ierr=talsh_tensor_init(T1,ZERO)
              ierr=talsh_tensor_contract("M(k,n)+=T(i,a,k)*V+(i,a,n)",T1,L1,L1)
              ierr=talsh_tensor_init(T2,ZERO)
              ierr=talsh_tensor_contract("N(k,n)+=T+(b,j,k)*V(b,j,n)",T2,L2,L2)
              ierr=talsh_tensor_contract("E()+=M(k,n)*N(k,n)",mp2_tensor,T1,T2)

              ierr=talsh_tensor_init(R1,ZERO)
              ierr=talsh_tensor_contract("M(i,j,k,n)+=T+(b,j,k)*V+(i,b,n)",R1,L2,L1)
              ierr=talsh_tensor_init(R2,ZERO)
              ierr=talsh_tensor_contract("N(i,j,k,n)+=T(i,a,k)*V(a,j,n)",R2,L1,L2)
              ierr=talsh_tensor_contract("E()+=M(i,j,k,n)*N(i,j,k,n)",mp2_tensor,R1,R2,scale=MINUS_ONE)

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
            ierr=talsh_tensor_destruct(temp)
        else
            ierr=talsh_tensor_destruct(T1)
            ierr=talsh_tensor_destruct(T2)
            ierr=talsh_tensor_destruct(R1)
            ierr=talsh_tensor_destruct(R2)
        end if
        ierr=talsh_tensor_destruct(mp2_tensor)
        !cholesky
        ierr=talsh_tensor_destruct(uvo_tensor)
        ierr=talsh_tensor_destruct(uov_tensor)
        ierr=talsh_tensor_destruct(L1)
        ierr=talsh_tensor_destruct(L2)
        !laplace
        call deallocate_laplace(laplace)
        deallocate(eps_occ)
        deallocate(eps_vir)
        deallocate(zero_chol)

        !shutdown TALSH
        ierr=talsh_stats()
        ierr = talsh_shutdown()

        call print_date('Leaving MP2_LAP Module')
      end subroutine talsh_mp2lap_driver

      !-------------------------------------------------------------------------------------------------------------------------------------
      subroutine mp2lap_get_DM (exa_input)
!       This module is used to construct the MP2 density matrix by talsh 
!       Written by Johann Pototschnig

        implicit none
   
        type(exacc_input), intent(in) :: exa_input

        !variables to initlize library
        integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT):: host_arg_max

        integer(C_INT)      :: ierr, rank

        !cholesky vector
        integer, dimension(3)           :: uov_dims, uvo_dims
        type(talsh_tens_t)              :: uov_tensor, uvo_tensor

        !Laplace denominator
        integer                                :: nocc, nvir, nchol, i, i1
        real(8), allocatable                   :: eps_occ(:),eps_vir(:),eps_occ_neg(:), zero_chol(:)
        type(laplace_info)                     :: laplace
        integer(INTD)                          :: dims4(1:4), dims2(1:2), dims1(1)
        type(talsh_tens_t)                     :: L1p,L2p,L1q,L2q

        !intermediates
        type(talsh_tens_t)                     :: temp, temp1, mp2_dens

        !avoiding nvir
        type(talsh_tens_t)                     :: T1,R1,R2,S1,one_mp2lap

        call print_date("Entering MP2_LAP Module")

        !Initialize TALSH
        buf_size=exa_input%talsh_buff*buf_size
        ierr=talsh_init(buf_size,host_arg_max)
        call print_date('Initialized talsh library')
        write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

        call interest_initialize(.true.)
        call print_date('Initialized interest library')

        nocc=exa_input%nocc
        nvir=exa_input%nvir

        !Create MP2 density matrix
        dims2=nvir
        ierr=talsh_tensor_construct(mp2_dens,C8,dims2,init_val=ZERO)

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
        laplace%bounds(1) = 2.0d0*(eps_vir(1)-eps_occ(nocc))
        laplace%bounds(2) = 2.0d0*(eps_vir(nvir)-eps_occ(1))
        if (exa_input%print_level>10) then
          write(*,*) '--- Laplace bounds :', laplace%bounds
          write(*,*) '--- occ orb. :', eps_occ
          write(*,*) '--- vir orb. :', eps_vir
        end if
        
!       get laplace points
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

        !create variables needed for loops
        ierr=talsh_tensor_construct(L1p,C8,uov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(L2p,C8,uvo_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(L1q,C8,uov_dims,init_val=ZERO)
        ierr=talsh_tensor_construct(L2q,C8,uvo_dims,init_val=ZERO)

        if (nocc<nchol) then
            dims4=nocc
            dims4(1:2)=nvir
            ierr=talsh_tensor_construct(temp,C8,dims4,init_val=ZERO)
            ierr=talsh_tensor_construct(temp1,C8,dims4,init_val=ZERO)  

        else
            dims4=nchol
            dims4(1:2)=nvir
            ierr=talsh_tensor_construct(R1,C8,dims4,init_val=ZERO)
            ierr=talsh_tensor_construct(R2,C8,dims4,init_val=ZERO)
            dims2=nchol
            ierr=talsh_tensor_construct(T1,C8,dims2,init_val=ZERO)
            dims2=nvir
            ierr=talsh_tensor_construct(S1,C8,dims2,init_val=ZERO)
            dims1(1) = 1
            ierr=talsh_tensor_construct(one_mp2lap,C8,dims1(1:0),init_val=ONE)
        end if

        do i=1,laplace%num

          ierr=talsh_tensor_init(L1p,ZERO)
          ierr=talsh_tensor_add("R(i,a,k)+=O(i,a,k)",L1p,uov_tensor)
          call laplace_mult(L1p,laplace%val(i)/1.0d0,eps_occ_neg,eps_vir,zero_chol)

          ierr=talsh_tensor_init(L2p,ZERO)
          ierr=talsh_tensor_add("R(a,i,k)+=O(a,i,k)",L2p,uvo_tensor)
          call laplace_mult(L2p,laplace%val(i)/1.0d0,eps_vir,eps_occ_neg,zero_chol)

          if (nocc<nchol) then
              !using nooc
              ierr=talsh_tensor_init(temp,ZERO)
              ierr=talsh_tensor_contract("D(a,b,i,j)+=T+(i,a,k)*V(b,j,k)",temp,L1p,L2p)
              ierr=talsh_tensor_contract("D(a,b,i,j)+=T+(i,b,k)*V(a,j,k)",temp,L1p,L2p,scale=MINUS_ONE)
          end if

          ! ---------------- MP2 density matrix
          do i1=1,laplace%num

            ierr=talsh_tensor_init(L1q,ZERO)
            ierr=talsh_tensor_add("R(i,a,k)+=O(i,a,k)",L1q,uov_tensor)
            call laplace_mult(L1q,laplace%val(i1)/1.0d0,eps_occ_neg,eps_vir,zero_chol)
                
            ierr=talsh_tensor_init(L2q,ZERO)
            ierr=talsh_tensor_add("R(a,i,k)+=O(a,i,k)",L2q,uvo_tensor)
            call laplace_mult(L2q,laplace%val(i1)/1.0d0,eps_vir,eps_occ_neg,zero_chol)

            if (nocc<nchol) then
              ierr=talsh_tensor_init(temp1,ZERO)
              ierr=talsh_tensor_contract("D(a,b,i,j)+=T+(i,a,k)*V(b,j,k)",temp1,L1q,L2q)
              ierr=talsh_tensor_contract("D(a,b,i,j)+=T+(i,b,k)*V(a,j,k)",temp1,L1q,L2q,scale=MINUS_ONE)
                
              ierr=talsh_tensor_contract("E(a,b)+=T+(a,c,i,j)*V(b,c,i,j)",mp2_dens,temp,temp1, &
                                          scale=ONE_HALF*laplace%wei(i)*laplace%wei(i1))
            else
              !using nchol
              ierr=talsh_tensor_init(R2,ZERO)
              ierr=talsh_tensor_contract("N(b,a,k,n)+=T+(b,j,k)*V(a,j,n)",R2,L2p,L2q)
              ierr=talsh_tensor_init(T1,ZERO)
              ierr=talsh_tensor_contract("M(k,n)+=T(i,c,k)*V+(i,c,n)",T1,L1p,L1q)
              ierr=talsh_tensor_init(S1,ZERO)
              ierr=talsh_tensor_contract("E(a,b)+=M(k,n)*N(b,a,k,n)",S1,T1,R2, &
                                          scale=ONE_HALF*laplace%wei(i)*laplace%wei(i1))
              ierr=talsh_tensor_add("D(a,b)+=E(a,b)",mp2_dens,S1)
              ierr=talsh_tensor_contract("D(a,b)+=E+(b,a)*K()",mp2_dens,S1,one_mp2lap)
          
              ierr=talsh_tensor_init(R1,ZERO)
              ierr=talsh_tensor_contract("M(c,a,k,n)+=V(i,c,k)*T+(i,a,n)",R1,L1p,L1q)
              ierr=talsh_tensor_init(R2,ZERO)
              ierr=talsh_tensor_contract("N(b,c,k,n)+=T+(b,j,k)*V(c,j,n)",R2,L2p,L2q)
              ierr=talsh_tensor_init(S1,ZERO)
              ierr=talsh_tensor_contract("E(a,b)+=M(c,a,k,n)*N(b,c,k,n)",S1,R1,R2, &
                                          scale=MINUS_ONE_HALF*laplace%wei(i)*laplace%wei(i1))
              ierr=talsh_tensor_add("D(a,b)+=E(a,b)",mp2_dens,S1)
              ierr=talsh_tensor_contract("D(a,b)+=E+(b,a)*K()",mp2_dens,S1,one_mp2lap)

            end if
          
          end do

        end do

        call write_talsh_matrix(nvir,mp2_dens,'DM_complex')

        !clean up
        if (nocc<nchol) then
            ierr=talsh_tensor_destruct(temp)
            ierr=talsh_tensor_destruct(temp1)
        else
            ierr=talsh_tensor_destruct(T1)
            ierr=talsh_tensor_destruct(R1)
            ierr=talsh_tensor_destruct(R2)
        end if
        ierr=talsh_tensor_destruct(mp2_dens)
        !cholesky
        ierr=talsh_tensor_destruct(uvo_tensor)
        ierr=talsh_tensor_destruct(uov_tensor)
        ierr=talsh_tensor_destruct(L1p)
        ierr=talsh_tensor_destruct(L2p)
        ierr=talsh_tensor_destruct(L1q)
        ierr=talsh_tensor_destruct(L2q)
        ierr=talsh_tensor_destruct(one_mp2lap)
        !laplace
        call deallocate_laplace(laplace)
        deallocate(eps_occ)
        deallocate(eps_vir)
        deallocate(zero_chol)

        !shutdown TALSH
        ierr=talsh_stats()
        ierr = talsh_shutdown()

        call print_date('Leaving MP2_LAP Module')
      end subroutine mp2lap_get_DM


end module talsh_mp2lap
