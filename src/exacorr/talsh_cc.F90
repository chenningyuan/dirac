module talsh_cc
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra
        use talsh
        use exacorr_datatypes
        use exacorr_utils
        use intermediates
        use talsh_common_routines
        use, intrinsic:: ISO_C_BINDING

        implicit none
        complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),MINUS_ONE=(-1.D0,0.D0), &
                                 ONE_QUARTER=(0.25D0,0.D0), MINUS_ONE_QUARTER=(-0.25D0,0.D0), &
                                 ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
        real(8), parameter    :: ONE_QUARTER_R=0.25D0, ONE_HALF_R=0.5D0

        type(talsh_tens_t) :: one_tensor

        private

        public talsh_cc_driver

!        interface scale_with_denominators
!           module procedure scale_with_denominators
!        end interface
!
        interface diis
           module procedure diis
        end interface

        interface get_subtensor
           module procedure get_subtensor
        end interface

       contains
!-------------------------------------------------------------------------------------------
        subroutine talsh_cc_driver (exa_input)

!This subroutine drives the ccsd energy calculation
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code

!        Written by Lucas Visscher, winter 2016/2017

         use exacorr_global
         use talsh_gradient
         use talsh_cholesky_integrals

         type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer(C_INT)       :: ierr
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max

!        Common tensors
         type(talsh_comm_tens) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens) :: int_t
!        Lambda tensors
         type(talsh_tens_t) :: l1_tensor, l2_tensor
!        intermediate tensors (same notation as in relccsd)
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims, oovv_dims, vvoo_dims, vovo_dims
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor
         integer(C_INT)     :: one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CCSD,CC2,... control and result variables
         real(8) :: e_core, scf_energy, mp2_energy, cc_energy, t_energy(3)
         real(8) :: t_target_precision
         real(8) :: t1diag
         integer :: ncycles

         logical :: CCD, lambda
         logical :: CC2 

         integer :: iff        ! loop variable to loop over field strengths
         integer :: fock_mode  ! choice of Fock build options
         
!        Make some noise so that we know we are here
         call print_date('Entered cc_driver routine')

!        Initialize libraries
         buf_size=exa_input%talsh_buff*buf_size
         ierr=talsh_init(buf_size,host_arg_max)
         call print_date('Initialized talsh library')
         write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max
         call interest_initialize(.true.)
         call print_date('Initialized interest library')

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         ierr=talsh_tensor_construct(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=talsh_tensor_construct(result_tensor,C8,one_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,one_dims)

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir
         t_target_precision = exa_input%t_econv
         CCD = exa_input%ccd
         CC2 = exa_input%cc2     !L
         lambda = exa_input%lambda
         ncycles = exa_input%ncycles

         call print_exainput(exa_input,get_nao())

         call print_date("Getting orbital energies")
!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("Retrieved orbital energies")

!        Get fixed two-body tensors
         call print_date('-Start- Integral Transformation')
         if (exa_input%moint_scheme.eq.41) then
           call get_CC_integrals_chol_vec (nocc,nvir,mo_occ,mo_vir,int_t,exa_input%t_cholesky,exa_input%print_level)
         else if (exa_input%moint_scheme.eq.42) then
           call get_CC_integrals_chol_all (nocc,nvir,mo_occ,mo_vir,int_t,exa_input%t_cholesky,exa_input%print_level)
         else
           call get_CC_integrals(nocc,nvir,mo_occ,mo_vir,int_t,exa_input%moint_scheme,exa_input%print_level)
         end if
         call print_date('-End- Integral Transformation')

!        Get fixed one-body tensors (Fock matrices)
         ! build using orbital energies supplied by SCF program unless MRCONEE is present
         ! third option (for open shells) still to be implemented
         if ( exist_one_el(e_core) ) then
            fock_mode = 2
         else
            fock_mode = 1
         end if
         call init_fock_talsh (fock_mode,mo_occ,mo_vir,eps_occ,eps_vir,int_t,scf_energy,exa_input%print_level)

!        create tensors needed in CC
         oo_dims = nocc
         vv_dims = nvir
         ierr=talsh_tensor_construct(comm_t%foo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%fvv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%goo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%hoo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%gvv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%hvv,C8,vv_dims,init_val=ZERO)
         if (.not.CCD) then
            ov_dims(1) = nocc
            ov_dims(2) = nvir
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=talsh_tensor_construct(comm_t%fov,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(comm_t%hov,C8,ov_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(comm_t%t1,C8,vo_dims,init_val=ZERO)
         end if
         oooo_dims = nocc
         vovo_dims = nvir
         vovo_dims(2) = nocc
         vovo_dims(4) = nocc
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(comm_t%a_int,C8,oooo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%h_int,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%tau,C8,vvoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%t2,C8,vvoo_dims,init_val=ZERO)
         if (ierr.eq.TALSH_SUCCESS) then
            call print_date('Succesfully allocated all tensors')
         else
            write(*,*) " ccdriver: error in tensor allocation, code=",ierr
            return
         end if

         ! copy Fock matrix to comm, add finite field to the Fock matrix
         call set_fock(exa_input,mo_occ,mo_vir,int_t,comm_t,-1)
         ! define a diagonal zeroth order Hamiltonian with H^1 = H - H^0 and H^0_pq = epsilon_p delta_pq
         call define_H0(scf_energy,eps_occ,eps_vir,comm_t,int_t%oooo, &
                        exa_input%level_shift,-1,fock_mode,exa_input%nff(1))

!        Initialize t1 amplitudes to fov
         if (.not.CCD) then 
            ierr=talsh_tensor_init(comm_t%t1)
            ierr=talsh_tensor_contract("T(a,i)+=V+(i,a)",comm_t%t1,comm_t%fov,one_tensor)
         end if

!        Initialize t2 amplitudes with the MP2 values
         ierr=talsh_tensor_init(comm_t%t2)
         ierr=talsh_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",comm_t%t2,int_t%oovv,one_tensor)

         if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)
         else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
         end if

!        Calculate and print the MP2 energy
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv) ! o(2)v(2)
         mp2_energy=result_tens(1)*ONE_QUARTER_R
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,i)*V(i,a)",result_tensor,comm_t%t1,comm_t%fov)
         mp2_energy=mp2_energy+result_tens(1)
         write(*,*) ""
         write(*,*) "    MP2 energy = ", mp2_energy
         write(*,*) ""
!=====================
!        Calculate the pair correlation energy for each pair of occupied spin-orbits
         call pair_correlation_energies(comm_t%t2, int_t%oovv, nocc, nvir)
!=====================

         ! begin finite field loop (nff(2) gives number of field strengths, nff(1) gives number of fields)
         ! in case no finite field is specified nff(2)=1 and nff(1)=0
         do iff=1,exa_input%nff(2)
         ! copy Fock matrix to comm, add finite field to the Fock matrix
         call set_fock(exa_input,mo_occ,mo_vir,int_t,comm_t,iff)
         ! define a diagonal zeroth order Hamiltonian with H^1 = H - H^0 and H^0_pq = epsilon_p delta_pq
         call define_H0 (scf_energy,eps_occ,eps_vir,comm_t,int_t%oooo, &
                         exa_input%level_shift,exa_input%print_level,fock_mode,exa_input%nff(1))

!        Start CCSD iterations
         call print_date('-Start- CC Iterations')

         if (.not.CCD) then
            call diis (comm_t%t2,comm_t%t1,mode=1)
         else
            call diis (comm_t%t2,mode=1)
         end if

         if (.not.CC2) then
            call solve_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                cc_energy,CCD,ncycles,exa_input%print_level)
         else if (CC2) then
            call solve_CC2_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                cc_energy,CCD,ncycles,exa_input%print_level)
         end if
         call print_date('-End- CC Iterations')

         if (.not.CCD) then
            ! calculate T1 diagnostic
            t1diag = tensor_norm2(comm_t%t1)
            t1diag = t1diag/sqrt(dble(nocc))
            if (exa_input%do_triples .and. .not.CC2) then
              !Calculate (T) energy correction
              call print_date('-Start- Triples')
              call t3_energy (t_energy,eps_occ,eps_vir, comm_t%t1,comm_t%t2,comm_t%fov, &
                              int_t%ooov,int_t%vovv,int_t%oovv,exa_input%print_level)
              if (exa_input%tripl_lap) then 
                call t3_energy_laplace (t_energy,comm_t,int_t,eps_occ,eps_vir,exa_input)
              end if
              call print_date('-End- Triples')
              call print_exaoutput(exa_input%talsh, CCD, CC2, scf_energy, mp2_energy, cc_energy, t1diag=t1diag, t_energy=t_energy)
            else
              call print_exaoutput(exa_input%talsh, CCD, CC2, scf_energy, mp2_energy, cc_energy, t1diag=t1diag)
            end if
         else
            call print_exaoutput(exa_input%talsh, CCD, CC2, scf_energy, mp2_energy, cc_energy)
         end if

         !Destruct diis tensors
         if (.not.CCD) then
            call diis (comm_t%t2,comm_t%t1,mode=2)
         else
            call diis (comm_t%t2,mode=2)
         end if

         ! end finite field loop
         end do

         if (lambda.and..not.CC2) then
          call print_date('-Start- Lambda Equations')
!         Construct lambda tensor
          if (.not.CCD) ierr=talsh_tensor_construct(l1_tensor,C8,ov_dims,init_val=ZERO)
          oovv_dims(1:2) = nocc
          oovv_dims(3:4) = nvir
          ierr=talsh_tensor_construct(l2_tensor,C8,oovv_dims,init_val=ZERO)

!         Initialize l1 and l2 with converged t amplitudes
          if (.not.CCD) ierr=talsh_tensor_contract("L(i,a)+=T+(a,i)*K()",l1_tensor,comm_t%t1,one_tensor)
          ierr=talsh_tensor_contract("L(i,j,a,b)+=T+(a,b,i,j)*K()",l2_tensor,comm_t%t2,one_tensor)

!         solve lambda-equations and caluclate density matrix in MO basis
          if (.not.CCD) then
            if (.not.CC2) then
                call solve_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir,eps_occ,eps_vir, &
                                    t_target_precision,ncycles,exa_input%print_level,l1_tensor)
             else  
                call solve_CC2_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir,eps_occ,eps_vir, &
                                    t_target_precision,ncycles,exa_input%print_level,l1_tensor)
            end if 
            call talsh_dm_driver(exa_input,comm_t%t2,l2_tensor,comm_t%t1,l1_tensor)
          else
            call solve_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir,eps_occ,eps_vir, &
                                        t_target_precision,ncycles,exa_input%print_level)
            call talsh_dm_driver(exa_input,comm_t%t2,l2_tensor)
          end if

          if (.not.CCD) ierr=talsh_tensor_destruct(l1_tensor)
          ierr=talsh_tensor_destruct(l2_tensor)
          call print_date('-End- Lambda Equations')
         else
            write(*,*) 'Lambda equations with CC2 wavefunctions not supported in this release' 
         end if

!        Destruct tensors and shut down library
         ierr=talsh_tensor_destruct(int_t%oooo)
         ierr=talsh_tensor_destruct(int_t%oovv)
         ierr=talsh_tensor_destruct(int_t%vovo)
         ierr=talsh_tensor_destruct(int_t%vvvv)
         ierr=talsh_tensor_destruct(int_t%ooov)
         ierr=talsh_tensor_destruct(int_t%vovv)
         ierr=talsh_tensor_destruct(comm_t%foo)
         ierr=talsh_tensor_destruct(comm_t%fov)
         ierr=talsh_tensor_destruct(comm_t%fvv)
         ierr=talsh_tensor_destruct(comm_t%tau)
         ierr=talsh_tensor_destruct(comm_t%hoo)
         if (.not.CCD) then
            ierr=talsh_tensor_destruct(comm_t%hov)
            ierr=talsh_tensor_destruct(comm_t%t1)
         end if
         ierr=talsh_tensor_destruct(comm_t%hvv)
         ierr=talsh_tensor_destruct(comm_t%goo)
         ierr=talsh_tensor_destruct(comm_t%gvv)
         ierr=talsh_tensor_destruct(comm_t%a_int)
         ierr=talsh_tensor_destruct(comm_t%h_int)
         ierr=talsh_tensor_destruct(result_tensor)
         ierr=talsh_tensor_destruct(one_tensor)
         ierr=talsh_tensor_destruct(comm_t%t2)

!        Get statistics
         ierr=talsh_stats()
         ierr = talsh_shutdown()
         deallocate (eps_occ)
         deallocate (eps_vir)
         deallocate (mo_occ)
         deallocate (mo_vir)

!        Free global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         write(*,*) ""
         call print_date('Leaving cc_driver routine')

        end subroutine talsh_cc_driver

        subroutine solve_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                ccsd_energy,CCD,ncycles,print_level)

!        Routine for solving the t amplitude equations

         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_comm_tens), intent(inout) :: comm_t
         type(talsh_intg_tens), intent(inout) :: int_t
!        intermediate tensors (same notation as in relccsd)
         !type(talsh_tens_t) :: hoo_int_tensor, hvv_int_tensor
         type(talsh_tens_t) :: tau2_tensor
         type(talsh_tens_t) :: cint_tensor
!        auxilliary tensors for t1 equation
         type(talsh_tens_t) :: hov_aux_tensor, hoo_aux_tensor
!        solution tensors
         type(talsh_tens_t) :: s1_tensor, s2_tensor
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor
         integer(INTD)      :: result_dims(1)
         !Tensor dimensions
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims
         integer(INTD), dimension(4) :: vvoo_dims, ovoo_dims
         !Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: iccsd
         integer, intent(in) :: ncycles
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         real(8), intent(out) :: ccsd_energy
         logical, intent(in)  :: CCD
         integer, intent(in)  :: print_level
         integer (INTD)       :: ierr

!        Initialize scalars that are to be used as tensors in contractions
         result_dims = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         if (.not.CCD) then
            oo_dims = nocc
            ierr=talsh_tensor_construct(hoo_aux_tensor,C8,oo_dims,init_val=ZERO)
            ov_dims(1) = nocc
            ov_dims(2) = nvir
            ierr=talsh_tensor_construct(hov_aux_tensor,C8,ov_dims,init_val=ZERO)
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=talsh_tensor_construct(s1_tensor,C8,vo_dims,init_val=ZERO)
         end if
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(s2_tensor,C8,vvoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
         ovoo_dims(1) = nocc
         ovoo_dims(2) = nvir
         ovoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(cint_tensor,C8,ovoo_dims,init_val=ZERO)

         do iccsd = 1, ncycles

          call get_tau(comm_t,CCD)

          if (print_level.gt.6) then
            if (.not.CCD) write(*,*) "t1    = ",tensor_norm2(comm_t%t1)
            write(*,*) "t2    = ",tensor_norm2(comm_t%t2)*ONE_QUARTER_R
            write(*,*) "Tau   = ",tensor_norm2(comm_t%tau)*ONE_QUARTER_R
          end if

          if (ierr.ne.TALSH_SUCCESS) then
             write(*,*) " ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

!NOTE: contractions marked with "validated" gave numerically significant contributions and could be validated by comparison
!to the reference implementation, the ones marked "moretests" were too small to decide whether it was correct or not.
!        calculate intermediaries

          ierr=talsh_tensor_init(s2_tensor)

!         The vvvv contribution
          ierr=talsh_tensor_contract("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,int_t%vvvv,comm_t%tau,scale=ONE_HALF) ! validated o(2)v(4)

!         vv ! validated
          call get_hvv(comm_t,int_t%oovv)
          call get_gvv(comm_t,int_t%vovv,CCD)

!         vovo ! moretests
          ierr=talsh_tensor_init(tau2_tensor)
          ierr=talsh_tensor_add("T(a,d,l,i)+=R(a,d,l,i)",tau2_tensor,comm_t%t2,scale=ONE_HALF)
          if (.not.CCD) ierr=talsh_tensor_contract("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,comm_t%t1,comm_t%t1)                       ! o(2)v(2)
          call get_h(comm_t,int_t,CCD)
          ierr=talsh_tensor_contract("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",comm_t%h_int,int_t%oovv,tau2_tensor)           ! o(3)v(3)

!         This term appears four times as we need to permute both i and j and a and b
          ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)

!         ovoo ! validated
          ierr=talsh_tensor_init(cint_tensor)
          ierr=talsh_tensor_contract("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,int_t%vovv,comm_t%tau,scale=MINUS_ONE_HALF) ! o(3)v(3)
          if (.not.CCD) then
            ierr=talsh_tensor_contract("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,int_t%vovo,comm_t%t1)                           ! o(3)v(2)
            ierr=talsh_tensor_contract("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,int_t%vovo,comm_t%t1,scale=MINUS_ONE)           ! o(3)v(2)
          end if
          ierr=talsh_tensor_contract("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,int_t%ooov,one_tensor)                            ! o(3)v(1)

!         oo ! validated
          call get_hoo(comm_t,int_t%oovv)
          call get_goo(comm_t,int_t%ooov,CCD)

!         ov ! validated
          if (.not.CCD) then
            call get_hov(comm_t,int_t%oovv)
            ierr=talsh_tensor_init(hov_aux_tensor)
            ierr=talsh_tensor_add("L(k,c)+=H(k,c)",hov_aux_tensor,comm_t%hov)                       ! copy common term
            ierr=talsh_tensor_add("H(k,c)+=F(k,c)",hov_aux_tensor,comm_t%fov,scale=MINUS_TWO)       ! subtract fock matrix twice
          end if

!         oooo ! validated
          call get_a(comm_t,int_t,CCD)

          if (ierr.ne.TALSH_SUCCESS) then
             write(*,*) " ccdriver: error in calculation of intermediates, code=",ierr
             return
          end if

          if (print_level.gt.8) then
             print*, "hoint  = ",tensor_norm2(comm_t%hoo)
             print*, "hvint  = ",tensor_norm2(comm_t%hvv)
             print*, "goint  = ",tensor_norm2(comm_t%goo)
             print*, "gvint  = ",tensor_norm2(comm_t%gvv)
             if (.not.CCD) print*, "hovint = ",tensor_norm2(comm_t%hov)
             print*, "a int  = ",tensor_norm2(comm_t%a_int)
             print*, "c int  = ",tensor_norm2(cint_tensor)
             print*, "h int  = ",tensor_norm2(comm_t%h_int)
          end if

!         This term appears twice as we need to permute i and j
          ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,comm_t%goo,comm_t%t2)  ! validated o(3)v(2)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,comm_t%goo,comm_t%t2)  ! validated o(3)v(2)
          if (.not.CCD) then
!           Again two times, also note the complex conjugation as we use vovv instead of vvvo
            ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,int_t%vovv,comm_t%t1)  ! validated o(2)v(3)
            ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,int_t%vovv,comm_t%t1,scale=MINUS_ONE)  ! validated o(2)v(3)
          end if

!         This term appears twice as we need to permute a and b
          ierr=talsh_tensor_contract("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! validated o(2)v(3)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! validated o(2)v(3)
          if (.not.CCD) then
!           Again two times, with the C intermediate (note: in RELCCSD this is done differently)
            ierr=talsh_tensor_contract("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
            ierr=talsh_tensor_contract("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,comm_t%t1) ! o(3)v(2)
          end if

!         The oooo contribution
          ierr=talsh_tensor_contract("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,comm_t%a_int,comm_t%tau,scale=ONE_HALF) !validated o(4)v(2)

!         ierr=talsh_tensor_add("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!         work-around because complex conjugation does not work in the above statement
          ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)    ! validated o(2)v(2)
!         end of work-around

!         calculate s1 tensor
          if (.not.CCD) then
            ierr=talsh_tensor_init(s1_tensor)
!           ierr=talsh_tensor_add("S(a,i)+=F+(i,a)",s1_tensor,fov_tensor)
!           work-around because complex conjugation does not work in the above statement
            ierr=talsh_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
!           end of work-around
            ierr=talsh_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,scale=ONE_HALF) !validated o(2)v(3)
            ierr=talsh_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                     !validated o(1)v(2)
            ierr=talsh_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                     !validated o(2)v(1)
            ierr=talsh_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)                     !validated o(2)v(2)
            ierr=talsh_tensor_init(hoo_aux_tensor)                                                                     !moretests
            ierr=talsh_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)                !moretests o(2)v(1)
            ierr=talsh_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)                     !moretests o(2)v(1)
            ierr=talsh_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,scale=ONE_HALF) !validated o(3)v(2)
            ierr=talsh_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,scale=MINUS_ONE)     !validated o(2)v(2)
          end if

          if (print_level.gt.6) then
            if (.not.CCD) write(*,*) "s1    = ",tensor_norm2(s1_tensor)
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)*ONE_QUARTER_R
          end if

          if (ierr.ne.TALSH_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

!         prepare for next iteration
          if (.not.CCD) then
            ierr=talsh_tensor_init(comm_t%t1)
            ierr=talsh_tensor_add("T1(a,i)+=S1(a,i)",comm_t%t1,s1_tensor)
          end if
          ierr=talsh_tensor_init(comm_t%t2)
          ierr=talsh_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",comm_t%t2,s2_tensor)

          if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
          end if

          if (.not.CCD) then
            call diis (comm_t%t2,comm_t%t1,t_convergence=t_convergence)
          else
            call diis (comm_t%t2,t_convergence=t_convergence)
          end if

          if (.not.CCD) then
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("E()+=T(a,i)*F(i,a)",result_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
            ccsd_energy=result_tens(1)
          else
            ccsd_energy = 0.D0
          end if

!         Update tau
          call get_tau(comm_t,CCD)

          ierr=talsh_tensor_init(result_tensor)
          ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER_R

          call print_iteration(iccsd,t_convergence,ccsd_energy,print_level)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         write(*,*) ""

         if (.not.CCD) ierr=talsh_tensor_destruct(s1_tensor)
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(tau2_tensor)
         ierr=talsh_tensor_destruct(cint_tensor)
         if (.not.CCD) then
            ierr=talsh_tensor_destruct(hoo_aux_tensor)
            ierr=talsh_tensor_destruct(hov_aux_tensor)
         end if
         ierr=talsh_tensor_destruct(result_tensor)

        end subroutine

        subroutine solve_CC2_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                ccsd_energy,CCD,ncycles,print_level)

!        Routine for solving the t amplitude equations


         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_comm_tens), intent(inout) :: comm_t
         type(talsh_intg_tens), intent(inout) :: int_t
!        intermediate tensors (same notation as in relccsd)
         type(talsh_tens_t) :: Wabef_tensor,Wmnij_tensor     !L
         type(talsh_tens_t) :: Wabej_tensor,Wmbij_tensor     !L

! Add tensor usefull to avoid tauCC2_tensor 
         type(talsh_tens_t) :: Hovvv_tensor,Hvovv_tensor     !L
         type(talsh_tens_t) :: Hvvov_tensor,Hvvvo_tensor     !L
         type(talsh_tens_t) :: Hooov_tensor,Hoovo_tensor
         type(talsh_tens_t) :: Hvooo_tensor,Hovoo_tensor
         
         type(talsh_tens_t) :: cint_ovoo_tensor
         type(talsh_tens_t) :: cint_vvvo_tensor
!        auxilliary tensors for t1 equation
         type(talsh_tens_t) :: hov_aux_tensor, hoo_aux_tensor
!        solution tensors
         type(talsh_tens_t) :: s1_tensor, s2_tensor
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor
         integer            :: result_dims(1)
         !Tensor dimensions
         integer, dimension(2) :: oo_dims, ov_dims, vo_dims
         integer, dimension(4) :: vvoo_dims, ovoo_dims

         integer, dimension(4) :: vvvv_dims, oooo_dims    !L
         integer, dimension(4) :: vvvo_dims    !L

! new dims for Hxxxx 
         integer, dimension(4) :: vovv_dims, ovvv_dims, vvov_dims !, vvov_dims
         integer, dimension(4) :: ooov_dims, oovo_dims, vooo_dims !, vvov_dims !, vvov_dims

         logical :: debug_local = .false.

         !Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: iccsd
         integer, intent(in) :: ncycles,print_level
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         real(8), intent(out) :: ccsd_energy
         logical, intent(in) :: CCD
         integer :: ierr

         if (print_level.gt.8) debug_local = .true.

!        Initialize scalars that are to be used as tensors in contractions
         result_dims = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,0,DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         
         oo_dims = nocc
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         vo_dims(1) = nvir
         vo_dims(2) = nocc

         if (.not.CCD) then
          ierr=talsh_tensor_construct(hoo_aux_tensor,C8,oo_dims,init_val=ZERO)
          ierr=talsh_tensor_construct(hov_aux_tensor,C8,ov_dims,init_val=ZERO)
         end if

         ierr=talsh_tensor_construct(s1_tensor,C8,vo_dims,init_val=ZERO)
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(s2_tensor,C8,vvoo_dims,init_val=ZERO)

         vvvv_dims(1:4) = nvir  !L
         oooo_dims(1:4) = nocc  !L
         ierr=talsh_tensor_construct(Wabef_tensor,C8,vvvv_dims,init_val=ZERO) !L
         ierr=talsh_tensor_construct(Wmnij_tensor,C8,oooo_dims,init_val=ZERO) !L

         vvvo_dims(1:3) = nvir !L
         vvvo_dims(4)   = nocc !L
         ierr=talsh_tensor_construct(cint_vvvo_tensor,C8,vvvo_dims,init_val=ZERO)

         ovoo_dims(1) = nocc
         ovoo_dims(2) = nvir
         ovoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(cint_ovoo_tensor,C8,ovoo_dims,init_val=ZERO)

         ierr=talsh_tensor_construct(Wabej_tensor,C8,vvvo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Wmbij_tensor,C8,ovoo_dims,init_val=ZERO)

! Creation of the 4 intermediaires tensor : Hvovv , Hovvv, Hvvov, Hvvvo
         vovv_dims(1:4) = nvir
         vovv_dims(2)   = nocc
         ovvv_dims(1:4) = nvir
         ovvv_dims(1)   = nocc
         vvov_dims(1:4) = nvir
         vvov_dims(3)   = nocc
         ierr=talsh_tensor_construct(Hovvv_tensor,C8,ovvv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Hvovv_tensor,C8,vovv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Hvvvo_tensor,C8,vvvo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Hvvov_tensor,C8,vvov_dims,init_val=ZERO)
         ooov_dims(1:4) = nocc
         ooov_dims(4)   = nvir
         oovo_dims(1:4) = nocc
         oovo_dims(3)   = nvir
         vooo_dims(1:4) = nocc
         vooo_dims(1)   = nvir
         ierr=talsh_tensor_construct(Hooov_tensor,C8,ooov_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Hoovo_tensor,C8,oovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Hovoo_tensor,C8,ovoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(Hvooo_tensor,C8,vooo_dims,init_val=ZERO)

! resorted integrals nedded below

         ierr=talsh_tensor_init(cint_vvvo_tensor)
         ierr=talsh_tensor_contract("C(b,c,a,i)+=V+(a,i,b,c)*K()",cint_vvvo_tensor,int_t%vovv,one_tensor) !Inversion, warning : no more + in +=

         ierr=talsh_tensor_init(cint_ovoo_tensor)
         ierr=talsh_tensor_contract("C(k,a,i,j)+=V+(i,j,k,a)*K()",cint_ovoo_tensor,int_t%ooov,one_tensor) !Inversion, warning : no more + in +=

         do iccsd = 1, ncycles


          ierr=talsh_tensor_init(s2_tensor)
          if (debug_local) print*, "norm s2 tensor at the start of CC2 t2 ",tensor_norm2(s2_tensor)

          if (iccsd==1) then
                 call get_tau(comm_t,CCD)
                 if (debug_local) print *, " call get_tau "
          end if

!  Start s2_tensor
!         Diagram 15 : s2_{ij}^{ab} += < ab || ij >
          ierr=talsh_tensor_contract("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)
          if (debug_local) print*, "norm s2 tensor after diagram 15 ",tensor_norm2(s2_tensor)

          ierr=talsh_tensor_init(Wabej_tensor,ZERO)
          ierr=talsh_tensor_init(Wmbij_tensor,ZERO)           


          ierr=talsh_tensor_add("W(a,b,e,j)+=C(a,b,e,j)",Wabej_tensor,cint_vvvo_tensor,scale=ONE)
          ierr=talsh_tensor_contract("W(a,b,e,j)+=T(a,m)*V(b,m,e,j)",Wabej_tensor,comm_t%t1,int_t%vovo,scale=ONE_HALF)
          ierr=talsh_tensor_contract("W(a,b,e,j)+=T(b,m)*V(a,m,e,j)",Wabej_tensor,comm_t%t1,int_t%vovo,scale=MINUS_ONE_HALF)
       
          ierr=talsh_tensor_contract("S(a,b,i,j)=T(e,i)*W(a,b,e,j)",s2_tensor,comm_t%t1,Wabej_tensor,scale=ONE) 
          ierr=talsh_tensor_contract("S(a,b,i,j)=T(e,j)*W(a,b,e,i)",s2_tensor,comm_t%t1,Wabej_tensor,scale=MINUS_ONE)
 
         
          ierr=talsh_tensor_add("W(m,b,i,j)+=C(m,b,i,j)",Wmbij_tensor,cint_ovoo_tensor,scale=ONE)
          ierr=talsh_tensor_contract("W(m,b,i,j)+=T(e,i)*V(b,m,e,j)",Wmbij_tensor,comm_t%t1,int_t%vovo,scale=MINUS_ONE_HALF)
          ierr=talsh_tensor_contract("W(m,b,i,j)+=T(e,j)*V(b,m,e,i)",Wmbij_tensor,comm_t%t1,int_t%vovo,scale=ONE_HALF)
       
          ierr=talsh_tensor_contract("S(a,b,i,j)=T(a,m)*W(m,b,i,j)",s2_tensor,comm_t%t1,Wmbij_tensor,scale=MINUS_ONE) 
          ierr=talsh_tensor_contract("S(a,b,i,j)=T(b,m)*W(m,a,i,j)",s2_tensor,comm_t%t1,Wmbij_tensor,scale=ONE)

!         additional CC2 W intermediates
!
!         a. W_{ij}^{mn} = V_{ij}^{mn} + P\_(ij)\sum_{e} t^e_i V_{ej}^{mn} + \frac{1}{4}P\_(ij)\sum_{ef} t^e_i t^f_j V_{ef}^{mn}
          ierr=talsh_tensor_init(Wmnij_tensor)
          ierr=talsh_tensor_contract("Wmnij(m,n,i,j)+=V(m,n,i,j)*K()",Wmnij_tensor,int_t%oooo,one_tensor)
          ierr=talsh_tensor_contract("Wmnij(m,n,i,j)+=V(m,n,i,e)*T(e,j)",Wmnij_tensor,int_t%ooov,comm_t%t1)
          ierr=talsh_tensor_contract("Wmnij(m,n,i,j)+=V(m,n,j,e)*T(e,i)",Wmnij_tensor,int_t%ooov,comm_t%t1,scale=MINUS_ONE)
          
          ierr=talsh_tensor_init(Hoovo_tensor)
          ierr=talsh_tensor_contract("Hoovo(m,n,e,j)+=V(m,n,e,f)*T(f,j)",& 
               &       Hoovo_tensor,int_t%oovv,comm_t%t1)
          ierr=talsh_tensor_contract("Wmnij(m,n,i,j)+=Hoovo(m,n,e,j)*T(e,i)",&
               &       Wmnij_tensor,Hoovo_tensor,comm_t%t1,scale=ONE_QUARTER)
          ierr=talsh_tensor_init(Hooov_tensor)
          ierr=talsh_tensor_contract("Hooov(m,n,j,f)+=V(m,n,e,f)*T(e,j)",& 
               &       Hooov_tensor,int_t%oovv,comm_t%t1)
          ierr=talsh_tensor_contract("Wmnij(m,n,i,j)+=Hooov(m,n,j,f)*T(f,i)",&
               &       Wmnij_tensor,Hooov_tensor,comm_t%t1,scale=MINUS_ONE_QUARTER)

!         b. W_{ij}^{ab} = V_{ef}^{ab} - P\_(ab)\sum_{m} t_m^b V_{ef}^{am} + \frac{1}{4}P\_(ab)\sum_{ef} t^e_i t^f_j V_{ef}^{mn}
          ierr=talsh_tensor_init(Wabef_tensor)
          ierr=talsh_tensor_contract("Wabef(a,b,e,f)+=V(a,b,e,f)*K()",Wabef_tensor,int_t%vvvv,one_tensor)
          ierr=talsh_tensor_contract("Wabef(a,b,e,f)+=V(a,m,e,f)*T(b,m)",Wabef_tensor,int_t%vovv,comm_t%t1,scale=MINUS_ONE)
          ierr=talsh_tensor_contract("Wabef(a,b,e,f)+=V(b,m,e,f)*T(a,m)",Wabef_tensor,int_t%vovv,comm_t%t1)

          ierr=talsh_tensor_init(Hvovv_tensor)
          ierr=talsh_tensor_init(Hovvv_tensor)
          ierr=talsh_tensor_contract("Hvovv(a,n,e,f)+=V(m,n,e,f)*T(a,m)",& 
               & Hvovv_tensor,int_t%oovv,comm_t%t1,scale=ONE)        !Operation cost : v^3*o^2
          ierr=talsh_tensor_contract("Wabef(a,b,e,f)+=Hvovv(a,n,e,f)*T(b,n)",& 
               & Wabef_tensor,Hvovv_tensor,comm_t%t1,scale=ONE_QUARTER)  !Operation cost : v^4*o
          ierr=talsh_tensor_contract("Hovvv(m,a,e,f)+=V(m,n,e,f)*T(a,n)",& 
               & Hovvv_tensor,int_t%oovv,comm_t%t1,scale=ONE)        ! Operation cost :v^3*o^2
          ierr=talsh_tensor_contract("Wabef(a,b,e,f)+=Hovvv(m,a,e,f)*T(b,m)",& 
               & Wabef_tensor,Hovvv_tensor,comm_t%t1,scale=MINUS_ONE_QUARTER)  ! Operation cost : v^4*o
          if (debug_local) print*, "norm Wabefi_CC2 new",tensor_norm2(Wabef_tensor)


!         Diagram 21 : s2_{ij}^{ab} += \frac{1}{2}\sum_{mn} tauCC2^{mn}_{ab} W_{ij}^{mn}

          ierr=talsh_tensor_init(Hvooo_tensor)
          ierr=talsh_tensor_contract("H(a,n,i,j)+=Wmnij(m,n,i,j)*T(a,m)",&
               & Hvooo_tensor,Wmnij_tensor,comm_t%t1,scale=ONE)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=H(a,n,i,j)*T(b,n)",&
               & s2_tensor,Hvooo_tensor,comm_t%t1,scale=ONE_HALF)
          ierr=talsh_tensor_init(Hovoo_tensor)
          ierr=talsh_tensor_contract("H(m,a,i,j)+=Wmnij(m,n,i,j)*T(a,n)",&
               & Hovoo_tensor,Wmnij_tensor,comm_t%t1,scale=ONE)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=H(m,a,i,j)*T(b,m)",&
               & s2_tensor,Hovoo_tensor,comm_t%t1,scale=MINUS_ONE_HALF)
          if (debug_local) print*, "norm s2 tensor after diagram 21 new ",tensor_norm2(s2_tensor)

!         Diagram 22 : s2_{ij}^{ab} += \frac{1}{2}\sum_{ef} tauCC2^{ij}_{ef} W_{ef}^{ab}  

          ierr=talsh_tensor_init(Hvvov_tensor)
          ierr=talsh_tensor_init(Hvvvo_tensor)
          ierr=talsh_tensor_contract("Hvvov(a,b,i,f)+=Wabef(a,b,e,f)*T(e,i)",& 
               & Hvvov_tensor,Wabef_tensor,comm_t%t1,scale=ONE)    ! Operationcost : v^4*o
          ierr=talsh_tensor_contract("S(a,b,i,j)+= Hvvov(a,b,i,f)*T(f,j)",& 
               & s2_tensor,Hvvov_tensor,comm_t%t1,scale=ONE_HALF)    !  Operation cost : v^3*o^2
          ierr=talsh_tensor_contract("Hvvvo(a,b,e,i)+=Wabef(a,b,e,f)*T(f,i)",& 
               & Hvvvo_tensor,Wabef_tensor,comm_t%t1,scale=ONE)                  !Operation cost : v^4*o
          ierr=talsh_tensor_contract("S(a,b,i,j)+= Hvvvo(a,b,e,i)*T(e,j)",& 
               & s2_tensor,Hvvvo_tensor,comm_t%t1,scale=MINUS_ONE_HALF)     ! Operation cost : v^3*o^2
          if (debug_local) print*, "norm s2 tensor after diagram 22 new ",tensor_norm2(s2_tensor)


!         Diagram 23 : s2_{ij}^{ab} +=  P\_(ab)\sum_{e} f_e^b t_{ij}^{ae} ! Use of the Fock matrix 
          ierr=talsh_tensor_contract("S(a,b,i,j)+=F(b,e)*T(a,e,i,j)",s2_tensor,comm_t%fvv,comm_t%t2)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=F(a,e)*T(b,e,i,j)",s2_tensor,comm_t%fvv,comm_t%t2,scale=MINUS_ONE)
          if (debug_local) print*, "norm s2 tensor after diagram 23 ",tensor_norm2(s2_tensor)

!         Diagram 24 : s2_{ij}^{ab} += -P\_(ij)\sum_{k} f_j^k t_{ik}^{ab} ! Use of the Fock matrix
          ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,j)*T(a,b,i,m)",s2_tensor,comm_t%foo,comm_t%t2,scale=MINUS_ONE)
          ierr=talsh_tensor_contract("S(a,b,i,j)+=F(m,i)*T(a,b,j,m)",s2_tensor,comm_t%foo,comm_t%t2)
          if (debug_local) print*, "norm s2 tensor after diagram 24 ",tensor_norm2(s2_tensor)

          if (debug_local) print*, "Norm s2 before s1",tensor_norm2(s2_tensor)

!        done with cc2 t2 equation, now cc2 t1 equation, which is exactly the same as the ccsd t1 equation

          ! Call F intermediates - used in \Ss ai
          call get_hoo(comm_t,int_t%oovv)
          call get_hvv(comm_t,int_t%oovv)
          call get_hov(comm_t,int_t%oovv) 
          if (.not.CCD) then
            ierr=talsh_tensor_init(hov_aux_tensor)
            ierr=talsh_tensor_add("L(k,c)+=H(k,c)",hov_aux_tensor,comm_t%hov)                       ! copy common term
            ierr=talsh_tensor_add("H(k,c)+=F(k,c)",hov_aux_tensor,comm_t%fov,scale=MINUS_TWO)       ! subtract fock matrix twice
            if (debug_local) print *, "In .not.CCD check for hov"
          end if

!         calculate s1 tensor
          ierr=talsh_tensor_init(s1_tensor)

!           Diagramme b1  :   \Ss ai= \f ai         
          ierr=talsh_tensor_contract("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
          if (debug_local) print*, "s1 - first step with \F ia ",tensor_norm2(s1_tensor)
!           
          ierr=talsh_tensor_contract("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,scale=ONE_HALF) !validated o(2)v(3)
          ierr=talsh_tensor_contract("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                     !validated o(1)v(2)
          ierr=talsh_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                     !validated o(2)v(1)
          ierr=talsh_tensor_contract("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)                     !validated o(2)v(2)
          ierr=talsh_tensor_init(hoo_aux_tensor)                                                                     !moretests
          ierr=talsh_tensor_contract("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)                !moretests o(2)v(1)
          ierr=talsh_tensor_contract("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)                     !moretests o(2)v(1)
          ierr=talsh_tensor_contract("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,scale=ONE_HALF) !validated o(3)v(2)
          ierr=talsh_tensor_contract("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,scale=MINUS_ONE)     !validated o(2)v(2)

          if (debug_local) print*, "Norm s1 before Energy",tensor_norm2(s1_tensor)
      
          if (ierr.ne.TALSH_SUCCESS) then
             print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

!         prepare for next iteration
          ierr=talsh_tensor_init(comm_t%t1)
          ierr=talsh_tensor_add("T1(a,i)+=S1(a,i)",comm_t%t1,s1_tensor)
          ierr=talsh_tensor_init(comm_t%t2)
          ierr=talsh_tensor_add("T2(a,b,i,j)+=S2(a,b,i,j)",comm_t%t2,s2_tensor)
          
          call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)

          call diis (comm_t%t2,comm_t%t1,t_convergence=t_convergence)

          ierr=talsh_tensor_init(result_tensor)
          ierr=talsh_tensor_contract("E()+=T(a,i)*F(i,a)",result_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
          ccsd_energy=result_tens(1)

!         Update tau
          call get_tau(comm_t,CCD)

          if (debug_local) print*, "Check for tau from get_tau", tensor_norm2(comm_t%tau)

          ierr=talsh_tensor_init(result_tensor)
          ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER_R

!         call print_iteration(iccsd,t_target_precision,t_convergence,ccsd_energy)
          call print_iteration(iccsd,t_convergence,ccsd_energy)
          if (t_convergence .lt. t_target_precision) exit

         end do


! Clean tensor
         ierr=talsh_tensor_destruct(s1_tensor)
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(cint_ovoo_tensor)
         ierr=talsh_tensor_destruct(cint_vvvo_tensor)

         ierr=talsh_tensor_destruct(Hvvov_tensor)
         ierr=talsh_tensor_destruct(Hvvvo_tensor)
         ierr=talsh_tensor_destruct(Hovvv_tensor)
         ierr=talsh_tensor_destruct(Hvovv_tensor)
         ierr=talsh_tensor_destruct(Hooov_tensor)
         ierr=talsh_tensor_destruct(Hoovo_tensor)
         ierr=talsh_tensor_destruct(Hovoo_tensor)
         ierr=talsh_tensor_destruct(Hvooo_tensor)

         ierr=talsh_tensor_destruct(Wabef_tensor)
         ierr=talsh_tensor_destruct(Wmnij_tensor) 

         ierr=talsh_tensor_destruct(Wabej_tensor)
         ierr=talsh_tensor_destruct(Wmbij_tensor)

         if (.not.CCD) then
            ierr=talsh_tensor_destruct(hoo_aux_tensor)
            ierr=talsh_tensor_destruct(hov_aux_tensor)
         end if
         ierr=talsh_tensor_destruct(result_tensor)

        end subroutine

        subroutine solve_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir,eps_occ,eps_vir, &
                                          t_target_precision,ncycles,print_level,l1_tensor)
!        Routine for solving the lambda equations

         integer, intent(in) :: nocc, nvir, print_level
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_comm_tens), intent(inout) :: comm_t
         type(talsh_intg_tens), intent(inout) :: int_t
         type(talsh_lambda_tens)              :: lambda_t
!        solution tensors
         type(talsh_tens_t), intent(inout)           :: l2_tensor
         type(talsh_tens_t), intent(inout), optional :: l1_tensor
         type(talsh_tens_t) :: s1_tensor, s2_tensor
!        Tensor dimensions
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims,ovoo_dims,ooov_dims,vovo_dims,oovv_dims, &
                                  vvvo_dims,vovv_dims,vvvv_dims
!        Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: ilambda
         integer, intent(in) :: ncycles
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence

         logical :: CCD = .false.

         integer(INTD) :: ierr

         write(*,*) ""
         write(*,'(A)') "*******************************************************************************"
         write(*,'(A)') "********************** Solving Left Eigenvector Equation **********************"
         write(*,'(A)') "*******************************************************************************"
         write(*,*) ""

         if (.not.present(l1_tensor)) CCD = .true.

!***************************************
!        construct fixed intermediates *
!***************************************

         oo_dims = nocc
         ierr=talsh_tensor_construct(lambda_t%fbar_oo,C8,oo_dims,init_val=ZERO)
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         ierr=talsh_tensor_construct(lambda_t%fbar_ov,C8,ov_dims,init_val=ZERO)
         vv_dims = nvir
         ierr=talsh_tensor_construct(lambda_t%fbar_vv,C8,vv_dims,init_val=ZERO)
         oooo_dims = nocc
         ierr=talsh_tensor_construct(lambda_t%w_oooo,C8,oooo_dims,init_val=ZERO)
         ovoo_dims(1)   = nocc
         ovoo_dims(2)   = nvir
         ovoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(lambda_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=talsh_tensor_construct(lambda_t%w_ooov,C8,ooov_dims,init_val=ZERO)
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=talsh_tensor_construct(lambda_t%w_vovo,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(lambda_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)
         vvvo_dims(1:3) = nvir
         vvvo_dims(4)   = nocc
         ierr=talsh_tensor_construct(lambda_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(lambda_t%w_vovv,C8,vovv_dims,init_val=ZERO)
         vvvv_dims = nvir
         ierr=talsh_tensor_construct(lambda_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)

!******************************************
!        obtain intermediates from CC(S)D *
!******************************************

         call get_tau(comm_t,CCD)
         call get_hvv(comm_t,int_t%oovv)
         call get_gvv(comm_t,int_t%vovv,CCD)
         call get_h(comm_t,int_t,CCD)
         call get_hoo(comm_t,int_t%oovv)
         call get_goo(comm_t,int_t%ooov,CCD)
         call get_hov(comm_t,int_t%oovv)
         call get_a(comm_t,int_t,CCD)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t)
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD)
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,nocc,nvir,CCD)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD)
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD)
         call print_date('Finished calculating fixed intermediates')

         if (print_level.gt.8) then
             if (.not.CCD) print*, "T1amp      = ", tensor_norm2(comm_t%t1)
             print*, "T2amp      = ", tensor_norm2(comm_t%t2)*ONE_QUARTER_R
             ! print*, "Tau        = ", tensor_norm2(comm_t%tau)*ONE_QUARTER_R
             ! print*, "Fbar_oo    = ", tensor_norm2(lambda_t%fbar_oo)
             print*, "Fbar_ov    = ", tensor_norm2(lambda_t%fbar_ov)
             ! print*, "Fbar_vv    = ", tensor_norm2(lambda_t%fbar_vv)
             ! print*, "w_vvvv     = ", tensor_norm2(lambda_t%w_vvvv)*ONE_QUARTER_R
             print*, "wbar_vovo  = ", tensor_norm2(lambda_t%wbar_vovo)
             print*, "w_vvvo     = ", tensor_norm2(lambda_t%w_vvvo)*0.5D0
             print*, "w_vovv     = ", tensor_norm2(lambda_t%w_vovv)*0.5D0
             print*, "w_vovo     = ", tensor_norm2(lambda_t%w_vovo)
             print*, "w_ovoo     = ", tensor_norm2(lambda_t%w_ovoo)*0.5D0
             print*, "w_ooov     = ", tensor_norm2(lambda_t%w_ooov)*0.5D0
             print*, "w_oooo     = ", tensor_norm2(lambda_t%w_oooo)*ONE_QUARTER_R
         end if

!************************************
!        construct solution tensors *
!************************************

         if (.not.CCD) ierr=talsh_tensor_construct(s1_tensor,C8,ov_dims,init_val=ZERO)
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(s2_tensor,C8,oovv_dims,init_val=ZERO)

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=talsh_tensor_construct(lambda_t%goo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(lambda_t%gvv,C8,vv_dims,init_val=ZERO)

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (l2_tensor,l1_tensor,mode=1)
         else
           call diis (l2_tensor,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)

!         debug
          if (print_level.gt.8) then
             if (.not.CCD) print*, "L1  = ", tensor_norm2(l1_tensor)
             print*, "L2  = ", tensor_norm2(l2_tensor)*ONE_QUARTER_R
             print*, "goo = ", tensor_norm2(lambda_t%goo)
             print*, "gvv = ", tensor_norm2(lambda_t%gvv)
          end if

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=talsh_tensor_init(s1_tensor)

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=talsh_tensor_add("S(i,a)+=F(i,a)",s1_tensor,lambda_t%fbar_ov)

!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,scale=MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,scale=MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,scale=ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,scale=ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=talsh_tensor_init(s2_tensor)

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=talsh_tensor_add("S(i,j,a,b)+=V(i,j,a,b)",s2_tensor,int_t%oovv)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |
!----------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t%fbar_vv)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |
!----------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |
!------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t%w_oooo,scale=ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |
!----------------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t%w_vovo)

!-------------------------------------------------------
!         term 6: S(ij,ab) += P(a,b) V(ij,ae) * G(e,b) |
!-------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,lambda_t%gvv)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,lambda_t%gvv)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(,j,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)
          end if

!-------------------------------------------------------
!         term 8: S(ij,ab) -= P(i,j) V(im,ab) * G(j,m) |
!-------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,lambda_t%goo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,lambda_t%goo)

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |
!-----------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t%w_vvvv,scale=ONE_HALF)

          if (ierr.ne.TALSH_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of lambda equations, code=",ierr
             return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=talsh_tensor_init(l1_tensor)
            ierr=talsh_tensor_add("L1(i,a)+=S1(i,a)",l1_tensor,s1_tensor)
          end if
          ierr=talsh_tensor_init(l2_tensor)
          ierr=talsh_tensor_add("L2(i,j,a,b)+=S2(i,j,a,b)",l2_tensor,s2_tensor)
          if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,l1_tensor,l2_tensor)
            call diis (l2_tensor,l1_tensor,t_convergence=t_convergence)
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=l2_tensor)
            call diis (l2_tensor,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", tensor_norm2(l1_tensor)
         write(*,'(A,ES22.16)') "  Final L2amp = ", tensor_norm2(l2_tensor)*ONE_QUARTER_R
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=talsh_tensor_destruct(s1_tensor)
           call diis (l2_tensor,l1_tensor,mode=2)
         else
           call diis (l2_tensor,mode=2)
         end if
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(lambda_t%fbar_oo)
         ierr=talsh_tensor_destruct(lambda_t%fbar_ov)
         ierr=talsh_tensor_destruct(lambda_t%fbar_vv)
         ierr=talsh_tensor_destruct(lambda_t%w_oooo)
         ierr=talsh_tensor_destruct(lambda_t%w_ooov)
         ierr=talsh_tensor_destruct(lambda_t%w_ovoo)
         ierr=talsh_tensor_destruct(lambda_t%w_vovo)
         ierr=talsh_tensor_destruct(lambda_t%w_vvvo)
         ierr=talsh_tensor_destruct(lambda_t%w_vovv)
         ierr=talsh_tensor_destruct(lambda_t%w_vvvv)
         ierr=talsh_tensor_destruct(lambda_t%wbar_vovo)
         ierr=talsh_tensor_destruct(lambda_t%goo)
         ierr=talsh_tensor_destruct(lambda_t%gvv)

        end subroutine solve_lambda_equations


        subroutine solve_CC2_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir, &
                                     eps_occ,eps_vir,t_target_precision,ncycles,print_level,l1_tensor)

!        Routine for solving the lambda equations

         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_comm_tens), intent(inout) :: comm_t
         type(talsh_comm_tens) :: comm_t_CC2
         type(talsh_intg_tens), intent(inout) :: int_t
         type(talsh_lambda_tens)              :: lambda_t,lambda_t_CC2
!        solution tensors
         type(talsh_tens_t), intent(inout)           :: l2_tensor
         type(talsh_tens_t), intent(inout), optional :: l1_tensor
         type(talsh_tens_t) :: s1_tensor, s2_tensor
!        Tensor dimensions
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vv_dims, vo_dims
         integer(INTD), dimension(4) :: oooo_dims,ovoo_dims,ooov_dims,vovo_dims,oovv_dims, &
                                  vvvo_dims,vovv_dims,vvvv_dims,vvoo_dims
!        Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: ilambda

         integer, intent(in) :: ncycles
         integer, intent(in) :: print_level
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence

         logical :: CCD = .false.

         integer(INTD) :: ierr

         write(*,*) ""
         write(*,'(A)') "*******************************************************************************"
         write(*,'(A)') "********************** Solving Left Eigenvector Equation **********************"
         write(*,'(A)') "*******************************************************************************"
         write(*,*) ""

         if (.not.present(l1_tensor)) CCD = .true.

!***************************************
!        construct fixed intermediates *
!***************************************

         oo_dims = nocc
         ierr=talsh_tensor_construct(lambda_t%fbar_oo,C8,oo_dims,init_val=ZERO)
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         ierr=talsh_tensor_construct(lambda_t%fbar_ov,C8,ov_dims,init_val=ZERO)
         vv_dims = nvir
         ierr=talsh_tensor_construct(lambda_t%fbar_vv,C8,vv_dims,init_val=ZERO)
         oooo_dims = nocc
         ierr=talsh_tensor_construct(lambda_t%w_oooo,C8,oooo_dims,init_val=ZERO)
         ovoo_dims(1)   = nocc
         ovoo_dims(2)   = nvir
         ovoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(lambda_t%w_ovoo,C8,ovoo_dims,init_val=ZERO)
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=talsh_tensor_construct(lambda_t%w_ooov,C8,ooov_dims,init_val=ZERO)
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=talsh_tensor_construct(lambda_t%w_vovo,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(lambda_t%wbar_vovo,C8,vovo_dims,init_val=ZERO)
         vvvo_dims(1:3) = nvir
         vvvo_dims(4)   = nocc
         ierr=talsh_tensor_construct(lambda_t%w_vvvo,C8,vvvo_dims,init_val=ZERO)
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(lambda_t%w_vovv,C8,vovv_dims,init_val=ZERO)
         vvvv_dims = nvir
         ierr=talsh_tensor_construct(lambda_t%w_vvvv,C8,vvvv_dims,init_val=ZERO)

!***
! List of NEW INTERMEDIATES for the CC2 Method
!***  
         ierr=talsh_tensor_construct(lambda_t_CC2%fbar_oo,C8,oo_dims,init_val=ZERO)   !cUSefull and different
         ierr=talsh_tensor_construct(lambda_t_CC2%fbar_ov,C8,ov_dims,init_val=ZERO)  ! Usefull but same as lambda1 and no difference with CCSD
         ierr=talsh_tensor_construct(lambda_t_CC2%fbar_vv,C8,vv_dims,init_val=ZERO)   ! Usefull and different
         ierr=talsh_tensor_construct(lambda_t_CC2%w_oooo,C8,oooo_dims,init_val=ZERO)   ! Usefull and different
!!!!!!  Inportant here :  The heart of w_vovo intermediate : h_int is the same as in 
!!!!!!          CCSD and for \Lambda_1 either :
         ierr=talsh_tensor_construct(lambda_t_CC2%w_vovo,C8,vovo_dims,init_val=ZERO) ! Usefull in CC2 - Same part as h_int
         ierr=talsh_tensor_construct(lambda_t_CC2%w_vvvv,C8,vvvv_dims,init_val=ZERO)   ! Usefull and different

!*
! End of list
!*
!******************************************
!        obtain intermediates from CC(S)D *
!******************************************

         call get_tau(comm_t,CCD)
         call get_hvv(comm_t,int_t%oovv)
         call get_gvv(comm_t,int_t%vovv,CCD)
         call get_h(comm_t,int_t,CCD)
         call get_hoo(comm_t,int_t%oovv)
         call get_goo(comm_t,int_t%ooov,CCD)
         call get_hov(comm_t,int_t%oovv)
         call get_a(comm_t,int_t,CCD)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t)
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD)
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,nocc,nvir,CCD)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD)
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD)
         call print_date('Finished calculating fixed intermediates')


!****************************************
!        get CC2 Specific Intermediates *
!****************************************

         vo_dims(1)=nvir
         vo_dims(2)=nocc
         vvoo_dims(1:2)=nvir
         vvoo_dims(3:4)=nocc
         ierr=talsh_tensor_construct(comm_t_CC2%t1,C8,vo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t_CC2%t2,C8,vvoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t_CC2%tau,c8,vvoo_dims,init_val=zero)
         ierr=talsh_tensor_construct(comm_t_CC2%h_int,c8,vovo_dims,init_val=zero)

         ierr=talsh_tensor_add("H(a,i)+=T(a,i)",comm_t_CC2%t1,comm_t%t1,scale=ONE)
         ierr=talsh_tensor_add("H(a,b,i,j)+=T(a,b,i,j)",comm_t_CC2%t2,comm_t%t2,scale=ONE)
         ierr=talsh_tensor_add("H(a,i)+=F(a,i)",comm_t_CC2%fov,comm_t%fov,scale=ONE)
         ierr=talsh_tensor_add("H(a,b)+=F(a,b)",comm_t_CC2%fvv,comm_t%fvv,scale=ONE)
         ierr=talsh_tensor_add("H(i,j)+=F(i,j)",comm_t_CC2%foo,comm_t%foo,scale=ONE)
         ierr=talsh_tensor_add("W(b,m,e,j)+=H(b,m,e,j)",comm_t_CC2%h_int,comm_t%h_int,scale=ONE)

         call get_tauCC2(comm_t_CC2)

         call get_fbar_vv_CC2(lambda_t_CC2%fbar_vv,comm_t_CC2,int_t) ! for CCSD : fbar_vv = gvv  (L! : change ?)
         call get_fbar_oo_CC2(lambda_t_CC2%fbar_oo,comm_t_CC2,int_t)
         call get_fbar_ov_CC2(lambda_t_CC2%fbar_ov,comm_t_CC2,int_t)
  
         call get_w_vovo_CC2(lambda_t_CC2%w_vovo,comm_t_CC2,int_t%oovv,nocc,nvir)
         call get_w_oooo_CC2(lambda_t_CC2%w_oooo,comm_t_CC2,int_t) 
         call get_w_vvvv_CC2(lambda_t_CC2%w_vvvv,comm_t_CC2,int_t) 

!************************************
!        construct solution tensors *
!************************************

         if (.not.CCD) ierr=talsh_tensor_construct(s1_tensor,C8,ov_dims,init_val=ZERO)
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(s2_tensor,C8,oovv_dims,init_val=ZERO)

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=talsh_tensor_construct(lambda_t%goo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(lambda_t%gvv,C8,vv_dims,init_val=ZERO)

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (l2_tensor,l1_tensor,mode=1)
         else
           call diis (l2_tensor,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=talsh_tensor_init(s1_tensor)

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=talsh_tensor_add("S(i,a)+=F(i,a)",s1_tensor,lambda_t%fbar_ov)

!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,l1_tensor,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,l1_tensor,lambda_t%fbar_oo,scale=MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,l1_tensor,lambda_t%w_vovo,scale=MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,l2_tensor,lambda_t%w_ovoo,scale=ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,l2_tensor,lambda_t%w_vvvo,scale=ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=talsh_tensor_contract("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=talsh_tensor_init(s2_tensor)

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=talsh_tensor_add("S(i,j,a,b)+=V(i,j,a,b)",s2_tensor,int_t%oovv)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |      Fbar_vv - CC2 
!----------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_vv)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |      Fbar_oo - CC2 
!----------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_oo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,l2_tensor,lambda_t_CC2%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |          W_oooo - CC2 
!------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,l2_tensor,lambda_t_CC2%w_oooo,scale=ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |       W vovo - CC2 
!----------------------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,l2_tensor,lambda_t_CC2%w_vovo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,l2_tensor,lambda_t_CC2%w_vovo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,l2_tensor,lambda_t_CC2%w_vovo)
          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,l2_tensor,lambda_t_CC2%w_vovo)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(ij,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,l1_tensor,lambda_t%w_ooov)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,l1_tensor,lambda_t%w_ooov)
          end if

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |          
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,l1_tensor,lambda_t%w_vovv)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,l1_tensor,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,l1_tensor,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=talsh_tensor_contract("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,l1_tensor,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |      W vvvv - CC2 
!-----------------------------------------------------

          ierr=talsh_tensor_contract("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,l2_tensor,lambda_t_CC2%w_vvvv,scale=ONE_HALF)

          if (ierr.ne.TALSH_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of lambda equations, code=",ierr
             return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=talsh_tensor_init(l1_tensor)
            ierr=talsh_tensor_add("L1(i,a)+=S1(i,a)",l1_tensor,s1_tensor)
          end if
          ierr=talsh_tensor_init(l2_tensor)
          ierr=talsh_tensor_add("L2(i,j,a,b)+=S2(i,j,a,b)",l2_tensor,s2_tensor)
          if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,l1_tensor,l2_tensor)
            call diis (l2_tensor,l1_tensor,t_convergence=t_convergence)
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=l2_tensor)
            call diis (l2_tensor,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", tensor_norm2(l1_tensor)
         write(*,'(A,ES22.16)') "  Final L2amp = ", tensor_norm2(l2_tensor)*ONE_QUARTER_R
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=talsh_tensor_destruct(s1_tensor)
           call diis (l2_tensor,l1_tensor,mode=2)
         else
           call diis (l2_tensor,mode=2)
         end if
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(lambda_t%fbar_oo)
         ierr=talsh_tensor_destruct(lambda_t%fbar_ov)
         ierr=talsh_tensor_destruct(lambda_t%fbar_vv)
         ierr=talsh_tensor_destruct(lambda_t%w_oooo)
         ierr=talsh_tensor_destruct(lambda_t%w_ooov)
         ierr=talsh_tensor_destruct(lambda_t%w_ovoo)
         ierr=talsh_tensor_destruct(lambda_t%w_vovo)
         ierr=talsh_tensor_destruct(lambda_t%w_vvvo)
         ierr=talsh_tensor_destruct(lambda_t%w_vovv)
         ierr=talsh_tensor_destruct(lambda_t%w_vvvv)
         ierr=talsh_tensor_destruct(lambda_t%wbar_vovo)
         ierr=talsh_tensor_destruct(lambda_t%goo)
         ierr=talsh_tensor_destruct(lambda_t%gvv)

        end subroutine solve_CC2_lambda_equations 


        subroutine init_fock_talsh (mode,mo_occ,mo_vir,eps_occ,eps_vir,int_t,scf_energy,print_level)

!        Written by Lucas Visscher, summer 2017
!        Constructs unrestricted Hartree-Fock matrix and computes corresponding UHF energy.

         use exacorr_global

         implicit none

         integer, intent(in)                  :: mode ! 1: Use precomputed orbital energies and HF energy from start file
                                                      ! 2: Build from scratch using one- and two-electron integrals
                                                      ! 3: Construct as F^{UHF} = F^{AvHF} + (F^{UHF} - F^{AvHF})
                                                      !    where F^{AvHF} is the average HF expression used in DIRAC
         integer, intent(in)                  :: mo_occ(:),mo_vir(:)   ! orbital indices
         real(8), intent(in)                  :: eps_occ(:),eps_vir(:) ! orbital energies (as read from start file)
         real(8), intent(out)                 :: scf_energy            ! computed (or read) single determinant UHF energy
         integer, intent(in)                  :: print_level
         
!        target tensors (2e integrals are input, UHF matrix sections to be filled here).
         type(talsh_intg_tens), intent(inout) :: int_t

         integer                         :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         type(talsh_tens_t)              :: id_mat_occ
         type(talsh_tens_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)
         integer(INTD)                   :: ierr
         logical                         :: one_el_exist
         real(8)                         :: e_core
         integer                         :: min_occ
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer                         :: i


         nocc = size(eps_occ)
         nvir = size(eps_vir)

         min_occ=minval(mo_occ)

!        Get foo tensor
         fo_dims(1:2) = nocc
         ierr=talsh_tensor_construct(int_t%foo,C8,fo_dims,init_val=ZERO)

!        Get fvv tensor
         fv_dims(1:2) = nvir
         ierr=talsh_tensor_construct(int_t%fvv,C8,fv_dims,init_val=ZERO)

!        Get fov tensor
         fov_dims(1) = nocc
         fov_dims(2) = nvir
         ierr=talsh_tensor_construct(int_t%fov,C8,fov_dims,init_val=ZERO)

         select case (mode)

         case (3)

!            construct Fock matrix from orbital energies and valence 2e integrals
!            to be written by Chima Chibueze

         case (2)

!           get one electron integrals
            one_el_exist = exist_one_el(e_core)

            if (.not. one_el_exist) stop "Requested Fock matrix construction, but MRCONEE is missing"

!           read the one electron integrals and copy them into the tensors

            ierr=talsh_tensor_get_body_access(int_t%foo,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_occ,nocc,min_occ)

            ierr=talsh_tensor_get_body_access(int_t%fvv,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_vir,nvir,mo_vir,nvir,min_occ)

            ierr=talsh_tensor_get_body_access(int_t%fov,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fov_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_vir,nvir,min_occ)

            call print_date('one electron integrals added to Fock Matrix')

            !create density (identity) matrix for the occupied space
            ierr=talsh_tensor_construct(id_mat_occ,C8,fo_dims,init_val=ZERO)
            ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, nocc
              h_tens(i,i) = ONE
            enddo

            ierr=talsh_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",int_t%foo,int_t%oooo,id_mat_occ, ONE)
            ierr=talsh_tensor_contract("F(i,a)+=V(k,i,l,a)*I(k,l)",int_t%fov,int_t%ooov,id_mat_occ, ONE)
            ierr=talsh_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",int_t%fvv,int_t%vovo,id_mat_occ, ONE)

            call print_date('two electron integrals added to Fock Matrix')

            if (print_level >10) then
              call print_tensor(int_t%foo, 1.0D-6, 'fo_tensor')
              call print_tensor(int_t%fov, 1.0D-6, 'fov_tensor')
              call print_tensor(int_t%fvv, 1.0D-6, 'fv_tensor')
            end if

            ierr=talsh_tensor_destruct(id_mat_occ)

         case (1)

            ! Fock matrix is defined as diagonal matrix with orbital energies on the diagonal
            ierr=talsh_tensor_get_body_access(int_t%foo,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            do i=1,fo_dims(1)
               h_tens(i,i)=eps_occ(i)
            end do

            ierr=talsh_tensor_get_body_access(int_t%fvv,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            do i=1,fv_dims(1)
               h_tens(i,i)=eps_vir(i)
            end do

            call get_scf ( scf_energy )
            write(*,*) ""
            write (*,'(A,F23.15)') "    SCF energy ( read ) = ", scf_energy
            write(*,*) ""

            if (print_level.gt.3) call print_orbital_energy(eps_occ,nocc,eps_vir,nvir,1)

         end select

         if (print_level >0) call print_date('Fock matrix initialized')

        end subroutine init_fock_talsh

        subroutine define_H0 (scf_energy,eps_occ,eps_vir,comm_t,oooo,level_shift,print_level,mode,nff)

!        Written by Lucas Visscher, spring 2021
!        Shifts diagonal part of the T1-equations to the LHS, defines these epsilons ("orbital energies") as H0 
!        and make foo and fvv consistent with this choice
!        recomputes scf_energy if mode > 1

         use exacorr_global

         implicit none

         integer, intent(in)                  :: mode ! for 1 the scf energy is not recomputed
         real(8), intent(inout)               :: eps_occ(:),eps_vir(:) ! (level-shifted) orbital energies, aka H0
         real(8), intent(in)                  :: level_shift
         integer, intent(in)                  :: print_level
         integer, intent(in)                  :: nff
         real(8), intent(inout)               :: scf_energy            ! recomputed single determinant UHF energy


!        common tensors: the foo and fvv parts get modfied here (diagonal is shifted to the left)
         type(talsh_comm_tens), intent(inout) :: comm_t
         type(talsh_tens_t)                   :: oooo

         integer                         :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         type(talsh_tens_t)              :: diag_mat_occ, diag_mat_vir
         type(talsh_tens_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)
         integer(INTD)                   :: ierr
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer                         :: i

         nocc = size(eps_occ)
         nvir = size(eps_vir)
         fo_dims(1:2) = nocc
         fv_dims(1:2) = nvir

         ! The epsilons are initialized from the diagonal of the UHF Fock matrix
         ierr=talsh_tensor_get_body_access(comm_t%foo,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
         do i=1,fo_dims(1)
            eps_occ(i)=h_tens(i,i)
         end do

         ierr=talsh_tensor_get_body_access(comm_t%fvv,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
         do i=1,fv_dims(1)
            eps_vir(i)=h_tens(i,i)
         end do

         !recompute SCF energy
         if (mode>1) then
            call compute_scf_energy(scf_energy,comm_t%foo,oooo,fo_dims,print_level)
         else
            if (nff.gt.0) stop "finite field requires core energy (MRCONEE)"
         end if

         !apply level shift
         call shift_orbital_energy(eps_vir,eps_occ,level_shift)
         if (print_level.gt.3) call print_orbital_energy(eps_occ,nocc,eps_vir,nvir,1)

         !put these orbital energies back at the diagonal to form H0
         ierr=talsh_tensor_construct(diag_mat_occ,C8,fo_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(diag_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
            h_tens(i,i) = eps_occ(i)
         end do

         ierr=talsh_tensor_construct(diag_mat_vir,C8,fv_dims,init_val=ZERO)
         ierr=talsh_tensor_get_body_access(diag_mat_vir,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
           h_tens(i,i) = eps_vir(i)
         end do

         ! create foo and fvv with H0 subtracted
         ierr=talsh_tensor_add("F(i,j)+=D(i,j)",comm_t%foo,diag_mat_occ,scale=MINUS_ONE)
         ierr=talsh_tensor_add("F(a,b)+=D(a,b)",comm_t%fvv,diag_mat_vir,scale=MINUS_ONE)
         if (print_level >0) then 
           write(*,*)
           write(*,*) "    Norm of Fock matrix (OO, no diagonal) = ", tensor_norm2(comm_t%foo)
           write(*,*) "    Norm of Fock matrix (OV, all        ) = ", tensor_norm2(comm_t%fov)
           write(*,*) "    Norm of Fock matrix (VV, no diagonal) = ", tensor_norm2(comm_t%fvv)
           write(*,*) 
         end if

         ierr=talsh_tensor_destruct(diag_mat_occ)
         ierr=talsh_tensor_destruct(diag_mat_vir)

        end subroutine define_H0

        subroutine compute_scf_energy(scf_energy,foo,oooo,foo_dims,print_level)

            use exacorr_global

!           Refactored by by Johann Pototschnig, spring 2021

            type(talsh_tens_t)                 :: foo
            type(talsh_tens_t)                 :: oooo
            real(8), intent(inout)             :: scf_energy
            integer(INTD)                      :: foo_dims(1:2)
            integer                            :: print_level

            type(talsh_tens_t)              :: id_mat_occ, h_tensor
            logical                         :: one_el_exist
            real(8)                         :: e_core
            integer(INTD)                   :: ierr
            type(C_PTR)                     :: body_p
            complex(8), pointer, contiguous :: h_tens(:,:)
            type(talsh_tens_t)              :: result_tensor
            integer(INTD)                   :: result_dims(1)
            complex(8), pointer             :: result_tens(:)
            integer                         :: i

!           get one electron integrals
            one_el_exist = exist_one_el(e_core)

            if (.not. one_el_exist) stop "Requested Fock matrix computation, but MRCONEE is missing"

            !Make a result tensor and arrange access to it via a Fortran pointer
            result_dims(1) = 1
            ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
            ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,result_tens,result_dims)

            !create density (identity) matrix for the occupied space
            ierr=talsh_tensor_construct(id_mat_occ,C8,foo_dims,init_val=ZERO)
            ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,foo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, foo_dims(1)
              h_tens(i,i) = ONE
            enddo
            ! recompute HF energy
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("R()=F(i,j)*I(i,j)",result_tensor,foo,id_mat_occ)
            ierr=talsh_tensor_construct(h_tensor,C8,foo_dims,init_val=ZERO)
            ierr=talsh_tensor_contract("H(i,j)+=G(i,k,j,l)*I(k,l)",h_tensor,oooo,id_mat_occ, MINUS_ONE_HALF)
            ierr=talsh_tensor_contract("R()+=H(i,j)*I(i,j)",result_tensor,h_tensor,id_mat_occ)
            ierr=talsh_tensor_destruct(h_tensor)

            scf_energy=real(result_tens(1),8) ! Casting into real(8)
            scf_energy=scf_energy+e_core
            if (print_level>-1) then
               write(*,*) "" 
               write (*,'(A,F23.15)') "   Recomputed energy (single determinant, calculated) = ", scf_energy
               write(*,*) ""
            end if

            ierr=talsh_tensor_destruct(result_tensor)
            ierr=talsh_tensor_destruct(id_mat_occ)

        end subroutine compute_scf_energy

        subroutine set_fock(exa_input,mo_occ,mo_vir,int_t,comm_t,iff)

!        Refactored by Lucas Visscher from routine by Johann Pototschnig, spring 2021
!        Fock matrix is copied from int_t to comm_t and finite field is added if available

         use exacorr_global

         implicit none

         type(exacc_input), intent(in)        :: exa_input
         integer, intent(in)                  :: mo_occ(:),mo_vir(:) ! the orbital indices
         integer, intent(in)                  :: iff !integer to loop over field strengths

!        integral tensors (for technical reasons defined as inout, but not to be changed here).
         type(talsh_intg_tens), intent(inout) :: int_t
!        target tensors to be augmented with the finite field
         type(talsh_comm_tens), intent(inout) :: comm_t

         integer                         :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         type(talsh_tens_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)
         integer(INTD)                   :: ierr
         integer                         :: i,min_occ
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)

!        copy the matrices
         ierr=talsh_tensor_init(comm_t%foo)
         ierr=talsh_tensor_init(comm_t%fvv)
         ierr=talsh_tensor_init(comm_t%fov)
         ierr=talsh_tensor_add("F(i,j)+=D(i,j)",comm_t%foo,int_t%foo)
         ierr=talsh_tensor_add("F(a,b)+=D(a,b)",comm_t%fvv,int_t%fvv)
         ierr=talsh_tensor_add("F(i,a)+=D(i,a)",comm_t%fov,int_t%fov)

!        if there is a finte field add it 
         if (exa_input%nff(1)>0 .and. iff>-1) then
            nocc = exa_input%nocc
            nvir = exa_input%nvir

            min_occ=minval(mo_occ)

!           Get foo tensor
            fo_dims(1:2) = nocc

!           Get fvv tensor
            fv_dims(1:2) = nvir

!           Get fov tensor
            fov_dims(1) = nocc
            fov_dims(2) = nvir

            ierr=talsh_tensor_get_body_access(comm_t%foo,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            call add_finitefield(h_tens,exa_input,iff,mo_occ,nocc,mo_occ,nocc,min_occ)

            ierr=talsh_tensor_get_body_access(comm_t%fvv,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            call add_finitefield(h_tens,exa_input,iff,mo_vir,nvir,mo_vir,nvir,min_occ)

            ierr=talsh_tensor_get_body_access(comm_t%fov,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fov_dims) ! to use as a regular Fortran array
            call add_finitefield(h_tens,exa_input,iff,mo_occ,nocc,mo_vir,nvir,min_occ)

            do i=1,exa_input%nff(1)
               write(*,*) ""
               write (*,'(3A,2E12.3)') "    Added finite field ",exa_input%ff_names(i), &
                                    " with strength ",exa_input%ff(i,iff)
               write(*,*) ""
            end do

            call print_date('Finite field added to Fock Matrix')
         end if

        end subroutine set_fock

        subroutine get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t,moint_scheme,print_level)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list.

!        Written by Lucas Visscher, summer 2017

         use exacorr_global
         use talsh_ao_to_mo

         integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(talsh_intg_tens),intent(inout)  :: int_t
         integer, intent(in)                  :: print_level
         integer, intent(in)                  :: moint_scheme

         integer(INTD) :: oovv_dims(1:4),vvvv_dims(1:4),oooo_dims(1:4),vovo_dims(1:4),ooov_dims(1:4),vovv_dims(1:4)

         integer(INTD)                      :: ierr
         complex(8), pointer, contiguous    :: ao_tens(:,:,:,:)
         type(C_PTR)                        :: body_p

         integer :: nao

!        auxilliary integral tensors (only needed inside this routine)
         type(talsh_tens_t) :: aoint_tensor, ovvo_tensor
         integer(INTD)      :: aoint_dims(1:4), ovvo_dims(1:4)

!        arrays needed to communicate with the MO integral generator
         integer:: nmo(4)
         integer, allocatable :: mo_list(:)

!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: minusone_tensor
         integer(INTD)      :: minusone_dims(1)

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get tensor with AO integrals
         aoint_dims      = nao
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
         if (ierr /= 0) then
            print*, "Could not construct ao integral tensor, increasing TALSH_BUFF may help"
            call quit ('error constructing aoint_tensor')
         end if

         call print_date('allocated AO integral tensor')
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens_complex=ao_tens)
         call print_date('finished ao integral calculation')

!        Get antisymmetric oovv tensor
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
!        whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
         if (moint_scheme<5) then
            allocate (mo_list(2*nocc+2*nvir))
            nmo(1) = nocc
            nmo(2) = nvir
            nmo(3) = nocc
            nmo(4) = nvir
            mo_list(1:nocc)                      = mo_occ(1:nocc)
            mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
            mo_list(nocc+nvir+1:2*nocc+nvir)     = mo_occ(1:nocc)
            mo_list(2*nocc+nvir+1:2*nocc+2*nvir) = mo_vir(1:nvir)
            call get_integral_tensor (aoint_tensor,int_t%oovv,12,nmo,mo_list)
            deallocate (mo_list)
         end if
         call print_date('finished oovv integrals')

!        Get antisymmetric oooo tensor
         oooo_dims = nocc
         ierr=talsh_tensor_construct(int_t%oooo,C8,oooo_dims,init_val=ZERO)
         if (moint_scheme<5) then
            allocate (mo_list(4*nocc))
            nmo = nocc
            mo_list(1:nocc)          = mo_occ(1:nocc)
            mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
            mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
            mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
            call get_integral_tensor (aoint_tensor,int_t%oooo,12,nmo,mo_list)
            deallocate (mo_list)
         end if
         call print_date('finished oooo integrals')

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=talsh_tensor_construct(int_t%vovo,C8,vovo_dims,init_val=ZERO)
         if (moint_scheme<5) then
            allocate (mo_list(2*nvir+2*nocc))
            nmo(1:2) = nvir
            nmo(3:4) = nocc
            mo_list(1:nvir)                      = mo_vir(1:nvir)
            mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
            mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
            mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
            call get_integral_tensor (aoint_tensor,int_t%vovo,0,nmo,mo_list)
            deallocate (mo_list)
!           Get ovvo tensor
            ovvo_dims(1) = nocc
            ovvo_dims(2) = nvir
            ovvo_dims(3) = nvir
            ovvo_dims(4) = nocc
            ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
            allocate (mo_list(2*nvir+2*nocc))
            nmo(1) = nocc
            nmo(2) = nvir
            nmo(3) = nvir
            nmo(4) = nocc
            mo_list(1:nocc)                      = mo_occ(1:nocc)
            mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
            mo_list(nvir+nocc+1:2*nvir+nocc)     = mo_vir(1:nvir)
            mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
            call get_integral_tensor (aoint_tensor,ovvo_tensor,0,nmo,mo_list)
            deallocate (mo_list)
!           Make anti-symmetrized vovo tensor
            minusone_dims(1) = 1
            ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
            ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,minusone_tensor)
            ierr=talsh_tensor_destruct(ovvo_tensor)
         end if
         call print_date('finished vovo integrals')

!        Get vvvv tensor
         vvvv_dims = nvir
         ierr=talsh_tensor_construct(int_t%vvvv,C8,vvvv_dims,init_val=ZERO)
         if (moint_scheme<5) then
            allocate (mo_list(4*nvir))
            nmo = nvir
            mo_list(1:nvir)          = mo_vir(1:nvir)
            mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
            mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
            mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
            call get_integral_tensor (aoint_tensor,int_t%vvvv,12,nmo,mo_list)
            deallocate (mo_list)
         end if
         call print_date('finished vvvv integrals')

!        Get ooov tensor
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=talsh_tensor_construct(int_t%ooov,C8,ooov_dims,init_val=ZERO)
         if (moint_scheme<5) then
            allocate (mo_list(3*nocc+nvir))
            nmo(1:3) = nocc
            nmo(4)   = nvir
            mo_list(1:nocc)               = mo_occ(1:nocc)
            mo_list(nocc+1:2*nocc)        = mo_occ(1:nocc)
            mo_list(2*nocc+1:3*nocc)      = mo_occ(1:nocc)
            mo_list(3*nocc+1:3*nocc+nvir) = mo_vir(1:nvir)
            call get_integral_tensor (aoint_tensor,int_t%ooov,12,nmo,mo_list)
            deallocate (mo_list)
         end if
         call print_date('finished ooov integrals')

!        Get vovv tensor
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(int_t%vovv,C8,vovv_dims,init_val=ZERO)
         if (moint_scheme<5) then
            allocate (mo_list(3*nvir+nocc))
            nmo(1)   = nvir
            nmo(2)   = nvir
            nmo(3)   = nocc
            nmo(4)   = nvir
            mo_list(1:nvir)                    = mo_vir(1:nvir)
            mo_list(nvir+1:2*nvir)             = mo_vir(1:nvir)
            mo_list(2*nvir+1:2*nvir+nocc)      = mo_occ(1:nocc)
            mo_list(2*nvir+nocc+1:3*nvir+nocc) = mo_vir(1:nvir)
            call get_integral_tensor (aoint_tensor,int_t%vovv,34,nmo,mo_list)
            deallocate (mo_list)
         end if
         call print_date('finished vovv integrals')

!        Clean-up
         ierr=talsh_tensor_destruct(aoint_tensor)

       end subroutine get_CC_integrals

       subroutine diis (t2_tensor,t1_tensor,t_convergence,mode)

         integer, parameter :: maxdim = 5 ! maximum number of tensors that will be stored
         integer, parameter :: mindim = 3 ! minimum number of tensors to allow extrapolation
         integer, intent(in), optional  :: mode ! 1:initialize, 2:clean-up
         real(8), intent(out), optional :: t_convergence ! convergence of the t1 and t2 amplitudes

         type(talsh_tens_t), intent(inout)  :: t2_tensor
         type(talsh_tens_t), intent(inout), optional  :: t1_tensor
         integer             :: i
         integer(INTD)       :: ierr, tens_rank
         integer(INTD), save :: t1_dims(1:2), t2_dims(1:4)
         logical             :: init, clean_up, CCD
         integer, save       :: itens, ntens
         complex(8), save    :: b_diis(maxdim,maxdim)
         complex(8)          :: bb(0:maxdim,0:maxdim),coeff(0:maxdim)
         real(8)             :: b_max, b_scale
         type(talsh_tens_t), save :: t1_inp_tensor, t2_inp_tensor
         type(talsh_tens_t), save :: t1_res_tensor(maxdim), t2_res_tensor(maxdim)
         type(talsh_tens_t), save :: t1_err_tensor(maxdim), t2_err_tensor(maxdim)
         integer             :: ipiv(0:maxdim), info  ! needed for ZHESV solver
         real(8)             :: work(maxdim,maxdim)   ! needed for ZHESV solver
         type(talsh_tens_t)  :: scalar_tensor
         integer(INTD)       :: scalar_dims(1)
         complex(8), pointer :: scalar_tens(:)
         type(C_PTR):: body_p

         init       = .false.
         clean_up   = .false.
         CCD = .false.

         if (present(mode)) then
            if (mode.eq.1) init       = .true.
            if (mode.eq.2) clean_up   = .true.
         end if

         if (.not.present(t1_tensor)) CCD = .true.

         if (init) then
!           get dimensions of the tensors
            if (.not.CCD) then
                ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
                if (ierr.ne.0 .or. tens_rank.ne.2) stop 'input t1-tensor corrupted in diis'
            end if
            ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
            if (ierr.ne.0 .or. tens_rank.ne.4) stop 'input t2-tensor corrupted in diis'
!           initialize DIIS matrix
            b_diis      = ZERO
!           create storage for input, result and error tensors (in memory for the time being)
            if (.not.CCD) ierr=talsh_tensor_construct(t1_inp_tensor,C8,t1_dims,init_val=ZERO)
            ierr=talsh_tensor_construct(t2_inp_tensor,C8,t2_dims,init_val=ZERO)
            do i = 1, maxdim
               if (.not.CCD) then
                    ierr=talsh_tensor_construct(t1_res_tensor(i),C8,t1_dims,init_val=ZERO)
                    ierr=talsh_tensor_construct(t1_err_tensor(i),C8,t1_dims,init_val=ZERO)
               end if
               ierr=talsh_tensor_construct(t2_res_tensor(i),C8,t2_dims,init_val=ZERO)
               ierr=talsh_tensor_construct(t2_err_tensor(i),C8,t2_dims,init_val=ZERO)
               if (ierr.ne.0) stop 'not enough space for diis'
            end do
            itens = 0
            ntens = 0
!           store the tensors in memory as current input tensors, needed to calculate the error vector later on
            if (.not.CCD) then
                ierr=talsh_tensor_init(t1_inp_tensor)
                ierr=talsh_tensor_add("T1O(a,i)+=T1N(a,i)",t1_inp_tensor,t1_tensor)
            end if
            ierr=talsh_tensor_init(t2_inp_tensor)
            ierr=talsh_tensor_add("T2O(a,b,i,j)+=T2N(a,b,i,j)",t2_inp_tensor,t2_tensor)
            return ! no extrapolation possible
         end if

         if (clean_up) then
!           de-allocate tensors
            if (.not.CCD) ierr=talsh_tensor_destruct(t1_inp_tensor)
            ierr=talsh_tensor_destruct(t2_inp_tensor)
            do i = 1, maxdim
               if (.not.CCD) then
                ierr=talsh_tensor_destruct(t1_res_tensor(i))
                ierr=talsh_tensor_destruct(t1_err_tensor(i))
               end if
               ierr=talsh_tensor_destruct(t2_res_tensor(i))
               ierr=talsh_tensor_destruct(t2_err_tensor(i))
            end do
            return
         end if

!        define scalar tensor needed to extract the error
         scalar_dims(1) = 1
         ierr=talsh_tensor_construct(scalar_tensor,C8,scalar_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(scalar_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,scalar_tens,scalar_dims)

!        update counters
         itens = itens + 1
         if (itens.gt.maxdim) itens = 1       ! cycle because we already have the maximum number of stored amplitudes, will overwrite the oldest one
         ntens = ntens + 1                    ! number of tensors stored after returning
         if (ntens.gt.maxdim) ntens = maxdim  ! we have reached the maximum number of tensors to be used in DIIS

!        store the result tensors
         if (.not.CCD) then
            ierr=talsh_tensor_init(t1_res_tensor(itens))
            ierr=talsh_tensor_add("T1R(a,i)+=T1(a,i)",t1_res_tensor(itens),t1_tensor)
         end if
         ierr=talsh_tensor_init(t2_res_tensor(itens))
         ierr=talsh_tensor_add("T2R(a,b,i,j)+=T2(a,b,i,j)",t2_res_tensor(itens),t2_tensor)

!        calculate the error tensors by comparing the result of the t-equations with the tensors used as input
         if (.not. CCD) then
            ierr=talsh_tensor_init(t1_err_tensor(itens))
            ierr=talsh_tensor_add("T1E(a,i)+=T1R(a,i)",t1_err_tensor(itens),t1_tensor)
            ierr=talsh_tensor_add("T1E(a,i)+=T1I(a,i)",t1_err_tensor(itens),t1_inp_tensor,scale=MINUS_ONE)
         end if
         ierr=talsh_tensor_init(t2_err_tensor(itens))
         ierr=talsh_tensor_add("T2E(a,b,i,j)+=T2R(a,b,i,j)",t2_err_tensor(itens),t2_tensor)
         ierr=talsh_tensor_add("T2E(a,b,i,j)+=T2I(a,b,i,j)",t2_err_tensor(itens),t2_inp_tensor,scale=MINUS_ONE)

!        update the b_diis array that holds the scalar products of the error tensors
         ierr=talsh_tensor_init(scalar_tensor)
         do i = 1, ntens
            ierr=talsh_tensor_init(scalar_tensor)
            if (.not.CCD) ierr=talsh_tensor_contract("E()+=T+(a,i)*T(a,i)",scalar_tensor,t1_err_tensor(itens),t1_err_tensor(i))
            ierr=talsh_tensor_contract("E()+=T+(a,b,i,j)*T(a,b,i,j)",scalar_tensor,t2_err_tensor(itens),t2_err_tensor(i))
            b_diis(itens,i) = scalar_tens(1)
            b_diis(i,itens) = dconjg(scalar_tens(1))
         end do

         if (present(t_convergence)) then
            t_convergence = dsqrt(dble(b_diis(itens,itens)))
         end if

!        transfer to a work array in the way requested by the solver and normalize the matrix
         b_max = maxval (abs(b_diis))
         if (b_max .gt. 1.D-30) then
             b_scale = 1.D0 / b_max
         else
             write(*,*) "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
             return
         end if
         bb = ZERO
         bb(1:ntens,0) = MINUS_ONE
         bb(0,1:ntens) = MINUS_ONE
         bb(1:ntens,1:ntens) = b_diis(1:ntens,1:ntens) * b_scale

!        solve the linear system of equations using the LAPACK solver ZHESV
         coeff(0)       = MINUS_ONE  ! on input this should be the right-hand-side of the matrix equation b coef = rhs
         coeff(1:ntens) = ZERO
         call ZHESV ('L',ntens+1,1,bb,maxdim+1,ipiv,coeff,maxdim+1,work,maxdim*maxdim,info)

!        make the extrapolated amplitudes if the space is large enough
         if (ntens.ge.mindim) then
           if (.not.CCD) ierr=talsh_tensor_init(t1_tensor)
           ierr=talsh_tensor_init(t2_tensor)
            do i = 1, ntens
               if (.not.CCD) ierr=talsh_tensor_add("T1(a,i)+=T1O(a,i)",t1_tensor,t1_res_tensor(i),scale=coeff(i))
               ierr=talsh_tensor_add("T2(a,b,i,j)+=T2O(a,b,i,j)",t2_tensor,t2_res_tensor(i),scale=coeff(i))
            end do
         end if

!        the extrapolated amplitudes are the new input tensors
         if (.not.CCD) then
            ierr=talsh_tensor_init(t1_inp_tensor)
            ierr=talsh_tensor_add("T1O(a,i)+=T1(a,i)",t1_inp_tensor,t1_tensor)
         end if
         ierr=talsh_tensor_init(t2_inp_tensor)
         ierr=talsh_tensor_add("T2O(a,b,i,j)+=T2(a,b,i,j)",t2_inp_tensor,t2_tensor)

         ierr=talsh_tensor_destruct(scalar_tensor)

       end subroutine diis

       subroutine t3_energy (t_energy,eps_occ,eps_vir,t1_tensor,t2_tensor,fov_tensor, &
                             ooov_tensor,vovv_tensor,oovv_tensor,print_level)

!        Evaluate +T, -T and (T) energy terms.
!        Definitions of terms: L. Visscher, T.J. Lee, K.G. Dyall, J Chem Phys. 105 (1996) 8769.
         
         type(talsh_tens_t), intent(inout)  :: t1_tensor, t2_tensor
         type(talsh_tens_t), intent(inout)  :: fov_tensor, ooov_tensor, vovv_tensor, oovv_tensor
         real(8), intent(in)                :: eps_occ(:),eps_vir(:)
         real(8), intent(out)               :: t_energy(3)
         integer, intent(in)                :: print_level

         integer(INTD)     :: t2_dims(4), abc_dims(3), ab_dims(2), ia_dims(2), abi_dims(3)
         type(talsh_tens_t):: w_tensor, y_tensor, we_tensor
         type(talsh_tens_t):: vovv_i_tensor, vovv_j_tensor, vovv_k_tensor
         type(talsh_tens_t):: ooov_ij_tensor, ooov_jk_tensor, ooov_ki_tensor
         type(talsh_tens_t):: t2_ij_tensor, t2_jk_tensor, t2_ki_tensor
         type(talsh_tens_t):: t2_i_tensor, t2_j_tensor, t2_k_tensor

!        scalars (need to be defined as tensor types)
         type(C_PTR)         :: body_p
         type(talsh_tens_t)  :: result_tensor
         integer(INTD)       :: result_dims(1)
         complex(8), pointer :: result_tens(:)

         real(8)        :: eps_ijk
         integer(INTD)  :: ierr, tens_rank
         integer        :: nocc, nvir, i, j, k
         complex(8), allocatable :: y(:,:,:) ! used in form_y_intermediate, but allocated here to have fewer alloc/dealloc calls

!        Get nocc and nvir from tensor information
         ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
         if (ierr.ne.0 .or. tens_rank.ne.4) stop 'error: t2 tensor corrupted'

         if (t2_dims(1).eq.t2_dims(2)) then
            nvir = t2_dims(1)
         else
            stop 'error: asymmetric tensor in t3_energy'
         end if

         if (t2_dims(3).eq.t2_dims(4)) then
            nocc = t2_dims(3)
         else
            stop 'error: asymmetric tensor in t3_energy'
         end if

!        Initialize scalars that are to be used as tensors in contractions
         result_dims(1) = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

!        Allocate tensors of dimension nvir**3 to hold intermediates
         abc_dims = nvir
         ierr=talsh_tensor_construct(w_tensor,C8,abc_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(y_tensor,C8,abc_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(we_tensor,C8,abc_dims,init_val=ZERO)
         allocate(y(nvir,nvir,nvir))

!        The following is a work-around, we need to implement contraction with restricted range for "inactive" indices
         ab_dims  = nvir
         ia_dims(1) = nocc
         ia_dims(2) = nvir
         abi_dims(1:2) = nvir
         abi_dims(3)   = nocc
         ierr=talsh_tensor_construct(vovv_i_tensor,C8,abc_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(vovv_j_tensor,C8,abc_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(vovv_k_tensor,C8,abc_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(t2_jk_tensor,C8,ab_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(t2_ki_tensor,C8,ab_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(t2_ij_tensor,C8,ab_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(ooov_ij_tensor,C8,ia_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(ooov_jk_tensor,C8,ia_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(ooov_ki_tensor,C8,ia_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(t2_i_tensor,C8,abi_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(t2_j_tensor,C8,abi_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(t2_k_tensor,C8,abi_dims,init_val=ZERO)

!        Get access to result
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         t_energy = ZERO
         do k = 3, nocc
            call get_subtensor (vovv_tensor,vovv_k_tensor,432,k)
            call get_subtensor (t2_tensor,t2_k_tensor,433,k)
            do j = 2, k-1
               call get_subtensor (vovv_tensor,vovv_j_tensor,432,j)
               call get_subtensor (t2_tensor,t2_j_tensor,433,j)
               call get_subtensor (t2_tensor,t2_jk_tensor,4334,j,k)
               call get_subtensor (ooov_tensor,ooov_jk_tensor,4312,j,k)
               do i = 1, j-1
                  call get_subtensor (vovv_tensor,vovv_i_tensor,432,i)
                  call get_subtensor (t2_tensor,t2_i_tensor,433,i)
                  call get_subtensor (t2_tensor,t2_ij_tensor,4334,i,j)
                  call get_subtensor (t2_tensor,t2_ki_tensor,4334,k,i)
                  call get_subtensor (ooov_tensor,ooov_ij_tensor,4312,i,j)
                  call get_subtensor (ooov_tensor,ooov_ki_tensor,4312,k,i)

                  ! Form w-intermediate, first without symmetrization and minus sign in the auxilliary tensor array we
                  ierr=talsh_tensor_init(we_tensor)
                  ! For this term we use the permutation symmetry <vv||ov> = - <vo||vv>*, with vvov(abe:i) = vovv*(eab:i)
                  ierr=talsh_tensor_contract("W(a,b,c)+=V+(e,a,b)*T(e,c)",we_tensor,vovv_i_tensor,t2_jk_tensor) ! o(3)v(4)
                  ierr=talsh_tensor_contract("W(a,b,c)+=V+(e,a,b)*T(e,c)",we_tensor,vovv_j_tensor,t2_ki_tensor) ! o(3)v(4)
                  ierr=talsh_tensor_contract("W(a,b,c)+=V+(e,a,b)*T(e,c)",we_tensor,vovv_k_tensor,t2_ij_tensor) ! o(3)v(4)
                  ! For this term we use the permutation symmetry <vo||oo> = - <oo||ov>*, with vooo(al:ij) = ooov*(la:ij)
                  ierr=talsh_tensor_contract("W(a,b,c)+=V+(l,a)*T(b,c,l)",we_tensor,ooov_ij_tensor,t2_k_tensor) ! o(4)v(3)
                  ierr=talsh_tensor_contract("W(a,b,c)+=V+(l,a)*T(b,c,l)",we_tensor,ooov_jk_tensor,t2_i_tensor) ! o(4)v(3)
                  ierr=talsh_tensor_contract("W(a,b,c)+=V+(l,a)*T(b,c,l)",we_tensor,ooov_ki_tensor,t2_j_tensor) ! o(4)v(3)
                  if(print_level > 1) write(*,*) "W[sym] =",tensor_norm2(we_tensor)

                  ! Copy into w array to symmetrize and add the minus sign omitted in the above contractions
                  ierr=talsh_tensor_init(w_tensor)
                  ! symmetrization via contraction as addition with permutation is not implemented
                  ierr=talsh_tensor_contract("W(a,b,c)+=W(a,b,c)",w_tensor,we_tensor,one_tensor,scale=MINUS_ONE)
                  ierr=talsh_tensor_contract("W(a,b,c)+=W(b,c,a)",w_tensor,we_tensor,one_tensor,scale=MINUS_ONE)
                  ierr=talsh_tensor_contract("W(a,b,c)+=W(c,a,b)",w_tensor,we_tensor,one_tensor,scale=MINUS_ONE)
                  if(print_level > 1) write(*,*) "W[+T ] =",tensor_norm2(w_tensor)

                  ! Form Y-intermediate and calculate 5th order -T energy contribution
                  call form_y3_intermediate (i,j,k,t1_tensor,t2_tensor,y_tensor,y)
                  if(print_level > 1) write(*,*) "Y[-T ] =",tensor_norm2(y_tensor)
                  ierr=talsh_tensor_init(result_tensor)
                  ierr=talsh_tensor_contract("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(3)=t_energy(3) + result_tens(1)

                  ! Form Y-intermediate and calculate 5th order (T) energy contribution
                  call form_y2_intermediate (i,j,k,t1_tensor,t2_tensor,fov_tensor,oovv_tensor,y_tensor,y)
                  if(print_level > 1) write(*,*) "Y[(T)] =",tensor_norm2(y_tensor)
                  eps_ijk = eps_occ(i) + eps_occ(j) + eps_occ(k)
                  call scale_with_denominators (eps_occ,eps_vir,t3_tensor=y_tensor,eps_ijk=eps_ijk)
                  ierr=talsh_tensor_init(result_tensor)
                  ierr=talsh_tensor_contract("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(2)=t_energy(2) + result_tens(1)

                  ! Copy w into y array and calculate 4th order +T energy contribution
                  ierr=talsh_tensor_init(y_tensor)
                  ierr=talsh_tensor_add("W(a,b,c)+=W(a,b,c)",y_tensor,w_tensor)
                  call scale_with_denominators (eps_occ,eps_vir,t3_tensor=y_tensor,eps_ijk=eps_ijk)
                  if(print_level > 1) write(*,*) " W[/D ] =",tensor_norm2(w_tensor)
                  ierr=talsh_tensor_init(result_tensor)
                  ierr=talsh_tensor_contract("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(1)=t_energy(1) + result_tens(1)
               end do
            end do
            if (print_level >1) write (*,*)  't3_energy:',t_energy
         end do

!        We looped over all (abc) combinations, scale with the factor of 1/6 needed to go to the unique ones.
         t_energy = t_energy / SIX

!        Destruct intermediate tensors
         deallocate(y)
         ierr=talsh_tensor_destruct(t2_i_tensor)
         ierr=talsh_tensor_destruct(t2_j_tensor)
         ierr=talsh_tensor_destruct(t2_k_tensor)
         ierr=talsh_tensor_destruct(t2_ij_tensor)
         ierr=talsh_tensor_destruct(t2_jk_tensor)
         ierr=talsh_tensor_destruct(t2_ki_tensor)
         ierr=talsh_tensor_destruct(vovv_i_tensor)
         ierr=talsh_tensor_destruct(vovv_j_tensor)
         ierr=talsh_tensor_destruct(vovv_k_tensor)
         ierr=talsh_tensor_destruct(ooov_ij_tensor)
         ierr=talsh_tensor_destruct(ooov_jk_tensor)
         ierr=talsh_tensor_destruct(ooov_ki_tensor)
         ierr=talsh_tensor_destruct(we_tensor)
         ierr=talsh_tensor_destruct(w_tensor)
         ierr=talsh_tensor_destruct(y_tensor)

       end subroutine t3_energy

       subroutine form_y2_intermediate (i,j,k,t1_tensor,t2_tensor,fov_tensor,oovv_tensor,y_tensor,y)

!        Get y intermediate for -T fifth order energy correction.

         type(talsh_tens_t), intent(inout) :: t1_tensor, t2_tensor, y_tensor
         type(talsh_tens_t), intent(inout) :: fov_tensor, oovv_tensor
         integer, intent(in)               :: i, j, k
         complex(8), intent(inout)         :: y(:,:,:)

         type(C_PTR):: body_p
         integer(INTD)        :: ierr, tens_rank
         complex(8), pointer  :: t1_tens(:,:), t2_tens(:,:,:,:), y_tens(:,:,:)
         complex(8), pointer  :: fov_tens(:,:),oovv_tens(:,:,:,:)
         integer(INTD)        :: t1_dims(2), t2_dims(4), y_dims(3), fov_dims(2), oovv_dims(4)
         integer              :: nvir, a, b, c

!        Determine source dimensions and get direct access to the tensor bodies

         tens_rank = talsh_tensor_rank(t1_tensor)
         if (tens_rank .ne. 2) stop 'Error: t1 tensor corrupted in form_y2_intermediate'
         ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
         ierr=talsh_tensor_get_body_access(t1_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,t1_tens,t1_dims)

         tens_rank = talsh_tensor_rank(t2_tensor)
         if (tens_rank .ne. 4) stop 'Error: t2 tensor corrupted in form_y2_intermediate'
         ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
         ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,t2_tens,t2_dims)

         tens_rank = talsh_tensor_rank(y_tensor)
         if (tens_rank .ne. 3) stop 'Error: y tensor corrupted in form_y2_intermediate'
         ierr = talsh_tensor_dimensions(y_tensor,tens_rank,y_dims)
         ierr=talsh_tensor_get_body_access(y_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,y_tens,y_dims)

         tens_rank = talsh_tensor_rank(fov_tensor)
         if (tens_rank .ne. 2) stop 'Error: fo tensor corrupted in form_y2_intermediate'
         ierr = talsh_tensor_dimensions(fov_tensor,tens_rank,fov_dims)
         ierr=talsh_tensor_get_body_access(fov_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,fov_tens,fov_dims)

         tens_rank = talsh_tensor_rank(oovv_tensor)
         if (tens_rank .ne. 4) stop 'Error: oovv tensor corrupted in form_y2_intermediate'
         ierr = talsh_tensor_dimensions(oovv_tensor,tens_rank,oovv_dims)
         ierr=talsh_tensor_get_body_access(oovv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         if (ierr.ne.0) stop 'error in y2_intermediate'
         call c_f_pointer(body_p,oovv_tens,oovv_dims)

!        Form the intermediate
         nvir = t1_dims(1)

!        Formula with only the P_ijk permutation operators
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y(a,b,c) =            dconjg(oovv_tens(i,j,a,b))*t1_tens(c,k) + t2_tens(a,b,i,j)*dconjg(fov_tens(k,c))
                  y(a,b,c) = y(a,b,c) + dconjg(oovv_tens(j,k,a,b))*t1_tens(c,i) + t2_tens(a,b,j,k)*dconjg(fov_tens(i,c))
                  y(a,b,c) = y(a,b,c) + dconjg(oovv_tens(k,i,a,b))*t1_tens(c,j) + t2_tens(a,b,k,i)*dconjg(fov_tens(j,c))
               end do
            end do
         end do

!        Do the P_abc permutation to fully symmetrize the expression
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y_tens(a,b,c)  = y(a,b,c) + y(c,a,b) + y(b,c,a)
               end do
            end do
         end do

       end subroutine form_y2_intermediate

       subroutine form_y3_intermediate (i,j,k,t1_tensor,t2_tensor,y_tensor,y)

!        Get y intermediate for -T fifth order energy correction.

         type(talsh_tens_t), intent(inout) :: t1_tensor, t2_tensor, y_tensor
         integer, intent(in)               :: i, j, k
         complex(8), intent(inout)         :: y(:,:,:)

         type(C_PTR):: body_p
         integer(INTD)        :: ierr, tens_rank
         complex(8), pointer  :: t1_tens(:,:), t2_tens(:,:,:,:), y_tens(:,:,:)
         integer(INTD)        :: t1_dims(2), t2_dims(4), y_dims(3)
         integer              :: nvir, a, b, c

!        Determine source dimensions and get direct access to the tensor bodies

         tens_rank = talsh_tensor_rank(t1_tensor)
         if (tens_rank .ne. 2) stop 'Error: t1 tensor corrupted in form_y_intermediate'
         ierr = talsh_tensor_dimensions(t1_tensor,tens_rank,t1_dims)
         ierr=talsh_tensor_get_body_access(t1_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,t1_tens,t1_dims)

         tens_rank = talsh_tensor_rank(t2_tensor)
         if (tens_rank .ne. 4) stop 'Error: t2 tensor corrupted in form_y_intermediate'
         ierr = talsh_tensor_dimensions(t2_tensor,tens_rank,t2_dims)
         ierr=talsh_tensor_get_body_access(t2_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,t2_tens,t2_dims)

         tens_rank = talsh_tensor_rank(y_tensor)
         if (tens_rank .ne. 3) stop 'Error: y tensor corrupted in form_y_intermediate'
         ierr = talsh_tensor_dimensions(y_tensor,tens_rank,y_dims)
         ierr=talsh_tensor_get_body_access(y_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,y_tens,y_dims)

!        Form the intermediate
         nvir = t1_dims(1)

!        Formula with only the P_ijk and P_ij permutation operators
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y(a,b,c) =            t1_tens(a,i)*t1_tens(b,j)*t1_tens(c,k) - t1_tens(a,j)*t1_tens(b,i)*t1_tens(c,k)
                  y(a,b,c) = y(a,b,c) + t1_tens(a,j)*t1_tens(b,k)*t1_tens(c,i) - t1_tens(a,k)*t1_tens(b,j)*t1_tens(c,i)
                  y(a,b,c) = y(a,b,c) + t1_tens(a,k)*t1_tens(b,i)*t1_tens(c,j) - t1_tens(a,i)*t1_tens(b,k)*t1_tens(c,j)
                  y(a,b,c) = y(a,b,c) / THREE
                  y(a,b,c) = y(a,b,c) + t2_tens(a,b,i,j)*t1_tens(c,k)+t2_tens(a,b,j,k)*t1_tens(c,i)+t2_tens(a,b,k,i)*t1_tens(c,j)
               end do
            end do
         end do

!        Do the P_abc permutation to fully symmetrize the expression
         do c = 1, nvir
            do b = 1, nvir
               do a = 1, nvir
                  y_tens(a,b,c)  = y(a,b,c) + y(c,a,b) + y(b,c,a)
               end do
            end do
         end do

       end subroutine form_y3_intermediate

       subroutine get_subtensor (tensor,subtensor,mode,i,j)

!        Severely limited implementation of extraction of a subtensor, as needed in t3_energy routine.
!        No shuffling of indices possible, only 4->3 (mode=1) and 4->2 (mode=2) extraction.
!        Type of extraction defined by 3 or 4 digit integer mode:
!        1st integer: dimension of source tensor
!        2nd integer: dimension of target tensor
!        3rd integer: passive index #1
!        4th integer: passive index #2

         type(talsh_tens_t), intent(inout) :: tensor, subtensor
         integer, intent(in)               :: mode,i
         integer, intent(in), optional     :: j

         type(C_PTR):: body_p
         integer(INTD)        :: ierr, tens_rank, subtens_rank
         complex(8), pointer  :: tens(:,:,:,:), subtens3(:,:,:), subtens2(:,:)
         integer(INTD)        :: tens_dims(4), subtens_dims(3)

!        Determine source dimensions and get direct access to the tensor body to enable copying
         tens_rank = talsh_tensor_rank(tensor)
         if (tens_rank .ne. 4) stop 'Error: Only extraction from 4-dimensional tensors possible'
         ierr = talsh_tensor_dimensions(tensor,tens_rank,tens_dims)
         ierr=talsh_tensor_get_body_access(tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,tens,tens_dims)

!        Determine target  dimensions and get direct access to the tensor body to enable copying
         subtens_rank = talsh_tensor_rank(subtensor)
         if (subtens_rank .eq. 3) then
            ierr = talsh_tensor_dimensions(subtensor,subtens_rank,subtens_dims)
            ierr=talsh_tensor_get_body_access(subtensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,subtens3,subtens_dims)
         elseif (subtens_rank .eq. 2) then
            ierr = talsh_tensor_dimensions(subtensor,subtens_rank,subtens_dims)
            ierr=talsh_tensor_get_body_access(subtensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,subtens2,subtens_dims(1:2))
         else
            stop 'Error: Only insertion into 3- or 2-dimensional tensors possible'
         end if

!        Copy the data into the subtensor
         if (mode .eq. 433) then
            subtens3(:,:,:) = tens(:,:,i,:)
         elseif (mode .eq. 432) then
            subtens3(:,:,:) = tens(:,i,:,:)
         elseif (mode .eq. 4312) then
            subtens2(:,:)   = tens(i,j,:,:)
         elseif (mode .eq. 4334) then
            subtens2(:,:)   = tens(:,:,i,j)
         else
            write (*,*) mode,' type of tensor extraction is not implemented'
            stop 'Error in get_subtensor'
         end if


       end subroutine get_subtensor

       subroutine t3_energy_laplace (t_energy,comm_t,int_t,eps_occ,eps_vir,exa_input)

!       Evaluate +T, -T and (T) energy terms.
!       Definitions of terms: L. Visscher, T.J. Lee, K.G. Dyall, J Chem Phys. 105 (1996) 8769.

        use type_laplace, only : laplace_info, deallocate_laplace

        real(8), intent(out)                 :: t_energy(3)
        type(talsh_comm_tens), intent(inout) :: comm_t
        type(talsh_intg_tens), intent(inout) :: int_t
        real(8), intent(in)                  :: eps_occ(:),eps_vir(:)
        type(exacc_input), intent(in)        :: exa_input

        !  variables to get laplace points
        type(laplace_info)    :: laplace

        real(8), allocatable  :: zero_occ(:),zero_vir(:),eps_occ_neg(:)
        real(8)               :: T(40),Res(40)
        integer(INTD)         :: ierr, tens_rank
        integer(INTD)         :: t2_dims(4)
        integer               :: nocc, nvir, i

!       Testing
        type(talsh_tens_t)  :: test_sum, test_h
        real(8)             :: test_val,test_ref
        integer(INTD)       :: test_ind(1:6), test_dims(6)
        real(8)             :: test_energy(3)

!       Get nocc and nvir from tensor information
        ierr = talsh_tensor_dimensions(comm_t%t2,tens_rank,t2_dims)
        if (ierr.ne.0 .or. tens_rank.ne.4) stop 'error: t2 tensor corrupted'

        if (t2_dims(1).eq.t2_dims(2)) then
          nvir = t2_dims(1)
        else
          stop 'error: asymmetric tensor in t3_energy'
        end if
        if (t2_dims(3).eq.t2_dims(4)) then
          nocc = t2_dims(3)
        else
          stop 'error: asymmetric tensor in t3_energy'
        end if

!       get a list of positive occupied orbitals
        allocate(eps_occ_neg(nocc))
        eps_occ_neg(1:nocc)=-eps_occ(1:nocc)
!       create list of empty orbital
        allocate(zero_occ(nocc))
        allocate(zero_vir(nvir))
        zero_occ=0.0
        zero_vir=0.0
!       determine Laplace quadrature parameters with the minimax algorithm
        laplace%bounds(1) = 3.0d0*(eps_vir(1)-eps_occ(nocc))
        laplace%bounds(2) = 3.0d0*(eps_vir(nvir)-eps_occ(1))

!       Testing
        test_dims=nocc 
        test_dims(1:3)=nvir
        ierr=talsh_tensor_construct(test_sum,C8,test_dims,init_val=ZERO)


        if (exa_input%print_level>1) write(*,*) laplace%bounds
        
        call set_laplace(exa_input%n_laplace,exa_input%print_level)
        call lap_driver(laplace)

!       start computation of triples
        t_energy=0
        test_energy=0
        Res=0
        do i=1,laplace%num
          write(*,*) 'laplace:',i,laplace%val(i),laplace%wei(i)
          call laplace_get_T(T,comm_t%t2,int_t%ooov,int_t%vovv,t2_dims, &
               eps_occ_neg,eps_vir,zero_occ,zero_vir,laplace%val(i)/2,exa_input%print_level)
          Res(1:19)=Res(1:19)+laplace%wei(i)*T(1:19)
          t_energy(1)=t_energy(1)+laplace%wei(i)*sum(T)

!          call laplace_total(t_t_energy,eps_occ,eps_vir,comm_t%t1,comm_t%t2,int_t%fov, &
!                             int_t%ooov,int_t%vovv,int_t%oovv,exa_input%print_level)
!          write(*,*) 'total laplace test:', t_t_energy

          ! Testing
          ierr=talsh_tensor_construct(test_h,C8,test_dims,init_val=ONE)
          call laplace_mult (test_h,laplace%val(i),eps_vir,eps_vir,eps_vir,eps_occ_neg,eps_occ_neg,eps_occ_neg)
          ierr=talsh_tensor_add("R(a,b,c,i,j,k)+=O(a,b,c,i,j,k)",test_sum,test_h,scale=DCMPLX(laplace%wei(i)))
          ierr=talsh_tensor_destruct(test_h)
           
        end do

        !Testing
        call laplace_testing(test_sum, eps_occ, eps_vir)
        ierr=talsh_tensor_destruct(test_sum)
        
        write(*,*) 'Res: '
        write(*,'(4e16.8)') Res
        write(*,*) 'Res: '
        write(*,'(5e16.8)') sum(Res(1:5)),sum(Res(6:14)),sum(Res(15:19)), &
                            sum(Res(1:2))+sum(Res(4:5))+Res(7)+sum(Res(10:16))+sum(Res(18:19))

!       clean up
        call deallocate_laplace(laplace)
        deallocate(zero_occ)
        deallocate(zero_vir)
        deallocate(eps_occ_neg)

    end subroutine t3_energy_laplace

    subroutine laplace_get_T(T,i_t2,i_ooov,i_vovv,t2_dims, &
               eps_occ,eps_vir,zero_occ,zero_vir,s, print_level)

      real(8), intent(out)              :: T(:)
      type(talsh_tens_t), intent(inout) :: i_t2,i_ooov,i_vovv
      real(8), intent(in)               :: eps_occ(:),eps_vir(:),zero_occ(:),zero_vir(:)
      real(8), intent(in)               :: s
      integer(INTD), intent(in)         :: t2_dims(4)
      integer, intent(in)               :: print_level

      type(talsh_tens_t)  :: t2_o,t2_v,ooov_o,vovv_v, I1, I2
      integer(INTD)       :: ierr
      integer(INTD)       :: dims4(4), dims2(2)
      integer             :: nvir, nocc

      ! scalars (need to be defined as tensor types)
      type(C_PTR)         :: body_p
      type(talsh_tens_t)  :: T_tensor
      integer(INTD)       :: T_dims(1)
      complex(8), pointer :: T_tens(:)

      if (print_level>1) call print_date('Entered laplace_get_T routine')
      if (print_level>1) write(*,*) '-------------'

      !Initialize scalars that are to be used as tensors in contractions
      T_dims(1) = 1
      ierr=talsh_tensor_construct(T_tensor,C8,T_dims(1:0),init_val=ZERO)
      ierr=talsh_tensor_get_body_access(T_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
      call c_f_pointer(body_p,T_tens,T_dims)
     
      nvir = t2_dims(1)
      nocc = t2_dims(3)
      T=0.0D0
      T_tens(1)=ZERO

      ! ------------------- compute first transformed integrals and amplitudes
      ierr=talsh_tensor_construct(t2_v,C8,t2_dims,init_val=ONE/SQRT(DCMPLX(nvir)))
      !ierr=talsh_tensor_add("R(e,a,i,j)+=O(e,a,i,j)",t2_v,i_t2,scale=ONE)
      call laplace_mult(t2_v,s,zero_vir,eps_vir,eps_occ,eps_occ)

      dims4=nvir
      dims4(2)=nocc
      ierr=talsh_tensor_construct(vovv_v,C8,dims4,init_val=ONE/SQRT(DCMPLX(nvir)))
      !ierr=talsh_tensor_add("R(e,i,a,b)=O(e,i,a,b)",vovv_v,i_vovv)
      !ierr=talsh_tensor_contract("R(e,i,a,b)+=O(e,i,a,b)",vovv_v,i_vovv,one_tensor,scale=MINUS_ONE)
      call laplace_mult(vovv_v,s,zero_vir,eps_occ,eps_vir,eps_vir)
      ! ------------------- compute first transformed integrals and amplitudes

      ! 0.25 T2^{fe}_{}G2^{fe}_{}
      dims2=nvir
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,e)+=X(f,c,j,k)*Y(e,c,j,k)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T2 wrong')
      ierr=talsh_tensor_contract("G(f,e)+=X+(f,i,b,a)*Y+(e,i,b,a)",I2,vovv_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G2 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,e)*G(f,e)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 1 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(1)=real(T_tens(1))
      

      ! 0.25 T5^{fe}_{ij}G4^{}_{jfie}
      dims4=nvir
      dims4(3:4)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,e,i,j)+=X(f,c,k,i)*Y(e,c,j,k)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T5 wrong')
      ierr=talsh_tensor_contract("G(f,e,i,j)+=X+(f,j,b,a)*Y+(e,i,b,a)",I2,vovv_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G4 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,e,i,j)*G(f,e,i,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 2 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(2)=real(T_tens(1))

      ! 0.25 T6^{fe}_{ik}G4^{}_{kfie}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,e,i,k)+=X(f,c,i,j)*Y(e,c,j,k)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T6 wrong')
      ierr=talsh_tensor_contract("G(f,e,i,k)+=X+(f,k,b,a)*Y+(e,i,b,a)",I2,vovv_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G4 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,e,i,k)*G(f,e,i,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 3 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(3)=real(T_tens(1))

      ! 0.25 T4^{faec}_{}G7^{ca}_{fe}
      dims4=nvir
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,a,e,c)+=X(f,a,j,k)*Y(e,c,j,k)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T4 wrong')
      ierr=talsh_tensor_contract("G(f,a,e,c)+=X+(f,i,c,b)*Y+(e,i,b,a)",I2,vovv_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G7 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,a,e,c)*G(f,a,e,c)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 4 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(4)=real(T_tens(1))

      ! 0.25 T4^{fbec}_{}G8^{cb}_{fe}
      dims4=nvir
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,b,e,c)+=X(f,b,j,k)*Y(e,c,j,k)",I1,t2_v,t2_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T4 wrong')
      ierr=talsh_tensor_contract("G(f,b,e,c)+=X+(f,i,a,c)*Y+(e,i,b,a)",I2,vovv_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G8 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,b,e,c)*G(f,b,e,c)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 5 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(5)=real(T_tens(1))

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I9^{fb}_{je}I13^{eb}_{jf}
      dims4=nvir
      dims4(4)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(f,b,e,j)+=X(f,a,i,j)*Y+(e,i,b,a)",I1,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I9 wrong')
      ierr=talsh_tensor_contract("N(f,b,e,j)+=X(e,c,j,k)*Y+(f,k,c,b)",I2,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I13 wrong')
      ierr=talsh_tensor_contract("R()+=M(f,b,e,j)*N(f,b,e,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 6 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(6)=real(T_tens(1))

      ! 0.25 I6^{eb}_{kf}I14^{fb}_{ke}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(f,b,e,k)+=X(f,a,k,i)*Y+(e,i,b,a)",I1,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I6 wrong')
      ierr=talsh_tensor_contract("N(f,b,e,k)+=X(e,c,j,k)*Y+(f,j,c,b)",I2,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I14 wrong')
      ierr=talsh_tensor_contract("R()+=M(f,b,e,k)*N(f,b,e,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 7 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(7)=real(T_tens(1))

      ! 0.25 I6^{fa}_{je}I14^{ea}_{jf}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(f,a,e,j)+=X(f,b,i,j)*Y+(e,i,b,a)",I1,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I6 wrong')
      ierr=talsh_tensor_contract("N(f,a,e,j)+=X(e,c,j,k)*Y+(f,k,a,c)",I2,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I14 wrong')
      ierr=talsh_tensor_contract("R()+=M(f,a,e,j)*N(f,a,e,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 8 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(8)=real(T_tens(1))

      ! 0.25 I9^{ea}_{kf}I13^{fa}_{ke}
      dims4=nvir
      dims4(4)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(e,a,f,k)+=X(f,b,k,i)*Y+(e,i,b,a)",I1,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I9 wrong')
      ierr=talsh_tensor_contract("N(e,a,f,k)+=X(e,c,j,k)*Y+(f,j,a,c)",I2,t2_v,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I13 wrong')
      ierr=talsh_tensor_contract("R()+=M(e,a,f,k)*N(e,a,f,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 6 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(9)=real(T_tens(1))

      ! ------------------- compute next transformed integrals and amplitudes
      ierr=talsh_tensor_construct(t2_o,C8,t2_dims,init_val=ONE/SQRT(DCMPLX(nocc)))
      !ierr=talsh_tensor_add("R(a,b,i,m)+=O(a,b,i,m)",t2_o,i_t2)
      call laplace_mult(t2_o,s,eps_vir,eps_vir,eps_occ,zero_occ)

      dims4=nocc
      dims4(4)=nvir
      ierr=talsh_tensor_construct(ooov_o,C8,dims4,init_val=ONE/SQRT(DCMPLX(nocc)))
      !ierr=talsh_tensor_add("R(i,j,m,a)=O(i,j,m,a)",ooov_o,i_ooov)
      !ierr=talsh_tensor_contract("R(i,j,m,a)+=O(i,j,m,a)",ooov_o,i_ooov,one_tensor,scale=MINUS_ONE)
      call laplace_mult(ooov_o,s,eps_occ,eps_occ,zero_occ,eps_vir)
      ! ------------------- compute next transformed integrals and amplitudes

      ! 0.25 T11^{ae}_{nk}G9^{na}_{ke}
      dims4=nvir
      dims4(3:4)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(a,e,n,k)+=X(e,c,j,k)*Y(c,a,j,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T11 wrong')
      ierr=talsh_tensor_contract("G(a,e,n,k)+=X+(e,i,b,a)*Y+(i,k,n,b)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G9 wrong')
      ierr=talsh_tensor_contract("R()+=T(a,e,n,k)*G(a,e,n,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 21 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(21)=real(T_tens(1))

      ! 0.25 T10^{ae}_{nj}G10^{na}_{ej}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(a,e,n,j)+=X(e,c,j,k)*Y(c,a,k,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T10 wrong')
      ierr=talsh_tensor_contract("G(a,e,n,j)+=X+(e,i,b,a)*Y+(j,i,n,b)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G10 wrong')
      ierr=talsh_tensor_contract("R()+=T(a,e,n,j)*G(a,e,n,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 22 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(22)=real(T_tens(1))

      ! 0.25 T12^{fc}_{im}G11^{cm}_{fi}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,c,i,m)+=X(f,b,k,i)*Y(b,c,k,m)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T12 wrong')
      ierr=talsh_tensor_contract("G(f,c,i,m)+=X+(f,j,a,c)*Y+(j,i,m,a)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G11 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,c,i,m)*G(f,c,i,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 23 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(23)=real(T_tens(1))

      ! 0.25 T9^{fc}_{im}G12^{cm}_{fi}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,c,j,m)+=X(f,b,j,k)*Y(b,c,k,m)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T9 wrong')
      ierr=talsh_tensor_contract("G(f,c,j,m)+=X+(f,i,a,c)*Y+(j,i,m,a)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G12 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,c,j,m)*G(f,c,j,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 24 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(24)=real(T_tens(1))

      ! 0.25 T13^{be}_{nj}G13^{nb}_{je}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(b,e,n,j)+=X(e,c,j,k)*Y(b,c,k,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T13 wrong')
      ierr=talsh_tensor_contract("G(b,e,n,j)+=X+(e,i,b,a)*Y+(j,i,n,a)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G13 wrong')
      ierr=talsh_tensor_contract("R()+=T(b,e,n,j)*G(b,e,n,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 25 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(25)=real(T_tens(1))

      ! 0.25 T15^{be}_{nk}G14^{nb}_{ke}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(b,e,n,k)+=X(e,c,j,k)*Y(b,c,j,n)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T15 wrong')
      ierr=talsh_tensor_contract("G(b,e,n,k)+=X+(e,i,b,a)*Y+(k,i,n,a)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G14 wrong')
      ierr=talsh_tensor_contract("R()+=T(b,e,n,k)*G(b,e,n,k)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 26 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(26)=real(T_tens(1))

      ! 0.25 T16^{fb}_{im}G15^{bm}_{fi}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,b,i,m)+=X(f,c,k,i)*Y(b,c,k,m)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T16 wrong')
      ierr=talsh_tensor_contract("G(f,b,i,m)+=X+(f,j,b,a)*Y+(i,j,m,a)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G15 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,b,i,m)*G(f,b,i,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 27 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(27)=real(T_tens(1))

      ! 0.25 T14^{fb}_{im}G16^{bm}_{fi}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(f,b,j,m)+=X(f,c,j,k)*Y(b,c,k,m)",I1,t2_v,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T14 wrong')
      ierr=talsh_tensor_contract("G(f,b,j,m)+=X+(f,i,b,a)*Y+(j,i,m,a)",I2,vovv_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G16 wrong')
      ierr=talsh_tensor_contract("R()+=T(f,b,j,m)*G(f,b,j,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 28 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(28)=real(T_tens(1))

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I3^{}_{ne}I4^{en}_{}
      dims2=nvir
      dims2(1)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_contract("M(n,e)+=X(a,b,i,n)*Y+(e,i,b,a)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I3 wrong')
      ierr=talsh_tensor_contract("N(n,e)+=X(e,c,j,k)*Y+(k,j,n,c)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I4 wrong')
      ierr=talsh_tensor_contract("R()+=M(n,e)*N(n,e)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 30 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(30)=real(T_tens(1))

      ! 0.25 I3^{}_{mf}I4^{mf}_{}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_contract("N(m,f)+=X(b,c,k,m)*Y+(f,k,c,b)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I3 wrong')
      ierr=talsh_tensor_contract("M(m,f)+=X(f,a,i,j)*Y+(j,i,m,a)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I4 wrong')
      ierr=talsh_tensor_contract("R()+=M(m,f)*N(m,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 30 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(34)=real(T_tens(1))

      ! 0.25 I1^{}_{knie}I16^{en}_{ki}
      dims4=nocc
      dims4(4)=nvir
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(k,n,i,e)+=X(a,b,k,n)*Y+(e,i,b,a)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=talsh_tensor_contract("N(k,n,i,e)+=X(e,c,j,k)*Y+(j,i,n,c)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I16 wrong')
      ierr=talsh_tensor_contract("F()+=M(k,n,i,e)*N(k,n,i,e)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 31 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(31)=real(T_tens(1))

      ! 0.25 I1^{}_{jnie}I7^{en}_{ij}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(j,n,i,e)+=X(a,b,j,n)*Y+(e,i,b,a)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=talsh_tensor_contract("N(j,n,i,e)+=X(e,c,j,k)*Y+(i,k,n,c)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I7 wrong')
      ierr=talsh_tensor_contract("F()+=M(j,n,i,e)*N(j,n,i,e)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 32 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(32)=real(T_tens(1))

      ! 0.25 I1^{}_{kmif}I16^{fm}_{ki}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(k,m,i,f)+=X(b,c,k,m)*Y+(f,i,c,b)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=talsh_tensor_contract("N(k,m,i,f)+=X(f,a,j,k)*Y+(j,i,m,a)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I16 wrong')
      ierr=talsh_tensor_contract("F()+=M(k,m,i,f)*N(k,m,i,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 36 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(36)=real(T_tens(1))

      ! 0.25 I1^{}_{kmjf}I7^{fm}_{kj}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(k,m,j,f)+=X(b,c,k,m)*Y+(f,j,c,b)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I1 wrong')
      ierr=talsh_tensor_contract("N(k,m,j,f)+=X(f,a,k,i)*Y+(j,i,m,a)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I7 wrong')
      ierr=talsh_tensor_contract("F()+=M(k,m,j,f)*N(k,m,j,f)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 37 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(37)=real(T_tens(1))


      ! 0.25 I2^{ecan}_{}I10^{ca}_{ne}
      dims4=nvir
      dims4(4)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(e,c,a,n)+=X(b,c,i,n)*Y+(e,i,b,a)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I10 wrong')
      ierr=talsh_tensor_contract("N(e,c,a,n)+=X(e,c,j,k)*Y+(j,k,n,a)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I2 wrong')
      ierr=talsh_tensor_contract("F()+=M(e,c,a,n)*N(e,c,a,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 33 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(33)=real(T_tens(1))

      ! 0.25 I2^{ecbn}_{}I11^{cb}_{ne}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(e,c,b,n)+=X(c,a,i,n)*Y+(e,i,b,a)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I11 wrong')
      ierr=talsh_tensor_contract("N(e,c,b,n)+=X(e,c,j,k)*Y+(j,k,n,b)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I2 wrong')
      ierr=talsh_tensor_contract("F()+=M(e,c,b,n)*N(e,c,b,n)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 35 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(35)=real(T_tens(1))

      ! 0.25 I2^{fcam}_{}I10^{ca}_{mf}
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(f,c,a,m))+=X(b,c,k,m)*Y+(f,k,b,a)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I10 wrong')
      ierr=talsh_tensor_contract("N(f,c,a,m)+=X(f,c,i,j)*Y+(j,i,m,a)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I2 wrong')
      ierr=talsh_tensor_contract("F()+=M(f,c,a,m)*N(f,c,a,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 38 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(38)=real(T_tens(1))

      ! 0.25 I2^{fbam}_{}I11^{ba}_{mf}
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(f,b,a,m))+=X(b,c,k,m)*Y+(f,k,a,c)",I1,t2_o,vovv_v)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I11 wrong')
      ierr=talsh_tensor_contract("N(f,b,a,m)+=X(f,b,i,j)*Y+(j,i,m,a)",I2,t2_v,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I2 wrong')
      ierr=talsh_tensor_contract("F()+=M(f,b,a,m)*N(f,b,a,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 39 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(39)=real(T_tens(1))

      ! ------------------- delete first transformed integrals and amplitudes
      ierr=talsh_tensor_destruct(t2_v)
      ierr=talsh_tensor_destruct(vovv_v)
      ! ------------------- delete first transformed integrals and amplitudes

      ! 0.25 T1^{}_{nm}G1^{}_{nm}
      dims2=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims2,init_val=ZERO)
      ierr=talsh_tensor_contract("T(n,m)+=X(b,c,k,n)*Y(b,c,k,m)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T1 wrong')
      ierr=talsh_tensor_contract("G(n,m)+=X+(j,i,n,a)*Y+(j,i,m,a)",I2,ooov_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G1 wrong')
      ierr=talsh_tensor_contract("R()+=T(n,m)*G(n,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 11 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(11)=real(T_tens(1))


      ! 0.25 T3^{}_{jnkm}G5^{nm}_{kj}
      dims4=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(j,n,k,m)+=X(b,c,j,n)*Y(b,c,k,m)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T3 wrong')
      ierr=talsh_tensor_contract("G(j,n,k,m)+=X+(i,k,n,a)*Y+(j,i,m,a)",I2,ooov_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G5 wrong')
      ierr=talsh_tensor_contract("R()+=T(j,n,k,m)*G(j,n,k,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 12 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(12)=real(T_tens(1))

      ! 0.25 T3^{}_{inkm}G6^{nm}_{ki}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(i,n,k,m)+=X(b,c,i,n)*Y(b,c,k,m)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T3 wrong')
      ierr=talsh_tensor_contract("G(i,n,k,m)+=X+(k,j,n,a)*Y+(j,i,m,a)",I2,ooov_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G6 wrong')
      ierr=talsh_tensor_contract("R()+=T(i,n,k,m)*G(i,n,k,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 13 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(13)=real(T_tens(1))

      ! 0.25 T7^{ab}_{nm}G3^{bnam}_{}
      dims4=nvir
      dims4(3:4)=nocc
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(a,b,n,m)+=Z(c,a,k,n)*Y(b,c,k,m)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T7 wrong')
      ierr=talsh_tensor_contract("G(a,b,n,m)+=X+(j,i,n,b)*Y+(j,i,m,a)",I2,ooov_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G3 wrong')
      ierr=talsh_tensor_contract("R()+=T(a,b,n,m)*G(a,b,n,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 14 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(14)=real(T_tens(1))

      ! 0.25 T8^{ac}_{nm}G3^{cnam}_{}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("T(a,c,n,m)+=Z(a,b,k,n)*Y(b,c,k,m)",I1,t2_o,t2_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction T8 wrong')
      ierr=talsh_tensor_contract("G(a,c,n,m)+=X+(j,i,n,c)*Y+(j,i,m,a)",I2,ooov_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction G3 wrong')
      ierr=talsh_tensor_contract("R()+=T(a,c,n,m)*G(a,c,n,m)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 15 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(15)=real(T_tens(1))

      !----------------------------------------------------------------------------------------------------

      ! 0.25 I5^{bm}_{nj}I15^{bn}_{mj}
      dims4=nocc
      dims4(1)=nvir
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(b,m,n,j)+=X(a,b,i,n)*Y+(j,i,m,a)",I1,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I5 wrong')
      ierr=talsh_tensor_contract("N(b,m,n,j)+=X(b,c,k,m)*Y+(k,j,n,c)",I2,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I15 wrong')
      ierr=talsh_tensor_contract("R()+=M(b,m,n,j)*N(b,m,n,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 16 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(16)=real(T_tens(1))

      ! 0.25 I8^{bm}_{ni}I12^{bn}_{mi}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(b,m,n,i)+=X(a,b,j,n)*Y+(j,i,m,a)",I1,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I8 wrong')
      ierr=talsh_tensor_contract("N(b,m,n,i)+=X(b,c,k,m)*Y+(i,k,n,c)",I2,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I12 wrong')
      ierr=talsh_tensor_contract("R()+=M(b,m,n,i)*N(b,m,n,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 17 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(17)=real(T_tens(1))

      ! 0.25 I8^{cn}_{mj}I12^{cm}_{nj}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(c,n,m,j)+=X(c,a,i,n)*Y+(j,i,m,a)",I1,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I8 wrong')
      ierr=talsh_tensor_contract("N(c,n,m,j)+=X(b,c,k,m)*Y+(k,j,n,b)",I2,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I12 wrong')
      ierr=talsh_tensor_contract("R()+=M(c,n,m,j)*N(c,n,m,j)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 18 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(18)=real(T_tens(1))

      ! 0.25 I15^{cm}_{ni}I5^{cm}_{nj}
      T_tens(1)=ZERO
      ierr=talsh_tensor_construct(I1,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_construct(I2,C8,dims4,init_val=ZERO)
      ierr=talsh_tensor_contract("M(c,m,n,i)+=X(c,a,j,n)*Y+(j,i,m,a)",I1,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I15 wrong')
      ierr=talsh_tensor_contract("N(c,m,n,i)+=X(b,c,k,m)*Y+(i,k,n,b)",I2,t2_o,ooov_o)
      if (ierr.ne.0) call quit('laplace_get_T: contraction I5 wrong')
      ierr=talsh_tensor_contract("R()+=M(c,m,n,i)*N(c,m,n,i)",T_tensor,I1,I2,scale=ONE_QUARTER)
      if (ierr.ne.0) call quit('laplace_get_T: contraction 19 wrong')
      ierr=talsh_tensor_destruct(I1)
      ierr=talsh_tensor_destruct(I2)
      T(19)=real(T_tens(1))

      write(*,*) 'Terms:'
      write(*,'(5e12.4)') T
      

      !clean up
      ierr=talsh_tensor_destruct(T_tensor)
      ierr=talsh_tensor_destruct(t2_o)
      ierr=talsh_tensor_destruct(ooov_o)

    end subroutine laplace_get_T


end module talsh_cc
