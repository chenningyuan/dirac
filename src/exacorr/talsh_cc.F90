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
                                 ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                                 MINUS_ONE_HALF=(-0.5D0,0.D0), TWO=(2.D0,0.D0), &
                                 MINUS_TWO=(-2.D0,0.D0), THREE=(3.0,0.0), SIX=(6.0,0.0)
        real(8), parameter    :: ONE_QUARTER=0.25D0

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
         use talsh_ao_to_mo, only       : tensor_norm2, compute_nelectrons
         use talsh_gradient
         use talsh_eom,      only       : talsh_eom_results_tens_t, talsh_eom_eigval_driver

         type(exacc_input), intent(in) :: exa_input

         integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer(INTD)       :: ierr
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
         integer(C_INT):: host_arg_max

!        Common tensors
         type(talsh_comm_tens_t) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_intg_tens_t) :: int_t
!        Lambda tensors
         type(talsh_tens_t) :: l1_tensor, l2_tensor
!        1- and 2-body cc intermediates 
         type(talsh_intermediates_tens_t) :: interm_t
!        eom results
         type(talsh_eom_results_tens_t) :: eom_results_t
!        intermediate tensors (same notation as in relccsd)
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims, oovv_dims, vvoo_dims, vovo_dims
!        scalars (need to be defined as tensor types)
         type(talsh_tens_t) :: result_tensor
         integer(INTD)      :: one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CCSD control and result variables
         real(8) :: mp2_energy, ccsd_energy, t_energy(3)
         real(8) :: t_target_precision
         integer :: ncycles, print_level

         logical :: CCD, lambda, eom
         character(len=2) :: eom_flavour
         logical :: debug_eom = .false.

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
         lambda = exa_input%lambda
         eom = exa_input%eom
         eom_flavour = exa_input%eom_flavour
         ncycles = exa_input%ncycles
         print_level = exa_input%print_level

         write(*,*) ""
         write(*,'(A)') " Characteristics of this calculation"
         write(*,'(A)') " -----------------------------------"
         write(*,'(A,L1)') " - Solve CCD equations :              ", exa_input%ccd
         write(*,'(A,L1)') " - Solve Lambda equations :           ", exa_input%lambda
         write(*,'(A,I3)') " - Number of atomic orbitals :      ", get_nao()
         write(*,'(A,I3)') " - Number of occupied spinors :     ", exa_input%nocc
         write(*,'(A,I3)') " - Number of virtual spinors :      ", exa_input%nvir
         write(*,'(A,ES10.1)') " - Target precision of amplitudes :", exa_input%t_econv
         write(*,'(A,I3)') " - Maximum number of cycles :       ", exa_input%ncycles
         write(*,'(A,I3)') " - Integral transformation scheme  ",exa_input%moint_scheme
         if (exa_input%moint_scheme.eq.42) then
            write(*,'(A,ES10.1)') " - Cholesky decomposition threshold ",exa_input%t_cholesky
         end if
         if (exa_input%level_shift.gt.1.0D-14) then
           write(*,'(A,ES10.1)') " - Level Shift :                     ",exa_input%level_shift
         end if
         write(*,'(A,I3)') " - Print Level :                    ",exa_input%print_level
         if (print_level > 1) then
           write(*,*) " - occupied spinors :                ", exa_input%mo_occ
           write(*,*) " - virtual spinors :                 ", exa_input%mo_vir
         end if
         write(*,*) " - buffer size: ",exa_input%talsh_buff,"GB"
         write(*,'(A)') ""

         call print_date("Getting orbital energies")
!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("Retrieved orbital energies")

         print *,'before get_CC_integrals'
         print *,' eps_occ:',eps_occ
         print *,' eps_vir:',eps_vir
!        Get fixed two-body tensors
         if (exa_input%moint_scheme.eq.42) then
           call get_CC_integrals_chol (nocc,nvir,mo_occ,mo_vir,int_t, exa_input%t_cholesky)
         else
           call get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t)
         end if

         print *,'before get_CC_fockmatrix'
         print *,' eps_occ:',eps_occ
         print *,' eps_vir:',eps_vir
!        Get fixed one-body tensors, needs to follow the two body tensors
         call get_CC_fockmatrix (nocc,nvir,mo_occ,mo_vir,eps_occ,eps_vir, &
                                 comm_t%foo,comm_t%fov,comm_t%fvv, int_t, print_level, exa_input%level_shift)

!        Initialize t1 amplitudes to zero (we assume canonical orbitals for the moment)
         if (.not.CCD) then
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=talsh_tensor_construct(comm_t%t1,C8,vo_dims,init_val=ZERO)
         end if

         print *,'before scale_with_denominators'
         print *,' eps_occ:',eps_occ
         print *,' eps_vir:',eps_vir
!        Initialize t2 amplitudes with the MP2 values
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(comm_t%t2,C8,vvoo_dims,init_val=ZERO)
!        ierr=talsh_tensor_add("T(a,b,i,j)+=V+(i,j,a,b)",t2_tensor,oovv_tensor)
!        work-around because complex conjugation does not work in the above statement
         ierr=talsh_tensor_contract("T(a,b,i,j)+=V+(i,j,a,b)",comm_t%t2,int_t%oovv,one_tensor) ! o(2)v(2)
!        end of work-around
         if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)
         else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
         end if

!        Calculate and print the MP2 energy
         ierr=talsh_tensor_init(result_tensor)
         ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv) ! o(2)v(2)
         mp2_energy=result_tens(1) * ONE_QUARTER
         print*,""
         print*," MP2 energy = ",mp2_energy
         print*,""

!=====================
!        Calculate the pair correlation energy for each pair of occupied spin-orbits
         call pair_correlation_energies(comm_t%t2, int_t%oovv, nocc, nvir)
!=====================

!        Start CCSD iterations

         oo_dims = nocc
         ierr=talsh_tensor_construct(comm_t%goo,C8,oo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%hoo,C8,oo_dims,init_val=ZERO)
         vv_dims = nvir
         ierr=talsh_tensor_construct(comm_t%gvv,C8,vv_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%hvv,C8,vv_dims,init_val=ZERO)
         if (.not.CCD) then
            ov_dims(1) = nocc
            ov_dims(2) = nvir
            ierr=talsh_tensor_construct(comm_t%hov,C8,ov_dims,init_val=ZERO)
         end if
         oooo_dims  = nocc
         ierr=talsh_tensor_construct(comm_t%a_int,C8,oooo_dims,init_val=ZERO)
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=talsh_tensor_construct(comm_t%h_int,C8,vovo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(comm_t%tau,C8,vvoo_dims,init_val=ZERO)
         if (.not.CCD) then
            call diis (comm_t%t2,comm_t%t1,mode=1)
         else
            call diis (comm_t%t2,mode=1)
         end if

         if (ierr.eq.TALSH_SUCCESS) then
            call print_date('Succesfully allocated all tensors')
         else
            print*," ccdriver: error in tensor allocation, code=",ierr
            return
         end if

         if (debug_exa) then
            print*, "Voovv = ",tensor_norm2(int_t%oovv)
         end if

         print *,'before solve_t_equations'
         print *,' eps_occ:',eps_occ
         print *,' eps_vir:',eps_vir
         call solve_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                ccsd_energy,CCD,ncycles)

         if (.not.CCD) then
!           Calculate (T) energy correction
            call t3_energy (t_energy,eps_occ,eps_vir, comm_t%t1,comm_t%t2,comm_t%fov,int_t%ooov,int_t%vovv,int_t%oovv)
            write (*,'(A,F20.15)') "  Final CCSD energy :               ",ccsd_energy
            write (*,'(A,F20.15)') "  4th order triples correction :    ",t_energy(1)
            write (*,'(A,F20.15)') "  5th order triples (T) correction :",t_energy(2)
            write (*,'(A,F20.15)') "  5th order triples -T  correction :",t_energy(3)
            write (*,'(A,F20.15)') "  CCSD+T  correlation energy :      ",ccsd_energy+t_energy(1)
            write (*,'(A,F20.15)') "  CCSD(T) correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(2)
            write (*,'(A,F20.15)') "  CCSD-T  correlation energy :      ",ccsd_energy+t_energy(1)+t_energy(3)
         else
            write (*,'(A,F20.15)') "  Final CCD energy     :            ",ccsd_energy
         end if

         !Destruct diis tensors
         if (.not.CCD) then
            call diis (comm_t%t2,comm_t%t1,mode=2)
         else
            call diis (comm_t%t2,mode=2)
         end if

         if (lambda) then

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
            call solve_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir, &
                                        eps_occ,eps_vir,t_target_precision,ncycles,l1_tensor)
            call talsh_dm_driver(exa_input,comm_t%t2,l2_tensor,comm_t%t1,l1_tensor)
          else
            call solve_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir, &
                                        eps_occ,eps_vir,t_target_precision,ncycles)
            call talsh_dm_driver(exa_input,comm_t%t2,l2_tensor)
          end if

         end if

!        Destruct tensors
         ierr=talsh_tensor_destruct(comm_t%tau)
         ierr=talsh_tensor_destruct(comm_t%hoo)
         if (.not.CCD) ierr=talsh_tensor_destruct(comm_t%hov)
         ierr=talsh_tensor_destruct(comm_t%hvv)
         ierr=talsh_tensor_destruct(comm_t%goo)
         ierr=talsh_tensor_destruct(comm_t%gvv)
         ierr=talsh_tensor_destruct(comm_t%a_int)
         ierr=talsh_tensor_destruct(comm_t%h_int)
         ierr=talsh_tensor_destruct(result_tensor)
         ierr=talsh_tensor_destruct(one_tensor)

!**************************
!        EOM calculations *
!**************************

         if (eom) then

!        Get fixed one-body tensors, needs to follow the two body tensors
            call get_CC_fockmatrix (nocc,nvir,mo_occ,mo_vir,eps_occ,eps_vir, &
                                 comm_t%foo,comm_t%fov,comm_t%fvv, int_t, print_level, exa_input%level_shift,keep_diagonal=.true.)
! create eom intermediates; actually, after the ground-state cc, and/or before lambda, we should create the intermediates and have them around (perhaps lambda and eom intermediates is too much to keep.
            call form_ccsd_intermediates(interm_t, exa_input%CCD, exa_input%nocc, exa_input%nvir, comm_t, int_t)
         
! with the ground-state intermediates done, and t amplitudes store

          if (debug_eom) then
             print*, "foo_tensor norm before eom #4", tensor_norm2(comm_t%foo)
             print*, "fov_tensor norm before eom #4", tensor_norm2(comm_t%fov)
             print*, "fvv_tensor norm before eom #4", tensor_norm2(comm_t%fvv)
          end if

           if (exa_input%eom_flavour.eq."EE" .and. exa_input%eom_properties .and. .not.lambda) then
             call quit("To do EOM-EE properties, you first need to calculate the Lambda amplitudes!")
           end if

! create eom intermediates

! solve eom eigenvalue equations; right now we do only one flavor but there's no reason not to handle different flavors in the same calculation,
! we just need different results instances
           call talsh_eom_eigval_driver(eom_results_t,exa_input,eps_occ,eps_vir,interm_t,int_t)
! TODO: if we want properties etc, now we call other drivers

         end if

         if (.not.CCD) then
           ierr=talsh_tensor_destruct(comm_t%t1)
           ierr=talsh_tensor_destruct(l1_tensor)
         end if
         ierr=talsh_tensor_destruct(comm_t%t2)
         ierr=talsh_tensor_destruct(l2_tensor)
! aspg
         ierr=talsh_tensor_destruct(comm_t%foo)
         ierr=talsh_tensor_destruct(comm_t%fov)
         ierr=talsh_tensor_destruct(comm_t%fvv)
! aspg
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
         print*, ""
         call print_date('Leaving cc_driver routine')

        end subroutine talsh_cc_driver

        subroutine solve_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                ccsd_energy,CCD,ncycles)

!        Routine for solving the t amplitude equations

         use talsh_ao_to_mo, only : tensor_norm2
         use exacorr_global, only : debug_exa

         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_comm_tens_t), intent(inout) :: comm_t
         type(talsh_intg_tens_t), intent(inout) :: int_t
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
         integer (INTD)       :: ierr

!        Initialize scalars that are to be used as tensors in contractions
         result_dims = 1
         ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
         ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,result_tens,result_dims)

         if (.not.CCD) then
            oo_dims = nocc
            ierr=talsh_tensor_construct(hoo_aux_tensor,C8,oo_dims,init_val=ZERO)
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=talsh_tensor_construct(s1_tensor,C8,vo_dims,init_val=ZERO)
         end if
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         ierr=talsh_tensor_construct(hov_aux_tensor,C8,ov_dims,init_val=ZERO)
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(s2_tensor,C8,vvoo_dims,init_val=ZERO)
         ierr=talsh_tensor_construct(tau2_tensor,C8,vvoo_dims,init_val=ZERO)
         ovoo_dims(1) = nocc
         ovoo_dims(2) = nvir
         ovoo_dims(3:4) = nocc
         ierr=talsh_tensor_construct(cint_tensor,C8,ovoo_dims,init_val=ZERO)

         do iccsd = 1, ncycles

          if (.not.CCD) then
              call get_tau(comm_t%tau,comm_t%t2,comm_t%t1)
          else
              call get_tau(comm_t%tau,comm_t%t2)
          end if

          if (ierr.ne.TALSH_SUCCESS) then
             print*," ccdriver: error in calculation of common blocks, code=",ierr
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
          call get_hov(comm_t,int_t%oovv,CCD)
          ierr=talsh_tensor_init(hov_aux_tensor)
          ierr=talsh_tensor_add("L(k,c)+=H(k,c)",hov_aux_tensor,comm_t%hov)                       ! copy common term
          ierr=talsh_tensor_add("H(k,c)+=F(k,c)",hov_aux_tensor,comm_t%fov,scale=MINUS_TWO)       ! subtract fock matrix twice

!         oooo ! validated
          call get_a(comm_t,int_t,CCD)

          if (ierr.ne.TALSH_SUCCESS) then
             print*," ccdriver: error in calculation of intermediates, code=",ierr
             return
          end if

          if (debug_exa) then
             if (.not.CCD) print*, "T1amp  = ",tensor_norm2(comm_t%t1)
             print*, "T2amp  = ",tensor_norm2(comm_t%t2)*ONE_QUARTER
             print*, "Tau    = ",tensor_norm2(comm_t%tau)*ONE_QUARTER
             print*, "hoint  = ",tensor_norm2(comm_t%hoo)
             print*, "goint  = ",tensor_norm2(comm_t%goo)
             print*, "hvint  = ",tensor_norm2(comm_t%hvv)
             print*, "gvint  = ",tensor_norm2(comm_t%gvv)
             if (.not.CCD) print*, "hovint = ",tensor_norm2(comm_t%hov)
             print*, "a int  = ",tensor_norm2(comm_t%a_int)*ONE_QUARTER
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

          if (ierr.ne.TALSH_SUCCESS) then
             print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
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
          if (.not.CCD) then
              call get_tau(comm_t%tau,comm_t%t2,comm_t%t1)
          else
              call get_tau(comm_t%tau,comm_t%t2)
          end if

          ierr=talsh_tensor_init(result_tensor)
          ierr=talsh_tensor_contract("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER

          call print_iteration(iccsd,t_convergence,ccsd_energy)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
         write(*,*) ""

         if (.not.CCD) ierr=talsh_tensor_destruct(s1_tensor)
         ierr=talsh_tensor_destruct(s2_tensor)
         ierr=talsh_tensor_destruct(tau2_tensor)
         ierr=talsh_tensor_destruct(cint_tensor)
         if (.not.CCD) then
            ierr=talsh_tensor_destruct(hoo_aux_tensor)
         end if
         ierr=talsh_tensor_destruct(hov_aux_tensor)
         ierr=talsh_tensor_destruct(result_tensor)

        end subroutine

        subroutine solve_lambda_equations(comm_t,int_t,l2_tensor,nocc,nvir, &
                                            eps_occ,eps_vir,t_target_precision,ncycles,l1_tensor)

!        Routine for solving the lambda equations

         use talsh_ao_to_mo, only : tensor_norm2
         use exacorr_global, only : debug_exa

         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(talsh_comm_tens_t), intent(inout) :: comm_t
         type(talsh_intg_tens_t), intent(inout) :: int_t
         type(talsh_lambda_tens_t)              :: lambda_t
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

         if (.not.CCD) then
              call get_tau(comm_t%tau,comm_t%t2,comm_t%t1)
         else
              call get_tau(comm_t%tau,comm_t%t2)
         end if
         call get_hvv(comm_t,int_t%oovv)
         call get_gvv(comm_t,int_t%vovv,CCD)
         call get_h(comm_t,int_t,CCD)
         call get_hoo(comm_t,int_t%oovv)
         call get_goo(comm_t,int_t%ooov,CCD)
         call get_hov(comm_t,int_t%oovv,CCD)
         call get_a(comm_t,int_t,CCD)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t)
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         if (.not.CCD) then
           call get_w_vovv(lambda_t%w_vovv,int_t,comm_t%t1)
         else
           call get_w_vovv(lambda_t%w_vovv,int_t)
         end if
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,nocc,nvir,CCD)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         if (.not.CCD) then
           call get_w_ooov(lambda_t%w_ooov,int_t,comm_t%t1)
         else
           call get_w_ooov(lambda_t%w_ooov,int_t)
         end if
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD)
         call print_date('Finished calculating fixed intermediates')

         if (debug_exa) then
             if (.not.CCD) print*, "T1amp      = ", tensor_norm2(comm_t%t1)
             print*, "T2amp      = ", tensor_norm2(comm_t%t2)*ONE_QUARTER
             ! print*, "Tau        = ", tensor_norm2(comm_t%tau)*ONE_QUARTER
             ! print*, "Fbar_oo    = ", tensor_norm2(lambda_t%fbar_oo)
             print*, "Fbar_ov    = ", tensor_norm2(lambda_t%fbar_ov)
             ! print*, "Fbar_vv    = ", tensor_norm2(lambda_t%fbar_vv)
             ! print*, "w_vvvv     = ", tensor_norm2(lambda_t%w_vvvv)*ONE_QUARTER
             print*, "wbar_vovo  = ", tensor_norm2(lambda_t%wbar_vovo)
             print*, "w_vvvo     = ", tensor_norm2(lambda_t%w_vvvo)*0.5D0
             print*, "w_vovv     = ", tensor_norm2(lambda_t%w_vovv)*0.5D0
             print*, "w_vovo     = ", tensor_norm2(lambda_t%w_vovo)
             print*, "w_ovoo     = ", tensor_norm2(lambda_t%w_ovoo)*0.5D0
             print*, "w_ooov     = ", tensor_norm2(lambda_t%w_ooov)*0.5D0
             print*, "w_oooo     = ", tensor_norm2(lambda_t%w_oooo)*ONE_QUARTER
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
          if (debug_exa) then
             if (.not.CCD) print*, "L1  = ", tensor_norm2(l1_tensor)
             print*, "L2  = ", tensor_norm2(l2_tensor)*ONE_QUARTER
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
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(ij,mb) |
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
             print*," ccdriver: error in evaluation of lambda equations, code=",ierr
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

         if (t_convergence .gt. t_target_precision) print*, "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", tensor_norm2(l1_tensor)
         write(*,'(A,ES22.16)') "  Final L2amp = ", tensor_norm2(l2_tensor)*ONE_QUARTER
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

        subroutine get_CC_fockmatrix (nocc,nvir,mo_occ,mo_vir,eps_occ,eps_vir, &
                                      fo_tensor,fov_tensor,fv_tensor, int_t, print_level, level_shift, keep_diagonal)

!        Written by Lucas Visscher, summer 2017

         use exacorr_global

         implicit none

         integer, intent(in)                  :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)                  :: print_level
         integer, intent(in)                  :: mo_occ(:),mo_vir(:) ! the orbital energies
         real(8), intent(inout)               :: eps_occ(:),eps_vir(:) ! orbital energies
         real(8), intent(in)                  :: level_shift
         logical, optional, intent(in)        :: keep_diagonal
!        target tensors to be filled and given back to the caller
         type(talsh_tens_t), intent(inout)    :: fo_tensor, fv_tensor, fov_tensor
         type(talsh_intg_tens_t), intent(inout) :: int_t

         type(talsh_tens_t)              :: id_mat_occ, diag_mat_occ, diag_mat_vir
         type(talsh_tens_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)
         integer(INTD)                   :: ierr
         logical                         :: one_el_exist
         type(one_el_t)                  :: one_int
         real(8)                         :: e_core, ev_energy, off_diag, hf_energy
         integer                         :: nkp, min_occ
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer                         :: i, j, iv, jv, imo, jmo
         type(talsh_tens_t)              :: result_tensor
         integer(INTD)                   :: result_dims(1)
         complex(8), pointer             :: result_tens(:)

         real(8)            :: e_lumo, e_homo
         integer            :: lumo
         real(8), parameter :: THRESHOLD=1.0D-8
         logical            :: keep_diagonal_local

         if (present(keep_diagonal)) then
           keep_diagonal_local = keep_diagonal
         else
           keep_diagonal_local = .false.
         end if

         min_occ=minval(mo_occ)
         result_dims(1)=1

!        Get foo tensor
         fo_dims(1:2) = nocc
         ierr=talsh_tensor_construct(fo_tensor,C8,fo_dims,init_val=ZERO)

!        Get fvv tensor
         fv_dims(1:2) = nvir
         ierr=talsh_tensor_construct(fv_tensor,C8,fv_dims,init_val=ZERO)

!        Get fov tensor
         fov_dims(1) = nocc
         fov_dims(2) = nvir
         ierr=talsh_tensor_construct(fov_tensor,C8,fov_dims,init_val=ZERO)

!        get one electron integrals
         one_el_exist = exist_one_el(e_core)

!        write one electron integrals to tensor
         if (one_el_exist) then

!           create density (identity) matrix for the occupied space
            ierr=talsh_tensor_construct(id_mat_occ,C8,fo_dims,init_val=ZERO)
            ierr=talsh_tensor_get_body_access(id_mat_occ,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, nocc
              h_tens(i,i) = ONE
            enddo

!           now write one electron integrals

            ierr=talsh_tensor_get_body_access(fo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_occ,nocc,min_occ)


            ierr=talsh_tensor_get_body_access(fv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_vir,nvir,mo_vir,nvir,min_occ)

            ierr=talsh_tensor_get_body_access(fov_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fov_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_vir,nvir,min_occ)

!           Get access to result
            result_dims(1) = 1
            ierr=talsh_tensor_construct(result_tensor,C8,result_dims(1:0),init_val=ZERO)
            ierr=talsh_tensor_get_body_access(result_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,result_tens,result_dims)

!           compute HF energy
            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("R()=F(i,j)*I(i,j)",result_tensor,fo_tensor,id_mat_occ)
            ierr=talsh_tensor_construct(h_tensor,C8,fo_dims,init_val=ZERO)
            ierr=talsh_tensor_contract("H(i,j)+=G(i,k,j,l)*I(k,l)",h_tensor,int_t%oooo,id_mat_occ, ONE_HALF)
            ierr=talsh_tensor_contract("R()+=H(i,j)*I(i,j)",result_tensor,h_tensor,id_mat_occ)
            ierr=talsh_tensor_destruct(h_tensor)

            hf_energy=real(result_tens(1),8) ! Casting into real(8)
            write (*,'(A,F23.15)') " Energy of reference wave function (HF single determinant) = ", hf_energy+e_core

            call print_date('one electron integrals added to Fock Matrix')

            ierr=talsh_tensor_contract("F(i,j)+=V(i,k,j,l)*I(k,l)",fo_tensor,int_t%oooo,id_mat_occ, ONE)

            ierr=talsh_tensor_contract("F(i,a)+=V(k,i,l,a)*I(k,l)",fov_tensor,int_t%ooov,id_mat_occ, ONE)

            ierr=talsh_tensor_contract("F(a,b)+=V(a,k,b,l)*I(k,l)",fv_tensor,int_t%vovo,id_mat_occ, ONE)
            call print_date('two electron integrals added to Fock Matrix')

            if (print_level >3) call print_tensor(fo_tensor, 1.0D-6, 'fo_tensor')
            if (print_level >3) call print_tensor(fov_tensor, 1.0D-6, 'fov_tensor')
            if (print_level >3) call print_tensor(fv_tensor, 1.0D-6, 'fv_tensor')

            !now we can get the correct epsilon values
            ierr=talsh_tensor_get_body_access(fo_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            do i=1,fo_dims(1)
              eps_occ(i)=h_tens(i,i)
            end do

            ierr=talsh_tensor_get_body_access(fv_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            do i=1,fv_dims(1)
              eps_vir(i)=h_tens(i,i)
            end do

            !apply level shit
            if (level_shift.gt.1.0D-14) then
              e_homo = maxval(eps_occ)
              e_lumo = minval(eps_vir)
              lumo   = minloc(eps_vir,1)

              do while(e_lumo-e_homo < THRESHOLD)
                if (e_lumo+level_shift-e_homo >THRESHOLD) then
                  eps_vir(lumo)=e_lumo+level_shift
                  print *, "lumo energy level shifted by ", level_shift," to ",e_lumo+level_shift
                else
                  call quit('Negative HOMO-LUMO gap')
                end if
                e_lumo = minval(eps_vir)
                lumo   = minloc(eps_vir, 1)
              end do
            end if

            !set diag_mat
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

            if (.not.keep_diagonal_local) then
              ierr=talsh_tensor_add("F(i,j)+=D(i,j)",fo_tensor,diag_mat_occ,scale=MINUS_ONE)
              ierr=talsh_tensor_add("F(a,b)+=D(a,b)",fv_tensor,diag_mat_vir,scale=MINUS_ONE)
              if (print_level >1) call print_date('diagonal substracted from Fock Matrix')
            end if

            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("R()+=F(i,j)*F(i,j)",result_tensor,fo_tensor,fo_tensor)
            off_diag=real(result_tens(1),8)
            if (print_level >0) print*," Norm of OO part Fock matrix (without diag) = ", off_diag

            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("R()+=F(i,a)*F(i,a)",result_tensor,fov_tensor,fov_tensor)
            off_diag=real(result_tens(1),8)
            if (print_level >0) print*," Norm of OV part Fock matrix =                ", off_diag

            ierr=talsh_tensor_init(result_tensor)
            ierr=talsh_tensor_contract("R()+=F(a,b)*F(a,b)",result_tensor,fv_tensor,fv_tensor)
            off_diag=real(result_tens(1),8)
            if (print_level >0) print*," Norm of VV part Fock matrix (without diag) = ", off_diag

            ierr=talsh_tensor_destruct(id_mat_occ)
            ierr=talsh_tensor_destruct(diag_mat_occ)
            ierr=talsh_tensor_destruct(diag_mat_vir)
            ierr=talsh_tensor_destruct(result_tensor)
         else
            if (level_shift.gt.1.0D-14) then
              e_homo = maxval(eps_occ)
              e_lumo = minval(eps_vir)
              lumo   = minloc(eps_vir,1)

              do while(e_lumo-e_homo < THRESHOLD)
                print*, "e_homo, e_lumo", e_homo, e_lumo
                if (e_lumo+level_shift-e_homo >THRESHOLD) then
                  eps_vir(lumo)=e_lumo+level_shift
                  print *, "lumo energy level shifted by ", level_shift," to ",e_lumo+level_shift
                else
                  call quit('Negative HOMO-LUMO gap')
                end if
                e_lumo = minval(eps_vir)
                lumo   = minloc(eps_vir, 1)
              end do
            end if
         end if
         if (print_level >0) call print_date('Fock matrix initialized')
        end subroutine get_CC_Fockmatrix

        subroutine get_CC_integrals (nocc,nvir,mo_occ,mo_vir,int_t)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list.

!        Written by Lucas Visscher, summer 2017

         use exacorr_global
         use talsh_ao_to_mo

         integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(talsh_intg_tens_t),intent(inout)    :: int_t
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

         call print_date('entering integral transformation, scheme = 2')

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get tensor with AO integrals
         aoint_dims      = nao
         ierr=talsh_tensor_construct(aoint_tensor,C8,aoint_dims,init_val=ZERO)
         call print_date('allocated AO integral tensor')
!        Initialize AO integral tensors (ideally via an initialization routine, for now via direct access)
         ierr=talsh_tensor_get_body_access(aoint_tensor,body_p,C8,int(0,C_INT),DEV_HOST)
         call c_f_pointer(body_p,ao_tens,aoint_dims(1:4)) ! to use <ao_tens> as a regular Fortran 4d array
         call compute_ao_integrals (nao,ao_tens)
         call print_date('finished ao integral calculation')

!        Get antisymmetric oovv tensor
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
!        whereas the dimensions for the mo lists in the index transform are in mullikens notation (13||24)
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
         call print_date('finished oovv integrals')

!        Get antisymmetric oooo tensor
         oooo_dims = nocc
         ierr=talsh_tensor_construct(int_t%oooo,C8,oooo_dims,init_val=ZERO)
         allocate (mo_list(4*nocc))
         nmo = nocc
         mo_list(1:nocc)          = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)   = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc) = mo_occ(1:nocc)
         mo_list(3*nocc+1:4*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,int_t%oooo,12,nmo,mo_list)
         deallocate (mo_list)
         call print_date('finished oooo integrals')

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
!        Get vovo tensor
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=talsh_tensor_construct(int_t%vovo,C8,vovo_dims,init_val=ZERO)
         allocate (mo_list(2*nvir+2*nocc))
         nmo(1:2) = nvir
         nmo(3:4) = nocc
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)               = mo_vir(1:nvir)
         mo_list(2*nvir+1:2*nvir+nocc)        = mo_occ(1:nocc)
         mo_list(2*nvir+nocc+1:2*nvir+2*nocc) = mo_occ(1:nocc)
         call get_integral_tensor (aoint_tensor,int_t%vovo,0,nmo,mo_list)
         deallocate (mo_list)
!        Get ovvo tensor
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
!        Make anti-symmetrized vovo tensor
         minusone_dims(1) = 1
         ierr=talsh_tensor_construct(minusone_tensor,C8,minusone_dims(1:0),init_val=MINUS_ONE)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,minusone_tensor)
         ierr=talsh_tensor_destruct(ovvo_tensor)
         call print_date('finished vovo integrals')

!        Get vvvv tensor
         vvvv_dims = nvir
         ierr=talsh_tensor_construct(int_t%vvvv,C8,vvvv_dims,init_val=ZERO)
         allocate (mo_list(4*nvir))
         nmo = nvir
         mo_list(1:nvir)          = mo_vir(1:nvir)
         mo_list(nvir+1:2*nvir)   = mo_vir(1:nvir)
         mo_list(2*nvir+1:3*nvir) = mo_vir(1:nvir)
         mo_list(3*nvir+1:4*nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,int_t%vvvv,12,nmo,mo_list)
         deallocate (mo_list)
         call print_date('finished vvvv integrals')

!        Get ooov tensor
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=talsh_tensor_construct(int_t%ooov,C8,ooov_dims,init_val=ZERO)
         allocate (mo_list(3*nocc+nvir))
         nmo(1:3) = nocc
         nmo(4)   = nvir
         mo_list(1:nocc)               = mo_occ(1:nocc)
         mo_list(nocc+1:2*nocc)        = mo_occ(1:nocc)
         mo_list(2*nocc+1:3*nocc)      = mo_occ(1:nocc)
         mo_list(3*nocc+1:3*nocc+nvir) = mo_vir(1:nvir)
         call get_integral_tensor (aoint_tensor,int_t%ooov,12,nmo,mo_list)
         deallocate (mo_list)
         call print_date('finished ooov integrals')

!        Get vovv tensor
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(int_t%vovv,C8,vovv_dims,init_val=ZERO)
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
         call print_date('finished vovv integrals')

!        Clean-up
         ierr=talsh_tensor_destruct(aoint_tensor)

       end subroutine get_CC_integrals


       subroutine get_CC_integrals_chol (nocc,nvir,mo_occ,mo_vir,int_t,t_cholesky)

!        Routine to get integrals needed in CC in the form of antisymmetric tensors.
!        In this implementation we take all AO integrals into memory, the result tensor themselves can be subsets of the
!        full integral list. We do assume that we deal with a closed shell (Kramers-restricted) case, so nocc and nvir
!        refer to Kramers pairs rather than spinors.

!        Written by Johann Pototschnig, summer 2019

         use exacorr_global
         use talsh_ao_to_mo
         use exacorr_cholesky

         implicit none

         integer, intent(in)    :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer, intent(in)    :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals
!        target integral tensors to be filled and given back to the caller
         type(talsh_intg_tens_t),intent(inout)    :: int_t

         real(8), intent(in)                  :: t_cholesky ! cholesky threshold

         integer, dimension(4)  :: oovv_dims,vvvv_dims,oooo_dims,vovo_dims,ooov_dims,vovv_dims

         integer                            :: ierr
         real(8), pointer, contiguous       :: ao_tens(:,:,:,:)
         type(C_PTR)                        :: body_p

         integer :: nao

!        auxilliary integral tensors (only needed inside this routine)
         type(talsh_tens_t)    :: aoint_tensor, ovvo_tensor
         integer, dimension(4) :: aoint_dims, ovvo_dims

!        arrays needed to communicate with the MO integral generator
         integer:: nmo(4)
         integer, allocatable :: mo_list(:)

!        cholesky
         real(8), allocatable            :: CD(:,:,:)
         complex(8), pointer, contiguous :: chol_tens(:,:,:)
         type(talsh_tens_t)              :: chol_tensor, test_tensor, err_tensor, CD0
         type(talsh_tens_t)              :: uoo_tensor, uov_tensor, uvo_tensor, uvv_tensor
         integer, dimension(3)           :: chol_dims, uoo_dims, uov_dims, uvo_dims, uvv_dims
         integer, dimension(4)           :: ovov_dims, vvov_dims
         type(talsh_tens_t)              :: result_tensor
         integer                         :: result_dims(1)
         complex(8), pointer             :: result_tens(:)
         integer                         :: i,j,m,ij

         call print_date('entering integral transformation, scheme = 42')
         print *, "Note: cholesky based ao to mo transformation fails on iMac16,1"

!        Retrieve basis set information
         nao     = get_nao()      ! number of basis functions

!        Get cholesky-vectors and copy them in to a talsh tensor
         call decompose_cholesky_talsh (chol_tensor, t_cholesky, nao)
         ierr = talsh_tensor_dimensions(chol_tensor,i,chol_dims)

!        perform ao to mo transformation for cholesky tensor
         uoo_dims(1:2)=nocc
         uoo_dims(3)=chol_dims(3)
         ierr=talsh_tensor_construct(uoo_tensor,C8,uoo_dims,init_val=ZERO)
         allocate(mo_list(2*nocc))
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nocc)            = mo_occ(1:nocc)
         call ao2mo_2ind (chol_tensor,uoo_tensor,uoo_dims(1:2),mo_list)
         deallocate(mo_list)

         uov_dims(1)=nocc
         uov_dims(2)=nvir
         uov_dims(3)=chol_dims(3)
         ierr=talsh_tensor_construct(uov_tensor,C8,uov_dims,init_val=ZERO)
         allocate(mo_list(nocc+nvir))
         mo_list(1:nocc)                      = mo_occ(1:nocc)
         mo_list(nocc+1:nocc+nvir)            = mo_vir(1:nvir)
         call ao2mo_2ind (chol_tensor,uov_tensor,uov_dims(1:2),mo_list)
         deallocate(mo_list)

         uvo_dims(1)=nvir
         uvo_dims(2)=nocc
         uvo_dims(3)=chol_dims(3)
         ierr=talsh_tensor_construct(uvo_tensor,C8,uvo_dims,init_val=ZERO)
         allocate(mo_list(nocc+nvir))
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:nvir+nocc)            = mo_occ(1:nocc)
         call ao2mo_2ind (chol_tensor,uvo_tensor,uvo_dims(1:2),mo_list)
         deallocate(mo_list)

         uvv_dims(1:2)=nvir
         uvv_dims(3)=chol_dims(3)
         ierr=talsh_tensor_construct(uvv_tensor,C8,uvv_dims,init_val=ZERO)
         allocate(mo_list(2*nvir))
         mo_list(1:nvir)                      = mo_vir(1:nvir)
         mo_list(nvir+1:nvir+nvir)            = mo_vir(1:nvir)
         call ao2mo_2ind (chol_tensor,uvv_tensor,uvv_dims(1:2),mo_list)
         deallocate(mo_list)

!        Get antisymmetric oovv tensor
!        Note that the dimensions for the tensor are given in physicists notation <12||34>
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(int_t%oovv,C8,oovv_dims,init_val=ZERO)
!use cholesky vector
         ovov_dims(1) = nocc
         ovov_dims(2) = nvir
         ovov_dims(3) = nocc
         ovov_dims(4) = nvir
         ierr=talsh_tensor_construct(test_tensor,C8,ovov_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uvo_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,int_t%oovv,12)
         ierr=talsh_tensor_destruct(test_tensor)

!        Get antisymmetric oooo tensor
         oooo_dims = nocc
         ierr=talsh_tensor_construct(int_t%oooo,C8,oooo_dims,init_val=ZERO)
!use cholesky vector
         ierr=talsh_tensor_construct(test_tensor,C8,oooo_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uoo_tensor,uoo_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,int_t%oooo,12)
         ierr=talsh_tensor_destruct(test_tensor)

!        Get anti-symmetrized vovo tensor, this needs to be done in two steps as two different classes contribute.
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=talsh_tensor_construct(int_t%vovo,C8,vovo_dims,init_val=ZERO)
         ovvo_dims(1) = nocc
         ovvo_dims(2) = nvir
         ovvo_dims(3) = nvir
         ovvo_dims(4) = nocc
         ierr=talsh_tensor_construct(ovvo_tensor,C8,ovvo_dims,init_val=ZERO)
!use cholesky vector
         ierr=talsh_tensor_construct(test_tensor,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uvv_tensor,uoo_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,int_t%vovo,0)
         ierr=talsh_tensor_destruct(test_tensor)
         ierr=talsh_tensor_construct(test_tensor,C8,(/nocc,nvir,nvir,nocc/),init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uvo_tensor,uvo_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,ovvo_tensor,0)
         ierr=talsh_tensor_destruct(test_tensor)
         ierr=talsh_tensor_contract("V(a,i,b,j)+=V(i,a,b,j)*C()",int_t%vovo,ovvo_tensor,one_tensor,scale=MINUS_ONE)
         ierr=talsh_tensor_destruct(ovvo_tensor)

!        Get vvvv tensor
         vvvv_dims = nvir
         ierr=talsh_tensor_construct(int_t%vvvv,C8,vvvv_dims,init_val=ZERO)
!use cholesky vector
         ierr=talsh_tensor_construct(test_tensor,C8,vvvv_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uvv_tensor,uvv_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,int_t%vvvv,12)
         ierr=talsh_tensor_destruct(test_tensor)

!        Get ooov tensor
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=talsh_tensor_construct(int_t%ooov,C8,ooov_dims,init_val=ZERO)
!use cholesky vector
         ierr=talsh_tensor_construct(test_tensor,C8,ooov_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uoo_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,int_t%ooov,12)
         ierr=talsh_tensor_destruct(test_tensor)

!        Get vovv tensor
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=talsh_tensor_construct(int_t%vovv,C8,vovv_dims,init_val=ZERO)
!use cholesky vector
         vvov_dims(1) = nvir
         vvov_dims(2) = nvir
         vvov_dims(3) = nocc
         vvov_dims(4) = nvir
         ierr=talsh_tensor_construct(test_tensor,C8,vvov_dims,init_val=ZERO)
         ierr=talsh_tensor_contract("T(p,q,r,s)+=K+(q,p,u)*L(r,s,u)",test_tensor,uvv_tensor,uov_tensor)
         if (ierr.ne.0) stop 'program error: cholesky not working'
         call mulliken_to_dirac_sort(test_tensor,int_t%vovv,34)
         ierr=talsh_tensor_destruct(test_tensor)

!        Clean-up
         ierr=talsh_tensor_destruct(aoint_tensor)
         ierr=talsh_tensor_destruct(chol_tensor)
         ierr=talsh_tensor_destruct(uoo_tensor)
         ierr=talsh_tensor_destruct(uov_tensor)
         ierr=talsh_tensor_destruct(uvo_tensor)
         ierr=talsh_tensor_destruct(uvv_tensor)

       end subroutine get_CC_integrals_chol

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
             print*, "WARNING: something is wrong in DIIS: found zero error for non-converged amplitudes"
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

       subroutine t3_energy (t_energy,eps_occ,eps_vir,t1_tensor,t2_tensor,fov_tensor,ooov_tensor,vovv_tensor,oovv_tensor)

!        Evaluate +T, -T and (T) energy terms.
!        Definitions of terms: L. Visscher, T.J. Lee, K.G. Dyall, J Chem Phys. 105 (1996) 8769.

         type(talsh_tens_t), intent(inout)  :: t1_tensor, t2_tensor
         type(talsh_tens_t), intent(inout)  :: fov_tensor, ooov_tensor, vovv_tensor, oovv_tensor
         real(8), intent(in)                :: eps_occ(:),eps_vir(:)
         real(8), intent(out)               :: t_energy(3)

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

                  ! Copy into w array to symmetrize and add the minus sign omitted in the above contractions
                  ierr=talsh_tensor_init(w_tensor)
                  ! symmetrization via contraction as addition with permutation is not implemented
                  ierr=talsh_tensor_contract("W(a,b,c)+=W(a,b,c)",w_tensor,we_tensor,one_tensor,scale=MINUS_ONE)
                  ierr=talsh_tensor_contract("W(a,b,c)+=W(b,c,a)",w_tensor,we_tensor,one_tensor,scale=MINUS_ONE)
                  ierr=talsh_tensor_contract("W(a,b,c)+=W(c,a,b)",w_tensor,we_tensor,one_tensor,scale=MINUS_ONE)

                  ! Form Y-intermediate and calculate 5th order -T energy contribution
                  call form_y3_intermediate (i,j,k,t1_tensor,t2_tensor,y_tensor,y)
                  ierr=talsh_tensor_init(result_tensor)
                  ierr=talsh_tensor_contract("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(3)=t_energy(3) + result_tens(1)

                  ! Form Y-intermediate and calculate 5th order (T) energy contribution
                  call form_y2_intermediate (i,j,k,t1_tensor,t2_tensor,fov_tensor,oovv_tensor,y_tensor,y)
                  eps_ijk = eps_occ(i) + eps_occ(j) + eps_occ(k)
                  call scale_with_denominators (eps_occ,eps_vir,t3_tensor=y_tensor,eps_ijk=eps_ijk)
                  ierr=talsh_tensor_init(result_tensor)
                  ierr=talsh_tensor_contract("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(2)=t_energy(2) + result_tens(1)

                  ! Copy w into y array and calculate 4th order +T energy contribution
                  ierr=talsh_tensor_init(y_tensor)
                  ierr=talsh_tensor_add("W(a,b,c)+=W(a,b,c)",y_tensor,w_tensor)
                  call scale_with_denominators (eps_occ,eps_vir,t3_tensor=y_tensor,eps_ijk=eps_ijk)
                  ierr=talsh_tensor_init(result_tensor)
                  ierr=talsh_tensor_contract("E()+=W+(a,b,c)*Y(a,b,c)",result_tensor,w_tensor,y_tensor) ! o(3)v(3)
                  t_energy(1)=t_energy(1) + result_tens(1)
               end do
            end do
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

    subroutine decompose_cholesky_talsh (CD, t_cholesky, nao)

!       Cholesky decomposition based on treating the shells seperately using talsh tensor operations
!       possible additions: 
!                          SVD (How to combine with tensor library)
!                          Additional screening in the computation of integrals

!       written by Johann Pototschnig, Summer 2019

        use exacorr_cholesky

        implicit none

        real(8), intent(in)                  :: t_cholesky ! cholesky threshold
        type(talsh_tens_t), intent(inout)    :: CD
        integer, intent(in)                  :: nao

        integer            :: ierr

        integer            :: lshell, kshell, lsize, ksize
        integer            :: l, k, l0, k0
        complex(8)         :: D_max, D_el, denom_fac
        type(talsh_tens_t) :: Mask, Diag, Res, Res2
        type(talsh_tens_t) :: Dlk, CD_v, CD_m, CD_int
        integer            :: dims2(1:2)
        integer            :: dims3(1:3), off3(1:3)
        integer            :: dims4(1:4), off4(1:4)

        integer            :: rank
        integer(4)         :: m

        ! accuracy for tensor printing
        real(8), parameter :: ACC=1.0E-9

        print *, "========= entered decompose_cholesky_talsh ============== "

        dims2=nao
        dims3=nao
        dims4=nao
        off4=0
        off3=0

        !Interest needs to be initalized for each thread, this function will initialize (or return immediately if it is).
        call interest_initialize(.false.)

        ierr=talsh_tensor_construct(Mask,R4,dims2,init_val=ONE)
        ierr=talsh_tensor_construct(Diag,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Res,C8,dims2,init_val=ZERO)
        ierr=talsh_tensor_construct(Res2,C8,dims2,init_val=ZERO)


        call talsh_get_diag_comp(Diag)
        print *, "========= diagonal computed ============== "

        call talsh_find_shells_comp(Diag, Mask, lshell, kshell, D_max)
        lsize=get_shellsize(lshell)
        ksize=get_shellsize(kshell)
        print *, "--- shell found: ", lshell, kshell, lsize, ksize, D_max

        m=0
        dims3(3)=1
        ierr=talsh_tensor_construct(CD,C8,dims3,init_val=ZERO)
        ierr=talsh_tensor_construct(CD_int,C8,dims3,init_val=ZERO)

        do while (abs(D_max)>t_cholesky)
          dims4(3)=lsize
          dims4(4)=ksize

          ierr=talsh_tensor_construct(Dlk,C8,dims4,init_val=ZERO)

          ierr=talsh_tensor_destruct(CD_int)
          dims3(3)=m+lsize*ksize
          ierr=talsh_tensor_construct(CD_int,C8,dims3,init_val=ZERO)
          if (m>0) then
            ierr = talsh_tensor_dimensions(CD,rank,dims3)
            off3=0
            ierr=talsh_tensor_insert(CD_int,CD,off3)
          end if

          !get the values for the current shell
          call talsh_get_shell_comp(Dlk, lshell, kshell)
          !call print_tensor(Dlk, ACC, "Dlk")

          call talsh_find_el_comp(Diag, Mask, lshell, kshell, l, k, D_el)
          !print *, '--- largest value found: ', l, k, D_el

          do while(abs(D_el)>abs(D_max/1000))
            call get_offset(l,l0)
            call get_offset(k,k0)

            !sum up previous CD vectors for the indices
            ierr=talsh_tensor_init(Res, ZERO)
            if (m>0) then
              dims3(1:2)=1
              dims3(3)=m
              ierr=talsh_tensor_construct(CD_v,C8,dims3,init_val=ZERO)
              dims3(1:2)=nao

              off3(1)=l-1
              off3(2)=k-1
              off3(3)=0
              ierr=talsh_tensor_slice(CD_v,CD_int,off3,0,DEV_HOST)
              ierr=talsh_tensor_reshape(CD_v,(/ m /))

              dims3(3)=m
              ierr=talsh_tensor_construct(CD_m,C8,dims3,init_val=ZERO)
              off3=0
              ierr=talsh_tensor_slice(CD_m,CD_int,off3,0,DEV_HOST)
              !call print_tensor(CD_m, ACC, "CD_m")

              ierr=talsh_tensor_contract("T(p,q)+=C(p,q,u)*D(u)",Res,CD_m,CD_v)

              ierr=talsh_tensor_destruct(CD_m)
              ierr=talsh_tensor_destruct(CD_v)
            end if

            !get Dlk for the current l and k
            ierr=talsh_tensor_init(Res2, ZERO)
            dims4(3)=1
            dims4(4)=1
            ierr=talsh_tensor_reshape(Res2,dims4)
            off4=0
            off4(3)=l-l0-1
            off4(4)=k-k0-1
            ierr=talsh_tensor_slice(Res2,Dlk,off4,0,DEV_HOST)
            ierr=talsh_tensor_reshape(Res2,dims2)

            !sum and denominate
            ierr=talsh_tensor_contract("T(p,q)+=V(p,q)",Res2,Res,one_tensor,scale=MINUS_ONE)
            ierr=talsh_tensor_init(Res, ZERO)
            denom_fac=ONE / zsqrt(D_el)
            ierr=talsh_tensor_contract("W(p,q)=V(p,q)",Res,Res2,one_tensor,scale=denom_fac)
            !call print_tensor(Res, ACC, "CD entry")

            !add Cholesky vector
            dims3(3)=1
            ierr=talsh_tensor_reshape(Res,dims3)
            off3(3)=m
            ierr=talsh_tensor_insert(CD_int,Res,off3)
            ierr=talsh_tensor_reshape(Res,dims2)

            !update diagonal
            call talsh_square_el(Res)
            ierr=talsh_tensor_contract("W(p,q)+=V(p,q)",Diag,Res,one_tensor,scale=MINUS_ONE)
            !call print_tensor(Diag, ACC, "Diag")

            !update mask
            call talsh_set_mask_done(Mask,l,k)

            m=m+1

            !update CD
            ierr=talsh_tensor_destruct(CD)
            dims3(3)=m
            ierr=talsh_tensor_construct(CD,C8,dims3,init_val=ZERO)
            off3=0
            ierr=talsh_tensor_slice(CD,CD_int,off3,0,DEV_HOST)
            !call print_tensor(CD, ACC, "CD")

            if (m>nao*nao) print *, 'Cholesky decomposition did not work'

            !find next element in this shell
            call talsh_find_el_comp(Diag, Mask, lshell, kshell, l, k, D_el)
            !print *, '--- largest value found: ', l, k, D_el

          end do

          !find shell with next largest element
          call talsh_find_shells_comp(Diag, Mask, lshell, kshell, D_max)
          lsize=get_shellsize(lshell)
          ksize=get_shellsize(kshell)
          write (*,'(a,i6,a,a,i6,i6,a,e12.2)') 'cholesky length =',m,' largest element:', &
          ' shell = ', lshell, kshell, ' value = ', abs(D_max)

          ierr=talsh_tensor_destruct(Dlk)

        end do

        print *, 'Cholesky decomposition done, used ',m,' elements'

        !clean up
        ierr=talsh_tensor_destruct(CD_int)

        ierr=talsh_tensor_destruct(Mask)
        ierr=talsh_tensor_destruct(Diag)
        ierr=talsh_tensor_destruct(Res)
        ierr=talsh_tensor_destruct(Res2)

        print *, "========= leaving decompose_cholesky_talsh ============== "

    end subroutine decompose_cholesky_talsh

end module talsh_cc
