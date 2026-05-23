module eti_cc
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_cfg, only: exacc_input
        use exacorr_utils
        use exacorr_printout
        use exacorr_scaling_parameters
        use intermediates_cc
        use intermediates_lambda
        use intermediates_eom_ccsd
        use eti_common_routines
        use eti_drivers_gs_ccsd
        use eti_drivers_gs_automatic
        use eti_diis
        use mfsolver_dirac_input
        use, intrinsic:: ISO_C_BINDING

        implicit none

        type(tensor_t) :: one_tensor

        private

        public eti_cc_driver
        public eti_init_fock

!        interface scale_with_denominators
!           module procedure scale_with_denominators
!        end interface
!

       contains
!-------------------------------------------------------------------------------------------
        subroutine eti_cc_driver (exa_input)

!This subroutine drives the ccsd energy calculation
!Currently everything is done in-core, purpose is to provide a working implementation that is suitable for debugging the production code

!        Written by Lucas Visscher, winter 2016/2017
         use os_utils
         use labeled_storage
         use exacorr_global
         use eti_gradient
         use eti_eom
         use eti_ccrsp,    only       : eti_CCRSP_driver
         use eti_QR_TPA,   only       : eti_QR_TPA_driver
         use eti_cholesky_integrals
         use eti_integrals
         use eti_triples

         type(exacc_input), intent(in) :: exa_input
         type(file_info_t)      :: exacc_file
         integer                :: nocc,nvir   ! the size of the mo basis for occupied and virtual spinors
         integer,  allocatable  :: mo_occ(:),mo_vir(:) ! the list of occupied and virtual orbitals

         integer(C_INT)       :: ierr
         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p

!        Common tensors
         type(eti_comm_tens_t) :: comm_t
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_intg_tens_t) :: int_t
!        Lambda tensors
!        1- and 2-body cc intermediates 
         type(eti_intermediates_tens_t) :: interm_t
!        eom results
         type(eti_eom_results_tens_t) :: eom_results_t
!        intermediate tensors (same notation as in relccsd)
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
         integer(INTD), dimension(4) :: oooo_dims, oovv_dims, vvoo_dims, vovo_dims
         integer(INTD), dimension(6) :: vvvooo_dims, ooovvv_dims
         integer(INTD), dimension(8) :: vvvvoooo_dims, oooovvvv_dims
         integer(INTD), dimension(10) :: vvvvvooooo_dims, ooooovvvvv_dims
!        scalars (need to be defined as tensor types)
         type(tensor_t) :: result_tensor
         integer(C_INT)     :: one_dims(1)
!        orbital energies
         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CCSD,CC2,... control and result variables
         real(8) :: e_core, scf_energy, mp2_energy, cc_energy, pert_energy(3)
         real(8) :: t_target_precision
         real(8) :: t1diag, t2diag
         integer :: ncycles
         logical :: CCD, lambda, eom, eom_lr, eom_qr, qr_tpa
         logical :: debug_eom = .false.
         logical :: CC2
         logical :: CCSDT, CCSDTQ, CCSDTQP
         logical :: do_perturbative
         integer :: iterative_level, iterative_scheme
         integer :: perturbative_level, perturbative_scheme
         logical :: do1rdm_gs, do2rdm_gs

         integer :: iff        ! loop variable to loop over field strengths
         integer :: fock_mode  ! choice of Fock build options
         integer :: rr         ! rank of cc-ci response

! this should be needed in the case of exatensor, but not for talsh i guess..
!        call mfsolver_dirac_input_sync_cw()

!        Make some noise so that we know we are here
         call print_date('Entered cc_driver routine')

!        The user may have modified the default executor via input
         if (exa_input%executor /= 0) then ! 0 means no input value is set
            ierr = set_default_executor(exa_input%executor)
            if (ierr /= 0) write (*,*) ' WARNING: incompatible executor specified in input, keeping default'
            select case(get_default_executor())
            case(1) 
                write(*,*) ' User set TAPP (executor 1) via input as the default for tensor operations'
            case(2)
                write(*,*) ' User set TALSH (executor 2) via input as the default for tensor operations'
            case(3)
                write(*,*) ' User set ExaTensor (executor 3) via input as the default for tensor operations'
            end select
         end if

!        The experimental Cholesky code works only with TALSH as it needs slicing and reshaping of tensors
         if (exa_input%moint_scheme > 40 .and. exa_input%executor /= 2) then
             write(*,*) ' Cholesky feature requested, thus therefore overruling default and forcing use of executor 2 (TALSH)'
             ierr = set_default_executor(2)
         end if
         
         call eti_executor_initialize(exa_input%talsh_buff,get_default_executor())

         call interest_initialize(.true.)
         call print_date('Initialized interest library')

!        Initialize scalars that are to be used as tensors in contractions
         one_dims(1) = 1
         ierr=construct_tensor_t(one_tensor,C8,one_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(result_tensor,C8,one_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,one_dims)

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nocc
         nvir = exa_input%nvir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ = exa_input%mo_occ
         mo_vir = exa_input%mo_vir
         t_target_precision = exa_input%t_econv

         do_perturbative = exa_input%do_perturbative
         CCD = exa_input%ccd
         CC2 = exa_input%cc2

         iterative_level     = exa_input%iterative_level
         iterative_scheme    = exa_input%iterative_scheme
         perturbative_level  = exa_input%perturbative_level
         perturbative_scheme = exa_input%perturbative_scheme

         CCSDT = .false.
         CCSDTQ = .false.
         CCSDTQP = .false.
         if (iterative_level.eq.3) CCSDT =.true.
         if (iterative_level.eq.4) CCSDTQ =.true.
         if (iterative_level.eq.5) CCSDTQP =.true.

         if (CCD.or.CC2.or.CCSDT.or.CCSDTQ.or.CCSDTQP) then
            do_perturbative = .false.
         end if

         lambda = exa_input%lambda
         eom = exa_input%do_eomcc
         eom_lr = exa_input%do_lr_cc
         eom_qr = exa_input%do_qr_cc
         qr_tpa = exa_input%qr_tpa
         ncycles = exa_input%ncycles
         do1rdm_gs = exa_input%do_1rdm_gs
         do2rdm_gs = exa_input%do_2rdm_gs
         if ((do1rdm_gs.or.do2rdm_gs).and..not.lambda) then
            call quit("specify .LAMBDA to generate 1RDM/2RDM")
         end if

         call print_date("MO occupation data to EXACC file")
         exacc_file%status = -1
         exacc_file%type = 2
         exacc_file%name = "EXACC.h5"
         call lab_write(exacc_file,'nocc',idata=nocc)
         call lab_write(exacc_file,'mo_occ',idata=mo_occ)

#if 0
         call print_exainput(exa_input,get_nao())
#endif

         call print_date("Getting orbital energies")
!        Get orbital energies from the DIRAC restart file
         allocate (eps_occ(nocc))
         allocate (eps_vir(nvir))
         call get_orbital_energies(nocc,mo_occ,eps_occ)
         call get_orbital_energies(nvir,mo_vir,eps_vir)
         call print_date("Retrieved orbital energies")

         ! print *,'before get_CC_integrals'
         ! print *,' eps_occ:',eps_occ
         ! print *,' eps_vir:',eps_vir
!        Get fixed two-body tensors
         call print_date('-Start- Integral Transformation')
         if (exa_input%moint_scheme.eq.42) then
           call get_CC_integrals_chol (nocc,nvir,mo_occ,mo_vir,int_t,exa_input%t_cholesky,exa_input%print_level)
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
         call eti_init_fock (fock_mode,mo_occ,mo_vir,eps_occ,eps_vir,int_t,scf_energy,exa_input%print_level)

!        create tensors needed in CC
         oo_dims = nocc
         vv_dims = nvir
         ierr=construct_tensor_t(comm_t%foo,C8,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%fvv,C8,vv_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%goo,C8,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%hoo,C8,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%gvv,C8,vv_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%hvv,C8,vv_dims,init_val=ZERO)
         if (.not.CCD) then
            ov_dims(1) = nocc
            ov_dims(2) = nvir
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=construct_tensor_t(comm_t%fov,C8,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(comm_t%hov,C8,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(comm_t%t1,C8,vo_dims,init_val=ZERO)
!           initialize t3 amplitudes to zero (will later update them to starting values based on mp2 amplitudes)
            if (iterative_level.ge.3) then
               vvvooo_dims(1:3) = nvir
               vvvooo_dims(4:6) = nocc
               ierr=construct_tensor_t(comm_t%t3,C8,vvvooo_dims,init_val=ZERO)
            end if
            if (iterative_level.ge.4) then
               vvvvoooo_dims(1:4) = nvir
               vvvvoooo_dims(5:8) = nocc
               ierr=construct_tensor_t(comm_t%t4,C8,vvvvoooo_dims,init_val=ZERO)
            end if
            if (iterative_level.ge.5) then
               vvvvvooooo_dims(1:5) = nvir
               vvvvvooooo_dims(6:10) = nocc
               ierr=construct_tensor_t(comm_t%t5,C8,vvvvvooooo_dims,init_val=ZERO)
            end if
         end if
         oooo_dims = nocc
         vovo_dims = nvir
         vovo_dims(2) = nocc
         vovo_dims(4) = nocc
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=construct_tensor_t(comm_t%a_int,C8,oooo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%h_int,C8,vovo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%tau,C8,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%t2,C8,vvoo_dims,init_val=ZERO)
         if (ierr.eq.ETI_SUCCESS) then
            call print_date('Succesfully allocated all tensors')
         else
            write(*,*) " ccdriver: error in tensor allocation, code=",ierr
            return
         end if

         ! copy Fock matrix to comm, add finite field to the Fock matrix
         call set_fock(exa_input,mo_occ,mo_vir,int_t,comm_t,-1)
         ! define a diagonal zeroth order Hamiltonian with H^1 = H - H^0 and H^0_pq = epsilon_p delta_pq
         call define_H0(scf_energy,eps_occ,eps_vir,comm_t,int_t%oooo, &
                        exa_input%level_shift,-1,fock_mode,exa_input%nff(1),CCD)

!        Initialize t1 amplitudes to fov
         if (.not.CCD) then 
            ierr=initialize_tensor_t(comm_t%t1)
            ierr=t_einsum("T(a,i)+=V+(i,a)",comm_t%t1,comm_t%fov,one_tensor)
         end if

!        Initialize t2 amplitudes with the MP2 values
         ierr=initialize_tensor_t(comm_t%t2)
         ierr=t_einsum("T(a,b,i,j)+=V+(i,j,a,b)",comm_t%t2,int_t%oovv,one_tensor)

         if (.not.CCD) then
            if (iterative_level.eq.5) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2, &
                              full_t3_tensor=comm_t%t3, full_t4_tensor=comm_t%t4, &
                              full_t5_tensor=comm_t%t5)
            else if (iterative_level.eq.4) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2, &
                              full_t3_tensor=comm_t%t3, full_t4_tensor=comm_t%t4)
            else if (iterative_level.eq.3) then
               call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2,full_t3_tensor=comm_t%t3)
            else ! CCSD and variants
               call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)
            end if
         else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
         end if

!        Calculate and print the MP2 energy
         ierr=initialize_tensor_t(result_tensor)
         ierr=t_einsum("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%t2,int_t%oovv) ! o(2)v(2)
         mp2_energy=result_tens(1)*ONE_QUARTER_R

         if (.not.CCD) then
           ierr=initialize_tensor_t(result_tensor)
           ierr=t_einsum("E()+=T(a,i)*V(i,a)",result_tensor,comm_t%t1,comm_t%fov)
           mp2_energy=mp2_energy+result_tens(1)
         end if

         write(*,*) ""
         write(*,*) "    MP2 energy = ", mp2_energy
         write(*,*) ""
         ! call print_tensor(comm_t%t2, 1.0d-18, "t2 mp2")
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
                         exa_input%level_shift,exa_input%print_level,fock_mode,exa_input%nff(1),CCD)

!        Start CCSD iterations
         call print_date('-Start- CC Iterations')

         call print_iterative_wavefunction_choice(iterative_scheme)

         if (.not.CCD) then
            if (iterative_level.eq.5) then
               call diis (t5_tensor=comm_t%t5,t4_tensor=comm_t%t4,t3_tensor=comm_t%t3, &
                          t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,mode=1)
            else if (iterative_level.eq.4) then
               call diis (t4_tensor=comm_t%t4,t2_tensor=comm_t%t2,t3_tensor=comm_t%t3,t1_tensor=comm_t%t1,mode=1)
            else if (iterative_level.eq.3) then
               call diis (t2_tensor=comm_t%t2,t3_tensor=comm_t%t3,t1_tensor=comm_t%t1,mode=1)
            else
               call diis (t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,mode=1)
            end if
         else
            call diis (t2_tensor=comm_t%t2,mode=1)
         end if

         if (exa_input%auto_kernels) then
            call solve_t_equations_automatic(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision,&
                                 cc_energy,CCD,iterative_level,iterative_scheme,exa_input%auto_scheme, &
                                 ncycles,exa_input%print_level)
         else if (CC2) then
            call solve_CC2_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                 cc_energy,CCD,ncycles,exa_input%print_level)
         else
            call solve_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                 cc_energy,CCD,iterative_level,iterative_scheme,ncycles,exa_input%print_level)
         end if
         call print_date('-End- CC Iterations')

         if (.not.CCD) then
            ! calculate T1 diagnostic
            t1diag = tensor_norm2(comm_t%t1)
            t1diag = t1diag/sqrt(dble(nocc))
         else
            t1diag = 0.0d0
         end if
            ! calculate T2 diagnostic
            t2diag = tensor_norm2(comm_t%t2)
            t2diag = t2diag/sqrt(dble(nocc))

         if (do_perturbative) then
           if (perturbative_scheme.eq.30) then
              !Calculate (T) energy correction
              call print_date('-Start- Triples')
              if (exa_input%tripl_lap) then
                call t3_energy_laplace (pert_energy,comm_t,int_t,eps_occ,eps_vir,exa_input)
              else
                call t3_energy (pert_energy,eps_occ,eps_vir, comm_t%t1,comm_t%t2,comm_t%fov, &
                           int_t%ooov,int_t%vovv,int_t%oovv,exa_input%print_level)
              end if
              call print_date('-End- Triples')
           end if
           call print_exaoutput(exa_input%talsh, iterative_scheme, perturbative_scheme, &
                                scf_energy, mp2_energy, cc_energy, t1diag, pert_energy)
         else
            call print_exaoutput(exa_input%talsh, iterative_scheme, perturbative_scheme,&
                                 scf_energy, mp2_energy, cc_energy, t1diag)
         end if

         !Destruct diis tensors
         if (.not.CCD) then
            if (iterative_level.eq.5) then
               call diis (t5_tensor=comm_t%t5,t4_tensor=comm_t%t4,t3_tensor=comm_t%t3,&
                          t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,mode=2)
            else if (iterative_level.eq.4) then
               call diis (t4_tensor=comm_t%t4,t3_tensor=comm_t%t3,t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,mode=2)
            else if (iterative_level.eq.3) then
               call diis (t3_tensor=comm_t%t3,t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,mode=2)
            else
               call diis (t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,mode=2)
            end if
         else
            call diis (t2_tensor=comm_t%t2,mode=2)
         end if

         ! end finite field loop
         end do

        if (exa_input%auto_kernels.and.lambda.and..not.CC2) then
           call print_date('-Start- Lambda Equations with generated code')
!          Construct and initialize lambda tensor
           if (.not.CCD) ierr=construct_tensor_t(comm_t%l1,C8,ov_dims,init_val=ZERO)
           oovv_dims(1:2) = nocc
           oovv_dims(3:4) = nvir
           ierr=construct_tensor_t(comm_t%l2,C8,oovv_dims,init_val=ZERO)
           if (.not.CCD) ierr=t_einsum("L(i,a)+=T+(a,i)*K()",comm_t%l1,comm_t%t1,one_tensor)
           ierr=t_einsum("L(i,j,a,b)+=T+(a,b,i,j)*K()",comm_t%l2,comm_t%t2,one_tensor)
           if (iterative_level.ge.3) then
              ooovvv_dims(1:3) = nocc
              ooovvv_dims(4:6) =nvir
              ierr=construct_tensor_t(comm_t%l3,C8,ooovvv_dims,init_val=ZERO)
              ierr=t_einsum("L(i,j,k,a,b,c)+=T+(a,b,c,i,j,k)*K()",comm_t%l3,&
                                                                comm_t%t3,one_tensor)
           end if 
           if (iterative_level.ge.4) then
              oooovvvv_dims(1:4) = nocc
              oooovvvv_dims(5:8) =nvir
              ierr=construct_tensor_t(comm_t%l4,C8,oooovvvv_dims,init_val=ZERO)
              ierr=t_einsum("L(i,j,k,l,a,b,c,d)+=T+(a,b,c,d,i,j,k,l)*K()",comm_t%l4,&
                                       comm_t%t4,one_tensor)
           end if
!          Solve the lambda-equations and calculate density matrix in MO basis
           call solve_lambda_equations_automatic(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,&
                            t_target_precision,CCD,iterative_level,iterative_scheme,exa_input%auto_scheme,&
                            ncycles,exa_input%print_level)
           call eti_dm_driver(exa_input, int_t, comm_t)
           call print_date('-End- Lambda Equations')
         
         else if (lambda.and..not.CC2) then
          call print_date('-Start- Lambda Equations')
!         Construct lambda tensor
          if (.not.CCD) ierr=construct_tensor_t(comm_t%l1,C8,ov_dims,init_val=ZERO)
          oovv_dims(1:2) = nocc
          oovv_dims(3:4) = nvir
          ierr=construct_tensor_t(comm_t%l2,C8,oovv_dims,init_val=ZERO)

!         Initialize l1 and l2 with converged t amplitudes
          if (.not.CCD) ierr=t_einsum("L(i,a)+=T+(a,i)*K()",comm_t%l1,comm_t%t1,one_tensor)
          ierr=t_einsum("L(i,j,a,b)+=T+(a,b,i,j)*K()",comm_t%l2,comm_t%t2,one_tensor)

!         solve lambda-equations and calculate density matrix in MO basis
          if (.not.CCD) then
            if (.not.CC2) then
                call solve_lambda_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir, &
                                    CCD,t_target_precision,ncycles,exa_input%print_level)
             else
                call solve_CC2_lambda_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir, &
                                    t_target_precision,ncycles,exa_input%print_level)
            end if 
            call eti_dm_driver(exa_input, int_t, comm_t)
          else
            call solve_lambda_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir, &
                                        CCD,t_target_precision,ncycles,exa_input%print_level)
            call eti_dm_driver(exa_input, int_t, comm_t)
          end if
          call print_date('-End- Lambda Equations')
         end if
!        Destruct tensors no longer needed for ground-state
!        note the amplitude tensor may still be reused in eom or response, so they need to live longer
         ierr=destroy_tensor_t(comm_t%tau)
         ierr=destroy_tensor_t(comm_t%hoo)
         if (.not.CCD) then
            ierr=destroy_tensor_t(comm_t%hov)
         end if
         ierr=destroy_tensor_t(comm_t%hvv)
         ierr=destroy_tensor_t(comm_t%goo)
         ierr=destroy_tensor_t(comm_t%gvv)
         ierr=destroy_tensor_t(comm_t%a_int)
         ierr=destroy_tensor_t(comm_t%h_int)
         ierr=destroy_tensor_t(result_tensor)
         ierr=destroy_tensor_t(one_tensor)
         
         ! Jan fixed for ETI
         ierr=destroy_tensor_t(int_t%foo)
         ierr=destroy_tensor_t(int_t%fvv)
         ierr=destroy_tensor_t(int_t%fov)

!**************************
!        EOM calculations *
!**************************

         if (eom.or.eom_lr.or.eom_qr.or.qr_tpa) then

!aspg: if requesting properties, check if lagrange multipliers are present. if these are not there
!      in EE-based calculations (IP and EA properties not yet implemented), stop 
            if ((exa_input%eom_input(1)%do_expval.or.eom_lr.or.eom_qr.or.qr_tpa).and..not.lambda) then
               call quit("To calculate properties, Lambda amplitudes are required ! Please check input.")
            end if

!aspg, this new handling of fock matrices has to be verified, we don't want any addition to finite field ere for example
!           Get fixed one-body tensors, needs to follow the two body tensors
            if ( exist_one_el(e_core) ) then
               fock_mode = 2
            else
               fock_mode = 1
            end if
            call eti_init_fock (fock_mode,mo_occ,mo_vir,eps_occ,eps_vir,int_t,scf_energy,exa_input%print_level)
            ! copy Fock matrix to comm, add finite field to the Fock matrix
            call set_fock(exa_input,mo_occ,mo_vir,int_t,comm_t,-1)
            ! define a diagonal zeroth order Hamiltonian with H^1 = H - H^0 and H^0_pq = epsilon_p delta_pq
            call define_H0(scf_energy,eps_occ,eps_vir,comm_t,int_t%oooo, &
                        exa_input%level_shift,-1,fock_mode,exa_input%nff(1),CCD,keep_diagonal=.true.)

! create eom intermediates
            call form_ccsd_intermediates(interm_t, exa_input%CCD, exa_input%nocc, exa_input%nvir, comm_t, int_t)

! solve eom eigenvalue equations; right now we do only one flavor but there's no reason not to handle different flavors in the same calculation,
! we just need different results instances
            if (eom_qr) then
               rr = 2
               call eti_CCRSP_driver(exa_input,interm_t,int_t,rr)
            else if (eom_lr) then
               rr = 1
               call eti_CCRSP_driver(exa_input,interm_t,int_t,rr)
! TODO: if we want properties etc, now we call other drivers
            else if (qr_tpa) then
               call eti_eom_results_create(eom_results_t,exa_input%eom_input(1)%nroots,.true.,.true.)
               call eti_eom_eigval_driver(eom_results_t,exa_input%eom_input(1),interm_t,int_t,&
                                          solve_right=.true.,solve_left=.false.)
 
               call eti_eom_eigval_driver(eom_results_t,exa_input%eom_input(1),interm_t,int_t,&
                                          solve_right=.false.,solve_left=.true.)
               call eti_QR_TPA_driver(exa_input,interm_t,int_t,eom_results_t)
               call eti_eom_results_destroy(eom_results_t)
            else
               call eti_eom_results_create(eom_results_t,exa_input%eom_input(1)%nroots,right=.true.,left=.false.)
               call eti_eom_eigval_driver(eom_results_t,exa_input%eom_input(1),interm_t,int_t,&
                             solve_right=.true.,solve_left=.false.)
               rr = 2
               ! call eti_CCRSP_driver(exa_input,interm_t,int_t,l1_tensor,l2_tensor, &
               !          rr,eom_results_t,eom_results_t_left)
               call eti_eom_results_destroy(eom_results_t)
            end if
         end if

         if (.not.CCD) then
           ierr=destroy_tensor_t(comm_t%t1)
           ierr=destroy_tensor_t(comm_t%l1)
           if (iterative_level.ge.3) ierr=destroy_tensor_t(comm_t%t3)
           if (iterative_level.ge.4) ierr=destroy_tensor_t(comm_t%t4)
           if (iterative_level.ge.5) ierr=destroy_tensor_t(comm_t%t5)
         end if
         ierr=destroy_tensor_t(comm_t%l2)
! aspg
         ierr=destroy_tensor_t(comm_t%foo)
         ierr=destroy_tensor_t(comm_t%fov)
         ierr=destroy_tensor_t(comm_t%fvv)
! aspg
         ierr=destroy_tensor_t(comm_t%t2)

!        shutdown tensor libraries
         call eti_executor_shutdown

         deallocate (eps_occ)
         deallocate (eps_vir)
         deallocate (mo_occ)
         deallocate (mo_vir)

!        Free global data used to interact with DIRAC/Interest
         call delete_global_data

!        Make some noise so that we know we are leaving
         write(*,*) ""
         call print_date('Leaving cc_driver routine')

        end subroutine eti_cc_driver


        subroutine eti_init_fock (mode,mo_occ,mo_vir,eps_occ,eps_vir,int_t,scf_energy,print_level)
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
         type(eti_intg_tens_t), intent(inout) :: int_t

         integer                         :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         type(tensor_t)              :: id_mat_occ
         type(tensor_t)              :: h_tensor
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
         ierr=construct_tensor_t(int_t%foo,C8,fo_dims,init_val=ZERO)

!        Get fvv tensor
         fv_dims(1:2) = nvir
         ierr=construct_tensor_t(int_t%fvv,C8,fv_dims,init_val=ZERO)

!        Get fov tensor
         fov_dims(1) = nocc
         fov_dims(2) = nvir
         ierr=construct_tensor_t(int_t%fov,C8,fov_dims,init_val=ZERO)

         select case (mode)

         case (3)

!            construct Fock matrix from orbital energies and valence 2e integrals
!            to be written by Chima Chibueze

         case (2)

!           get one electron integrals
            one_el_exist = exist_one_el(e_core)

            if (.not. one_el_exist) stop "Requested Fock matrix construction, but MRCONEE is missing"

!           read the one electron integrals and copy them into the tensors

            ierr=direct_data_pointer(int_t%foo,body_p)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_occ,nocc,min_occ)

            ierr=direct_data_pointer(int_t%fvv,body_p)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_vir,nvir,mo_vir,nvir,min_occ)

            ierr=direct_data_pointer(int_t%fov,body_p)
            call c_f_pointer(body_p,h_tens,fov_dims) ! to use as a regular Fortran array
            call get_one_el(h_tens,mo_occ,nocc,mo_vir,nvir,min_occ)

            call print_date('one electron integrals added to Fock Matrix')

            !create density (identity) matrix for the occupied space
            ierr=construct_tensor_t(id_mat_occ,C8,fo_dims,init_val=ZERO)
            ierr=direct_data_pointer(id_mat_occ,body_p)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, nocc
              h_tens(i,i) = ONE
            enddo

            ierr=t_einsum("F(i,j)+=V(i,k,j,l)*I(k,l)",int_t%foo,int_t%oooo,id_mat_occ, scale=ONE)
            ierr=t_einsum("F(i,a)+=V(k,i,l,a)*I(k,l)",int_t%fov,int_t%ooov,id_mat_occ, scale=ONE)
            ierr=t_einsum("F(a,b)+=V(a,k,b,l)*I(k,l)",int_t%fvv,int_t%vovo,id_mat_occ, scale=ONE)

            call print_date('two electron integrals added to Fock Matrix')

            if (print_level >10) then
              call print_tensor(int_t%foo, 1.0D-6, 'fo_tensor')
              call print_tensor(int_t%fov, 1.0D-6, 'fov_tensor')
              call print_tensor(int_t%fvv, 1.0D-6, 'fv_tensor')
            end if

            ierr=destroy_tensor_t(id_mat_occ)

         case (1)

            ! Fock matrix is defined as diagonal matrix with orbital energies on the diagonal
            ierr=direct_data_pointer(int_t%foo,body_p)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            do i=1,fo_dims(1)
               h_tens(i,i)=eps_occ(i)
            end do

            ierr=direct_data_pointer(int_t%fvv,body_p)
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

        end subroutine eti_init_fock

        subroutine define_H0 (scf_energy,eps_occ,eps_vir,comm_t,oooo,level_shift,print_level,mode,nff,CCD,keep_diagonal)

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
         logical, intent(in)                  :: CCD


!        common tensors: the foo and fvv parts get modfied here (diagonal is shifted to the left)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(tensor_t)                   :: oooo

         integer                         :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         type(tensor_t)              :: diag_mat_occ, diag_mat_vir
         type(tensor_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)
         integer(INTD)                   :: ierr
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)
         integer                         :: i

         logical, intent(in), optional:: keep_diagonal
         logical            :: keep_diagonal_local

         if (present(keep_diagonal)) then
           keep_diagonal_local = keep_diagonal
         else
           keep_diagonal_local = .false.
         end if

         nocc = size(eps_occ)
         nvir = size(eps_vir)
         fo_dims(1:2) = nocc
         fv_dims(1:2) = nvir

         ! The epsilons are initialized from the diagonal of the UHF Fock matrix
         ierr=direct_data_pointer(comm_t%foo,body_p)
         call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
         do i=1,fo_dims(1)
            eps_occ(i)=h_tens(i,i)
         end do

         ierr=direct_data_pointer(comm_t%fvv,body_p)
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
         ierr=construct_tensor_t(diag_mat_occ,C8,fo_dims,init_val=ZERO)
         ierr=direct_data_pointer(diag_mat_occ,body_p)
         call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nocc
            h_tens(i,i) = eps_occ(i)
         end do

         ierr=construct_tensor_t(diag_mat_vir,C8,fv_dims,init_val=ZERO)
         ierr=direct_data_pointer(diag_mat_vir,body_p)
         call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
         h_tens = ZERO
         do i = 1, nvir
           h_tens(i,i) = eps_vir(i)
         end do

         if (.not.keep_diagonal_local) then
            ! create foo and fvv with H0 subtracted
            ierr=t_einsum("F(i,j)+=D(i,j)",comm_t%foo,diag_mat_occ,scale=MINUS_ONE)
            ierr=t_einsum("F(a,b)+=D(a,b)",comm_t%fvv,diag_mat_vir,scale=MINUS_ONE)
         end if
         if (print_level >0) then
           write(*,*)
           if (keep_diagonal_local) then
              write(*,*) "    Norm of Fock matrix (OO, all        ) = ", tensor_norm2(comm_t%foo)
           else
              write(*,*) "    Norm of Fock matrix (OO, no diagonal) = ", tensor_norm2(comm_t%foo)
           end if
           ! for CCD comm_t%fov is not constructed
           if (.not.CCD) write(*,*) "    Norm of Fock matrix (OV, all        ) = ", tensor_norm2(comm_t%fov)
           if (keep_diagonal_local) then
              write(*,*) "    Norm of Fock matrix (VV, all        ) = ", tensor_norm2(comm_t%fvv)
           else
              write(*,*) "    Norm of Fock matrix (VV, no diagonal) = ", tensor_norm2(comm_t%fvv)
           end if
           write(*,*)
         end if

         ierr=destroy_tensor_t(diag_mat_occ)
         ierr=destroy_tensor_t(diag_mat_vir)

        end subroutine define_H0

        subroutine compute_scf_energy(scf_energy,foo,oooo,foo_dims,print_level)

            use exacorr_global

!           Refactored by by Johann Pototschnig, spring 2021

            type(tensor_t)                 :: foo
            type(tensor_t)                 :: oooo
            real(8), intent(inout)             :: scf_energy
            integer(INTD)                      :: foo_dims(1:2)
            integer                            :: print_level

            type(tensor_t)              :: id_mat_occ, h_tensor
            logical                         :: one_el_exist
            real(8)                         :: e_core
            integer(INTD)                   :: ierr
            type(C_PTR)                     :: body_p
            complex(8), pointer, contiguous :: h_tens(:,:)
            type(tensor_t)              :: result_tensor
            integer(INTD)                   :: result_dims(1)
            complex(8), pointer             :: result_tens(:)
            integer                         :: i

!           get one electron integrals
            one_el_exist = exist_one_el(e_core)

            if (.not. one_el_exist) stop "Requested Fock matrix computation, but MRCONEE is missing"

            !Make a result tensor and arrange access to it via a Fortran pointer
            result_dims(1) = 1
            ierr=construct_tensor_t(result_tensor,C8,result_dims(1:0),init_val=ZERO)
            ierr=direct_data_pointer(result_tensor,body_p)
            call c_f_pointer(body_p,result_tens,result_dims)

            !create density (identity) matrix for the occupied space
            ierr=construct_tensor_t(id_mat_occ,C8,foo_dims,init_val=ZERO)
            ierr=direct_data_pointer(id_mat_occ,body_p)
            call c_f_pointer(body_p,h_tens,foo_dims) ! to use as a regular Fortran array
            h_tens = ZERO
            do i = 1, foo_dims(1)
              h_tens(i,i) = ONE
            enddo
            ! recompute HF energy
            ierr=initialize_tensor_t(result_tensor)
            ierr=t_einsum("R()=F(i,j)*I(i,j)",result_tensor,foo,id_mat_occ)
            ierr=construct_tensor_t(h_tensor,C8,foo_dims,init_val=ZERO)
            ierr=t_einsum("H(i,j)+=G(i,k,j,l)*I(k,l)",h_tensor,oooo,id_mat_occ, scale=MINUS_ONE_HALF)
            ierr=t_einsum("R()+=H(i,j)*I(i,j)",result_tensor,h_tensor,id_mat_occ)
            ierr=destroy_tensor_t(h_tensor)

            scf_energy=real(result_tens(1),8) ! Casting into real(8)
            scf_energy=scf_energy+e_core
            if (print_level>-1) then
               write(*,*) ""
               write (*,'(A,F23.15)') "   Recomputed energy (single determinant, calculated) = ", scf_energy
               write(*,*) ""
            end if

            ierr=destroy_tensor_t(result_tensor)
            ierr=destroy_tensor_t(id_mat_occ)

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
         type(eti_intg_tens_t), intent(inout) :: int_t
!        target tensors to be augmented with the finite field
         type(eti_comm_tens_t), intent(inout) :: comm_t

         integer                         :: nocc,nvir  ! the size of the mo basis for occupied and virtual spinors
         type(tensor_t)              :: h_tensor
         integer(INTD)                   :: fo_dims(1:2),fv_dims(1:2),fov_dims(1:2)
         integer(INTD)                   :: ierr
         integer                         :: i,min_occ
         type(C_PTR)                     :: body_p
         complex(8), pointer, contiguous :: h_tens(:,:)

!        copy the matrices
         ierr=initialize_tensor_t(comm_t%foo)
         ierr=initialize_tensor_t(comm_t%fvv)
         if (.not.exa_input%ccd) ierr=initialize_tensor_t(comm_t%fov)
         ierr=t_einsum("F(i,j)+=D(i,j)",comm_t%foo,int_t%foo)
         ierr=t_einsum("F(a,b)+=D(a,b)",comm_t%fvv,int_t%fvv)
         if (.not.exa_input%ccd) ierr=t_einsum("F(i,a)+=D(i,a)",comm_t%fov,int_t%fov)

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

            ierr=direct_data_pointer(comm_t%foo,body_p)
            call c_f_pointer(body_p,h_tens,fo_dims) ! to use as a regular Fortran array
            call add_finitefield(h_tens,exa_input,iff,mo_occ,nocc,mo_occ,nocc,min_occ)

            ierr=direct_data_pointer(comm_t%fvv,body_p)
            call c_f_pointer(body_p,h_tens,fv_dims) ! to use as a regular Fortran array
            call add_finitefield(h_tens,exa_input,iff,mo_vir,nvir,mo_vir,nvir,min_occ)

            if (.not.exa_input%ccd) ierr=direct_data_pointer(comm_t%fov,body_p)
            if (.not.exa_input%ccd) call c_f_pointer(body_p,h_tens,fov_dims) ! to use as a regular Fortran array
            if (.not.exa_input%ccd) call add_finitefield(h_tens,exa_input,iff,mo_occ,nocc,mo_vir,nvir,min_occ)

            do i=1,exa_input%nff(1)
               write(*,*) ""
               write (*,'(3A,2E12.3)') "    Added finite field ",exa_input%ff_names(i), &
                                    " with strength ",exa_input%ff(i,iff)
               write(*,*) ""
            end do

            call print_date('Finite field added to Fock Matrix')
         end if

        end subroutine set_fock

end module eti_cc
