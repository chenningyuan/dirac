program test_teqs_auto_eti

   use tensor_algebra
   use exacorr_tensor_interface
   use exacorr_datatypes
   use exacorr_printout
   use os_utils
   use exacorr_utils
   use eti_common_routines
   use eti_sigma_eomccsd
   use exacorr_scaling_parameters
   use eti_drivers_gs_automatic

   implicit none

   type(eti_intg_tens_t) :: int_t
   type(eti_intermediates_tens_t) :: interm_t
   type(eti_sa_tens_t) :: r_tensor, sigma_tensor
   type(eti_comm_tens_t) :: comm_t

   integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims, vv_dims
   integer(INTD), dimension(4) :: oooo_dims, oovv_dims, vvoo_dims, vovo_dims
   integer(INTD), dimension(6) :: vvvooo_dims, ooovvv_dims
   integer(INTD), dimension(8) :: vvvvoooo_dims, oooovvvv_dims
   integer(INTD), dimension(10) :: vvvvvooooo_dims, ooooovvvvv_dims

   integer :: nocc, nvir 
   integer :: ierr, ncall, i, j
   integer :: num_gpus = 0, igpu, num_gpus_env = -1
   integer :: size_buf = 60
   integer :: auto_scheme = 2
   integer :: print_level = 1

   integer :: executor
   integer :: buf_size = 1024_4*1024_4*1024_4 !desired Host argument buffer size in bytes

   integer :: tensor_datatype = 3 ! possible values are 1: real double; 2: real single; 3: complex double; 4: complex single
   character*2 :: datatype_label(4) = (/ 'R8', 'R4', 'C8', 'C4' /)


         real(8), allocatable :: eps_occ(:),eps_vir(:)
!        CCSD,CC2,... control and result variables
         real(8) :: e_core, scf_energy, mp2_energy, cc_energy, pert_energy(3)
         real(8) :: t_target_precision
         real(8) :: t1diag, t2diag
         integer :: ncycles = 1
         logical :: CCD, lambda, eom, eom_lr, eom_qr, qr_tpa
         logical :: CC2
         logical :: CCSDT, CCSDTQ, CCSDTQP
         logical :: do_perturbative
         integer :: iterative_level = 2, iterative_scheme = 21
         integer :: perturbative_level = 3 , perturbative_scheme = 30
         logical :: do1rdm_gs, do2rdm_gs
   
   nocc = 6
   nvir = 12
   ncall = 1

!  Make some noise so that we know we are here
   call print_date('Entered test_teq_auto code')
   call print_iterative_wavefunction_choice(iterative_scheme)

   do j = 3, 3, 2
   tensor_datatype = j
   call print_date('Using data type '//datatype_label(j)//' for tensors')
!  Initialize libraries
   executor = 1 ! tapp, can toggle to other executors as needed
   write(*,*) "    executor: ", executor
   call eti_executor_initialize(buf_size,executor)

   allocate (eps_occ(nocc))
   allocate (eps_vir(nvir))
   eps_occ = -0.5d0
   eps_vir = 0.5d0
   call print_date("Retrieved orbital energies")
!        create tensors needed in CC
         oo_dims = nocc
         vv_dims = nvir
         ierr=construct_tensor_t(comm_t%foo,tensor_datatype,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%fvv,tensor_datatype,vv_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%goo,tensor_datatype,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%hoo,tensor_datatype,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%gvv,tensor_datatype,vv_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%hvv,tensor_datatype,vv_dims,init_val=ZERO)
         if (.not.CCD) then
            ov_dims(1) = nocc
            ov_dims(2) = nvir
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=construct_tensor_t(comm_t%fov,tensor_datatype,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(comm_t%hov,tensor_datatype,ov_dims,init_val=ZERO)
            ierr=construct_tensor_t(comm_t%t1,tensor_datatype,vo_dims,init_val=ZERO)
!           initialize t3 amplitudes to zero (will later update them to starting values based on mp2 amplitudes)
            if (iterative_level.ge.3) then
               vvvooo_dims(1:3) = nvir
               vvvooo_dims(4:6) = nocc
               ierr=construct_tensor_t(comm_t%t3,tensor_datatype,vvvooo_dims,init_val=ZERO)
            end if
            if (iterative_level.ge.4) then
               vvvvoooo_dims(1:4) = nvir
               vvvvoooo_dims(5:8) = nocc
               ierr=construct_tensor_t(comm_t%t4,tensor_datatype,vvvvoooo_dims,init_val=ZERO)
            end if
            if (iterative_level.ge.5) then
               vvvvvooooo_dims(1:5) = nvir
               vvvvvooooo_dims(6:10) = nocc
               ierr=construct_tensor_t(comm_t%t5,tensor_datatype,vvvvvooooo_dims,init_val=ZERO)
            end if
         end if
         oooo_dims = nocc
         vovo_dims = nvir
         vovo_dims(2) = nocc
         vovo_dims(4) = nocc
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         vvoo_dims(3:4) = nocc
         ierr=construct_tensor_t(comm_t%a_int,tensor_datatype,oooo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%h_int,tensor_datatype,vovo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%tau,tensor_datatype,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t%t2,tensor_datatype,vvoo_dims,init_val=ZERO)
         if (ierr.eq.TALSH_SUCCESS) then
            call print_date('Succesfully allocated all tensors')
         else
            write(*,*) " ccdriver: error in tensor allocation, code=",ierr
            return
         end if

! creating integral tensors
   ierr=construct_tensor_t(int_t%foo,tensor_datatype,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"xx1 ",ierr
   ierr=construct_tensor_t(int_t%fov,tensor_datatype,(/nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx2 ",ierr
   ierr=construct_tensor_t(int_t%fvv,tensor_datatype,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx3 ",ierr
   ierr=construct_tensor_t(int_t%oooo,tensor_datatype,(/nocc,nocc,nocc,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"xx4 ",ierr
   ierr=construct_tensor_t(int_t%ooov,tensor_datatype,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx5 ",ierr
   ierr=construct_tensor_t(int_t%vovo,tensor_datatype,(/nvir,nocc,nvir,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"xx6 ",ierr
   ierr=construct_tensor_t(int_t%oovv,tensor_datatype,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx7 ",ierr
   ierr=construct_tensor_t(int_t%vovv,tensor_datatype,(/nvir,nocc,nvir,nvir/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"xx8 ",ierr
   ierr=construct_tensor_t(int_t%vvvv,tensor_datatype,(/nvir,nvir,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx9 ",ierr
   call print_date('Created 1/2-body Hamiltonian tensors')
   print *, '  #occupied',nocc 
   print *, '  #virtuals',nvir

   do i = 1, ncall
   call print_date('Starting t equations')
   ncycles = 1
   CCD=.false.
   t_target_precision = 1.0e-5
   ierr=initialize_tensor_t(comm_t%t1)
   ierr=initialize_tensor_t(comm_t%t2)
   ierr=initialize_tensor_t(comm_t%tau)
   if (iterative_level.ge.3) ierr=initialize_tensor_t(comm_t%t3)
   if (iterative_level.ge.4) ierr=initialize_tensor_t(comm_t%t4)
   if (iterative_level.ge.5) ierr=initialize_tensor_t(comm_t%t5)
   call solve_t_equations_automatic(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision,&
                                 cc_energy,CCD,iterative_level,iterative_scheme,auto_scheme, &
                                 ncycles,print_level,benchmarking_skip_diis=.true.) ! ,tensor_datatype)
   end do
                         
   call print_date('Done with t equations')

   call eti_executor_shutdown()
   call print_date('Done for data type '//datatype_label(j)//' for tensors')
   end do
                           
!  Make some noise so that we know we are leaving
   write(*,*) ""
   call print_date('Leaving test_sigma program')

end program 
