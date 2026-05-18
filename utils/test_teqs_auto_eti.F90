program test_eti_sigma_ee

   use tensor_algebra
   use exacorr_tensor_interface
   use exacorr_datatypes
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

   integer :: iterative_level, iterative_scheme
   integer :: perturbative_level, perturbative_scheme


   integer :: nocc, nvir 
   integer :: ierr, ncall, i, j
   integer :: num_gpus = 0, igpu, num_gpus_env = -1
   integer :: size_buf = 60

   integer :: executor
   integer :: buf_size = 1024_4*1024_4*1024_4 !desired Host argument buffer size in bytes

   integer :: tensor_datatype = 2 ! possible values are 1: real double; 2: real single; 3: complex double; 4: complex single
   character*2 :: datatype_label(4) = (/ 'R8', 'R4', 'tensor_datatype', 'C4' /)

   
   nocc = 4
   nvir = 20
   ncall = 10

!  Make some noise so that we know we are here
   call print_date('Entered test_sigma code')

   do j = 3, 3
   tensor_datatype = j
   call print_date('Using data type '//datatype_label(j)//' for tensors')
!  Initialize libraries
   executor = 1 ! tal-sh, can toggle to other executors as needed
   write(*,*) "    executor: ", executor
   call eti_executor_initialize(buf_size,executor)

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

! creating intermediate tensors
   ierr=construct_tensor_t(interm_t%fbar_oo,tensor_datatype,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x10 ",ierr
   ierr=construct_tensor_t(interm_t%goo,tensor_datatype,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x11 ",ierr

   ierr=construct_tensor_t(interm_t%fbar_ov,tensor_datatype,(/nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x12 ",ierr
   ierr=construct_tensor_t(interm_t%fbar_vv,tensor_datatype,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x13 ",ierr
   ierr=construct_tensor_t(interm_t%gvv,tensor_datatype,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x14 ",ierr

   ierr=construct_tensor_t(interm_t%w_oooo,tensor_datatype,(/nocc,nocc,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x15 ",ierr
   ierr=construct_tensor_t(interm_t%w_ovoo,tensor_datatype,(/nocc,nvir,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x16 ",ierr
   ierr=construct_tensor_t(interm_t%w_ooov,tensor_datatype,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x17 ",ierr
   ierr=construct_tensor_t(interm_t%w_vovo,tensor_datatype,(/nvir,nocc,nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x18 ",ierr
   ierr=construct_tensor_t(interm_t%wbar_vovo,tensor_datatype,(/nvir,nocc,nvir,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"x19 ",ierr
   ierr=construct_tensor_t(interm_t%w_vvvo,tensor_datatype,(/nvir,nvir,nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x20 ",ierr
   ierr=construct_tensor_t(interm_t%w_vovv,tensor_datatype,(/nvir,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x21 ",ierr
   ierr=construct_tensor_t(interm_t%w_vvvv,tensor_datatype,(/nvir,nvir,nvir,nvir/),init_val=ONE)
   if (ierr.ne.0) print *,"x22 ",ierr

   ierr=construct_tensor_t(interm_t%eta_oovv,tensor_datatype,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x23 ",ierr
!
   interm_t%nocc = nocc
   interm_t%nvir = nvir

   ierr=construct_tensor_t(interm_t%gs%r1,tensor_datatype,(/nvir,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x24 ",ierr
   ierr=construct_tensor_t(interm_t%gs%r2,tensor_datatype,(/nvir,nvir,nocc,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x25 ",ierr

   interm_t%sigma_rhs_p => eti_sigma_eomccsd_rhs_ee
   call solve_t_equations_automatic(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision,&
                                 cc_energy,CCD,iterative_level,iterative_scheme,exa_input%auto_scheme, &
                                 ncycles,exa_input%print_level)


   call print_date('Created EOMCC intermediates tensors')
!!!!!!
! setting up trial vector
   ierr=construct_tensor_t(r_tensor%s1,tensor_datatype,(/nvir,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x26 ",ierr
   ierr=construct_tensor_t(sigma_tensor%s1,tensor_datatype,(/nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x27 ",ierr
   ierr=construct_tensor_t(r_tensor%s2,tensor_datatype,(/nvir,nvir,nocc,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"x28 ",ierr
   ierr=construct_tensor_t(sigma_tensor%s2,tensor_datatype,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x29 ",ierr
   call print_date('Created EOMCC sigma and trial vector tensors')

!!!!!
! evaluating the sigma vector a number of times
   call print_date('Starting sigma vector evaluations')
   print *, "  Will carry out ",ncall," sigma vector evaluations" 
   do i = 1, ncall
      print *,"  Sigma vector evaluation #",i
      call interm_t%sigma_rhs_p(int_t,interm_t,r_tensor,sigma_tensor)
      print *,"    info, norm(sigma):",tensor_norm2(r_tensor%s1)+tensor_norm2(r_tensor%s2)
      call print_date('After sigma vector evaluation')
   end do
   call print_date('End of sigma vector evaluations')
   print *, "  Number of calls to sigma vector routine",ncall

   print *, 'Cleaning up tensors'
   ierr=destroy_tensor_t(sigma_tensor%s1)
   ierr=destroy_tensor_t(sigma_tensor%s2)
   ierr=destroy_tensor_t(r_tensor%s1)
   ierr=destroy_tensor_t(r_tensor%s1)
   ierr=destroy_tensor_t(interm_t%gs%r1)
   ierr=destroy_tensor_t(interm_t%gs%r2)

   ierr=destroy_tensor_t(interm_t%w_oooo)
   ierr=destroy_tensor_t(interm_t%w_ovoo)
   ierr=destroy_tensor_t(interm_t%w_ooov)
   ierr=destroy_tensor_t(interm_t%w_vovo)
   ierr=destroy_tensor_t(interm_t%wbar_vovo)
   ierr=destroy_tensor_t(interm_t%w_vvvo)
   ierr=destroy_tensor_t(interm_t%w_vovv)
   ierr=destroy_tensor_t(interm_t%w_vvvv)
   ierr=destroy_tensor_t(interm_t%eta_oovv)


   ierr=destroy_tensor_t(int_t%foo) 
   ierr=destroy_tensor_t(int_t%fov) 
   ierr=destroy_tensor_t(int_t%fvv) 
   ierr=destroy_tensor_t(int_t%oooo) 
   ierr=destroy_tensor_t(int_t%ooov) 
   ierr=destroy_tensor_t(int_t%vovo) 
   ierr=destroy_tensor_t(int_t%oovv) 
   ierr=destroy_tensor_t(int_t%vovv) 
   ierr=destroy_tensor_t(int_t%vvvv) 
   ierr=destroy_tensor_t(interm_t%fbar_oo) 
   ierr=destroy_tensor_t(interm_t%goo)
   ierr=destroy_tensor_t(interm_t%fbar_ov)
   ierr=destroy_tensor_t(interm_t%fbar_vv)
   ierr=destroy_tensor_t(interm_t%gvv) 

   call eti_executor_shutdown()
   call print_date('Done for data type '//datatype_label(j)//' for tensors')
   end do
                           
!  Make some noise so that we know we are leaving
   write(*,*) ""
   call print_date('Leaving test_sigma program')

end program 
