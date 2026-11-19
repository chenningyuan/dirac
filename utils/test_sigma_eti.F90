program test_eti_sigma_ee

   use tensor_algebra
   use exacorr_tensor_interface
   use exacorr_datatypes
   use os_utils
   use exacorr_utils
   use eti_common_routines
   use eti_sigma_eomccsd

   implicit none

   complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),TWO=(2.D0,0.D0), THREE=(3.0,0.0),   &
                            SIX=(6.0,0.0), MINUS_ONE=(-1.D0,0.D0), MINUS_TWO=(-2.D0,0.D0),       &
                            ONE_HALF=(0.5D0,0.D0), ONE_QUARTER=(0.25D0,0.D0),                    &
                            ONE_TWELFTH=(1.D0,0.D0)/12, MINUS_ONE_HALF=(-0.5D0,0.D0),            &
                            MINUS_ONE_QUARTER=(-0.25D0,0.D0), MINUS_ONE_TWELFTH=(-1.D0,0.D0)/12
   real(8), parameter    :: ONE_QUARTER_R=0.25D0, ONE_HALF_R=0.5D0, ONE_SIXTH_R=(1.D0/6.D0),     &
                            ONE_36_R=(1.D0/36.D0)

   type(eti_intg_tens_t) :: int_t
   type(eti_intermediates_tens_t) :: interm_t
   type(eti_sa_tens_t) :: r_tensor, sigma_tensor

   integer :: nocc, nvir 
   integer :: ierr, ncall, i
   integer :: num_gpus = 0, igpu, num_gpus_env = -1
   integer :: size_buf = 60

   integer :: executor

   integer :: buf_size = 1024_4*1024_4*1024_4 !desired Host argument buffer size in bytes

   nocc = 30
   nvir = 160
   ncall = 10

!  Make some noise so that we know we are here
   call print_date('Entered test_sigma code')

!  Initialize libraries
   executor = 2 ! tal-sh, can toggle to other executors as needed
   write(*,*) "    executor: ", executor
   call eti_executor_initialize(buf_size,executor)

! creating integral tensors
   ierr=construct_tensor_t(int_t%foo,C8,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"xx1 ",ierr
   ierr=construct_tensor_t(int_t%fov,C8,(/nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx2 ",ierr
   ierr=construct_tensor_t(int_t%fvv,C8,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx3 ",ierr
   ierr=construct_tensor_t(int_t%oooo,C8,(/nocc,nocc,nocc,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"xx4 ",ierr
   ierr=construct_tensor_t(int_t%ooov,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx5 ",ierr
   ierr=construct_tensor_t(int_t%vovo,C8,(/nvir,nocc,nvir,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"xx6 ",ierr
   ierr=construct_tensor_t(int_t%oovv,C8,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx7 ",ierr
   ierr=construct_tensor_t(int_t%vovv,C8,(/nvir,nocc,nvir,nvir/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"xx8 ",ierr
   ierr=construct_tensor_t(int_t%vvvv,C8,(/nvir,nvir,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx9 ",ierr
   call print_date('Created 1/2-body Hamiltonian tensors')
   print *, '  #occupied',nocc 
   print *, '  #virtuals',nvir

! creating intermediate tensors
   ierr=construct_tensor_t(interm_t%fbar_oo,C8,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x10 ",ierr
   ierr=construct_tensor_t(interm_t%goo,C8,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x11 ",ierr

   ierr=construct_tensor_t(interm_t%fbar_ov,C8,(/nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x12 ",ierr
   ierr=construct_tensor_t(interm_t%fbar_vv,C8,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x13 ",ierr
   ierr=construct_tensor_t(interm_t%gvv,C8,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x14 ",ierr

   ierr=construct_tensor_t(interm_t%w_oooo,C8,(/nocc,nocc,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x15 ",ierr
   ierr=construct_tensor_t(interm_t%w_ovoo,C8,(/nocc,nvir,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x16 ",ierr
   ierr=construct_tensor_t(interm_t%w_ooov,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x17 ",ierr
   ierr=construct_tensor_t(interm_t%w_vovo,C8,(/nvir,nocc,nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x18 ",ierr
   ierr=construct_tensor_t(interm_t%wbar_vovo,C8,(/nvir,nocc,nvir,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"x19 ",ierr
   ierr=construct_tensor_t(interm_t%w_vvvo,C8,(/nvir,nvir,nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x20 ",ierr
   ierr=construct_tensor_t(interm_t%w_vovv,C8,(/nvir,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x21 ",ierr
   ierr=construct_tensor_t(interm_t%w_vvvv,C8,(/nvir,nvir,nvir,nvir/),init_val=ONE)
   if (ierr.ne.0) print *,"x22 ",ierr

   ierr=construct_tensor_t(interm_t%eta_oovv,C8,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x23 ",ierr
!
   interm_t%nocc = nocc
   interm_t%nvir = nvir

   ierr=construct_tensor_t(interm_t%gs%r1,C8,(/nvir,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x24 ",ierr
   ierr=construct_tensor_t(interm_t%gs%r2,C8,(/nvir,nvir,nocc,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x25 ",ierr

   interm_t%sigma_rhs_p => eti_sigma_eomccsd_rhs_ee

   call print_date('Created EOMCC intermediates tensors')
!!!!!!
! setting up trial vector
   ierr=construct_tensor_t(r_tensor%s1,C8,(/nvir,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x26 ",ierr
   ierr=construct_tensor_t(sigma_tensor%s1,C8,(/nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x27 ",ierr
   ierr=construct_tensor_t(r_tensor%s2,C8,(/nvir,nvir,nocc,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"x28 ",ierr
   ierr=construct_tensor_t(sigma_tensor%s2,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
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

   call eti_executor_shutdown()

                           
!  Make some noise so that we know we are leaving
   write(*,*) ""
   call print_date('Leaving test_sigma program')

end program 
