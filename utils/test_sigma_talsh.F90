program test_talsh_sigma_ee

   use tensor_algebra
   use talsh
   use exacorr_datatypes
!  use intermediates
!  use intermediates_eom_ccsd
   use os_utils
   use exacorr_utils
   use talsh_common_routines
!  use mfsolver_dirac_input
   use talsh_sigma_eomccsd

   implicit none

   complex(8), parameter :: ZERO=(0.D0,0.D0),ONE=(1.D0,0.D0),TWO=(2.D0,0.D0), THREE=(3.0,0.0),   &
                            SIX=(6.0,0.0), MINUS_ONE=(-1.D0,0.D0), MINUS_TWO=(-2.D0,0.D0),       &
                            ONE_HALF=(0.5D0,0.D0), ONE_QUARTER=(0.25D0,0.D0),                    &
                            ONE_TWELFTH=(1.D0,0.D0)/12, MINUS_ONE_HALF=(-0.5D0,0.D0),            &
                            MINUS_ONE_QUARTER=(-0.25D0,0.D0), MINUS_ONE_TWELFTH=(-1.D0,0.D0)/12
   real(8), parameter    :: ONE_QUARTER_R=0.25D0, ONE_HALF_R=0.5D0, ONE_SIXTH_R=(1.D0/6.D0),     &
                            ONE_36_R=(1.D0/36.D0)

   type(talsh_intg_tens_t) :: int_t
   type(talsh_intermediates_tens_t) :: interm_t
   type(talsh_tens_t) :: r1_tensor, r2_tensor, sigma1_tensor, sigma2_tensor

   integer :: nocc, nvir 
   integer :: ierr, ncall, i
   integer :: num_gpus = 0, igpu, num_gpus_env = -1
   integer :: size_buf = 60

   integer(C_SIZE_T):: buf_size=1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
   integer(C_INT):: host_arg_max

   nocc = 30
   nvir = 160
   ncall = 10

! splash of talsh

!  Make some noise so that we know we are here
   call print_date('Entered test_sigma code')

   call talsh_get_number_gpus(num_gpus)
   num_gpus_env = get_environment_integer('TALSH_GPUS',-1)
   if (num_gpus_env.eq.-1) then
      print *,' Environment variable TALSH_GPUS undefined'
      num_gpus = max(num_gpus,num_gpus_env)
   else
      print *,' Environment variable TALSH_GPUS defined, set to:',num_gpus_env
      num_gpus = min(num_gpus,num_gpus_env)
   end if
   print *,' GPUs in use:',num_gpus
   
!  Initialize libraries
   buf_size=buf_size*size_buf
   ierr=talsh_init(buf_size,host_arg_max,gpu_list=(/(igpu,igpu=0,num_gpus-1)/))
   call print_date('Initialized talsh library')
   write(*,'("  Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,buf_size,host_arg_max

! creating integral tensors
   ierr=talsh_tensor_construct(int_t%foo,C8,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"xx1 ",ierr
   ierr=talsh_tensor_construct(int_t%fov,C8,(/nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx2 ",ierr
   ierr=talsh_tensor_construct(int_t%fvv,C8,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx3 ",ierr
   ierr=talsh_tensor_construct(int_t%oooo,C8,(/nocc,nocc,nocc,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"xx4 ",ierr
   ierr=talsh_tensor_construct(int_t%ooov,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx5 ",ierr
   ierr=talsh_tensor_construct(int_t%vovo,C8,(/nvir,nocc,nvir,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"xx6 ",ierr
   ierr=talsh_tensor_construct(int_t%oovv,C8,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx7 ",ierr
   ierr=talsh_tensor_construct(int_t%vovv,C8,(/nvir,nocc,nvir,nvir/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"xx8 ",ierr
   ierr=talsh_tensor_construct(int_t%vvvv,C8,(/nvir,nvir,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"xx9 ",ierr
   call print_date('Created 1/2-body Hamiltonian tensors')
   print *, '  #occupied',nocc 
   print *, '  #virtuals',nvir

! creating intermediate tensors
   ierr=talsh_tensor_construct(interm_t%fbar_oo,C8,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x10 ",ierr
   ierr=talsh_tensor_construct(interm_t%goo,C8,(/nocc,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x11 ",ierr

   ierr=talsh_tensor_construct(interm_t%fbar_ov,C8,(/nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x12 ",ierr
   ierr=talsh_tensor_construct(interm_t%fbar_vv,C8,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x13 ",ierr
   ierr=talsh_tensor_construct(interm_t%gvv,C8,(/nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x14 ",ierr

   ierr=talsh_tensor_construct(interm_t%w_oooo,C8,(/nocc,nocc,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x15 ",ierr
   ierr=talsh_tensor_construct(interm_t%w_ovoo,C8,(/nocc,nvir,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x16 ",ierr
   ierr=talsh_tensor_construct(interm_t%w_ooov,C8,(/nocc,nocc,nocc,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x17 ",ierr
   ierr=talsh_tensor_construct(interm_t%w_vovo,C8,(/nvir,nocc,nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x18 ",ierr
   ierr=talsh_tensor_construct(interm_t%wbar_vovo,C8,(/nvir,nocc,nvir,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"x19 ",ierr
   ierr=talsh_tensor_construct(interm_t%w_vvvo,C8,(/nvir,nvir,nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x20 ",ierr
   ierr=talsh_tensor_construct(interm_t%w_vovv,C8,(/nvir,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x21 ",ierr
   ierr=talsh_tensor_construct(interm_t%w_vvvv,C8,(/nvir,nvir,nvir,nvir/),init_val=ONE)
   if (ierr.ne.0) print *,"x22 ",ierr

   ierr=talsh_tensor_construct(interm_t%eta_oovv,C8,(/nocc,nocc,nvir,nvir/),init_val=ZERO)
   if (ierr.ne.0) print *,"x23 ",ierr
!
   interm_t%nocc = nocc
   interm_t%nvir = nvir

   ierr=talsh_tensor_construct(interm_t%t1,C8,(/nvir,nocc/),init_val=MINUS_ONE)
   if (ierr.ne.0) print *,"x24 ",ierr
   ierr=talsh_tensor_construct(interm_t%t2,C8,(/nvir,nvir,nocc,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x25 ",ierr

   interm_t%sigma_rhs_p => talsh_sigma_eomccsd_rhs_ee

   call print_date('Created EOMCC intermediates tensors')
!!!!!!
! setting up trial vector
   ierr=talsh_tensor_construct(r1_tensor,C8,(/nvir,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x26 ",ierr
   ierr=talsh_tensor_construct(sigma1_tensor,C8,(/nvir,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x27 ",ierr
   ierr=talsh_tensor_construct(r2_tensor,C8,(/nvir,nvir,nocc,nocc/),init_val=ONE_QUARTER)
   if (ierr.ne.0) print *,"x28 ",ierr
   ierr=talsh_tensor_construct(sigma2_tensor,C8,(/nvir,nvir,nocc,nocc/),init_val=ZERO)
   if (ierr.ne.0) print *,"x29 ",ierr
   call print_date('Created EOMCC sigma and trial vector tensors')

!!!!!
! evaluating the sigma vector a number of times
   call print_date('Starting sigma vector evaluations')
   print *, "  Will carry out ",ncall," sigma vector evaluations" 
   do i = 1, ncall
      print *,"  Sigma vector evaluation #",i
      call interm_t%sigma_rhs_p(int_t,interm_t,r1_tensor,r2_tensor,sigma1_tensor,sigma2_tensor)
      print *,"    info, norm(sigma):",tensor_norm2(r1_tensor)+tensor_norm2(r2_tensor)
      call print_date('After sigma vector evaluation')
   end do
   call print_date('End of sigma vector evaluations')
   print *, "  Number of calls to sigma vector routine",ncall

   ierr=talsh_stats()
   if (ierr.ne.0) print *,"x30 ",ierr
   ierr=talsh_shutdown()
   if (ierr.ne.0) print *,"x31 ",ierr
                           
!  Make some noise so that we know we are leaving
   write(*,*) ""
   call print_date('Leaving test_sigma program')

end program 
