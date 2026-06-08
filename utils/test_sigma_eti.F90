program test_eti_sigma_ee

   use tensor_algebra
   use exacorr_tensor_interface
   use exacorr_datatypes
   use os_utils
   use exacorr_utils
   use eti_common_routines
   use eti_sigma_eomccsd
   use exacorr_scaling_parameters

   implicit none

   type(eti_intg_tens_t) :: int_t
   type(eti_intermediates_tens_t) :: interm_t
   type(eti_sa_tens_t) :: r_tensor, sigma_tensor

   integer :: nocc, nvir 
   integer :: ierr, ncall, i, j
   integer :: num_gpus = 0, igpu, num_gpus_env = -1
   integer :: size_buf = 60

   integer :: executor
   integer :: buf_size = 1024_4*1024_4*1024_4 !desired Host argument buffer size in bytes

   integer :: tensor_datatype = 2 ! possible values are 1: real double; 2: real single; 3: complex double; 4: complex single
   character*2 :: datatype_label(4) = (/ 'R4', 'R8', 'C4', 'C8' /)

   
   nocc = 4
   nvir = 20
   ncall = 1

!  Make some noise so that we know we are here
   call print_date('Entered test_sigma code')

   !  Initialize libraries
   executor = 1 ! tapp, can toggle to other executors as needed
   write(*,*) "    executor: ", executor
   call eti_executor_initialize(buf_size,executor)

   do j = 1, 4 
   tensor_datatype = j-1
   print *,''
   call print_date('Using data type '//datatype_label(j)//' for tensors')
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

   call print_date('Created EOMCC intermediates tensors')
!!!!!!
! setting up trial vector
   ierr=construct_tensor_t(r_tensor%s1,tensor_datatype,(/nvir,nocc/),init_val=ONE)
   if (ierr.ne.0) print *,"x26 ",ierr
   print *,'constructing s1'
   ierr=construct_tensor_t(sigma_tensor%s1,tensor_datatype,(/nvir,nocc/),init_val=MINUS_ONE)
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
      print *,"    info, datatype s1:",r_tensor%s1%datatype
      print *,"    info, norm(sigma1):",tensor_norm2(r_tensor%s1)
      print *,"    info, norm(sigma2):",tensor_norm2(r_tensor%s2)
      call print_date('After sigma vector evaluation')
   end do
   call print_date('End of sigma vector evaluations')
   print *, "  Number of calls to sigma vector routine",ncall

   print *, 'Cleaning up tensors'
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
   ierr=destroy_tensor_t(interm_t%w_oooo)
   ierr=destroy_tensor_t(interm_t%w_ovoo)
   ierr=destroy_tensor_t(interm_t%w_ooov)
   ierr=destroy_tensor_t(interm_t%w_vovo)
   ierr=destroy_tensor_t(interm_t%wbar_vovo)
   ierr=destroy_tensor_t(interm_t%w_vvvo)
   ierr=destroy_tensor_t(interm_t%w_vovv)
   ierr=destroy_tensor_t(interm_t%w_vvvv)
   ierr=destroy_tensor_t(interm_t%eta_oovv)
   ierr=destroy_tensor_t(interm_t%gs%r1)
   ierr=destroy_tensor_t(interm_t%gs%r2)
   ierr=destroy_tensor_t(r_tensor%s1)
   ierr=destroy_tensor_t(sigma_tensor%s1)
   ierr=destroy_tensor_t(r_tensor%s2)
   ierr=destroy_tensor_t(sigma_tensor%s2)


   call print_date('Done for data type '//datatype_label(j)//' for tensors')
   end do
   call eti_executor_shutdown()
                           
!  Make some noise so that we know we are leaving
   write(*,*) ""
   call print_date('Leaving test_sigma program')

end program 
