module eti_drivers_gs_ccsd
!This module contains the routines needed to calculate ae coupled cluster wave function

        use tensor_algebra
        use exacorr_tensor_interface
        use exacorr_datatypes
        use exacorr_cfg, only: exacc_input
        use exacorr_utils
        use intermediates_cc
        use intermediates_lambda
        use eti_common_routines
        use exacorr_tensor_interface
        use exacorr_printout
        use exacorr_scaling_parameters
        use eti_diis
        use, intrinsic:: ISO_C_BINDING

        implicit none

        type(tensor_t) :: one_tensor

        private

        
        public solve_t_equations
        public solve_CC2_t_equations
        public solve_lambda_equations
        public solve_CC2_lambda_equations


       contains
!-------------------------------------------------------------------------------------------

        subroutine solve_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                ccsd_energy,CCD,iterative_level,iterative_scheme,ncycles,print_level)

!        Routine for solving the t amplitude equations

         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(eti_intg_tens_t), intent(inout) :: int_t
!        intermediate tensors (same notation as in relccsd)
         !type(tensor_t) :: hoo_int_tensor, hvv_int_tensor
         type(tensor_t) :: tau2_tensor
         type(tensor_t) :: cint_tensor
!        auxilliary tensors for t1 equation
         type(tensor_t) :: hov_aux_tensor, hoo_aux_tensor
!        solution tensors
         type(tensor_t) :: s1_tensor, s2_tensor
         type(tensor_t) :: vvvo_aux_tensor
!        scalars (need to be defined as tensor types)
         type(tensor_t) :: result_tensor
         integer(INTD)      :: result_dims(1)
         !Tensor dimensions
         integer(INTD), dimension(2) :: oo_dims, ov_dims, vo_dims
         integer(INTD), dimension(4) :: vvoo_dims, ovoo_dims, vvvo_dims
         !Orbital energies for scaling
         real(8), intent(in)   :: eps_occ(:), eps_vir(:)
         integer :: iccsd
         integer, intent(in) :: ncycles
         real(8), intent(in) :: t_target_precision
         real(8) :: t_convergence
         real(8), intent(out) :: ccsd_energy
         logical, intent(in)  :: CCD
         integer, intent(in)  :: iterative_level
         integer, intent(in)  :: iterative_scheme
         integer, intent(in)  :: print_level
         integer (INTD)       :: ierr
         type(eti_precision_t):: eti_precision

         eti_precision%storage_t = comm_t%t2%datatype

!        Initialize scalars that are to be used as tensors in contractions
         result_dims = 1
         ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,result_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(result_tensor,eti_precision%storage_t,result_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,result_dims)

         if (.not.CCD) then
            oo_dims = nocc
            ierr=construct_tensor_t(hoo_aux_tensor,eti_precision%storage_t,oo_dims,init_val=ZERO)
            vo_dims(1) = nvir
            vo_dims(2) = nocc
            ierr=construct_tensor_t(s1_tensor,eti_precision%storage_t,vo_dims,init_val=ZERO)
         end if
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         ierr=construct_tensor_t(hov_aux_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=construct_tensor_t(s2_tensor,eti_precision%storage_t,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(tau2_tensor,eti_precision%storage_t,vvoo_dims,init_val=ZERO)
         ovoo_dims(1) = nocc
         ovoo_dims(2) = nvir
         ovoo_dims(3:4) = nocc
         ierr=construct_tensor_t(cint_tensor,eti_precision%storage_t,ovoo_dims,init_val=ZERO)

         do iccsd = 1, ncycles

          call get_tau(comm_t,CCD)

! aspg, TODO: harmonize t1,t2,t3,t4 etc with general excitation rank setup
          if (print_level.gt.6) then
            if (.not.CCD) write(*,*) "t1    = ",tensor_norm2(comm_t%t1)
            write(*,*) "t2    = ",tensor_norm2(comm_t%t2)*ONE_QUARTER_R
            write(*,*) "Tau   = ",tensor_norm2(comm_t%tau)*ONE_QUARTER_R
          end if

          if (print_level.gt.2) then
            if (.not.CCD) write(*,*) "s1    = ",tensor_norm2(s1_tensor)
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)*ONE_QUARTER_R
          end if

          if (ierr.ne.ETI_SUCCESS) then
             write(*,*) " ccdriver: error in calculation of common blocks, code=",ierr
             return
          end if

!NOTE: contractions marked with "validated" gave numerically significant contributions and could be validated by comparison
!to the reference implementation, the ones marked "moretests" were too small to decide whether it was correct or not.
!        calculate intermediaries

          ierr=initialize_tensor_t(s2_tensor)

!         The vvvv contribution
          ierr=t_einsum("S(a,b,i,j)+=V(a,b,c,d)*T(c,d,i,j)",s2_tensor,int_t%vvvv,comm_t%tau,scale=ONE_HALF) ! validated o(2)v(4)

!         vv ! validated
          call get_hvv(comm_t,int_t%oovv)
          call get_gvv(comm_t,int_t%vovv,CCD)

!         vovo ! moretests
          ierr=initialize_tensor_t(tau2_tensor)
          ierr=t_einsum("T(a,d,l,i)+=R(a,d,l,i)",tau2_tensor,comm_t%t2,scale=ONE_HALF)
          if (.not.CCD) ierr=t_einsum("T(a,d,l,i)+=T(a,l)*T(d,i)",tau2_tensor,comm_t%t1,comm_t%t1)                       ! o(2)v(2)
          call get_h(comm_t,int_t,CCD)
          ierr=t_einsum("H(a,k,c,i)+=V(k,l,c,d)*T(a,d,l,i)",comm_t%h_int,int_t%oovv,tau2_tensor)           ! o(3)v(3)

!         This term appears four times as we need to permute both i and j and a and b
          ierr=t_einsum("S(a,b,i,j)+=H(a,k,c,i)*T(c,b,j,k)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)
          ierr=t_einsum("S(a,b,i,j)+=H(a,k,c,j)*T(c,b,k,i)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)
          ierr=t_einsum("S(a,b,i,j)+=H(b,k,c,i)*T(a,c,j,k)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)
          ierr=t_einsum("S(a,b,i,j)+=H(b,k,c,j)*T(a,c,k,i)",s2_tensor,comm_t%h_int,comm_t%t2) ! validated o(3)v(3)

!         ovoo ! validated
          ierr=initialize_tensor_t(cint_tensor)
          ierr=t_einsum("C(k,b,i,j)+=V(b,k,c,d)*T(c,d,i,j)",cint_tensor,int_t%vovv,comm_t%tau,scale=MINUS_ONE_HALF) ! o(3)v(3)
          if (.not.CCD) then
            ierr=t_einsum("C(k,b,i,j)+=V(b,k,c,i)*T(c,j)",cint_tensor,int_t%vovo,comm_t%t1)                           ! o(3)v(2)
            ierr=t_einsum("C(k,b,i,j)+=V(b,k,c,j)*T(c,i)",cint_tensor,int_t%vovo,comm_t%t1,scale=MINUS_ONE)           ! o(3)v(2)
          end if
          ierr=t_einsum("C(k,b,i,j)+=V+(i,j,k,b)*K()",cint_tensor,int_t%ooov,one_tensor)                            ! o(3)v(1)

!         oo ! validated
          call get_hoo(comm_t,int_t%oovv)
          call get_goo(comm_t,int_t%ooov,CCD)

!         ov ! validated
          if (.not.CCD) then
            call get_hov(comm_t,int_t%oovv)
            ierr=initialize_tensor_t(hov_aux_tensor)
            ierr=t_einsum("L(k,c)+=H(k,c)",hov_aux_tensor,comm_t%hov)                       ! copy common term
            ierr=t_einsum("H(k,c)+=F(k,c)",hov_aux_tensor,comm_t%fov,scale=MINUS_TWO)       ! subtract fock matrix twice
          end if

!         oooo ! validated
          call get_a(comm_t,int_t,CCD)

          if (ierr.ne.ETI_SUCCESS) then
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
          ierr=t_einsum("S(a,b,i,j)+=G(k,j)*T(a,b,i,k)",s2_tensor,comm_t%goo,comm_t%t2)  ! validated o(3)v(2)
          ierr=t_einsum("S(a,b,i,j)+=G(k,i)*T(a,b,k,j)",s2_tensor,comm_t%goo,comm_t%t2)  ! validated o(3)v(2)
          if (.not.CCD) then
!           Again two times, also note the complex conjugation as we use vovv instead of vvvo
            ierr=t_einsum("S(a,b,i,j)+=V+(c,j,a,b)*T(c,i)",s2_tensor,int_t%vovv,comm_t%t1)  ! validated o(2)v(3)
            ierr=t_einsum("S(a,b,i,j)+=V+(c,i,a,b)*T(c,j)",s2_tensor,int_t%vovv,comm_t%t1,scale=MINUS_ONE)  ! validated o(2)v(3)
          end if

!         This term appears twice as we need to permute a and b
          ierr=t_einsum("S(a,b,i,j)+=G(a,c)*T(c,b,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! validated o(2)v(3)
          ierr=t_einsum("S(a,b,i,j)+=G(b,c)*T(a,c,i,j)",s2_tensor,comm_t%gvv,comm_t%t2)  ! validated o(2)v(3)
          if (.not.CCD) then
!           Again two times, with the C intermediate (note: in RELCCSD this is done differently)
            ierr=t_einsum("S(a,b,i,j)+=C(k,b,i,j)*T(a,k)",s2_tensor,cint_tensor,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
            ierr=t_einsum("S(a,b,i,j)+=C(k,a,i,j)*T(b,k)",s2_tensor,cint_tensor,comm_t%t1) ! o(3)v(2)
          end if

!         The oooo contribution
          ierr=t_einsum("S(a,b,i,j)+=A(k,l,i,j)*T(a,b,k,l)",s2_tensor,comm_t%a_int,comm_t%tau,scale=ONE_HALF) !validated o(4)v(2)

!         ierr=t_einsum("S(a,b,i,j)+=V+(i,j,a,b)",s2_tensor,oovv_tensor)
!         work-around because complex conjugation does not work in the above statement
          ierr=t_einsum("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)    ! validated o(2)v(2)
!         end of work-around

!         calculate s1 tensor
          if (.not.CCD) then
            ierr=initialize_tensor_t(s1_tensor)
!           ierr=t_einsum("S(a,i)+=F+(i,a)",s1_tensor,fov_tensor)
!           work-around because complex conjugation does not work in the above statement
            ierr=t_einsum("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
!           end of work-around
            ierr=t_einsum("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,scale=ONE_HALF) !validated o(2)v(3)
            ierr=t_einsum("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                     !validated o(1)v(2)
            ierr=t_einsum("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                     !validated o(2)v(1)
            ierr=t_einsum("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)                     !validated o(2)v(2)
            ierr=initialize_tensor_t(hoo_aux_tensor)                                                                     !moretests
            ierr=t_einsum("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)                !moretests o(2)v(1)
            ierr=t_einsum("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)                     !moretests o(2)v(1)
            ierr=t_einsum("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,scale=ONE_HALF) !validated o(3)v(2)
            ierr=t_einsum("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,scale=MINUS_ONE)     !validated o(2)v(2)
          end if

          if (print_level.gt.6) then
            if (.not.CCD) write(*,*) "s1    = ",tensor_norm2(s1_tensor)
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)*ONE_QUARTER_R
          end if

          if (ierr.ne.ETI_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

          if (print_level.gt.2) then
            if (.not.CCD) write(*,*) "s1    = ",tensor_norm2(s1_tensor)
            write(*,*) "s2    = ",tensor_norm2(s2_tensor)*ONE_QUARTER_R
          end if

!         prepare for next iteration
          if (.not.CCD) then
            ierr=initialize_tensor_t(comm_t%t1)
            ierr=t_einsum("T1(a,i)+=S1(a,i)",comm_t%t1,s1_tensor)
          end if
          ierr=initialize_tensor_t(comm_t%t2)
          ierr=t_einsum("T2(a,b,i,j)+=S2(a,b,i,j)",comm_t%t2,s2_tensor)

          if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,t1_tensor=comm_t%t1,t2_tensor=comm_t%t2)
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%t2)
          end if
         !  call print_tensor(comm_t%t2, 1.0d-18, "t2 in solve_t_eqs before diis")

          if (.not.CCD) then
            call diis (t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,t_convergence=t_convergence)
          else
            call diis (t2_tensor=comm_t%t2,t_convergence=t_convergence)
          end if

         ! call print_tensor(comm_t%t2, 1.0d-18, "t2 in solve_t_eqs, after diis")

          if (.not.CCD) then
            ierr=initialize_tensor_t(result_tensor)
            ierr=t_einsum("E()+=T(a,i)*F(i,a)",result_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
            ccsd_energy=result_tens(1)
          else
            ccsd_energy = 0.D0
          end if

!         Update tau
          call get_tau(comm_t,CCD)

          ierr=initialize_tensor_t(result_tensor)
          ierr=t_einsum("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER_R

          call print_iteration(iccsd,t_convergence,ccsd_energy,print_level)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "-----------------------------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         write(*,*) ""

         if (.not.CCD) then
            ierr=destroy_tensor_t(s1_tensor)
         end if
         ierr=destroy_tensor_t(s2_tensor)
         ierr=destroy_tensor_t(tau2_tensor)
         ierr=destroy_tensor_t(cint_tensor)
         if (.not.CCD) then
            ierr=destroy_tensor_t(hoo_aux_tensor)
            ierr=destroy_tensor_t(hov_aux_tensor)
         end if
         ierr=destroy_tensor_t(result_tensor)
         ierr=destroy_tensor_t(one_tensor)
         ! call print_tensor(comm_t%t2, 1.0d-18, "t2 in solve_t_eqs, leaving")

        end subroutine

        subroutine solve_CC2_t_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir,t_target_precision, &
                                ccsd_energy,CCD,ncycles,print_level)

!        Routine for solving the t amplitude equations


         complex(8), pointer :: result_tens(:)
         type(C_PTR):: body_p
         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(eti_intg_tens_t), intent(inout) :: int_t
!        intermediate tensors (same notation as in relccsd)
         type(tensor_t) :: Wabef_tensor,Wmnij_tensor     !L
         type(tensor_t) :: Wabej_tensor,Wmbij_tensor     !L

! Add tensor usefull to avoid tauCC2_tensor 
         type(tensor_t) :: Hovvv_tensor,Hvovv_tensor     !L
         type(tensor_t) :: Hvvov_tensor,Hvvvo_tensor     !L
         type(tensor_t) :: Hooov_tensor,Hoovo_tensor
         type(tensor_t) :: Hvooo_tensor,Hovoo_tensor
         
         type(tensor_t) :: cint_ovoo_tensor
         type(tensor_t) :: cint_vvvo_tensor
!        auxilliary tensors for t1 equation
         type(tensor_t) :: hov_aux_tensor, hoo_aux_tensor
!        solution tensors
         type(tensor_t) :: s1_tensor, s2_tensor
!        scalars (need to be defined as tensor types)
         type(tensor_t) :: result_tensor
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
         type(eti_precision_t):: eti_precision

         eti_precision%storage_t = comm_t%t2%datatype


         if (print_level.gt.8) debug_local = .true.

!        Initialize scalars that are to be used as tensors in contractions
         result_dims = 1
         ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,result_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(result_tensor,eti_precision%storage_t,result_dims(1:0),init_val=ZERO)
         ierr=direct_data_pointer(result_tensor,body_p)
         call c_f_pointer(body_p,result_tens,result_dims)


         oo_dims = nocc
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         vo_dims(1) = nvir
         vo_dims(2) = nocc

         if (.not.CCD) then
          ierr=construct_tensor_t(hoo_aux_tensor,eti_precision%storage_t,oo_dims,init_val=ZERO)
          ierr=construct_tensor_t(hov_aux_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
         end if

         ierr=construct_tensor_t(s1_tensor,eti_precision%storage_t,vo_dims,init_val=ZERO)
         vvoo_dims(1:2) = nvir
         vvoo_dims(3:4) = nocc
         ierr=construct_tensor_t(s2_tensor,eti_precision%storage_t,vvoo_dims,init_val=ZERO)

         vvvv_dims(1:4) = nvir  !L
         oooo_dims(1:4) = nocc  !L
         ierr=construct_tensor_t(Wabef_tensor,eti_precision%storage_t,vvvv_dims,init_val=ZERO) !L
         ierr=construct_tensor_t(Wmnij_tensor,eti_precision%storage_t,oooo_dims,init_val=ZERO) !L

         vvvo_dims(1:3) = nvir !L
         vvvo_dims(4)   = nocc !L
         ierr=construct_tensor_t(cint_vvvo_tensor,eti_precision%storage_t,vvvo_dims,init_val=ZERO)

         ovoo_dims(1) = nocc
         ovoo_dims(2) = nvir
         ovoo_dims(3:4) = nocc
         ierr=construct_tensor_t(cint_ovoo_tensor,eti_precision%storage_t,ovoo_dims,init_val=ZERO)

         ierr=construct_tensor_t(Wabej_tensor,eti_precision%storage_t,vvvo_dims,init_val=ZERO)
         ierr=construct_tensor_t(Wmbij_tensor,eti_precision%storage_t,ovoo_dims,init_val=ZERO)

! Creation of the 4 intermediaires tensor : Hvovv , Hovvv, Hvvov, Hvvvo
         vovv_dims(1:4) = nvir
         vovv_dims(2)   = nocc
         ovvv_dims(1:4) = nvir
         ovvv_dims(1)   = nocc
         vvov_dims(1:4) = nvir
         vvov_dims(3)   = nocc
         ierr=construct_tensor_t(Hovvv_tensor,eti_precision%storage_t,ovvv_dims,init_val=ZERO)
         ierr=construct_tensor_t(Hvovv_tensor,eti_precision%storage_t,vovv_dims,init_val=ZERO)
         ierr=construct_tensor_t(Hvvvo_tensor,eti_precision%storage_t,vvvo_dims,init_val=ZERO)
         ierr=construct_tensor_t(Hvvov_tensor,eti_precision%storage_t,vvov_dims,init_val=ZERO)
         ooov_dims(1:4) = nocc
         ooov_dims(4)   = nvir
         oovo_dims(1:4) = nocc
         oovo_dims(3)   = nvir
         vooo_dims(1:4) = nocc
         vooo_dims(1)   = nvir
         ierr=construct_tensor_t(Hooov_tensor,eti_precision%storage_t,ooov_dims,init_val=ZERO)
         ierr=construct_tensor_t(Hoovo_tensor,eti_precision%storage_t,oovo_dims,init_val=ZERO)
         ierr=construct_tensor_t(Hovoo_tensor,eti_precision%storage_t,ovoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(Hvooo_tensor,eti_precision%storage_t,vooo_dims,init_val=ZERO)

! resorted integrals nedded below

         ierr=initialize_tensor_t(cint_vvvo_tensor)
         ierr=t_einsum("C(b,c,a,i)+=V+(a,i,b,c)*K()",cint_vvvo_tensor,int_t%vovv,one_tensor) !Inversion, warning : no more + in +=

         ierr=initialize_tensor_t(cint_ovoo_tensor)
         ierr=t_einsum("C(k,a,i,j)+=V+(i,j,k,a)*K()",cint_ovoo_tensor,int_t%ooov,one_tensor) !Inversion, warning : no more + in +=

         do iccsd = 1, ncycles


          ierr=initialize_tensor_t(s2_tensor)
          if (debug_local) print*, "norm s2 tensor at the start of CC2 t2 ",tensor_norm2(s2_tensor)

          if (iccsd==1) then
                 call get_tau(comm_t,CCD)
                 if (debug_local) print *, " call get_tau "
          end if

!  Start s2_tensor
!         Diagram 15 : s2_{ij}^{ab} += < ab || ij >
          ierr=t_einsum("S(a,b,i,j)+=V+(i,j,a,b)*K()",s2_tensor,int_t%oovv,one_tensor)
          if (debug_local) print*, "norm s2 tensor after diagram 15 ",tensor_norm2(s2_tensor)

          ierr=initialize_tensor_t(Wabej_tensor,ZERO)
          ierr=initialize_tensor_t(Wmbij_tensor,ZERO)           


          ierr=t_einsum("W(a,b,e,j)+=C(a,b,e,j)",Wabej_tensor,cint_vvvo_tensor,scale=ONE)
          ierr=t_einsum("W(a,b,e,j)+=T(a,m)*V(b,m,e,j)",Wabej_tensor,comm_t%t1,int_t%vovo,scale=ONE_HALF)
          ierr=t_einsum("W(a,b,e,j)+=T(b,m)*V(a,m,e,j)",Wabej_tensor,comm_t%t1,int_t%vovo,scale=MINUS_ONE_HALF)
       
          ierr=t_einsum("S(a,b,i,j)=T(e,i)*W(a,b,e,j)",s2_tensor,comm_t%t1,Wabej_tensor,scale=ONE) 
          ierr=t_einsum("S(a,b,i,j)=T(e,j)*W(a,b,e,i)",s2_tensor,comm_t%t1,Wabej_tensor,scale=MINUS_ONE)
 
         
          ierr=t_einsum("W(m,b,i,j)+=C(m,b,i,j)",Wmbij_tensor,cint_ovoo_tensor,scale=ONE)
          ierr=t_einsum("W(m,b,i,j)+=T(e,i)*V(b,m,e,j)",Wmbij_tensor,comm_t%t1,int_t%vovo,scale=MINUS_ONE_HALF)
          ierr=t_einsum("W(m,b,i,j)+=T(e,j)*V(b,m,e,i)",Wmbij_tensor,comm_t%t1,int_t%vovo,scale=ONE_HALF)
       
          ierr=t_einsum("S(a,b,i,j)=T(a,m)*W(m,b,i,j)",s2_tensor,comm_t%t1,Wmbij_tensor,scale=MINUS_ONE) 
          ierr=t_einsum("S(a,b,i,j)=T(b,m)*W(m,a,i,j)",s2_tensor,comm_t%t1,Wmbij_tensor,scale=ONE)

!         additional CC2 W intermediates
!
!         a. W_{ij}^{mn} = V_{ij}^{mn} + P\_(ij)\sum_{e} t^e_i V_{ej}^{mn} + \frac{1}{4}P\_(ij)\sum_{ef} t^e_i t^f_j V_{ef}^{mn}
          ierr=initialize_tensor_t(Wmnij_tensor)
          ierr=t_einsum("Wmnij(m,n,i,j)+=V(m,n,i,j)*K()",Wmnij_tensor,int_t%oooo,one_tensor)
          ierr=t_einsum("Wmnij(m,n,i,j)+=V(m,n,i,e)*T(e,j)",Wmnij_tensor,int_t%ooov,comm_t%t1)
          ierr=t_einsum("Wmnij(m,n,i,j)+=V(m,n,j,e)*T(e,i)",Wmnij_tensor,int_t%ooov,comm_t%t1,scale=MINUS_ONE)
          
          ierr=initialize_tensor_t(Hoovo_tensor)
          ierr=t_einsum("Hoovo(m,n,e,j)+=V(m,n,e,f)*T(f,j)",& 
               &       Hoovo_tensor,int_t%oovv,comm_t%t1)
          ierr=t_einsum("Wmnij(m,n,i,j)+=Hoovo(m,n,e,j)*T(e,i)",&
               &       Wmnij_tensor,Hoovo_tensor,comm_t%t1,scale=ONE_QUARTER)
          ierr=initialize_tensor_t(Hooov_tensor)
          ierr=t_einsum("Hooov(m,n,j,f)+=V(m,n,e,f)*T(e,j)",& 
               &       Hooov_tensor,int_t%oovv,comm_t%t1)
          ierr=t_einsum("Wmnij(m,n,i,j)+=Hooov(m,n,j,f)*T(f,i)",&
               &       Wmnij_tensor,Hooov_tensor,comm_t%t1,scale=MINUS_ONE_QUARTER)

!         b. W_{ij}^{ab} = V_{ef}^{ab} - P\_(ab)\sum_{m} t_m^b V_{ef}^{am} + \frac{1}{4}P\_(ab)\sum_{ef} t^e_i t^f_j V_{ef}^{mn}
          ierr=initialize_tensor_t(Wabef_tensor)
          ierr=t_einsum("Wabef(a,b,e,f)+=V(a,b,e,f)*K()",Wabef_tensor,int_t%vvvv,one_tensor)
          ierr=t_einsum("Wabef(a,b,e,f)+=V(a,m,e,f)*T(b,m)",Wabef_tensor,int_t%vovv,comm_t%t1,scale=MINUS_ONE)
          ierr=t_einsum("Wabef(a,b,e,f)+=V(b,m,e,f)*T(a,m)",Wabef_tensor,int_t%vovv,comm_t%t1)

          ierr=initialize_tensor_t(Hvovv_tensor)
          ierr=initialize_tensor_t(Hovvv_tensor)
          ierr=t_einsum("Hvovv(a,n,e,f)+=V(m,n,e,f)*T(a,m)",& 
               & Hvovv_tensor,int_t%oovv,comm_t%t1,scale=ONE)        !Operation cost : v^3*o^2
          ierr=t_einsum("Wabef(a,b,e,f)+=Hvovv(a,n,e,f)*T(b,n)",& 
               & Wabef_tensor,Hvovv_tensor,comm_t%t1,scale=ONE_QUARTER)  !Operation cost : v^4*o
          ierr=t_einsum("Hovvv(m,a,e,f)+=V(m,n,e,f)*T(a,n)",& 
               & Hovvv_tensor,int_t%oovv,comm_t%t1,scale=ONE)        ! Operation cost :v^3*o^2
          ierr=t_einsum("Wabef(a,b,e,f)+=Hovvv(m,a,e,f)*T(b,m)",& 
               & Wabef_tensor,Hovvv_tensor,comm_t%t1,scale=MINUS_ONE_QUARTER)  ! Operation cost : v^4*o
          if (debug_local) print*, "norm Wabefi_CC2 new",tensor_norm2(Wabef_tensor)


!         Diagram 21 : s2_{ij}^{ab} += \frac{1}{2}\sum_{mn} tauCC2^{mn}_{ab} W_{ij}^{mn}

          ierr=initialize_tensor_t(Hvooo_tensor)
          ierr=t_einsum("H(a,n,i,j)+=Wmnij(m,n,i,j)*T(a,m)",&
               & Hvooo_tensor,Wmnij_tensor,comm_t%t1,scale=ONE)
          ierr=t_einsum("S(a,b,i,j)+=H(a,n,i,j)*T(b,n)",&
               & s2_tensor,Hvooo_tensor,comm_t%t1,scale=ONE_HALF)
          ierr=initialize_tensor_t(Hovoo_tensor)
          ierr=t_einsum("H(m,a,i,j)+=Wmnij(m,n,i,j)*T(a,n)",&
               & Hovoo_tensor,Wmnij_tensor,comm_t%t1,scale=ONE)
          ierr=t_einsum("S(a,b,i,j)+=H(m,a,i,j)*T(b,m)",&
               & s2_tensor,Hovoo_tensor,comm_t%t1,scale=MINUS_ONE_HALF)
          if (debug_local) print*, "norm s2 tensor after diagram 21 new ",tensor_norm2(s2_tensor)

!         Diagram 22 : s2_{ij}^{ab} += \frac{1}{2}\sum_{ef} tauCC2^{ij}_{ef} W_{ef}^{ab}

          ierr=initialize_tensor_t(Hvvov_tensor)
          ierr=initialize_tensor_t(Hvvvo_tensor)
          ierr=t_einsum("Hvvov(a,b,i,f)+=Wabef(a,b,e,f)*T(e,i)",& 
               & Hvvov_tensor,Wabef_tensor,comm_t%t1,scale=ONE)    ! Operationcost : v^4*o
          ierr=t_einsum("S(a,b,i,j)+= Hvvov(a,b,i,f)*T(f,j)",& 
               & s2_tensor,Hvvov_tensor,comm_t%t1,scale=ONE_HALF)    !  Operation cost : v^3*o^2
          ierr=t_einsum("Hvvvo(a,b,e,i)+=Wabef(a,b,e,f)*T(f,i)",& 
               & Hvvvo_tensor,Wabef_tensor,comm_t%t1,scale=ONE)                  !Operation cost : v^4*o
          ierr=t_einsum("S(a,b,i,j)+= Hvvvo(a,b,e,i)*T(e,j)",& 
               & s2_tensor,Hvvvo_tensor,comm_t%t1,scale=MINUS_ONE_HALF)     ! Operation cost : v^3*o^2
          if (debug_local) print*, "norm s2 tensor after diagram 22 new ",tensor_norm2(s2_tensor)


!         Diagram 23 : s2_{ij}^{ab} +=  P\_(ab)\sum_{e} f_e^b t_{ij}^{ae} ! Use of the Fock matrix 
          ierr=t_einsum("S(a,b,i,j)+=F(b,e)*T(a,e,i,j)",s2_tensor,comm_t%fvv,comm_t%t2)
          ierr=t_einsum("S(a,b,i,j)+=F(a,e)*T(b,e,i,j)",s2_tensor,comm_t%fvv,comm_t%t2,scale=MINUS_ONE)
          if (debug_local) print*, "norm s2 tensor after diagram 23 ",tensor_norm2(s2_tensor)

!         Diagram 24 : s2_{ij}^{ab} += -P\_(ij)\sum_{k} f_j^k t_{ik}^{ab} ! Use of the Fock matrix
          ierr=t_einsum("S(a,b,i,j)+=F(m,j)*T(a,b,i,m)",s2_tensor,comm_t%foo,comm_t%t2,scale=MINUS_ONE)
          ierr=t_einsum("S(a,b,i,j)+=F(m,i)*T(a,b,j,m)",s2_tensor,comm_t%foo,comm_t%t2)
          if (debug_local) print*, "norm s2 tensor after diagram 24 ",tensor_norm2(s2_tensor)

          if (debug_local) print*, "Norm s2 before s1",tensor_norm2(s2_tensor)

!        done with cc2 t2 equation, now cc2 t1 equation, which is exactly the same as the ccsd t1 equation

          ! Call F intermediates - used in \Ss ai
          call get_hoo(comm_t,int_t%oovv)
          call get_hvv(comm_t,int_t%oovv)
          if (.not.CCD) then
            call get_hov(comm_t,int_t%oovv) 
            ierr=initialize_tensor_t(hov_aux_tensor)
            ierr=t_einsum("L(k,c)+=H(k,c)",hov_aux_tensor,comm_t%hov)                       ! copy common term
            ierr=t_einsum("H(k,c)+=F(k,c)",hov_aux_tensor,comm_t%fov,scale=MINUS_TWO)       ! subtract fock matrix twice
            if (debug_local) print *, "In .not.CCD check for hov"
          end if

!         calculate s1 tensor
          ierr=initialize_tensor_t(s1_tensor)

!           Diagramme b1  :   \Ss ai= \f ai         
          ierr=t_einsum("S(a,i)+=F+(i,a)*K()",s1_tensor,comm_t%fov,one_tensor)
          if (debug_local) print*, "s1 - first step with \F ia ",tensor_norm2(s1_tensor)
!           
          ierr=t_einsum("S(a,i)+=V(a,k,c,d)*T(c,d,i,k)",s1_tensor,int_t%vovv,comm_t%tau,scale=ONE_HALF) !validated o(2)v(3)
          ierr=t_einsum("S(a,i)+=H(a,c)*T(c,i)",s1_tensor,comm_t%hvv,comm_t%t1)                     !validated o(1)v(2)
          ierr=t_einsum("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,comm_t%hoo,comm_t%t1)                     !validated o(2)v(1)
          ierr=t_einsum("S(a,i)+=H(k,c)*T(a,c,i,k)",s1_tensor,comm_t%hov,comm_t%t2)                     !validated o(2)v(2)
          ierr=initialize_tensor_t(hoo_aux_tensor)                                                                     !moretests
          ierr=t_einsum("L(k,i)+=H(k,c)*T(c,i)",hoo_aux_tensor,hov_aux_tensor,comm_t%t1)                !moretests o(2)v(1)
          ierr=t_einsum("S(a,i)+=H(k,i)*T(a,k)",s1_tensor,hoo_aux_tensor,comm_t%t1)                     !moretests o(2)v(1)
          ierr=t_einsum("S(a,i)+=V(k,l,i,c)*T(c,a,k,l)",s1_tensor,int_t%ooov,comm_t%tau,scale=ONE_HALF) !validated o(3)v(2)
          ierr=t_einsum("S(a,i)+=V(a,k,c,i)*T(c,k)",s1_tensor,int_t%vovo,comm_t%t1,scale=MINUS_ONE)     !validated o(2)v(2)

          if (debug_local) print*, "Norm s1 before Energy",tensor_norm2(s1_tensor)

          if (ierr.ne.ETI_SUCCESS) then
             print*," ccdriver: error in evaluation of amplitude equations, code=",ierr
             return
          end if

!         prepare for next iteration
          ierr=initialize_tensor_t(comm_t%t1)
          ierr=t_einsum("T1(a,i)+=S1(a,i)",comm_t%t1,s1_tensor)
          ierr=initialize_tensor_t(comm_t%t2)
          ierr=t_einsum("T2(a,b,i,j)+=S2(a,b,i,j)",comm_t%t2,s2_tensor)
          
          call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%t1,comm_t%t2)

          if (.not.CCD) then
            call diis (t2_tensor=comm_t%t2,t1_tensor=comm_t%t1,t_convergence=t_convergence)
          else
            call diis (t2_tensor=comm_t%t2,t_convergence=t_convergence)
          end if

          ierr=initialize_tensor_t(result_tensor)
          ierr=t_einsum("E()+=T(a,i)*F(i,a)",result_tensor,comm_t%t1,comm_t%fov) ! o(1)v(1)
          ccsd_energy=result_tens(1)

!         Update tau
          call get_tau(comm_t,CCD)

          if (debug_local) print*, "Check for tau from get_tau", tensor_norm2(comm_t%tau)

          ierr=initialize_tensor_t(result_tensor)
          ierr=t_einsum("E()+=T(a,b,i,j)*V(i,j,a,b)",result_tensor,comm_t%tau,int_t%oovv) ! o(2)v(2)
          ccsd_energy=ccsd_energy + result_tens(1) * ONE_QUARTER_R

!         call print_iteration(iccsd,t_target_precision,t_convergence,ccsd_energy)
          call print_iteration(iccsd,t_convergence,ccsd_energy)
          if (t_convergence .lt. t_target_precision) exit

         end do


! Clean tensor
         ierr=destroy_tensor_t(s1_tensor)
         ierr=destroy_tensor_t(s2_tensor)
         ierr=destroy_tensor_t(cint_ovoo_tensor)
         ierr=destroy_tensor_t(cint_vvvo_tensor)

         ierr=destroy_tensor_t(Hvvov_tensor)
         ierr=destroy_tensor_t(Hvvvo_tensor)
         ierr=destroy_tensor_t(Hovvv_tensor)
         ierr=destroy_tensor_t(Hvovv_tensor)
         ierr=destroy_tensor_t(Hooov_tensor)
         ierr=destroy_tensor_t(Hoovo_tensor)
         ierr=destroy_tensor_t(Hovoo_tensor)
         ierr=destroy_tensor_t(Hvooo_tensor)

         ierr=destroy_tensor_t(Wabef_tensor)
         ierr=destroy_tensor_t(Wmnij_tensor) 

         ierr=destroy_tensor_t(Wabej_tensor)
         ierr=destroy_tensor_t(Wmbij_tensor)

         if (.not.CCD) then
            ierr=destroy_tensor_t(hoo_aux_tensor)
            ierr=destroy_tensor_t(hov_aux_tensor)
         end if
         ierr=destroy_tensor_t(result_tensor)
         ierr=destroy_tensor_t(one_tensor)

        end subroutine

        subroutine solve_lambda_equations(comm_t,int_t,nocc,nvir,eps_occ,eps_vir, &
                                          CCD,t_target_precision,ncycles,print_level)
!        Routine for solving the lambda equations

         logical, intent(in) :: CCD
         integer, intent(in) :: nocc, nvir, print_level
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(eti_intg_tens_t), intent(inout) :: int_t
         type(eti_lambda_tens_t)              :: lambda_t
!        solution tensors
         type(tensor_t) :: s1_tensor, s2_tensor
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
         integer(INTD)      :: result_dims(1)
         integer(INTD) :: ierr
         type(eti_precision_t):: eti_precision

         eti_precision%storage_t = comm_t%t2%datatype

         write(*,*) ""
         write(*,'(A)') "*******************************************************************************"
         write(*,'(A)') "**************** Solving for Ground-State Lagrange Multipliers ****************"
         write(*,'(A)') "*******************************************************************************"
         write(*,*) ""

!***************************************
!        construct fixed intermediates *
!***************************************
         result_dims=1
         ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,result_dims(1:0),init_val=ONE)

         oo_dims = nocc
         ierr=construct_tensor_t(lambda_t%fbar_oo,eti_precision%storage_t,oo_dims,init_val=ZERO)
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         ierr=construct_tensor_t(lambda_t%fbar_ov,eti_precision%storage_t,ov_dims,init_val=ZERO)
         vv_dims = nvir
         ierr=construct_tensor_t(lambda_t%fbar_vv,eti_precision%storage_t,vv_dims,init_val=ZERO)
         oooo_dims = nocc
         ierr=construct_tensor_t(lambda_t%w_oooo,eti_precision%storage_t,oooo_dims,init_val=ZERO)
         ovoo_dims(1)   = nocc
         ovoo_dims(2)   = nvir
         ovoo_dims(3:4) = nocc
         ierr=construct_tensor_t(lambda_t%w_ovoo,eti_precision%storage_t,ovoo_dims,init_val=ZERO)
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=construct_tensor_t(lambda_t%w_ooov,eti_precision%storage_t,ooov_dims,init_val=ZERO)
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=construct_tensor_t(lambda_t%w_vovo,eti_precision%storage_t,vovo_dims,init_val=ZERO)
         ierr=construct_tensor_t(lambda_t%wbar_vovo,eti_precision%storage_t,vovo_dims,init_val=ZERO)
         vvvo_dims(1:3) = nvir
         vvvo_dims(4)   = nocc
         ierr=construct_tensor_t(lambda_t%w_vvvo,eti_precision%storage_t,vvvo_dims,init_val=ZERO)
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=construct_tensor_t(lambda_t%w_vovv,eti_precision%storage_t,vovv_dims,init_val=ZERO)
         vvvv_dims = nvir
         ierr=construct_tensor_t(lambda_t%w_vvvv,eti_precision%storage_t,vvvv_dims,init_val=ZERO)

!******************************************
!        obtain intermediates from CC(S)D *
!******************************************

         call get_tau(comm_t,CCD)
         call get_hvv(comm_t,int_t%oovv)
         call get_gvv(comm_t,int_t%vovv,CCD)
         call get_h(comm_t,int_t,CCD)
         call get_hoo(comm_t,int_t%oovv)
         call get_goo(comm_t,int_t%ooov,CCD)
         if (.not.CCD) call get_hov(comm_t,int_t%oovv)
         call get_a(comm_t,int_t,CCD)
         call print_date("Finished calculating common blocks")

!*********************************
!        get fixed intermediates *
!*********************************

         call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD)
         call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t)
         if (.not.CCD) call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov)
         call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD)
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,nocc,nvir,CCD)
         call get_w_oooo(lambda_t%w_oooo,comm_t%a_int)
         call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
         call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD)
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD)
         call print_date('Finished calculating fixed intermediates')

         if (print_level>0) then
             if (.not.CCD) print*, "T1amp      = ", tensor_norm2(comm_t%t1)
             print*, "T2amp      = ", tensor_norm2(comm_t%t2)*ONE_QUARTER_R
             ! print*, "Tau        = ", tensor_norm2(comm_t%tau)*ONE_QUARTER_R
              print*, "Fbar_oo    = ", tensor_norm2(lambda_t%fbar_oo)
             print*, "Fbar_ov    = ", tensor_norm2(lambda_t%fbar_ov)
             print*, "Fbar_vv    = ", tensor_norm2(lambda_t%fbar_vv)
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

         if (.not.CCD) ierr=construct_tensor_t(s1_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=construct_tensor_t(s2_tensor,eti_precision%storage_t,oovv_dims,init_val=ZERO)

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=construct_tensor_t(lambda_t%goo,eti_precision%storage_t,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(lambda_t%gvv,eti_precision%storage_t,vv_dims,init_val=ZERO)

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=1)
         else
           call diis (t2_tensor=comm_t%l2,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,comm_t%l2,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,comm_t%l2,comm_t%t2)

!         debug
          if (print_level.gt.8) then
             if (.not.CCD) print*, "L1  = ", tensor_norm2(comm_t%l1)
             print*, "L2  = ", tensor_norm2(comm_t%l2)*ONE_QUARTER_R
             print*, "goo = ", tensor_norm2(lambda_t%goo)
             print*, "gvv = ", tensor_norm2(lambda_t%gvv)
          end if

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=initialize_tensor_t(s1_tensor)

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=t_einsum("S(i,a)+=F(i,a)",s1_tensor,lambda_t%fbar_ov)

!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,comm_t%l1,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,comm_t%l1,lambda_t%fbar_oo,scale=MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=t_einsum("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,comm_t%l1,lambda_t%w_vovo,scale=MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,comm_t%l2,lambda_t%w_ovoo,scale=ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,comm_t%l2,lambda_t%w_vvvo,scale=ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=t_einsum("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=initialize_tensor_t(s2_tensor)

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,b)",s2_tensor,int_t%oovv)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |
!----------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,comm_t%l2,lambda_t%fbar_vv)
          ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,comm_t%l2,lambda_t%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |
!----------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,comm_t%l2,lambda_t%fbar_oo)
          ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,comm_t%l2,lambda_t%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |
!------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,comm_t%l2,lambda_t%w_oooo,scale=ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |
!----------------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,comm_t%l2,lambda_t%w_vovo)
          ierr=t_einsum("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,comm_t%l2,lambda_t%w_vovo)
          ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,comm_t%l2,lambda_t%w_vovo)
          ierr=t_einsum("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,comm_t%l2,lambda_t%w_vovo)

!-------------------------------------------------------
!         term 6: S(ij,ab) += P(a,b) V(ij,ae) * G(e,b) |
!-------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,e)*G(e,b)",s2_tensor,int_t%oovv,lambda_t%gvv)
          ierr=t_einsum("S(i,j,a,b)+=V(i,j,e,b)*G(e,a)",s2_tensor,int_t%oovv,lambda_t%gvv)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(,j,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=t_einsum("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,comm_t%l1,lambda_t%w_ooov)
            ierr=t_einsum("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,comm_t%l1,lambda_t%w_ooov)
          end if

!-------------------------------------------------------
!         term 8: S(ij,ab) -= P(i,j) V(im,ab) * G(j,m) |
!-------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=V(i,m,b,a)*G(j,m)",s2_tensor,int_t%oovv,lambda_t%goo)
          ierr=t_einsum("S(i,j,a,b)+=V(j,m,a,b)*G(i,m)",s2_tensor,int_t%oovv,lambda_t%goo)

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=t_einsum("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,comm_t%l1,lambda_t%w_vovv)
            ierr=t_einsum("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,comm_t%l1,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=t_einsum("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,comm_t%l1,lambda_t%fbar_ov)
            ierr=t_einsum("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,comm_t%l1,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=t_einsum("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,comm_t%l1,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=t_einsum("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,comm_t%l1,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |
!-----------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,comm_t%l2,lambda_t%w_vvvv,scale=ONE_HALF)

          if (ierr.ne.ETI_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of lambda equations, code=",ierr
             return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=initialize_tensor_t(comm_t%l1)
            ierr=t_einsum("L1(i,a)+=S1(i,a)",comm_t%l1,s1_tensor)
          end if
          ierr=initialize_tensor_t(comm_t%l2)
          ierr=t_einsum("L2(i,j,a,b)+=S2(i,j,a,b)",comm_t%l2,s2_tensor)
          if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%l1,comm_t%l2)
            call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,t_convergence=t_convergence)
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%l2)
            call diis (t2_tensor=comm_t%l2,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", tensor_norm2(comm_t%l1)
         write(*,'(A,ES22.16)') "  Final L2amp = ", tensor_norm2(comm_t%l2)*ONE_QUARTER_R
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=destroy_tensor_t(s1_tensor)
           call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=2)
         else
           call diis (t2_tensor=comm_t%l2,mode=2)
         end if
         ierr=destroy_tensor_t(s2_tensor)
         ierr=destroy_tensor_t(one_tensor)
         ierr=destroy_tensor_t(lambda_t%fbar_oo)
         ierr=destroy_tensor_t(lambda_t%fbar_ov)
         ierr=destroy_tensor_t(lambda_t%fbar_vv)
         ierr=destroy_tensor_t(lambda_t%w_oooo)
         ierr=destroy_tensor_t(lambda_t%w_ooov)
         ierr=destroy_tensor_t(lambda_t%w_ovoo)
         ierr=destroy_tensor_t(lambda_t%w_vovo)
         ierr=destroy_tensor_t(lambda_t%w_vvvo)
         ierr=destroy_tensor_t(lambda_t%w_vovv)
         ierr=destroy_tensor_t(lambda_t%w_vvvv)
         ierr=destroy_tensor_t(lambda_t%wbar_vovo)
         ierr=destroy_tensor_t(lambda_t%goo)
         ierr=destroy_tensor_t(lambda_t%gvv)

        end subroutine solve_lambda_equations


        subroutine solve_CC2_lambda_equations(comm_t,int_t,nocc,nvir, &
                                     eps_occ,eps_vir,t_target_precision,ncycles,print_level)

!        Routine for solving the lambda equations
         use intermediates_CC2

         integer, intent(in) :: nocc, nvir
!        fixed 1- and 2-body tensors (fock matrix elements and two-electron integrals)
         type(eti_comm_tens_t), intent(inout) :: comm_t
         type(eti_comm_tens_t) :: comm_t_CC2
         type(eti_intg_tens_t), intent(inout) :: int_t
         type(eti_lambda_tens_t)              :: lambda_t,lambda_t_CC2
!        solution tensors
         type(tensor_t) :: s1_tensor, s2_tensor
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
         integer(INTD) :: result_dims(1)
         type(eti_precision_t):: eti_precision

         eti_precision%storage_t = comm_t%t2%datatype


         write(*,*) ""
         write(*,'(A)') "*******************************************************************************"
         write(*,'(A)') "**************** Solving for Ground-State Lagrange Multipliers ****************"
         write(*,'(A)') "*******************************************************************************"
         write(*,*) ""

!***************************************
!        construct fixed intermediates *
!***************************************
         result_dims=1
         ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,result_dims(1:0),init_val=ONE)

         oo_dims = nocc
         ierr=construct_tensor_t(lambda_t%fbar_oo,eti_precision%storage_t,oo_dims,init_val=ZERO)
         ov_dims(1) = nocc
         ov_dims(2) = nvir
         ierr=construct_tensor_t(lambda_t%fbar_ov,eti_precision%storage_t,ov_dims,init_val=ZERO)
         vv_dims = nvir
         ierr=construct_tensor_t(lambda_t%fbar_vv,eti_precision%storage_t,vv_dims,init_val=ZERO)
         oooo_dims = nocc
         ierr=construct_tensor_t(lambda_t%w_oooo,eti_precision%storage_t,oooo_dims,init_val=ZERO)
         ovoo_dims(1)   = nocc
         ovoo_dims(2)   = nvir
         ovoo_dims(3:4) = nocc
         ierr=construct_tensor_t(lambda_t%w_ovoo,eti_precision%storage_t,ovoo_dims,init_val=ZERO)
         ooov_dims(1:3) = nocc
         ooov_dims(4)   = nvir
         ierr=construct_tensor_t(lambda_t%w_ooov,eti_precision%storage_t,ooov_dims,init_val=ZERO)
         vovo_dims(1) = nvir
         vovo_dims(2) = nocc
         vovo_dims(3) = nvir
         vovo_dims(4) = nocc
         ierr=construct_tensor_t(lambda_t%w_vovo,eti_precision%storage_t,vovo_dims,init_val=ZERO)
         ierr=construct_tensor_t(lambda_t%wbar_vovo,eti_precision%storage_t,vovo_dims,init_val=ZERO)
         vvvo_dims(1:3) = nvir
         vvvo_dims(4)   = nocc
         ierr=construct_tensor_t(lambda_t%w_vvvo,eti_precision%storage_t,vvvo_dims,init_val=ZERO)
         vovv_dims(1)   = nvir
         vovv_dims(2)   = nocc
         vovv_dims(3:4) = nvir
         ierr=construct_tensor_t(lambda_t%w_vovv,eti_precision%storage_t,vovv_dims,init_val=ZERO)
         vvvv_dims = nvir
         ierr=construct_tensor_t(lambda_t%w_vvvv,eti_precision%storage_t,vvvv_dims,init_val=ZERO)

!***
! List of NEW INTERMEDIATES for the CC2 Method
!***  
         ierr=construct_tensor_t(lambda_t_CC2%fbar_oo,eti_precision%storage_t,oo_dims,init_val=ZERO)   !cUSefull and different
         ierr=construct_tensor_t(lambda_t_CC2%fbar_ov,eti_precision%storage_t,ov_dims,init_val=ZERO)  ! Usefull but same as lambda1 and no difference with CCSD
         ierr=construct_tensor_t(lambda_t_CC2%fbar_vv,eti_precision%storage_t,vv_dims,init_val=ZERO)   ! Usefull and different
         ierr=construct_tensor_t(lambda_t_CC2%w_oooo,eti_precision%storage_t,oooo_dims,init_val=ZERO)   ! Usefull and different
!!!!!!  Inportant here :  The heart of w_vovo intermediate : h_int is the same as in 
!!!!!!          CCSD and for \Lambda_1 either :
         ierr=construct_tensor_t(lambda_t_CC2%w_vovo,eti_precision%storage_t,vovo_dims,init_val=ZERO) ! Usefull in CC2 - Same part as h_int
         ierr=construct_tensor_t(lambda_t_CC2%w_vvvv,eti_precision%storage_t,vvvv_dims,init_val=ZERO)   ! Usefull and different

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
         if (.not.CCD) call get_hov(comm_t,int_t%oovv)
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
         ierr=construct_tensor_t(comm_t_CC2%t1,eti_precision%storage_t,vo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t_CC2%t2,eti_precision%storage_t,vvoo_dims,init_val=ZERO)
         ierr=construct_tensor_t(comm_t_CC2%tau,c8,vvoo_dims,init_val=zero)
         ierr=construct_tensor_t(comm_t_CC2%h_int,c8,vovo_dims,init_val=zero)

         ierr=t_einsum("H(a,i)+=T(a,i)",comm_t_CC2%t1,comm_t%t1,scale=ONE)
         ierr=t_einsum("H(a,b,i,j)+=T(a,b,i,j)",comm_t_CC2%t2,comm_t%t2,scale=ONE)
         ierr=t_einsum("H(a,i)+=F(a,i)",comm_t_CC2%fov,comm_t%fov,scale=ONE)
         ierr=t_einsum("H(a,b)+=F(a,b)",comm_t_CC2%fvv,comm_t%fvv,scale=ONE)
         ierr=t_einsum("H(i,j)+=F(i,j)",comm_t_CC2%foo,comm_t%foo,scale=ONE)
         ierr=t_einsum("W(b,m,e,j)+=H(b,m,e,j)",comm_t_CC2%h_int,comm_t%h_int,scale=ONE)

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

         if (.not.CCD) ierr=construct_tensor_t(s1_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
         oovv_dims(1:2) = nocc
         oovv_dims(3:4) = nvir
         ierr=construct_tensor_t(s2_tensor,eti_precision%storage_t,oovv_dims,init_val=ZERO)

!********************************************
!        construct mixed type intermediates *
!********************************************

         ierr=construct_tensor_t(lambda_t%goo,eti_precision%storage_t,oo_dims,init_val=ZERO)
         ierr=construct_tensor_t(lambda_t%gvv,eti_precision%storage_t,vv_dims,init_val=ZERO)

!*****************************************************************
!        re-initialize diis, so it can be used in lambda routine *
!*****************************************************************

         if (.not.CCD) then
           call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=1)
         else
           call diis (t2_tensor=comm_t%l2,mode=1)
         end if

!*****************
!        iterate *
!*****************

         do ilambda = 1, ncycles

!***********************************
!         form mixed intermediates *
!***********************************

          call get_gvv_lambda(lambda_t%gvv,comm_t%l2,comm_t%t2)
          call get_goo_lambda(lambda_t%goo,comm_t%l2,comm_t%t2)

!******************************
!         calculate s1 tensor *
!******************************

          if (.not.CCD) then
            ierr=initialize_tensor_t(s1_tensor)

!----------------------------------------
!           term 1: S(i,a) += Fbar(i,a) |
!----------------------------------------

            ierr=t_einsum("S(i,a)+=F(i,a)",s1_tensor,lambda_t%fbar_ov)

!-------------------------------------------------
!           term 2: S(i,a) += L(i,e) * Fbar(e,a) |
!-------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(i,e)*F(e,a)",s1_tensor,comm_t%l1,lambda_t%fbar_vv)

!-------------------------------------------------
!           term 3: S(i,a) -= L(m,a) * Fbar(i,m) |
!-------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(m,a)*F(i,m)",s1_tensor,comm_t%l1,lambda_t%fbar_oo,scale=MINUS_ONE)

!------------------------------------------------
!           term 6: S(i,a) -= G(n,m) * W(mi,na) |
!------------------------------------------------

            ierr=t_einsum("S(i,a)+=G(n,m)*W(i,m,n,a)",s1_tensor,lambda_t%goo,lambda_t%w_ooov)

!------------------------------------------------
!           term 7: S(i,a) += L(m,e) * W(ie,am) |
!------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(m,e)*W(e,i,a,m)",s1_tensor,comm_t%l1,lambda_t%w_vovo,scale=MINUS_ONE)

!---------------------------------------------------------
!           term 8: S(i,a) -= 1/2 * L(mn,ae)  * W(ie,mn) |
!---------------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(m,n,a,e)*W(i,e,n,m)",s1_tensor,comm_t%l2,lambda_t%w_ovoo,scale=ONE_HALF)

!------------------------------------------------------
!           term 4: S(i,a) += 1/2 L(im,ef) * W(ef,am) |
!------------------------------------------------------

            ierr=t_einsum("S(i,a)+=L(i,m,e,f)*W(e,f,a,m)",s1_tensor,comm_t%l2,lambda_t%w_vvvo,scale=ONE_HALF)

!------------------------------------------------
!           term 5: S(i,a) -= G(f,e) * W(ei,fa) |
!------------------------------------------------

            ierr=t_einsum("S(i,a)+=G(f,e)*W(e,i,a,f)",s1_tensor,lambda_t%gvv,lambda_t%w_vovv)
          end if

!******************************
!         calculate s2 tensor *
!******************************

          ierr=initialize_tensor_t(s2_tensor)

!---------------------------------------
!         term 1: S(ij,ab) += V(ij,ab) |
!---------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=V(i,j,a,b)",s2_tensor,int_t%oovv)

!----------------------------------------------------------
!         term 2: S(ij,ab) += P(a,b) L(ij,ae) * Fbar(e,b) |      Fbar_vv - CC2
!----------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,j,a,e)*F(e,b)",s2_tensor,comm_t%l2,lambda_t_CC2%fbar_vv)
          ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,b)*F(e,a)",s2_tensor,comm_t%l2,lambda_t_CC2%fbar_vv)

!----------------------------------------------------------
!         term 3: S(ij,ab) -= P(i,j) L(im,ab) * Fbar(j,m) |      Fbar_oo - CC2
!----------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,m,b,a)*F(j,m)",s2_tensor,comm_t%l2,lambda_t_CC2%fbar_oo)
          ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,b)*F(i,m)",s2_tensor,comm_t%l2,lambda_t_CC2%fbar_oo)

!------------------------------------------------------
!         term 4: S(ij,ab) += 1/2 L(mn,ab) * W(ij,mn) |          W_oooo - CC2
!------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(m,n,a,b)*W(i,j,m,n)",s2_tensor,comm_t%l2,lambda_t_CC2%w_oooo,scale=ONE_HALF)

!----------------------------------------------------------------
!         term 5: S(ij,ab) += P(i,j) P(a,b) L(im,ae) * W(je,bm) |       W vovo - CC2
!----------------------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,m,e,a)*W(e,j,b,m)",s2_tensor,comm_t%l2,lambda_t_CC2%w_vovo)
          ierr=t_einsum("S(i,j,a,b)+=L(i,m,b,e)*W(e,j,a,m)",s2_tensor,comm_t%l2,lambda_t_CC2%w_vovo)
          ierr=t_einsum("S(i,j,a,b)+=L(j,m,a,e)*W(e,i,b,m)",s2_tensor,comm_t%l2,lambda_t_CC2%w_vovo)
          ierr=t_einsum("S(i,j,a,b)+=L(j,m,e,b)*W(e,i,a,m)",s2_tensor,comm_t%l2,lambda_t_CC2%w_vovo)

!-------------------------------------------------------
!         term 7: S(ij,ab) -= P(a,b) L(m,a) * W(ij,mb) |
!-------------------------------------------------------

          if (.not.CCD) then
            ierr=t_einsum("S(i,j,a,b)+=L(m,a)*W(j,i,m,b)",s2_tensor,comm_t%l1,lambda_t%w_ooov)
            ierr=t_einsum("S(i,j,a,b)+=L(m,b)*W(i,j,m,a)",s2_tensor,comm_t%l1,lambda_t%w_ooov)
          end if

!----------------------------------------------------------
!         term 9+10: S(ij,ab) += P(i,j) L(i,e) * W(ej,ab) |
!----------------------------------------------------------

          if (.not.CCD) then
            ierr=t_einsum("S(i,j,a,b)+=L(i,e)*W(e,j,a,b)",s2_tensor,comm_t%l1,lambda_t%w_vovv)
            ierr=t_einsum("S(i,j,a,b)+=L(j,e)*W(e,i,b,a)",s2_tensor,comm_t%l1,lambda_t%w_vovv)
          end if

!----------------------------------------------------------------
!         term 11: S(ij,ab) += P(i,j) P(a,b) L(i,a) * Fbar(j,b) |
!----------------------------------------------------------------

          if (.not.CCD) then
            ierr=t_einsum("S(i,j,a,b)+=L(i,a)*F(j,b)",s2_tensor,comm_t%l1,lambda_t%fbar_ov)
            ierr=t_einsum("S(i,j,a,b)+=L(i,b)*F(j,a)",s2_tensor,comm_t%l1,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=t_einsum("S(i,j,a,b)+=L(j,a)*F(i,b)",s2_tensor,comm_t%l1,lambda_t%fbar_ov,scale=MINUS_ONE)
            ierr=t_einsum("S(i,j,a,b)+=L(j,b)*F(i,a)",s2_tensor,comm_t%l1,lambda_t%fbar_ov)
          end if

!-----------------------------------------------------
!         term 12: S(ij,ab) += 1/2 L(ij,ef) W(ef,ab) |      W vvvv - CC2
!-----------------------------------------------------

          ierr=t_einsum("S(i,j,a,b)+=L(i,j,e,f)*W(e,f,a,b)",s2_tensor,comm_t%l2,lambda_t_CC2%w_vvvv,scale=ONE_HALF)

          if (ierr.ne.ETI_SUCCESS) then
             write(*,*) " ccdriver: error in evaluation of lambda equations, code=",ierr
             return
          end if

!*************************************
!         prepare for next iteration *
!*************************************

          if (.not.CCD) then
            ierr=initialize_tensor_t(comm_t%l1)
            ierr=t_einsum("L1(i,a)+=S1(i,a)",comm_t%l1,s1_tensor)
          end if
          ierr=initialize_tensor_t(comm_t%l2)
          ierr=t_einsum("L2(i,j,a,b)+=S2(i,j,a,b)",comm_t%l2,s2_tensor)
          if (.not.CCD) then
            call scale_with_denominators (eps_occ,eps_vir,nocc,comm_t%l1,comm_t%l2)
            call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,t_convergence=t_convergence)
          else
            call scale_with_denominators (eps_occ,eps_vir,nocc,t2_tensor=comm_t%l2)
            call diis (t2_tensor=comm_t%l2,t_convergence=t_convergence)
          end if

          call print_iteration(ilambda,t_convergence)
          if (t_convergence .lt. t_target_precision) exit

         end do

         write(*,*) "------------------------"
         write(*,*) ""

         if (t_convergence .gt. t_target_precision) write(*,*) "WARNING: Non-converged amplitudes!"
         if (.not.CCD) write(*,'(A,ES22.16)') "  Final L1amp = ", tensor_norm2(comm_t%l1)
         write(*,'(A,ES22.16)') "  Final L2amp = ", tensor_norm2(comm_t%l2)*ONE_QUARTER_R
         write(*,*) ""

!*****************
!        cleanup *
!*****************

         if (.not.CCD) then
           ierr=destroy_tensor_t(s1_tensor)
           call diis (t2_tensor=comm_t%l2,t1_tensor=comm_t%l1,mode=2)
         else
           call diis (t2_tensor=comm_t%l2,mode=2)
         end if
         ierr=destroy_tensor_t(s2_tensor)
         ierr=destroy_tensor_t(lambda_t%fbar_oo)
         ierr=destroy_tensor_t(lambda_t%fbar_ov)
         ierr=destroy_tensor_t(lambda_t%fbar_vv)
         ierr=destroy_tensor_t(lambda_t%w_oooo)
         ierr=destroy_tensor_t(lambda_t%w_ooov)
         ierr=destroy_tensor_t(lambda_t%w_ovoo)
         ierr=destroy_tensor_t(lambda_t%w_vovo)
         ierr=destroy_tensor_t(lambda_t%w_vvvo)
         ierr=destroy_tensor_t(lambda_t%w_vovv)
         ierr=destroy_tensor_t(lambda_t%w_vvvv)
         ierr=destroy_tensor_t(lambda_t%wbar_vovo)
         ierr=destroy_tensor_t(lambda_t%goo)
         ierr=destroy_tensor_t(lambda_t%gvv)

        end subroutine solve_CC2_lambda_equations

end module 
