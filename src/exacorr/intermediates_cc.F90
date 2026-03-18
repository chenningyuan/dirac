module intermediates_CC

!Module to make intermediates needed in CC and Lambda equations
     use exacorr_scaling_parameters

     implicit none
     private

     interface get_tau
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_tau_exatns
#endif
        module procedure get_tau_eti
     end interface get_tau

     interface get_hoo
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_hoo_exatns
#endif
        module procedure get_hoo_eti
     end interface get_hoo

     interface get_hov
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_hov_exatns
#endif
        module procedure get_hov_eti
     end interface get_hov

     interface get_hvv
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_hvv_exatns
#endif
        module procedure get_hvv_eti
     end interface get_hvv

     interface get_goo
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_goo_exatns
#endif
        module procedure get_goo_eti
     end interface get_goo

     interface get_gvv
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_gvv_exatns
#endif
        module procedure get_gvv_eti
     end interface get_gvv

     interface get_a
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_a_exatns
#endif
        module procedure get_a_eti
     end interface get_a

     interface get_h
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_h_exatns
#endif
        module procedure get_h_eti
     end interface get_h

!    CCSD intermediates
     public get_tau
     public get_hoo
     public get_hov
     public get_hvv
     public get_goo
     public get_gvv
     public get_a
     public get_h

     contains

       subroutine get_tau_exatns(comm_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       make tau ! validated
        ierr=exatns_tensor_init(comm_t%tau,'ZERO')
        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,b,i,j)*K()",comm_t%tau,comm_t%t2,one_tensor)
        if (.not.CCD) then
            ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1) ! o(2)v(2)
            ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,MINUS_ONE) ! o(2)v(2)
        end if

       end subroutine get_tau_exatns


       subroutine get_hoo_exatns(comm_t,oovv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%hoo,'ZERO')
        ierr=exatns_tensor_contract("H(k,i)+=F(k,i)*K()",comm_t%hoo,comm_t%foo,one_tensor, MINUS_ONE)
        ierr=exatns_tensor_contract("H(k,i)+=V(k,l,c,d)*T(c,d,l,i)",comm_t%hoo,oovv_tensor,comm_t%tau,ONE_HALF) ! o(3)v(2)

       end subroutine get_hoo_exatns

       subroutine get_hov_exatns(comm_t,oovv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%hov,'ZERO')
        ierr=exatns_tensor_contract("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,oovv_tensor,comm_t%t1) ! o(2)v(2)
        ierr=exatns_tensor_contract("H(k,c)+=F(k,c)*K()",comm_t%hov,comm_t%fov,one_tensor)           ! add fock matrix

       end subroutine get_hov_exatns

       subroutine get_hvv_exatns(comm_t,oovv_tensor,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%hvv,'ZERO')
        ierr=exatns_tensor_contract("H(a,c)+=F(a,c)*K()",comm_t%hvv,comm_t%fvv,one_tensor)
        ierr=exatns_tensor_contract("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%hvv,oovv_tensor,comm_t%tau,ONE_HALF) ! o(2)v(3)

       end subroutine get_hvv_exatns

       subroutine get_goo_exatns(comm_t,ooov_tensor,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: ooov_tensor
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%goo,'ZERO')
        ierr=exatns_tensor_contract("G(k,i)+=H(k,i)*K()",comm_t%goo,comm_t%hoo,one_tensor)
        if (.not.CCD) then 
            ierr=exatns_tensor_contract("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,ooov_tensor,comm_t%t1,MINUS_ONE) ! o(3)v(1)
            ierr=exatns_tensor_contract("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1,MINUS_ONE) ! o(2)v(1)
        end if

       end subroutine get_goo_exatns

       subroutine get_gvv_exatns(comm_t,vovv_tensor,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: vovv_tensor
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%gvv,'ZERO')
        ierr=exatns_tensor_contract("G(a,c)+=H(a,c)*K()",comm_t%gvv,comm_t%hvv,one_tensor)
        if (.not.CCD) then
            ierr=exatns_tensor_contract("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,vovv_tensor,comm_t%t1) ! o(1)v(3)
            ierr=exatns_tensor_contract("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,MINUS_ONE) ! o(1)v(2)
        end if

       end subroutine get_gvv_exatns

       subroutine get_a_exatns(comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t,exatns_intg_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(exatns_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%a_int,'ZERO')
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,j)*K()",comm_t%a_int,int_t%oooo,one_tensor)
        if (.not.CCD) then 
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
            ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1,MINUS_ONE) ! o(4)v(1)
        end if
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%a_int,int_t%oovv,comm_t%tau,ONE_HALF) ! o(4)v(2)

       end subroutine get_a_exatns

       subroutine get_h_exatns(comm_t,int_t,CCD,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t,exatns_intg_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(exatns_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

        ierr=exatns_tensor_init(comm_t%h_int,'ZERO')
        if (.not.CCD) then
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1) ! o(2)v(3)
            ierr=exatns_tensor_contract("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,MINUS_ONE) ! o(3)v(2)
        end if
        ierr=exatns_tensor_contract("H(a,k,c,i)+=V(a,k,c,i)*K()",comm_t%h_int,int_t%vovo,one_tensor)

       end subroutine get_h_exatns

!!!!!!!!!!!!!!!!!!!!!!!

       subroutine get_tau_eti(comm_t,CCD)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        logical, intent(in) :: CCD
        integer :: ierr

!       tau ! validated
        ierr=initialize_tensor_t(comm_t%tau)
        ierr=t_einsum("T(a,b,i,j)+=R(a,b,i,j)",comm_t%tau,comm_t%t2)
        if (.not.CCD) then
            ierr=t_einsum("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1)                  ! o(2)v(2)
            ierr=t_einsum("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,scale=MINUS_ONE)  ! o(2)v(2)
        end if

       end subroutine get_tau_eti

       subroutine get_hoo_eti(comm_t,oovv_tensor)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(tensor_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       oo ! validated
        ierr=initialize_tensor_t(comm_t%hoo)
        ierr=t_einsum("H(k,i)+=F(k,i)",comm_t%hoo,comm_t%foo,scale=MINUS_ONE)
        ierr=t_einsum("H(k,i)+=V(k,l,c,d)*T(c,d,l,i)",comm_t%hoo,oovv_tensor,comm_t%tau,scale=ONE_HALF) ! o(3)v(2)

       end subroutine get_hoo_eti

       subroutine get_hov_eti(comm_t,oovv_tensor)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(tensor_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       ov ! validated
        ierr=initialize_tensor_t(comm_t%hov)
        ierr=t_einsum("H(k,c)+=V(k,l,c,d)*T(d,l)",comm_t%hov,oovv_tensor,comm_t%t1) ! o(2)v(2)
        ierr=t_einsum("H(k,c)+=F(k,c)",comm_t%hov,comm_t%fov)

       end subroutine get_hov_eti

       subroutine get_hvv_eti(comm_t,oovv_tensor)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(tensor_t),intent(inout) :: oovv_tensor
        integer :: ierr

!       vv ! validated
        ierr=initialize_tensor_t(comm_t%hvv)
        ierr=t_einsum("H(a,c)+=F(a,c)",comm_t%hvv,comm_t%fvv)
        ierr=t_einsum("H(a,c)+=V(k,l,c,d)*T(d,a,k,l)",comm_t%hvv,oovv_tensor,comm_t%tau,scale=ONE_HALF) ! o(2)v(3)

       end subroutine get_hvv_eti

       subroutine get_goo_eti(comm_t,ooov_tensor,CCD)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(tensor_t),intent(inout) :: ooov_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       oo ! validated
        ierr=initialize_tensor_t(comm_t%goo)
        ierr=t_einsum("G(k,i)+=H(k,i)",comm_t%goo,comm_t%hoo)
        if (.not.CCD) then
            ierr=t_einsum("G(k,i)+=V(k,l,i,c)*T(c,l)",comm_t%goo,ooov_tensor,comm_t%t1,scale=MINUS_ONE)     ! o(3)v(1)
            ierr=t_einsum("G(k,i)+=F(k,c)*T(c,i)",comm_t%goo,comm_t%fov,comm_t%t1,scale=MINUS_ONE)                         ! o(2)v(1)
        end if

       end subroutine get_goo_eti

       subroutine get_gvv_eti(comm_t,vovv_tensor,CCD)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(tensor_t),intent(inout) :: vovv_tensor
        logical, intent(in) :: CCD
        integer :: ierr

!       vv ! validated
        ierr=initialize_tensor_t(comm_t%gvv)
        ierr=t_einsum("G(a,c)+=H(a,c)",comm_t%gvv,comm_t%hvv)
        if (.not.CCD) then
            ierr=t_einsum("G(a,c)+=V(a,k,c,d)*T(d,k)",comm_t%gvv,vovv_tensor,comm_t%t1)             ! o(1)v(3)
            ierr=t_einsum("G(a,c)+=F(k,c)*T(a,k)",comm_t%gvv,comm_t%fov,comm_t%t1,scale=MINUS_ONE) ! o(1)v(2)
        end if

       end subroutine get_gvv_eti

       subroutine get_a_eti(comm_t,int_t,CCD)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(eti_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

!       oooo ! validated
        ierr=initialize_tensor_t(comm_t%a_int)
        ierr=t_einsum("A(k,l,i,j)+=V(k,l,i,j)",comm_t%a_int,int_t%oooo)
        if (.not.CCD) then
            ierr=t_einsum("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",comm_t%a_int,int_t%ooov,comm_t%t1) ! o(4)v(1)
            ierr=t_einsum("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",comm_t%a_int,int_t%ooov,comm_t%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        end if
        ierr=t_einsum("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",comm_t%a_int,int_t%oovv,comm_t%tau, &
                                    scale=ONE_HALF)

       end subroutine get_a_eti

       subroutine get_h_eti(comm_t,int_t,CCD)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t
        type(eti_intg_tens_t),intent(inout) :: int_t
        logical, intent(in) :: CCD
        integer :: ierr

!       vovo ! more tests
        ierr=initialize_tensor_t(comm_t%h_int)
        if (.not.CCD) then
            ierr=t_einsum("H(a,k,c,i)+=V(a,k,c,d)*T(d,i)",comm_t%h_int,int_t%vovv,comm_t%t1)                 ! o(2)v(3)
            ierr=t_einsum("H(a,k,c,i)+=V(k,l,i,c)*T(a,l)",comm_t%h_int,int_t%ooov,comm_t%t1,scale=MINUS_ONE) ! o(3)v(2)
        end if
        ierr=t_einsum("H(a,k,c,i)+=V(a,k,c,i)",comm_t%h_int,int_t%vovo)

       end subroutine get_h_eti

end module intermediates_CC
