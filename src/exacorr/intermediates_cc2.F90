module intermediates_cc2

!Module to make intermediates needed in Lambda equations

     implicit none
     private

     complex(8), parameter :: ZERO=(0.D0,0.D0), MINUS_ONE=(-1.D0,0.D0), ONE=(1.D0,0.D0), &
                              ONE_QUARTER_C=(0.25D0,0.D0), ONE_HALF=(0.5D0,0.D0), &
                              MINUS_ONE_HALF=(-0.5D0,0.D0), &
                              MINUS_TWO=(-2.D0,0.D0), TWO=(2.D0,0.D0)
    real(8), parameter     :: ONE_QUARTER=0.25D0

!! Start CC2
     interface get_tauCC2 !L
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_tauCC2_exatns
#endif
        module procedure get_tauCC2_eti
     end interface get_tauCC2

     interface get_fbar_oo_CC2
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_fbar_oo_CC2_exatns
#endif
        module procedure get_fbar_oo_CC2_eti
     end interface get_fbar_oo_CC2

     interface get_fbar_vv_CC2
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_fbar_vv_CC2_exatns
#endif
        module procedure get_fbar_vv_CC2_eti
     end interface get_fbar_vv_CC2

!fbar_ov_CC2 Not usefulll same as CCSD but here to check equations
     interface get_fbar_ov_CC2
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_fbar_ov_CC2_exatns
#endif
        module procedure get_fbar_ov_CC2_eti
     end interface get_fbar_ov_CC2


     interface get_w_vvvv_CC2
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_w_vvvv_CC2_exatns
#endif
        module procedure get_w_vvvv_CC2_eti
     end interface get_w_vvvv_CC2

     interface get_w_oooo_CC2
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_w_oooo_CC2_exatns
#endif
        module procedure get_w_oooo_CC2_eti
     end interface get_w_oooo_CC2

     interface get_w_vovo_CC2
#if (defined (VAR_MPI) && !defined(EXA_TALSH_ONLY))
        module procedure get_w_vovo_CC2_exatns
#endif
        module procedure get_w_vovo_CC2_eti
     end interface get_w_vovo_CC2

!! End CC2

!    Lambda CC2 Intermediates
     public get_tauCC2  !L
     public get_fbar_oo_CC2
     public get_fbar_ov_CC2 ! Not usefull here just to check equations
     public get_fbar_vv_CC2
     public get_w_vovo_CC2
     public get_w_oooo_CC2
     public get_w_vvvv_CC2


     contains

       subroutine get_tauCC2_exatns(comm_t,one_tensor) !L

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(exatns_comm_tens_t),intent(inout) :: comm_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       make tau ! validated
!L        ierr=exatns_tensor_init(comm_t%tau)
!!        ierr=exatns_tensor_contract("T(a,b,i,j)+=R(a,b,i,j)*K()",comm_t%tau,comm_t%t2,one_tensor)
        ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t%tau,comm_t%t1,comm_t%t1) ! o(2)v(2)
        ierr=exatns_tensor_contract("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t%tau,comm_t%t1,comm_t%t1,MINUS_ONE) ! o(2)v(2)
 
       end subroutine get_tauCC2_exatns

       subroutine get_fbar_vv_CC2_exatns(fbar_vv,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens_t,exatns_comm_tens_t

!        type(exatns_lambda_tens_t),intent(inout) :: lambda_t_CC2
        type(tens_rcrsv_t),intent(inout) :: fbar_vv

        type(exatns_comm_tens_t),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens_t),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       \Ft ea =\f ea -\SUMs \f ma \ts em - \SUMd mf \V mefa \ts fm - \undemi
!       \Vmnaf \taup efmn (\tau in CC2)  (Eq. A5)

        ierr=exatns_tensor_init(fbar_vv,'ZERO')
        ierr=exatns_tensor_contract("H(e,a)+=F(e,a)*K()",fbar_vv,comm_t_CC2%fvv,one_tensor)
        ierr=exatns_tensor_contract("H(e,a)+=F(m,a)*T(e,m)",fbar_vv,comm_t_CC2%fov,comm_t_CC2%t1,MINUS_ONE)
        ierr=exatns_tensor_contract("H(e,a)+=V(m,e,f,a)*T(f,m)",fbar_vv,int_t%vovv,comm_t_CC2%t1,ONE)
        ierr=exatns_tensor_contract("H(e,a)+=V(m,n,a,f)*T(e,f,m,n)",fbar_vv,int_t%oovv,comm_t_CC2%tau,MINUS_ONE_HALF)

       end subroutine get_fbar_vv_CC2_exatns

       subroutine get_fbar_oo_CC2_exatns(fbar_oo,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens_t,exatns_comm_tens_t

        type(tens_rcrsv_t),intent(inout) :: fbar_oo
        type(exatns_comm_tens_t),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens_t),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!      \Ft im =\f im +\SUMs e \f ie \ts em +\SUMd en \V inme \ts en +\ undemi
!      \SUMt mef \V inef \taup efmn  (\tau in CC2)  (Eq. A4)

        ierr=exatns_tensor_init(fbar_oo,'ZERO')
        ierr=exatns_tensor_contract("H(i,m)+=F(i,m)*K()",fbar_oo,comm_t_CC2%foo,one_tensor)
        ierr=exatns_tensor_contract("H(i,m)+=F(i,e)*T(e,m)",fbar_oo,comm_t_CC2%fov,comm_t_CC2%t1,ONE)
        ierr=exatns_tensor_contract("H(i,m)+=V(i,n,m,e)*T(e,n)",fbar_oo,int_t%ooov,comm_t_CC2%t1,ONE)
        ierr=exatns_tensor_contract("H(i,m)+=V(i,n,e,f)*T(e,f,m,n)",fbar_oo,int_t%oovv,comm_t_CC2%tau,ONE_HALF)

       end subroutine get_fbar_oo_CC2_exatns

       subroutine get_fbar_ov_CC2_exatns(lambda_t_CC2,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens_t,exatns_comm_tens_t,exatns_lambda_tens_t

        type(exatns_lambda_tens_t),intent(inout) :: lambda_t_CC2
        type(exatns_comm_tens_t),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens_t),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!      \Ft me =\f me +\SUMd nf \V mnef \ts fn (Eq. A6)
        ierr=exatns_tensor_init(lambda_t_CC2%fbar_ov,'ZERO')
        ierr=exatns_tensor_contract("H(m,e)+=F(m,e)*K()",lambda_t_CC2%fbar_ov,comm_t_CC2%fov,one_tensor)
        ierr=exatns_tensor_contract("H(m,e)+=V(m,n,e,f)*T(f,n)",lambda_t_CC2%fbar_ov,int_t%oovv,comm_t_CC2%t1,ONE)

       end subroutine get_fbar_ov_CC2_exatns

       subroutine get_w_vovo_CC2_exatns(w_vovo,comm_t_CC2,oovv_tensor,vvoo_id,vvoo_root,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_comm_tens_t

        type(tens_rcrsv_t),intent(inout) :: w_vovo
        type(exatns_comm_tens_t),intent(inout) :: comm_t_CC2
        type(tens_rcrsv_t),intent(inout) :: oovv_tensor
        type(tens_rcrsv_t) :: t1t1_tensor
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer(INTD) :: ierr
        integer(INTD), dimension(4) :: vvoo_id
        integer(INTL), dimension(4) :: vvoo_root
!!!!!!!!!!!!!       Auxiliary tensor
         ierr=exatns_tensor_create(t1t1_tensor,"t1t1",vvoo_id,vvoo_root,EXA_DATA_KIND_C8)


!        \W bmje = \V bmje +\SUMs f \V mbef \ts fj - \SUMs n \V mnej \ts bn +
!                  +\SUMd nf \V mnef (\td fbjn - \ts fj \ts bn)

        ierr=exatns_tensor_contract("W(b,m,e,j)+=H(b,m,e,j)*K()",w_vovo,comm_t_CC2%h_int,one_tensor)
        ierr=exatns_tensor_init(t1t1_tensor,'ZERO')
        ierr=exatns_tensor_contract("T(f,b,j,n)+=T(f,j)*T(b,n)",t1t1_tensor,comm_t_CC2%t1,comm_t_CC2%t1,MINUS_ONE)
!ReDo these Lines for Check CC2=CCSD if t2 \neg 0 - That's why t1t1 is negative        !
        ierr=exatns_tensor_contract("T(f,b,j,n)+=R(f,b,j,n)*K()",t1t1_tensor,comm_t_CC2%t2,one_tensor)
!End of
        ierr=exatns_tensor_contract("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,t1t1_tensor)

       end subroutine get_w_vovo_CC2_exatns

       subroutine get_w_oooo_CC2_exatns(w_oooo_CC2,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens_t,exatns_comm_tens_t

        type(tens_rcrsv_t),intent(inout) :: w_oooo_CC2
        type(exatns_comm_tens_t),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens_t),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

!       oooo ! validated
!        ierr=exatns_tensor_init(w_oooo_CC2)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,j)*K()",w_oooo_CC2,int_t%oooo,one_tensor)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1) ! o(4)v(1)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1, &
                                    MINUS_ONE) ! o(4)v(1)
        ierr=exatns_tensor_contract("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",w_oooo_CC2,int_t%oovv,comm_t_CC2%tau, &
                                    ONE_HALF)

       end subroutine get_w_oooo_CC2_exatns

       subroutine get_w_vvvv_CC2_exatns(w_vvvv_CC2,comm_t_CC2,int_t,one_tensor)

        use exatensor
        use exacorr_datatypes, only : exatns_intg_tens_t,exatns_comm_tens_t
        
        type(tens_rcrsv_t),intent(inout) :: w_vvvv_CC2
        type(exatns_comm_tens_t),intent(inout) :: comm_t_CC2
        type(exatns_intg_tens_t),intent(inout) :: int_t
        type(tens_rcrsv_t),intent(inout) :: one_tensor
        integer :: ierr

! \W efab=\V efab -\Pm ef \Vmfab \ts em + \undemi \V mnab \taup efmn
!        ierr=exatns_tensor_init(w_vvvv_CC2)
        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(e,f,a,b)*K()",w_vvvv_CC2,int_t%vvvv,one_tensor)

        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1)
        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1,MINUS_ONE)

        ierr=exatns_tensor_contract("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv_CC2,int_t%oovv,comm_t_CC2%tau,ONE_HALF)

       end subroutine get_w_vvvv_CC2_exatns
 
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  End CC2 Exatens
!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!
!!!!!!!    Start CC2 Talsh
!!!!!!!
       subroutine get_tauCC2_eti(comm_t_CC2)  !L

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        integer :: ierr

        logical :: shunt_CC2=.false.

!       tau ! validated
        ierr=initialize_tensor_t(comm_t_CC2%tau)


!       Here I do a 'usual tau' to check
        if(shunt_CC2) then
                ierr=t_einsum("T(a,b,i,j)=H(a,b,i,j)",comm_t_CC2%tau,comm_t_CC2%t2)
                print*, " "
                print*, " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING SHUNT IN TAU_CC2  "
                print*, " "
                print*, "t2 comm CC2 - Subroutine   get_tauCC2 access to t2   "
                print*, " "
        end if

        ierr=t_einsum("T(a,b,i,j)+=T(a,i)*T(b,j)",comm_t_CC2%tau,comm_t_CC2%t1,comm_t_CC2%t1)                  ! o(2)v(2)
        ierr=t_einsum("T(a,b,i,j)+=T(a,j)*T(b,i)",comm_t_CC2%tau,comm_t_CC2%t1,comm_t_CC2%t1,scale=MINUS_ONE)  ! o(2)v(2)

       end subroutine get_tauCC2_eti

       subroutine get_fbar_oo_CC2_eti(fbar_oo,comm_t_CC2,int_t)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t
      
        type(tensor_t),intent(inout)   :: fbar_oo
        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        type(eti_intg_tens_t),intent(inout) :: int_t
        integer :: ierr

!      \Ft im =\f im +\SUMs e \f ie \ts em +\SUMd en \V inme \ts en +\ undemi
!      \SUMt mef \V inef \taup efmn  (\tau in CC2)  (Eq. A4)

        ierr=initialize_tensor_t(fbar_oo)
        ierr=t_einsum("H(i,m)+=F(i,m)",fbar_oo,comm_t_CC2%foo)
        ierr=t_einsum("H(i,m)+=F(i,e)*T(e,m)",fbar_oo,comm_t_CC2%fov,comm_t_CC2%t1,scale=ONE)
        ierr=t_einsum("H(i,m)+=V(i,n,m,e)*T(e,n)",fbar_oo,int_t%ooov,comm_t_CC2%t1,scale=ONE)
        ierr=t_einsum("H(i,m)+=V(i,n,e,f)*T(e,f,m,n)",fbar_oo,int_t%oovv,comm_t_CC2%tau,scale=ONE_HALF)

       end subroutine get_fbar_oo_CC2_eti

       subroutine get_fbar_vv_CC2_eti(fbar_vv,comm_t_CC2,int_t)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t

!        type(eti_lambda_tens_t),intent(inout) :: lambda_t_CC2
        type(tensor_t),intent(inout)   :: fbar_vv
        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        type(eti_intg_tens_t),intent(inout) :: int_t
        integer :: ierr

!       \Ft ea =\f ea -\SUMs \f ma \ts em - \SUMd mf \V mefa \ts fm - \undemi
!       \Vmnaf \taup efmn (\tau in CC2)  (Eq. A5)

        ierr=initialize_tensor_t(fbar_vv)
        ierr=t_einsum("H(e,a)+=F(e,a)",fbar_vv,comm_t_CC2%fvv)
        ierr=t_einsum("H(e,a)+=F(m,a)*T(e,m)",fbar_vv,comm_t_CC2%fov,comm_t_CC2%t1,scale=MINUS_ONE)
        ierr=t_einsum("H(e,a)+=V(e,m,f,a)*T(f,m)",fbar_vv,int_t%vovv,comm_t_CC2%t1,scale=MINUS_ONE)
        ierr=t_einsum("H(e,a)+=V(m,n,a,f)*T(e,f,m,n)",fbar_vv,int_t%oovv,comm_t_CC2%tau,scale=MINUS_ONE_HALF)

       end subroutine get_fbar_vv_CC2_eti

       subroutine get_fbar_ov_CC2_eti(fbar_ov,comm_t_CC2,int_t)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t

        type(tensor_t),intent(inout)   :: fbar_ov
        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        type(eti_intg_tens_t),intent(inout) :: int_t
        integer :: ierr

!      \Ft me =\f me +\SUMd nf \V mnef \ts fn (Eq. A6)
        ierr=initialize_tensor_t(fbar_ov)
        ierr=t_einsum("H(m,e)+=F(m,e)",fbar_ov,comm_t_CC2%fov)
        ierr=t_einsum("H(m,e)+=V(m,n,e,f)*T(f,n)",fbar_ov,int_t%oovv,comm_t_CC2%t1,scale=ONE)

       end subroutine get_fbar_ov_CC2_eti

       subroutine get_w_vovo_CC2_eti(w_vovo,comm_t_CC2,oovv_tensor,nocc,nvir)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_comm_tens_t

        type(tensor_t),intent(inout) :: w_vovo
        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        type(tensor_t),intent(inout) :: oovv_tensor
        integer,intent(in) :: nocc, nvir
        type(tensor_t) :: t1t1_aux_tensor
        integer(INTD) :: vvoo_dims(4)
        integer(INTD) :: ierr

        logical :: shunt_CC2 =.false.

!       Auxiliary tensor
        vvoo_dims(1:2) = nvir
        vvoo_dims(3:4) = nocc
        ierr=construct_tensor_t(t1t1_aux_tensor,C8,vvoo_dims,init_val=ZERO)

!        \W bmje = \V bmje +\SUMs f \V mbef \ts fj - \SUMs n \V mnej \ts bn +
!                  +\SUMd nf \V mnef (\td fbjn - \ts fj \ts bn)

        ierr=t_einsum("W(b,m,e,j)+=H(b,m,e,j)",w_vovo,comm_t_CC2%h_int)
        ierr=initialize_tensor_t(t1t1_aux_tensor)
        ierr=t_einsum("T(f,b,j,n)+=T(f,j)*T(b,n)",t1t1_aux_tensor,comm_t_CC2%t1,comm_t_CC2%t1,scale=ONE)


!ReDo these Lines for Check CC2=CCSD if t2 \neg 0 - That's why t1t1 is negative        !
        if(shunt_CC2) then
                print*, " "
                print*, " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING SHUNT IN W VOVO_CC2  "
                print*, " "
                ierr=t_einsum("T(f,b,j,n)+=R(f,b,j,n)",t1t1_aux_tensor,comm_t_CC2%t2,scale=ONE)
        end if
!End of

        ierr=t_einsum("W(b,m,e,j)+=V(m,n,e,f)*T(f,b,j,n)",w_vovo,oovv_tensor,t1t1_aux_tensor,scale=ONE)



       end subroutine get_w_vovo_CC2_eti

       subroutine get_w_oooo_CC2_eti(w_oooo_CC2,comm_t_CC2,int_t)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t

        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        type(tensor_t),intent(inout) :: w_oooo_CC2
        type(eti_intg_tens_t),intent(inout) :: int_t
        integer :: ierr

!       \W mnij = \V mnij + \Pm mn \SUMs e \V ijen \ts em + \undemi \SUMd ef \V ijef \taup efmn

        ierr=t_einsum("A(k,l,i,j)+=V(k,l,i,j)",w_oooo_CC2,int_t%oooo)
        ierr=t_einsum("A(k,l,i,j)+=V(k,l,i,c)*T(c,j)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1) ! o(4)v(1)
        ierr=t_einsum("A(k,l,i,j)+=V(k,l,j,c)*T(c,i)",w_oooo_CC2,int_t%ooov,comm_t_CC2%t1, &
                                    scale=MINUS_ONE) ! o(4)v(1)
        ierr=t_einsum("A(k,l,i,j)+=V(k,l,c,d)*T(c,d,i,j)",w_oooo_CC2,int_t%oovv,comm_t_CC2%tau, &
                                    scale=ONE_HALF)

       end subroutine get_w_oooo_CC2_eti

       subroutine get_w_vvvv_CC2_eti(w_vvvv_CC2,comm_t_CC2,int_t)

        use exacorr_tensor_interface
        use exacorr_datatypes, only : eti_intg_tens_t,eti_comm_tens_t
        
        type(tensor_t),intent(inout) :: w_vvvv_CC2
        type(eti_comm_tens_t),intent(inout) :: comm_t_CC2
        type(eti_intg_tens_t),intent(inout) :: int_t
        integer :: ierr

! \W efab=\V efab -\Pm ef \Vmfab \ts em + \undemi \V mnab \taup efmn
        ierr=t_einsum("W(e,f,a,b)+=V(e,f,a,b)",w_vvvv_CC2,int_t%vvvv)

        ierr=t_einsum("W(e,f,a,b)+=V(f,m,a,b)*T(e,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1)
        ierr=t_einsum("W(e,f,a,b)+=V(e,m,a,b)*T(f,m)",w_vvvv_CC2,int_t%vovv,comm_t_CC2%t1,scale=MINUS_ONE)

        ierr=t_einsum("W(e,f,a,b)+=V(m,n,a,b)*T(e,f,m,n)",w_vvvv_CC2,int_t%oovv,comm_t_CC2%tau,scale=ONE_HALF)

       end subroutine get_w_vvvv_CC2_eti
       
end module intermediates_cc2
