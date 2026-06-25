module intermediates_eom_rsp

!Module to make intermediates needed in CC-CI(EOM-CC) response equations

     use exacorr_tensor_interface
     use tensor_algebra
     use exacorr_datatypes
     use exacorr_scaling_parameters
     use eti_common_routines

     implicit none
     private

     interface get_epsilon_prop_vo
        module procedure get_epsilon_prop_vo_eti
     end interface get_epsilon_prop_vo

     interface get_epsilon_prop_vvoo
        module procedure get_epsilon_prop_vvoo_eti
     end interface get_epsilon_prop_vvoo

     interface get_O_vo
        module procedure get_O_vo_eti
     end interface get_O_vo

     interface get_eta_prop_ov
        module procedure get_eta_prop_ov_eti
     end interface get_eta_prop_ov

     interface get_eta_prop_ov_eom
        module procedure get_eta_prop_ov_eom_eti
     end interface get_eta_prop_ov_eom

     interface get_eta_prop_oovv
        module procedure get_eta_prop_oovv_eti
     end interface get_eta_prop_oovv

     interface get_prop_bar_oo
        module procedure get_prop_bar_oo_eti
     end interface get_prop_bar_oo

     interface get_prop_bar_vv
        module procedure get_prop_bar_vv_eti
     end interface get_prop_bar_vv

     interface get_prop_bar_ov
        module procedure get_prop_bar_ov_eti
     end interface get_prop_bar_ov

     interface get_Q_ovoo
        module procedure get_Q_ovoo_eti
     end interface get_Q_ovoo

     interface get_Q_vvvo
        module procedure get_Q_vvvo_eti
     end interface get_Q_vvvo


!      interface form_ccsd_intermediates
!         module procedure form_ccsd_intermediates_eti
!      end interface form_ccsd_intermediates

   !   interface destroy_eomccsd_rsp_intermediates
   !      module procedure destroy_eomccsd_rsp_intermediates_eti
   !   end interface destroy_eomccsd_rsp_intermediates

!    top-level driver for forming intermediates with CCSD wavefunctions
!      public form_ccsd_intermediates
   !   public destroy_eomccsd_rsp_intermediates
     public get_epsilon_prop_vo, get_epsilon_prop_vvoo, get_O_vo, get_eta_prop_ov, get_eta_prop_oovv
     public get_prop_bar_oo, get_prop_bar_vv, get_prop_bar_ov, get_Q_ovoo, get_Q_vvvo
     public get_eta_prop_ov_eom

     contains


!-----------------------------------------------------------------------------------------------------------------------
     subroutine get_epsilon_prop_vo_eti(epsilon_prop_vo_tensor,orb_int_tensor,t2_tensor,t1_tensor,LR_Aux_oo,vo_dims)

         type(tensor_t), intent(inout) :: epsilon_prop_vo_tensor
         type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
         type(tensor_t), intent(inout) :: t2_tensor,t1_tensor
         type(tensor_t), intent(inout)            :: LR_Aux_oo
         type(tensor_t)            :: Transpose_vo
         type(tensor_t) :: one_tensor
         integer(INTD), dimension(2)   :: vo_dims
         integer(C_INT)     :: one_dims(1)
         integer :: ierr
         type(tensor_t) :: uncon_tensor
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = t2_tensor%datatype

         one_dims(1) = 1
         ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ONE)
         ierr=construct_tensor_t(epsilon_prop_vo_tensor,eti_precision%storage_t,vo_dims,init_val=ZERO)
         ! if(print_level.ge.2) print *,'LH_eti LR_Aux_vo=Y ai ', ierr
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(epsilon_prop_vo_tensor),epsilon_prop_vo_tensor%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(orb_int_tensor%prop_vo),orb_int_tensor%prop_vo%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(t1_tensor),t1_tensor%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(t2_tensor),t2_tensor%datatype
         print *,'debug in get_epsilon_prop_vo_eti',eti_precision%storage_t
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(orb_int_tensor%prop_vv),orb_int_tensor%prop_vv%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(orb_int_tensor%prop_oo),orb_int_tensor%prop_oo%datatype

!LH     \Chis ai = +\Y ai
         ierr=t_einsum('C(a,i)+=Y(a,i)*K()',epsilon_prop_vo_tensor,orb_int_tensor%prop_vo,one_tensor)
         ! if(print_level.ge.2) print *,'LH_eti epsilon_prop_vo_tensor ierr ', ierr

!LH     \Chis ai += +\Y ae \ts ei 
         ierr=t_einsum('C(a,i)+=Y(a,e)*T(e,i)',epsilon_prop_vo_tensor,orb_int_tensor%prop_vv,t1_tensor,scale=ONE)
         ! if(print_level.ge.2) print *,'LH_eti epsilon_prop_vo_tensor ierr ', ierr

!LH     \Chis ai += -\Y mi \ts am
         ierr=t_einsum('C(a,i)+=Y(m,i)*T(a,m)',epsilon_prop_vo_tensor,orb_int_tensor%prop_oo,t1_tensor,scale=MINUS_ONE)
         ! if(print_level.ge.2) print *,'LH_eti epsilon_prop_vo_tensor ierr ', ierr

!LH     \Chis ai += - \Y me \ts ei \ts am
         ierr=initialize_tensor_t(LR_Aux_oo)
         ierr=t_einsum('A(m,i)+=Y(m,e)*T(e,i)',LR_Aux_oo,orb_int_tensor%prop_ov,t1_tensor)
         ierr=t_einsum('C(a,i)+=A(m,i)*T(a,m)',epsilon_prop_vo_tensor,LR_Aux_oo,t1_tensor,scale=MINUS_ONE)
         ! if(print_level.ge.2) print *,'LH_eti epsilon_prop_vo_tensor ierr ', ierr

!LH     \Chis ai += +\Y me \td aeim
         ierr=t_einsum('C(a,i)+=Y(m,e)*T(a,e,i,m)',epsilon_prop_vo_tensor,orb_int_tensor%prop_ov,t2_tensor,scale=ONE)
         ! if(print_level.ge.2) print *,'LH_eti epsilon_prop_vo_tensor ierr ', ierr
         ierr=destroy_tensor_t(one_tensor)

! !XY     \Chis ai += +\ts ai \Ys me \ts em  (Unconnected)
!          ierr=t_einsum('A()+=Y(m,e)*T(e,m)',uncon_tensor,orb_int_tensor%prop_ov,t1_tensor)
!          ierr=t_einsum('C(a,i)+=A()*T(a,i)',epsilon_prop_vo_tensor,uncon_tensor,t1_tensor,scale=ONE)

!          ! if(print_level.ge.2) print *,"tensor_norm2(epsilon_prop_vo_tensor)",tensor_norm2(epsilon_prop_vo_tensor)

     end subroutine get_epsilon_prop_vo_eti



     subroutine get_epsilon_prop_vvoo_eti(epsilon_prop_vvoo_tensor,orb_int_tensor, &
                                    t2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,vvoo_dims)
             
               type(tensor_t), intent(inout) :: epsilon_prop_vvoo_tensor
               type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
               type(tensor_t), intent(inout) :: t2_tensor,t1_tensor
               type(tensor_t), intent(inout)            :: LR_Aux_oo,LR_Aux_vv
               integer(INTD), dimension(4)   :: vvoo_dims
               integer :: ierr
               logical :: debug_local = .false.
               type(tensor_t) :: uncon_tensor
               integer(C_INT)     :: one_dims(1)
               type(eti_precision_t) :: eti_precision

           eti_precision%storage_t = t2_tensor%datatype


               
           ierr=construct_tensor_t(epsilon_prop_vvoo_tensor,eti_precision%storage_t,vvoo_dims,init_val=ZERO)

         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(epsilon_prop_vvoo_tensor),epsilon_prop_vvoo_tensor%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(orb_int_tensor%prop_vo),orb_int_tensor%prop_vo%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(t1_tensor),t1_tensor%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(t2_tensor),t2_tensor%datatype
         print *,'debug in get_epsilon_prop_vo_eti',eti_precision%storage_t
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(orb_int_tensor%prop_vv),orb_int_tensor%prop_vv%datatype
         print *,'debug in get_epsilon_prop_vo_eti',tensor_norm2(orb_int_tensor%prop_oo),orb_int_tensor%prop_oo%datatype



!LH     \Chid abij = +\Pm ab \Y bf \td afij
         !   ierr=t_einsum('C(a,b,i,j)+=Y(b,f)*T(a,f,i,j)',epsilon_prop_vvoo_tensor,orb_int_tensor%prop_vv,t2_tensor,scale=ONE)
           ierr=t_einsum('C(a,b,i,j)+=Y(b,f)*T(a,f,i,j)', &
                                    epsilon_prop_vvoo_tensor,orb_int_tensor%prop_vv,t2_tensor,scale=ONE)

           ierr=t_einsum('C(a,b,i,j)+=Y(a,f)*T(b,f,i,j)', & 
                                    epsilon_prop_vvoo_tensor,orb_int_tensor%prop_vv,t2_tensor,scale=MINUS_ONE)
           if(debug_local) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr

!LH     \Chid abij += - \Pm ij \Y mj \td abim
           ierr=t_einsum('C(a,b,i,j)+=Y(m,j)*T(a,b,i,m)', &
                                    epsilon_prop_vvoo_tensor,orb_int_tensor%prop_oo,t2_tensor,scale=MINUS_ONE)
         !   if(print_level.ge.2) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr
           ierr=t_einsum('C(a,b,i,j)+=Y(m,i)*T(a,b,j,m)', & 
                                       epsilon_prop_vvoo_tensor,orb_int_tensor%prop_oo,t2_tensor,scale=ONE)        
           if(debug_local) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr


!LH     \Chid abij += - \Pm ij \Y me \ts ei \td abmj
           ierr=initialize_tensor_t(LR_Aux_oo)
           ierr=t_einsum('A(m,i)+=Y(m,e)*T(e,i)', &
                                    LR_Aux_oo,orb_int_tensor%prop_ov,t1_tensor,scale=ONE)
           ierr=t_einsum('C(a,b,i,j)+=A(m,i)*T(a,b,m,j)', &
                                    epsilon_prop_vvoo_tensor,LR_Aux_oo,t2_tensor,scale=MINUS_ONE)
           ierr=t_einsum('C(a,b,i,j)+=A(m,j)*T(a,b,m,i)', &
                                    epsilon_prop_vvoo_tensor,LR_Aux_oo,t2_tensor,scale=ONE)
           if(debug_local) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr

!LH     \Chid abij += - \Pm ab \Y mf \ts am \td fbij
           ierr=initialize_tensor_t(LR_Aux_vv)
           ierr=t_einsum('A(a,f)+=Y(m,f)*T(a,m)', &
                                    LR_Aux_vv,orb_int_tensor%prop_ov,t1_tensor,scale=ONE)
           ierr=t_einsum('C(a,b,i,j)+=A(a,f)*T(f,b,i,j)', &
                                    epsilon_prop_vvoo_tensor,LR_Aux_vv,t2_tensor,scale=MINUS_ONE)
           ierr=t_einsum('C(a,b,i,j)+=A(b,f)*T(f,a,i,j)', &
                                    epsilon_prop_vvoo_tensor,LR_Aux_vv,t2_tensor,scale=ONE)
           if(debug_local) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr

! !XY    \Chid abij += + \td abij \Ys me \ts em (Unconnected)
!            ierr=t_einsum('A()+=Y(m,e)*T(e,m)', &
!                                     uncon_tensor,orb_int_tensor%prop_ov,t1_tensor,scale=ONE)
!            ierr=t_einsum('C(a,b,i,j)+=A()*T(a,b,i,j)', &
!                                     epsilon_prop_vvoo_tensor,uncon_tensor,t2_tensor,scale=ONE)
!            if(debug_local) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr

         !   if(debug_local) print *,"tensor_norm2(epsilon_prop_vvoo_tensor)",tensor_norm2(epsilon_prop_vvoo_tensor)

         !  call print_tensor(epsilon_prop_vvoo_tensor,1.0d-10,"T2 part of epison_prop")

     end subroutine get_epsilon_prop_vvoo_eti



     subroutine get_eta_prop_ov_eti(eta_prop_ov_tensor,orb_int_tensor,l1_tensor,l2_tensor,&
                           t2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,LR_Aux_vo,ov_dims)
               
         type(tensor_t), intent(inout) :: eta_prop_ov_tensor
         type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
         type(tensor_t), intent(inout) :: t2_tensor,t1_tensor
         type(tensor_t), intent(inout) :: l2_tensor,l1_tensor
         type(tensor_t), intent(inout)            :: LR_Aux_oo,LR_Aux_vv,LR_Aux_vo
         type(tensor_t) :: one_tensor
         integer(INTD), dimension(2)   :: ov_dims
         integer :: ierr

         type(tensor_t)            :: Transpose_vo
         integer(INTD), dimension(2)   :: vo_dims
         integer(C_INT)     :: one_dims(1)
         type(tensor_t) :: uncon_tensor
         logical :: debug_local = .false.
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = t2_tensor%datatype


         ierr=construct_tensor_t(eta_prop_ov_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_tran ', ierr

         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(eta_prop_ov_tensor),eta_prop_ov_tensor%datatype
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(orb_int_tensor%prop_vo),orb_int_tensor%prop_vo%datatype
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(t1_tensor),t1_tensor%datatype
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(t2_tensor),t2_tensor%datatype
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(l1_tensor),l1_tensor%datatype
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(l2_tensor),l2_tensor%datatype
         print *,'debug in get_eta_prop_ov_eti',eti_precision%storage_t
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(orb_int_tensor%prop_vv),orb_int_tensor%prop_vv%datatype
         print *,'debug in get_eta_prop_ov_eti',tensor_norm2(orb_int_tensor%prop_oo),orb_int_tensor%prop_oo%datatype


!LH Eat1
!LH     \Eta_1 ia=+\Ys ia
         ierr=t_einsum('E(i,a)+=Y(i,a)',eta_prop_ov_tensor,orb_int_tensor%prop_ov,scale=ONE)                    
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o1 ', ierr

!LH     \Eta_1 ia=+\Ys ea \lambdas ie
         ierr=t_einsum('E(i,a)+=Y(e,a)*L(i,e)',eta_prop_ov_tensor,orb_int_tensor%prop_vv,l1_tensor,scale=ONE)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o3 ', ierr

!LH     \Eta_1 ia=-\Ys im \lambdas ma
         ierr=t_einsum('E(i,a)+=Y(i,m)*L(m,a)',eta_prop_ov_tensor,orb_int_tensor%prop_oo,l1_tensor,scale=MINUS_ONE)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o4 ', ierr

!LH     \Eta_1=-\Ys ma \ts em \lambdas ie
         ierr=initialize_tensor_t(LR_Aux_oo)       
         ierr=t_einsum('H(i,m)+=T(e,m)*L(i,e)',LR_Aux_oo,t1_tensor,l1_tensor,scale=ONE)
         ierr=t_einsum('E(i,a)+=H(i,m)*Y(m,a)',eta_prop_ov_tensor,LR_Aux_oo,orb_int_tensor%prop_ov,scale=MINUS_ONE)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o5 ', ierr

!LH     \Eta_1=-\Ys ie \ts em \llambdas ma
         ierr=initialize_tensor_t(LR_Aux_oo)
         ierr=t_einsum('H(i,m)+=Y(i,e)*T(e,m)',Lr_Aux_oo,orb_int_tensor%prop_ov,t1_tensor,scale=ONE)
         ierr=t_einsum('E(i,a)+=H(i,m)*L(m,a)',eta_prop_ov_tensor,LR_Aux_oo,l1_tensor,scale=MINUS_ONE)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o6 ', ierr

!******************************RED********************************************************************************************** 
! !LH     \Eta_1 = + \ts fm \lambdad imae \Ys ef 
!          ierr=initialize_tensor_t(LR_Aux_vo)
!          ierr=t_einsum('H(e,m)+=T(f,m)*Y(e,f)',LR_Aux_vo,t1_tensor,orb_int_tensor%prop_vv)
!          ierr=t_einsum('E(i,a)+=H(e,m)*L(i,m,a,e)',eta_prop_ov_tensor,LR_Aux_vo,l2_tensor,scale=ONE)
!          if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o18 ', ierr

! !LH      \Eta_1 = - \ts en \lambdad imae \Ys nm 
!          ierr=initialize_tensor_t(LR_Aux_vo)
!          ierr=t_einsum('H(e,m)+=T(e,n)*Y(n,m)',LR_Aux_vo,t1_tensor,orb_int_tensor%prop_oo)
!          ierr=t_einsum('E(i,a)+=H(e,m)*L(i,m,a,e)',eta_prop_ov_tensor,LR_Aux_vo,l2_tensor,scale=MINUS_ONE)
!          if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o19 ', ierr

         
! !LH     \Eta_1 =+\lambdad imae \Ys em 
!          ierr=t_einsum('E(i,a)+=L(i,m,a,e)*A(e,m)',eta_prop_ov_tensor,l2_tensor,Transpose_vo,scale=ONE)
!          if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o20 ', ierr
! !******************************************************************************************************************************* 



!LH     \Eta_1= -one_half \td femn \lambdad mife \Ys na 
         ierr=initialize_tensor_t(LR_Aux_oo)
         ierr=t_einsum('H(i,n)+=T(f,e,m,n)*L(m,i,f,e)',LR_Aux_oo,t2_tensor,l2_tensor)
         ierr=t_einsum('E(i,a)+=H(i,n)*Y(n,a)',& 
                                 eta_prop_ov_tensor,LR_Aux_oo,orb_int_tensor%prop_ov,scale=MINUS_ONE_HALF)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o28 ', ierr

!LH     \Eta_1=-one_half \td  fenm \lambdad nmfa \Ys ie 
         ierr=initialize_tensor_t(LR_Aux_vv)
         ierr=t_einsum('H(e,a)+=T(f,e,n,m)*L(n,m,f,a)',LR_Aux_vv,t2_tensor,l2_tensor)
         ierr=t_einsum('E(i,a)+=H(e,a)*Y(i,e)',&
                                 eta_prop_ov_tensor,LR_Aux_vv,orb_int_tensor%prop_ov,scale=MINUS_ONE_HALF)
         if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o29 ', ierr
         ! ierr=destroy_tensor_t(one_tensor)


! !******************************RED**********************************************************************************************                                 
! ! !LH     \Eta_1=\td efmn \lambdad imae \Ys nf
!          ierr=initialize_tensor_t(LR_Aux_vo)
!          ierr=t_einsum('H(e,m)+=T(e,f,m,n)*Y(n,f)',LR_Aux_vo,t2_tensor,orb_int_tensor%prop_ov)
!          ierr=t_einsum('E(i,a)+=H(e,m)*L(i,m,a,e)',eta_prop_ov_tensor,LR_Aux_vo,l2_tensor,scale=ONE)    
!          if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o30 ', ierr

! ! !LH     \Eta_1=-\ts en \ts fm \lambdad imae \Ys nf 

!          ierr=initialize_tensor_t(LR_Aux_oo)
!          ierr=t_einsum('H(n,m)+=T(f,m)*Y(n,f)',LR_Aux_oo,t1_tensor,orb_int_tensor%prop_ov)
!          ierr=initialize_tensor_t(LR_Aux_vo)
!          ierr=t_einsum('H(e,m)+=T(e,n)*Y(n,m)',LR_Aux_vo,t1_tensor,LR_Aux_oo)
!          ! ! ! ! ! ! ierr=t_einsum('E(i,a)+=H(e,m)*L(i,m,a,e)',eta_prop_ov_tensor,LR_Aux_vo,l2_tensor,scale=ONE)
!          ierr=t_einsum('E(i,a)+=H(e,m)*L(i,m,a,e)',eta_prop_ov_tensor,LR_Aux_vo,l2_tensor,scale=MINUS_ONE) !XY
!          if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr_o39 ', ierr
! !*******************************************************************************************************************************  


! ! !XY     \Eta_1 += +\lambdas kc \Ys ia \ts ai  (Unconnected)
!          ierr=t_einsum('A()+=Y(i,a)*T(a,i)',uncon_tensor,orb_int_tensor%prop_ov,t1_tensor)
!          ierr=t_einsum('E(i,a)+=A()*L(i,a)',eta_prop_ov_tensor,uncon_tensor,l1_tensor,scale=ONE)

!          if(debug_local) print *,"tensor_norm2(epsilon_prop_vo_tensor)",tensor_norm2(eta_prop_ov_tensor)
! !          ierr=destroy_tensor_t(Transpose_vo)


         ! ierr=initialize_tensor_t(LR_Aux_vvoo)
         ! ierr=t_einsum('H(e,f,n,m)+=T(e,n)*T(f,m)',LR_Aux_vvoo,t1_tensor,t1_tensor)
         ! ierr=initialize_tensor_t(LR_Aux_vo)
         ! ierr=t_einsum('H(e,m)+=T(e,f,m,n)*Y(i,e)',LR_Aux_vo,LR_Aux_vvoo,orb_int_tensor%prop_ov)
         ! ierr=t_einsum('E(i,a)+=H(e,m)*L(i,m,a,e)',eta_prop_ov_tensor,LR_Aux_vo,l2_tensor,scale=ONE)
         ! ierr=destroy_tensor_t(Transpose_vo)


     end subroutine get_eta_prop_ov_eti



     subroutine get_eta_prop_oovv_eti(eta_prop_oovv_tensor,orb_int_tensor,&
                        l1_tensor,l2_tensor,t1_tensor,LR_Aux_oo,LR_Aux_vv,oovv_dims)

      type(tensor_t), intent(inout) :: eta_prop_oovv_tensor
      type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
      type(tensor_t), intent(inout) :: t1_tensor
      type(tensor_t), intent(inout) :: l2_tensor,l1_tensor
      type(tensor_t), intent(inout)            :: LR_Aux_oo,LR_Aux_vv
      integer(INTD), dimension(4)   :: oovv_dims
      integer :: ierr
      logical :: debug_local = .false.
      integer(C_INT)     :: one_dims(1)
      type(tensor_t) :: uncon_tensor
      type(eti_precision_t) :: eti_precision

      eti_precision%storage_t = l2_tensor%datatype


      ierr=construct_tensor_t(eta_prop_oovv_tensor,eti_precision%storage_t,oovv_dims,init_val=ZERO)

!LH     \Eta_2=+ \Pmm ij \Pmm ab \lambdas ia \Ys jb
      ierr=t_einsum('E(i,j,a,b)+=L(i,a)*Y(j,b)',eta_prop_oovv_tensor,l1_tensor,orb_int_tensor%prop_ov,scale=ONE)
      ierr=t_einsum('E(i,j,a,b)+=L(j,a)*Y(i,b)',eta_prop_oovv_tensor,l1_tensor,orb_int_tensor%prop_ov,scale=MINUS_ONE)
      ierr=t_einsum('E(i,j,a,b)+=L(i,b)*Y(j,a)',eta_prop_oovv_tensor,l1_tensor,orb_int_tensor%prop_ov,scale=MINUS_ONE)
      ierr=t_einsum('E(i,j,a,b)+=L(j,b)*Y(i,a)',eta_prop_oovv_tensor,l1_tensor,orb_int_tensor%prop_ov,scale=ONE)
      if(debug_local) print *,'LH_eti eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=+\Pmm ab \lambdad ijae \Ys eb 
      ierr=t_einsum('E(i,j,a,b)+=L(i,j,a,e)*Y(e,b)', &
                                 eta_prop_oovv_tensor,l2_tensor,orb_int_tensor%prop_vv,scale=ONE)
      ierr=t_einsum('E(i,j,a,b)+=L(i,j,b,e)*Y(e,a)', & 
                                 eta_prop_oovv_tensor,l2_tensor,orb_int_tensor%prop_vv,scale=MINUS_ONE)
      if(debug_local) print *,'LH_eti eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2= - \Pmm ij \lambdad imab \Ys jm 
      ierr=t_einsum('E(i,j,a,b)+=L(i,m,a,b)*Y(j,m)', &
                                 eta_prop_oovv_tensor,l2_tensor,orb_int_tensor%prop_oo,scale=MINUS_ONE)
      ierr=t_einsum('E(i,j,a,b)+=L(j,m,a,b)*Y(i,m)', &
                                 eta_prop_oovv_tensor,l2_tensor,orb_int_tensor%prop_oo,scale=ONE)
      if(debug_local) print *,'LH_eti eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=-\Pmm ab \ts em \lambdad ijae \Ys mb 
      ierr=initialize_tensor_t(LR_Aux_vv)
      ierr=t_einsum('H(e,b)+=T(e,m)*Y(m,b)',LR_Aux_vv,t1_tensor,orb_int_tensor%prop_ov)
      ierr=t_einsum('E(i,j,a,b)+=H(e,b)*L(i,j,a,e)',eta_prop_oovv_tensor,LR_Aux_vv,l2_tensor,scale=MINUS_ONE)         
      ierr=t_einsum('E(i,j,a,b)+=H(e,a)*L(i,j,b,e)',eta_prop_oovv_tensor,LR_Aux_vv,l2_tensor,scale=ONE)         
      if(debug_local) print *,'LH_eti eta_prop_oovv_tensor ierr ', ierr

!LH     \Eta_2=-\Pmm ij \ts em \lambdad imab\Ys je 
      ierr=initialize_tensor_t(LR_Aux_oo)
      ierr=t_einsum('H(j,m)+=T(e,m)*Y(j,e)',LR_Aux_oo,t1_tensor,orb_int_tensor%prop_ov)
      ierr=t_einsum('E(i,j,a,b)+=H(j,m)*L(i,m,a,b)',eta_prop_oovv_tensor,LR_Aux_oo,l2_tensor,scale=MINUS_ONE)
      ierr=t_einsum('E(i,j,a,b)+=H(i,m)*L(j,m,a,b)',eta_prop_oovv_tensor,LR_Aux_oo,l2_tensor,scale=ONE)
      if(debug_local) print *,'LH_eti eta_prop_ov_tensor ierr ', ierr


! !XY    \Eta_2 += + \lambda klcd \Ys ia \ts ai (Unconnected)
!            ierr=t_einsum('A()+=Y(i,a)*T(a,i)', &
!                                     uncon_tensor,orb_int_tensor%prop_ov,t1_tensor,scale=ONE)
!            ierr=t_einsum('E(i,j,a,b)+=A()*L(i,j,a,b)', &
!                                     eta_prop_oovv_tensor,uncon_tensor,l2_tensor,scale=ONE)
!            if(debug_local) print *,'LH_eti epsilon_prop_vvoo_tensor ierr ', ierr

      ! if(debug_local) print*,"tensor_norm2(eta_prop_oovv_tensor)",tensor_norm2(eta_prop_oovv_tensor)

      ! call print_tensor(eta_prop_oovv_tensor,1.0d-10,"T2 part of eta_prop3_after")
     end subroutine get_eta_prop_oovv_eti


     subroutine get_eta_prop_ov_eom_eti(eta_prop_ov_eom_tensor,eta_prop_ov_tensor,epsilon_prop_vo_tensor,&
                                    l2_tensor,ov_dims)

         type(tensor_t), intent(inout) :: eta_prop_ov_eom_tensor
         type(tensor_t), intent(inout) :: eta_prop_ov_tensor
         type(tensor_t), intent(inout) :: epsilon_prop_vo_tensor
         type(tensor_t), intent(inout) :: l2_tensor
         integer(INTD), dimension(2)   :: ov_dims
         type(C_PTR) :: body_p, body_p1
         complex(8), pointer, contiguous:: etaX_tens(:,:), etaXeom_tens(:,:)
         integer            :: ierr
         type(eti_precision_t) :: eti_precision

         eti_precision%storage_t = l2_tensor%datatype

         ierr=construct_tensor_t(eta_prop_ov_eom_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
         ! eta_prop_ov_eom_tensor = eta_prop_ov_tensor

         ! ierr=direct_data_pointer(eta_prop_ov_eom_tensor,body_p)
         ! ierr=direct_data_pointer(eta_prop_ov_tensor,body_p1)

         ! call c_f_pointer(body_p,etaXeom_tens,ov_dims)
         ! call c_f_pointer(body_p1,etaX_tens,ov_dims)

         ! etaXeom_tens = etaX_tens

         ierr=t_einsum('nX(i,a)+=n(i,a)',eta_prop_ov_eom_tensor,eta_prop_ov_tensor,scale=ONE)
         ierr=t_einsum("nX(i,a)+=Tl(i,j,a,b)*Epis(b,j)",eta_prop_ov_eom_tensor,& 
                        l2_tensor,epsilon_prop_vo_tensor)

         ! print *, "ierr: ", ierr
         ! call print_tensor(eta_prop_ov_tensor,1.0d-18,"eta_prop_ov in etaeom ") 
         ! call print_tensor(epsilon_prop_vo_tensor,1.0d-18,"T1 part of epsilon_prop_vo in etarom")  
         ! call print_tensor(l2_tensor,1.0d-18,"T2 part of lambda in etaeom")   
         ! call print_tensor(eta_prop_ov_eom_tensor,1.0d-18,"eta_ov_eom")

     end subroutine get_eta_prop_ov_eom_eti


     subroutine get_O_vo_eti(O_vo_tensor,orb_int_tensor,vo_dims)
        
        type(tensor_t), intent(inout) :: O_vo_tensor
        type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
        integer(INTD), dimension(2)   :: vo_dims
        integer :: ierr

        type(tensor_t) :: one_tensor
        integer(C_INT)     :: one_dims(1)
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = orb_int_tensor%prop_ov%datatype

        one_dims(1) = 1
        ierr=construct_tensor_t(one_tensor,eti_precision%storage_t,one_dims(1:0),init_val=ONE)

        ierr=construct_tensor_t(O_vo_tensor,eti_precision%storage_t,vo_dims,init_val=ZERO)
        ierr=t_einsum("O(a,i)+=X+(i,a)*K()",O_vo_tensor,orb_int_tensor%prop_ov,one_tensor)

     end subroutine get_O_vo_eti


     subroutine get_prop_bar_oo_eti(prop_bar_oo_tensor,orb_int_tensor,t1_tensor,oo_dims)

        type(tensor_t), intent(inout) :: prop_bar_oo_tensor
        type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
        type(tensor_t), intent(inout) :: t1_tensor
        integer(INTD), dimension(2)   :: oo_dims
        integer :: ierr
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = t1_tensor%datatype


        ierr=construct_tensor_t(prop_bar_oo_tensor,eti_precision%storage_t,oo_dims,init_val=ZERO)
        ierr=t_einsum("Xbar(m,i)+=X(m,i)",prop_bar_oo_tensor,orb_int_tensor%prop_oo)
      !   ierr=t_einsum("Xbar(m,i)+=X(m,e)*t1(e,i)",prop_bar_oo_tensor,orb_int_tensor%prop_ov,t1_tensor)

     end subroutine get_prop_bar_oo_eti


     subroutine get_prop_bar_vv_eti(prop_bar_vv_tensor,orb_int_tensor,t1_tensor,vv_dims)

        type(tensor_t), intent(inout) :: prop_bar_vv_tensor
        type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
        type(tensor_t), intent(inout) :: t1_tensor
        integer(INTD), dimension(2)   :: vv_dims
        integer :: ierr
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = t1_tensor%datatype


        ierr=construct_tensor_t(prop_bar_vv_tensor,eti_precision%storage_t,vv_dims,init_val=ZERO)
        ierr=t_einsum("Xbar(a,e)+=X(a,e)",prop_bar_vv_tensor,orb_int_tensor%prop_vv)
      !   ierr=t_einsum("Xbar(a,e)+=X(m,e)*t1(a,m)",prop_bar_vv_tensor,orb_int_tensor%prop_ov,t1_tensor,scale=MINUS_ONE)

     end subroutine get_prop_bar_vv_eti


     subroutine get_prop_bar_ov_eti(prop_bar_ov_tensor,orb_int_tensor,ov_dims)

        type(tensor_t), intent(inout) :: prop_bar_ov_tensor
        type(eti_orb_prop_integ_tens_t), intent(inout) :: orb_int_tensor
        integer(INTD), dimension(2)   :: ov_dims
        integer :: ierr
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = orb_int_tensor%prop_ov%datatype

        ierr=construct_tensor_t(prop_bar_ov_tensor,eti_precision%storage_t,ov_dims,init_val=ZERO)
        ierr=t_einsum("Xbar(m,e)+=X(m,e)",prop_bar_ov_tensor,orb_int_tensor%prop_ov)

     end subroutine get_prop_bar_ov_eti


     subroutine get_Q_ovoo_eti(Q_ovoo_tensor,prop_bar_ov_tensor,t2_tensor,ovoo_dims)

        type(tensor_t), intent(inout) :: Q_ovoo_tensor
        type(tensor_t), intent(inout) :: prop_bar_ov_tensor
        type(tensor_t), intent(inout) :: t2_tensor
        integer(INTD), dimension(4)   :: ovoo_dims
        integer :: ierr
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = t2_tensor%datatype

        ierr=construct_tensor_t(Q_ovoo_tensor,eti_precision%storage_t,ovoo_dims,init_val=ZERO)
        ierr=t_einsum("Q(m,b,i,j)=Xbar(m,f)*t2(b,f,i,j)",Q_ovoo_tensor,prop_bar_ov_tensor,t2_tensor)

     end subroutine 


     subroutine get_Q_vvvo_eti(Q_vvvo_tensor,prop_bar_ov_tensor,t2_tensor,vvvo_dims)

        type(tensor_t), intent(inout) :: Q_vvvo_tensor
        type(tensor_t), intent(inout) :: prop_bar_ov_tensor
        type(tensor_t), intent(inout) :: t2_tensor
        integer(INTD), dimension(4)   :: vvvo_dims
        integer :: ierr
        type(eti_precision_t) :: eti_precision

        eti_precision%storage_t = t2_tensor%datatype

        ierr=construct_tensor_t(Q_vvvo_tensor,eti_precision%storage_t,vvvo_dims,init_val=ZERO)
        ierr=t_einsum("Q(a,b,e,j)=Xbar(n,e)*t2(a,b,j,n)",Q_vvvo_tensor,prop_bar_ov_tensor,t2_tensor)
      !   ierr=t_einsum("Q(a,b,e,j)=Xbar(n,e)*t2(a,b,j,n)",Q_vvvo_tensor,prop_bar_ov_tensor,t2_tensor,scale=MINUS_ONE)

     end subroutine 
       
end module
