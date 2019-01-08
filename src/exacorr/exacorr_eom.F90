module exacorr_eom

#if defined (VAR_MPI)

        use exatensor
        use exacorr_datatypes
        use intermediates

        implicit none
        complex(8), parameter :: MINUS_ONE=(-1.D0,0.D0),ONE_HALF=(0.5D0,0.D0)

        private

        public exacorr_eom_driver

       contains

        subroutine exacorr_eom_driver (exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor
!        Type of EOM
         character(len=2) :: flavour

!        tensor ids + roots
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(2)       :: oo_id, ov_id, vv_id
         integer(INTL),dimension(2) :: oo_root, ov_root, vv_root
         integer,dimension(4)       :: oooo_id, ovoo_id, ooov_id, vovo_id, vvvo_id, vovv_id, vvvv_id
         integer(INTL),dimension(4) :: oooo_root, ovoo_root, ooov_root, vovo_root, vvvo_root, vovv_root, vvvv_root

         logical :: CCD
         integer :: ierr

!        Make some noise so that we know we are here
         call interface_mpi_comm_rank(global_communicator,my_MPI_rank)
         print*, " Starting EOM-CC calculation with exatensor"

         CCD = exa_input%ccd

!        Create all ids and roots
         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oooo_id = (/occ_space_id, occ_space_id, occ_space_id, occ_space_id/)
         ovoo_id = (/occ_space_id, vir_space_id, occ_space_id, occ_space_id/)
         ooov_id = (/occ_space_id, occ_space_id, occ_space_id, vir_space_id/)
         vovo_id = (/vir_space_id, occ_space_id, vir_space_id, occ_space_id/)
         vvvo_id = (/vir_space_id, vir_space_id, vir_space_id, occ_space_id/)
         vovv_id = (/vir_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvvv_id = (/vir_space_id, vir_space_id, vir_space_id, vir_space_id/)
         oooo_root = (/occ_space_root, occ_space_root, occ_space_root, occ_space_root/)
         ovoo_root = (/occ_space_root, vir_space_root, occ_space_root, occ_space_root/)
         ooov_root = (/occ_space_root, occ_space_root, occ_space_root, vir_space_root/)
         vovo_root = (/vir_space_root, occ_space_root, vir_space_root, occ_space_root/)
         vvvo_root = (/vir_space_root, vir_space_root, vir_space_root, occ_space_root/)
         vovv_root = (/vir_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvvv_root = (/vir_space_root, vir_space_root, vir_space_root, vir_space_root/)

!        Construct fixed intermediates
         ierr=exatns_tensor_create(lambda_t%fbar_oo,"fbar_oo",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%fbar_ov,"fbar_ov",ov_id,ov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%fbar_vv,"fbar_vv",vv_id,vv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_oooo,"w_oooo",oooo_id,oooo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_ovoo,"w_ovoo",ovoo_id,ovoo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_ooov,"w_ooov",ooov_id,ooov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vovo,"w_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%wbar_vovo,"wbar_vovo",vovo_id,vovo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vvvo,"w_vvvo",vvvo_id,vvvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vovv,"w_vovv",vovv_id,vovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%w_vvvv,"w_vvvv",vvvv_id,vvvv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%goo,"goo_lambda",oo_id,oo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(lambda_t%gvv,"gvv_lambda",vv_id,vv_root,EXA_DATA_KIND_C8)

!        Get fixed intermediates
         call print_date('Calculating fixed intermediates')
         if (flavour.eq."EA".or.flavour.eq."EE") then
            call get_w_vvvv(lambda_t%w_vvvv,comm_t,int_t,CCD,one_tensor)
            call get_wbar_vovo(lambda_t%wbar_vovo,comm_t%t2,int_t,one_tensor)
         end if
         call get_fbar_ov(lambda_t%fbar_ov,comm_t%hov,one_tensor)
         if (flavour.eq."EA".or.flavour.eq."EE") then
            call get_w_vvvo(lambda_t,comm_t,int_t,CCD,one_tensor)
            call get_w_vovv(lambda_t%w_vovv,comm_t%t1,int_t,CCD,one_tensor)
         end if
         call get_w_vovo(lambda_t%w_vovo,comm_t,int_t%oovv,occ_space_id,occ_space_root, &
                         vir_space_id,vir_space_root,CCD,one_tensor)
         if (flavour.eq."IP".or.flavour.eq."EE") then
            call get_w_oooo(lambda_t%w_oooo,comm_t%a_int,one_tensor)
            call get_w_ovoo(lambda_t,comm_t,int_t,CCD,one_tensor)
            call get_w_ooov(lambda_t%w_ooov,comm_t%t1,int_t,CCD,one_tensor)
         end if
         call get_fbar_vv(lambda_t%fbar_vv,comm_t%gvv,one_tensor)
         call get_fbar_oo(lambda_t%fbar_oo,comm_t,CCD,one_tensor)
         if (flavour.eq."EE") then
            call get_gvv_lambda(lambda_t%gvv,l2_tensor,comm_t%t2)
            call get_goo_lambda(lambda_t%goo,l2_tensor,comm_t%t2)
         end if
         call print_date('Finished calculating fixed intermediates')

         if (flavour.eq."IP") then
            call eom_ip(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
         else if (flavour.eq."EE") then
            call eom_ee(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
         else if (flavour.eq."EA") then
            call eom_ea(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)
         else
            print*, "You should pick a valid EOM type."
         end if

         ierr=exatns_tensor_destroy(lambda_t%fbar_oo)
         ierr=exatns_tensor_destroy(lambda_t%fbar_ov)
         ierr=exatns_tensor_destroy(lambda_t%fbar_vv)
         ierr=exatns_tensor_destroy(lambda_t%w_oooo)
         ierr=exatns_tensor_destroy(lambda_t%w_ooov)
         ierr=exatns_tensor_destroy(lambda_t%w_ovoo)
         ierr=exatns_tensor_destroy(lambda_t%w_vovo)
         ierr=exatns_tensor_destroy(lambda_t%w_vvvo)
         ierr=exatns_tensor_destroy(lambda_t%w_vovv)
         ierr=exatns_tensor_destroy(lambda_t%w_vvvv)
         ierr=exatns_tensor_destroy(lambda_t%wbar_vovo)
         ierr=exatns_tensor_destroy(lambda_t%goo)
         ierr=exatns_tensor_destroy(lambda_t%gvv)

!        Make some noise so that we know we are leaving
         print*, " Leaving eom_driver routine"
         return

        end subroutine exacorr_eom_driver

        subroutine eom_ip(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor

         integer                    :: nocc, nvir
         integer, allocatable       :: mo_occ, mo_vir
         !logical                   :: CCD
         real(8)                    :: t_target_precision
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(3)       :: voo_id, oov_id
         integer(INTL),dimension(3) :: voo_root, oov_root

!        Solution tensors
         type(tens_rcrsv_t) :: r1_coef, r2_coef,l1_coef,l2_coef
         type(tens_rvrsv_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!        Auxiliary tensor
         type(tens_rcrsv_t) :: v_aux

         integer :: neom=30
         integer :: ierr, i

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         !CCD=exa_input%ccd

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root

!        Create all id and root arrays
         voo_id = (/vir_space_id, occ_space_id, occ_space_id/)
         oov_id = (/occ_space_id, occ_space_id, vir_space_id/)
         voo_root = (/vir_space_root, occ_space_root, occ_space_root/)
         oov_root = (/occ_space_root, occ_space_root, vir_space_root/)

!        Create all coefficient and sigma tensors
         ierr=exatns_tensor_create(r1_coef,"r1_coef",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_coef,"l1_coef",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_coef,"r2_coef",voo_id,voo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_coef,"l2_coef",oov_id,oov_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r1_sigma,"r1_sigma",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_sigma,"l1_sigma",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_sigma,"r2_sigma",voo_id,voo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_sigma,"l2_sigma",oov_id,oov_root,EXA_DATA_KIND_C8)

!        Create aux tensor
         ierr=exatns_create(v_aux,"v_aux",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)

         do i = 1, neom

!         r1_sigma
          ierr=exatns_tensor_init(r1_sigma)
          ierr=exatns_tensor_contract("S(i)+=F(m,i)*R(m)",r1_sigma,lambda_t%fbar_oo,r1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i)+=F(m,e)*R(e,i,m)",r1_sigma,lambda_t%fbar_ov,r2_coef)
          ierr=exatns_tensor_contract("S(i)+=W(m,n,i,e)*R(e,m,n)",r1_sigma,lambda_t%w_ooov,r2_coef,ONE_HALF)

!         r2_sigma
          ierr=exatns_tensor_init(r2_sigma)
          ierr=exatns_tensor_contract("S(a,i,j)+=W(m,a,j,i)*R(m)",r2_sigma,lambda_t%w_ovoo,r1_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_create("S(a,i,j)+=F(m,i)*r(a,j,m)",r2_sigma,lambda_t%fbar_oo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=F(m,j)*r(a,m,i)",r2_sigma,lambda_t%fbar_oo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=W(m,n,i,j)*r(a,m,n)",r2_sigma,lambda_t%w_oooo,r2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_create("S(a,i,j)+=W(a,m,e,i)*R(e,j,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=W(a,m,e,j)*R(e,m,i)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,i,j)+=F(a,e)*R(e,i,j)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_init(v_aux)
          !ierr=exatns_tensor_create("U(e)+=W(m,n,e,f)*R(f,m,n)",v_aux,w_oovv,r2_coef,ONE_HALF) ==> NO W_OOVV DEFINED ==> V_OOVV??
          ierr=exatns_tensor_create("S(a,i,j)+=U(e)*T(e,a,i,j)",r2_sigma,v_aux,comm_t%t2)

!         l1_sigma
          ierr=exatns_tensor_init(l1_sigma)
          ierr=exatns_tensor_contract("S(i)+=F(i,m)*L(m)",l1_sigma,lambda_t%fbar_oo,l1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(i)+=W(i,e,m,n)*L(m,n,e)",l1_sigma,lambda_t%w_ovoo,l2_coef,ONE_HALF)

!         l2_sigma
          ierr=exatns_tensor_init(l2_sigma)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(j,i,m,a)*L(m)",l2_sigma,lambda_t%w_ooov,l1_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(i,e,m,n)*L(m,n,e)",l2_sigma,lambda_t%w_ovoo,l2_coef,ONE_HALF)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a)+=F(j,m)*L(m,i,a)",l2_sigma,lambda_t%fbar_oo,l2_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=F(i,m)*L(j,m,a)",l2_sigma,lambda_t%fbar_oo,l2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a)+=W(e,j,a,m)*L(m,i,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(e,i,a,m)*L(j,m,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,j,a)+=W(i,j,m,n)*L(m,n,a)",l2_sigma,lambda_t%w_oooo,l2_coef,ONE_HALF)
          ierr=exatns_tensor_init(v_aux)
          ierr=exatns_tensor_contract("U(e)+=L(m,n,f)*T(e,f,m,n)",v_aux,l2_coef,comm_t%t2,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,j,a)+=V(i,j,e,a)*U(e)",l2_sigma,int_t%oovv,v_aux)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,j,a)+=L(i)*F(j,a)",l2_sigma,l1_coef,lambda_t%fbar_ov)
          ierr=exatns_tensor_contract("S(i,j,a)+=L(j)*F(i,a)",l2_sigma,l1_coef,lambda_t%fbar_ov,MINUS_ONE)

         end do

         deallocate(mo_occ)
         deallocate(mo_vir)

         ierr=exatns_tensor_destroy(r1_coef)
         ierr=exatns_tensor_destroy(l1_coef)
         ierr=exatns_tensor_destroy(r2_coef)
         ierr=exatns_tensor_destroy(l2_coef)
         ierr=exatns_tensor_destroy(r1_sigma)
         ierr=exatns_tensor_destroy(l1_sigma)
         ierr=exatns_tensor_destroy(r2_sigma)
         ierr=exatns_tensor_destroy(l2_sigma)
         ierr=exatns_tensor_destroy(v_aux)

        end subroutine eom_ip

        subroutine eom_ea(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor

         integer                    :: nocc, nvir
         integer, allocatable       :: mo_occ, mo_vir
         !logical                   :: CCD
         real(8)                    :: t_target_precision
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(3)       :: vvo_id, ovv_id
         integer(INTL),dimension(3) :: vvo_root, ovv_root

!        Solution tensors
         type(tens_rcrsv_t) :: r1_coef, r2_coef,l1_coef,l2_coef
         type(tens_rvrsv_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!        Auxiliary tensor
         type(tens_rcrsv_t) :: o_aux

         integer :: neom=30
         integer :: ierr, i

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         !CCD=exa_input%ccd

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root

         !Create all id and root arrays, neccessary for exatns_tensor_create
         vvo_id = (/vir_space_id, vir_space_id, occ_space_id/)
         ovv_id = (/occ_space_id, vir_space_id, vir_space_id/)
         vvo_root = (/vir_space_root, vir_space_root, occ_space_root/)
         ovv_root = (/occ_space_root, vir_space_root, vir_space_root/)

!        Create all coefficient and sigma tensors
         ierr=exatns_tensor_create(r1_coef,"r1_coef",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_coef,"l1_coef",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_coef,"r2_coef",vvo_id,vvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_coef,"l2_coef",ovv_id,ovv_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r1_sigma,"r1_sigma",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l1_sigma,"l1_sigma",(/vir_space_id/),(/vir_space_root/),EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(r2_sigma,"r2_sigma",vvo_id,vvo_root,EXA_DATA_KIND_C8)
         ierr=exatns_tensor_create(l2_sigma,"l2_sigma",ovv_id,ovv_root,EXA_DATA_KIND_C8)

!        Create aux tensor
         ierr=exatns_create(o_aux,"o_aux",(/occ_space_id/),(/occ_space_root/),EXA_DATA_KIND_C8)

         do i = 1, neom

!         r1_sigma
          ierr=exatns_tensor_init(r1_sigma)
          ierr=exatns_tensor_contract("S(a)+=F(a,e)*R(e)",r1_sigma,lambda_t%fbar_vv,r1_coef)
          ierr=exatns_tensor_contract("S(a)+=F(m,e)*R(e,a,m)",r1_sigma,lambda_t%fbar_ov,r2_coef)
          ierr=exatns_tensor_contract("S(a)+=W(a,m,f,e)*R(e,f,n)",r1_sigma,lambda_t%w_vovv,r2_coef,ONE_HALF)

!         r2_sigma
          ierr=exatns_tensor_init(r2_sigma)
          ierr=exatns_tensor_contract("S(a,b,i)+=W(a,b,e,i)*R(e)",r2_sigma,lambda_t%w_vvvo,r1_coef)
!         Permute a and b, so term appears twice
          ierr=exatns_tensor_create("S(a,b,i)+=F(a,e)*r(e,b,i)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=F(b,e)*r(a,e,i)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=W(m,n,i,j)*r(a,m,n)",r2_sigma,lambda_t%w_oooo,r2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_create("S(a,b,i)+=W(a,m,e,i)*R(e,j,m)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=W(a,m,e,j)*R(e,m,i)",r2_sigma,lambda_t%w_vovo,r2_coef)
          ierr=exatns_tensor_create("S(a,b,i)+=F(a,e)*R(e,i,j)",r2_sigma,lambda_t%fbar_vv,r2_coef)
          ierr=exatns_tensor_init(v_aux)
          !ierr=exatns_tensor_create("U(m)+=W(m,n,e,f)*R(f,m,n)",v_aux,w_oovv,r2_coef,ONE_HALF) ==> NO W_OOVV DEFINED ==> V_OOVV??
          ierr=exatns_tensor_create("S(a,b,i)+=U(e)*T(e,a,i,j)",r2_sigma,v_aux,comm_t%t2)

!         l1_sigma
          ierr=exatns_tensor_init(l1_sigma)
          ierr=exatns_tensor_contract("S(a)+=F(i,m)*L(m)",l1_sigma,lambda_t%fbar_oo,l1_coef,MINUS_ONE)
          ierr=exatns_tensor_contract("S(a)+=W(i,e,m,n)*L(m,n,e)",l1_sigma,lambda_t%w_ovoo,l2_coef,ONE_HALF)

!         l2_sigma
          ierr=exatns_tensor_init(l2_sigma)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(j,i,m,a)*L(m)",l2_sigma,lambda_t%w_ooov,l1_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(i,e,m,n)*L(m,n,e)",l2_sigma,lambda_t%w_ovoo,l2_coef,ONE_HALF)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,a,b)+=F(j,m)*L(m,i,a)",l2_sigma,lambda_t%fbar_oo,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=F(i,m)*L(j,m,a)",l2_sigma,lambda_t%fbar_oo,l2_coef)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,a,b)+=W(e,j,a,m)*L(m,i,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(e,i,a,m)*L(j,m,e)",l2_sigma,lambda_t%w_vovo,l2_coef)
          ierr=exatns_tensor_contract("S(i,a,b)+=W(i,j,m,n)*L(m,n,a)",l2_sigma,lambda_t%w_oooo,l2_coef,ONE_HALF)
          ierr=exatns_tensor_init(v_aux)
          ierr=exatns_tensor_contract("U(e)+=L(m,n,f)*T(e,f,m,n)",v_aux,l2_coef,comm_t%t2,ONE_HALF)
          ierr=exatns_tensor_contract("S(i,a,b)+=V(i,j,e,a)*U(e)",l2_sigma,int_t%oovv,v_aux)
!         Permute i and j, so term appears twice
          ierr=exatns_tensor_contract("S(i,a,b)+=L(i)*F(j,a)",l2_sigma,l1_coef,lambda_t%fbar_ov)
          ierr=exatns_tensor_contract("S(i,a,b)+=L(j)*F(i,a)",l2_sigma,l1_coef,lambda_t%fbar_ov,MINUS_ONE)

         end do

         deallocate(mo_occ)
         deallocate(mo_vir)

         ierr=exatns_tensor_destroy(r1_coef)
         ierr=exatns_tensor_destroy(l1_coef)
         ierr=exatns_tensor_destroy(r2_coef)
         ierr=exatns_tensor_destroy(l2_coef)
         ierr=exatns_tensor_destroy(r1_sigma)
         ierr=exatns_tensor_destroy(l1_sigma)
         ierr=exatns_tensor_destroy(r2_sigma)
         ierr=exatns_tensor_destroy(l2_sigma)
         ierr=exatns_tensor_destroy(o_aux)

        end subroutine eom_ea

        subroutine eom_ee(exa_input,comm_t,lambda_t,int_t,l1_tensor,l2_tensor,dims,one_tensor)

         type(exacc_input), intent(in) :: exa_input
!        Common tensors
         type(exatns_comm_tens),intent(inout)   :: comm_t
         type(exatns_lambda_tens,intent(inout)  :: lambda_t
!        fixed 2-body tensors (two-electron integrals)
         type(exatns_intg_tens),intent(inout) :: int_t
!        Lambda tensors
         type(tens_rcrsv_t),intent(inout) :: l1_tensor, l2_tensor
!        space ids and roots
         type(space_dims),intent(in) :: dims
!        scalar tensor
         type(tens_rcrsv_t), intent(inout) :: one_tensor

         integer                    :: nocc, nvir
         integer, allocatable       :: mo_occ, mo_vir
         !logical                   :: CCD
         real(8)                    :: t_target_precision
         integer                    :: occ_space_id, vir_space_id
         integer(INTL)              :: occ_space_root, vir_space_root
         integer,dimension(2)       :: oo_id, ov_id, vo_id, vv_id
         integer,dimension(4)       :: oovv_id, vvoo_id
         integer(INTL),dimension(2) :: oo_root, ov_root, vo_root, vv_root
         integer(INTL),dimension(4) :: oovv_root, vvoo_root

!        Solution tensors
         type(tens_rcrsv_t) :: r1_coef, r2_coef,l1_coef,l2_coef
         type(tens_rvrsv_t) :: r1_sigma, r2_sigma, l1_sigma, l2_sigma
!        Auxiliary tensor
         type(tens_rcrsv_t) :: oo_aux, vv_aux

         integer :: neom=30
         integer :: ierr, i

         !DO NOT CHANGE THESE VARIABLES
         nocc = exa_input%nkr_occ
         nvir = exa_input%nkr_vir
         allocate(mo_occ(nocc))
         allocate(mo_vir(nvir))
         mo_occ(1:nocc) = exa_input%mokr_occ(1:exa_input%nkr_occ)
         mo_vir(1:nvir) = exa_input%mokr_vir(1:exa_input%nkr_vir)
         t_target_precision = exa_input%t_econv
         !CCD=exa_input%ccd

         occ_space_id=dims%occ_space_id
         vir_space_id=dims%vir_space_id
         occ_space_root=dims%occ_space_root
         vir_space_root=dims%vir_space_root

         !Create all id and root arrays, neccessary for exatns_tensor_create
         oo_id = (/occ_space_id, occ_space_id/)
         ov_id = (/occ_space_id, vir_space_id/)
         vo_id = (/vir_space_id, occ_space_id/)
         vv_id = (/vir_space_id, vir_space_id/)
         oo_root = (/occ_space_root, occ_space_root/)
         ov_root = (/occ_space_root, vir_space_root/)
         vo_root = (/vir_space_root, occ_space_root/)
         vv_root = (/vir_space_root, vir_space_root/)
         oovv_id = (/occ_space_id, occ_space_id, vir_space_id, vir_space_id/)
         vvoo_id = (/vir_space_id, vir_space_id, occ_space_id, occ_space_id/)
         oovv_root = (/occ_space_root, occ_space_root, vir_space_root, vir_space_root/)
         vvoo_root = (/vir_space_root, vir_space_root, occ_space_root, occ_space_root/)

         do i = 1, neom
         end do

         deallocate(mo_occ)
         deallocate(mo_vir)

        end subroutine eom_ee

end module exacorr_eom
